#! /usr/bin/env python3
# Author: Izaak Neutelings (November 2025)
# Description: This script studies the correlation between VSjet/mu/e discriminators.
#              It attempts to find WPs in 2D (e.g. VSjet vs. VSe)
# Instructions:
#   python3 plotTauROCCurves2D.py output/tauAlgoEffs_2024_DY.root -t _DY -v1
import os, sys
import numpy as np
from math import log10
sys.path.append("python")
from utils import loadconfig
from rocutils import addCommonROCArgs, getROCFiles
from histutils import getHistFromFile
from plotutils import setera, addLogitAxis, getcanvas2D, setBinaryColorPalette
from ROOT import gROOT, TCanvas, TBox, TText,\
                 kRed, kWhite, kOrange, kDashed
import cmsstyle as CMS
CONFIG = loadconfig() # settings like cuts, binning, labels, ...


def makeEffHist2D(hist,name=None,title=None,scale=100.,zmin=1e-9,verb=0):
  """
  Convert 2D histogram into a map of efficiency for applying
  lower cuts on x and y (i.e. sum over bins with x>=x_cut and y>=x_cut).
  Includes overflow bins, but not underflow.
  """
  if verb>=1:
    print(f">>> makeEffHist2D: {hist.GetName()!r}")
  if name==None:
    name = hist.GetName()+"_eff"
  if title==None:
    xtit  = hist.GetXaxis().GetTitle()
    ytit  = hist.GetYaxis().GetTitle()
    title = f"Efficiency;Lower cut on {xtit};Lower cut on {ytit};Efficiency [%]"
  h_eff = hist.Clone(name)
  h_eff.SetTitle(title)
  h_eff.Reset()
  nxbins = hist.GetXaxis().GetNbins()
  nybins = hist.GetYaxis().GetNbins()
  zlast  = hist.GetBinContent(nxbins+1,nybins+1)
  if zlast<=0: # force drawing zero/empty bin
    hist.SetBinContent(nxbins+1,nybins+1,zmin)
  for iy in range(nybins+1,0,-1): # from overflow (nybins+1) down to first bin (iy=1)
    for ix in range(nxbins+1,0,-1): # from overflow (nxbins+1) down to first bin (ix=1)
      zval = hist.GetBinContent(ix,iy)
      if ix<=nxbins: # add cumulative value from right neighbor
        zval += h_eff.GetBinContent(ix+1,iy)
      if iy<=nybins: # add cumulative value from top neighbor (this adds overcount)
        zval += h_eff.GetBinContent(ix,iy+1)
        if ix<=nxbins: # subtract top right diagonal neighbor (to remove overcount)
          zval -= h_eff.GetBinContent(ix+1,iy+1)
      h_eff.SetBinContent(ix,iy,zval) # add cumulative
      h_eff.SetBinError(ix,iy,1e10)
  ntot1 = h_eff.GetBinContent(1,1) # this should correspond to 100%
  ###ntot0 = hist.Integral(1,nxbins+1,1,nybins+1)
  ###assert ntot0==ntot1, f"{ntot0=:.7g} (from input hist) != {ntot1=:.7g} (from eff. map)" # consistency check
  if ntot1!=0: # convert to efficiency in units of %
    h_eff.Scale(scale/ntot1)
  ###print(h_eff.GetBinContent(nxbins+1,nybins+1),h_eff.GetBinContent(nxbins,nybins))
  return h_eff
  

def makeCorrHist2D(hist,name=None,title=None,rebin=-1,verb=0):
  """
  Convert 2D histogram into a scan of the correlation factor between raw scores.
  The scan is performed over lower cuts on these raw output scores.
  """
  if verb>=1:
    print(f">>> makeCorrHist2D: {hist.GetName()!r}, {rebin=}")
  if name==None:
    name = hist.GetName()+"_corr"
  if title==None:
    xtit  = hist.GetXaxis().GetTitle()
    ytit  = hist.GetYaxis().GetTitle()
    title = f"Correlation;Lower cut on {xtit};Lower cut on {ytit};Correlation [%]"
  if rebin>=2:
    #oldhist = hist
    hist = hist.Rebin2D(rebin,rebin,hist.GetName()+"_rebin")
    #print(f"makeCorrHist2D: {oldhist=!r} -> {hist=!r}")
  h_cor = hist.Clone(name)
  h_cor.Reset() # reset bin contents
  nxbins = hist.GetXaxis().GetNbins()
  nybins = hist.GetYaxis().GetNbins()
  if verb>=1:
    print(f">>> makeCorrHist2D: Computing correlations for nxbins x nybins = {nxbins} x {nybins}...")
  for ix in range(1,nxbins+1):
    hist.GetXaxis().SetRange(ix,nxbins+1) # limit the range to compute the correlation
    for iy in range(1,nybins+1):
      hist.GetYaxis().SetRange(iy,nybins+1) # limit the range to compute the correlation
      zval = hist.GetCorrelationFactor() # e.g., between VSe and VSjet
      h_cor.SetBinContent(ix,iy,zval)
    if verb>=2:
      print(f">>> makeCorrHist2D: Done ix / nxbins = {ix} / {nxbins} = {100.*ix/nxbins}%...")
  hist.GetXaxis().SetRange(1,nxbins+1) # reset range
  hist.GetYaxis().SetRange(1,nxbins+1) # reset range
  if verb>=1:
    print(f">>> makeCorrHist2D: Done computing correlations...")
  return h_cor
  

def drawconts(canvas, hist, zconts=None, wpconts=None, graphs=None, tsize=0.027, xfrac=0.12, verb=1):
  """Draw WP contour labels along y axis."""
  if verb>=1:
    print(f">>> drawconts: {hist=}, {wpconts=}, {xfrac=}")
  objects = [ ]
  
  # COMPUTE & DRAW CONTOURS
  if wpconts and not graphs:
    if verb>=1:
      print(f">>> drawconts: Creating contours...")
    graphs = [ ]
    canvas2 = TCanvas('canvas_tmp','tmp',500,500)
    wplevels = np.array(wpconts,dtype='float64')
    wphist = hist.Clone(hist.GetName()+'_wpconts')
    wphist.SetLineWidth(2)
    wphist.SetLineColor(kRed)
    wphist.SetContour(len(wpconts),wplevels)
    wphist.Draw('CONT Z LIST')
    canvas2.Update() # create list
    contours = gROOT.GetListOfSpecials().FindObject('contours')
    for level, level_list in zip(wpconts,contours):
      for graph in level_list:
        graph.SetTitle(str(level))
        graph.SetLineColor(kRed)
        graph.SetFillStyle(0)
        graph.SetLineWidth(2)
        graph = graph.Clone() # clone to avoid loss when list is deleted
        graphs.append(graph)
    canvas2.Close()
  
  # DRAW CONTOURS
  canvas.cd()
  for graph in graphs:
    graph.Draw('L SAME')
  
  # GET CONTOUR LEVELS of z axis
  canvas.Update() # ensure contours & palette
  if zconts is not None:
    nconts = len(zconts)
  else:
    nconts = hist.GetContour()
    zarray = np.zeros(nconts, dtype='float64')
    hist.GetContour(zarray) # fill zarray
    zconts = list(zarray)
  
  # Z VALUE COLOR help function
  palette = hist.GetListOfFunctions().FindObject('palette')
  def _getzcol(x,y):
    zval = hist.GetBinContent(hist.FindBin(x,y))
    if zval<zmin:
      return kWhite
    iz   = next((i for i, z in enumerate(zconts) if z>=zval),None)
    if iz is not None:
      if canvas.GetLogz():
        zlog = log10(zmin)+(log10(zmax)-log10(zmin))*(iz-0.5)/nconts
        zval = 10**zlog
      else:
        zval = zmin+(zmax-zmin)*(iz-0.5)/nconts
    return palette.GetValueColor(zval)
  
  # DRAW LABELS
  frame  = canvas.GetPrimitive('hframe')
  xmin   = frame.GetXaxis().GetXmin()
  xmax   = frame.GetXaxis().GetXmax()
  ###ymin   = frame.GetYaxis().GetXmin()
  ymax   = frame.GetYaxis().GetXmax()
  zmin   = hist.GetMinimum()
  zmax   = hist.GetMaximum()
  x      = xmin+(xmax-xmin)*xfrac
  dx     = 0.62*tsize*(canvas.GetUxmax()-canvas.GetUxmin()) # box width scale
  dy     = 0.002*(canvas.GetUymax()-canvas.GetUymin()) # box height scale (~line width)
  #graphs = graphs[::-1] # large to small contours
  graphs = sorted(graphs,key=lambda g: g.Eval(x)) # large to small contours
  for i, graph in enumerate(graphs):
    y   = graph.Eval(x)
    zwp = float(graph.GetTitle().replace('%',''))
    if zwp==int(zwp):
      zwp = int(zwp)
    zstr = f"{zwp}%"
    wstr = dx*(len(zstr)-0.6*zstr.count('.')+0.5*zstr.count('%'))
    graph.SetTitle(zstr)
    
    # AVOID OUT OF FRAME (scan to higher x values)
    getxrange = lambda x: ((x-0.8*wstr),(x+0.6*wstr))
    x1, x2 = getxrange(x)
    y1 = graph.Eval(x1)
    while xfrac<0.95 and y1>ymax and (y+5*dy)>ymax:
      xfrac += 0.5*wstr/(canvas.GetUxmax()-canvas.GetUxmin()) # small step
      x = xmin+(xmax-xmin)*xfrac # overwrite previous x to avoid out-of-frame
      y = graph.Eval(x)
      y1 = graph.Eval(x1)
      x1, x2 = getxrange(x)
    if verb>=1:
      print(f">>> drawconts: {i=:2} -> {repr(zstr):>7s}, (x,y)=({x:6.3f},{y:6.3f}), (x1,x2,y1)=({x1:6.3f},{x2:6.3f},{y1:6.3f}), "
            f"{xfrac=:8.6g}, {(y+5*dy)=:6.3f}, {ymax=:5.3g}")
    
    # REMOVE OVERLAP (scan to higher x values)
    while xfrac<0.95 and len(ographs:=[ 
            g for g in graphs[:i] if any(abs(g.Eval(x_)-graph.Eval(x_))<3*dy for x_ in [x1,x2])
          ])>=1:
      for ograph in ographs:
        if verb>=2:
          print(f">>> drawconts: overlap old: {x=:.4g}, {y=:.4g}, {xfrac=:.4g} between {zwp=} and {ograph.GetTitle()}")
        xfrac += 1.4*wstr/(canvas.GetUxmax()-canvas.GetUxmin())
        x = xmin+(xmax-xmin)*xfrac # overwrite previous x to avoid overlap
        y = graph.Eval(x)
        x1, x2 = getxrange(x)
        if verb>=2:
          print(f">>> drawconts: overlap new: {x=:.4g}, {y=:.4g}, {xfrac=:.4g}")
    
    # DRAW BOXES overlapping with line
    x1, x2 = (x-wstr/2), (x+wstr/2)
    if abs(graph.Eval(x1)-graph.Eval(x2))<8*dy: # horizontal
      for dy_ in [-dy,+dy]:
        box = TBox(x1,y,x2,y+dy_)
        box.SetLineWidth(0)
        box.SetFillColor(_getzcol(x,y+dy_)) #4 if dy_>0 else 6)
        box.Draw()
        objects.append(box)
    else: # vertical
      hstr = tsize*(canvas.GetUymax()-canvas.GetUymin())
      y1, y2 = (y-0.47*hstr), (y+0.47*hstr)
      for dx_ in [-dx,+dx]:
        box = TBox(x,y1,x+dx_,y2)
        box.SetLineWidth(0)
        box.SetFillColor(_getzcol(x+dx_,y)) #4 if dy_>0 else 6)
        box.Draw()
        objects.append(box)
    
    # DRAW CONTOUR LABEL with box
    text = TText(x,y,zstr)
    text.SetTextFont(42)
    text.SetTextAlign(22)
    text.SetTextColor(kRed)
    text.SetTextSize(tsize)
    text.Draw()
    objects.append(text)
  return graphs, objects
  

def plotFracHist2D(infile,algo,xdisc,ydisc,gm='real',prefix='qlty',outdir='plots/2d',tag='',
                   rebin=10,merge='',threshold=0.0004,graph_dict=None,verb=0,**kwargs):
  """Prepare & plot of fraction of candididates discriminants."""
  if isinstance(infile,dict):
    infile = infile[gm]
  h_raw = getHistFromFile(f"raw/{algo}/raw{algo}{xdisc}_{ydisc}_logit_{prefix}_{gm}",
                          infile,scale=100,merge=merge,threshold=threshold,verb=verb)
  if rebin: # reduce number of x bins by evenly rebinning
    h_raw.Rebin2D(rebin,rebin)
  xtitle = f"Raw {algo}{xdisc} score"
  ytitle = f"Raw {algo}{ydisc} score"
  if gm=='real':
    ztitle = "Fraction of #tau_{#lower[-0.2]{h}} candidates [%]"
  else:
    GM_dict = CONFIG['text']['GM']
    ztitle = f"Fraction of {GM_dict[gm]} mis-ID [%]"
  pname  = f"frac_{algo}{xdisc}-{ydisc}_{gm}{tag}"
  xmax   = 7.9
  ###if graph_dict:
  ###  graphs = [g for g in graph_dict.get('real',[]) if any(str(d) in g.GetTitle() for d in [60,90])]
  graphs = graph_dict.get(gm,None)
  kwargs2 = dict( # defauly key-word arguments
    extra_graphs=graphs,xtitle=xtitle,ytitle=ytitle,ztitle=ztitle,logz=True,
    xmin=-xmax,ymin=-xmax,xmax=xmax,ymax=xmax,zmin=2e-4,zmax=2.0
  )
  kwargs2.update(kwargs) # override with user settings
  return plotFracHist2DCanvas(h_raw,outdir=outdir,verb=verb,name=pname,**kwargs2)
  

def plotFracHist2DCanvas(hist,tag="",outdir='plots/2d',verb=0,**kwargs):
  """Plot 2D efficiency map."""
  if verb>=1:
    print(f">>> plotFracHist2DCanvas: {hist.GetName()!r}")
  exts    = kwargs.pop('ext',      ['.png','.pdf']  )
  pname   = kwargs.pop('name',     hist.GetName()   )
  zconts  = kwargs.pop('zconts',   None             )
  #graphs  = kwargs.pop('graphs',   None             )
  extra_graphs = kwargs.pop('extra_graphs',  None   )
  
  # SET CANVAS
  CMS.SetCMSPalette()
  obj_dict = getcanvas2D(hist,verb=verb,**kwargs)
  canvas = obj_dict['canvas']
  
  # SET CONTOUR LEVELS
  if zconts:
    levels = np.array(zconts,dtype='float64')
    hist.SetContour(len(zconts),levels)
  else:
    hist.SetContour(40)
  
  # DRAW
  hist.Draw('COLZ SAME')
  
  # ADD CONTOUR VALUE ALONG Y
  if extra_graphs:
    graphs, objects = drawconts(canvas,hist,graphs=extra_graphs,zconts=zconts)
    #objects.extend(graphs+objects) # prevent garbage collection
  
  #### DRAW EXTRA GRAPHS
  ###if extra_graphs:
  ###  for graph in extra_graphs:
  ###    graph.Draw('LSAME')
  
  # CUSTOM AXIS for logit
  lsize = obj_dict['frame'].GetXaxis().GetLabelSize()
  obj_dict['xaxis'] = addLogitAxis(canvas,axis='x',tsize=lsize,offset=0.3)
  obj_dict['yaxis'] = addLogitAxis(canvas,axis='y',tsize=lsize,offset=0.2)
  
  # WRITE
  CMS.CMS_lumi(canvas,0)
  #CMS.UpdatePad(canvas)
  canvas.Modified()
  canvas.Update()
  canvas.RedrawAxis()
  if tag:
    pname += tag
  if outdir:
    pname = os.path.join(outdir,pname)
  for ext in exts:
    canvas.SaveAs(pname+ext)
  canvas.Close()
  #return graphs
  

def plotEffHist2D(infile,algo,xdisc,ydisc,gm='real',prefix='qlty',outdir='plots/2d',tag='',
                  merge='xy',threshold=0.0004,graph_dict=None,verb=0,**kwargs):
  """Prepare & plot efficiency map."""
  if isinstance(infile,dict):
    infile = infile[gm]
  h_raw = getHistFromFile(f"raw/{algo}/raw{algo}{xdisc}_{ydisc}_logit_{prefix}_{gm}",
                          infile,merge=merge,threshold=threshold,verb=verb)
  h_eff = makeEffHist2D(h_raw) # convert raw score histogram to efficiency map by integrating
  xtitle = f"Lower cut on {algo}{xdisc}"
  ytitle = f"Lower cut on {algo}{ydisc}"
  if gm=='real':
    ztitle = "Efficiency of #tau_{#lower[-0.2]{h}} candidates [%]"
  else:
    GM_dict = CONFIG['text']['GM']
    ztitle = f"{GM_dict[gm]} mis-ID rate [%]"
  pname  = f"eff_{algo}{xdisc}-{ydisc}_{gm}{tag}"
  xmax   = 7.6
  ymax   = 7.9
  zmin   = 0
  logz   = False
  graphs = None
  if gm=='real':
    zconts  = [z for z in range(0,40,10)]+[z for z in range(40,100,5)]+[98,99,99.5,100]
    wpconts = [50,60,70,80,90,95,98,99,99.5]
  else: # backgrounds
    logz    = True
    zconts  = [0.001,0.01,0.1,1,2,4,10,20,30,50,80,90,100]
    wpconts = [0.01,1,2,4,10,20,50]
    zmin    = zconts[0]
    if graph_dict:
      graphs = [g for g in graph_dict.get('real',[]) if any(str(d) in g.GetTitle() for d in [60,90])]
  kwargs2 = dict( # defauly key-word arguments
    zconts=zconts,wpconts=wpconts,extra_graphs=graphs,
    xtitle=xtitle,ytitle=ytitle,ztitle=ztitle,logz=logz,
    xmin=-xmax,ymin=-ymax,xmax=xmax,ymax=ymax,zmin=zmin
  )
  kwargs2.update(kwargs) # override with user settings
  return plotEffHist2DCanvas(h_eff,outdir=outdir,verb=verb,name=pname,**kwargs2)
  

def plotEffHist2DCanvas(hist,tag="",outdir='plots/2d',verb=0,**kwargs):
  """Plot 2D efficiency map."""
  if verb>=1:
    print(f">>> plotEffHist2DCanvas: {hist.GetName()!r}")
  logz    = kwargs.get('logz',     False            )
  exts    = kwargs.pop('ext',      ['.png','.pdf']  )
  pname   = kwargs.pop('name',     hist.GetName()   )
  zconts  = kwargs.pop('zconts',   None             )
  wpconts = kwargs.pop('wpconts',  None             )
  extra_graphs = kwargs.pop('extra_graphs',  None   )
  kwargs.setdefault('ztitle', "Efficiency [%]" )
  kwargs.setdefault('zmin',   1e-3 if logz else 0 )
  kwargs.setdefault('zmax',   100 )
  kwargs.setdefault('zoffset', 1.0 if logz else 0.95 )
  
  # SET CANVAS
  CMS.SetCMSPalette()
  obj_dict = getcanvas2D(hist,verb=verb,**kwargs)
  canvas = obj_dict['canvas']
  
  # SET CONTOUR LEVELS
  if verb>=1:
    print(f">>> plotEffHist2DCanvas: Creating contours...")
  if zconts:
    levels = np.array(zconts,dtype='float64')
    hist.SetContour(len(zconts),levels)
  
  # DRAW
  hist.Draw('COLZ SAME')
  
  # ADD CONTOUR VALUE ALONG Y
  if zconts:
    graphs, objects = drawconts(canvas,hist,zconts=zconts,wpconts=wpconts)
    #objects.extend(graphs+objects) # prevent garbage collection
  
  # DRAW EXTRA GRAPHS
  if extra_graphs:
    for graph in extra_graphs:
      graph.SetLineColor(kOrange+1)
      graph.SetLineStyle(kDashed)
      graph.Draw('L SAME')
  
  # CUSTOM AXIS for logit
  lsize = obj_dict['frame'].GetXaxis().GetLabelSize()
  obj_dict['xaxis'] = addLogitAxis(canvas,axis='x',tsize=lsize,offset=0.3)
  obj_dict['yaxis'] = addLogitAxis(canvas,axis='y',tsize=lsize,offset=0.2)
  
  # WRITE
  CMS.CMS_lumi(canvas,0)
  #CMS.UpdatePad(canvas)
  canvas.Modified()
  canvas.Update()
  canvas.RedrawAxis()
  if tag:
    pname += tag
  if outdir:
    pname = os.path.join(outdir,pname)
  for ext in exts:
    canvas.SaveAs(pname+ext)
  canvas.Close()
  return graphs
  

def plotCorrHist2D(infile,algo,xdisc,ydisc,gm='real',prefix='qlty',outdir='plots/2d',tag='',
                   merge='',threshold=0.0004,graph_dict=None,verb=0,**kwargs):
  """Prepare & plot efficiency map."""
  if isinstance(infile,dict):
    infile = infile[gm]
  h_raw = getHistFromFile(f"raw/{algo}/raw{algo}{xdisc}_{ydisc}_logit_{prefix}_{gm}",
                          infile,merge=merge,threshold=threshold,verb=verb)
  h_cor = makeCorrHist2D(h_raw,rebin=20) # convert raw score histogram to efficiency map by integrating
  xtitle = f"Lower cut on {algo}{xdisc}"
  ytitle = f"Lower cut on {algo}{ydisc}"
  if gm=='real':
    ztitle = "Correlation for #tau_{#lower[-0.2]{h}} candidates"
  else:
    gm_dict = CONFIG['text']['GM']
    ztitle = f"Correlation for {gm_dict[gm]} mis-ID"
  pname  = f"corr_{algo}{xdisc}-{ydisc}_{gm}{tag}"
  xmax   = 7.6
  ymax   = 7.9
  graphs = graph_dict.get(gm,None)
  kwargs2 = dict( # defauly key-word arguments
    xtitle=xtitle,ytitle=ytitle,ztitle=ztitle,extra_graphs=graphs,
    xmin=-xmax,ymin=-ymax,xmax=xmax,ymax=ymax,zmin=-1,zmax=+1
  )
  kwargs2.update(kwargs) # override with user settings
  return plotCorrHist2DCanvas(h_cor,outdir=outdir,verb=verb,name=pname,**kwargs2)
  

def plotCorrHist2DCanvas(hist,tag="",outdir='plots/2d',verb=0,**kwargs):
  """Plot 2D correlation map for lower cuts on the x and y axis."""
  if verb>=1:
    print(f">>> plotCorrHist2DCanvas: {hist.GetName()!r}")
  exts    = kwargs.pop('ext',      ['.png','.pdf']   )
  pname   = kwargs.pop('name',     hist.GetName()    )
  zconts  = kwargs.pop('zconts',   None              )
  extra_graphs = kwargs.pop('extra_graphs',  None   )
  
  # SET CANVAS
  setBinaryColorPalette()
  obj_dict = getcanvas2D(hist,verb=verb,**kwargs)
  canvas = obj_dict['canvas']
  
  # DRAW
  if verb>=1:
    print(f">>> plotCorrHist2DCanvas: Creating contours...")
  #hist.SetContour(40)
  hist.Draw('COLZ SAME')
  
  # ADD CONTOUR VALUE ALONG Y
  if extra_graphs:
    graphs, objects = drawconts(canvas,hist,graphs=extra_graphs,zconts=zconts)
    #objects.extend(graphs+objects) # prevent garbage collection
  
  # CUSTOM AXIS for logit
  lsize = obj_dict['frame'].GetXaxis().GetLabelSize()
  obj_dict['xaxis'] = addLogitAxis(canvas,axis='x',tsize=lsize,offset=0.3)
  obj_dict['yaxis'] = addLogitAxis(canvas,axis='y',tsize=lsize,offset=0.2)
  
  # WRITE
  CMS.CMS_lumi(canvas,0)
  canvas.Modified()
  canvas.Update()
  canvas.RedrawAxis()
  if tag:
    pname += tag
  if outdir:
    pname = os.path.join(outdir,pname)
  for ext in exts:
    canvas.SaveAs(pname+ext)
  canvas.Close()
  

def main(args):
  
  # GLOBAL SETTINGS from config.yaml
  genmatches    = CONFIG['genmatches']
  tauAlgos      = CONFIG['tauAlgos']
  
  # USER SETTINGS
  verbosity = args.verbosity
  infnames  = args.inputFile
  outdir    = args.outdir
  tag       = args.tag
  era       = args.era
  exts      = ['.png','.pdf']
  outdir    = os.path.join(outdir,'2d')
  os.makedirs(outdir,exist_ok=True)
  
  # INPUT FILES
  infile0, infiles = getROCFiles(infnames,tag=tag,genmatches=genmatches)
  
  # ERA/LUMI TEXT
  era, etag = setera(era,infnames[0])
  
  # EFFICIENCY CURVES
  xdisc = 'VSjet'
  ydisc = 'VSe'
  tags  = ['qlty'] # for DPS note
  tauAlgos = ['DeepTau2018v2p5']
  genmatches = ['real','jet','elec']
  genmatches = ['jet']
  #genmatches = ['elec']
  graph_dict = { } # dictionary of contour graphs per gm
  merge = 'xy' # reduce number of bins by merging (reduces PDF size)
  if 'nomerge' in tag:
    merge = ''
  for dtag in tags: # numerator
    for algo in tauAlgos:
      for gm in genmatches:
        prefix = dtag
        ptag   = f"{etag}_{dtag}{tag}"
        kwargs = dict( # common key=word arguments
          gm=gm,prefix=prefix,tag=ptag,graph_dict=graph_dict,
          ext=exts,outdir=outdir,verb=verbosity
        )
        graphs = plotEffHist2D(infiles,algo,xdisc,ydisc,merge=merge,**kwargs)
        graph_dict[gm] = graphs # contours
        plotFracHist2D(infiles,algo,xdisc,ydisc,merge='',**kwargs)
        plotCorrHist2D(infiles,algo,xdisc,ydisc,merge='',**kwargs)
  
  # CLOSE
  infile0.Close()
  for infile in infiles.values():
    infile.Close()
  

if __name__ == '__main__':
  from argparse import ArgumentParser
  description = """This script computes & plots efficiencies and ROC curves of tau ID algorithms."""
  parser  = ArgumentParser(description=description,epilog="Good luck!")
  methods = ['all','raw','qConf','eff','roc','comp']
  addCommonROCArgs(parser, outdir='plots', methods=methods)
  args = parser.parse_args()
  main(args)
  print(">>> Done!")
  
