#! /usr/bin/env python3
# Author: Izaak Neutelings (August 2025)
# Description: Plot the ROC curves from output from makeTauHists.py
# Instructions:
#   python3 plotTauROCCurves.py output/tauAlgoEffs_2024_DY.root -t _DY -v1
#   python3 plotTauROCCurves.py output/tauAlgoEffs_2024_HTT.root fake=output/tauAlgoEffs_2024_DY.root -t _Hto2Tau-DY -m roc -v1
import os, sys, re
import numpy as np
from math import log10
sys.path.append("python")
from utils import loadconfig, warning
from histutils import getHistFromFile, getHistSliceFromFile, getbinrange, divideByBinWidths
from plotutils import setera, getgraphratio, plotHistComparison, drawtext, addLogitAxis,\
                      kRawLineColor, kMyDashed, lcolors, mcolors, mstyles
from rocutils import addCommonROCArgs, getROCFiles, getProcFromFileName, getwpcuts
from computeTauIDWPs import computeEff, getCutForEff, getMisIDForEff, makeEffGraph
from ROOT import gStyle, TGraph,\
                 kWhite, kGray, kRed, kBlue, kFullDotLarge, kFullDiamond, kDashed
import cmsstyle as CMS

CONFIG  = loadconfig() # settings like cuts, binning, labels, ...
REXP_vVS = re.compile(r"(\s+v[1-9\.]+)(VS[a-z]+)") # swap version number and VSjet/mu/e


class ROCCurve:
  """
  Container class for bundling related TGraph objects:
  1. ROC curve of identification efficiency (solid line)
  2. ROC curve of reconstruction efficiency (dashed line)
  3. WP points (solid circles)
  """
  def __init__(self,g_rocid,g_rocreco=None,corr_eff=1,verb=0):
    self.g_rocreco = g_rocreco # dashed
    self.g_rocid   = g_rocid   # solid
    self.g_rocwps  = None      # circles
    self.g_roc     = None      # g_rocreco + g_rocid
    self.g_ratreco = None      # ratio g_rocreco/den
    self.g_ratid   = None      # ratio g_rocid/den
    self.g_ratwps  = None      # ratio g_rocwps/den
    self.corr_eff  = corr_eff  # correct efficiency for common denominator
    
  def __getattr__(self, name):
    print(f"ROCCurve.__getattr__({name!r})")
    return getattr(self.g_rocid,name) # fallback
    
  def addWPs(self,wp2eff_dict,logy=True,verb=0):
    """Compute WPs (eff,midID) and create graph."""
    self.g_rocwps = TGraph(len(wp2eff_dict))
    for ig, (wp, eff) in enumerate(wp2eff_dict.items()):
      eff  *= 0.01*self.corr_eff # correct for common denominator
      misid = getMisIDForEff(eff,self.g_rocid,logy=logy,verb=verb)
      #misid = self.g_rocid.Eval(eff)
      self.g_rocwps.SetPoint(ig,eff,misid)
    return self.g_rocwps
    
  def draw(self,lcolor=None,mcolor=None,lstyle=1,i=None):
    """Draw all ROC curves."""
    if i is not None:
      lcolor = lcolors[i%len(lcolors)]
      mcolor = mcolors[i%len(mcolors)] # darker
      mstyle = mstyles[i%len(mstyles)]
      msize  = 1.0 if mstyle==kFullDotLarge else 2.0 if mstyle==kFullDiamond else 1.4
    if mcolor is None:
      mcolor = lcolor
      mstyle = kFullDotLarge
    if self.g_rocreco: # reconstruction ROC (dashed)
      self.g_rocreco.SetFillColor(0)
      self.g_rocreco.SetLineStyle(kMyDashed) # dashed
      self.g_rocreco.SetLineWidth(2)
      self.g_rocreco.SetLineColor(lcolor)
      self.g_rocreco.SetMarkerSize(0)
      if self.g_rocreco.GetN()>0:
        self.g_rocreco.Draw('SAME L')
    if self.g_rocid: # identification ROC (solid)
      self.g_rocid.SetFillColor(0)
      self.g_rocid.SetLineStyle(lstyle) # solid (default)
      self.g_rocid.SetLineWidth(2)
      self.g_rocid.SetLineColor(lcolor)
      self.g_rocid.SetMarkerSize(msize)
      self.g_rocid.SetMarkerColor(mcolor)
      self.g_rocid.SetMarkerStyle(mstyle)
      if self.g_rocid.GetN()>0:
        #self.g_rocid.Draw('SAME LP0') # line+marker
        self.g_rocid.Draw('SAME L') # line only
    if self.g_rocwps: # working points (circles)
      self.g_rocwps.SetLineWidth(0)
      self.g_rocwps.SetMarkerColor(mcolor)
      self.g_rocwps.SetMarkerStyle(mstyle)
      self.g_rocwps.SetMarkerSize(msize)
      if self.g_rocwps.GetN()>0:
        self.g_rocwps.Draw('SAME P') # marker only
    
  def drawratio(self,g_den):
    """Compute ratio w.r.t. g_den and draw."""
    y_force = 1. if (g_den==self) else None # avoid noise deviating from ratio=1
    if isinstance(g_den,ROCCurve):
      g_den = g_den.getsinglegraph()
    if self.g_rocreco is not None:
      self.g_ratreco = getgraphratio(self.g_rocreco,g_den,y_force=y_force)
      self.g_ratreco.Draw('SAME L')
    self.g_ratid = getgraphratio(self.g_rocid,g_den,y_force=y_force)
    #self.g_ratid.Draw('SAME PL ][') # marker only
    self.g_ratid.Draw('SAME L') # line only
    if self.g_rocwps and self.g_ratid.GetN()>0:
      self.g_ratwp = self.g_rocwps.Clone(self.g_rocwps.GetName()+'_num')
      xmid = max(self.g_ratid.GetX()) # presumed point between ID & reco
      for i in range(self.g_ratwp.GetN()):
        x = self.g_ratwp.GetPointX(i) # efficiency
        if (x<=xmid or not self.g_ratreco):
          y = self.g_ratid.Eval(x) # misid
        else:
          y = self.g_ratreco.Eval(x) # misid
        self.g_ratwp.SetPoint(i,x,y)
      self.g_ratwp.Draw(f'SAME P ][')
    return (self.g_ratreco,self.g_ratid,self.g_ratwps)
    
  def getsinglegraph(self):
    """Merge g_rocreco and g_rocid into one TGraph."""
    if self.g_roc is None:
      self.g_roc = self.g_rocid.Clone(self.g_rocid.GetName()+'_full')
      if self.g_rocreco:
        for i, (x, y) in enumerate(zip(self.g_rocreco.GetX(),self.g_rocreco.GetY())):
          self.g_roc.InsertPointBefore(i,x,y) # preserve order: g_rocreco before g_rocid
    return self.g_roc
  

def makeROCGraphFromHists(h_real,h_fake,algo,gm_real='real',gm_fake='jet',den_real=None,den_fake=None,verb=0):
  """Prepare ROC curves as ROCCurve objects."""
  GM_dict = CONFIG['text']['GM']
  g_rocid = TGraph()
  g_rocid.SetName(f"roc_{algo}_{gm_fake}_vs_{gm_real}")
  g_rocid.SetTitle(f"{algo} ROC curve;{algo} efficiency;{GM_dict[gm_fake]} mis-ID rate")
  g_rocreco = g_rocid.Clone(g_rocid.GetName()+"_reco") if den_real else None
  iter = zip(computeEff(h_real,den=den_real,verb=verb),computeEff(h_fake,den=den_fake,verb=verb))
  x_last, y_last = None, None
  nbins = max(h_real.GetXaxis().GetNbins(),h_fake.GetXaxis().GetNbins())
  istep = 1 if (verb>=3 or nbins<50) else 20 if nbins>500 else 10 if nbins>200 else 5 if nbins>100 else 2
  for i, (real, fake) in enumerate(iter): # from 100% to 0%
    cut_real, eff_real, _ = real
    cut_fake, eff_fake, _ = fake
    if verb>=2 and (verb>=3 or i%istep==0 or i>=nbins):
      print(f">>> makeROCGraphFromHists: {i=}: (cut,eff)_real=({cut_real:.4g},{eff_real:.4g}), (cut,eff)_fake=({cut_fake:.4g},{eff_fake:.4g})")
    assert abs(cut_real-cut_fake)<1e-5, (f"Cuts on scores are inexpectedly different: cut_real={cut_real} and "
                                         f"cut_fake={cut_fake} for {h_real.GetName()!r} and {h_fake.GetName()!r}")
    if i<=1 and g_rocreco: # reconstruction ROC with common denominator (two points)
      g_rocreco.AddPoint(eff_real,eff_fake)
    if i>=1 or not g_rocreco: # identification ROC (after reconstruction ROC)
      if x_last is not None:
        if x_last==eff_real and y_last==eff_fake:
          continue # do not repeat previous point
        while x_last-eff_real>0.10: # check for large distance with previous eff.
          # insert linearly interpolated point to avoid wrong interpolation in ROC curve plot
          # (with logarithmic y scale) that mismatches the ratio panel (with linear y scale)
          x = x_last-0.10 # new point to insert
          y = y_last+(eff_fake-y_last)/(eff_real-x_last)*(x-x_last) # interpolate for x
          g_rocid.AddPoint(x,y) # insert (x,y)
          x_last, y_last = x, y
      g_rocid.AddPoint(eff_real,eff_fake)
      x_last, y_last = eff_real, eff_fake
  corr_eff = 1.
  if den_real:
    ntot = h_real.Integral(1,h_real.GetXaxis().GetNbins()+1) # include overflow
    corr_eff = ntot/den_real
  graph = ROCCurve(g_rocid,g_rocreco,corr_eff=corr_eff)
  return graph
  

def makeROCGraphsFromFiles(infiles,algos,disc,gm_real,gm_fake,prefix,wp2eff_dict=None,
                           wptag='',y=None,slice=None,tag='',verb=0):
  """Get raw output scores from input files and prepare ROC curves as ROCCurve objects."""
  
  # GET INPUT FILES
  assert gm_real in infiles, f"makeROCGraphsFromFiles: {gm_real=} not in {infiles}!"
  assert gm_fake in infiles, f"makeROCGraphsFromFiles: {gm_fake=} not in {infiles}!"
  infile_real = infiles[gm_real]
  infile_fake = infiles[gm_fake]
  if verb>=1:
    print(f">>> makeROCGraphsFromFiles: Constructing ROC graphs for "
          f"{gm_real}={infile_real.GetName()} and {gm_fake}={infile_fake.GetName()}")
  
  # GET COMMON DENOMINATOR
  den_real, den_fake = None, None
  if 'common' in tag or ('hps' in tag and 'hps3' not in tag):
    if 'hps3' in tag or 'all2' in tag:
      prefix_ = 'all2'
    elif 'hps' in tag:
      prefix_ = prefix.replace('logit_','')
    else:
      prefix_ = 'all'
    if slice and 'VS' not in y: # integrate ntot in pt/dm distribution
      h_real = getHistFromFile(f"{y}/{y}_{prefix_}_{gm_real}",infile_real,verb=verb)
      h_fake = getHistFromFile(f"{y}/{y}_{prefix_}_{gm_fake}",infile_fake,verb=verb)
      den_real = h_real.Integral(*getbinrange(h_real,*slice,verb=verb))
      den_fake = h_fake.Integral(*getbinrange(h_fake,*slice,verb=verb))
      ###print(f"{y}/{y}_{prefix_}_{gm_real}",f"{y}/{y}_{prefix_}_{gm_fake}",den_real,den_fake)
    else: # 1D histogram
      h_real = getHistFromFile(f"ntau/ntau_{prefix_}_{gm_real}",infile_real,verb=verb)
      h_fake = getHistFromFile(f"ntau/ntau_{prefix_}_{gm_fake}",infile_fake,verb=verb)
      den_real = h_real.GetBinContent(2)
      den_fake = h_fake.GetBinContent(2)
  elif slice and 'VS' in y and len(algos)==1: # normalize all ROCs to "No WP"
    algo   = algos[0]
    xvar   = f"raw{algo}" + (f"{y}_{disc}" if y=='VSjet' else f"{disc}_{y}")
    h_real = getHistFromFile(f"raw/{algo}/{xvar}_{prefix}_{gm_real}",infile_real,verb=verb)
    h_fake = getHistFromFile(f"raw/{algo}/{xvar}_{prefix}_{gm_fake}",infile_fake,verb=verb)
    den_real = h_real.Integral(0,h_real.GetXaxis().GetNbins(),0,h_real.GetYaxis().GetNbins())
    den_fake = h_fake.Integral(0,h_fake.GetXaxis().GetNbins(),0,h_fake.GetYaxis().GetNbins())
  
  # GET ROC CURVES by computing efficiencies from raw scores
  graphs = [ ]
  if 'VS' in y: # create ROCs per WP slice from TH2
    slice = tuple(('Reco' if wp=='None' else wp) for wp in slice) # map None/Reco
    if y=='VSjet':
      axis = 'x'   # VSjet is always x axis
      xdisc = y    # 'VSjet'
      ydisc = disc # 'VSe' or 'VSmu'
    else:
      axis = 'y'   # axis to be integrated over slice
      xdisc = disc # 'VSjet'
      ydisc = y    # 'VSe' or 'VSmu'
    algos_seen = set() # ensure unique elements, preserving order
    for algo in algos:
      if algo in algos_seen: # only do once
        continue
      else:
        algos_seen.add(algo)
      wp2cut = getwpcuts(algo,outdir='output/wps',tag=wptag,discriminants=[y],wps=slice,logit=True)
      assert wp2cut, f"Could not retrieve WP for {algo=}, {wptag=}. You may need to run computeTauIDWPs.py with the correct tag first!"
      for wp in slice:
        wp_    = 'Reco' if wp=='None' else wp
        wpcut  = wp2cut[y][wp_]['cut']
        slice_ = (wpcut,'inf') # lower cut (should align with bin edge)
        xvar   = f"raw{algo}{xdisc}_{ydisc}_{prefix}"
        h_real = getHistSliceFromFile(f"raw/{algo}/{xvar}_{gm_real}",infile_real,*slice_,axis=axis,verb=verb)
        h_fake = getHistSliceFromFile(f"raw/{algo}/{xvar}_{gm_fake}",infile_fake,*slice_,axis=axis,verb=verb)
        ###print(den_real,h_real.Integral(),h_real,repr(h_real))
        ###print(den_fake,h_fake.Integral(),h_fake,repr(h_fake))
        g_roc  = makeROCGraphFromHists(h_real,h_fake,algo,gm_real=gm_real,gm_fake=gm_fake,
                                       den_real=den_real,den_fake=den_fake,verb=verb)
        if wp2eff_dict:
          g_roc.addWPs(wp2eff_dict,logy=True,verb=verb)
        graphs.append(g_roc)
  else: # inclusive ROCs from TH1 or pT/DM-dependent ROCs from TH2D
    for algo in algos:
      if 'qConf' in disc:
        xvar = f"{disc}{algo}_{prefix}"
      else: # raw score
        xvar = f"raw{algo}{disc}_{prefix}"
      if slice: # get TH1 slice from TH2
        h_real = getHistSliceFromFile(f"{y}/{algo}/raw{algo}{disc}_{y}_{prefix}_{gm_real}",infile_real,*slice,axis='y',verb=verb)
        h_fake = getHistSliceFromFile(f"{y}/{algo}/raw{algo}{disc}_{y}_{prefix}_{gm_fake}",infile_fake,*slice,axis='y',verb=verb)
        #print(h_real.Integral(),h_fake.Integral())
      else: # get TH1 histogram from file
        hdir = f"charge/{algo}" if 'qConf' in disc else f"raw/{algo}"
        h_real = getHistFromFile(f"{hdir}/{xvar}_{gm_real}",infile_real)
        h_fake = getHistFromFile(f"{hdir}/{xvar}_{gm_fake}",infile_fake)
      g_roc  = makeROCGraphFromHists(h_real,h_fake,algo,gm_real=gm_real,gm_fake=gm_fake,
                                    den_real=den_real,den_fake=den_fake,verb=verb)
      if wp2eff_dict:
        g_roc.addWPs(wp2eff_dict,logy=True,verb=verb)
      graphs.append(g_roc)
  
  if not graphs:
    warning(f"Returning graphs...",pre='makeROCGraphsFromFiles: ')
  return graphs
  

def plotRawScores(infile,algo,disc,genmatches,prefix='qlty',norm=True,verb=0,**kwargs):
  """Plot raw output scores."""
  gm_dict   = CONFIG['text']['gm']
  
  # GET EFFICIENCIES
  hists = [ ]
  ifill = [0]  # fill 'real' histogram with transparent color
  xbins = None # variable binning shared between all histograms
  xvar  = f"{disc}{algo}" if 'qConf' in disc else f"raw{algo}{disc}"
  for i, gm in enumerate(genmatches):
    
    # GET HISTOGRAM
    hdir = f"charge/{algo}" if 'qConf' in disc else f"raw/{algo}"
    hist = getHistFromFile(f"{hdir}/{xvar}_{prefix}_{gm}",infile)
    if gm=='muon' and ('mu' in disc) or (gm=='jet' and 'jet' in disc) or\
      (gm=='elec' and 'e' in disc and 'jet' not in disc):
      ifill.append(i) # fill target background to highlight
    
    # REBIN
    if norm: # only rebin if also normalized (use non-normalized for debugging only)
      if 'logit' in prefix and hist.GetXaxis().GetNbins()>=800:
        hist = hist.Rebin(4)
      elif 'raw' in xvar: # rebin raw output score in the middle
        bw = hist.GetXaxis().GetBinWidth(1)
        if xbins is None: # initialize only once
          e = 0.01*bw
          xbins = np.concatenate((np.arange(0.0,0.2-e,bw),np.arange(0.2,0.4-e,2*bw),np.arange(0.4,0.6-e,4*bw),
                                  np.arange(0.6,0.8-e,2*bw),np.arange(0.8,1.0+e,bw)))
        hist = hist.Rebin(len(xbins)-1,hist.GetName()+"_rebin",xbins)
        divideByBinWidths(hist)
    
    hists.append(hist)
  
  # PLOT HISTOGRAM COMPARISON
  kwargs2 = dict( # defauly key-word arguments
    logy=True, logr=True, rmin=5e-2, rmax=2e1, xleg=0.32,
    ymarg=(1.38 if 'Deep' in algo else 1.32 if 'VSjet' in disc else 1.24),
    logyrange=(5e6 if 'PNetCHSVSmu' in xvar else 5e5),
    xtitle=f"{algo} charge assignment probability" if 'qConf' in xvar else f"Raw {algo}{disc} score",
    ytitle="A.U." if norm else "Reconstructed#kern[-0.3]{ }#tau_{h} candidates",
    entries=[gm_dict[gm] for gm in genmatches],
  )
  kwargs2.update(kwargs) # override with user settings
  plotHistComparison(hists,xvar,ifill=ifill,norm=norm,**kwargs2)
  

def plotEfficiencies(infile,algo,disc,gm="real",prefix='logit_qlty_',wp2eff_dict=None,verb=0,**kwargs):
  """Plot tau ID efficiencies."""
  
  # GET EFFICIENCIES
  xvar  = f"raw{algo}{disc}_{prefix.strip('_')}"
  h_raw = getHistFromFile(f"raw/{algo}/{xvar}_{gm}",infile)
  g_eff = makeEffGraph(h_raw,algo+disc,gm=gm,verb=verb)
  
  # USER SETTINGS
  xtitle     = f"Lower cut on raw{algo}{disc}"
  ytitle     = f"{algo}{disc} efficiency [%]"
  xmin, xmax = h_raw.GetXaxis().GetXmin(), h_raw.GetXaxis().GetXmax()
  ymin, ymax = 0, 110.
  outdir     = kwargs.pop('outdir',   ""       )
  xmax       = kwargs.pop('xmax',     xmax     )
  xmin       = kwargs.pop('xmin',     xmin     )
  logy       = kwargs.pop('logy',     False    )
  write      = kwargs.pop('write',    False    )
  cname      = kwargs.pop('cname',    'canvas' )
  tag        = kwargs.pop('tag',       ""      )
  exts       = kwargs.pop('ext',      ['.png'] ) #,'.pdf']:
  
  # PREPARE CANVAS
  #CMS.SetLumi(lumi)
  canvas = CMS.cmsCanvas(cname,xmin,xmax,ymin,ymax,xtitle,ytitle,square=False,
                         extraSpace=0.03,yTitOffset=1.05,iPos=0)
  gStyle.SetErrorX(0)       # size of horizontal error bars
  gStyle.SetEndErrorSize(5) # size of bars at the end of error bars
  #canvas.SetLeftMargin(0.11)
  canvas.SetLogy(logy)
  frame = canvas.GetPrimitive('hframe')
  frame.GetXaxis().SetTitleOffset(1.05)
  if (nchars:=len(ytitle))>=32:
    ysize = 32/nchars*frame.GetYaxis().GetTitleSize()
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitleSize(ysize)
    #frame.GetXaxis().SetTitleOffset(1.06)
  
  # DRAW RAW SCORE
  hmax = h_raw.GetMaximum()
  if hmax>0:
    for i in range(1,h_raw.GetXaxis().GetNbins()+1): # convert to log10
      y = h_raw.GetBinContent(i)
      h_raw.SetBinContent(i,log10(y) if y>0 else 0)
    h_raw.Scale(100.0/h_raw.GetMaximum())
    h_raw.SetFillColorAlpha(kRawLineColor,0.06)
    h_raw.SetLineColor(kRawLineColor)
    h_raw.Draw('SAME HIST')
  
  # DRAW EFFICIENCY
  g_eff.SetLineColor(kRed)
  g_eff.SetMarkerColor(kRed+1)
  g_eff.SetMarkerStyle(kFullDotLarge) # kFullDotLarge=8, scaleable dot
  g_eff.SetMarkerSize(0.25)
  #g_eff.Draw('SAME L P0') # no errors
  g_eff.Draw('SAME L PE0')
  
  # ADD WPs
  # https://root.cern.ch/doc/master/classTAttMarker.html#ATTMARKER21
  objs_wp = [ ] # to keep objects in memory
  if wp2eff_dict and gm=='real':
    g_wps = TGraph(len(wp2eff_dict))
    g_wps.SetLineWidth(0)
    g_wps.SetMarkerColor(kBlue)
    g_wps.SetMarkerStyle(51) # thick ROOT.kStar
    g_wps.SetMarkerSize(1)
    for i, (wp, eff) in enumerate(wp2eff_dict.items()):
      cut = getCutForEff(eff,g_eff,verb=verb)
      g_wps.SetPoint(i,cut,eff)
      g_line = TGraph(3,np.array([xmin,cut,cut],'d'),np.array([eff,eff,0],'d'))
      g_line.SetLineWidth(1)
      g_line.SetLineStyle(7)
      g_line.SetLineColor(kBlue-9)
      g_line.Draw('SAME L')
      objs_wp.append(g_line)
    g_wps.Draw('SAME P')
    objs_wp.append(g_wps)
  
  # CUSTOM AXIS for logit
  ticks = [ ]
  if 'logit' in prefix:
    ticks.append(addLogitAxis(canvas))
  
  # WRITE to PNG
  canvas.RedrawAxis()
  pname = os.path.join(outdir,f"eff_{algo}{disc}_{gm}{tag}")
  for ext in exts:
    canvas.SaveAs(pname+ext)
  
  # WRITE TO FILE
  if write:
    write.cd()
    gStyle.Write('style',TCanvas.kOverwrite)
    g_eff.Write(pname,TCanvas.kOverwrite)
    canvas.Write(cname,TCanvas.kOverwrite)
  canvas.Close()
  

def plotROCCurves(infiles,algos,disc,gm_real='real',gm_fake=None,prefix='logit',wp2eff_dict=None,
                  wptag='',y=None,slice=None,infile_comp=None,text=None,tag='',verb=0,**kwargs):
  """Prepare ROC curves as ROCCurve objects and plot comparison."""
  gm_dict   = CONFIG['text']['gm']
  GM_dict   = CONFIG['text']['GM']
  algo_dict = CONFIG['text']['algo']
  
  # DEFAULTS
  xtitle = "#tau_{h} ID efficiency (w.r.t. reconstruction)"
  if 'common' in tag:
    xtitle = "#tau_{h} efficiency (w.r.t. HPS+CHS+PUPPI)"
  elif 'hps3' in tag or 'hps3' in prefix:
    pass
  elif 'hps2' in tag or 'hps2' in prefix:
    xtitle = "#tau_{h} efficiency (w.r.t. HPS, excl. 2-prong)"
  elif 'hps' in tag or 'hps' in prefix:
    xtitle = "#tau_{h} efficiency (w.r.t. HPS)"
  
  # FAKE GENMATCH
  ymarg = 1.18 if 'common' in tag else 1.1
  xmin, ymin = 2e-1, 6e-5
  if gm_fake is None:
    if 'jet' in disc or 'qConf' in disc: # VSjet
      ymin = 6e-4
      ymarg = 1.05
      gm_fake = 'jet'
    elif 'mu' in disc: # VSmu
      if 'common' in tag:
        ymarg = 1.17 if len(text)<=1 else 1.25
        xmin = 0.4
      else:
        xmin = 0.90 if ('hps2' in tag and 'hps3' not in tag) else 0.98
      gm_fake = 'muon'
    elif 'e' in disc: # VSe
      if 'common' in tag:
        ymarg = 1.17 if len(text)<=1 else 1.25
        xmin = 0.30
      gm_fake = 'elec'
    else:
      raise IOError(f"Unknown {disc=!r} discriminator")
  if y and 'VS' in y:
    if 'mu' not in y:
      xmin = 0.2
    ymin  *= 0.1
    ymarg *= 1.2
  
  # GET GRAPHS
  graphs_comp = None
  graphs = makeROCGraphsFromFiles(infiles,algos,disc,gm_real,gm_fake,prefix,wp2eff_dict=wp2eff_dict,
                                  wptag=wptag,y=y,slice=slice,tag=tag,verb=verb)
  ###if infile_comp: # inputs to compare
  ###  graphs_comp = makeROCGraphsFromFiles(infile_comp,algos,disc,gm_real,gm_fake,prefix,wp2eff_dict=wp2eff_dict,
  ###                                      y=y,slice=slice,tag=tag,verb=verb)
  
  # PROCES TITLES
  
  xtext2, ytext2 = None, None
  ###if infiles[gm_real]==infiles[gm_fake]:
  ###  proctitles = [ f"{GM_dict[gm_real]} and {gm_dict[gm_fake]} from {getProcFromFileName(infiles[gm_real])}", ]
  ###else:
  proctitles = [
    f"{GM_dict[gm_real]} from {getProcFromFileName(infiles[gm_real])}",
    f"{GM_dict[gm_fake]} from {getProcFromFileName(infiles[gm_fake])}",
  ]
  if any(s in tag for s in ['all2_common','hps3_common']) and '_DY' in tag: # disable for DPS note
    proctitles = None
  if 'VS' in y:
    xtext2 = 0.49
    ytext2 = 0.17 # bottom right
  
  # PLOT ROC CURVES
  entries = [REXP_vVS.sub(r"\2\1",algo_dict.get(a,a)+disc) for a in algos]
  entries = kwargs.get('entries',entries)
  rtitle  = "Algo. / DeepTau  "
  if 'VS' in y and len(algos)==1:
    rtitle  = f"Algo. / No WP " #{entries[0]}"
  kwargs2 = dict( # defauly key-word arguments (update with optional user kwargs)
    xtitle=xtitle,xmin=xmin,ymin=ymin,
    ymarg=ymarg,
    ytitle=gm_dict[gm_fake].replace('jet','Jet')+" mis-ID rate",
    entries=entries,
    text=text,
    text2=proctitles, # extra text below legend
    xtext2=xtext2,ytext2=ytext2, # custom position
    dolegbkg=('VSmu' in disc),
    rtitle=rtitle,
  )
  kwargs2.update(kwargs) # override with user settings
  fname = f"roc_{disc}{tag}" #_{gm_real}-{gm_fake}
  plotROCCurvesCanvas(fname,graphs,graphs_comp,tag=tag,**kwargs2)
  

def plotROCCurvesCanvas(fname,graphs,graphs_comp=None,tag='',**kwargs):
  """Plot tau ID ROC curves as ROCCurve objects."""
  
  # USER SETTINGS
  xtitle   = kwargs.pop('xtitle',   "Efficiency"  )
  ytitle   = kwargs.pop('ytitle',   "Mis-ID rate" )
  rtitle   = kwargs.pop('rtitle',   "Ratio"       )
  logy     = kwargs.pop('logy',     True          )
  logx     = kwargs.pop('logx',     False         )
  xmin     = kwargs.pop('xmin', 1e-4 if logx else 0 )
  xmax     = kwargs.pop('xmax',     1             )
  ymin     = kwargs.pop('ymin', 1e-1 if logy else 0 )
  ymax     = kwargs.pop('ymax',     None          )
  ymarg    = kwargs.pop('ymarg',    1.1           ) # padding inside plot
  rmin     = kwargs.pop('rmin',     0.0           )
  rmax     = kwargs.pop('rmax',     2.0           )
  text     = kwargs.pop('text',     None          ) # extra text in top corner inside frame
  text2    = kwargs.pop('text2',    None          ) # extra text below legend
  xtext2   = kwargs.pop('xtext2',   None          ) # x coordinate of extra text
  ytext2   = kwargs.pop('ytext2',   None          ) # x coordinate of extra text
  cname    = kwargs.pop('cname',    'canvas'      )
  header   = kwargs.pop('header',   None          )
  entries  = kwargs.pop('entries',  None          ) or [ ]
  xleg     = kwargs.pop('xleg',     0.20          )
  dy       = 0.05*max(2,1 if text==None else len(text) if isinstance(text,list) else text.count('\n'))
  yleg     = kwargs.pop('yleg',     0.84-dy       )
  dolegbkg = kwargs.pop('dolegbkg', True          )
  outdir   = kwargs.pop('outdir',   None          )
  exts     = kwargs.pop('ext',      ['.png']      ) #,'.pdf']
  if ymax is None:
    ymax = 1. # 100% efficiency
    if logy:
      span = ymarg*abs(log10(ymax/ymin))
      ymax = ymin*10.**(span)
    else:
      ymax = ymin+ymarg*(ymax-ymin)
  
  # PREPARE CANVAS
  canvas = CMS.cmsDiCanvas(cname,xmin,xmax,ymin,ymax,rmin,rmax,xtitle,ytitle,rtitle,
                           extraSpace=0.03,iPos=0,square=True)
  gStyle.SetErrorX(0)       # size of horizontal error bars
  gStyle.SetEndErrorSize(5) # size of bars at the end of error bars
  gStyle.SetGridColor(kGray)
  
  # EXTRA TEXT
  if text:
    #H, W = gPad.GetWh()*gPad.GetHNDC(), gPad.GetWw()*gPad.GetWNDC()
    #gStyle.SetExponentOffset(-0.12*H/W,0.015,'y')
    drawtext(text)
  
  # PREPARE LEGEND
  nrows  = len(graphs) + int(bool(header))
  tsize  = 0.045
  hleg   = 0.052*nrows
  wleg   = 0.41 #if 'VSe' in fname
  legend = CMS.cmsLeg(xleg,yleg,xleg+wleg,yleg-hleg,textSize=tsize)
  legend.SetMargin(0.06/wleg)
  if header:
    legend.SetHeader(header)
  if dolegbkg:
    legend.SetFillColorAlpha(kWhite,0.60)
    legend.SetFillStyle(1001)
  
  # DRAW PREPARE STACK
  mainpad  = canvas.cd(1)
  lowerpad = canvas.GetPad(2)
  mainpad.SetLogx(logx)
  mainpad.SetLogy(logy)
  mainpad.SetGrid(1,1)
  mframe = mainpad.GetPrimitive('hframe')
  rframe = lowerpad.GetPrimitive('hframe')
  rframe.GetXaxis().SetTitleOffset(1.02)
  mframe.GetYaxis().SetTitleOffset(1.07)
  rframe.GetYaxis().CenterTitle(True)
  if len(rtitle)>14:
    rsize = rframe.GetXaxis().GetTitleSize()
    rframe.GetYaxis().SetTitleSize(0.95*rsize)
    rframe.GetYaxis().SetTitleOffset(0.5)
  if not logx and abs(xmax-xmin)<=0.06:
    for pad in [mainpad,lowerpad]:
      for obj in pad.GetListOfPrimitives():
        if 'hframe' in obj.GetName():
          obj.GetXaxis().SetNdivisions(505,False)
  
  # DRAW ROC CURVES (solid)
  for i, graph in enumerate(graphs):
    graph.draw(i=i) # draw all graphs (ID, reco, WPs)
    entry = entries[i] if i<len(entries) else graph.g_rocid.GetTitle()
    legend.AddEntry(graph.g_rocid,entry,'lp')
  
  # DRAW EXTRA ROC CURVES (dashed)
  if graphs_comp:
    for i, graph in enumerate(graphs_comp):
      graph.draw(i=i,lstyle=kDashed) # draw all graphs (ID, reco, WPs)
  
  # DRAW LEGEND
  legend.Draw()
  CMS.UpdatePad(mainpad)
  
  # EXTRA TEXT
  if text2:
    if (xtext2,ytext2)!=(None,None): # custom
      drawtext(text2,x=xtext2,y=ytext2)
    else: # below legend
      drawtext(text2,below=legend)
  
  # DRAW RATIO
  lowerpad = canvas.cd(2) # lower panel
  lowerpad.SetLogx(logx)
  lowerpad.SetGrid(1,1)
  lowerpad.RedrawAxis('G') # force grid
  dengraph = graphs[0]
  for numgraph in graphs:
    numgraph.drawratio(dengraph) # draw all graphs (ID, reco, WPs)
  
  # DRAW RATIO of EXTRA ROC CURVES (dashed)
  if graphs_comp:
    for numgraph in graphs_comp:
      numgraph.drawratio(dengraph) # draw all graphs (ID, reco, WPs)
  
  # WRITE to PNG
  CMS.UpdatePad(lowerpad)
  pname = os.path.join(outdir,fname).replace('$TAG',tag)
  for ext in exts:
    canvas.SaveAs(pname+ext)
  canvas.Close()
  

def main(args):
  
  # GLOBAL SETTINGS from config.yaml
  genmatches    = CONFIG['genmatches']
  tauAlgos      = CONFIG['tauAlgos']
  discriminants = CONFIG['discriminants']
  algo_dict     = CONFIG['text']['algo']
  pfix_dict     = CONFIG['text']['selection']
  wp2eff_dicts  = CONFIG['wp2eff_dicts']
  
  # USER SETTINGS
  verbosity = args.verbosity
  methods   = args.methods
  infnames  = args.inputFile
  outdir    = args.outdir
  tag       = args.tag
  era       = args.era
  exts      = ['.png','.pdf']
  effdir    = os.path.join(outdir,'eff')
  rawdir    = os.path.join(outdir,'raw')
  rocdir    = os.path.join(outdir,'roc')
  wpdir     = os.path.join(outdir,'wp')
  for dir in [outdir,effdir,rawdir,rocdir,wpdir]:
    if dir:
      os.makedirs(dir,exist_ok=True)
  gm_real = 'real'
  
  # INPUT FILES
  infile0, infiles = getROCFiles(infnames,tag=tag,genmatches=genmatches)
  #infile_comp = TFile(infnames[1],'READ') if len(infnames)>=2 else None
  
  # ERA/LUMI TEXT
  era, etag = setera(era,infnames[0])
  
  
  ########################
  #   RAW SCORE CURVES   #
  ########################
  
  if 'all' in methods or 'raw' in methods:
    prefixes = ['hps3',] #'qlty']
    for prefix in prefixes:
      for algo in tauAlgos:
        for disc in discriminants: # 'VSjet', 'VSmu', 'VSe'
          tkey = 'qlty_hps' if 'DeepTau' in algo else prefix if 'hps' in prefix else 'qlty_jet'
          kwargs = dict(ymin=5e-7,ymax=60,ext=exts,outdir=rawdir,text=pfix_dict[tkey],verb=verbosity)
          if prefix=='qlty':
            plotRawScores(infile0,algo,disc,genmatches,prefix=prefix,norm=False,tag=f"{etag}_{prefix}_abs{tag}",**kwargs)
          plotRawScores(infile0,algo,disc,genmatches,prefix=prefix,norm=True,tag=f"{etag}_{prefix}{tag}",**kwargs)
          plotRawScores(infile0,algo,disc,genmatches,prefix='logit_'+prefix,norm=True,tag=f"{etag}_{prefix}_logit{tag}",
                        xmin=-9,xmax=9,**kwargs)
  
  
  #############
  #   QCONF   #
  #############
  
  if 'all' in methods or 'qConf' in methods:
    prefixes = ['all','hps','hps2','hps3','hps_hps3','rec_hps3']
    for prefix in prefixes:
      for algo in tauAlgos:
        if 'DeepTau' in algo:
          continue
        text = pfix_dict[prefix]
        text[-1] = text[-1].replace("/ |qConf| > 0.2","").replace(", |qConf| > 0.2","")
        kwargs = dict(xmin=-0.5,xmax=0.5,ymin=6e-5,ymax=2.5,logyrange=5e3,text=text,
                      prefix=prefix,ext=exts,outdir=rawdir,verb=verbosity)
        ###plotRawScores(infile0,algo,'qConf',genmatches,norm=False,tag=etag+"_{prefix}"+tag,**kwargs)
        plotRawScores(infile0,algo,'qConf',genmatches,norm=True,tag=etag+f"_{prefix}_norm"+tag,**kwargs)
  
  
  #########################
  #   EFFICIENCY CURVES   #
  #########################
  
  if 'all' in methods or 'eff' in methods:
    prefixes = ['hps3',] #'qlty']
    for prefix in prefixes:
      for algo in tauAlgos:
        for disc in discriminants:
          wp2eff_dict = wp2eff_dicts.get(disc,None)
          tkey = 'qlty_hps' if 'DeepTau' in algo else prefix if 'hps' in prefix else 'qlty_jet'
          kwargs = dict(gm=gm_real,wp2eff_dict=wp2eff_dict,text=pfix_dict[tkey],ext=exts,outdir=effdir,verb=verbosity)
          plotEfficiencies(infile0,algo,disc,prefix=prefix,tag=f"{etag}_{prefix}{tag}",**kwargs)
          plotEfficiencies(infile0,algo,disc,prefix='logit_'+prefix,tag=f"{etag}_{prefix}_logit{tag}",xmin=-9,xmax=9,**kwargs)
  
  
  #################################
  #   ROC CURVES for raw scores   #
  #################################
  
  if 'all' in methods or 'roc' in methods:
    ###discriminants = ['VSmu']
    for disc in discriminants: # 'VSjet', 'VSmu', 'VSe'
      wp2eff_dict = wp2eff_dicts.get(disc,None)
      #prefixes = ['all2_common','hps3_common'] # for DPS note
      prefixes = ['all_common','all2_common','hps','hps2','hps3','hps3_common','qlty','qlty_common']
      for dtag in prefixes: # numerator
        if dtag in ['all2_common','hps3_common'] and '_DY' in tag:
          setera(None,extratext="Simulation Preliminary") # for DPS note
        if 'qlty' in dtag:
          tkey = 'qlty'
          #tkey = 'qlty_hps' if 'DeepTau' in algo else 'qlty_jet'
        else:
          tkey = (dtag if 'hps' in dtag else 'all').replace('_common','')
        prefix = 'logit_'+dtag.replace('_common','')
        ptag   = f"{etag}_{dtag}{tag}"
        texts  = pfix_dict[tkey]
        ###plotROCCurves(infiles,tauAlgos,disc,gm_real=gm_real,wp2eff_dict=wp2eff_dict,text=texts,
        ###              prefix=prefix,ext=exts,outdir=rocdir,tag=ptag,verb=verbosity)
        ###setera(era) # reset for internal documentation
        
        #### SPLIT BY DM or PT bin
        ###if dtag=='hps3_common' and disc=='VSjet':
        ###  
        ###  for dm in [0,1,10,11]:
        ###    ptag_ = f"_dm{dm}{ptag}"
        ###    texts_ = texts.copy()
        ###    texts_[0] += f", DM {dm_dict.get(dm,dm)}"
        ###    plotROCCurves(infiles,tauAlgos,disc,gm_real=gm_real,wp2eff_dict=wp2eff_dict,text=texts_,
        ###                  y='dm',slice=(dm,),prefix=prefix,ext=exts,outdir=rocdir,tag=ptag_,verb=verbosity)
        ###  
        ###  for ptbin in [(20,30),(30,50),(50,'inf')]:
        ###    texts_ = texts.copy()
        ###    ptag   = f"{ptbin[0]}to{ptbin[1]}"
        ###    ptag_  = f"_pt{ptag}{ptag}"
        ###    ptstr  = "p_{#lower[-0.15]{T}}"
        ###    if ptbin[-1]=='inf':
        ###      texts_[0] = f"{ptstr} > {ptbin[0]} GeV"
        ###    else:
        ###      texts_[0] = f"{ptbin[0]} < {ptstr}#kern[-.4]{{ }}#leq {ptbin[1]} GeV"
        ###    texts_[0] += ", |#eta| < 2.5, |d_{#lower[-0.15]{z}}| < 0.2 cm"
        ###    plotROCCurves(infiles,tauAlgos,disc,gm_real=gm_real,wp2eff_dict=wp2eff_dict,text=texts_,
        ###                  y='dm',slice=ptbin,prefix=prefix,ext=exts,outdir=rocdir,tag=ptag_,verb=verbosity)
        
        # FOR DIFFERENT WPs
        if disc in ['VSjet','VSe'] and 'qlty' in dtag: #any(s in dtag for s in ['qlty']):
          for wpdisc in ['VSe','VSmu','VSjet']:
          #for wpdisc in ['VSe','VSjet']:
            if disc==wpdisc or (disc=='VSe' and wpdisc=='VSmu'): continue
            wp0    = 'VLoose' if wpdisc=='VSmu' else 'VVVLoose' # loosest WP
            ###wps    = ['None','Loose','Tight']
            wps    = ['None',wp0,'Loose','Tight']
            wptag  = ptag.replace('_common',"") #prefix.replace('logit','')+tag
            
            # COMPARE ALGORITHM for a given WP
            if 'common' in dtag: # plot only with common denominator
              for wp in wps:
                wpstr  = 'No' if 'No' in wp else wp
                texts_ = texts+[f"{wpstr} {wpdisc} WP"]
                ptag_  = f"-{wpdisc}{wp}{ptag}"
                plotROCCurves(infiles,tauAlgos,disc,gm_real=gm_real,wp2eff_dict=wp2eff_dict,text=texts_,wptag=wptag,
                              y=wpdisc,slice=[wp],prefix=prefix,ext=exts,outdir=wpdir,tag=ptag_,verb=verbosity)
            
            # COMPARE WPs for a given tau algorithm
            for algo in tauAlgos:
              texts_  = [t for t in texts if ( # filter text
                          ('DeepTau' in algo or 'HPS' not in t) and
                          ('DeepTau' not in algo or 'qConf' not in t)
                        )] + [REXP_vVS.sub(r"\2\1",algo_dict.get(algo,algo)+disc)]
              ptag_   = f"-{wpdisc}_{algo}{ptag}"
              wptag   = ptag.replace('_common',"") #prefix.replace('logit','')+tag
              entries = [(f"No {wpdisc} WP" if 'No' in wp else f"{wp} {wpdisc}") for wp in wps]
              plotROCCurves(infiles,[algo],disc,gm_real=gm_real,wp2eff_dict=wp2eff_dict,text=texts_,entries=entries,
                            wptag=wptag,y=wpdisc,slice=wps,prefix=prefix,ext=exts,outdir=wpdir,tag=ptag_,verb=verbosity)
  
  
  #### ROC CURVES comparisons
  ###if ('all' in methods or 'comp' in methods) and infile_comp:
  ###  for disc in discriminants:
  ###    wp2eff_dict = wp2eff_dicts.get(disc,None)
  ###    tags = ['hps3','hps3_common']
  ###    for dtag in tags: # numerator
  ###      tkey   = dtag.replace('_common','')
  ###      prefix = 'logit_'+dtag.replace('_common','')
  ###      dtag   = f"{etag}_{dtag}{tag}"
  ###      texts  = pfix_dict[tkey]
  ###      plotROCCurves(infile0,tauAlgos,disc,gm_real=gm_real,wp2eff_dict=wp2eff_dict,text=texts,
  ###                    infile_comp=infile_comp,prefix=prefix,ext=exts,outdir=rocdir,tag=dtag,verb=verbosity)
  
  
  ############################
  #   ROC CURVES for qConf   #
  ############################
  
  if 'all' in methods or ('roc' in methods and 'qConf' in methods):
    jetAlgos = [a for a in tauAlgos if 'DeepTau' not in a]
    for dtag in ['all']:
      prefix = 'logit_'+dtag
      dtag   = f"{etag}_{dtag}{tag}"
      plotROCCurves(infiles,jetAlgos,'qConf',gm_real=gm_real,gm_fake='jet',text=pfix_dict['all'],
                    prefix=prefix,ext=exts,outdir=rocdir,tag=dtag,verb=verbosity)
  
  
  # CLOSE
  infile0.Close()
  for infile in infiles.values():
    infile.Close()
  

if __name__ == '__main__':
  from argparse import ArgumentParser, Action
  description = """This script computes & plots efficiencies and ROC curves of tau ID algorithms."""
  parser  = ArgumentParser(description=description,epilog="Good luck!")
  methods = ['all','raw','qConf','eff','roc','comp']
  addCommonROCArgs(parser, outdir='plots', methods=methods)
  args = parser.parse_args()
  main(args)
  print(">>> Done!")
  
