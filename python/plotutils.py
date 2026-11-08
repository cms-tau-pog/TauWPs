# Author: Izaak Neutelings (August 2025)
import os, re
from math import log, ceil
from utils import warning
import numpy as np
from histutils import normalize, getxlimits, getylimits
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True # to avoid conflict with argparse
from ROOT import gROOT, gStyle, gPad, TAttAxis, TCanvas, TLatex, TH1, TH2, TH3, TGraph, THStack, TLine,\
                 kFullDotLarge, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross,\
                 TAttMarker, TAttLine, TColor, kBlack, kGray, kWhite, kSolid # kRed, kBlue, kMagenta, kOrange
gROOT.SetBatch(True) # do not open GUI windows
import cmsstyle as CMS
CMS.setCMSStyle()

# COLORS
def darken(c,f=0.8):
  c = gROOT.GetColor(c)
  return TColor.GetColor(f*c.GetRed(),f*c.GetGreen(),f*c.GetBlue())
kRawFillColor = TColor.GetColor(245,245,255)
kRawLineColor = TColor.GetColor(230,220,255) #TColor.GetColor(230,230,255)
kMyDashed     = 11; gStyle.SetLineStyleString(kMyDashed,"40 15")
kDashDoubleDotted = 8
#lcolors = [ kBlack, kRed, kBlue, kMagenta, kOrange ]
lcolors = [ TColor.GetColor(c) for c in CMS.petroff_6 ] #[::-1]
mcolors = [ darken(c,0.85) for c in lcolors ]
mstyles = [ kFullDotLarge, kFullTriangleUp, kFullTriangleDown, kFullDiamond, kFullCross ]


# CUSTOM COLO BINARY PALETTE, see core/base/src/TColor.cxx
kMyTemperature = None # array of colors
def setBinaryColorPalette(ncolors=None,alpha=1):
  """Create custom red-white-blue color palette to display binary values
  like correlations between -1 and 0."""
  global kMyTemperature
  if ncolors is None: # set default
    ncolors = 100 if kMyTemperature is None else len(kMyTemperature)
  if kMyTemperature is None or ncolors!=len(kMyTemperature):
    red   = np.array([ 0.80, 0.90, 1.00, 0.60, 0.02, ], dtype='float64')
    green = np.array([ 0.20, 0.80, 1.00, 0.80, 0.20, ], dtype='float64')
    blue  = np.array([ 0.10, 0.60, 1.00, 0.90, 0.65, ], dtype='float64')
    stops = np.array([i/(len(red)-1.) for i in range(len(red))], dtype='float64')
    ifirst = TColor.CreateGradientColorTable(len(red), stops, red, green, blue, ncolors)
    kMyTemperature = np.array([ifirst+i for i in range(ncolors)], dtype='int32')
  gStyle.SetPalette(ncolors,kMyTemperature,alpha)
  return kMyTemperature
  

def setera(era,string=None,extratext="Simulation Internal"):
  if era is None and string:
    matches = re.findall(r"20[12][0-9]",string)
    era = matches[0] if matches else ""
    etag = f"_{era}"
  else:
    etag = ""
  CMS.cmsstyle.cms_lumi = (era or "")
  #CMS.cms_lumi = era
  #CMS.SetLumi(lumi=None,run=era)
  CMS.SetEnergy(13.6, unit='TeV')
  CMS.SetExtraText(extratext)
  return era, etag
  

def getcanvas2D(hist,verb=0,**kwargs):
  """Plot 2D correlation map for lower cuts on the x and y axis."""
  if verb>=1:
    print(f">>> getcanvas2D: {hist.GetName()!r}")
  cname   = kwargs.pop('cname',    'canvas'          )
  logz    = kwargs.pop('logz',     False             )
  xtitle  = kwargs.pop('xtitle',   None              )
  ytitle  = kwargs.pop('ytitle',   None              )
  ztitle  = kwargs.pop('ztitle',   None              )
  zoffset = kwargs.pop('zoffset', 1.02 if logz else 1.00 )
  xmin    = kwargs.pop('xmin', hist.GetXaxis().GetXmin() )
  xmax    = kwargs.pop('xmax', hist.GetXaxis().GetXmax() )
  ymin    = kwargs.pop('ymin', hist.GetYaxis().GetXmin() )
  ymax    = kwargs.pop('ymax', hist.GetYaxis().GetXmax() )
  zmin    = kwargs.pop('zmin',     None              )
  zmax    = kwargs.pop('zmax',     None              )
  cwid    = kwargs.pop('cwid',     1000              )
  chgt    = kwargs.pop('chgt',     800               )
  tsize   = kwargs.pop('tsize',    0.050             )
  lsize   = kwargs.pop('lsize',    0.82*tsize        )
  frame   = kwargs.pop('frame',    None              )
  logit   = kwargs.pop('logit',    ''                )
  
  # SET CANVAS
  obj_dict = { }
  if verb>=2:
    print(f">>> getcanvas2D: Creating canvas (xmin,ymin,xmax,ymax)=({xmin:.5g},{ymin:.5g},{xmax:.5g},{ymax:.5g})...")
  lmarg, rmarg = 0.13, 0.17
  bmarg, tmarg = 0.12, 0.073
  canvas = TCanvas(cname,cname,cwid,chgt)
  canvas.SetMargin(lmarg,rmarg,bmarg,tmarg) # LRBT
  canvas.SetTicks(1,1)
  canvas.SetLogz(logz)
  ###gStyle.SetPaintTextFormat(".1f") # for marker text
  obj_dict['canvas'] = canvas
  
  # SET AXES
  if not frame:
    frame = canvas.DrawFrame(xmin,ymin,xmax,ymax)
  frame.GetXaxis().SetTitleSize(tsize)
  frame.GetYaxis().SetTitleSize(tsize)
  frame.GetZaxis().SetTitleSize(tsize)
  frame.GetXaxis().SetLabelSize(lsize)
  frame.GetYaxis().SetLabelSize(lsize)
  frame.GetZaxis().SetLabelSize(lsize)
  frame.GetXaxis().SetTitleOffset(1.07)
  frame.GetYaxis().SetTitleOffset(1.23)
  frame.GetZaxis().SetTitleOffset(zoffset)
  frame.GetXaxis().SetLabelOffset(0.0044)
  frame.GetYaxis().SetLabelOffset(0.0025)
  frame.GetZaxis().SetLabelOffset(-0.003 if logz else 0.005)
  if logz:
    frame.GetZaxis().SetNoExponent(True)
    hist.GetZaxis().SetNoExponent(True)
  ###hist.GetXaxis().SetRangeUser(xmin,xmax)
  ###hist.GetYaxis().SetRangeUser(ymin,ymax)
  if zmin is not None:
    hist.SetMinimum(zmin)
  if zmax is not None:
    hist.SetMaximum(zmax)
  if xtitle is not None:
    frame.GetXaxis().SetTitle(xtitle)
  if ytitle is not None:
    frame.GetYaxis().SetTitle(ytitle)
  if ztitle is not None:
    frame.GetZaxis().SetTitle(ztitle)
    hist.GetZaxis().SetTitle(ztitle)
    TAttAxis.Copy(frame.GetZaxis(),hist.GetZaxis()) # copy all atrributes
  frame.Draw('AXIS')
  obj_dict['frame'] = frame
  
  # CUSTOM AXIS for logit
  if 'x' in logit:
    obj_dict['xaxis'] = addLogitAxis(canvas,axis='x',tsize=lsize,offset=0.3)
  if 'y' in logit:
    obj_dict['yaxis'] = addLogitAxis(canvas,axis='y',tsize=lsize,offset=0.2)
  
  return obj_dict
  

def getlegend(xleg,yleg,entries,objs=None,header=None,ncols=2,tsize=0.045,dobkg=False):
  if objs:
    nobjs = len(objs)
  elif isinstance(entries,list):
    nobjs = len(entries)
  else: # assume number
    nobjs = entries
  nrows  = ceil(nobjs/float(ncols)) + int(bool(header))
  hleg   = 0.052*nrows
  wleg   = 0.24*(0.8*ncols if ncols>=1 else 0)
  legend = CMS.cmsLeg(xleg,yleg,xleg+wleg,yleg-hleg,textSize=tsize,columns=ncols)
  if ncols>=2:
    legend.SetMargin(0.24)
    legend.SetColumnSeparation(0.05)
  if header:
    legend.SetHeader(header)
  if dobkg:
    legend.SetFillColorAlpha(kWhite,0.50)
  for i, obj in enumerate(objs):
    if hasattr(obj,'GetTitle'):
      entry = obj.GetTitle()
      style = 'l'
    else: # assume string
      entry = obj
      obj = 0 # dummy
      style = ''
    if entries and i<len(entries):
      entry = entries[i]
    legend.AddEntry(obj,entry,style)
  return legend
  

def drawtext(lines,font=42,align=13,x=None,y=None,tsize=0.046,color=kBlack,below=None):
  """Draw TLaTeX text in the corner."""
  latex = TLatex()
  latex.SetTextSize(tsize)
  latex.SetTextAlign(align if align>0 else 13 if below else 11)
  latex.SetTextFont(font)
  latex.SetTextColor(color)
  latex.SetNDC(True)
  L, R = gPad.GetLeftMargin(), gPad.GetRightMargin()
  T, B = gPad.GetTopMargin(),  gPad.GetBottomMargin()
  if isinstance(lines,str):
    lines = lines.split('\n')
  if align==0: # outside frame
    x_ = L
    y_ = (1.015-T)
  elif below:
    if y is None:
      y = 0.02 # extra space
    x_ = below.GetX1NDC()
    y_ = below.GetY1NDC()-y # extra space from legend
  else:
    if x is None:
      x = 0.06
    if y is None:
      y = 0.95
    x_ = L+(1-L-R)*x # convert frame to canvas coordinates
    y_ = B+(1-T-B)*y # convert frame to canvas coordinates - i*1.2*tsize
  for i, line in enumerate(lines):
    y__ = y_ - i*1.2*tsize
    if i==0 and line.count('{gen}')>=2: # small correction
      y__ += 0.017
    latex.DrawLatex(x_,y__,line)
  return latex
  

def gethistratio(numhist,denhist,name=None,tag="_ratio"):
  """Make ratio between two histograms."""
  if isinstance(denhist,THStack): # convert to hist
    denhist = denhist.GetStack().Last()
  if name is None:
    name = numhist.GetName()
  rathist = numhist.Clone(name+tag) # numerator
  rathist.SetBinErrorOption(denhist.kNormal)
  rathist.Divide(denhist) # numerator / denominator = Obs. / Exp.
  for ibin in range(0,denhist.GetXaxis().GetNbins()+2): # include under-/overflow
    if denhist.GetBinContent(ibin)==0==numhist.GetBinContent(ibin):
      rathist.SetBinContent(ibin,1) # set manually
  return rathist
  

def getgraphratio(numgraph,dengraph,incl_zero=True,y_force=None,verb=0):
  """Create a TGraph that is the ratio of numgraph to dengraph using Eval from dengraph."""
  if verb>=2:
    print(f">>> getgraphratio: Compute ratio between {numgraph.GetName()} and  {dengraph.GetName()}...")
  npoints = numgraph.GetN()
  ratgraph = TGraph(npoints)
  TAttLine.Copy(numgraph,ratgraph)
  TAttMarker.Copy(numgraph,ratgraph)
  ratgraph.SetName(f"ratio_{numgraph.GetName()}_over_{dengraph.GetName()}")
  ratgraph.SetTitle(f"Ratio {numgraph.GetTitle()} / {dengraph.GetTitle()}")
  
  # LOOP over numerator
  ig = 0
  for i in range(npoints):
    x = numgraph.GetPointX(i)
    y_num = numgraph.GetPointY(i)
    if y_force is not None: # force to user's number
      ratio = y_force
    else: # compute ratio numerically
      y_den = dengraph.Eval(x)
      if y_num==y_den and (y_den>0 or (y_den==0 and incl_zero)):
        ratio = 1
      elif y_den==0:
        return
      else:
        ratio = 1.0 if abs(y_num-y_den)/y_den<1e-5 else (y_num / y_den)
      if verb>=3:
        print(f">>> getgraphratio: (ig,x,y_num,y_den)=({ig},{x:.4g},{y_num:.4g},{y_den:.4g}) => ratio={ratio:.4g}")
    ratgraph.SetPoint(ig,x,ratio)
    ig += 1
  
  return ratgraph
  

def addLogitAxis(pad=gPad,frame=None,hidelabels=False,axis='x',tsize=None,offset=0.5):
  """Add custom transformed x axis for a logit histogram to the frame."""
  ticks = [ ]
  labels = [ ]

  # GET FRAME
  if frame is None:
    frames = [o for o in pad.GetListOfPrimitives() if 'hframe' in o.GetName() and hasattr(o,'GetXaxis')]
    if not frames:
      frames = [o for o in pad.GetListOfPrimitives() if hasattr(o,'GetXaxis')][:1] # assume first is frame
  else:
    frames = [frame]
  if not frames:
    warning(f"Could find a frame to make a custom logit axis...",pre="addLogitAxis: ")
    return ticks, labels
  oldPad = gPad # just in case user-given pad!=gPad
  pad.cd()
  pad.SetTickx(0)
  xmin, xmax, ltick = None, None, None
  for frame in frames:
    taxis = frame.GetXaxis() if axis=='x' else frame.GetYaxis()
    if xmin is None:
      xmin, xmax = getxlimits(taxis)
      ltick = taxis.GetTickLength()
      if tsize is None:
        tsize = 0.86*taxis.GetLabelSize()
    taxis.SetTickLength(0)
    taxis.SetLabelSize(0)
  
  # CREATE CUSTOM AXIS (TGAxis did not work well...)
  if axis=='x':
    xminmarg, xmaxmarg = pad.GetLeftMargin(), pad.GetRightMargin()
    yminmarg, ymaxmarg = pad.GetBottomMargin(), pad.GetTopMargin()
  else:
    xminmarg, xmaxmarg = pad.GetBottomMargin(), pad.GetTopMargin()
    yminmarg, ymaxmarg = pad.GetLeftMargin(), pad.GetRightMargin()
  for y in [yminmarg,1-ymaxmarg]:
    fwid = 1-xminmarg-xmaxmarg # frame width
    dy = ltick*(fwid+0.02)*(1 if y==yminmarg else -1)
    for rx in [0.0002,0.002,0.01,0.05,0.2,0.5,0.8,0.95,0.99,0.998,0.9998]:
      lx = log(rx/(1-rx)) # convert logit(x) back to x
      if lx<xmin or lx>xmax: continue # out of range
      x = xminmarg+fwid*(lx-xmin)/(xmax-xmin)
      line = TLine()
      line.SetLineWidth(1)
      line.SetLineColor(kBlack)
      if axis=='x':
        line.DrawLineNDC(x,y,x,y+dy)
      else:
        line.DrawLineNDC(y,x,y+dy,x)
      if y==yminmarg and not hidelabels:
        text = TLatex()
        text.SetTextSize(tsize)
        text.SetTextFont(42)
        text.SetTextColor(kBlack)
        text.SetNDC(True)
        if axis=='x':
          tx = x + (-0.008 if rx<0.001 else -0.010 if rx==0.002 else
                    0.008 if rx>0.999 else  0.010 if rx==0.998 else 0)
          text.SetTextAlign(23) # center top
          text.DrawLatex(tx,y-offset*ltick,f"{rx:.4g}")
        else:
          text.SetTextAlign(32) # right center
          text.DrawLatex(y-offset*ltick,x,f"{rx:.4g}")
        labels.append(text)
      ticks.append(line)
  oldPad.cd()

  return ticks, labels
  

def plotHistComparison(hists,fname,ratio=True,norm=True,verb=0,**kwargs):
  """Plot histogram comparisons."""
  
  # NORMALIZE
  legobjset = [hists] # objects (TH1, str, ...) for legend
  hists = [h for h in hists if isinstance(h,TH1)]
  if norm:
    normalize(hists)
  
  # USER SETTINGS
  xmin, xmax = hists[0].GetXaxis().GetXmin(), hists[0].GetXaxis().GetXmax()
  ytitle     = "A.U." if norm else "Reconstructed#kern[-0.3]{ }#tau_{h} candidates"
  xtitle     = kwargs.pop('xtitle',    "p_{T} [GeV]"  )
  ytitle     = kwargs.pop('ytitle',    ytitle         )
  rtitle     = kwargs.pop('rtitle',    "Ratio"        )
  logy       = kwargs.pop('logy',      False          )
  logx       = kwargs.pop('logx',      False          )
  logr       = kwargs.pop('logr',      False          )
  xmin       = kwargs.pop('xmin',      1e-4 if logx else xmin )
  xmax       = kwargs.pop('xmax',      xmax           )
  ymarg      = kwargs.pop('ymarg',     1.20           ) # padding inside plot
  #ybmarg     = kwargs.pop('ymarg',     1.1 if logy else None ) # padding inside plot
  logyrange  = kwargs.pop('logyrange', 5e5            ) # ymin = ymax/logyrange
  ymin, ymax = getylimits(hists,logyrange=logyrange,ymarg=ymarg,xmin=xmin,xmax=xmax,logy=logy)
  ymin       = kwargs.pop('ymin',      ymin           )
  ymax       = kwargs.pop('ymax',      ymax           )
  rmin       = kwargs.pop('rmin',      5e-2 if logr else 0 )
  rmax       = kwargs.pop('rmax',      2e0            )
  text       = kwargs.pop('text',      None           )
  bmarg      = kwargs.pop('bmarg',     0.305 if ratio else 0.145 )
  grid       = kwargs.pop('grid',      False          )
  rgrid      = kwargs.pop('rgrid',     grid           )
  cname      = kwargs.pop('cname',     'canvas'       )
  headers    = kwargs.pop('header',    None           )
  lwidths    = kwargs.pop('lwidth',    [2]            )
  lstyles    = kwargs.pop('lstyle',    [kSolid]       )
  entries    = kwargs.pop('entries',   [ ]            )
  legobjset  = kwargs.pop('legobjs',   legobjset      ) # list of list of objects for legend 
  xleg       = kwargs.pop('xleg',      0.32           )
  dy         = 0.05*max(2,1 if text==None else len(text) if isinstance(text,list) else text.count('\n'))
  ylegs      = kwargs.pop('yleg',      0.86-dy        )
  ncols      = kwargs.pop('ncols',     2              )
  ifill      = kwargs.pop('ifill',     [ ]            ) # index of histograms to fill
  icol0      = kwargs.pop('icol0',     0              ) # shift colors indices
  outdir     = kwargs.pop('outdir',    None           )
  tag        = kwargs.pop('tag',       ""             )
  exts       = kwargs.pop('ext',       ['.png']       ) #,'.pdf']
  #write      = kwargs.pop('write',     None           ) # write to ROOT file
  logit      = kwargs.pop('logit',     'logit' in tag ) # custom logit(x) axis & ticks
  
  # LIST OF OBJECTS to avoid garbage collection
  ticks    = [ ]
  legends  = [ ]
  rathists = [ ]
  
  # PREPARE CANVAS
  ###CMS.setCMSStyle(True)
  ###CMS.cmsGrid(True)
  if ratio:
    canvas = CMS.cmsDiCanvas(cname,xmin,xmax,ymin,ymax,rmin,rmax,xtitle,ytitle,rtitle,
                             extraSpace=0.03,iPos=0,square=True)
  else:
    canvas = CMS.cmsCanvas(cname,xmin,xmax,ymin,ymax,xtitle,ytitle,
                           extraSpace=0.03,iPos=0,square=False)
  
  # SET STYLE
  gStyle.SetErrorX(0)        # size of horizontal error bars
  gStyle.SetEndErrorSize(5)  # size of bars at the end of error bars
  if grid or rgrid:
    gStyle.SetGridColor(kGray) # grid color
    gStyle.SetGridWidth(1)     # grid width
    gStyle.SetGridStyle(3)     # grid style
  
  # PREPARE FRAME
  mainpad = canvas.cd(1) # = canvas for !ratio
  if not ratio:
    mainpad.SetBottomMargin(bmarg)
  mainpad.SetLogx(logx)
  mainpad.SetLogy(logy)
  if grid:
    mainpad.SetGrid(1,1)
    mainpad.RedrawAxis('G') # force grid
  if logit: # custom logit(x) axis & ticks
    ticks.append(addLogitAxis(mainpad,hidelabels=ratio))
  mframe = mainpad.GetPrimitive('hframe')
  if ratio:
    mframe.GetXaxis().SetLabelSize(0)
    mframe.GetYaxis().SetLabelOffset(0.007)
    mframe.GetXaxis().SetTitleOffset(1.02)
    mframe.GetYaxis().SetTitleOffset(1.04)
  else:
    mframe.GetXaxis().SetLabelOffset(0.005)
    mframe.GetYaxis().SetLabelOffset(0.006)
    mframe.GetXaxis().SetTitleOffset(1.00)
    mframe.GetYaxis().SetTitleOffset(1.14 if ymax<1e-2 else 1.10)
  if verb>=2:
    for obj in list(mainpad.GetListOfPrimitives()):
      print(obj,repr(obj),obj.GetTitle())
  
  # DRAW HISTOGRAMS
  for i, hist in enumerate(hists):
    lcolor = lcolors[(icol0+i)%len(lcolors)]
    lstyle = lstyles[i if i<len(lstyles) else -1]
    lwidth = lwidths[i if i<len(lwidths) else -1]
    if i in ifill: # fill (transparent)
      hist.SetFillColorAlpha(lcolor,0.05 if i==0 else 0.04)
    else: # no fill
      hist.SetFillColor(0)
    hist.SetLineWidth(lwidth)
    hist.SetLineStyle(lstyle)
    hist.SetLineColor(lcolor)
    hist.Draw('SAME HIST')
  
  # EXTRA TEXT
  if text:
    drawtext(text)
  
  # DRAW LEGEND(s)
  if legobjset and not isinstance(legobjset[0],list):
    legobjset = [legobjset]
  for i, legobjs in enumerate(legobjset):
    header = headers[i] if isinstance(headers,list) else headers
    if isinstance(ylegs,list):
      yleg = ylegs[i] 
    else:
      yleg = ylegs
      if i>0:
        yleg -= 0.02 + abs(legends[i-1].GetY2()-legends[i-1].GetY1())
    legend = getlegend(xleg,yleg,entries,objs=legobjs,header=header,ncols=ncols,tsize=0.045)
    legend.Draw()
    legends.append(legend)
  CMS.UpdatePad(mainpad)
  
  # DRAW RATIO
  if ratio:
    lowerpad = canvas.cd(2) # lower panel with ratio plot
    rframe = lowerpad.GetPrimitive('hframe')
    lowerpad.SetLogx(logx)
    lowerpad.SetLogy(logr)
    lowerpad.SetBottomMargin(bmarg)
    if rgrid:
      lowerpad.SetGrid(1,1)
      lowerpad.RedrawAxis('G') # force grid
    if logit: # custom logit(x) axis & ticks
      ticks.append(addLogitAxis(lowerpad))
    rframe.GetXaxis().SetLabelOffset(0.007)
    rframe.GetYaxis().SetLabelOffset(0.007)
    rframe.GetYaxis().CenterTitle(True)
    denhist = hists[0]
    for numhist in hists:
      if isinstance(hist,str):
        continue
      rathist = gethistratio(numhist,denhist)
      rathist.Draw('SAME HIST ][')
      rathists.append(rathist)
    CMS.UpdatePad(lowerpad)
  
  # WRITE to PNG
  pname = os.path.join(outdir,f"{fname}{tag}")
  for ext in exts: #+['.C']:
    canvas.SaveAs(pname+ext)
  
  canvas.Close()
  
