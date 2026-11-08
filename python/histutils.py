#! /usr/bin/env python3
# Author: Izaak Neutelings (August 2025)
import os
from math import log10
import numpy as np
from utils import warning
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True # to avoid conflict with argparse
from ROOT import gROOT, TNamed, TH2, TH3, TH2D #, TGraph
gROOT.SetBatch(True) # do not open GUI windows

# more meaningful representation for debugging:
TNamed.__repr__ = lambda s: f"<{s.__class__.__name__}({s.GetName()!r}) at {hex(id(s))}>"


def getHistFromFile(hname,infile,nonzero=True,norm=False,scale=None,merge='',threshold=0.0004,verb=0):
  """Help function to get nonzero histogram from file."""
  if verb>=1:
    print(f">>> getHistFromFile: {infile.GetPath()}{hname}")
  ntot = None
  hist = infile.Get(hname)
  if not hist:
    warning(f"Could not retrieve {hname} in {infile.GetPath()}...",pre='getHistFromFile: ')
    hdir = infile.Get(os.path.dirname(hname))
    if hdir:
      hdir.ls()
    else:
      infile.ls()
    return None
  hist.SetDirectory(0) # detach from file
  if nonzero: # assert nonzero histogram
    ntot = getintegral(hist,of=True)
    if ntot<=0:
      warning(f"{hname} has ntot={ntot}. Ignoring...",pre='getHistFromFile: ')
      return None
  if norm or (scale is not None):
    if ntot is None:
      ntot = getintegral(hist,of=True)
    if ntot==0:
      warning(f"{hname} has ntot={ntot}. Cannot normalize...",pre='getHistFromFile: ')
    else:
      scale = 1. if scale is None else float(scale)
      hist.Scale(scale/ntot)
  if 'x' in merge: # reduce number of x bins
    hist = mergebins(hist,threshold=threshold/len(merge),axis='x',verb=verb) # reduce number of bins
  if 'y' in merge: # reduce number of y bins
    hist = mergebins(hist,threshold=threshold,axis='y',verb=verb) # reduce number of bins
  return hist
  

def getHistSliceFromFile(hname,infile,xmin,xmax=None,axis='x',nonzero=True,verb=0):
  """Help function to get slice of nonzero histogram from file."""
  hist2d = getHistFromFile(hname,infile,nonzero=nonzero,verb=verb)
  if xmin is not None or xmax is not None:
    hist1d = getslice(hist2d,xmin,xmax=xmax,axis=axis,verb=verb)
  if verb+2>=2:
    print(f">>> getHistSliceFromFile: {hist2d.GetName()!r} (ntot={hist2d.Integral()}) "
          f"-> {hist1d.GetName()!r} (ntot={hist1d.Integral()}) for ({axis}min,{axis}max)=({xmin},{xmax})")
  return hist1d
  

def getintegral(hist,of=True):
  """Get integral. Include overflow if of==True"""
  assert hasattr(hist,'Integral'), f"{hist} does not have an Integral method..."
  if not of: # do not include overflow
    return hist.Integral()
  if isinstance(hist,TH3):
    return hist.Integral(1,hist.GetXaxis().GetNbins()+1,1,hist.GetYaxis().GetNbins()+1,
                         1,hist.GetZaxis().GetNbins()+1)
  if isinstance(hist,TH2):
    return hist.Integral(1,hist.GetXaxis().GetNbins()+1,1,hist.GetYaxis().GetNbins()+1)
  return hist.Integral(1,hist.GetXaxis().GetNbins()+1)
  

def normalize(hists,norm=1.,of=False):
  if not isinstance(hists,(list,tuple)):
    hists = [ hists ]
  for hist in hists:
    if hasattr(hist,'Integral'):
      sum = getintegral(hist,of=of)
      if sum!=0:
        hist.Scale(norm/sum)
  return hists
  

def getslice(hist2d,xmin,xmax=None,axis='x',name=None,verb=0):
  """
  Get TH1 slice from TH2, e.g. x=[xmin,xmax].
  Assume given xmin and xmax align with bin edges.
  """
  if isinstance(xmin,(list,tuple)):
    assert xmax is None, f"{xmin=} and {xmax=}, but expected numbers and/or xmax=None..."
    assert len(xmin)==2, f"{len(xmin)=}, but expected tuple with exactly 2 numbers..."
    xmin, xmax = xmin
  hname2d = hist2d.GetName()
  ixmin, ixmax = getbinrange(hist2d,xmin,xmax=xmax,axis=axis,verb=verb)
  if name is None:
    if xmax is None:
      name = f"{hname2d}_{xmin}".replace('.','p')
    else:
      name = f"{hname2d}_{xmin}-{xmax}".replace('.','p')
  if axis=='x': # integrate bins [ixmin,ixmax] to convert TH2 -> TH1
    hist1d = hist2d.ProjectionY(name,ixmin,ixmax)
  else:
    hist1d = hist2d.ProjectionX(name,ixmin,ixmax)
  return hist1d
  

def getbinrange(hist,xmin,xmax=None,axis='x',verb=0):
  """Get range of bin indices (ixmin, ixmax) that contain x=[xmin,xmax].
  If xmax==None, it will be assumed only one (integer) bin is needed.
  If xmax=='max', xmax will be set to the x axis maximum.
  If xmax=='inf', xmax will be set to the maximum including overflow.
  """
  xaxis = hist.GetXaxis() if axis=='x' else hist.GetYaxis()
  ixmin = xaxis.FindBin(xmin)
  ilast = xaxis.GetNbins() # maximum index
  sxmax = str(xmax)
  xmin_ = xaxis.GetBinLowEdge(ixmin)
  if xmin!=xmin_:
    warning(f"{xmin=} does not align with {xmin_=}, the lower bin edge of bin {ixmin} of {hist!r}",
            pre='getbinrange: ')
  if xmax is None: # assume integer
    ixmax = ixmin
  elif xmax=='max':
    ixmax = ilast # maximum, excluding overflow
  elif xmax=='inf':
    ixmax = ilast+1 # maximum, including overflow
  elif not isinstance(xmax,(int,float)):
    raise IOError()
  else: # assume (finite) float
    sxmax = f"{xmax:.4g}" # float notation
    ixmax = xaxis.FindBin(xmax)
    if xmax==xaxis.GetBinLowEdge(ixmax):
      ixmax -= 1
    else: # xmax must be an upper edge
      warning(f"xmax={xmax} does not align with upper bin edge of bin {ixmax} of {hist!r}",
              pre='getbinrange: ')
  if verb>=2:
    xmin_ = xaxis.GetBinLowEdge(ixmin)
    if ixmax<=ilast:
      of = "" # overflow text
      xmax_ = xaxis.GetBinUpEdge(ixmax)
    else:
      of = " + overflow" # overflow text
      xmax_ = xaxis.GetBinUpEdge(ilast)
    x = 'y' if axis=='x' else 'x'
    print(f">>> getbinrange: ({x}min,{x}max)=({xmin:.4g},{sxmax}) "
          f"=> (i{x}min,i{x}max)=({ixmin},{ixmax}) "
          f" with edges [{xmin_},{xmax_}]{of} of {hist.GetName()!r}")
  return (ixmin,ixmax)
  

def getxlimits(axis):
  """
  Compute axis limits [xmin,xmax].
  In case of a histogram, it has to align with bin edges.
  """
  if hasattr(axis,'GetXaxis'):
    axis = axis.GetXaxis()
  if axis.TestBit(axis.kAxisRange): # if SetRangeUser was used...
    xmin = axis.GetBinLowEdge(axis.GetFirst())
    xmax = axis.GetBinUpEdge(axis.GetLast())
  else:
    xmin = axis.GetXmin()
    xmax = axis.GetXmax()
  return xmin, xmax
  

def getylimits(hists,**kwargs):
  """
  Compute y axis limits from list of histograms.
  Add some margin on top via ymarg
  """
  logy      = kwargs.pop('logy',      False )
  ymin      = kwargs.pop('ymin',      None  )
  ymax      = kwargs.pop('ymax',      None  )
  xmin      = kwargs.pop('xmin',      None  )
  xmax      = kwargs.pop('xmax',      None  )
  ymarg     = kwargs.pop('ymarg',     1.25  ) # top margin
  ybmarg    = kwargs.pop('ybmarg',    None  ) # bottom margin
  logyrange = kwargs.pop('logyrange', 5e5   )
  if isinstance(hists,(list,tuple)):
    hists = [ h for h in hists  if h]
  else:
    hists = [ hists ]
  if xmin is not None or xmax is not None: # restrict range to compute y limits
    if xmin is None: xmin = hists[0].GetXaxis().GetXmin()
    if xmax is None: xmax = hists[0].GetXaxis().GetXmax()
    for hist in hists:
      hist.GetXaxis().SetRangeUser(xmin,xmax)
  #hmax = max(1,max(h.GetMaximum() for h in hists))
  hmax = max(h.GetMaximum() for h in hists)
  if ymin is None:
    if ybmarg:
      hmin = min(h.GetMinimum(0) for h in hists)
      if logy:
        span = ybmarg*abs(log10(hmax/hmin))
        ymin = hmin*10.**(-span)
      else:
        ymin = ybmarg*hmax
    if logy and (ymin is None or ymin<=0):
      ymin = (ymax or hmax or logyrange)/logyrange
    elif ymin is None:
      ymin = 0
  if ymax is None:
    if logy:
      #print(ymarg,hmax,ymax,ymin,ymax/ymin)
      span = ymarg*abs(log10(hmax/ymin))
      ymax = ymin*10.**(span)
    else:
      ymax = ymarg*hmax
  #print(ymin, ymax, ymax, ymarg, hmax)
  return ymin, ymax
  

def divideByBinWidths(hist):
  """Loop over bins and divide by bin width."""
  for i in range(1,hist.GetXaxis().GetNbins()+1): # exclude under-/overflow
    y = hist.GetBinContent(i)
    w = hist.GetXaxis().GetBinWidth(i)
    if y>0 and w>0:
      hist.SetBinContent(i,y/w)
  return hist
  

def mergebins(hold, threshold=0.002, name=None, title=None, axis='x',
              nmaxzeros=3, tolerance=1e-10, verb=0):
  """
  Merge consecutive bins if below some threshold.
  """
  axis1, axis2 = getaxes(hold,axis=axis)
  nbins1    = axis1.GetNbins() # old
  nbins2    = axis2.GetNbins() # old
  ztot      = hold.Integral() # (1,nxbins,1,nybins)
  zthres    = threshold*ztot # merge threshold
  if verb>=1:
    min1, max1 = axis1.GetXmax(), axis1.GetXmin()
    print(f">>> mergebins: Merging {axis} bins ({nbins1},{min1},{max1}) with {threshold=}")
  
  # HELP FUNCTION
  if not isinstance(hold,TH2): # assume TH1
    def getzval(ix):
      return hold.GetBinContent(ix)
  elif axis=='x': # assume TH2
    def getzval(ix):
      return hold.Integral(ix,ix,1,nbins2)
  else: # assume TH2
    def getzval(iy):
      return hold.Integral(1,nbins2,iy,iy)
  
  # LOOP OVER BINS
  zsum   = 0 # cumulative sum/intergral of bins to be merged
  nsum   = 1 # number of bins to be merged
  zsums  = [ ] # for debugging
  edges = [axis1.GetBinLowEdge(1)]
  for ix in range(1,nbins1+1):
    zval = getzval(ix)
    
    # merge bin with previous if cumulative sum is below treshold
    merge = (zsum+zval<zthres)
    # merge single bin if it's above the threshold
    if (zsum==0 and nsum==1 and zval>=zthres):
      merge = True
    # do not merge if several previous bins in a row were all empty (assume nonnegative)
    elif zsum==0 and zval!=0 and nsum>=nmaxzeros:
      merge = False
    # do not merge if several next bins in a row are all empty (assume nonnegative)
    elif zsum!=0 and zval==0 and all(getzval(ix+i)==0 for i in range(1,nmaxzeros)):
      merge = False
    if merge: # merge bins with previous
      zsum += zval
      nsum += 1
      continue
    
    # close previous bin by registering new edge, and reset sum
    edge = axis1.GetBinLowEdge(ix)
    if verb>=2:
      print(f">>> mergebins: Merging {nsum:3} bins [{edges[-1]:6.4g},{edge:6.4g}] with {zsum=:5}")
    edges.append(edge)
    zsums.append(zsum)
    zsum = zval
    nsum = 1
  
  # CLOSE FINAL BIN
  edge = axis1.GetBinUpEdge(nbins1)
  if verb>=2:
    print(f">>> mergebins: Merging {nsum:3} bins [{edges[-1]:6.4g},{edge:6.4g}] with zsum={zsum:5} (final bin!)")
  zsums.append(zsum)
  edges.append(edge) # close final bin edge
  if verb>=2:
    edgestr = [f"{x:.5g}" for x in edges]
    print(f">>> mergebins: Old bins {nbins1} -> {len(edges)-1} with edges {edgestr}")
  elif verb>=1:
    print(f">>> mergebins: Old bins {nbins1} -> {len(edges)-1}")
  zsum = sum(zsums)
  zden = max(ztot,zsum)
  if zden!=0:
    diff = abs(zsum-ztot)/zden # error (should be 0)
    assert diff<tolerance, (f"mergebins: {zsum=} != {ztot=} for {diff=} and {tolerance=} "
                            "=> There may be hole in the algorithms' logic")
  
  # CREATE NEW HISTOGRAM
  hnew = rebinhist(hold,edges,name=name,title=title,axis=axis)
  return hnew
  

def getaxes(hist,axis='x'):
  """Return axes. Swap if axis!='x'."""
  axis1 = hist.GetXaxis()
  axis2 = hist.GetYaxis()
  return (axis1, axis2) if axis=='x' else (axis2, axis1)
  

def rebinhist(hold, edges1, name=None, title=None, axis='x', verb=0):
  """
  Rebin a 2D histogram along one axis ('x' or 'y')
  using given edges aligned with original binning.
  """
  if verb>=1:
    print(f">>> rebinhist: Rebinning {hold.GetName()!r} with {edges1}")
  if name is None:
    name = hold.GetName()+"_rebin"+axis
  if title is None:
    title = hold.GetTitle()
  
  # 1D HISTOGRAM (use standard TH1::Rebin routine)
  if not isinstance(hold,TH2):
    xbins = np.array(edges1,dtype='float64')
    hnew  = hold.Rebin(len(edges1)-1,name,xbins)
    hnew.SetTitle(title)
    return hnew
  
  # SET AXES (1=to be rebinned, 2=not rebinned)
  axis1, axis2 = getaxes(hold,axis=axis)
  nbins1 = axis1.GetNbins() # old
  nbins2 = axis2.GetNbins() # old
  nbins1_new = len(edges1)-1
  
  # NEW BINNING
  bins1 = (nbins1_new,np.array(edges1,dtype='float64'))
  if axis2.IsVariableBinSize():
    edges2 = [axis2.GetBinLowEdge(i) for i in range(1,nbins2+2)]
    edges2 = np.array(edges2,dtype='float64')
    bins2 = (nbins2, np.array(edges2,dtype='float64'))
  else:
    bins2 = (nbins2, axis2.GetXmin(), axis2.GetXmax())
  
  # CREATE NEW HISTOGRAM
  if axis=='x':
    hnew = TH2D(name, title, *bins1, *bins2)
  else:
    hnew = TH2D(name, title, *bins2, *bins1)
  
  # FILL NEW HISTOGRAM
  minbw = 0.1*min(axis1.GetBinWidth(i) for i in range(1,nbins1+1))
  for i1_new in range(1,nbins1_new+1):
    i1_min = axis1.FindBin(edges1[i1_new-1]+minbw)
    i1_max = axis1.FindBin(edges1[i1_new]  -minbw)
    for i2 in range(1,nbins2+1):
      if axis=='x': # integrate x axis
        zval = hold.Integral(i1_min,i1_max,i2,i2)
        hnew.SetBinContent(i1_new,i2,zval)
      else: # integrate y axis
        zval = hold.Integral(i2,i2,i1_min,i1_max)
        hnew.SetBinContent(i2,i1_new,zval)
  hnew.SetEntries(hold.GetEntries())
  
  return hnew
  
