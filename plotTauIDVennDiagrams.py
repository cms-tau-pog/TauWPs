#! /usr/bin/env python3
# Author: Izaak Neutelings (August 2025)
# Description: Study overlap and complementarity/recovery
# Instructions:
#   python3 plotTauIDVennDiagrams.py output/tauAlgoEffs_2024.root -v1
# References:
#   https://repository.cern/records/ysxvz-waj51
import os, sys
from datetime import datetime
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True # to avoid conflict with argparse
from ROOT import gROOT, gStyle, TCanvas, TFile, TH2F, TLine, TLatex,\
                 kBlack, kGray, kRed
sys.path.append("python")
from utils import warning, loadconfig
from histutils import getHistFromFile
from plotutils import drawtext, setera, lcolors
from vennutils import TH2Venn, VennBins, drawvenncircles, getvbincenters #, polar2xy
import cmsstyle as CMS
CMS.SetCMSPalette()
gROOT.SetBatch(True) # do not open GUI windows
CONFIG = loadconfig()

# HISTOGRAM BINS
ibins_dict = { # (binx,biny,binz) -> VennBins enums
  (2,2,2): VennBins.AnBnC, # intersection of A, B, and C
  (2,2,1): VennBins.AnB,   # intersection of A and B
  (2,1,2): VennBins.AnC,   # intersection of A and C
  (1,2,2): VennBins.BnC,   # intersection of B and C
  (2,1,1): VennBins.A,     # A without B or C
  (1,2,1): VennBins.B,     # B without A or C
  (1,1,2): VennBins.C,     # C without A or B
}


def getVennHist3DFromFile(infile,gm,wp=None,norm=True,verb=0):
  """Return TH3 Venn histogram from file."""
  postfix = '' if wp is None else f"_{wp}"
  hname   = f"venn/ntau_venn{postfix}_{gm}"
  hist    = getHistFromFile(hname,infile)
  algos   = [a.GetTitle() for a in [hist.GetXaxis(),hist.GetYaxis(),hist.GetZaxis()]]
  ntot    = hist.Integral()-hist.GetBinContent(1,1,1)
  if hist.GetBinContent(1,1,1)!=0: # sanity check
    warning(f"(1,1,1) bin has yval={yval:.4g} (expected 0)",pre='getVennHist3DFromFile')
  if norm and ntot>0: # normalize
    hist.Scale(100./ntot) # in percentage
    ntot = 1.
  return hist, algos
  

def getVennDataFromFile(infile,gm,wp,norm=True,verb=0):
  """Get data from TH3 object from file.
  Returns dictionary of count per partition."""
  algo_dict = CONFIG['text']
  h_venn, algos = getVennHist3DFromFile(infile,gm,wp,norm=norm,verb=verb)
  part_dict = { } # count per partition
  for ibins in ibins_dict:
    yval = h_venn.GetBinContent(*ibins)
    ckey = [ ] # collection key in part_dict
    for i, algo in zip(ibins,algos):
      if i==2: # algo==True
        algo = algo_dict.get(algo,algo)
        ckey.append(algo)
    if ckey: # e.g. ('DeepTau','PNetCHS') or ('DeepTau','PNetCHS','UParTPUPPI')
      part_dict[tuple(ckey)] = yval
  return part_dict, algos
  

def getVennHist(infile,gm,wp=None,norm=True,verb=0,**kwargs):
  """Convert TH3 object to TH2Venn (=TH2Poly) object."""
  h_3d, algos = getVennHist3DFromFile(infile,gm,wp,norm=norm,verb=verb)
  hname = h_3d.GetName()
  kwargs.setdefault('name',hname)
  h_venn = TH2Venn(**kwargs)
  for ibins, vbin in ibins_dict.items():
    yval = h_3d.GetBinContent(*ibins)
    h_venn.SetBinContent(vbin,yval)
  return h_venn, algos
  

def makeVennDiagram(infile,genmatches,wp,tag='',plotdir='plots/venn',verb=0):
  """Create Venn diagram histogram, compute (relative) fractions & plot."""
  algo_dict = CONFIG['text']['algo']
  gm_dict   = CONFIG['text']['gm']
  
  # GET DATA
  algos = None
  for gm in genmatches:
    hist, algos_ = getVennHist(infile,gm,wp,norm=True,verb=verb,marg=0.35)
    if algos is None:
      algos = algos_
    elif algos!=algos_:
      warning(f"Different algorithms extracted from axes, algos={algos}!={algos_}",pre="plotVennTables")
    
    # SUM
    sublabs = [ ]
    relfracs = [ ]
    for set in 'ABC':
      ntot = 0
      for vbin in VennBins:
        if set in vbin.name:
          ntot += hist.GetBinContent(vbin)
      sublabs.append(f"({ntot:.1f}%)")
      relfracs.append(ntot)
    
    # PLOT
    zmax   = max(hist.GetMaximum(),74 if gm=='real' else 40)
    labels = [algo_dict.get(a,a) for a in algos]
    #print(labels,algos)
    ztitle = "Partition fraction [%]"
    text   = [gm_dict[gm].replace('jet','q/g').replace('real','Real').replace('gen','Gen').replace('muon','Muon')]
    if gm!='real':
      text[0] += " mis-ID"
    etag = f"_{gm}"
    if wp is None:
      text.append("No WPs")
    else:
      etag += f"_{wp}"
      text.append(f"{wp} VSjet") #(VVVL VSe, VL VSmu)
    plotVennDiagram(hist,labels=labels,sublabels=sublabs,drawopt='COLZ TEXT',
                    zmax=zmax,ztitle=ztitle,text=text,outdir=plotdir,tag="_fraction"+etag+tag)
    plotVennDiagram(hist,labels=labels,sublabels=sublabs,drawopt='COLZ',relfracs=relfracs,
                    zmax=zmax,ztitle=ztitle,text=text,outdir=plotdir,tag="_relfracs"+etag+tag)
  

def plotVennDiagram(hist, drawopt='COL TEXT L', tag="", outdir='plots/venn', **kwargs):
  """Plot Venn diagram."""
  logz     = kwargs.pop('logz',      False           )
  ztitle   = kwargs.pop('ztitle',    None            )
  text     = kwargs.pop('text',      None            )
  format   = kwargs.pop('format',    ".1f"           )
  relfracs = kwargs.pop('relfracs',  None            ) # indicate relative fractions
  labels   = kwargs.pop('labels',    None            )
  sublabs  = kwargs.pop('sublabels', None            )
  zmax     = kwargs.pop('zmax',      None            )
  exts     = kwargs.pop('ext',       ['.png','.pdf'] )
  
  # CANVAS
  cname = 'canvas'
  lmarg, rmarg = 0, 0
  bmarg, tmarg = 0.02, 0.02
  width, height = 1100, 1000
  if 'COLZ' in drawopt:
    width = 1300 # extra space for color scale bar
    rmarg = 0.18
  canvas = TCanvas(cname, cname, width, height)
  
  # DRAW
  tsize = 0.050
  gStyle.SetPaintTextFormat(format)
  hist.SetMarkerSize(2.5)
  hist.SetMarkerColor(kRed)
  canvas.SetMargin(lmarg,rmarg,bmarg,tmarg) # LRBT
  canvas.SetFrameLineColor(0)
  canvas.SetFrameFillStyle(0)
  canvas.SetTicks(False,False)
  canvas.SetLogz(logz)
  if ztitle:
    hist.GetZaxis().SetLabelOffset(-0.003 if logz else 0.006)
    hist.GetZaxis().SetTitleSize(1.05*tsize)
    hist.GetZaxis().SetTitleOffset(0.99)
    hist.GetZaxis().SetTitle(ztitle)
  if zmax is not None:
    hist.SetMaximum(zmax)
  for axis in [hist.GetXaxis(),hist.GetYaxis()]:
    axis.SetLabelSize(0)
    axis.SetTickLength(0)
    axis.SetAxisColor(0,0) # transparent
  canvas.Update()
  
  # DRAW
  hist.Draw(drawopt)
  circles = drawvenncircles(colors=lcolors)
  
  # INDICATE RELATIVE FRACTIONS
  latex = TLatex()
  latex.SetTextSize(tsize)
  latex.SetTextFont(42)
  if relfracs:
    bincenters = getvbincenters()
    for vbin, (x,y) in bincenters.items():
      bname = vbin.name.replace('n','')
      nsets = len(bname)
      for iset, set in enumerate(bname):
        ialgo = 'ABC'.index(set) # algo index
        ntot  = relfracs[ialgo]
        zval  = hist.GetBinContent(vbin)
        label = f"{100.0*zval/ntot:.1f}" if ntot>0 else ''
        if nsets==1: # A, B, or C
          y_ = y
          align = 22
        elif nsets==2: # AnB, AnC, or BnC
          y_ = y+0.05*(0.5-iset)
          align = 21+2*iset
        else: # AnBnC
          y_ = y+0.105*(0.68-iset)
          align = 21+iset
        latex.SetTextAlign(align)
        latex.SetTextColor(kGray+1)
        latex.SetTextFont(62) # bold (as shadow)
        latex.DrawLatex(x,y_,label)
        latex.SetTextColor(lcolors[ialgo])
        latex.SetTextFont(42) # normal (undo bold)
        latex.DrawLatex(x,y_,label)
  
  # SET LABELS of CIRCLES (=sets)
  latex.SetTextSize(tsize)
  if labels is not None:
    aligns = [13,33,31] #aligns = [31,11,13]
    xys = [(0.02,0.985-tmarg),(0.98-1.02*rmarg,0.985-tmarg),(0.98-1.02*rmarg,0.08)]
    #polars = ((1.52,105),(1.62,64),(1.7,-70))
    for i, label in enumerate(labels):
      x, y = xys[i] #*polar2xy(*polars[i])
      if i==2 and len(label)<14: # tweak
        x -= 0.004*max(0,15-len(label))
      latex.SetTextFont(62) # bold
      latex.SetTextAlign(aligns[i])
      latex.SetTextColor(lcolors[i])
      latex.DrawLatexNDC(x,y,label)
      if sublabs is not None: # (percentage)
        latex.SetTextFont(42) # normal
        latex.DrawLatexNDC(x,y-1.05*tsize,sublabs[i])
  
  # TEXT (tau type, selection)
  if text:
    nrows = len(text) if isinstance(text,list) else text.count('\n')
    y = 0.01+1.1*(nrows-1)*tsize
    drawtext(text,x=0.01,y=y,tsize=tsize,align=11)
  
  # WRITE
  pname = os.path.join(outdir,f"venn_diagram{tag}")
  for ext in exts:
    canvas.SaveAs(pname+ext)
  

def makeVennTables(infile,genmatches,wps,tag='',outdir='output',plotdir='plots/venn',verb=0):
  """
  Compute Venn diagram of three tau IDs, print & plot table.
  See https://repository.cern/records/ysxvz-waj51
  """
  
  # GET DATA
  algos = None
  part_dict = { } # count per partition
  for wp in wps:
    part_dict[wp] = { }
    for gm in genmatches:
      fracs, algos_ = getVennDataFromFile(infile,gm,wp,norm=True,verb=verb)
      part_dict[wp][gm] = fracs
      if algos is None:
        algos = algos_
      elif algos!=algos_:
        warning(f"Different algorithms extracted from axes, algos={algos}!={algos_}",pre="plotVennTables")
  
  # ROW & COLUMN LABELS
  wp0, gm0 = wps[0], genmatches[0]
  parts    = list(part_dict[wp0][gm0].keys()) # algo partitions
  lab_dict = { p: (f"{p[0]} only" if len(p)==1 else "All three" if len(p)==3 else ' & '.join(p))
                for p in parts } # labels per partition
  wp_dict  = { wp: ("No WPs" if wp is None else (wp+" VSjet")) for wp in wps }
  
  # HEADERS
  gm_real  = 'real'
  rows_frac = [ "Partition fraction" ]
  rows_prty = [ "Purity = eff / mis-ID" ]
  rows_cntm = [ "Contamination = mis-ID / eff" ]
  labwid    = max(len(l) for l in lab_dict.values())
  gms_prty  = [gm for gm in genmatches if gm!=gm_real]
  for rows, gms in [(rows_frac,genmatches),(rows_prty,gms_prty),(rows_cntm,gms_prty)]:
    colwid = 7*len(gms)-2
    wphead = '|'.join('  '+wp_dict[wp].ljust(colwid) for wp in wps)+'|'
    gmhead = (' '.join(f"{gm:>6}" for gm in gms)+" |")*len(wps)+' '*(labwid+2)
    rows.extend([wphead,gmhead])
  
  # CONSTRUCT HISTOGRAMS
  nrows       = len(parts)
  ncols_frac  = len(wps)*len(genmatches)
  ncols_prty  = len(wps)*len(gms_prty) # without real
  ztitle_prty = "Purity = f_{#lower[-0.35]{genuine}} / f_{#lower[-0.25]{mis-ID}}"
  ztitle_cntm = "Contamination = f_{#lower[-0.35]{mis-ID}} / f_{#lower[-0.25]{genuine}}"
  hist_frac   = TH2F('fraction',"Partition fraction [%]",ncols_frac,0,ncols_frac,nrows,0,nrows)
  hist_prty   = TH2F('purity',ztitle_prty,ncols_prty,0,ncols_prty,nrows,0,nrows)
  hist_cntm   = TH2F('contam',ztitle_cntm,ncols_prty,0,ncols_prty,nrows,0,nrows)
  
  # FILL DATA into tables and histograms
  for iy, part in enumerate(parts,1): # loop over partitions
    row_frac, row_prty, row_cntm = "", "", ""
    ix_frac, ix_prty = 0, 0
    for wp in wps:
      eff_real = part_dict[wp][gm_real][part]
      for gm in genmatches:
        ix_frac += 1
        eff = part_dict[wp][gm][part]
        row_frac += f"{eff:6.1f} "
        hist_frac.SetBinContent(ix_frac,iy,eff)
        if gm!=gm_real: # purity = efficiency / mid-ID
          ix_prty += 1
          if eff>0: #and eff_real>0:
            purity = (eff_real/eff) #100.0*
            contam = (eff/eff_real) if eff>0 else 999.
            row_prty += f"{purity:6.1f} "
            row_cntm += f"{contam:6.1f} "
            hist_prty.SetBinContent(ix_prty,iy,purity)
            hist_cntm.SetBinContent(ix_prty,iy,contam)
          else:
            rows_prty += ' '*6
            rows_cntm += ' '*6
      row_frac += '|'
      row_prty += '|'
      row_cntm += '|'
    rows_frac.append(row_frac+' '+lab_dict[part])
    rows_prty.append(row_prty+' '+lab_dict[part])
    rows_cntm.append(row_cntm+' '+lab_dict[part])
  
  # PRINT
  print(">>> ")
  for rows in [rows_frac,rows_prty,rows_cntm]:
    for i, row in enumerate(rows):
      if i in [1,2]:
        row = "\033[4m"+row+"\033[0m"
      print(f">>>   "+row)
    print(">>> ")
  
  # WRITE TO TXT
  outfname = os.path.join(outdir,f"venn_frac{tag}.txt")
  with open(outfname,'w') as outfile:
    outfile.write(f"# Generated by {__file__}.plotVennTables on {datetime.now():%Y-%m-%d %H:%M:%S}\n\n")
    for rows in [rows_frac,rows_prty]:
      outfile.write('\n')
      for i, row in enumerate(rows):
        if i in [1,2]:
          row = row.replace(' ','_')
        outfile.write(row+'\n')
      outfile.write('\n')
  
  # PLOT TABLE
  hist_prty.SetMaximum(25)
  hist_cntm.SetMaximum(25)
  xlabs_frac = [CONFIG['text']['gm_short'][gm] for gm in genmatches*len(wps)]
  xlabs_prty = [CONFIG['text']['gm_short'][gm] for gm in gms_prty*len(wps)]
  ylabs      = [lab_dict[p] for p in parts]
  wplabs     = [wp_dict[wp] for wp in wps]
  plotVennTables(hist_frac,xlabs_frac,ylabs,wplabs,outdir=plotdir)
  plotVennTables(hist_frac,xlabs_frac,ylabs,wplabs,outdir=plotdir,tag="_log",logz=True)
  plotVennTables(hist_prty,xlabs_prty,ylabs,wplabs,outdir=plotdir)
  plotVennTables(hist_prty,xlabs_prty,ylabs,wplabs,outdir=plotdir,tag="_log",logz=True)
  plotVennTables(hist_cntm,xlabs_prty,ylabs,wplabs,outdir=plotdir)
  plotVennTables(hist_cntm,xlabs_prty,ylabs,wplabs,outdir=plotdir,tag="_log",logz=True)
  

def plotVennTables(hist,xlabs,ylabs,wplabs,tag="",outdir='plots/venn',**kwargs):
  cname = 'canvas'
  ztitle = kwargs.pop('ztitle', hist.GetTitle() )
  logz   = kwargs.pop('logz',   False           )
  exts   = kwargs.pop('ext',    ['.png','.pdf'] )
  xmin, xmax = hist.GetXaxis().GetXmin(), int(hist.GetXaxis().GetXmax())
  ymin, ymax = hist.GetXaxis().GetXmin(), int(hist.GetYaxis().GetXmax())
  
  # SET CANVAS
  lmarg, rmarg = 0.36, 0.13
  bmarg, tmarg = 0.12, 0.073
  cw, ch = 1400, 700
  canvas = TCanvas(cname,cname,cw,ch)
  canvas.SetMargin(lmarg,rmarg,bmarg,tmarg) # LRBT
  canvas.SetTicks(1,1)
  canvas.SetLogz(logz)
  gStyle.SetPaintTextFormat(".1f") # for marker text
  
  # SET AXES
  tsize = 0.050
  lsize = 1.16*tsize
  objects = [ ] # to prevent python's garbage collection
  hist.GetXaxis().SetLabelOffset(0.0044)
  hist.GetYaxis().SetLabelOffset(0.0025)
  hist.GetZaxis().SetLabelOffset(-0.004 if logz else 0.004)
  hist.GetXaxis().SetLabelSize(1.2*lsize)
  hist.GetYaxis().SetLabelSize(lsize)
  hist.GetYaxis().SetTitleSize(tsize)
  hist.GetZaxis().SetTitleSize(tsize)
  hist.GetXaxis().SetTitleOffset(1.05)
  hist.GetZaxis().SetTitleOffset(0.80 if logz else 0.76)
  hist.GetXaxis().SetNdivisions(6)
  hist.GetYaxis().SetNdivisions(6)
  hist.GetXaxis().SetTickLength(0.06)
  hist.GetZaxis().SetTitle(ztitle)
  hist.SetMarkerSize(1.8 if xmax==12 else 1.9)
  hist.SetMarkerColor(kRed)
  
  # DRAW
  #print(frame.SetBins(xmax,0,xmax,ymax,0,ymax))
  #print(xmax,ymax,frame.GetXaxis().GetNbins(),frame.GetYaxis().GetNbins())
  hist.Draw('COLZ TEXT')
  
  # COLOR BAR / PALETTE (reduce width)
  canvas.Update()
  palette = hist.GetListOfFunctions().FindObject('palette')
  if palette:
    x1, x2 = palette.GetX1NDC(), palette.GetX2NDC()
    pw = abs(x2-x1) # original width
    palette.SetX1NDC(x1+0.08*pw) # left edge
    palette.SetX2NDC(x1+0.90*pw) # right edge
  
  # Y AXIS LABELS
  latex = TLatex()
  latex.SetTextSize(0.75*lsize)
  latex.SetTextAlign(32)
  latex.SetTextFont(42)
  latex.SetTextColor(kBlack)
  objects.append(latex)
  for i, label in enumerate(ylabs,1):
    ###hist.GetYaxis().SetBinLabel(i,label)
    ###hist.GetYaxis().ChangeLabel(i,0,-1,32,-1,-1,label)
    y = bmarg+(i-0.5)/len(ylabs)*(1-bmarg-tmarg)
    if '&' in label:
      #parts = label.split('&')
      latex.DrawLatexNDC(lmarg-0.005,y+0.22*lsize,label)
      latex.DrawLatexNDC(lmarg-0.005,y-0.48*lsize,"only")
    else:
      latex.DrawLatexNDC(lmarg-0.005,y,label)
  
  # X AXIS LABELS
  for i, label in enumerate(xlabs,1):
    if label=='e': # tweak
      label = f"#lower[0.09]{{{label}}}"
    elif label=='#mu': # tweak
      label = f"#lower[0.10]{{{label}}}"
    elif '#tau' in label: # tweak
      label = f"#lower[0.09]{{{label}}}"
    hist.GetXaxis().SetBinLabel(i,label)
    #hist.GetXaxis().ChangeLabel(i,-10,0.75*lsize)#,-1,-1,-1,label)
  hist.GetYaxis().SetLabelSize(0)
  
  # X AXIS SUBDIVISIONS (WP)
  ncols = len(xlabs)/len(wplabs)
  ndc2y = bmarg*(ymax-ymin)/(1-bmarg-tmarg)
  latex.SetTextSize(0.97*tsize)
  latex.SetTextAlign(21)
  for i, wplab in enumerate(wplabs):
    #latex.DrawLatex(ncols*(i+0.52),ymin-0.84*ndc2y,wplab)
    latex.DrawLatex(ncols*(i+0.5),ymin-0.84*ndc2y,wplab)
    if i>0:
      x = ncols*i
      line = TLine(x,ymin-0.9*ndc2y,x,ymax)
      line.SetLineWidth(1)
      line.SetLineColor(kBlack)
      line.Draw('SAME')
      objects.append(line)
  for y in [1,4]:
    line = TLine(xmin,y,xmax,y)
    line.SetLineWidth(1)
    line.SetLineColor(kBlack)
    line.Draw('SAME')
    objects.append(line)
  
  # WRITE
  CMS.CMS_lumi(canvas,0)
  #CMS.UpdatePad(canvas)
  canvas.Modified()
  canvas.Update()
  canvas.RedrawAxis()
  pname = os.path.join(outdir,f"venn_table_{hist.GetName()}{tag}")
  for ext in exts:
    canvas.SaveAs(pname+ext)
  

def main(args):
  
  # GLOBAL SETTINGS from config.yaml
  genmatches    = CONFIG['genmatches']
  vennWPs       = CONFIG['venn']['wps']
  
  # USER SETTINGS
  verbosity = args.verbosity
  infnames  = args.inputFile
  outdir    = args.outdir
  plotdir   = args.plotdir
  tag       = args.tag
  era       = args.era
  exts      = ['.png','.pdf']
  for dir in [plotdir,outdir]:
    if dir:
      os.makedirs(dir,exist_ok=True)
  
  # INPUT
  infname = infnames[0].replace('$TAG',tag)
  infile = TFile(infname,'READ')
  
  # ERA/LUMI TEXT
  era, etag = setera(era,infname)
  
  # VENN TABLE
  wps = [None]+vennWPs
  setera(None,extratext="Simulation Preliminary") # for DPS note
  makeVennTables(infile,genmatches,wps,tag=etag+tag,outdir=outdir,verb=verbosity)
  setera(era) # reset
  
  # VENN DIAGRAM
  for wp in wps:
    makeVennDiagram(infile,genmatches,wp,tag=etag+tag,plotdir=plotdir,verb=verbosity)
  
  # CLOSE
  infile.Close()
  

if __name__ == '__main__':
  from argparse import ArgumentParser
  description = """This script compares tau ID algorithms in tables and Venn diagrams."""
  parser  = ArgumentParser(description=description,epilog="Good luck!")
  parser.add_argument('inputFile',           nargs='+',
                      metavar='INPUT_FILE',  help="input ROOT file from makeTauHists.py" )
  parser.add_argument('-o', '--plotdir',     default='plots/venn',
                                             help="output directory for plots, default=%(default)r" )
  parser.add_argument('-O', '--outdir',      default='output',
                                             help="output directory for data (txt & yaml), default=%(default)r" )
  parser.add_argument('-e', '--era',         help="era for lumi text in plots (e.g. '2024')" )
  parser.add_argument('-t', '--tag',         default="",
                                             help="tag for default output ROOT file name" )
  parser.add_argument('-v', '--verbose',     dest='verbosity', type=int, nargs='?', const=1, default=0,
                                             help="set level of verbosity, default=%(default)s" )
  args = parser.parse_args()
  main(args)
  print(">>> Done!")
  
