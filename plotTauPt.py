#! /usr/bin/env python3
# Author: Izaak Neutelings (September 2025)
# Instructions:
#   python3 plotTauPt.py output/tauAlgoEffs_2024_DY.root -v1
import os, sys
sys.path.append("python")
from utils import loadconfig
from histutils import getHistFromFile
from plotutils import setera, plotHistComparison, kDashDoubleDotted
from ROOT import TFile, kSolid, kDashed #, kDotted
CONFIG = loadconfig() # settings like cuts, binning, labels, ...


def main(args):
  
  # GLOBAL SETTINGS from config.yaml
  tauAlgos  = CONFIG['tauAlgos']
  algo_dict = CONFIG['text']['algo']
  text_dict = CONFIG['text']['selection']
  
  # USER SETTINGS
  verbosity = args.verbosity
  plotdir   = args.plotdir
  tag       = args.tag
  era       = args.era
  infname   = args.inputFile.replace('$TAG',tag)
  exts      = ['.png','.pdf']
  for dir in [plotdir]:
    if dir:
      os.makedirs(dir,exist_ok=True)
  
  # ERA/LUMI TEXT
  infile = TFile.Open(infname,'READ')
  era, etag = setera(era,infname)
  
  # PT RESPONSE PLOTS
  xvar = 'ptres'
  prefixes = [ 'qlty', 'hps_hps3', 'rec_hps3', 'hps_hps3_VVVLoose', 'rec_hps3_VVVLoose' ]
  for prefix in prefixes:
    if 'hps3_VVVLoose' in prefix:
      setera(None,extratext="Simulation Preliminary") # for DPS note
    tkey     = prefix.replace('_VVVLoose','') #'qlty_hps' if 'DeepTau' in algo else 'qlty_jet'
    texts    = text_dict.get(tkey,text_dict['all'])
    ###texts[0] = cut_dict['all_gen'] # add generator superscript
    if 'VVVLoose' in prefix:
      texts.append("VVVLoose VSjet, VVVLoose VSe, Loose VSmu")
    hists    = [ ]
    
    # EXTRA LEGEND
    legobjs1 = [ ]
    legobjs2 = [ ]
    lwidths  = [ ]
    lstyles  = [ ]
    for algo in tauAlgos:
      entry = algo_dict.get(algo,algo)
      if 'DeepTau' in algo:
        if ('jpt' in xvar or 'rec' in prefix):
          continue
        if 'Loose' not in prefix:
          entry = "HPS"
      hname = f"{xvar}{algo}_{prefix}"
      hist  = getHistFromFile(f"pt/{algo}/{hname}",infile)
      hist.SetTitle(entry)
      if 'DeepTau' in algo:
        legobjs1.append(hist) # first legend
        lwidths.append(3)
        lstyles.append(kSolid)
      else: # PNet or UParT
        legobjs2.append(hist) # second legend
        lwidths.append(2)
        lstyles.append(kSolid if 'CHS' in algo else
                       kDashed if 'PNet' in algo else kDashDoubleDotted)
        #lstyle = kDashed if 'PNetPUPPI' in algo else kSolid
      hists.append(hist)
    headers = ["Regressed to p_{#lower[-0.2]{T}}^{#lower[0.35]{gen}}:"]
    if legobjs1: # two legends
      icol0 = 0
      yleg  = 0.84-0.05*len(texts)
      headers.insert(0,"Reconstructed:")
      legobjs = [legobjs1,legobjs2]
    else: # one legend
      icol0 = 1
      yleg  = 0.65
      legobjs = [legobjs2]
    
    ymarg  = 1.20+0.055*len(texts)
    fname  = f"{xvar}{etag}_{prefix}{tag}"
    xtitle = "Response p_{#lower[-0.2]{T}}^{#lower[0.2]{reco}} / p_{#lower[-0.2]{T}}^{#lower[0.25]{gen}}"
    plotHistComparison(hists,fname,outdir=plotdir,ext=exts,ratio=False,grid=True,
                       xtitle=xtitle,icol0=icol0,lwidth=lwidths,lstyle=lstyles,
                       header=headers,legobjs=legobjs,ncols=1,xleg=0.175,yleg=yleg,
                       xmin=0.33,xmax=1.70,ymarg=ymarg,text=texts,verb=verbosity)
    setera(era) # reset
  
  #### PT PLOTS
  ###gm_real  = 'real'
  ###prefixes = [ 'qlty' ]
  ###for xvar, prefix in prefixes:
  ###  tkey = prefix #'qlty_hps' if 'DeepTau' in algo else 'qlty_jet'
  ###  text = text_dict.get(tkey)
  ###  hists  = [ ]
  ###  for algo in tauAlgos:
  ###    if 'DeepTau' in algo and 'jpt' in xvar:
  ###      continue
  ###    hname = f"{xvar}{algo}_{prefix}_{gm_real}"
  ###    hist  = getHistFromFile(f"{algo}/{hname}",infile)
  ###    hists.append(hist)
  ###  fname = f"{xvar}_{gm_real}{etag}_{prefix}{tag}"
  ###  plotHistComparison(hists,fname,outdir=plotdir,
  ###                     ext=exts,text=text,verb=verbosity)
  
  infile.Close()
  

if __name__ == '__main__':
  from argparse import ArgumentParser
  description = """This script plots pT and pT response."""
  parser  = ArgumentParser(description=description,epilog="Good luck!")
  parser.add_argument('inputFile',          #nargs='+',
                      metavar='INPUT_FILE', help="input ROOT file from makeTauHists.py" )
  parser.add_argument('-o', '--plotdir',    default='plots/pt',
                                            help="output directory for plots, default=%(default)r" )
  parser.add_argument('-e', '--era',        help="era for lumi text in plots (e.g. '2024')" )
  parser.add_argument('-t', '--tag',        default="",
                                            help="tag for default output ROOT file name" )
  parser.add_argument('-v', '--verbose',    dest='verbosity', type=int, nargs='?', const=1, default=0,
                                            help="set level of verbosity, default=%(default)s" )
  args = parser.parse_args()
  main(args)
  print(">>> Done!")
  
