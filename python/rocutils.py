# Author: Izaak Neutelings (November 2026)
import os
import yaml
from math import atanh
from utils import warning, loadconfig
from ROOT import TFile
CONFIG = loadconfig() # settings like cuts, binning, labels, ...


def getwpcuts(algo,outdir=CONFIG['wp']['outdir'],tag="",discriminants=None,
              wps=None,defcut=CONFIG['wp']['cut'],logit=False):
  """Help function to read YAML file with WPs cuts."""
  fname = os.path.join(outdir,CONFIG['wp']['file'].replace('$ALGO',algo).replace('$TAG',tag))
  wp2cut_dict = { }
  if not os.path.isfile(fname):
    warning(f"Could not find WP cuts file {fname}! Returning empty dict...",pre='getwpcuts: ')
  else:
    with open(fname,'r') as infile:
      wp2cut_dict = yaml.safe_load(infile)
  ###if wps:
  ###  wps = tuple('Reco' if wp=='None' else wp for wp in wps)
  if not wp2cut_dict:
    warning(f"Could not load WP cut for {algo}! Use plotTauROCCurves.py to create a YAML. Using default cut={defcut}...",pre="getwpcuts: ")
    if discriminants and wps:
      wp2cut_dict = { d: { wp: { 'cut': defcut } for wp in wps } for d in discriminants }
  if discriminants and wps: # sanity checks
    for disc in discriminants:
      if disc not in wp2cut_dict:
        warning(f"Did not find {disc!r} in wp2cut dictionary! Using default cut={defcut}...",pre="getwpcuts: ")
        wp2cut_dict[disc] = { wp: { 'cut': defcut } for wp in wps } # defaults
      for wp in wps:
        if disc=='VSmu' and any(wp.startswith(s) for s in ['VVL','VT']):
          continue # impossibly WPs for VSmu
        if wp not in wp2cut_dict[disc]:
          warning(f"Did not find WP {wp!r} for {disc!r} in wp2cut dictionary! Using default cut={defcut}...",pre="getwpcuts: ")
  if logit: # transform cut on (linear) score to cut on logit(score)
    for disc in wp2cut_dict:
      for wp, cutinfo in wp2cut_dict[disc].items():
        cut = cutinfo['cut']
        cutinfo['cut'] = 2*atanh(2*cut-1) if cut>0 else -99
  return wp2cut_dict
  

def getROCFiles(infnames,tag='',genmatches=CONFIG['genmatches']):
  """
  Parse input files from command line. Allow keys like
    'muon=tauAlgoEffs_DY.root elec=tauAlgoEffs_DY.root jet=tauAlgoEffs_TT.root'
  """
  allowedkeys = CONFIG['rockeys']
  
  # MAIN INPUT
  infname0 = infnames[0].replace('$TAG',tag)
  infile0  = TFile.Open(infname0,'READ')
  
  # GET ALTERNATIVE FILES
  infiles = { k: infile0 for k in genmatches } # split tau candidates by genmatch
  for infname in infnames:
    nkeys = infname.count('=')
    if nkeys==1:
      key, infname = infname.split('=')
      assert key in allowedkeys, f"Unrecognized {key=}! Allowed: {allowedkeys}"
    elif nkeys==0:
      key = 'real'
    else:
      raise IOError(f"{infname} does contains {nkeys} '='; expected exaclty 0 or 1!")
    infile_ = TFile.Open(infname,'READ')
    if key=='fake': # use for elec, muon, jet
      for gm in genmatches:
        if gm!='real':
          infiles[gm] = infile_
    else:
      infiles[key] = infile_
  
  return infile0, infiles
  

def getProcFromFileName(fname):
  """Guess the physics process from the file name."""
  if isinstance(fname,TFile):
    fname = os.path.basename(fname.GetName())
  proc_dict = CONFIG.get('text', {}).get('process', {})
  for prockey in proc_dict:
    if prockey in fname:
      return proc_dict[prockey]
  return 'unkown'
  

def addCommonROCArgs(parser, outdir='plots', methods=None): #addwptag=False
  #from argparse import Action
  allowedkeys = CONFIG['rockeys']
  parser.add_argument('inputFile',          nargs='+', action='store', #action=KeyFileAction,
                      metavar='INPUT_FILE', help="input ROOT file(s) from makeTauHists.py; add other files with a key "
                                                 "to get the fakes from other sources, like 'muon=tauAlgoEffs_DY.root elec=tauAlgoEffs_DY.root jet=tauAlgoEffs_TT.root'"
                                                 f", allowed keys: {', '.join(allowedkeys)}" )
  parser.add_argument('-o', '--outdir',     default=outdir,
                                            help="output directory for plots, default=%(default)r" )
  parser.add_argument('-e', '--era',        help="era for lumi text in plots, e.g. '2024'" )
  parser.add_argument('-t', '--tag',        default="",
                                            help="tag for default output ROOT file name" )
  ###if addwptag:
  ### parser.add_argument('-w', '--wptag',     default=CONFIG['wp']['tag'], #"_2024_hps3",
  ###                                          help="tag for WP YAML, e.g. $output/wps_$algo$wptag.yaml" )
  if methods:
    parser.add_argument('-m', '--method',     dest='methods', nargs='+', choices=methods, default='all',
                                              help="methods to preform; default=%(default)r" )
  parser.add_argument('-v', '--verbose',    dest='verbosity', type=int, nargs='?', const=1, default=0,
                                            help="set level of verbosity, default=%(default)s" )
  
