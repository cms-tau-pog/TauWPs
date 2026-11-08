#! /usr/bin/env python3
# Author: Izaak Neutelings (August 2025)
# Instructions to install for ROOT v6.32:
#   . /cvmfs/sft.cern.ch/lcg/views/LCG_106c/x86_64-el9-gcc13-opt/setup.sh
#   pip3 install --user dask_lxplus
# Instructions to setup the environment with each clean session:
#   kinit && aklog # Kerberos authentication for HTCondor
#   voms-proxy-init -voms cms -valid 200:0 # GRID authentication
#   export PYTHONPATH=~/.local/lib/python3.11/site-packages:$PYTHONPATH
# Instructions to run locally:
#   python3 makeTauHists.py /store/mc/RunIII2024Summer24NanoAODv15/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/100000/1059354e-54e8-4beb-bb83-742e23d89675.root
#   python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --nthreads 16 -t _2024
#   python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --nthreads 16 -t _2024_HTT -f HTo2Tau
#   python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --nthreads 16 -t _2024_DY -f DY
#   python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --nthreads 16 -t _2024_TT -f TTto
# Instructions to submit jobs to HTCondor via Dask:
#   tmux # use tmux screen or whatever to wait for jobs to finish (write down lxplus node)
#   python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --dask -t _2024 # runs all
#   python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --dask -t _2024_HTT -f HTo2Tau
#   python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --dask -t _2024_DY -f DY
#   python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --dask -t _2024_TT -f TTto
#   # wait, or deattach the tmux screen to check condor_q
# Sources:
#   https://root.cern/doc/master/distrdf002__dask__connection_8py.html
#   https://batchdocs.web.cern.ch/specialpayload/dask.html
#   https://root-forum.cern.ch/t/distributed-rdataframe-fromspec/58123
#   https://root-forum.cern.ch/t/using-rdataframe-with-dask-and-rdatasetspec/58230
import time; TIME0 = time.perf_counter()
import os, sys, re, glob
import functools # for @functools.cache
import yaml # type: ignore
import numpy as np
from datetime import datetime
from fnmatch import fnmatch
sys.path.append("python")
from utils import loadconfig, warning, elapsed, BASE_DIR
from rocutils import getwpcuts
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True # to avoid conflict with argparse
from ROOT import gROOT, RDF, TFile
CONFIG = loadconfig(verb=0) # settings like cuts, binning, labels, ...

# CACHES
COL_DICT  = { # global dictionary to bookkeep all column definitions
  # colname: {
  #   'definition': ... # actual definition in RDataFrame
  #   'expanded':   ... # definition with expanded definitions
  # }
}
MASK_DICT = { # global dictionary of mask definitions (selections) for tau algorithms
  # algo: { # individual cuts for given tau type & ID algorithm
  #   'pt':   ... # require algo-specific pt cut (regressed for PNet/UParT)
  #   'id':   ... # require ID (and consequently corresponding reconstructed object)
  #   'all':  ... # require ID (and consequently corresponding reconstructed object)
  #   ...
  # }
}


def iter_gm(genmatches):
  """Help function to iterate over genmatches, return key and text label."""
  for gm in genmatches:
    yield gm, CONFIG['text']['gm'][gm]
  

@functools.cache
def getregexp(key):
  """Create regular expression. Cached to avoid recompiling and speed up."""
  regexp  = rf"\b{key}\b(?!\[)"
  regexp1 = re.compile(regexp)
  regexp2 = re.compile(rf'(\|\s*{regexp}|{regexp}\s*\|)')
  return regexp1, regexp2
  

def add2dict(colname,coldef,verb=0):
  """Add column definition to dictionary, including expanded form.
  This will be used for debugging the various mask definitons."""
  global COL_DICT
  colexp = coldef # column expression that will be fully expanded
  if colname.endswith('_elec') or colname.endswith('_jet'):
    return # ignore these fakes (_muon is enough)
  for key, value in COL_DICT.items():
    if key not in colexp:
      continue # speed up
    expval = value['expanded']
    regexp1, regexp2 = getregexp(key)
    if '|' in expval or regexp2.search(colexp):
      expval = f"({expval})"
    colexp = regexp1.sub(expval,colexp) # full expand expression
    if not any(s in colexp for s in ['|',' (',') ','[']):
      parts = [ ] # preserve order
      for cond in colexp.split(" & "):
        if cond not in parts: # deduplicate
          parts.append(cond)
      colexp = " & ".join(parts)
  if verb>=3:
    print(f">>> add2dict: Adding column {colname!r}={coldef!r}")
  COL_DICT[colname] = { 'definition': coldef, 'expanded': colexp }
  

def joinmasks_or(mask,algos,postfix=''):
  """Join masks as an OR."""
  if '$ALGO' not in mask:
    mask += '_$ALGO'
  if postfix:
    mask += postfix
  return (' | '.join(mask.replace('$ALGO',a) for a in algos))
  

def getcollectiontext(algo,postfix,extra=""):
  """Help function to create histogram model."""
  text = ""
  if 'hps' in postfix:
    key = 'hps3' if 'hps3' in postfix else 'hps2' if 'hps2' in postfix else 'hps'
    text = CONFIG['text']['collection'][key]
  elif 'qlty' in postfix:
    key = 'qlty_hps' if 'DeepTau' in algo else 'qlty_jet'
    text = CONFIG['text']['collection'][key]
  text = ', '.join(s for s in [text,extra] if s)
  return text
  

def gethistmodel(xvar,algo='',gm=None,postfix='',logit=None,xkey=None,xextra=None):
  """
  Help function to create histogram model for RDataFrame::Histo1D.
  Returns
  - xname:  full name of the RDataFrame column
  - htitle: histogram title
  - xtitle: x axis title
  - xbins:  binning, (nbins,xmin,xmax) for fixed, (nedges,array) for variable
  """
  
  # LOGIT
  if logit is None:
    isScore = any(t in xvar for t in ['VS','qConf'])
    logit = ('logit' in postfix and isScore)
  
  # CONSTRUCT htitle
  gmtext = CONFIG['text']['gm'].get(gm,"")
  htitle = ', '.join(x for x in [algo,gmtext] if x)
  
  # CONSTRUCT xname, xtitle, xbins
  xname  = f'{xvar}{algo}' # e.g. "ptPNetCHS_all_real"
  xkey   = xvar if xkey==None else xkey
  xtitle = CONFIG['title'].get(xkey,None)
  xbins  = CONFIG['binning'].get(xkey,None)
  if xvar=='ntau':
    xtitle = getcollectiontext(algo,postfix,xextra)
    xbins  = CONFIG['binning'].get('ntau',(2,0,2))
  elif xvar=='pt':
    ttype  = CONFIG['text']['regressAlgo'].get(algo,"")
    xtitle = f"{ttype} p_{{T}}"
    xbins  = CONFIG['binning'].get('pt',(50,0,250))
  elif 'ptres' in xvar: # ptresp or jptresp
    ###ttype = ('HPS' if 'DeepTau' in algo else 'CHS' if 'CHS' in algo else
    ###        'PUPPI' if 'PUPPI' in algo else algo)
    if 'jpt' in xvar:
      ttype = CONFIG['text']['recoAlgo'].get(algo,"")
      if 'HPS' not in ttype:
        ttype += " jet"
    else:
      ttype = CONFIG['text']['regressAlgo'].get(algo,"")
    xtitle = ttype + " " + CONFIG['title'].get('ptres',"p_{T} response")
    xbins  = CONFIG['binning'].get('ptresp',(80,0,2))
  elif 'VS' in xvar: # assume raw score
    xname  = f"raw{algo}{xvar}" # e.g. "rawPNetCHSVSjet_all_real"
    xtitle = f"Raw {algo}{xvar} score"
    xbins  = CONFIG['binning'].get('raw',(200,0,1)) # binning for linear raw output score in [0,1]
  if not xtitle or not xbins:
    raise IOError(f"gethistmodel: Did not find xtitle or xbins for {xvar=!r}, {algo=!r}, {xkey=!r}")
  
  # EXTRA TAGS
  if postfix:
    xname += "_"+postfix
  if gm:
    xname += '_'+gm
  if logit:
    xtitle = xtitle[0].lower()+xtitle[1:]
    xtitle = f"Logit of {xtitle}"
    xbins  = CONFIG['binning'].get('logit',(1000,-10,10))
  if xextra and xvar not in ['ntau']:
    xtitle = f"{xtitle} ({xextra})"
  
  # VARIABLE BINS
  if (nbins:=len(xbins))!=3 or isinstance(xbins,list):
    xbins = (nbins-1,np.array(xbins,'d'))
  
  return xname, htitle, xtitle.strip(), xbins
  

def createConnection(njobs=5,outdir='logs',nthreads=None,verb=0):
  """Create connection to setup Dask client."""
  from dask.distributed import Client
  from dask_lxplus import CernCluster
  import sys, socket
  if verb>=1:
    print(">>> createConnection: Setting up CernCluster...")
  #user = os.getenv('USER')
  #logdir = f"/afs/cern.ch/user/{user[0]}/{user}/test"
  uselcg = not all(k in os.environ for k in ['CMSSW_BASE','CMSSW_BASE'])
  cluster = CernCluster(
    cores=(1 if nthreads==None else nthreads),
    memory='2000MB', # NOTE: 1 CPU slot has <4GiB
    disk='1000MB',
    death_timeout='60',
    lcg=uselcg,
    nanny=False,
    container_runtime='none',
    log_directory=outdir, #os.path.abspath(outdir)
    scheduler_options={
      'port': 8786, # only one per lxplus node available
      'host': socket.gethostname(),
    },
    job_extra={
      '+MaxRuntime': 600, # seconds = 10 min.
      #'+JobFlavour': '"espresso"', # 20 min.
      #'+JobFlavour': '"microcentury"', # 60 min.
      '+AccountingGroup': '"group_u_CMST3"', # priority for EP-CMG members
    },
    #worker_command='distributed.cli.dask_worker', # prevent `KeyError: 'worker-command'`
    python=sys.executable, # prevent `KeyError: 'python'`
    extra=['--worker-port 10000:10100']
  )
  cluster.adapt(
    minimum=int(njobs),
    maximum=int(max(njobs,500))
  )
  if verb>=1:
    print(">>> createConnection: Script:")
    print(cluster.job_script())
  client = Client(cluster)
  return cluster, client
  

def compileMacro(verb=0):
  """Compile macros to load custom help functions like logit, reco2gm, gen2recodm, ptres"""
  macro = os.path.join(BASE_DIR,"macros/maps.cpp")
  if verb>=1:
    print(f">>> compileMacro: Compiling custom macro {macro}...")
  status = gROOT.LoadMacro(macro+'+O')
  assert status==0, f"Failed to load macro {macro}"
  

def getRDataFrame(infiles,nthreads=None,dask=False,verb=0):
  """
  Initiate RDataFrame.
  
  infiles:   list of nanoAOD ROOT files
  nthreads:  number of maximum number of threads for RDataFrame
  dask:      boolean, if True, Dask is used by RDataFrame to submit jobs to lxplus HTCondor
  verb:      int, verbosity level
  """
  if verb>=1:
    print(f">>> getRDataFrame: Creating RDataFrame ({nthreads=!r}, {dask=!r})...")
  if dask: # distributed (HTCondor via Dask)
    cluster, client = createConnection(njobs=len(infiles)/2,nthreads=nthreads,verb=verb)
    RDataFrame = RDF.Experimental.Distributed.Dask.RDataFrame
    dframe = RDataFrame('Events',infiles,daskclient=client)
  else: # run locally
    from ROOT import RDataFrame
    ROOT.EnableImplicitMT(8 if nthreads==None else nthreads) # set before creating RDataFrame
    dframe = RDataFrame('Events',infiles)
    if hasattr(RDF.Experimental,'AddProgressBar'):
      RDF.Experimental.AddProgressBar(dframe)
  return dframe
  

def getmaskdict(tauAlgos,wptag='',verb=0):
  """Construct dictionary of common masks (selections)."""
  discriminants = CONFIG['discriminants']
  vennWPs       = CONFIG['venn']['wps']  # VSjet WPs
  dmWPs         = CONFIG['dm']['wps']    # VSjet WPs
  allWPs        = set(vennWPs+dmWPs)     # VSjet WPs
  pt_cut        = CONFIG['cut']['pt']    # min. pt
  cut_qConf     = CONFIG['cut']['qConf'] # min. abs(qConf) for PNet/UParT
  
  # TODO? rename:
  #   all  -> allAlgoPtCut
  #   all2 -> allSamePtCut
  #   qlty -> qltyCharge
  #   hps  -> hpsIncl2Pr
  #   hps2 -> hpsExcl2Pr
  #   hps3 -> hpsExcl2PrPlusRecov
  
  # MASK dictionary defining individual and common cuts
  masks_dict = {
    # algo: { # individual cuts for given tau type & ID algorithm
    #   'pt':   ... # require algo-specific pt cut (regressed for PNet/UParT)
    #   'id':   ... # require ID (and consequently corresponding reconstructed object)
    #   'all':  ... # require ID (and consequently corresponding reconstructed object)
    #   'hps3': ... # HPS pT cut and no 2-prong DMs if HPS-matched, otherwise |qConf|>0.2 for recovery
    #   'qlty': ... # require ID + some quality cuts on prongness/charge
    #   wp:     ... # require WP
    # }
  }
  
  # INDIVDUAL CUTS per tau type (HPS/CHS/PUPPI), and tau ID algorithm (DeepTau/PNet/UParT)
  # Note: "quality" charge for PNet/UParT mean |qConf| > 0.2 to remove many jet fakes,
  #       while for DeepTau it means excluding 2-prong DMs
  for algo in tauAlgos:
    mask_pt    = f"pt{algo}>{pt_cut}" # pt depends on algorithm
    mask_id    = f"Tau_raw{algo}VSjet>-0.5" # require ID algo (and consequently tau type) exist
    mask_ptid  = f"{mask_pt} & {mask_id}"
    mask_no2p  = "mask_no2prong" #"Tau_decayMode!=5 & Tau_decayMode!=6"
    if 'DeepTau' in algo: # HPS+DeepTau
      mask_qlty  = f"{mask_id} & {mask_no2p}" # high-quality charge = excl. 2-prong DMs
      mask_hps3  = f"{mask_id} & {mask_no2p}" # excl. 2-prong DMs
    else: # AK4+PNet or AK4+UParT
      mask_qlty  = f"{mask_id} & abs(Tau_qConf{algo})>{cut_qConf}" # high-quality charge
      mask_hps3  = f"{mask_id} & {mask_no2p} & (Tau_idDecayModeNewDMs | abs(Tau_qConf{algo})>{cut_qConf})" # excl. 2-prong DMs
    masks_dict[algo] = {
      'pt':    mask_pt,   # algo-specific pT cut
      'id':    mask_id,   # ID, excl. pT cut
      'ptid':  mask_ptid, # ID, incl pT cut
      'qlty':  mask_qlty, # ID. + charge quality (excl. 2-prong DM for HPS, |qConf|>0.2 only for recovered taus)
      'hps3':  mask_hps3, # excl. 2-prong HPS DMs for all taggers, apply |qConf|>0.2 only for recovered taus
    }
    if algo in tauAlgos: # apply some reasonable WP
      wps = set(['VVVLoose','VLoose','Loose']) | allWPs # union
      wp2cut = getwpcuts(algo,outdir='output/wps',tag=wptag,discriminants=discriminants,wps=wps,defcut=0.05)
      lepcut = (f"Tau_raw{algo}VSmu>{wp2cut['VSmu']['Loose']['cut']:.5g} & " # second loosest VSmu WPs
                f"Tau_raw{algo}VSe>{wp2cut['VSe']['VVVLoose']['cut']:.5g}") # loosest VSe WPs
      for wp in allWPs: # e.g. 'VVVLoose' & 'Medium' VSjet
        wpcut = wp2cut['VSjet'][wp]['cut']
        masks_dict[algo][wp] = mask_id.replace("VSjet>-0.5",f"VSjet>{wpcut:.5g} & {lepcut}")
  if verb>=2:
    print(f">>> getmaskdict: {masks_dict=!r}")
  
  return masks_dict
  

def addColumnsMasks(dframe,tauAlgos,wptag='',verb=0):
  """
  Add mask columns for tau selection requirements.
  
  Returns dframe with all mask-related columns defined:
    - base masks (mask_gen, mask_eta, mask_real, mask_elec, etc.)
    - algo-specific selections (mask_all_*, mask_hps_*, etc.)
    - genmatch masks
    - combined cuts (mask_all, mask_qlty, mask_venn)
    - working point masks
  """
  global CONFIG, MASK_DICT
  genmatches    = CONFIG['genmatches']
  vennAlgos     = CONFIG['venn']['algos']
  vennWPs       = CONFIG['venn']['wps']
  dmWPs         = CONFIG['dm']['wps']
  pt_cut        = CONFIG['cut']['pt']
  eta_cut       = CONFIG['cut']['eta']
  
  # get algo-specific pt/charge/WP cuts
  MASK_DICT = getmaskdict(tauAlgos,wptag=wptag,verb=verb)
  
  def Define(colname,coldef):
    nonlocal dframe # get latest, local dframe
    add2dict(colname,coldef,verb=verb)
    dframe = dframe.Define(colname,coldef)
  
  if verb>=1:
    print(">>> addColumnsMasks: Defining mask columns...")
  
  # define pT by algo (HPS uses Tau_pt directly, PNet/UParT use regressed pT)
  for algo in tauAlgos:
    if 'DeepTau' in algo:
      Define(f"pt{algo}","Tau_pt")
    else:
      Define(f"pt{algo}", f"Tau_pt*Tau_ptCorr{algo}")
  
  # define common masks
  Define('mask_gen',     f"GenVisTau_pt>{pt_cut} & abs(GenVisTau_eta)<{eta_cut}")
  Define('mask_eta',     f"abs(Tau_eta)<{eta_cut} & abs(Tau_dz)<0.2")
  Define('mask_no2prong',"Tau_decayMode!=5 & Tau_decayMode!=6")
  Define('mask_all2',    f"mask_eta & Tau_pt>{pt_cut}")
  Define('mask_hps',     f"mask_all2 & Tau_idDecayModeNewDMs")
  Define('mask_hps2',    f"mask_hps & mask_no2prong")
  Define('mask_real',    "Tau_genPartFlav==5") # genuine taus
  Define('mask_elec',    "Tau_genPartFlav==1 | Tau_genPartFlav==3") # electron fakes
  Define('mask_muon',    "Tau_genPartFlav==2 | Tau_genPartFlav==4") # muons fakes
  Define('mask_jet',     "Tau_genPartFlav==0") # q/g jet fakes
  Define('gm',           "reco2gm(Tau_genPartFlav, Tau_genPartIdx, GenVisTau_status)") # map to []
  
  # define algo-specific masks
  for algo in tauAlgos:
    Define(f'mask_all_{algo}',  f"mask_eta & {MASK_DICT[algo]['ptid']}")
    Define(f'mask_all2_{algo}', f"mask_all2 & {MASK_DICT[algo]['id']}")
    Define(f'mask_hps_{algo}',  f"mask_hps & {MASK_DICT[algo]['id']}")
    Define(f'mask_hps2_{algo}', f"mask_hps2 & {MASK_DICT[algo]['qlty']}")
    Define(f'mask_hps3_{algo}', f"mask_all2 & {MASK_DICT[algo]['hps3']}")
    Define(f'mask_qlty_{algo}', f"mask_all_{algo} & {MASK_DICT[algo]['qlty']}") # charge quality
    
    # genmatch combinations
    for gm in genmatches:
      Define(f'mask_all_{algo}_{gm}',  f"mask_all_{algo} & mask_{gm}")
      Define(f'mask_all2_{algo}_{gm}', f"mask_all2_{algo} & mask_{gm}")
      Define(f'mask_hps_{algo}_{gm}',  f"mask_hps_{algo} & mask_{gm}")
      Define(f'mask_hps2_{algo}_{gm}', f"mask_hps2_{algo} & mask_{gm}")
      Define(f'mask_hps3_{algo}_{gm}', f"mask_hps3_{algo} & mask_{gm}")
      Define(f'mask_qlty_{algo}_{gm}', f"mask_qlty_{algo} & mask_{gm}")
  
  # combined cuts (OR across algos)
  # Note: Needed for common denomenator.
  Define('mask_all',  joinmasks_or('mask_all_$ALGO',tauAlgos))
  Define('mask_qlty', joinmasks_or('mask_qlty_$ALGO',tauAlgos))
  
  # Venn diagram base masks
  # Note:
  #   mask_venn(_gm) is the union of the three algos in vennAlgos
  #   Needed for common denomenator.
  Define('mask_venn', joinmasks_or('mask_hps3_$ALGO',vennAlgos))
  for gm in genmatches:
    Define(f'mask_venn_{gm}', f"mask_venn & mask_{gm}")
  
  # common WP masks
  for wp in set(vennWPs+dmWPs):
    for algo in tauAlgos:
      Define(f'mask_{algo}_{wp}',      MASK_DICT[algo][wp])
      Define(f'mask_hps3_{algo}_{wp}', f"mask_hps3_{algo} & mask_{algo}_{wp}")
  
  # WP masks for just the Venn diagrams
  for wp in vennWPs:
    Define(f'mask_venn_{wp}',          joinmasks_or(f'mask_hps3_$ALGO_{wp}',vennAlgos))
    for gm in genmatches:
      Define(f'mask_venn_{wp}_{gm}',   f"mask_venn_{wp} & mask_{gm}")
  
  # control variables
  for gm in genmatches:
    
    # number of tau candidates per type & genmatch
    Define(f'ntau_all_{gm}',   f"mask_all[mask_{gm}]")
    Define(f'ntau_all2_{gm}',  f"mask_all2[mask_{gm}]")
    Define(f'ntau_hps_{gm}',   f"mask_hps[mask_{gm}]")
    Define(f'ntau_hps2_{gm}',  f"mask_hps2[mask_{gm}]")
    Define(f'ntau_qlty_{gm}',  f"mask_qlty[mask_{gm}]")
    
    # reconstructed pT and decay mode (inclusive)
    Define(f'pt_all2_{gm}',    f"Tau_pt[mask_all2 & mask_{gm}]") # inclusive pT
    Define(f'recodm_all2_{gm}',f"Tau_decayMode[mask_all2 & mask_{gm}]")
  
  return dframe
  

def addColumnsRawScores(dframe,tauAlgos,verb=0):
  """
  Add raw discriminant score columns for all algorithms.
  
  Includes linear raw scores and logit-transformed versions for fine binning.
  """
  if verb>=1:
    print(">>> addColumnsRawScores: Defining raw score columns...")
  discriminants = CONFIG['discriminants']
  genmatches    = CONFIG['genmatches']
  xclip         = CONFIG['binning']['logit'][-1] # clip logit, e.g. between [-9,+9]
  
  def Define(colname,coldef):
    nonlocal dframe # get latest, local dframe
    add2dict(colname,coldef,verb=verb)
    dframe = dframe.Define(colname,coldef)
  
  for algo in tauAlgos:
    for disc in discriminants:
      for gm in genmatches:
        xvar   = f"raw{algo}{disc}"
        
        # linear scores [0,1]
        Define(f'{xvar}_all_{gm}',        f"Tau_raw{algo}{disc}[mask_all_{algo}_{gm}]")
        Define(f'{xvar}_all2_{gm}',       f"Tau_raw{algo}{disc}[mask_all2_{algo}_{gm}]")
        Define(f'{xvar}_hps_{gm}',        f"Tau_raw{algo}{disc}[mask_hps_{algo}_{gm}]")
        Define(f'{xvar}_hps2_{gm}',       f"Tau_raw{algo}{disc}[mask_hps2_{algo}_{gm}]")
        Define(f'{xvar}_hps3_{gm}',       f"Tau_raw{algo}{disc}[mask_hps3_{algo}_{gm}]")
        Define(f'{xvar}_qlty_{gm}',       f"Tau_raw{algo}{disc}[mask_qlty_{algo}_{gm}]")
        
        # logit-transformed (for zooming into peaks at 0 and 1)
        Define(f'{xvar}_logit_all_{gm}',  f"logit({xvar}_all_{gm}, {xclip})")
        Define(f'{xvar}_logit_all2_{gm}', f"logit({xvar}_all2_{gm},{xclip})")
        Define(f'{xvar}_logit_hps_{gm}',  f"logit({xvar}_hps_{gm}, {xclip})")
        Define(f'{xvar}_logit_hps2_{gm}', f"logit({xvar}_hps2_{gm},{xclip})")
        Define(f'{xvar}_logit_hps3_{gm}', f"logit({xvar}_hps3_{gm},{xclip})")
        Define(f'{xvar}_logit_qlty_{gm}', f"logit({xvar}_qlty_{gm},{xclip})")
  
  return dframe
  

def addColumnsKinematics(dframe,tauAlgos,verb=0):
  """Add kinematic columns: pT, charge probability (qConf), decay modes, and pT response."""
  if verb>=1:
    print(">>> addColumnsKinematics: Defining kinematic columns...")
  genmatches    = CONFIG['genmatches']
  dmWPs         = CONFIG['dm']['wps']
  dmPostfixes   = CONFIG['dm']['postfixes']
  dmPostfixesWP = CONFIG['dm']['postfixes_wp']
  dmPostfixesGM = CONFIG['dm']['postfixes_gm']
  vennWPs       = CONFIG['venn']['wps']
  xclip         = CONFIG['binning']['logit'][-1]
  
  def Define(colname,coldef):
    nonlocal dframe # get latest, local dframe
    add2dict(colname,coldef,verb=verb)
    dframe = dframe.Define(colname,coldef)
  
  for algo in tauAlgos:
    for gm in genmatches:
      
      # charge probability (qConf)
      # - qConf     = inclusive (HPS+recovery)
      # - qConf_hps = HPS only
      # - qConf_rec = recovery only (recovered by PNet/UParT w.r.t. HPS)
      # Note:
      #   Remove the qConf cuts from "hps3" to get the full qConf spectrum, but still exclude 2 prongs
      if 'DeepTau' not in algo:
        mask_all2_noQConf = f"mask_all2_{algo}_{gm} & mask_no2prong" # remove qConf cut
        Define(f'qConf{algo}_all_{gm}',       f"Tau_qConf{algo}[mask_all_{algo}_{gm}]")
        Define(f'qConf{algo}_all2_{gm}',      f"Tau_qConf{algo}[mask_all2_{algo}_{gm}]")
        Define(f'qConf{algo}_hps_{gm}',       f"Tau_qConf{algo}[mask_hps & mask_{gm}]")
        Define(f'qConf{algo}_hps2_{gm}',      f"Tau_qConf{algo}[mask_hps2 & mask_{gm}]")
        Define(f'qConf{algo}_hps3_{gm}',      f"Tau_qConf{algo}[{mask_all2_noQConf}]")
        Define(f'qConf{algo}_hps_hps3_{gm}',  f"Tau_qConf{algo}[{mask_all2_noQConf} & Tau_idDecayModeNewDMs]")
        Define(f'qConf{algo}_rec_hps3_{gm}',  f"Tau_qConf{algo}[{mask_all2_noQConf} & !Tau_idDecayModeNewDMs]")
        Define(f'qConf{algo}_logit_all_{gm}', f"logit(0.5+qConf{algo}_all_{gm},{xclip})")
      
      # pT (momentum)
      # Note:
      #   Remove pT cut for pt{algo}_all to get full spectrum before cut
      mask_all_noPT = f"mask_eta & {MASK_DICT[algo]['id']} & mask_{gm}" # remove pT cut
      Define(f'pt{algo}_all_{gm}',  f"pt{algo}[{mask_all_noPT}]")       # regressed pT
      Define(f'pt{algo}_all2_{gm}', f"Tau_pt[mask_all2_{algo}_{gm}]")   # HPS pT
      Define(f'pt{algo}_hps_{gm}',  f"Tau_pt[mask_hps_{algo}_{gm}]")    # HPS pT
      Define(f'pt{algo}_hps3_{gm}', f"Tau_pt[mask_hps3_{algo}_{gm}]")   # HPS pT
      Define(f'pt{algo}_qlty_{gm}', f"pt{algo}[mask_qlty_{algo}_{gm}]") # regressed pT
    
    # pT response (regressed pT / gen visible pT) for genuine taus
    # - ptres     = inclusive (HPS+recovery)
    # - ptres_hps = HPS only
    # - ptres_rec = recovery only (recovered by PNet/UParT w.r.t. HPS)
    wp = vennWPs[0] # loosest
    ptres  = lambda m: f"ptres(GenVisTau_pt, pt{algo}[{m}], Tau_genPartIdx[{m}])"
    Define(f'ptres{algo}_qlty',          ptres(f'mask_qlty_{algo}')+"[mask_gen]")
    Define(f'ptres{algo}_hps_hps3',      ptres(f'mask_hps3_{algo}_real & Tau_idDecayModeNewDMs')+"[mask_gen]")
    Define(f'ptres{algo}_hps_hps3_{wp}', ptres(f'mask_hps3_{algo}_{wp} & Tau_idDecayModeNewDMs & mask_real')+"[mask_gen]")
    if 'DeepTau' not in algo:
      Define(f'ptres{algo}_rec_hps3',      ptres(f'mask_hps3_{algo}_real & !Tau_idDecayModeNewDMs')+"[mask_gen]")
      Define(f'ptres{algo}_rec_hps3_{wp}', ptres(f'mask_hps3_{algo}_{wp} & !Tau_idDecayModeNewDMs & mask_real')+"[mask_gen]")
    
    # decay mode vs. genmatch
    # Note:
    # - recodm, reco2gm, hpsdm have dimension of reconstructed tau candidates
    # - gendm, gen2recodm      have dimension of genuine (hadronic) taus
    dm = "Tau_decayMode"
    if 'DeepTau' not in algo:
      dm += algo # e.g. Tau_decayModePNetCHS
    postfixes = [(None,pf) for pf in dmPostfixes[algo]] +\
                [(wp,pf) for pf in dmPostfixesWP[algo] for wp in dmWPs] # add WPs
    for wp, postfix in postfixes:
      if wp is None: # inclusive
        mask = f"mask_{postfix}_{algo}"
      else: # apply WP
        mask = f"mask_{algo}_{wp} & mask_{postfix}_{algo}"
        postfix = f"{postfix}_{wp}"
      gen2recodm = f"gen2recodm(nGenVisTau, Tau_genPartIdx[{mask}], recodm{algo}_{postfix})"
      Define(f'recodm{algo}_{postfix}',     f"{dm}[{mask}]")
      Define(f'reco2gm{algo}_{postfix}',    f"gm[{mask}]")
      Define(f'gendm{algo}_{postfix}',      f"GenVisTau_status[mask_gen]")
      Define(f'gen2recodm{algo}_{postfix}', f"{gen2recodm}[mask_gen]")
      
      # split recodm by genmatch
      if postfix in dmPostfixesGM[algo]:
        for gm in genmatches:
          hpsdm = f'recodm{algo}_{postfix}_{gm}' if 'DeepTau' in algo else f"Tau_decayMode[{mask}_{gm}]"
          Define(f'recodm{algo}_{postfix}_{gm}',f"{dm}[{mask}_{gm}]")
          Define(f'hpsdm{algo}_{postfix}_{gm}', hpsdm) # use HPS Tau_decayMode explicitly
  
  return dframe
  

def addColumnsVenn(dframe,verb=0):
  """Add Venn diagram mask columns for studying algorithm overlap and complementarity."""
  if verb>=1:
    print(">>> addColumnsVenn: Defining Venn diagram columns...")
  genmatches = CONFIG['genmatches']
  vennAlgos  = CONFIG['venn']['algos']
  vennWPs    = CONFIG['venn']['wps']
  
  def Define(colname,coldef):
    nonlocal dframe # get latest, local dframe
    add2dict(colname,coldef,verb=verb)
    dframe = dframe.Define(colname,coldef)
  
  for algo in vennAlgos:
    for gm in genmatches:
      Define(f'ntau_venn_{algo}_{gm}', f"mask_hps3_{algo}[mask_venn_{gm}]")
      for wp in vennWPs:
        Define(f'ntau_venn_{algo}_{wp}_{gm}', f"mask_hps3_{algo}_{wp}[mask_venn_{wp}_{gm}]")
  
  return dframe
  

def addColumns(dframe,tauAlgos,wptag='',maskfname='masks.yaml',verb=0):
  """
  Initiate RDataFrame, & add useful columns to RDataFrame.
  
  dframe:        RDataFrame instance
  tauAlgos:      list of tau ID algorithms, e.g. ['DeepTau','PNet']
  discriminants: list of distriminants, e.g. ['VSjet','VSmu','VSe']
  genmatches:    list of genmatches, e.g. ['real','jet']
  wptag:         string, tag for YAML file with cuts for WPs ('{outdir}/wps_{algo}{wptag}.yaml"')
  maskfname:     string, file name for YAML file with masks used by RDataFrame
  verb:          int, verbosity level
  """
  if verb>=1:
    print(">>> addColumns: Adding columns in semantic modules...")
  
  # 0. Load macros like logit, reco2gm, gen2recodm, ptres
  compileMacro(verb=verb)
  
  # 1. Add all mask columns (selections, genmatch combinations, WPs)
  dframe = addColumnsMasks(dframe,tauAlgos,wptag=wptag,verb=verb)
  
  # 2. Add raw discriminant scores (linear and logit-transformed)
  dframe = addColumnsRawScores(dframe,tauAlgos,verb=verb)
  
  # 3. Add kinematics (pT, charge, decay mode, pT response)
  dframe = addColumnsKinematics(dframe,tauAlgos,verb=verb)
  
  # 4. Add Venn diagram columns (algorithm overlap study)
  dframe = addColumnsVenn(dframe,verb=verb)
  
  # WRITE masks to YAML
  if maskfname:
    global COL_DICT, MASK_DICT
    if verb>=1:
      nmasks = sum(len(v) for v in MASK_DICT.values())
      print(f">>> addColumns: Writing {nmasks} masks & {len(COL_DICT)} columns to {maskfname}...")
    MASK_DICT['common'] = { p: COL_DICT[f'mask_{p}']['expanded'] for p in ['all','all2','hps','hps2','qlty','venn'] if f'mask_{p}' in COL_DICT }
    MASK_DICT['columns_masks'] = { k: v for k, v in COL_DICT.items() if k.startswith('mask') }
    MASK_DICT['columns_vars']  = { k: v for k, v in COL_DICT.items() if not k.startswith('mask') }
    with open(maskfname,'w') as outfile:
      outfile.write(f"# Generated by {__file__}.addColumns on {datetime.now():%Y-%m-%d %H:%M:%S}\n")
      yaml.dump(MASK_DICT, outfile, default_flow_style=False, sort_keys=False)
  
  return dframe
  

def addHists(dframe,tauAlgos,verb=0):
  """
  Add histograms to RDataFrame.
  """
  dmWPs         = CONFIG['dm']['wps']
  dmPostfixes   = CONFIG['dm']['postfixes']
  dmPostfixesWP = CONFIG['dm']['postfixes_wp']
  dmPostfixesGM = CONFIG['dm']['postfixes_gm']
  vennAlgos     = CONFIG['venn']['algos']
  vennWPs       = CONFIG['venn']['wps']
  discriminants = CONFIG['discriminants']
  genmatches    = CONFIG['genmatches']
  recocuts      = CONFIG['text']['collection']['all'] # baseline cuts
  ytitle        = CONFIG['title']['ngentau']
  ntitle        = CONFIG['title']['ntau'] # number of tau candidates
  results = { } # dict of histogram -> directory name
  
  # HELP FUNCTION
  def addHist(dname,hmodel,*vars):
    if verb>=1: print(f">>>   Adding {hmodel!r}, {vars=}")
    method = getattr(dframe,f"Histo{len(vars)}D") # dframe.Histo1D/2D/3D
    results.setdefault(dname,[ ]).append(method(hmodel,*vars))
  
  # DEFINE 1D HISTOGRAMS: counting gen or reco taus
  # - 'all'   = everything passing basic reco cuts, algo-dependent pt cut
  # - 'all2'  = everything passing basic reco cuts, common pt cut
  # - 'qlty'  = with quality cuts (no 2-prong for HPS, |qConf| > 0.2 for PNet/UParT)
  # - 'hps'   = HPS incl. 2-prong DMs, w/o qConf cut
  # - 'hps2'  = HPS excl. 2-prong DMs, w/ qConf cut
  # - 'hps3'  = HPS excl. 2-prong HPS DMs, and |qConf| > 0.2 if recovered by PNet/UParT (no HPS match)
  # - 'logit' = 2*tanh(2x-1) transform to zoom into peaks at 0 and 1 (plus fine binning)
  htitle        = f";{recocuts};{ntitle}"
  addHist('ntau',('ntau_gen',htitle,2,0,2),'mask_gen')
  for postfix in ['all','all2','hps','hps2','qlty']:
    for gm in genmatches:
      xvar, htitle, xtitle, xbins = gethistmodel('ntau','',gm,postfix,xextra=recocuts)
      htitle = f"{htitle};{xtitle};{ntitle}"
      addHist('ntau',(xvar,htitle,*xbins),xvar)
  
  # DEFINE 1D HISTOGRAMS: pt (for common denominator)
  for postfix in ['all2']:
    for gm in genmatches:
      xvar, htitle, xtitle, xbins = gethistmodel('pt','',gm,postfix)
      htitle = f"{htitle};{xtitle};{ntitle}"
      addHist('pt',(xvar,htitle,*xbins),xvar)
  
  # DEFINE 1D HISTOGRAMS: DM (for common denominator)
  for postfix in ['all2']:
    for gm in genmatches:
      xvar, htitle, xtitle, xbins = gethistmodel('recodm','',gm,postfix,xkey='dm')
      hname  = xvar.replace('recodm','dm')
      htitle = f"{htitle};{xtitle};{ntitle}"
      addHist('dm',(hname,htitle,*xbins),xvar)
  
  for algo in tauAlgos:
    rawdir = 'raw/'+algo
    qcdir  = 'charge/'+algo
    ptdir  = 'pt/'+algo
    dmdir  = 'dm/'+algo
    
    # DEFINE 1D HISTOGRAMS: number of tau candidates
    for postfix in ['all','hps','hps2','qlty']:
      xvar, htitle, xtitle, xbins = gethistmodel('ntau',algo,None,postfix,xextra=recocuts)
      htitle = f"{htitle};{xtitle};{ntitle}"
      addHist('ntau',(xvar,htitle,2,0,2),f"mask_{postfix}_{algo}")
    
    # DEFINE 1D HISTOGRAMS: pt
    for postfix in ['all','all2','qlty','hps']:
      for gm in genmatches:
        xvar, htitle, xtitle, xbins = gethistmodel('pt',algo,gm,postfix,xextra=recocuts)
        htitle = f"{htitle};{xtitle};{ntitle}"
        addHist(ptdir,(xvar,htitle,*xbins),xvar)
    
    # DEFINE 1D HISTOGRAMS: pt response
    wp = vennWPs[0] # take loosest (assume first WP in list)
    for postfix in ['qlty','hps_hps3',f'hps_hps3_{wp}','rec_hps3',f'rec_hps3_{wp}']:
      if 'DeepTau' in algo and 'rec' in postfix:
        continue
      xvar, htitle, xtitle, xbins = gethistmodel('ptres',algo,None,postfix,xextra=recocuts)
      htitle = f"{htitle};{xtitle};{ntitle}"
      addHist(ptdir,(xvar,htitle,*xbins),xvar)
    
    # DEFINE 1D HISTOGRAMS: raw output scores
    postfixes = [ 'logit_all','logit_all2','hps','logit_hps','logit_hps2',
                 'hps3','logit_hps3','qlty','logit_qlty' ]
    for postfix in postfixes:
      for disc in discriminants:
        for gm in genmatches:
          xvar, htitle, xtitle, xbins = gethistmodel(disc,algo,gm,postfix)
          htitle = f"{htitle};{xtitle};{ntitle}"
          hmodel = (xvar,htitle,*xbins)
          if verb>=1:
            print(f">>>   Adding {hmodel!r}")
          addHist(rawdir,hmodel,xvar)
    
    #### DEFINE 1D HISTOGRAMS: jet pt response
    ###if 'DeepTau' not in algo:
    ###  for postfix in ['qlty','hps_hps3','rec_hps3']:
    ###    xvar, htitle, xtitle, xbins = gethistmodel('jptres',algo,postfix=postfix,xextra=recocuts)
    ###    htitle = f"{algo};{xtitle};Number of reconstructed tau candidates"
    ###    addHist(ptdir,(xvar,htitle,*xbins),xvar))
    
    # DEFINE 1D HISTOGRAMS: charge quality probability (qConf)
    if 'DeepTau' not in algo:
      postfixes = ['all','all2','logit_all','hps','hps2','hps3','hps_hps3','rec_hps3']
      for postfix in postfixes:
        for gm in genmatches:
          xvar, htitle, xtitle, xbins = gethistmodel('qConf',algo,gm,postfix,xextra=recocuts)
          htitle = f"{htitle};{xtitle};{ntitle}"
          addHist(qcdir,(xvar,htitle,*xbins),xvar)
    
    # DEFINE 2D HISTOGRAMS: VSmu/e vs. VSjet to look for correlations
    xdisc = 'VSjet'
    postfix = 'logit_qlty' # arctanh transform, with quality cuts
    if xdisc in discriminants:
      for ydisc in discriminants:
        for gm in genmatches:
          if xdisc==ydisc:
            continue
          xvar, htitle, xtitle, xbins = gethistmodel(xdisc,algo,gm,postfix)
          yvar, htitle, ytitle, ybins = gethistmodel(ydisc,algo,gm,postfix)
          hname  = f"raw{algo}{xdisc}_{ydisc}_{postfix}_{gm}"
          htitle = f"{htitle};{xtitle};{ytitle};{ntitle}"
          hmodel = (hname,htitle,*xbins,*ybins)
          addHist(rawdir,hmodel,xvar,yvar)
    
    #### DEFINE 3D HISTOGRAMS: VSe vs. VSmu vs. VSjet to look for correlations
    #### WARNING: This may lead to crashes, likely due to large number of bins (1000^3).
    ####          Might be fixed in ROOT >= v6.32 ?
    ###xdisc, ydisc, zdisc = 'VSjet', 'VSmu', 'VSe'
    ###for postfix in ['logit_hps3','logit_qlty']: # arctanh transform, with quality cuts
    ###  if xdisc in discriminants:
    ###    for gm in genmatches:
    ###      xvar, htitle, xtitle, xbins = gethistmodel(xdisc,algo,gm,postfix)
    ###      yvar, htitle, ytitle, ybins = gethistmodel(ydisc,algo,gm,postfix)
    ###      zvar, htitle, ztitle, zbins = gethistmodel(zdisc,algo,gm,postfix)
    ###      hname  = f"{postfix}_raw{algo}{xdisc}_{ydisc}_{zdisc}_{gm}"
    ###      htitle = f"{gm};{xtitle};{ytitle};{ztitle}"
    ###      hmodel = (hname,htitle,*xbins,*ybins,*zbins)
    ###      addHist(rawdir,hmodel,xvar,yvar,zvar)
    
    # DEFINE 2D HISTOGRAMS: qConf vs. VSjet to look for correlations
    postfix = 'logit_all' # arctanh transform, w/o quality cuts
    if xdisc in discriminants and 'DeepTau' not in algo:
      for gm in genmatches:
        xvar, htitle, xtitle, xbins = gethistmodel(xdisc,  algo,gm,postfix)
        yvar, htitle, ytitle, ybins = gethistmodel('qConf',algo,gm,postfix)
        hname  = f"raw{algo}{xdisc}_qConf_{postfix}_{gm}"
        htitle = f"{htitle};{xtitle};{ytitle};{ntitle}"
        hmodel = (hname,htitle,*xbins,*ybins)
        addHist(qcdir,hmodel,xvar,yvar)
    
    # DEFINE 2D HISTOGRAMS: pT vs. VSjet
    for ypostfix in ['hps3','qlty']:
      xpostfix = 'logit_'+ypostfix
      for xdisc in discriminants:
        for gm in genmatches:
          xvar, htitle, xtitle, xbins = gethistmodel(xdisc,algo,gm,xpostfix)
          yvar, htitle, ytitle, ybins = gethistmodel('pt', algo,gm,ypostfix,xkey='pt2')
          hname  = f"raw{algo}{xdisc}_pt_{xpostfix}_{gm}"
          htitle = f"{htitle};{xtitle};{ytitle};{ntitle}"
          hmodel = (hname,htitle,*xbins,*ybins)
          addHist(ptdir,hmodel,xvar,yvar)
    
    # DEFINE 2D HISTOGRAMS: reco DM vs. gen DM/GENMATCH
    postfixes = dmPostfixes[algo] +\
                [f"{pf}_{wp}" for pf in dmPostfixesWP[algo] for wp in dmWPs] # add WPs
    for postfix in postfixes:
      
      # DEFINE 2D HISTOGRAMS: reco vs. gen decay mode to study migration
      # Note: recodm, reco2gm, hpsdm have dimension of reconstructed tau candidates
      xvar, htitle, xtitle, xbins = gethistmodel('reco2gm',algo,None,postfix,xkey='gm')
      yvar, htitle, ytitle, ybins = gethistmodel('recodm', algo,None,postfix,xkey='dm')
      hname  = f"dm{algo}_gm_{postfix}"
      htitle = f"{htitle};{xtitle};{ytitle};{ntitle}"
      hmodel = (hname,htitle,*xbins,*ybins)
      addHist(dmdir,hmodel,xvar,yvar)
      
      # DEFINE 2D HISTOGRAMS: reco DM vs. genmatch to study migration
      # Note: gendm, gen2recodm have dimension of genuine (hadronic) taus
      xvar, htitle, xtitle, xbins = gethistmodel('gendm',algo,None,postfix,xkey='gendm')
      yvar, htitle, ytitle, ybins = gethistmodel('gen2recodm',algo,None,postfix,xkey='dm')
      hname  = f"dm{algo}_gendm_{postfix}"
      htitle = f"{htitle};{xtitle};{ytitle};{ntitle}"
      hmodel = (hname,htitle,*xbins,*ybins)
      addHist(dmdir,hmodel,xvar,yvar)
    
    # DEFINE 2D HISTOGRAMS: HPS DM vs. VSjet
    for ypostfix in dmPostfixesGM[algo]:
      xpostfix = f"logit_{ypostfix}"
      for xdisc in discriminants:
        for gm in genmatches:
          xvar, htitle, xtitle, xbins = gethistmodel(xdisc,  algo,gm,xpostfix)
          yvar, htitle, ytitle, ybins = gethistmodel('hpsdm',algo,gm,ypostfix,xkey='hpsdm')
          hname  = f"raw{algo}{xdisc}_dm_{xpostfix}_{gm}"
          htitle = f"{htitle};{xtitle};{ytitle};{ntitle}"
          hmodel = (hname,htitle,*xbins,*ybins)
          addHist(dmdir,hmodel,xvar,yvar)
  
  # DEFINE 3D VENN HISTOGRAMS
  if len(vennAlgos)==3:
    for gm in genmatches:
      for wp in ['no WP']+vennWPs:
        _wp    = '' if wp=='no WP' else ('_'+wp)
        hname  = f"ntau_venn{_wp}_{gm}"
        gmtext = CONFIG['text']['gm'][gm]
        htitle = f"Venn histogram {wp} ({gmtext});{';'.join(vennAlgos)}" #;{ntitle}
        hmodel = (hname,htitle,*CONFIG['binning']['venn'])
        xvars  = (f"ntau_venn_{algo}{_wp}_{gm}" for algo in vennAlgos)
        addHist('venn',hmodel,*xvars)
  
  return dframe, results
  

def parseInFiles(infname,url="root://xrootd-cms.infn.it/",filters=None,verb=0):
  """
  Parse user's list of input files.
  If ASCII txt file, parse file's content to list of ROOT files.
  """
  infnames = [ ]
  if isinstance(infname,list):
    for fname in infname:
      infnames.extend(parseInFiles(fname,filters=filters))
    return infnames
  assert isinstance(infname,str), f"Expected string: infname={infname!r}"
  if infname.lower().endswith('.root'):
    return glob.glob(infname) if '*' in infname else [infname]
  with open(infname, 'r') as txtfile:
    for line in txtfile:
      line = line.strip() # remove leading or trailing spaces
      if line and line[0]!='#' and line.lower().endswith('.root'):
        if filters and not any(fnmatch(line,f) for f in filters):
          continue # no match
        if line.startswith('/store/mc/'):
          line = url+line # add url
        if '*' in line: # expand wild card
          infnames.extend(glob.glob(line))
        else: # add single file
          infnames.append(line)
      elif verb>=1:
        print(f">>> parseInFiles: Ignoring line={line}...")
  return infnames
  

def writeOutput(outfname,res_dict,verb=0):
  """Write histograms to an output ROOT file."""
  outfile = TFile.Open(outfname,'RECREATE')
  for dname, hists in res_dict.items():
    tdir = outfile
    if dname!='.': # subdirectory
      for subdir in dname.split('/'):
        tdir = tdir.mkdir(subdir,subdir,True)
    tdir.cd()
    if verb>=1:
      print(f">>>   Writing {len(hists)} histograms to {outfname}:/{dname}...")
    for hist in hists: # write
      if verb>=2:
        print(f">>>     Writing {outfname}:/{dname}/{hist.GetName()}...")
      hist.GetValue().Write(hist.GetName())
  outfile.Close()
  

def main(args):
  global CONFIG
  verbosity = args.verbosity
  filters   = ['*'+f+'*' for f in args.filters] if args.filters else None
  infnames  = parseInFiles(args.inputFile,filters=filters)
  tauAlgos  = args.tauAlgos
  outdir    = args.outdir
  tag       = args.tag
  wptag     = args.wptag
  outfname  = os.path.join(outdir,args.outputFile.replace('$TAG',tag))
  maskfname = os.path.splitext(outfname)[0]+'.yaml'
  usedask   = args.dask
  dryrun    = args.dryrun
  vennAlgos = CONFIG['venn']['algos']
  if outdir:
    os.makedirs(outdir,exist_ok=True)
  if any(a not in tauAlgos for a in vennAlgos):
    raise IOError(f"{vennAlgos=!r} not contained in {tauAlgos=!r}")
  
  # CREATE RDATAFRAME
  print(f">>> Initiating RDataFrame for {len(infnames)} files...")
  dframe = getRDataFrame(infnames,dask=usedask,verb=verbosity)
  
  # DEFINE COLUMNS
  print(">>> Defining columns...")
  dframe = addColumns(dframe,tauAlgos,wptag=wptag,maskfname=maskfname,verb=verbosity)
  
  # DEFINE HISTOGRAMS
  print(">>> Defining histograms...")
  dframe, res_dict = addHists(dframe,tauAlgos,verb=verbosity)
  
  # RUN EVENT LOOP
  if dryrun:
    return
  print(">>> Start RDataFrame event loop...")
  res_list = [x for l in res_dict.values() for x in l]
  time1 = time.perf_counter()
  if usedask:
    if res_list:
      res_list[0].GetValue() # trigger event loop
  else:
    RDF.RunGraphs(res_list)
  print(f">>> RDataFrame finished in {elapsed(time1)}")
  
  # WRITE TO FILE
  print(f">>> Writing {len(res_list)} histograms to {outfname} ...")
  writeOutput(outfname,res_dict,verb=verbosity)
  if verbosity>=1:
    print(f">>> Written {len(res_list)} histograms to {outfname} !")
  

if __name__ == '__main__':
  from argparse import ArgumentParser
  description = """This script compute efficiencies of tau ID algorithms."""
  parser  = ArgumentParser(description=description,epilog="Good luck!")
  parser.add_argument('inputFile',           nargs='+',
                      metavar='INPUT_FILE',  help="input ROOT file(s) or txt file with list of ROOT files" )
  parser.add_argument('-o', '--output',      dest='outputFile', default="tauAlgoEffs$TAG.root",
                      metavar='OUTPUT_FILE', help="output ROOT file, default=%(default)r" )
  parser.add_argument('-O', '--outdir',      default="output",
                                             help="output directory, default=%(default)r" )
  parser.add_argument('-A', '--tauAlgos',    nargs='+', choices=CONFIG['tauAlgos'], default=CONFIG['tauAlgos'],
                                             help="list of tau algos, default=%(default)s" )
  parser.add_argument('-t', '--tag',         default="",
                                             help="tag for default output ROOT file name" )
  parser.add_argument('-w', '--wptag',       default=CONFIG['wp']['tag'], #"_2024_hps3",
                                             help="tag for WP YAML, e.g. $output/wps_$algo$wptag.yaml" )
  parser.add_argument('-f', '--filter',      dest='filters', nargs='+',
                                             help="glob pattern to filter ROOT file name from input list" )
  parser.add_argument('-d', '--dask',        action='store_true',
                                             help="use Dask with lxplus HTCondor to submit jobs" )
  parser.add_argument('-D', '--dryrun',      action='store_true',
                                             help="do not run RDataFrame (for debuggin)" )
  parser.add_argument('-j', '--nthreads',    default=None,
                                             help="number of threads, default: 8 locally, 1 for Dask" )
  parser.add_argument('-v', '--verbose',     dest='verbosity', type=int, nargs='?', const=1, default=0,
                                             help="set level of verbosity, default=%(default)s" )
  args = parser.parse_args()
  main(args)
  print(f">>> Done in {elapsed(TIME0)}")
  
