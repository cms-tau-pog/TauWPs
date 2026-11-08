#! /usr/bin/env python3
# Author: Izaak Neutelings (August 2025)
# Description:
#   The following was an attempt to analyze miniAOD with RDataFrame.
#   It works for the PAT Tau collection after jumping through some hoops.
#   However, I did not succeeed in getting it to work for the PAT Jet
#   or GenParticle collections, whose branches have a high splitting level...
#   Conclusion: It is easier to create custom nanoAOD without a cut
#   on the VSjet score to obtain an unbiased sample of tau candidates.
# Instructions
#   SCRAM_ARCH=el9_amd64_gcc12
#   cmsrel -n CMSSW_15_0_4_TauWPs CMSSW_15_0_4 # match to miniAOD version !
#   cd CMSSW_15_0_4_TauWPs/src/
#   cmsenv
#   python3 -m pip install -r requirements.txt # cmsstyle & dask_lxplus
#   cd $CMSSW_BASE/src/
#   git clone -b RDataFrame git@github.com:cms-tau-pog/TauWPs.git TauWPs
import ROOT; ROOT.PyConfig.IgnoreCommandLineOptions = True # to avoid conflict with argparse
from ROOT import RDF, RDataFrame, TFile


def getRDataFrame(infiles,tauAlgos,discriminants,nthreads=None,dask=False,verb=0):
  """
  Initiate RDataFrame, & add useful columns to RDataFrame. 
  """
  scores = ["raw"+a+d for a in tauAlgos for d in discriminants]
  
  print(">>> Loading libraries...")
  ROOT.gSystem.Load("libFWCoreFWLite")
  ROOT.gROOT.ProcessLine('FWLiteEnabler::enable();')
  ROOT.gSystem.Load("libDataFormatsPatCandidates")
  
  # DATA MODEL
  # The code below transforms RVeC of PAT objects to a RVeCs of primitive datatypes
  # that RDataFrame can more easily manipulate
  # Unfortunately this does not work as easily for "split" branches like Jets...
  print(">>> Declaring data model...")
  datamodel = (r"""
    #include "DataFormats/PatCandidates/interface/Tau.h"
    #include "DataFormats/PatCandidates/interface/PackedCandidate.h"
    #include <tuple>
    #include <vector>
    using RVecTau = ROOT::VecOps::RVec<pat::Tau>;
    using RVecJet = ROOT::VecOps::RVec<pat::Jet>;
    struct TauVecs {
      ROOT::RVecF pt;
      ROOT::RVecF eta;
      ROOT::RVecF dz;
  """ # add raw output scores
  + "\n  ".join(f"    ROOT::RVecF {s};" for s in scores) +
  r"""
    };
    TauVecs makeTauVecs(const RVecTau& taus) {
      TauVecs out;
      out.pt.reserve(taus.size());
      out.eta.reserve(taus.size());
  """ # reserve raw output scores
  + "\n  ".join(f"    out.{s}.reserve(taus.size());" for s in scores) +
  r"""
      for (const auto &t : taus) {
        out.pt.emplace_back(t.pt());
        out.eta.emplace_back(t.eta());
        out.dz.emplace_back(t.leadChargedHadrCand().isNonnull()?
          //dynamic_cast<const pat::PackedCandidate*>(&*t.leadChargedHadrCand())->dz():0);
          dynamic_cast<const pat::PackedCandidate*>(t.leadChargedHadrCand().get())->dz():0);
  """
  + "\n  ".join(f'      out.{s}.emplace_back(t.isTauIDAvailable("{s}")?t.tauID("{s}"):-1);' for s in scores) +
  r"""
      }
      return out;
    };
  """)
  print(datamodel)
  ROOT.gInterpreter.Declare(datamodel)
  
  print(f">>> Creating RDataFrame...")
  ROOT.EnableImplicitMT(8 if nthreads==None else nthreads) # set before creating RDataFrame
  dframe = RDataFrame('Events',infiles)
  if hasattr(RDF.Experimental,'AddProgressBar'):
    RDF.Experimental.AddProgressBar(dframe)
  
  # ADD MASKS for different generator components (genmatch)
  print(f">>> Defining columns...")
  dframe = (dframe
              .Define("taus",      "makeTauVecs(patTaus_slimmedTaus__PAT.obj)")
              .Define('mask_base', "taus.pt>20 & abs(taus.eta)<2.5 & abs(taus.dz)<0.2")\
              .Define("pt",        "taus.pt")
              .Define('pt_cut',    "pt[mask_base]")
  )
  
  # ADD MASKS for different tau ID algorithms
  for algo in tauAlgos:
    dframe = dframe.Define(f"mask_{algo}", f"mask_base & taus.raw{algo}VSjet>=0")
    for disc in discriminants:
      dframe = dframe.Define(f"raw{algo}{disc}", f"taus.raw{algo}{disc}[mask_{algo}]")
  
  return dframe
  

def main():
  # https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/interface/Tau.h
  # https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/src/Tau.cc
  # https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATTauHybridProducer.cc
  # https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/taus_cff.py
  infiles = [
    "ff14e667-2d71-4d34-b9b2-1cbc2ac700df.root",
    #"root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2810000/ff14e667-2d71-4d34-b9b2-1cbc2ac700df.root"
  ]
  tauAlgos = ['DeepTau2018v2p5','PNet','UParT']
  discriminants = ['VSjet','VSmu','VSe']
  dframe = getRDataFrame(infiles,tauAlgos,discriminants,nthreads=1,verb=2)
  
  print(f">>> Adding histograms...")
  results = [
    dframe.Histo1D(('pt','pt',100,0,200),'pt'),
    dframe.Histo1D(('pt_cut','pt_cut',100,0,200),'pt_cut')
  ]
  for algo in tauAlgos:
    for disc in discriminants:
      xvar = f"raw{algo}{disc}"
      hist = dframe.Histo1D((xvar,xvar,100,0,1),xvar)
      results.append(hist)
  
  print(f">>> Triggering event loop...")
  #RDF.RunGraphs(results)
  results[0].GetValue()
  
  # WRITE TO FILE
  print(f">>> Writing...")
  outfile = TFile('test.root','RECREATE')
  for hist in results:
    hist.GetValue().Write(hist.GetName())
  outfile.Close()
  

if __name__ == '__main__':
  main()
  print(">>> Done!")

