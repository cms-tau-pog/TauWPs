# TauWPs

The tools in this repository compute & study efficiencies of tau identification algorithms in CMS.
The scripts rely on [ROOT's RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html)
for processing nanoAOD files, and [`cmsstyle`](https://cmsstyle.readthedocs.io) for plotting.


## Installation

<!--### For nanoAOD: LCG's ROOT-->

For nanoAOD files, we do not need CMSSW, but can use the more lightweight
**[LCG environment](https://lcgdocs.web.cern.ch/lcgdocs/lcgreleases/introduction/)**.
Please use ROOT v6.32 or some useful RDataFrame features, e.g.
```bash
. /cvmfs/sft.cern.ch/lcg/views/LCG_106c/x86_64-el9-gcc13-opt/setup.sh
```
Install CMS Style and Dask for lxplus and the `rdataframe` branch of this repository:
```bash
python3 -m pip install -r requirements.txt # cmsstyle & dask_lxplus
git clone -b rdataframe git@github.com:cms-tau-pog/TauWPs.git TauWPs
```
With **each new shell session**, set up Kerberos and GRID authentication (if needed),
and add your local python path:
```bash
kinit && aklog # Kerberos authentication for HTCondor
voms-proxy-init -voms cms -valid 200:0 # GRID authentication
export PYTHONPATH=~/.local/lib/python3.11/site-packages:$PYTHONPATH
cd TauWPs/
```

<!--### For miniAOD: CMSSW

> [!WARNING]
> This part does not work so well. Better to use NanoAOD.

To process miniAOD files, you need a matching CMSSW version to access the objects correctly.
Look up the fragment of your sample on [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/),
or on PmdV's pages for [Run 2 UL](https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis)
and [Run 3](https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun3Analysis).

For 2024, MC was produced with `CMSSW_15_0_4`:
```bash
SCRAM_ARCH=el9_amd64_gcc12
cmsrel -n CMSSW_15_0_4_TauWPs CMSSW_15_0_4 # match to miniAOD version !
cd CMSSW_15_0_4_TauWPs/src/
cmsenv
```
Install CMS Style and Dask for lxplus and this repository:
```bash
python3 -m pip install -r requirements.txt # cmsstyle & dask_lxplus
cd $CMSSW_BASE/src/
git clone -b RDataFrame git@github.com:cms-tau-pog/TauWPs.git TauWPs
```
With each new shell session, set up Kerberos and GRID authentication (if needed),
setup CMSSW, and add your local python path:
```bash
kinit && aklog # Kerberos authentication for HTCondor
voms-proxy-init -voms cms -valid 200:0 # GRID authentication
cd CMSSW_15_0_4_TauWPs/src/
cmsenv
export PYTHONPATH=~/.local/lib/python3.11/site-packages:$PYTHONPATH
cd $CMSSW_BASE/src/TauWPs
```-->


## Creating histograms from nanoAOD

There are several scripts: 
1. `makeTauHists.py`: Create histograms from nanoAOD.
1. `plotTauROCCurves.py`: Plot the histograms created by `makeTauHists.py`, compute efficiencies & misidentification probability, plot efficiency curves, ROC curves, etc.
1. `plotTauPt.py`: Plot the pT response (ptreco/ptgen).
1. `plotTauIDVennDiagrams.py`: Plot Venn diagrams and tables to compare overlap of the different reco/IDs.

### Process nanoAOD locally
Run RDataFrame multithreaded on one or more nanoAOD files as follows
```bash
python3 makeTauHists.py --nthreads 8 -t _2024 /store/mc/RunIII2024Summer24NanoAODv15/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/100000/1059354e-54e8-4beb-bb83-742e23d89675.root
```
Passing a tag via `-t` is optional, but helps organize the output files.

You can also organize larger file lists into ASCII files
(see [`samples/samples*txt`](blob/master/samples/)).
```bash
python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --nthreads 8 -t _2024_HTT -f HTo2Tau
python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --nthreads 8 -t _2024_DY -f DY
python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --nthreads 8 -t _2024_TT -f TTto
```
Filter ROOT files with the `-f`/`--filter` flags (e.g. if you want to use Higgs samples as the source of real taus, but ttbar for jet fakes):
```bash
python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --dask -t _2024_HTT -f HTo2Tau
python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --dask -t _2024_DY -f DY
python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --dask -t _2024_TT -f TTto
```

### Process nanoAOD on HTCondor via Dask
Running RDataFrame locally is quite fast. The RDataFrame defined in `makeTauHists.py` can process
roughly 500k events per minute using 16 threads and nanoAOD samples on EOS.
However, for a huge amount of files, especially if they stored on the GRID, one can submit to HTCondor via Dask with the `--dask` flag.
In case HTCondor takes long to start and return jobs,
it might be a good idea to create a persistent `tmux` or `screen` session in which you start `makeTauHists.py`.
Please find instructions for CERN's lxplus on [this](https://cern.service-now.com/service-portal?id=kb_article&n=KB0008111)
and [this page](https://hsf-training.github.io/analysis-essentials/shell-extras/persistent-screen.html).
```bash
systemctl --user start tmux.service # to have tmux screen persist on lxplus9
tmux # use tmux screen or whatever to wait for jobs to finish (write down lxplus node)
python3 makeTauHists.py samples/samples_nanos_2024_passthrough.txt --dask -t _2024
# wait, or deattach the tmux screen with Ctrl+b and then d to check condor_q
```

> [!WARNING]
> Currently, only lxplus's HTCondor is supported, but one can easily configure other clusters for Dask.

### Histogram outputs

The output histograms from `makeTauHists.py` are organized into the following directories in a single ROOT file:
1. `ntau/`: Control plots of total number of generated or reconstructed taus passing some selection (for example, `ntau_all_real`).
1. `pt/` For pT spectra, pT response (ptreco/ptgen), and pT vs. raw scores (for pT-binned ROC curves).
1. `dm/` For DM vs. raw scorces (for DM-binned ROC curves).
1. `raw/`: Raw discriminator scores, with one subdirectory per reco/ID collection.
1. `charge/`: Charge-assignment probability (`qConf`) histograms, with one subdirectory per reco/ID collection.
1. `venn/`: 3D histograms comparing a set of three tau reco/ID collections (e.g. HPS+DeepTau vs. CHS+PNet vs. PUPPI+UParT).

Histogram names put the variable and algorithm first, followed by the selection and transformation postfixes and, where applicable, the genmatch.
```bash
${dir}/${xvar}${algorithm}_${postfix}_${genmatch}
```
For example, `raw/PNetCHS/rawPNetCHSVSjet_logit_qlty_real` is the logit-transformed VSjet score for genuine taus passing the `qlty` selection. A count histogram is named `ntau/ntau_all_real`. Some multidimensional histograms have custom names; for example, Venn histograms use `venn/ntau_venn[_WP]_{genmatch}`.

Selection and transformation postfixes indicate:
1. `all`: Baseline tau selection with algorithm-specific pT cuts, e.g. `abs(Tau_eta)<2.5 & abs(Tau_dz)<0.2`.
1. `all2`: Baseline tau selection with a common pT cut, e.g. `Tau_pt>20 & abs(Tau_eta)<2.5 & abs(Tau_dz)<0.2`, where `Tau_pt` is HPS-reconstructed pT if the tau candidate has an HPS match, otherwise the raw, uncorrected pT of the CHS or PUPPI jet.
1. `qlty`: Require that the target tau ID exists (>-0.5), plus quality cuts (excluding 2-prong decay modes for HPS+DeepTau, and requiring |qConf| > 0.2 for PNet or UParT).
1. `hps`: Require the tau candidate is a HPS object (including reconstructed 2-prong decay modes).
1. `hps2`: Same as `hps`, except exclude 2-prong decay modes.
1. `hps3`: All HPS-matched tau candidates excluding 2-prong DMs, _plus_ "recovery" tau candidates (i.e., CHS/PUPPI without a HPS match) that pass |qConf| > 0.2. This is used as standard for the comparison in the DPS note.
1. `logit`: Logit transformation, logit(x) = 2*atanh(2*x-1), to "zoom" more into the densely distributed raw score around 0 and 1, plus finer binning. This makes computing efficiencies more accurate.

> [!IMPORTANT]
> To compare HPS+DeepTau tau candidates to tau candidates reconstructed/recovered with AK4 CHS or PNet,
> one can define the following "common denominator" that is the "union" of "HPS+CHS+PUPPI" candidates.
> For `all2`, a common pT definition is used:
> ```cpp
> abs(Tau_eta)<2.5 & abs(Tau_dz)<0.2 & Tau_pt>20.0
> ```
> Note that `Tau_pt` is the pT of HPS tau if it is reconstructed as an HPS tau, otherwise, it is the value of the AK4 jet pT.
> For `all`, the pT cut depends on the algorithm (HPS pT for HPS taus, regressed pT for PNet or UParT):
> ```cpp
> abs(Tau_eta)<2.5 & abs(Tau_dz)<0.2 & (
>   (Tau_pt>20.0 & Tau_rawDeepTau2018v2p5VSjet>-0.5)
>   | (Tau_pt*Tau_ptCorrPNetCHS>20.0 & Tau_rawPNetCHSVSjet>-0.5)
>   | (Tau_pt*Tau_ptCorrPNetPUPPI>20.0 & Tau_rawPNetPUPPIVSjet>-0.5)
>   | (Tau_pt*Tau_ptCorrUParTPUPPI>20.0 & Tau_rawUParTPUPPIVSjet>-0.5)
> )
> ```
> where 2-prong HPS candidates are exlcuded for HPS+DeepTau,
> and `|qConf|>0.2` require a high charge quality for PNet and UParT.

### YAML output

Unfortunately, the code that sets up the RDataFrame and many parallel definitions of selections (as "masks")
can get a bit confusing and opaque, so the `makeTauHists.py` script also creates a YAML file
that contains a dictionary with the definition of all variables and masks (=selection strings).


## Analyzing results

The ouput from `makeTauHists.py` can be used as inputs to the `plot*.py` scripts,
which in turn create plots in `plots/`

### Compute the WPs

First compute the WPs, which will be used by `plotTauROCCurves.py`:
```bash
python3 computeTauIDWPs.py output/tauAlgoEffs_2024.root -t _2024 -v1
```

### Compute & plot efficiencies

To create the basic ROC curves, run
```bash
python3 plotTauROCCurves.py output/tauAlgoEffs_2024.root --era 2024 -t _2024
```

You can run different "methods" with `plotTauROCCurves.py`,
by passing the right keyword to `-m`/`--method`:

- `raw`:   Plot raw output scores, comparing different tau components (real, muon, electron, jet).
- `qConf`: Plot charge assignment probability (for PNet and UParT), comparing different tau components.
- `roc`:   Plot ROC curves, starting from raw output scores.
- `eff`:   Plot efficiency curves (incl. WPs).
- `all` (default): Run all the plotting methods above.

Compute working-point cuts separately with `computeTauIDWPs.py`; `wp` is not a `plotTauROCCurves.py` method. The script writes its text and YAML tables to `plots/wps` by default. Use `-o output/wps` to write the YAML files to the directory used by `makeTauHists.py` and WP-sliced ROC plots. The WP tag in the generated filenames must match the tag expected by the consumer.

If you have produced separate histograms for different samples (e.g. split by physics process),
you can run
```bash
python3 computeTauIDWPs.py output/tauAlgoEffs_2024_HTT.root fake=output/tauAlgoEffs_2024_DY.root -t _HTT-DY -v1
python3 plotTauROCCurves.py output/tauAlgoEffs_2024_HTT.root fake=output/tauAlgoEffs_2024_DY.root -t _HTT-DY -m roc -v1
```
where the genuine taus will be taken from `output/tauAlgoEffs_2024_HTT.root`,
and the fake ones (electrons, muons, or jets) from `output/tauAlgoEffs_2024_DY.root`.


### Plot pt response

Use `plotTauPt.py` to plot pT response.
```bash
python3 plotTauPt.py output/tauAlgoEffs_2024.root --era 2024 -t _2024
```

### Plot Venn diagrams

Use `plotTauIDVennDiagrams.py` to compute the Venn diagram of HPS+DeepTau, CHS+PNet, PUPPI+UParT, to study where the overlap and which reco/ID recovers what.
```bash
python3 plotTauIDVennDiagrams.py output/tauAlgoEffs_2024.root --era 2024 -t _2024
```

## :bulb: Pro-tips
Here are some quick tips 'n tricks to understand the data you are looking at.


### DAS
In a clean environment, you query DAS from the command line.
```bash
# Look for datasets with pattern
dasgoclient --query="/GluGlu*HTo2Tau*/*/NANO*"
dasgoclient --query="/TTtoLNu2Q*/*/NANO*"

# Get file lists (use full dataset path)
dasgoclient --query="/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM file"
```

### NanoAOD
Check the contents in the ROOT command line as
```bash
root -l root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/2520000/004eb958-f042-48a7-b7a3-2f1cd5185215.root

Events->Print("*Tau*")
```

Quickly draw the raw output score for real tau candidates, e.g., in H &rarr; &tau;&tau; samples:
```bash
root -l root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/100000/1059354e-54e8-4beb-bb83-742e23d89675.root

Events->Draw("Tau_rawDeepTau2018v2p5VSjet >> h(100,0,1)","Tau_genPartFlav==5")
Events->Draw("Tau_rawPNetVSjet >> h(100,0,1)","Tau_genPartFlav==5")
Events->Draw("Tau_rawUParTVSjet >> h(100,0,1)","Tau_genPartFlav==5")
```

Quickly draw the raw output score for jet fakes, e.g., in ttbar samples:
```bash
root -l root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/2520000/004eb958-f042-48a7-b7a3-2f1cd5185215.root

Events->Draw("Tau_rawDeepTau2018v2p5VSjet >> h(100,0,1)","Tau_genPartFlav==0")
Events->Draw("Tau_rawPNetVSjet >> h(100,0,1)","Tau_genPartFlav==0")
Events->Draw("Tau_rawUParTVSjet >> h(100,0,1)","Tau_genPartFlav==0")
```

### RDataFrame
RDataFrame is also easy to use in C++ and/or the ROOT command line:
```cpp
from ROOT import RDataFrame
df = RDataFrame("Events","root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/100000/1059354e-54e8-4beb-bb83-742e23d89675.root")
df = df.Define('pt', "Tau_pt[abs(Tau_eta)<2.5 & Tau_genPartFlav==5"])
h = df.Histo1D(('pt','pt',100,0,100),'pt')
h->Draw() // event loop is triggered here
```
Or in python:
```python
from ROOT import RDataFrame
df = RDataFrame('Events',"root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24NanoAODv15/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v2-v2/100000/1059354e-54e8-4beb-bb83-742e23d89675.root")
df = df.Define('pt', "Tau_pt[abs(Tau_eta)<2.5 & Tau_genPartFlav==5]") # the cuts are applied as a "mask", a vector of booleans
h = df.Histo1D(('pt','pt',100,0,100),'pt')
h.GetValue().Draw() # event loop is triggered here
```
You can count the number of all taus that pass some requirement as follows:
```python
df = df.Define('mask_den', "Tau_pt>20 & abs(Tau_eta)<2.5 & Tau_genPartFlav==5 & Tau_rawDeepTau2018v2p5VSjet>=0")
df = df.Define('mask_num', "mask_den & Tau_rawDeepTau2018v2p5VSjet>0.5")
den = df.Sum('mask_den') # denominator, counting all True values of all tau objects in all events 
num = df.Sum('mask_num') # numerator
ratio = num.GetValue() / den.GetValue() # event loop is triggered here
```
You can match reconstructed tau objects with generated tau objects using `Take`
(ROOT v32 or newer to use a default for if the index is -1):
```python
df = df.Define('mask',   "Tau_pt>20 & Tau_genPartFlav==5") # selection on reco tau
df = df.Define('recodm', "Tau_decayMode[mask]")
df = df.Define('gendm',  "Take<int>(GenVisTau_status,Tau_genPartIdx[mask],-1)") # defauly -1 if Tau_genPartIdx is out of GenVisTau_status's range
h1 = df.Histo1D(('gendm','Generator DM',20,-1,19),'gendm')
h2 = df.Histo2D(('reco_vs_gendm','Reconstructed vs. Generator DM',20,-1,19,20,-1,19),'gendm','recodm')
h2.Draw('COLZ') # event loop is triggered here for both h1 and h2
h1.Scale(16./h1.GetMaximum()); h1.SetLineColor(2); h1.SetLineWidth(2) # thick red
h1.Draw('SAME HIST')
```
Note that if the files are available locally (e.g. on EOS), you can also use a wildcard:
```python
indir = "/eos/home-b/bluj/data/NanoAOD/DYto2Tau_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8/DYto2Tau_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15TauRelaxed-150X_mcRun3_2024_realistic_v2-v1/250709_161740/0000"
df = RDataFrame('Events',f"{indir}/*.root")
```

### MiniAOD

Quickly plot tau branches in the ROOT command line with
```bash
root -l root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2810000/ff14e667-2d71-4d34-b9b2-1cbc2ac700df.root

Events->Draw("patTaus_slimmedTaus__PAT.obj.pt() >> h(100,0,100)")
```

Use the `edmDumpEventContent` command in a CMSSW environment to check the contents of a miniAOD file.
```bash
# Dump content of Events tree to find desired (tau) object
edmDumpEventContent root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2810000/ff14e667-2d71-4d34-b9b2-1cbc2ac700df.root | grep -iC2 tau

# Dump full path of Tau PAT object, e.g. patTaus_slimmedTaus__PAT
edmDumpEventContent --all --regex '[Tt][Aa][Uu]' root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2810000/ff14e667-2d71-4d34-b9b2-1cbc2ac700df.root
```
Use the `edmProvDump` command to dump the "provenance" of tracked configuration parameters,
which is useful to know with which settings the sample was produced:
```bash
# Dump configuration of tau-related configuration
edmProvDump root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2810000/ff14e667-2d71-4d34-b9b2-1cbc2ac700df.root | grep -iC2 'tau\|VSjet\|chargeAssignmentProbMin'

# Dump configuration of just slimmedTau collections
edmProvDump --findMatch 'slimmedTaus' root://xrootd-cms.infn.it//store/mc/RunIII2024Summer24MiniAODv6/GluGluH-HTo2Tau_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/2810000/ff14e667-2d71-4d34-b9b2-1cbc2ac700df.root
```

These should also work for nanoAOD files, depending on how they were processed.
More information can be found on [this page](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookEdmUtilities)
and [this page](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookEdmInfoOnDataFile).

> [!IMPORTANT]  
> For _unbiased_ computation of efficiencies, you need to make sure the input nanoAOD samples
> have no "hidden" cuts on the "slimmed" tau collection, which default central sample have.
> For this exercise, you probably need custom nanoAOD with these cuts removed.
> With `edmProvDump` you can check that there are no lower cuts on VSjet.
> For example, central nanoAODs for 2024 were produced with the following configuration for hybrid taus:
> ```python
> Module: slimmedTaus PAT
>  parameters: {
>   @module_edm_type: string tracked  = 'EDProducer'
>   @module_label: string tracked  = 'slimmedTaus'
>   @module_type: string tracked  = 'PATTauHybridProducer'
>   UTagLabel: string tracked  = 'pfUnifiedParticleTransformerAK4JetTags'
>   UTagScoreNames: vstring tracked  = {'pfUnifiedParticleTransformerAK4JetTags:probb','pfUnifiedParticleTransformerAK4JetTags:probbb','pfUnifiedParticleTransformerAK4JetTags:problepb','pfUnifiedParticleTransformerAK4JetTags:probc','pfUnifiedParticleTransformerAK4JetTags:probs','pfUnifiedParticleTransformerAK4JetTags:probu','pfUnifiedParticleTransformerAK4JetTags:probd','pfUnifiedParticleTransformerAK4JetTags:probg','pfUnifiedParticleTransformerAK4JetTags:probele','pfUnifiedParticleTransformerAK4JetTags:probmu','pfUnifiedParticleTransformerAK4JetTags:probtaup1h0p','pfUnifiedParticleTransformerAK4JetTags:probtaup1h1p','pfUnifiedParticleTransformerAK4JetTags:probtaup1h2p','pfUnifiedParticleTransformerAK4JetTags:probtaup3h0p','pfUnifiedParticleTransformerAK4JetTags:probtaup3h1p','pfUnifiedParticleTransformerAK4JetTags:probtaum1h0p','pfUnifiedParticleTransformerAK4JetTags:probtaum1h1p','pfUnifiedParticleTransformerAK4JetTags:probtaum1h2p','pfUnifiedParticleTransformerAK4JetTags:probtaum3h0p','pfUnifiedParticleTransformerAK4JetTags:probtaum3h1p','pfUnifiedParticleTransformerAK4JetTags:ptcorr','pfUnifiedParticleTransformerAK4JetTags:ptreshigh','pfUnifiedParticleTransformerAK4JetTags:ptreslow','pfUnifiedParticleTransformerAK4JetTags:ptnu','pfUnifiedParticleTransformerAK4JetTags:probemudata','pfUnifiedParticleTransformerAK4JetTags:probemumc','pfUnifiedParticleTransformerAK4JetTags:probdimudata','pfUnifiedParticleTransformerAK4JetTags:probdimumc','pfUnifiedParticleTransformerAK4JetTags:probmutaudata','pfUnifiedParticleTransformerAK4JetTags:probmutaumc'}
>   UtagPtCorrName: string tracked  = 'ptcorr'
>   addGenJetMatch: bool tracked  = true
>   chargeAssignmentProbMin: double tracked  = 0.20000000000000001 # this choice affect the effiency computation
>   checkTauScoreIsBest: bool tracked  = false
>   dRMax: double tracked  = 0.40000000000000002
>   genJetMatch: InputTag tracked  = tauGenJetMatchPUPPIJet
>   jetEtaMax: double tracked  = 2.5
>   jetPtMin: double tracked  = 15
>   jetSource: InputTag tracked  = slimmedJetsPuppi
>   src: InputTag tracked  = slimmedTausWithUTagCHS
>   tagPrefix: string tracked  = 'byUTagPUPPI'
>   tauScoreMin: double tracked  = -1
>   usePFLeptonsAsChargedHadrons: bool tracked  = true
>   vsJetMin: double tracked  = 0.050000000000000003   # this can bias effiency computation
>  }
> ```


### Source code

You can dig into the CMSSW source code on GitHub to know how the miniAOD or nanoAOD were produced. Note you need to select the correct branch corresponding to the CMSSW release that was used for production (e.g. [`CMSSW_15_0_X`](https://github.com/cms-sw/cmssw/blob/CMSSW_15_0_X/DataFormats/PatCandidates/src/Tau.cc)).
- [`DataFormats/PatCandidates/src/Tau.cc`](https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/src/Tau.cc) and [`DataFormats/PatCandidates/interface/Tau.h`](https://github.com/cms-sw/cmssw/blob/master/DataFormats/PatCandidates/interface/Tau.h) define the PAT Tau objects in miniAOD.
- [`PhysicsTools/PatAlgos/plugins/PATTauHybridProducer.cc`](https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATTauHybridProducer.cc) is the "HybridTau" producer that matches to HPS taus and "recovers" tau candidates from AK4 CHS (+PNet) and AK4 PUPPI (+PNet or UParT).
- [`PhysicsTools/NanoAOD/python/taus_cff.py`](https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/taus_cff.py) are used to configure the Tau collection in nanoAOD.

The configuration for these producers for a given miniAOD or nanoAOD file can be extracted with `edmProvDump`, see the previous section.


## References

- Michał Bluj (4/2025) https://indico.cern.ch/event/1535282/#2-unified-tau-taggers-performa
- [CMS-DP-2025-081](https://cds.cern.ch/record/2948917) [[plots](https://btv-public.docs.cern.ch/DPNotesRun3/DP-2025-081/DP-2025-081/)]: "Flavour tagging performance of the updated Unified Particle Transformer algorithm with the CMS experiment at sqrt(s)= 13.6 TeV"
- [CMS-DP-2025-073](https://cds.cern.ch/record/2946445) [[plots](https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsDP25073)]: "Comparison of the performance of tau reconstruction and identification algorithms in Run 3"
- [CMS-DP-2025-074](https://cds.cern.ch/record/2946446) [[plots](https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsDP25074)]: "Validation of tau identification algorithms in Run 3 data and simulation"
