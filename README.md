# HH NanoAOD
Instructions for HH analysis using NanoAOD


## Quick start

```bash
cmsrel CMSSW_10_6_4
cd CMSSW_10_6_4/src
cmsenv

mkdir PhysicsTools/
git clone https://github.com/vivannguyen/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone https://github.com/vivannguyen/NUNanoAOD.git PhysicsTools/MonoZ

cd $CMSSW_BASE/src/PhysicsTools/NanoAODTools
git checkout remotes/origin/dev-nanoAODv7
cd $CMSSW_BASE/src
scram b -j 10

## Example: running over a signal MC file
To run locally, edit condor_Run2_proc.py so that the xsections yaml is from the data directory.
Uncomment:
   with open(os.path.dirname(__file__) +'../data/xsections_{}.yaml'.format(options.era)) as file:
Comment out: 
   with open(os.path.dirname(__file__) +'xsections_{}.yaml'.format(options.era)) as file:

cd $CMSSW_BASE/src/PhysicsTools/MonoZ/condor/
voms-proxy-init -voms cms --valid 168:00
python condor_Run2_proc.py --isMC=1 --doSyst=0 --era=2016 --nevt=1000 --infile=root://cms-xrd-global.cern.ch//store/mc/RunIISummer16NanoAODv7/GluGluToRadionToHHTo2B2ZTo2L2J_M-600_narrow_13TeV-madgraph-v2/NANOAODSIM/PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8-v1/60000/A760EFF8-EC29-5B48-B918-E38D12296512.root

## Example: submitting jobs to condor
To run on condor, edit condor_Run2_proc.py 
Comment out:
   with open(os.path.dirname(__file__) +'../data/xsections_{}.yaml'.format(options.era)) as file:
Uncomment out: 
   with open(os.path.dirname(__file__) +'xsections_{}.yaml'.format(options.era)) as file:
Can also turn off systematics by changing --doSyst default value

python run.py --isMC=1 --era=2016 --tag=[tagname] --input=../data/list_2016_MC.txt 
python run.py --isMC=0 --era=2016 --tag=[tagname] --input=../data/list_2016_Data_Nano02Apr2020.txt


```

## EFT reweighting
Code is based on framework from dihiggs group. 
```
https://gitlab.cern.ch/jlidrych/eftinterfaceforhhtobbzzllqq
```
### Electron scale factor
Tools provided by CMS EGamma POG
```
https://github.com/jlidrych/egm_tnp_analysis
```
### Muon scale factor
Tools provided by CMS Muon POG
```
https://gitlab.cern.ch/cms-muonPOG/spark_tnp
```

