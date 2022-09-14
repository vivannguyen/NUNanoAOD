import os, sys, re
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
# import PSet
import yaml
#Import the NanoAOD-tools that we will need
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme import jetRecalib
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSysColl import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.lepSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
#Import the MonoZ analysis tools
#from PhysicsTools.MonoZ.MonoZProducer import *
#from PhysicsTools.MonoZ.HHProducer import *
from PhysicsTools.MonoZ.ScaleFactorProducer import *
from PhysicsTools.MonoZ.GenWeightProducer import *
#from PhysicsTools.MonoZ.EWProducer import *
#from PhysicsTools.MonoZ.ADDProducer import *
#from PhysicsTools.MonoZ.NvtxPUreweight import *
from PhysicsTools.MonoZ.PhiXYCorrection import *
from PhysicsTools.MonoZ.BtagEventWeightProducer import *
from PhysicsTools.MonoZ.TriggerSFProducerForHH import *
from PhysicsTools.MonoZ.AngularVariablesProducerForHH import *
#from PhysicsTools.MonoZ.GenMonoZProducer import *
import argparse

parser = argparse.ArgumentParser("")
parser.add_argument('-isMC'   , '--isMC'   , type=int, default=1     , help="")
parser.add_argument('-jobNum' , '--jobNum' , type=int, default=1     , help="")
parser.add_argument('-era'    , '--era'    , type=str, default="2018", help="")
parser.add_argument('-doSyst' , '--doSyst' , type=int, default=0     , help="")
parser.add_argument('-infile' , '--infile' , type=str, default=None  , help="")
parser.add_argument('-dataset', '--dataset', type=str, default="X"   , help="")
parser.add_argument('-nevt'   , '--nevt'   , type=str, default=-1    , help="")
parser.add_argument('-json'   , '--json'   , type=str, default=None  , help="")
parser.add_argument('-period' , '--period' , type=str, default="Run2016B"  , help="")
options  = parser.parse_args()

def inputfile(nanofile):
   tested   = False
   forceaaa = False
   pfn=os.popen("edmFileUtil -d %s"%(nanofile)).read()
   pfn=re.sub("\n","",pfn)
   print nanofile," -> ",pfn

   if (os.getenv("GLIDECLIENT_Group","") != "overflow" and
       os.getenv("GLIDECLIENT_Group","") != "overflow_conservative" and not
       forceaaa ):
      if not tested:
         print "Testing file open"
         testfile=ROOT.TFile.Open(pfn)
         if testfile and testfile.IsOpen() :
            print "Test OK"
            nanofile=pfn
            testfile.Close()
         else:
            print 'pozor'
            if "root://cms-xrd-global.cern.ch/" not in nanofile:
               nanofile = "root://cms-xrd-global.cern.ch/" + nanofile
            forceaaa=True
      else:
         print 'nothing'
 #        nanofile = pfn
   else:
#      print 'nic'
       if "root://cms-xrd-global.cern.ch/" not in nanofile:
          nanofile = "root://cms-xrd-global.cern.ch/" + nanofile
   print 'final file is: ',nanofile
   return nanofile

options.infile = inputfile(options.infile)

print "---------------------------"
print " -- options  = ", options
print " -- is MC    = ", options.isMC
print " -- jobNum   = ", options.jobNum
print " -- era      = ", options.era
print " -- in file  = ", options.infile
print " -- dataset  = ", options.dataset
print "---------------------------"

xsection = 1.0
#nevents = 1
if options.isMC:
   condtag_ = "NANOAODSIM"
   if options.dataset == "X":
      options.dataset = options.infile
      options.dataset = options.dataset.split('/store')[1].split("/")
      condtag_ = options.dataset[5]
      options.dataset = options.dataset[3]
      runperiod_ = "None"
   print "[check] condtag_ == ", condtag_
   print "[check] dataset  == ", options.dataset
else:
   if options.dataset == "X":
      options.dataset = options.infile
      options.dataset = options.dataset.split('/store')[1].split("/")
      condtag_ = options.dataset[2]
      options.dataset = options.dataset[3]
      runperiod_ = condtag_.split(options.era)[1][:1]
   else:
      options.dataset = options.dataset.split("/")
      condtag_ = options.dataset[2]
      options.dataset = options.dataset[1]
      runperiod_ = condtag_.split(options.era)[1][:1]
   print "[check] condtag_ == ", condtag_
   print "[check] dataset  == ", options.dataset

dataset = options.dataset 
#condtag_ = "SUEP_QCD"

if options.isMC:
#   with open(os.path.dirname(__file__) +'../data/xsections_{}.yaml'.format(options.era)) as file:
   with open(os.path.dirname(__file__) +'xsections_{}.yaml'.format(options.era)) as file:
       #MC_xsecs = yaml.full_load(file)
       MC_xsecs = yaml.safe_load(file)
   try:
       xsection *= MC_xsecs[dataset]["xsec"]
       xsection *= MC_xsecs[dataset]["kr"]
       xsection *= MC_xsecs[dataset]["br"]
       print("The dataset name: ", dataset)
       print("The xsection: ", xsection)
   except:
       print("WARNING: I did not find the xsection for that MC sample. Check the dataset name and the relevant yaml file")

pre_selection = "(Sum$(Electron_pt>8 && abs(Electron_eta)<2.5) >= 2 || Sum$(Muon_pt>8 && abs(Muon_eta)<2.5) >= 2 ) "
pre_selection += "&& Flag_METFilters"

if float(options.nevt) > 0:
   print " passing this cut and : ", options.nevt
   pre_selection += ' && (Entry$ < {})'.format(options.nevt)

#modules_era   = [
#    GenWeightProducer(
#       isMC = options.isMC,
#       dopdf = False if ("ADD" in options.dataset or "Unpart" in options.dataset) else True
#    )
#]

if "QCD" in dataset:
    modules_era = [ GenWeightProducer( isMC = options.isMC, xsec = xsection, dopdf =  False, do_xsecscale = True) ]
else:
    modules_era = [ GenWeightProducer( isMC = options.isMC, xsec = xsection, dopdf =  True, do_xsecscale = True) ]

# modules_era = [ GenWeightProducer( isMC = options.isMC, xsec = xsection, dopdf =  True) ]

pro_syst = [ "ElectronEn", "MuonEn", "jesTotal", "jer"]
ext_syst = [ "puWeight", "PDF", "MuonSF", "ElecronSF", "EWK", "nvtxWeight","TriggerSFWeight","btagEventWeight", "QCDScale0w", "QCDScale1w", "QCDScale2w"]

jetmetCorrector = createJMECorrector(isMC = options.isMC, dataYear = options.era, runPeriod= runperiod_ , jesUncert="Total") # add METfixEE2017

if options.isMC:
   print "sample : ", options.dataset, " candtag : ", condtag_
   options.period = condtag_
   try:
      combineHLT = yaml.safe_load(open("combineHLT_Run2_MET.yaml"))
   except yaml.YAMLError as exc:
      print(exc)
   if options.era=="2016":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2016All.MC", "") + ")"
   if options.era=="2017":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2017All.MC", "") + ")"
   if options.era=="2018":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2018All.MC", "") + ")"

   if options.era=="2016":
      modules_era.append(jetmetCorrector())
      modules_era.append(puAutoWeight_2016())
      modules_era.append(PrefCorr())
      modules_era.append(muonScaleRes2016())
      modules_era.append(lepSF_2016())
      ext_syst.append("PrefireWeight")
   if options.era=="2017":
      modules_era.append(jetmetCorrector())
      modules_era.append(puAutoWeight_2017())
      modules_era.append(PrefCorr())
      modules_era.append(muonScaleRes2017())
      modules_era.append(lepSF_2017())
      ext_syst.append("PrefireWeight")
   if options.era=="2018":
      modules_era.append(jetmetCorrector())
      modules_era.append(puAutoWeight_2018())
      modules_era.append(muonScaleRes2018())
      modules_era.append(lepSF_2018())


   modules_era.append(ScaleFactorProducer(isMC=options.isMC, era=str(options.era), period=str(options.period), do_syst=0, syst_var=''))

   if options.era=="2016":
      modules_era.append(TriggerSF_2016())
   if options.era=="2017":
      modules_era.append(TriggerSF_2017())
   if options.era=="2018":
      modules_era.append(TriggerSF_2018())

   # for shift-based systematics
   if options.doSyst:
      for sys in pro_syst:
         for var in ["Up", "Down"]:
	    # if "jesTotal" in sys and options.doSyst==1: modules_era.append(PhiXYCorrection(era=options.era,isMC=options.isMC,sys=sys+var))
	    # if "jer" in sys and options.doSyst==1: modules_era.append(PhiXYCorrection(era=options.era,isMC=options.isMC,sys=sys+var))
            modules_era.append(ScaleFactorProducer(options.isMC, str(options.era), period=str(options.period),do_syst=options.doSyst, syst_var=sys+var))

else:
   if options.era=="2016":
      modules_era.append(muonScaleRes2016())
   if options.era=="2017":
      modules_era.append(muonScaleRes2017())
   if options.era=="2018":
      modules_era.append(muonScaleRes2018())

   options.period = condtag_
   print "sample : ", options.dataset, " candtag : ", condtag_
   try:
      combineHLT = yaml.safe_load(open("combineHLT_Run2_MET.yaml"))
   except yaml.YAMLError as exc:
      print(exc)
   if options.era=="2016":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2016All.%s" % options.dataset, "") + ")"
   if options.era=="2017":
      if 'Run2017B' in condtag_:
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2017B.%s" % options.dataset, "") + ")"
      else:
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2017All.%s" % options.dataset, "") + ")"
   if options.era=="2018":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2018All.%s" % options.dataset, "") + ")"

   print " -- era : ",
   modules_era.append(jetmetCorrector())   
   modules_era.append(ScaleFactorProducer(isMC=options.isMC, era=str(options.era), period=str(options.period),do_syst=0, syst_var=''))

   if options.era=="2016":
       options.json = "Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt"
   if options.era=="2017":
       options.json = "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"
   if options.era=="2018":
       options.json = "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"
   print "---- JSON used is : ", options.json


for i in modules_era:
   print "modules : ", i

print "Selection : ", pre_selection

p = PostProcessor(
   ".", [options.infile],
   cut=pre_selection,
   branchsel="keep_and_drop.txt",
   outputbranchsel="keep_and_drop_post.txt",
   haddFileName="tree_%s.root" % str(options.jobNum),
   modules=modules_era,
   provenance=True,
   noOut=False,
   fwkJobReport=True,
   jsonInput=options.json
)
p.run()
