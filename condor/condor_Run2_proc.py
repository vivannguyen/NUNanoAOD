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
from PhysicsTools.MonoZ.DiHiggsProducer import *
from PhysicsTools.MonoZ.GenWeightProducer import *
from PhysicsTools.MonoZ.GenTopProducer import *
#from PhysicsTools.MonoZ.EWProducer import *
#from PhysicsTools.MonoZ.ADDProducer import *
#from PhysicsTools.MonoZ.NvtxPUreweight import *
from PhysicsTools.MonoZ.PhiXYCorrection import *
from PhysicsTools.MonoZ.BtagEventWeightProducer import *
from PhysicsTools.MonoZ.TriggerSFProducerForHH import *
from PhysicsTools.MonoZ.AngularVariablesProducerForHH import *
#from PhysicsTools.MonoZ.BDTdiscriminantProducerForHH import *
from PhysicsTools.MonoZ.StoreEventsProducerForHH import *
#from PhysicsTools.MonoZ.GenMonoZProducer import *
import argparse

parser = argparse.ArgumentParser("")
parser.add_argument('-isMC'   , '--isMC'   , type=int, default=1     , help="")
parser.add_argument('-jobNum' , '--jobNum' , type=int, default=1     , help="")
parser.add_argument('-era'    , '--era'    , type=str, default="2018", help="")
parser.add_argument('-doSyst' , '--doSyst' , type=int, default=1     , help="")
parser.add_argument('-signal' , '--signal' , type=str, default='GF'  , help="")
parser.add_argument('-infile' , '--infile' , type=str, default=None  , help="")
parser.add_argument('-dataset', '--dataset', type=str, default="X"   , help="")
parser.add_argument('-nevt'   , '--nevt'   , type=str, default=-1    , help="")
parser.add_argument('-json'   , '--json'   , type=str, default=None  , help="")

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
            if "root://cms-xrd-global.cern.ch/" not in nanofile:
               nanofile = "root://cms-xrd-global.cern.ch/" + nanofile
            forceaaa=True
      else:
         nanofile = pfn
   else:
#      print 'nic'
       if "root://cms-xrd-global.cern.ch/" not in nanofile:
          nanofile = "root://cms-xrd-global.cern.ch/" + nanofile
   return nanofile

options.infile = inputfile(options.infile)

# Use EE noise mitigation for 2017
metBranchName = "METFixEE2017" if options.era=="2017" else "MET"

print "---------------------------"
print " -- options  = ", options
print " -- is MC    = ", options.isMC
print " -- jobNum   = ", options.jobNum
print " -- era      = ", options.era
print " -- in file  = ", options.infile
print " -- dataset  = ", options.dataset
print " -- MET      = ", metBranchName
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

# for VBF
if options.signal=='VBF':
    pre_selection = "(Sum$(Electron_pt>10 && abs(Electron_eta)<2.5) >= 2 || Sum$(Muon_pt>8 && abs(Muon_eta)<2.5) >= 2 ) && Sum$(Jet_pt>15)>= 4 "
# for EMU control region
elif options.signal=='EMU':
    pre_selection = "(Sum$(Electron_pt>10 && abs(Electron_eta)<2.5) >= 1 || Sum$(Muon_pt>8 && abs(Muon_eta)<2.5) >= 1 ) && Sum$(Jet_pt>15)>= 4 "
# for GF
else:
    pre_selection = "(Sum$(Electron_pt>10 && abs(Electron_eta)<2.5) >= 2 || Sum$(Muon_pt>8 && abs(Muon_eta)<2.5) >= 2 ) && Sum$(Jet_pt>15 && abs(Jet_eta) < 2.4 )>= 4 "
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

pro_syst = [ "ElectronEn", "MuonEn", "jer", "UnclusteredEn", "bRegScale", "bRegSmear"]
ext_syst = [ "puWeight", "PDF", "MuonSF", "ElecronSF","TriggerSFWeight","btagEventWeight", "QCDScale0w", "QCDScale1w", "QCDScale2w"]

jetmetCorrector = createJMECorrector(isMC = options.isMC, dataYear = options.era, runPeriod= runperiod_ , jesUncert="Merged", metBranchName=metBranchName) # add METfixEE2017

if options.isMC:
   if options.era=="2016":
      for ix in ['jesAbsolute', 'jesAbsolute_2016', 'jesBBEC1', 'jesBBEC1_2016', 'jesEC2', 'jesEC2_2016', 'jesFlavorQCD', 'jesHF', 'jesHF_2016', 'jesRelativeBal', 'jesRelativeSample_2016']:
         pro_syst.append(ix)
   if options.era=="2017":
      for ix in ['jesAbsolute', 'jesAbsolute_2017', 'jesBBEC1', 'jesBBEC1_2017', 'jesEC2', 'jesEC2_2017', 'jesFlavorQCD', 'jesHF', 'jesHF_2017', 'jesRelativeBal', 'jesRelativeSample_2017']:
         pro_syst.append(ix)
   if options.era=="2018":
      for ix in ['jesAbsolute', 'jesAbsolute_2018', 'jesBBEC1', 'jesBBEC1_2018', 'jesEC2', 'jesEC2_2018', 'jesFlavorQCD', 'jesHF', 'jesHF_2018', 'jesRelativeBal', 'jesRelativeSample_2018','jesHEMIssue']:
         pro_syst.append(ix)

if options.isMC:
   print "sample : ", options.dataset, " candtag : ", condtag_
   try:
      combineHLT = yaml.safe_load(open("combineHLT_Run2.yaml"))
   except yaml.YAMLError as exc:
      print(exc)

   if options.era=="2016":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2016.MC", "") + ")"
   if options.era=="2017":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2017.MC", "") + ")"
   if options.era=="2018":
      pre_selection = pre_selection + " && (" + combineHLT.get("Run2018.MC", "") + ")"

   if options.era=="2016":
      modules_era.append(puAutoWeight_2016())
      modules_era.append(PrefCorr())
      modules_era.append(jetmetCorrector())
      modules_era.append(btagSFProducer("Legacy2016", "deepjet",['shape_corr'],jesSystsForShape=['jesAbsolute', 'jesAbsolute_2016', 'jesBBEC1', 'jesBBEC1_2016', 'jesEC2', 'jesEC2_2016', 'jesFlavorQCD', 'jesHF', 'jesHF_2016', 'jesRelativeBal', 'jesRelativeSample_2016']))
      modules_era.append(muonScaleRes2016())
      modules_era.append(lepSF_2016())
      # modules_era.append(nvtxWeight_2016())
      ext_syst.append("PrefireWeight")
   if options.era=="2017":
      modules_era.append(puAutoWeight_2017())
      modules_era.append(PrefCorr())
      modules_era.append(jetmetCorrector())
      modules_era.append(btagSFProducer("2017", "deepjet",['shape_corr'],jesSystsForShape=['jesAbsolute', 'jesAbsolute_2017', 'jesBBEC1', 'jesBBEC1_2017', 'jesEC2', 'jesEC2_2017', 'jesFlavorQCD', 'jesHF', 'jesHF_2017', 'jesRelativeBal', 'jesRelativeSample_2017']))
      modules_era.append(muonScaleRes2017())
      modules_era.append(lepSF_2017())
      # modules_era.append(nvtxWeight_2017())
      ext_syst.append("PrefireWeight")
   if options.era=="2018":
      modules_era.append(puAutoWeight_2018())
      modules_era.append(jetmetCorrector())
      modules_era.append(btagSFProducer("2018","deepjet",['shape_corr'],jesSystsForShape=['jesAbsolute', 'jesAbsolute_2018', 'jesBBEC1', 'jesBBEC1_2018', 'jesEC2', 'jesEC2_2018', 'jesFlavorQCD', 'jesHF', 'jesHF_2018', 'jesRelativeBal', 'jesRelativeSample_2018','jesHEMIssue']))
      modules_era.append(muonScaleRes2018())
      modules_era.append(lepSF_2018())
      # modules_era.append(nvtxWeight_2018())

   # modules_era.append(PhiXYCorrection(era=options.era,isMC=options.isMC,sys=''))
   modules_era.append(DiHiggsProducer(isMC=options.isMC, era=str(options.era), do_syst=0, signal=options.signal, syst_var=''))

   modules_era.append(AngularVariablesProducerForHH(signal=options.signal))

   if "TTTo2L2Nu" in options.dataset or "TTToHadronic" in options.dataset or "TTToSemiLeptonic" in options.dataset or "TTJets" in options.dataset or "TT_" in options.dataset:
      modules_era.append(GenTopProducer())

   if options.era=="2016":
      modules_era.append(TriggerSF_2016(syst=''))
   if options.era=="2017":
      modules_era.append(TriggerSF_2017(syst=''))
   if options.era=="2018":
      modules_era.append(TriggerSF_2018(syst=''))



   # for shift-based systematics
   if options.doSyst:
      for sys in pro_syst:
         for var in ["Up", "Down"]:
	    # if "jesTotal" in sys and options.doSyst==1: modules_era.append(PhiXYCorrection(era=options.era,isMC=options.isMC,sys=sys+var))
	    # if "jer" in sys and options.doSyst==1: modules_era.append(PhiXYCorrection(era=options.era,isMC=options.isMC,sys=sys+var))
            modules_era.append(DiHiggsProducer(options.isMC, str(options.era), do_syst=options.doSyst, signal=options.signal, syst_var=sys+var))
            modules_era.append(AngularVariablesProducerForHH(do_syst=options.doSyst, signal=options.signal, syst_var=sys+var))

            if options.era=="2016":
               modules_era.append(TriggerSF_2016(syst=sys+var))
            if options.era=="2017":
               modules_era.append(TriggerSF_2017(syst=sys+var))
            if options.era=="2018":
               modules_era.append(TriggerSF_2018(syst=sys+var))
#            modules_era.append(BDTdiscriminantProducerForHH(era=str(options.era), do_syst=options.doSyst, syst_var=sys+var))
   modules_era.append(StoreEventsProducerForHH(str(options.era),options.doSyst))

else:
   if options.era=="2016":
      modules_era.append(muonScaleRes2016())
   if options.era=="2017":
      modules_era.append(muonScaleRes2017())
   if options.era=="2018":
      modules_era.append(muonScaleRes2018())

   print "sample : ", options.dataset, " candtag : ", condtag_
   try:
      combineHLT = yaml.safe_load(open("combineHLT_Run2.yaml"))
   except yaml.YAMLError as exc:
      print(exc)
   if options.era=="2016":
      if 'Run2016H' in condtag_:
	pre_selection = pre_selection + " && (" + combineHLT.get("Run2016H.%s" % options.dataset, "") + ")"
      else:
        pre_selection = pre_selection + " && (" + combineHLT.get("Run2016All.%s" % options.dataset, "") + ")"
   if options.era=="2017":
      if 'Run2017B' in condtag_:
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2017B.%s" % options.dataset, "") + ")"
      elif 'Run2017C' in condtag_:
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2017C.%s" % options.dataset, "") + ")"
      else:
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2017DF.%s" % options.dataset, "") + ")"
   if options.era=="2018":
      if ('Run2018A' in condtag_) or ('Run2018B' in condtag_):
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2018AB.%s" % options.dataset, "") + ")"
      else:
         pre_selection = pre_selection + " && (" + combineHLT.get("Run2018CD.%s" % options.dataset, "") + ")"

   print " -- era : ",

   runPeriod = condtag_.split(options.era)[1][:1]
   jmeCorrections = createJMECorrector(isMC=False, 
                                       dataYear=options.era, 
                                       runPeriod=runPeriod, 
                                       jesUncert="Total", 
                                       metBranchName=metBranchName)   

   modules_era.append(jetmetCorrector())
  
#   if options.era=="2016":
#      if 'Run2016G' or 'Run2016H' in condtag_:
#         modules_era.append(getattr(jetRecalib, 'jetRecalib2016GH')() )
#      elif 'Run2016E' or 'Run2016F' in condtag_:
#         modules_era.append(getattr(jetRecalib, 'jetRecalib2016EF')() )
#      else:
#         modules_era.append(getattr(jetRecalib('AK4PFchs'), 'jetRecalib2016BCD')() )
#   if options.era=="2017":
#      if 'Run2017D' or 'Run2017E' in condtag_:
#         modules_era.append(getattr(jetRecalib, 'jetRecalib2017DE')() )
#      else:
#         modules_era.append(getattr(jetRecalib, 'jetRecalib2017%s' % condtag_.split(options.era)[1])() )
#   if options.era=="2018":
#      modules_era.append(getattr(jetRecalib, 'jetRecalib2018%s' % condtag_.split(options.era)[1][:1])() )

#   modules_era.append(PhiXYCorrection(era=options.era,isMC=options.isMC,sys=''))
   modules_era.append(DiHiggsProducer  (isMC=options.isMC, era=str(options.era), do_syst=0, signal=options.signal, syst_var=''))
   modules_era.append(AngularVariablesProducerForHH(signal=options.signal))
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
