import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class TriggerSFProducerForHH(Module):
    def __init__(self,targetfile,name="TriggerSFWeight",norm=True,verbose=False,doSysVar=True,syst_var=''):
        self.targetfile = targetfile
        self.name = name
        self.norm = norm
        self.verbose = verbose
        self.doSysVar = doSysVar
        self.syst_var = syst_var
        self.do_syst = True
        if self.syst_var !='':
          self.syst_suffix = '_sys_' + self.syst_var if self.do_syst else ''
        else:
          self.syst_suffix = syst_var


    def loadHisto(self,filename,hname):
        tf = ROOT.TFile.Open(filename)
        hist = tf.Get(hname)
        hist.SetDirectory(None)
        tf.Close()
        return hist
    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
#        print "begin...", self.name
        self.out = wrappedOutputTree
        self.out.branch(self.name+"{}".format(self.syst_suffix), "F")
        self.out.branch(self.name+"{}_Up".format(self.syst_suffix), "F")
        self.out.branch(self.name+"{}_Down".format(self.syst_suffix), "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
 #       print "here we are.."
        """process event, return True (go to next module) or False (fail, go to next event)"""
        lep_cat = 0
        ev_cat = 0
        good_event = 0
        if hasattr(event,"lep_category{}".format(self.syst_suffix)):
            lep_cat = int(getattr(event,"lep_category{}".format(self.syst_suffix)))
        if hasattr(event,"event_category{}".format(self.syst_suffix)):
            ev_cat = int(getattr(event,"event_category{}".format(self.syst_suffix)))
        if hasattr(event,"good_event{}".format(self.syst_suffix)):
            good_event = int(getattr(event,"good_event{}".format(self.syst_suffix)))
        if lep_cat < 1 :
            l1_pt = 0
            l1_eta = 0
            l2_pt = 0
            l2_eta = 0
            l1_flavor = 0
        else:
            l1_pt = float(getattr(event,"leading_lep_pt{}".format(self.syst_suffix)))
            l1_eta = abs(float(getattr(event,"leading_lep_eta{}".format(self.syst_suffix))))
            l2_pt = float(getattr(event,"trailing_lep_pt{}".format(self.syst_suffix)))
            l2_eta = abs(float(getattr(event,"trailing_lep_eta{}".format(self.syst_suffix))))
            l1_flavor = int(getattr(event,"leading_lep_flavor{}".format(self.syst_suffix)))

        weight = 1
        weightError = 0

        if lep_cat==1:# or lep_cat==5 or lep_cat==7 : #these are MM. MML and MMLL lepton categories
            hist = self.loadHisto(self.targetfile,"muon-AA") # trigger SF are not eta dependent because of low stats
#	    if l1_eta <= 1.5 and l2_eta <= 1.5:
#	    	hist = self.loadHisto(self.targetfile,"muon-BB")
#            elif l1_eta >= 1.5 and l2_eta <= 1.5:
#                hist = self.loadHisto(self.targetfile,"muon-EB")
#            elif l1_eta <= 1.5 and l2_eta >= 1.5:
#                hist = self.loadHisto(self.targetfile,"muon-BE")
#            elif l1_eta >= 1.5 and l2_eta >= 1.5:
#                hist = self.loadHisto(self.targetfile,"muon-EE")
        elif lep_cat==2:# or lep_cat==4 or lep_cat==6 : #these are EE. EEL and EELL lepton categories
            hist = self.loadHisto(self.targetfile,"electron-AA")
#            if l1_eta <= 1.5 and l2_eta <= 1.5:
#                hist = self.loadHisto(self.targetfile,"electron-BB")
#            elif l1_eta >= 1.5 and l2_eta <= 1.5:
#                hist = self.loadHisto(self.targetfile,"electron-EB")
#            elif l1_eta <= 1.5 and l2_eta >= 1.5:
#                hist = self.loadHisto(self.targetfile,"electron-BE")
#            elif l1_eta >= 1.5 and l2_eta >= 1.5:
#                hist = self.loadHisto(self.targetfile,"electron-EE")

        if lep_cat==1:
            nxBins = 8
            nyBins = 10
        else:
            nxBins = 7
            nyBins = 9

	searchbinx = -1
	searchbiny = -1
        if good_event == 1:
            for xbin in range(1,nxBins+1):
                if l1_pt > hist.GetXaxis().GetBinLowEdge(nxBins) + hist.GetXaxis().GetBinWidth(nxBins):
                    if lep_cat==1:
                        searchbinx = 8
                    else:
                        searchbinx = 7
                    break
                if l1_pt > hist.GetXaxis().GetBinLowEdge(xbin) and l1_pt < hist.GetXaxis().GetBinLowEdge(xbin) + hist.GetXaxis().GetBinWidth(xbin) :
                    searchbinx = xbin
            for ybin in range(1,nyBins+1):
                if l2_pt > hist.GetYaxis().GetBinLowEdge(nyBins) + hist.GetYaxis().GetBinWidth(nyBins):
                    if lep_cat==1:
                        searchbiny = 7
                    else:
                        searchbiny = 9
                    break
                if l2_pt > hist.GetYaxis().GetBinLowEdge(ybin) and l2_pt < hist.GetYaxis().GetBinLowEdge(ybin) + hist.GetYaxis().GetBinWidth(ybin) :
                    searchbiny = ybin

            weight = hist.GetBinContent(searchbinx,searchbiny)
            self.out.fillBranch(self.name+"{}".format(self.syst_suffix), weight)
            if self.doSysVar:
                weightError = hist.GetBinErrorUp(searchbinx,searchbiny)
                self.out.fillBranch(self.name+"{}_Up".format(self.syst_suffix), weight+weightError)
                weightError = hist.GetBinErrorLow(searchbinx,searchbiny)
                self.out.fillBranch(self.name+"{}_Down".format(self.syst_suffix), weight-weightError)
        else:
            self.out.fillBranch(self.name+"{}".format(self.syst_suffix), -99)
            self.out.fillBranch(self.name+"{}_Up".format(self.syst_suffix), -99)
            self.out.fillBranch(self.name+"{}_Down".format(self.syst_suffix), -99)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
TrigSF_2016 = "%s/src/PhysicsTools/MonoZ/data/TriggerSFs/trigger-2016-miniISO.root" % os.environ['CMSSW_BASE']
TriggerSF_2016 = lambda syst : TriggerSFProducerForHH(TrigSF_2016,verbose=False, doSysVar=True , syst_var=syst)

TrigSF_2017 = "%s/src/PhysicsTools/MonoZ/data/TriggerSFs/trigger-2017-miniISO.root" % os.environ['CMSSW_BASE']
TriggerSF_2017 = lambda syst : TriggerSFProducerForHH(TrigSF_2017,verbose=False, doSysVar=True, syst_var=syst)

TrigSF_2018 = "%s/src/PhysicsTools/MonoZ/data/TriggerSFs/trigger-2018-miniISO.root" % os.environ['CMSSW_BASE']
TriggerSF_2018 = lambda syst : TriggerSFProducerForHH(TrigSF_2018,verbose=False, doSysVar=True, syst_var=syst)

