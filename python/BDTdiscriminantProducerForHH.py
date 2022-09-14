import ROOT
import os
import numpy as np
import array
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class BDTdiscriminantProducerForHH(Module):
    def __init__(self, era,name="BDTdiscrim", verbose=False, do_syst=False, syst_var=''):
        self.era = era
        self.targetfileMuon = "%s/src/PhysicsTools/MonoZ/data/TMVA/TMVAClassification_BDT_M260_mm_%s.weights.xml"% (os.environ['CMSSW_BASE'], self.era)
        self.targetfileElectron = "%s/src/PhysicsTools/MonoZ/data/TMVA/TMVAClassification_BDT_M260_ee_%s.weights.xml"% (os.environ['CMSSW_BASE'], self.era)
        self.name = name
        self.verbose = verbose
        self.syst_var = syst_var
        self.do_syst = do_syst
        if self.syst_var !='':
          self.syst_suffix = '_sys_' + self.syst_var if self.do_syst else ''
        else:
          self.syst_suffix = syst_var


        
        self._bdtvars = ['Mbb_H','DR_bb_H','Mjj_Z','DR_jj_Z','M_uu','DR_muon1muon2','DR_uu_bb_H','DR_u1Hj1','DR_u1Hj2','DR_u2Hj1','DR_u2Hj2','DR_uu_jj_Z','DR_u1Zj1','DR_u1Zj2','DR_u2Zj1','DR_u2Zj2','abs(cosThetaStarMu)','abs(cosTheta_hbb_uu)','abs(cosTheta_zuu_hzz)','abs(DPhi_muon1met)','abs(phi1_uu)','abs(phi1_zjj_uu)']
        self._bdtvarnames = {}
#        self._bdtvarnamesElectron = {}
#        self._testvars = ['lep_category']
#        self._testvarnames = {}
#        for ix in self._testvars:
#            self._testvarnames[ix] = array.array('f',[0])

        for ix in self._bdtvars:
            self._bdtvarnames[ix] = array.array('f',[0])        

        self.BDTreader_Muon = self.loadBDT(self.targetfileMuon,"mm")

#        self.BDTreader_Electron = self.loadBDT(self.targetfileElectron,"ee")
# as placeholder let use 
        self.BDTreader_Electron = self.loadBDT(self.targetfileMuon,"ee")
        
    def loadBDT(self,filename,channel):
        reader = ROOT.TMVA.Reader("!Color")
  #      _bdtvars = ['Mbb_H','DR_bb_H','Mjj_Z','DR_jj_Z','M_uu','DR_muon1muon2','DR_uu_bb_H','DR_u1Hj1','DR_u1Hj2','DR_u2Hj1','DR_u2Hj2','DR_uu_jj_Z','DR_u1Zj1','DR_u1Zj2','DR_u2Zj1','DR_u2Zj2','abs(cosThetaStarMu)','abs(cosTheta_hbb_uu)','abs(cosTheta_zuu_hzz)','abs(DPhi_muon1met)','abs(phi1_uu)','abs(phi1_zjj_uu)']
  #      _bdtvarnames = {}
        for ix in self._bdtvars:
            reader.AddVariable(ix, self._bdtvarnames[ix])
        reader.BookMVA("BDT_classifier_%s_%s"%(channel,self.syst_suffix),filename)
        return reader

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("%s%s"%(self.name,self.syst_suffix),"F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        lep_cat = 0
        if hasattr(event,"lep_category"):
            lep_cat = int(getattr(event,"lep_category"))

        
 #       for ix in self._testvars:
 #           self._testvarnames[ix][0] = float(getattr(event,ix))
 #           print'test pozor: ',self._testvarnames[ix][0]
   
        if lep_cat == 1:
            BDTreader = self.BDTreader_Muon
            channel = "mm"
        if lep_cat == 2:
            BDTreader = self.BDTreader_Electron
            channel = "ee"


        self._bdtvarnames['Mbb_H'][0] = -99.0
        self._bdtvarnames['DR_bb_H'][0] = -99.0
        self._bdtvarnames['Mjj_Z'][0] = -99.0
        self._bdtvarnames['DR_jj_Z'][0] = -99.0
        self._bdtvarnames['M_uu'][0] = -99.0
        self._bdtvarnames['DR_muon1muon2'][0] = -99.0
        self._bdtvarnames['DR_uu_bb_H'][0] = -99.0
        self._bdtvarnames['DR_u1Hj1'][0] = -99.0
        self._bdtvarnames['DR_u1Hj2'][0] = -99.0
        self._bdtvarnames['DR_u2Hj1'][0] = -99.0
        self._bdtvarnames['DR_u2Hj2'][0] = -99.0
        self._bdtvarnames['DR_uu_jj_Z'][0] = -99.0
        self._bdtvarnames['DR_u1Zj1'][0] = -99.0
        self._bdtvarnames['DR_u1Zj2'][0] = -99.0
        self._bdtvarnames['DR_u2Zj1'][0] = -99.0
        self._bdtvarnames['DR_u2Zj2'][0] = -99.0
        self._bdtvarnames['abs(cosThetaStarMu)'][0] = -99.0
        self._bdtvarnames['abs(cosTheta_hbb_uu)'][0] = -99.0
        self._bdtvarnames['abs(cosTheta_zuu_hzz)'][0] = -99.0
        self._bdtvarnames['abs(DPhi_muon1met)'][0] = -99.0
        self._bdtvarnames['abs(phi1_uu)'][0] = -99.0
        self._bdtvarnames['abs(phi1_zjj_uu)'][0] = -99.0

        bdtdisc = BDTreader.EvaluateMVA("BDT_classifier_%s_%s"%(channel,self.syst_suffix))
        self.out.fillBranch("%s%s"%(self.name,self.syst_suffix), bdtdisc)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
#BDTdiscrim_2016 = lambda : BDTdiscriminantProducerForHH(era="2016", verbose=False)
#BDTdiscrim_2017 = lambda : BDTdiscriminantProducerForHH(era="2017", verbose=False)
#BDTdiscrim_2018 = lambda : BDTdiscriminantProducerForHH(era="2018", verbose=False)

