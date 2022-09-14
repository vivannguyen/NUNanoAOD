import ROOT
import sys
import numpy as np
import math
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import getValueReader

ROOT.PyConfig.IgnoreCommandLineOptions = True


class GenHHProducer(Module):
    def __init__(self, isMC, era, do_syst=False, syst_var=''):
        self.isMC = isMC
        self.era = era
        self.do_syst = do_syst
        self.syst_var = syst_var
        self.zmass = 91.1873
        self.Hmass = 125.10
        if self.syst_var !='':
          self.syst_suffix = '_sys_' + self.syst_var if self.do_syst else ''
        else:
          self.syst_suffix = syst_var

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def getCosThetaStar_CS(self, h1, h2):
        #cos theta star angle in the Collins Soper frame
        p1, p2, hh = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()
        p1.SetPxPyPzE(0, 0,  6500., 6500)
        p2.SetPxPyPzE(0, 0, -6500., 6500)
        hh = h1 + h2
        boost = ROOT.TVector3(- hh.BoostVector())
        #boost = - hh.BoostVector();
        p1.Boost(boost)
        p2.Boost(boost)
        h1.Boost(boost)
        #TVector3 CSaxis = p1.Vect().Unit() - p2.Vect().Unit()
        CSaxis = ROOT.TVector3(p1.Vect().Unit() - p2.Vect().Unit())
        CSaxis.Unit()

        return abs(math.cos(   CSaxis.Angle( h1.Vect().Unit() )    ))

    def getCosThetaStar_CS_new(self, h1, h2):
        hh_lor = h1 + h2;
        hh = h1 + h2;

        h1_lor = h1;
        h_1 = h1

        h_1.Boost(-hh.BoostVector());   
        
        return abs(h_1.CosTheta());

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("gen_HH_pt", "F")
        self.out.branch("gen_HH_eta", "F")
        self.out.branch("gen_HH_phi", "F")
        self.out.branch("gen_HH_mass", "F")
        self.out.branch("gen_HH_cos_CS", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        # only valid in the MC samples

        gen_part = Collection(event, "GenPart")
        
        genHiggs = []
        for part in gen_part:
            if(abs(part.pdgId)==25 and part.genPartIdxMother==0):
                genHiggs.append(part.p4())
        
        genHH = ROOT.TLorentzVector()
        genHH = genHiggs[0] + genHiggs[1]
        
        self.out.fillBranch("gen_HH_pt", genHH.Pt())
        self.out.fillBranch("gen_HH_mass", genHH.M())
        self.out.fillBranch("gen_HH_cos_CS", self.getCosThetaStar_CS_new(genHiggs[0],genHiggs[1]))
        self.out.fillBranch("gen_HH_phi", genHH.Phi())
        self.out.fillBranch("gen_HH_eta", genHH.Eta())


        return True
