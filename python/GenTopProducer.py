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


class GenTopProducer(Module):
    def __init__(self):
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("ttbarweight_nominal", "F")
        self.out.branch("ttbarweight_up", "F")
        self.out.branch("ttbarweight_down", "F")
 
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        # only valid in the MC samples

        gen_part = Collection(event, "GenPart")
        p4Gen_top = ROOT.TLorentzVector()
        p4Gen_antitop = ROOT.TLorentzVector()
        

        for part in gen_part:
            if (part.statusFlags & 256) == 0:#isLastCopy
                continue
            if (part.statusFlags & 8192) == 0:#fromHardProcess
                continue
            if (part.pdgId) == 6:
                p4Gen_top = part.p4()
            if (part.pdgId) == -6:
                p4Gen_antitop = part.p4()

        top_weight = math.exp(-2.02274e-01 + 1.09734e-04*p4Gen_top.Pt() + -1.30088e-07*p4Gen_top.Pt()**2 + (5.83494e+01/(p4Gen_top.Pt()+1.96252e+02)))
        antitop_weight = math.exp(-2.02274e-01 + 1.09734e-04*p4Gen_antitop.Pt()+ -1.30088e-07*p4Gen_antitop.Pt()**2 + (5.83494e+01/(p4Gen_antitop.Pt()+1.96252e+02)))

        if top_weight*antitop_weight > 0:
            reweight = math.sqrt(top_weight*antitop_weight)
        else:
            reweight = 1

        self.out.fillBranch("ttbarweight_nominal", reweight) # apply once
        self.out.fillBranch("ttbarweight_up", reweight*reweight) ## apply twice
        self.out.fillBranch("ttbarweight_down", 1) # no apply

        return True
