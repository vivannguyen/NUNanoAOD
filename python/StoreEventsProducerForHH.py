import ROOT
import sys
import numpy as np
import math
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import getValueReader
import PhysicsTools.NanoAODTools.postprocessing.tools as tk

ROOT.PyConfig.IgnoreCommandLineOptions = True

class StoreEventsProducerForHH(Module):
    def __init__(self, era, do_syst = False):
        self.era = era
        self.do_syst = do_syst

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        goodEvent = int(getattr(event, "good_event"))
        if goodEvent==1:
            return True

        if self.do_syst==False and goodEvent==0:
            return False

        pro_syst = [ "ElectronEn", "MuonEn", "jer"]        

        if self.era=="2016":
            for ix in ['jesAbsolute', 'jesAbsolute_2016', 'jesBBEC1', 'jesBBEC1_2016', 'jesEC2', 'jesEC2_2016', 'jesFlavorQCD', 'jesHF', 'jesHF_2016', 'jesRelativeBal', 'jesRelativeSample_2016']:
                pro_syst.append(ix)
        if self.era=="2017":
            for ix in ['jesAbsolute', 'jesAbsolute_2017', 'jesBBEC1', 'jesBBEC1_2017', 'jesEC2', 'jesEC2_2017', 'jesFlavorQCD', 'jesHF', 'jesHF_2017', 'jesRelativeBal', 'jesRelativeSample_2017']:
                pro_syst.append(ix)
        if self.era=="2018":
            for ix in ['jesAbsolute', 'jesAbsolute_2018', 'jesBBEC1', 'jesBBEC1_2018', 'jesEC2', 'jesEC2_2018', 'jesFlavorQCD', 'jesHF', 'jesHF_2018', 'jesRelativeBal', 'jesRelativeSample_2018','jesHEMIssue']:
                pro_syst.append(ix)

        for sys in pro_syst:
            for var in ["Up", "Down"]:
                syst_name = sys+var
                eventGood = int(getattr(event, "good_event_sys_{}".format(syst_name)))

                if eventGood==1:
                    return True

        return False



