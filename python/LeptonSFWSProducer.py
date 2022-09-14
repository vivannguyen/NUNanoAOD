import ROOT
from importlib import import_module
import numpy as np
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import array as ar

ROOT.PyConfig.IgnoreCommandLineOptions = True

class LeptonSFWSProducer(Module):
    def __init__(self, isMC, era=2017, sample="DY", do_syst=False, syst_var='', weight_syst=False):
        self.isMC = isMC
        self.era = era
        self.sample = sample
        self.do_syst = do_syst
        self.syst_var = syst_var
        self.weight_syst = weight_syst
        if self.syst_var !='':
            self.syst_suffix = '_sys_' + syst_var if do_syst else ''
        else:
            self.syst_suffix = syst_var
        self.writeHistFile=True

    def passbut(self, event, excut=None, cat="signal"):
        event_pass = True
        for cut in self.selection[cat]:
            if self.weight_syst:
                cut = cut.format(sys="")
            else:
                cut = cut.format(sys=self.syst_suffix)
            if excut is not None:
                if excut in cut:
                    continue
            if eval(cut) is False:
                event_pass = False
        return event_pass

    def beginJob(self, histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)
        self.cats = {
            1  : "cat_BB" ,
            2  : "cat_BE" ,
            3  : "cat_EB" ,
            4  : "cat_EE" ,
            5  : "cat_BB_num",
            6  : "cat_BE_num",
            7  : "cat_EB_num",
            8  : "cat_EE_num",
        }
        self.selection = {
            "signal" : [
                "event.met_pt > 40",
                "event.event_category == 1",
                "abs(event.DY_mass - 91) > 15",
                "event.DY_mass > 15 "
            ]
        }

        self.h_DYmass = ROOT.TH1F(
            'DY_mass{}{}'.format("_" + self.sample, self.syst_suffix),
            'DY_mass{}{}'.format("_" + self.sample, self.syst_suffix),
            25, 0, 250
        )
        self.h_muon = {}
        self.h_electron = {}
        for i,cat in self.cats.items():
            self.h_muon[i] = ROOT.TH2F(
                'h_muon{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_muon{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                8, ar.array('d',[20,25,30,35,40,50,60,70,90]),
                10, ar.array('d',[10,15,20,25,30,35,40,50,60,70,90])
            )
            self.h_electron[i] = ROOT.TH2F(
                'h_electron{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_electron{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                7, ar.array('d',[25,30,35,40,50,60,70,90]),
                9, ar.array('d',[15,20,25,30,35,40,50,60,70,90])
            )

#        self.h_met = {}
#        self.h_mT = {}
#        for i,cat in self.cats.items():
#	    self.h_mT[i] = ROOT.TH1F(
#                'MT{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
#                'MT{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
#		13, ar.array('d', [0,100,200,250,300,350,400,500,600,700,800,1000,1200,2000])
#            )
            # different binning for different regions
#            if cat == 'catNRB' or cat=="catTOP" or cat=="catDY":
#                self.h_met[i] = ROOT.TH1F(
#                    'measMET{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
#                    'measMET{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
#                    5, ar.array('d', [50,60,70,80,90,100])
#                )
#            else:
#                self.h_met[i] = ROOT.TH1F(
#                    'measMET{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
#                    'measMET{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
#                    12, ar.array('d', [50,100,125,150,175,200,250,300,350,400,500,600,1000])
#                )

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        prevdir = ROOT.gDirectory
        outputFile.cd()
        for i,cat in self.cats.items():
            self.h_electron[i].Write()
            self.h_muon[i].Write()
        self.h_DYmass.Write()
        prevdir.cd()

    def analyze(self, event):
        # only valid in the MC samples
        try:
            lep_category = getattr(event, "lep_category{}".format(self.syst_suffix))
            event_category = getattr(event, "event_category{}".format(self.syst_suffix))
            leading_lep_pt = getattr(event,"leading_lep_pt{}".format(self.syst_suffix))
            leading_lep_eta = getattr(event,"leading_lep_eta{}".format(self.syst_suffix))
            trailing_lep_pt = getattr(event,"trailing_lep_pt{}".format(self.syst_suffix))
            trailing_lep_eta = getattr(event,"trailing_lep_eta{}".format(self.syst_suffix))
            electron_trigger = getattr(event,"electron_trigger{}".format(self.syst_suffix))
            muon_trigger = getattr(event,"muon_trigger{}".format(self.syst_suffix))
        except:
            lep_category = getattr(event, "lep_category")
            event_category = getattr(event, "event_category")
            leading_lep_pt = getattr(event,"leading_lep_pt")
            leading_lep_eta = getattr(event,"leading_lep_eta")
            trailing_lep_pt = getattr(event,"trailing_lep_pt")
            trailing_lep_eta = getattr(event,"trailing_lep_eta")
            electron_trigger = getattr(event,"electron_trigger")
            muon_trigger = getattr(event,"muon_trigger")
        if self.weight_syst:
            meas_MET 	 = getattr(event, "met_pt")
        else:
            meas_MET 	 = getattr(event, "met_pt{}".format(self.syst_suffix))

        # cross-section
        weight = 1.0
        try:
            weight = getattr(event, "xsecscale")
        except:
            return "ERROR: weight branch doesn't exist"


        #weight *= event.ADDWeight #for ADD samples only (EFT weights)
        # pu uncertainty
        if self.isMC:

            try:
                if "ADD" in self.syst_suffix:
		    weight *= event.ADDWeight #for ADD samples only (EFT weights)
            except:
                pass
            # weight *= event.genWeight
            if "puWeight" in self.syst_suffix:
                if "Up" in self.syst_suffix:
                    weight *= event.puWeightUp
                else:
                    weight *= event.puWeightDown
            else:
            	weight *= event.puWeight
                  
            if "MuonSF" in self.syst_suffix:
                if "Up" in self.syst_suffix:
                    weight *= event.w_muon_SFUp
                else:
                    weight *= event.w_muon_SFDown
            else:
                weight *= event.w_muon_SF
            # Electron SF
            if "ElecronSF" in self.syst_suffix:
                if "Up" in self.syst_suffix:
                    weight *= event.w_electron_SFUp
                else:
                    weight *= event.w_electron_SFDown
            else:
                weight *= event.w_electron_SF
	    #Prefire Weight
            try:
                if "PrefireWeight" in self.syst_suffix:
                    if "Up" in self.syst_suffix:
                        weight *= event.PrefireWeight_Up
                    else:
                        weight *= event.PrefireWeight_Down
                else:
                    weight *= event.PrefireWeight
            except:
                pass
        if (weight == 0.):
            print ' problem'
        if ( (lep_category == 1) and self.passbut(event, cat="signal") ):
            if (abs(leading_lep_eta)<1.5 and abs(trailing_lep_eta) < 1.5):
                self.h_muon[1].Fill(leading_lep_pt,trailing_lep_pt,weight)
                if(muon_trigger == 1):
                    self.h_muon[5].Fill(leading_lep_pt,trailing_lep_pt,weight)
            if (abs(leading_lep_eta)<1.5 and abs(trailing_lep_eta) > 1.5):
                self.h_muon[2].Fill(leading_lep_pt,trailing_lep_pt,weight)
                if(muon_trigger == 1):
                    self.h_muon[6].Fill(leading_lep_pt,trailing_lep_pt,weight)
            if (abs(leading_lep_eta)>1.5 and abs(trailing_lep_eta) < 1.5):
                self.h_muon[3].Fill(leading_lep_pt,trailing_lep_pt,weight)
                if(muon_trigger == 1):
                    self.h_muon[7].Fill(leading_lep_pt,trailing_lep_pt,weight)
            if (abs(leading_lep_eta)>1.5 and abs(trailing_lep_eta) > 1.5):
                self.h_muon[4].Fill(leading_lep_pt,trailing_lep_pt,weight)
                if(muon_trigger == 1):
                    self.h_muon[8].Fill(leading_lep_pt,trailing_lep_pt,weight)

        if ( (lep_category == 2) and self.passbut(event, cat="signal") ):
            if (abs(leading_lep_eta)<1.5 and abs(trailing_lep_eta) < 1.5):
                self.h_electron[1].Fill(leading_lep_pt,trailing_lep_pt,weight)
                if(electron_trigger == 1):
                    self.h_electron[5].Fill(leading_lep_pt,trailing_lep_pt,weight)
            if (abs(leading_lep_eta)<1.5 and abs(trailing_lep_eta) > 1.5):
                self.h_electron[2].Fill(leading_lep_pt,trailing_lep_pt,weight)
                if(electron_trigger == 1):
                    self.h_electron[6].Fill(leading_lep_pt,trailing_lep_pt,weight)
            if (abs(leading_lep_eta)>1.5 and abs(trailing_lep_eta) < 1.5):
                self.h_electron[3].Fill(leading_lep_pt,trailing_lep_pt,weight)
                if(electron_trigger == 1):
                    self.h_electron[7].Fill(leading_lep_pt,trailing_lep_pt,weight)
            if (abs(leading_lep_eta)>1.5 and abs(trailing_lep_eta) > 1.5):
                self.h_electron[4].Fill(leading_lep_pt,trailing_lep_pt,weight)
                if(electron_trigger == 1):
                    self.h_electron[8].Fill(leading_lep_pt,trailing_lep_pt,weight)
        return True
