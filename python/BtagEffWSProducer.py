import ROOT
from importlib import import_module
import numpy as np
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import array as ar

ROOT.PyConfig.IgnoreCommandLineOptions = True

class BtagEffWSProducer(Module):
    def __init__(self, isMC, era, sample="DY", do_syst=False, syst_var='', weight_syst=False):
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


    def btagDeepJet_id(self, wp):
        # ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        if (self.era == "2016" and wp == "loose"):
            return 0.0614
        elif (self.era == "2016" and wp == "medium"):
            return 0.3093
        elif (self.era == "2016" and wp == "tight"):
            return 0.7221
        elif (self.era == "2017" and wp == "loose"):
            return 0.0521
        elif (self.era == "2017" and wp == "medium"):
            return 0.3033
        elif (self.era == "2017" and wp == "tight"):
            return 0.7489
        elif (self.era == "2018" and wp == "loose"):
            return 0.0494
        elif (self.era == "2018" and wp == "medium"):
            return 0.2770
        elif (self.era == "2018" and wp == "tight"):
            return 0.7264

    def btag_id(self, wp):
        # ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        if (self.era == "2016" and wp == "loose"):
            return 0.2219
        elif (self.era == "2016" and wp == "medium"):
            return 0.6324
        elif (self.era == "2016" and wp == "tight"):
            return 0.8958
        elif (self.era == "2017" and wp == "loose"):
            return 0.1522
        elif (self.era == "2017" and wp == "medium"):
            return 0.4941
        elif (self.era == "2017" and wp == "tight"):
            return 0.8001
        elif (self.era == "2018" and wp == "loose"):
            return 0.1241
        elif (self.era == "2018" and wp == "medium"):
            return 0.4184
        elif (self.era == "2018" and wp == "tight"):
            return 0.7527


    def beginJob(self, histFile=None,histDirName=None):
        Module.beginJob(self,histFile,histDirName)
        self.cats = {
            1  : "cat_Bottom" ,
            2  : "cat_Charm" ,
            3  : "cat_Light" ,
        }
        self.selection = {
            "signal" : [
                "event.event_category == 1",
                "event.Zlep_cand_mass > 15 ",
                "event.leading_lep_pt > 30",
                "event.trailing_lep_pt > 20"
#                "event.Zlep_cand_pt > 15 "
            ]
        }

        self.h_muon_num = {}
        self.h_electron_num = {}
        self.h_muon_den = {}
        self.h_electron_den = {}

        self.h_muon_num_csv = {}
        self.h_electron_num_csv = {}
        self.h_muon_den_csv = {}
        self.h_electron_den_csv = {}

        for i,cat in self.cats.items():
            self.h_muon_num[i] = ROOT.TH2F(
                'h_muon_num{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_muon_num{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                9, ar.array('d',[20,30,50,70,100,140,200,300,600,1000]),
                5, ar.array('d',[0.,0.5,1.,1.5,2.0,2.4])
            )
            self.h_electron_num[i] = ROOT.TH2F(
                'h_electron_num{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_electron_num{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                9, ar.array('d',[20,30,50,70,100,140,200,300,600,1000]),
                5, ar.array('d',[0.,0.5,1.,1.5,2.0,2.4])
            )
            self.h_muon_den[i] = ROOT.TH2F(
                'h_muon_den{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_muon_den{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                9, ar.array('d',[20,30,50,70,100,140,200,300,600,1000]),
                5, ar.array('d',[0.,0.5,1.,1.5,2.0,2.4])
            )
            self.h_electron_den[i] = ROOT.TH2F(
                'h_electron_den{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_electron_den{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                9, ar.array('d',[20,30,50,70,100,140,200,300,600,1000]),
                5, ar.array('d',[0.,0.5,1.,1.5,2.0,2.4])
            )

            self.h_muon_num_csv[i] = ROOT.TH2F(
                'h_muon_num_csv{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_muon_num_csv{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                9, ar.array('d',[20,30,50,70,100,140,200,300,600,1000]),
                5, ar.array('d',[0.,0.5,1.,1.5,2.0,2.4])
            )
            self.h_electron_num_csv[i] = ROOT.TH2F(
                'h_electron_num_csv{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_electron_num_csv{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                9, ar.array('d',[20,30,50,70,100,140,200,300,600,1000]),
                5, ar.array('d',[0.,0.5,1.,1.5,2.0,2.4])
            )
            self.h_muon_den_csv[i] = ROOT.TH2F(
                'h_muon_den_csv{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_muon_den_csv{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                9, ar.array('d',[20,30,50,70,100,140,200,300,600,1000]),
                5, ar.array('d',[0.,0.5,1.,1.5,2.0,2.4])
            )
            self.h_electron_den_csv[i] = ROOT.TH2F(
                'h_electron_den_csv{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                'h_electron_den_csv{}{}{}'.format("_" + self.sample, "_" + cat, self.syst_suffix),
                9, ar.array('d',[20,30,50,70,100,140,200,300,600,1000]),
                5, ar.array('d',[0.,0.5,1.,1.5,2.0,2.4])
            )


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        prevdir = ROOT.gDirectory
        outputFile.cd()
        for i,cat in self.cats.items():
            self.h_electron_num[i].Write()
            self.h_muon_num[i].Write()
            self.h_electron_den[i].Write()
            self.h_muon_den[i].Write()

            self.h_electron_num_csv[i].Write()
            self.h_muon_num_csv[i].Write()
            self.h_electron_den_csv[i].Write()
            self.h_muon_den_csv[i].Write()

        prevdir.cd()

    def analyze(self, event):
        try:
            lep_category = getattr(event, "lep_category{}".format(self.syst_suffix))
            event_category = getattr(event, "event_category{}".format(self.syst_suffix))
            leading_lep_pt = getattr(event,"leading_lep_pt{}".format(self.syst_suffix))
            leading_lep_eta = getattr(event,"leading_lep_eta{}".format(self.syst_suffix))
            trailing_lep_pt = getattr(event,"trailing_lep_pt{}".format(self.syst_suffix))
            trailing_lep_eta = getattr(event,"trailing_lep_eta{}".format(self.syst_suffix))
            Zpt = getattr(event, "Zlep_cand_pt{}".format(self.syst_suffix))
        except:
            lep_category = getattr(event, "lep_category")
            event_category = getattr(event, "event_category")
            leading_lep_pt = getattr(event,"leading_lep_pt")
            leading_lep_eta = getattr(event,"leading_lep_eta")
            trailing_lep_pt = getattr(event,"trailing_lep_pt")
            trailing_lep_eta = getattr(event,"trailing_lep_eta")
            Zpt = getattr(event, "Zlep_cand_pt")
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
            print ' pozor'
        # only valid in the MC samples
#        print "Zpt: ",Zpt
#        print "channels: ",lep_category
#        print "l lepton: ",leading_lep_pt,leading_lep_eta
#        print "t lepton: ",trailing_lep_pt,trailing_lep_eta

        jets = list(Collection(event, "Jet"))
        for j in jets:
            flavour =  j.hadronFlavour
            tagged = -99
            taggedCSV = -99
            if j.pt < 20.0 or abs(j.eta) > 2.4:
                continue
            if not j.jetId :
                continue
            if j.btagDeepFlavB > self.btagDeepJet_id("loose"):
                tagged = 1
            else:
                tagged = 0
            if j.btagDeepB > self.btag_id("loose"):
                taggedCSV = 1
            else:
                taggedCSV = 0
            if ( (lep_category == 1) and self.passbut(event, cat="signal") ):
                if(flavour == 5):
                    self.h_muon_den[1].Fill(j.pt,abs(j.eta),weight)
                    self.h_muon_den_csv[1].Fill(j.pt,abs(j.eta),weight)
                    if(tagged == 1):
                        self.h_muon_num[1].Fill(j.pt,abs(j.eta),weight)
                    if(taggedCSV == 1):
                        self.h_muon_num_csv[1].Fill(j.pt,abs(j.eta),weight)
                if(flavour == 4):
                    self.h_muon_den[2].Fill(j.pt,abs(j.eta),weight)
                    self.h_muon_den_csv[2].Fill(j.pt,abs(j.eta),weight)
                    if(tagged == 1):
                        self.h_muon_num[2].Fill(j.pt,abs(j.eta),weight)
                    if(taggedCSV == 1):
                        self.h_muon_num_csv[2].Fill(j.pt,abs(j.eta),weight)
                if(flavour == 0):
                    self.h_muon_den[3].Fill(j.pt,abs(j.eta),weight)
                    self.h_muon_den_csv[3].Fill(j.pt,abs(j.eta),weight)
                    if(tagged == 1):
                        self.h_muon_num[3].Fill(j.pt,abs(j.eta),weight)
                    if(taggedCSV == 1):
                        self.h_muon_num_csv[3].Fill(j.pt,abs(j.eta),weight)

            if ( (lep_category == 2) and self.passbut(event, cat="signal") ):
                if(flavour == 5):
                    self.h_electron_den[1].Fill(j.pt,abs(j.eta),weight)
                    self.h_electron_den_csv[1].Fill(j.pt,abs(j.eta),weight)
                    if(tagged == 1):
                        self.h_electron_num[1].Fill(j.pt,abs(j.eta),weight)
                    if(taggedCSV == 1):
                        self.h_electron_num_csv[1].Fill(j.pt,abs(j.eta),weight)
                if(flavour == 4):
                    self.h_electron_den[2].Fill(j.pt,abs(j.eta),weight)
                    self.h_electron_den_csv[2].Fill(j.pt,abs(j.eta),weight)
                    if(tagged == 1):
                        self.h_electron_num[2].Fill(j.pt,abs(j.eta),weight)
                    if(taggedCSV == 1):
                        self.h_electron_num_csv[2].Fill(j.pt,abs(j.eta),weight)
                if(flavour == 0):
                    self.h_electron_den[3].Fill(j.pt,abs(j.eta),weight)
                    self.h_electron_den_csv[3].Fill(j.pt,abs(j.eta),weight)
                    if(tagged == 1):
                        self.h_electron_num[3].Fill(j.pt,abs(j.eta),weight)
                    if(taggedCSV == 1):
                        self.h_electron_num_csv[3].Fill(j.pt,abs(j.eta),weight)
        return True
