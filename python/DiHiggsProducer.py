import ROOT
import sys, os
import numpy as np
import math
from importlib import import_module
import itertools
from copy import deepcopy
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import PhysicsTools.NanoAODTools.postprocessing.tools as tk

ROOT.PyConfig.IgnoreCommandLineOptions = True


class DiHiggsProducer(Module):
    def __init__(self, isMC, era, do_syst = False, syst_var=''):
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

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("good_event{}".format(self.syst_suffix), "I")
        self.out.branch("met_pt{}".format(self.syst_suffix), "F")
        self.out.branch("met_phi{}".format(self.syst_suffix), "F")
        self.out.branch("ngood_leptons{}".format(self.syst_suffix), "I")
        self.out.branch("nextra_leptons{}".format(self.syst_suffix), "I")
        self.out.branch("lep_category{}".format(self.syst_suffix), "I") 
        # 1 = dimuon channel, 2 = dielectron channel
        self.out.branch("event_category{}".format(self.syst_suffix), "I") 
        # used for ABCD method
        # 1 = A OS+ISO
        # 2 = B OS+non-ISO
        # 3 = C SS+ISO
        # 4 = D SS+non-ISO
        self.out.branch("dR_l1l2{}".format(self.syst_suffix), "F")
        self.out.branch("event_category_30{}".format(self.syst_suffix), "I")
        self.out.branch("event_category_35{}".format(self.syst_suffix), "I")
        self.out.branch("event_category_40{}".format(self.syst_suffix), "I")

        self.out.branch("event_category_mini30{}".format(self.syst_suffix), "I")
        self.out.branch("event_category_mini35{}".format(self.syst_suffix), "I")

        self.out.branch("leading_Hbb_pt{}".format(self.syst_suffix), "F")
        self.out.branch("leading_Hbb_eta{}".format(self.syst_suffix), "F")
        self.out.branch("leading_Hbb_phi{}".format(self.syst_suffix), "F")
        self.out.branch("leading_Hbb_btag{}".format(self.syst_suffix), "F")


        self.out.branch("trailing_Hbb_pt{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_Hbb_eta{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_Hbb_phi{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_Hbb_btag{}".format(self.syst_suffix), "F")

        self.out.branch("leading_lep_pt{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_eta{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_phi{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_iso{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_iso_chg{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_pt{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_eta{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_phi{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_iso{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_iso_chg{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_flavor{}".format(self.syst_suffix), "I")
        self.out.branch("trailing_lep_flavor{}".format(self.syst_suffix), "I")

        self.out.branch("leading_jet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("leading_jet_eta{}".format(self.syst_suffix), "F")
        self.out.branch("leading_jet_phi{}".format(self.syst_suffix), "F")
        self.out.branch("leading_jet_qgl{}".format(self.syst_suffix), "F")

        self.out.branch("trailing_jet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_jet_eta{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_jet_phi{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_jet_qgl{}".format(self.syst_suffix), "F")

        self.out.branch("met_filter{}".format(self.syst_suffix), "I")

        #Adding in the Higgs boson candidate variables
        self.out.branch("Higgsbb_cand_pt{}".format(self.syst_suffix), "F")
        self.out.branch("Higgsbb_cand_eta{}".format(self.syst_suffix), "F")
        self.out.branch("Higgsbb_cand_phi{}".format(self.syst_suffix), "F")
        self.out.branch("Higgsbb_cand_mass{}".format(self.syst_suffix), "F")

        self.out.branch("HiggsZZ_cand_pt{}".format(self.syst_suffix), "F")
        self.out.branch("HiggsZZ_cand_eta{}".format(self.syst_suffix), "F")
        self.out.branch("HiggsZZ_cand_phi{}".format(self.syst_suffix), "F")
        self.out.branch("HiggsZZ_cand_mass{}".format(self.syst_suffix), "F")

        self.out.branch("Zlep_cand_pt{}".format(self.syst_suffix), "F")
        self.out.branch("Zlep_cand_eta{}".format(self.syst_suffix), "F")
        self.out.branch("Zlep_cand_phi{}".format(self.syst_suffix), "F")
        self.out.branch("Zlep_cand_mass{}".format(self.syst_suffix), "F")

        self.out.branch("Zjet_cand_pt{}".format(self.syst_suffix), "F")
        self.out.branch("Zjet_cand_eta{}".format(self.syst_suffix), "F")
        self.out.branch("Zjet_cand_phi{}".format(self.syst_suffix), "F")
        self.out.branch("Zjet_cand_mass{}".format(self.syst_suffix), "F")

        self.out.branch("HH_cand_pt{}".format(self.syst_suffix), "F")
        self.out.branch("HH_cand_eta{}".format(self.syst_suffix), "F")
        self.out.branch("HH_cand_phi{}".format(self.syst_suffix), "F")
        self.out.branch("HH_cand_mass{}".format(self.syst_suffix), "F")

        self.out.branch("jetHT{}".format(self.syst_suffix), "F")
        self.out.branch("ngood_jets{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_jets_noHbb{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_bjets{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_bjetsM{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_bjetsT{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_bjetsC{}".format(self.syst_suffix), "I")
        self.out.branch("lead_jet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("lead_bjet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("delta_phi_j_met{}".format(self.syst_suffix), "F")
        self.out.branch("leptoncharge{}".format(self.syst_suffix), "I")

        self.out.branch("nhad_taus{}".format(self.syst_suffix), "I")
        self.out.branch("lead_tau_pt{}".format(self.syst_suffix), "F")

        if self.isMC:
            self.out.branch("w_muon_SF{}".format(self.syst_suffix), "F")
            self.out.branch("w_muon_SFUp{}".format(self.syst_suffix), "F")
            self.out.branch("w_muon_SFDown{}".format(self.syst_suffix), "F")
            self.out.branch("w_electron_SF{}".format(self.syst_suffix), "F")
            self.out.branch("w_electron_SFUp{}".format(self.syst_suffix), "F")
            self.out.branch("w_electron_SFDown{}".format(self.syst_suffix), "F")
        if self.isMC:
            self.out.branch("w_btag_SF{}".format(self.syst_suffix),"F")
        if self.isMC and len(self.syst_suffix)==0:
            self.out.branch("w_btag_SF_sys_up_hf".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_up_lf".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_up_cferr1".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_up_cferr2".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_up_hfstats1".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_up_hfstats2".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_up_lfstats1".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_up_lfstats2".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_down_hf".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_down_lf".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_down_cferr1".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_down_cferr2".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_down_hfstats1".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_down_hfstats2".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_down_lfstats1".format(self.syst_suffix),"F")
            self.out.branch("w_btag_SF_sys_down_lfstats2".format(self.syst_suffix),"F")
        if self.isMC and len(self.syst_suffix)==0:
            self.out.branch("met_ptUnclustEnUp{}".format(self.syst_suffix), "F")
            self.out.branch("met_phiUnclustEnUp{}".format(self.syst_suffix), "F")
            self.out.branch("met_ptUnclustEnDown{}".format(self.syst_suffix), "F")
            self.out.branch("met_phiUnclustEnDown{}".format(self.syst_suffix), "F")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def electron_id(self, electron, wp):
        pass_id = 0
        if ( wp == "80"):
            try:
                pass_id = electron.mvaFall17V2Iso_WP80
            except:
                try:
                    pass_id = electron.mvaFall17V1Iso_WP80
                except:
                    try:
                        pass_id = electron.mvaFall17Iso_WP80
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "

            return pass_id
        elif ( wp == "90"):
            try:
                pass_id = electron.mvaFall17V2Iso_WP90
            except:
                try:
                    pass_id = electron.mvaFall17V1Iso_WP90
                except:
                    try:
                        pass_id = electron.mvaFall17Iso_WP90
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "

            return pass_id
        elif ( wp == "WPL"):
            try:
                pass_id = electron.mvaFall17V2Iso_WPL
            except:
                try:
                    pass_id = electron.mvaFall17V1Iso_WPL
                except:
                    try:
                        pass_id = electron.mvaFall17Iso_WPL
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "

            return pass_id

    def electron_id_noiso(self, electron, wp):
        pass_id = 0
        if ( wp == "80"):
            try:
                pass_id = electron.mvaFall17V2noIso_WP80
            except:
                try:
                    pass_id = electron.mvaFall17V1noIso_WP80
                except:
                    try:
                        pass_id = electron.mvaFall17noIso_WP80
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "

            return pass_id
        elif (wp == "90"):
            try:
                pass_id = electron.mvaFall17V2noIso_WP90
            except:
                try:
                    pass_id = electron.mvaFall17V1noIso_WP90
                except:
                    try:
                        pass_id = electron.mvaFall17noIso_WP90
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "

            return pass_id
        elif ( wp == "WPL"):
 #           print "wp loose"
            try:          
                pass_id = electron.mvaFall17V2noIso_WPL
            except:
                try:
                    pass_id = electron.mvaFall17V1noIso_WPL
                except:
                    try:
                        pass_id = electron.mvaFall17noIso_WPL
                    except ValueError:
                        print "[error] not mvaFall17 electron id found ... "
            return pass_id


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

    def met_filter(self, flag, filter_mask=True):
        return filter_mask and (
              (flag.HBHENoiseFilter)
           and (flag.HBHENoiseIsoFilter)
           and (flag.EcalDeadCellTriggerPrimitiveFilter)
           and (flag.goodVertices)
           and (flag.eeBadScFilter)
           and (flag.globalTightHalo2016Filter)
#           and (flag.BadChargedCandidateFilter)
           and (flag.BadPFMuonFilter)
        )

    def duplicate_removal(self):
        """
        For data, same event could come from different datasets
        FIXME: need to be implemented check the source from
        the old MonoZ code
        https://github.com/NEUAnalyses/monoZ_Analysis/blob/master/src/MonoZSelector.cc#L463
        """
        pass

    def analyze(self, event):
        """
        process event, return True (go to next module)
        or False (fail, go to next event)
        """
        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))
        jets = list(Collection(event, "Jet"))
        taus = list(Collection(event, "Tau"))
        flag = Object(event, "Flag")
        met = Object(event, "MET")

        # in case of systematic take the shifted values are default
        # For the central values, need to include jetMetTool all the time
        # Jet systematics
        if self.syst_var == "":
            syst_var = "nom"
        else:
            syst_var = self.syst_var
        # checking something
#        print 'syst var is: ',syst_var
        try:
            var_jet_pts = getattr(event,  "Jet_pt_{}".format(syst_var), None)
            if var_jet_pts:
                for i,jet in enumerate(jets):
                    jet.pt = var_jet_pts[i]
            else:
                print 'WARNING: jet pts with variation {}'
                'not available, using the nominal value'.format(syst_var)
        except:
 #           print 'WARNING: jet pts with variation {} not available, using the nominal value'.format(syst_var)
            var_jet_pts = getattr(event,  "Jet_pt_nom", None)
            for i,jet in enumerate(jets):
                jet.pt = var_jet_pts[i]

        ## btag SF factor
        if self.syst_var == "":
            syst_var_btag = ""
        else:
  #          print "sys name origin: ",self.syst_var
            syst_var_tmp = self.syst_var
            if "Up" in syst_var_tmp:
                syst_var_tmp1 = syst_var_tmp.replace("Up", "")
                syst_var_tmp2 = syst_var_tmp1.replace("sys", "")
                syst_var_btag = "_up_jes"+syst_var_tmp2
 #               print "test sys name: ",sys_var_btag
            if "Down" in syst_var_tmp:
                syst_var_tmp1 = syst_var_tmp.replace("Down", "")
                syst_var_tmp2 = syst_var_tmp1.replace("sys", "")
                syst_var_btag = "_down_jes"+syst_var_tmp2
#                print "test sys name: ",sys_var_btag

#jet.btagSF_deepjet_shape
        try:
            var_jet_b = getattr(event,  "Jet_btagSF_deepjet_shape{}".format(syst_var_btag), None)
            if var_jet_b:
                for i,jet in enumerate(jets):
                    jet.btagSF_deepjet_shape = var_jet_b[i]
            else:
                print 'WARNING 0: jet btagSF with variation {} not available, using the nominal value'.format(syst_var_btag)
        except:
#            print 'WARNING 1: jet btagSF with variation {} not available, using the nominal value'.format(syst_var_btag)
            var_jet_b = getattr(event,  "Jet_btagSF_deepjet_shape", None)
            for i,jet in enumerate(jets):
                jet.btagSF_deepjet_shape = var_jet_b[i]

        if self.era == "2017":
            try:
                var_met_pt  = getattr(event,  "METFixEE2017_pt_{}".format(syst_var), None)
                var_met_phi = getattr(event, "METFixEE2017_phi_{}".format(syst_var), None)
                if var_met_pt:
                    met.pt = var_met_pt
                else:
                    print 'WARNING: MET pt with variation '
                    '{} not available, using the nominal value'.format(syst_var)
                if var_met_phi:
                    met.phi = var_met_phi
                else:
                    print 'WARNING: MET phi with variation {}'
                    'not available, using the nominal value'.format(syst_var)
            except:
                var_met_pt  = getattr(event,  "METFixEE2017_pt", None)
                var_met_phi = getattr(event, "METFixEE2017_phi", None)
                if var_met_pt:
                    met.pt = var_met_pt
                if var_met_phi:
                    met.phi = var_met_phi
        else:
            try:
                var_met_pt  = getattr(event,  "MET_pt_{}".format(syst_var), None)
                var_met_phi = getattr(event, "MET_phi_{}".format(syst_var), None)
                if var_met_pt:
                    met.pt = var_met_pt
                else:
                    print 'WARNING: MET pt with variation '
                    '{} not available, using the nominal value'.format(syst_var)
                if var_met_phi:
                    met.phi = var_met_phi
                else:
                    print 'WARNING: MET phi with variation {}'
                    'not available, using the nominal value'.format(syst_var)
            except:
                var_met_pt  = getattr(event,  "MET_T1_pt", None)
                var_met_phi = getattr(event, "MET_T1_phi", None)
                if var_met_pt:
                    met.pt = var_met_pt
                if var_met_phi:
                    met.phi = var_met_phi

        met_p4 = ROOT.TLorentzVector()
        met_p4.SetPtEtaPhiM(met.pt,0.0,met.phi, 0.0)

        if self.isMC and len(self.syst_suffix)==0:
            var_met_pt_Up  = 0
            var_met_phi_Up = 0
            var_met_pt_Down  = 0
            var_met_phi_Down = 0
            try:
                if self.era == "2017":
                    var_met_pt_Up = getattr(event,  "METFixEE2017_pt_unclustEnUp", None)
                else:
                    var_met_pt_Up = getattr(event,  "MET_pt_unclustEnUp", None)
            except:
                var_met_pt_Up = -99
            try:
                if self.era == "2017":
                    var_met_phi_Up = getattr(event, "METFixEE2017_phi_unclustEnUp", None)
                else:
                    var_met_phi_Up = getattr(event, "MET_phi_unclustEnUp", None)
            except:
                var_met_phi_Up = -99
            try:
                if self.era == "2017":
                    var_met_pt_Down = getattr(event,  "METFixEE2017_pt_unclustEnDown", None)
                else:
                    var_met_pt_Down = getattr(event,  "MET_pt_unclustEnDown", None)
            except:
                var_met_pt_Down = -99
            try:
                if self.era == "2017":
                    var_met_phi_Down = getattr(event, "METFixEE2017_phi_unclustEnDown", None)
                else:
                    var_met_phi_Down = getattr(event, "MET_phi_unclustEnDown", None)
            except:
                var_met_phi_Down = -99

            self.out.fillBranch("met_ptUnclustEnUp{}".format(self.syst_suffix),var_met_pt_Up)
            self.out.fillBranch("met_phiUnclustEnUp{}".format(self.syst_suffix), var_met_phi_Up)
            self.out.fillBranch("met_ptUnclustEnDown{}".format(self.syst_suffix),var_met_pt_Down)
            self.out.fillBranch("met_phiUnclustEnDown{}".format(self.syst_suffix),var_met_phi_Down)

        # Electrons Energy
        if "ElectronEn" in self.syst_var:
            (met_px, met_py) = ( met.pt*np.cos(met.phi), met.pt*np.sin(met.phi) )
            if "Up" in self.syst_var:
                for i, elec in enumerate(electrons):
                    met_px = met_px + (elec.energyErr)*np.cos(elec.phi)/math.cosh(elec.eta)
                    met_py = met_py + (elec.energyErr)*np.sin(elec.phi)/math.cosh(elec.eta)
                    elec.pt = elec.pt + elec.energyErr/math.cosh(elec.eta)
            else:
                for i, elec in enumerate(electrons):
                    met_px = met_px - (elec.energyErr)*np.cos(elec.phi)/math.cosh(elec.eta)
                    met_py = met_py - (elec.energyErr)*np.sin(elec.phi)/math.cosh(elec.eta)
                    elec.pt = elec.pt - elec.energyErr/math.cosh(elec.eta)
            met.pt  = math.sqrt(met_px**2 + met_py**2)
            met.phi = math.atan2(met_py, met_px)

        # Muons Energy
        #if self.isMC:
        muons_pts = getattr(event, "Muon_corrected_pt")
        for i, muon in enumerate(muons):
            muon.pt = muons_pts[i]

        if "MuonEn" in self.syst_var:
            (met_px, met_py) = ( met.pt*np.cos(met.phi), met.pt*np.sin(met.phi) )
            if "Up" in self.syst_var:
                muons_pts = getattr(event, "Muon_correctedUp_pt")
                for i, muon in enumerate(muons):
                    met_px = met_px - (muons_pts[i] - muon.pt)*np.cos(muon.phi)
                    met_py = met_py - (muons_pts[i] - muon.pt)*np.sin(muon.phi)
                    muon.pt = muons_pts[i]
            else:
                muons_pts = getattr(event, "Muon_correctedDown_pt")
                for i, muon in enumerate(muons):
                    met_px =met_px - (muons_pts[i] - muon.pt)*np.cos(muon.phi)
                    met_py =met_py - (muons_pts[i] - muon.pt)*np.sin(muon.phi)
                    muon.pt = muons_pts[i]
            met.pt  = math.sqrt(met_px**2 + met_py**2)
            met.phi = math.atan2(met_py, met_px)
            
        # filling and contructing the event categorisation
        self.out.fillBranch("met_pt{}".format(self.syst_suffix), met.pt)
        self.out.fillBranch("met_phi{}".format(self.syst_suffix), met.phi)

        pass_met_filter = self.met_filter(flag, True)
        self.out.fillBranch("met_filter{}".format(self.syst_suffix), pass_met_filter)

        # count electrons and muons
        good_leptons = []
        good_muons = []
        good_electrons = []
        lep_category = -1
        event_category = -1

        muons.sort(key=lambda muon: muon.pt, reverse=True)
        electrons.sort(key=lambda el: el.pt, reverse=True)
        # Choose loose/medium-quality e/mu for event categorization
        for idx,mu in enumerate(muons):
            pass_ips = abs(mu.dxy) < 0.02 and abs(mu.dz) < 0.1
            pass_fid = abs(mu.eta) < 2.4 and mu.pt >= (20 if idx==0 else 10)
            pass_ids = mu.mediumId
            if pass_fid and pass_ids and pass_ips:
                good_muons.append(mu)
        for idy,el in enumerate(electrons):
            pass_ips_dxy = abs(el.dxy) < 0.10 #if abs(el.eta) < 1.442 else 0.10)
            pass_ips_dz =  abs(el.dz ) < 0.20 #if abs(el.eta) < 1.442 else 0.20)
            pass_ids = self.electron_id_noiso(el,"WPL")
            pass_fid = el.pt >= (25 if idy==0 else 15) and abs(el.eta) <= 2.5
            # changing to MVA based ID :
            if pass_fid and pass_ids and pass_ips_dz and pass_ips_dxy:
                good_electrons.append(el)

        # let sort the muons in pt
        good_muons.sort(key=lambda x: x.pt, reverse=True)
        good_electrons.sort(key=lambda x: x.pt, reverse=True)

        # Find any remaining e/mu that pass looser selection
        extra_leptons = []
        for mu in muons:
            isoLep   = mu.pfRelIso04_all
            pass_ids = mu.looseId and isoLep <= 0.25             
            pass_ips = abs(mu.dxy) < 0.02 and abs(mu.dz) < 0.1
            pass_fid = abs(mu.eta) < 2.4 and mu.pt >= 10
            if tk.closest(mu, good_muons)[1] < 0.01:
                continue
            if pass_fid and pass_ids and pass_ips:
                extra_leptons.append(mu)

        for el in electrons:
            pass_fid = abs(el.eta) < 2.5 and el.pt >= 10
            if tk.closest(el, good_electrons)[1] < 0.01:
                continue
            pass_ips_dxy = abs(el.dxy) < (0.05 if abs(el.eta) < 1.442 else 0.10)
            pass_ips_dz =  abs(el.dz ) < (0.10 if abs(el.eta) < 1.442 else 0.20)
            pass_ids = el.miniPFRelIso_all < 0.4
            if pass_fid and pass_ids and pass_ips_dxy and pass_ips_dz and self.electron_id_noiso(el, "WPL"):
                extra_leptons.append(el)



        # find categories
        z_candidate = []
        zcand_p4 = ROOT.TLorentzVector()
        emulated_met = ROOT.TLorentzVector()
        all_lepton_p4 = ROOT.TLorentzVector()
        rem_lepton_p4 = ROOT.TLorentzVector()

        good_leptons = good_electrons + good_muons
        good_leptons.sort(key=lambda x: x.pt, reverse=True)

        _lead_lep_pt = good_leptons[0].pt if len(good_leptons) else 0.0
        _lead_lep_eta = good_leptons[0].eta if len(good_leptons) else 0.0
        _trail_lep_pt = good_leptons[1].pt if len(good_leptons) >= 2 else 0.0
        _trail_lep_eta = good_leptons[1].eta if len(good_leptons) >= 2 else 0.0
        _lead_lep_pdgId = good_leptons[0].pdgId if len(good_leptons)  else 0.0
        _trail_lep_pdgId = good_leptons[1].pdgId if len(good_leptons) >= 2 else 0.0


	_leading_lep_flavor = 0
	if len(good_muons) and len(good_electrons):
		if good_muons[0].pt > good_electrons[0].pt: _leading_lep_flavor = 1

        self.out.fillBranch("leading_lep_pt{}".format(self.syst_suffix), _lead_lep_pt)
        self.out.fillBranch("leading_lep_eta{}".format(self.syst_suffix), _lead_lep_eta)
        self.out.fillBranch("trailing_lep_pt{}".format(self.syst_suffix), _trail_lep_pt)
        self.out.fillBranch("trailing_lep_eta{}".format(self.syst_suffix), _trail_lep_eta)
        self.out.fillBranch("leading_lep_flavor{}".format(self.syst_suffix), _leading_lep_flavor)

        ngood_leptons = len(good_leptons)
        nextra_leptons = len(extra_leptons)
        ngood_muons = len(good_muons)
        ngood_electrons = len(good_electrons)
#        print "ngood_electrons: ",ngood_electrons
        if False:
            print "number of leptons [all, good, extra]: ", ngood_leptons, " : ", nextra_leptons
            print "        CBId electrons : ", [e.cutBased for e in good_electrons]
            print "        WP90 electrons : ", [e.mvaFall17Iso_WP90 for e in good_electrons]
            print "             muons     : ", [e.tightId for e in good_muons]
            print "        lepton pts     : ", [e.pt for e in good_leptons]

        self.out.fillBranch("ngood_leptons{}".format(self.syst_suffix), ngood_leptons)
        self.out.fillBranch("nextra_leptons{}".format(self.syst_suffix), nextra_leptons)

        # Leptons efficiency/Trigger/Isolation Scale factors
        # These are applied only of the first 2 leading leptons
        if self.isMC:
            w_muon_SF     = w_electron_SF     = 1.0
            w_muon_SFUp   = w_electron_SFUp   = 1.0
            w_muon_SFDown = w_electron_SFDown = 1.0
            if ngood_leptons >= 2:
                if abs(good_leptons[0].pdgId) == 11:
                    w_electron_SF     *=  good_leptons[0].SF
                    w_electron_SFUp   *= (good_leptons[0].SF + good_leptons[0].SFErr)
                    w_electron_SFDown *= (good_leptons[0].SF - good_leptons[0].SFErr)
                if abs(good_leptons[1].pdgId) == 11:
                    w_electron_SF     *=  good_leptons[1].SF
                    w_electron_SFUp   *= (good_leptons[1].SF + good_leptons[1].SFErr)
                    w_electron_SFDown *= (good_leptons[1].SF - good_leptons[1].SFErr)
                if abs(good_leptons[0].pdgId) == 13:
                    w_muon_SF     *=  good_leptons[0].SF
                    w_muon_SFUp   *= (good_leptons[0].SF + good_leptons[0].SFErr)
                    w_muon_SFDown *= (good_leptons[0].SF - good_leptons[0].SFErr)
                if abs(good_leptons[1].pdgId) == 13:
                    w_muon_SF     *=  good_leptons[1].SF
                    w_muon_SFUp   *= (good_leptons[1].SF + good_leptons[1].SFErr)
                    w_muon_SFDown *= (good_leptons[1].SF - good_leptons[1].SFErr)
            self.out.fillBranch("w_muon_SF"        , w_muon_SF        )
            self.out.fillBranch("w_muon_SFUp"      , w_muon_SFUp      )
            self.out.fillBranch("w_muon_SFDown"    , w_muon_SFDown    )
            self.out.fillBranch("w_electron_SF"    , w_electron_SF    )
            self.out.fillBranch("w_electron_SFUp"  , w_electron_SFUp  )
            self.out.fillBranch("w_electron_SFDown", w_electron_SFDown)

        dR_l1l2 = -99
        leptoncharge_test = 0
        lep_category = 0
        event_category = 0
        
        event_category_35 = 0
        event_category_30 = 0
        event_category_40 = 0
        event_category_mini30 = 0
        event_category_mini35 = 0
        
        if ngood_leptons < 2:
            lep_category = -1
            event_category = -1
            
            event_category_30 = -1
            event_category_35 = -1
            event_category_40 = -1
            event_category_mini30 = -1
            event_category_mini35 = -1
            
        if ngood_muons == 2 and ngood_electrons == 0 and nextra_leptons == 0:
            lep_category = 1 #dimuon channel
            isoLep0   = good_muons[0].miniPFRelIso_all
            isoLep1   = good_muons[1].miniPFRelIso_all
            dR_l1l2 = tk.deltaR(good_muons[0].eta, good_muons[0].phi,good_muons[1].eta, good_muons[1].phi,)
            
            normalisoLep0   = good_muons[0].pfRelIso04_all
            normalisoLep1   = good_muons[1].pfRelIso04_all
            '''
            miniisoLep0   = good_muons[0].miniPFRelIso_all
            miniisoLep1   = good_muons[1].miniPFRelIso_all
            '''
            if (good_muons[0].pdgId * good_muons[1].pdgId == -13*13): # opposite sign: A or B region
                if (isoLep0 < 0.4 and isoLep1 < 0.4):
                    event_category = 1
                else:
                    event_category = 2
                
                if (isoLep0 < 0.30 and isoLep1 < 0.30):
                    event_category_mini30 = 1
                else:
                    event_category_mini30 = 2
                if (isoLep0 < 0.35 and isoLep1 < 0.35):
                    event_category_mini35 = 1
                else:
                    event_category_mini35 = 2
               
                if (normalisoLep0 < 0.40 and normalisoLep1 < 0.40):
                    event_category_40 = 1
                else:
                    event_category_40 = 2
                if (normalisoLep0 < 0.30 and normalisoLep1 < 0.30):
                    event_category_30 = 1
                else:
                    event_category_30 = 2
                if (normalisoLep0 < 0.35 and normalisoLep1 < 0.35):
                    event_category_35 = 1
                else:
                    event_category_35 = 2
                
            else: #same sign: C or D region
                if (isoLep0 < 0.4 and isoLep1 < 0.4):
                    event_category = 3
                else:
                    event_category = 4
                
                if (isoLep0 < 0.30 and isoLep1 < 0.30):
                    event_category_mini30 = 3
                else:
                    event_category_mini30 = 4
                if (isoLep0 < 0.35 and isoLep1 < 0.35):
                    event_category_mini35 = 3
                else:
                    event_category_mini35 = 4
                if (normalisoLep0 < 0.40 and normalisoLep1 < 0.40):
                    event_category_40 = 3
                else:
                    event_category_40 = 4
                if (normalisoLep0 < 0.30 and normalisoLep1 < 0.30):
                    event_category_30 = 3
                else:
                    event_category_30 = 4
                if (normalisoLep0 < 0.35 and normalisoLep1 < 0.35):
                    event_category_35 = 3
                else:
                    event_category_35 = 4
            
            if tk.deltaR(good_muons[0].eta, good_muons[0].phi,good_muons[1].eta, good_muons[1].phi,) < 0.3:
                event_category +=10 
                event_category_30 +=10
                event_category_35 +=10
                event_category_40 +=10
                event_category_mini30 +=10
                event_category_mini35 +=10
                
        elif ngood_electrons == 2 and ngood_muons == 0 and nextra_leptons == 0:
            lep_category = 2 #dielectron channel
            dR_l1l2 = tk.deltaR(good_electrons[0].eta, good_electrons[0].phi,good_electrons[1].eta, good_electrons[1].phi,)
            if (good_electrons[0].pdgId * good_electrons[1].pdgId == -11*11)==(good_electrons[0].charge * good_electrons[1].charge < 0):
                leptoncharge_test = 0
            else:
                leptoncharge_test = 1

            isoLep0   = good_electrons[0].miniPFRelIso_all
            isoLep1   = good_electrons[1].miniPFRelIso_all
            normalisoLep0   = good_electrons[0].pfRelIso03_all
            normalisoLep1   = good_electrons[1].pfRelIso03_all
            '''
            miniisoLep0   = good_electrons[0].miniPFRelIso_all
            miniisoLep1   = good_electrons[1].miniPFRelIso_all
            '''
            if (good_electrons[0].charge * good_electrons[1].charge < 0): # opposite sign: A or B region
                if (isoLep0 < 0.4 and isoLep1 < 0.4):
                    event_category = 1
                else:
                    event_category = 2
                if (isoLep0 < 0.30 and isoLep1 < 0.30):
                    event_category_mini30 = 1
                else:
                    event_category_mini30 = 2
                if (isoLep0 < 0.35 and isoLep1 < 0.35):
                    event_category_mini35 = 1
                else:
                    event_category_mini35 = 2
                if (normalisoLep0 < 0.40 and normalisoLep1 < 0.40):
                    event_category_40 = 1
                else:
                    event_category_40 = 2
                if (normalisoLep0 < 0.30 and normalisoLep1 < 0.30):
                    event_category_30 = 1
                else:
                    event_category_30 = 2
                if (normalisoLep0 < 0.35 and normalisoLep1 < 0.35):
                    event_category_35 = 1
                else:
                    event_category_35 = 2
            else: #same sign: C or D region
                if (isoLep0 < 0.4 and isoLep1 < 0.4):
                    event_category = 3
                else:
                    event_category = 4
                if (isoLep0 < 0.30 and isoLep1 < 0.30):
                    event_category_mini30 = 3
                else:
                    event_category_mini30 = 4
                if (isoLep0 < 0.35 and isoLep1 < 0.35):
                    event_category_mini35 = 3
                else:
                    event_category_mini35 = 4
                if (normalisoLep0 < 0.40 and normalisoLep1 < 0.40):
                    event_category_40 = 3
                else:
                    event_category_40 = 4
                if (normalisoLep0 < 0.30 and normalisoLep1 < 0.30):
                    event_category_30 = 3
                else:
                    event_category_30 = 4
                if (normalisoLep0 < 0.35 and normalisoLep1 < 0.35):
                    event_category_35 = 3
                else:
                    event_category_35 = 4

            if tk.deltaR(good_electrons[0].eta, good_electrons[0].phi,good_electrons[1].eta, good_electrons[1].phi,) < 0.3:
                event_category +=10 
                
                event_category_30 +=10
                event_category_35 +=10
                event_category_40 +=10
                event_category_mini30 +=10
                event_category_mini35 +=10
                
        self.out.fillBranch("leptoncharge{}".format(self.syst_suffix), leptoncharge_test)
        self.out.fillBranch("ngood_leptons{}".format(self.syst_suffix), ngood_leptons)
        self.out.fillBranch("nextra_leptons{}".format(self.syst_suffix), nextra_leptons)
        self.out.fillBranch("lep_category{}".format(self.syst_suffix), lep_category)
        self.out.fillBranch("event_category{}".format(self.syst_suffix), event_category)
        self.out.fillBranch("dR_l1l2{}".format(self.syst_suffix), dR_l1l2)
        
        self.out.fillBranch("event_category_30{}".format(self.syst_suffix),event_category_30 )
        self.out.fillBranch("event_category_35{}".format(self.syst_suffix),event_category_35 )
        self.out.fillBranch("event_category_40{}".format(self.syst_suffix), event_category_40)
        self.out.fillBranch("event_category_mini30{}".format(self.syst_suffix),event_category_mini30 )
        self.out.fillBranch("event_category_mini35{}".format(self.syst_suffix),event_category_mini35 )
        

        # process jet
        good_jets  = []
        good_jetsNoPULoose = []
        good_jetsNoPUMedium = []
        good_jetsNoPUTight = []
        good_bjets = []
        good_bjetsM = []
        good_bjetsT = []
        good_bjetsC = []
        btagSF = 1

        if self.isMC and len(self.syst_suffix)==0:
            btagSF_up_hf = 1
            btagSF_up_lf = 1
            btagSF_up_cferr1 = 1
            btagSF_up_cferr2 = 1
            btagSF_up_hfstats1 = 1
            btagSF_up_hfstats2 = 1
            btagSF_up_lfstats1 = 1
            btagSF_up_lfstats2 = 1
            btagSF_down_hf = 1
            btagSF_down_lf = 1
            btagSF_down_cferr1 = 1
            btagSF_down_cferr2 = 1
            btagSF_down_hfstats1 = 1
            btagSF_down_hfstats2 = 1
            btagSF_down_lfstats1 = 1
            btagSF_down_lfstats2 = 1

        jetHT = 0
        for jet in jets:
            if jet.btagDeepFlavB > self.btagDeepJet_id("loose"):
                jet.pt = jet.pt*jet.bRegCorr
            if jet.pt < 30.0 or abs(jet.eta) > 2.4:
                continue
            if jet.pt < 50.0 and jet.puId < 4:
                continue
            if not jet.jetId :
                continue
            if tk.closest(jet, good_leptons)[1] < 0.4:
                continue
            good_jets.append(jet)
            jetHT += jet.pt
            #https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration
            if self.isMC:
                btagSF *=jet.btagSF_deepjet_shape
                if self.isMC and len(self.syst_suffix)==0:
                    btagSF_up_hf *=jet.btagSF_deepjet_shape_up_hf
                    btagSF_up_lf *=jet.btagSF_deepjet_shape_up_lf
                    btagSF_up_cferr1 *=jet.btagSF_deepjet_shape_up_cferr1
                    btagSF_up_cferr2 *=jet.btagSF_deepjet_shape_up_cferr2
                    btagSF_up_hfstats1 *=jet.btagSF_deepjet_shape_up_hfstats1
                    btagSF_up_hfstats2 *=jet.btagSF_deepjet_shape_up_hfstats2
                    btagSF_up_lfstats1 *=jet.btagSF_deepjet_shape_up_lfstats1
                    btagSF_up_lfstats2 *=jet.btagSF_deepjet_shape_up_lfstats2
                    btagSF_down_hf *=jet.btagSF_deepjet_shape_down_hf
                    btagSF_down_lf *=jet.btagSF_deepjet_shape_down_lf
                    btagSF_down_cferr1 *=jet.btagSF_deepjet_shape_down_cferr1
                    btagSF_down_cferr2 *=jet.btagSF_deepjet_shape_down_cferr2
                    btagSF_down_hfstats1 *=jet.btagSF_deepjet_shape_down_hfstats1
                    btagSF_down_hfstats2 *=jet.btagSF_deepjet_shape_down_hfstats2
                    btagSF_down_lfstats1 *=jet.btagSF_deepjet_shape_down_lfstats1
                    btagSF_down_lfstats2 *=jet.btagSF_deepjet_shape_down_lfstats2
 
            # Count b-tag with loose WP DeepCSV
            # ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
            if jet.btagDeepFlavB > self.btagDeepJet_id("loose"):
#            if jet.btagDeepB > self.btag_id("loose"):
                good_bjets.append(jet)
            if jet.btagDeepFlavB > self.btagDeepJet_id("medium"):
                good_bjetsM.append(jet)
            if jet.btagDeepFlavB > self.btagDeepJet_id("tight"):
                good_bjetsT.append(jet)
            if jet.btagDeepB > self.btag_id("loose"):
                good_bjetsC.append(jet)

        good_jets.sort(key=lambda jet: jet.pt, reverse=True)
        good_bjets.sort(key=lambda jet: jet.btagDeepFlavB, reverse=True)
        if self.isMC:
            self.out.fillBranch("w_btag_SF{}".format(self.syst_suffix),btagSF)

        if self.isMC and len(self.syst_suffix)==0:
            self.out.fillBranch("w_btag_SF_sys_up_hf".format(self.syst_suffix),btagSF_up_hf)
            self.out.fillBranch("w_btag_SF_sys_up_lf".format(self.syst_suffix),btagSF_up_lf)
            self.out.fillBranch("w_btag_SF_sys_up_cferr1".format(self.syst_suffix),btagSF_up_cferr1)
            self.out.fillBranch("w_btag_SF_sys_up_cferr2".format(self.syst_suffix),btagSF_up_cferr2)
            self.out.fillBranch("w_btag_SF_sys_up_hfstats1".format(self.syst_suffix),btagSF_up_hfstats1)
            self.out.fillBranch("w_btag_SF_sys_up_hfstats2".format(self.syst_suffix),btagSF_up_hfstats2)
            self.out.fillBranch("w_btag_SF_sys_up_lfstats1".format(self.syst_suffix),btagSF_up_lfstats1)
            self.out.fillBranch("w_btag_SF_sys_up_lfstats2".format(self.syst_suffix),btagSF_up_lfstats2)
            self.out.fillBranch("w_btag_SF_sys_down_hf".format(self.syst_suffix),btagSF_down_hf)
            self.out.fillBranch("w_btag_SF_sys_down_lf".format(self.syst_suffix),btagSF_down_lf)
            self.out.fillBranch("w_btag_SF_sys_down_cferr1".format(self.syst_suffix),btagSF_down_cferr1)
            self.out.fillBranch("w_btag_SF_sys_down_cferr2".format(self.syst_suffix),btagSF_down_cferr2)
            self.out.fillBranch("w_btag_SF_sys_down_hfstats1".format(self.syst_suffix),btagSF_down_hfstats1)
            self.out.fillBranch("w_btag_SF_sys_down_hfstats2".format(self.syst_suffix),btagSF_down_hfstats2)
            self.out.fillBranch("w_btag_SF_sys_down_lfstats1".format(self.syst_suffix),btagSF_down_lfstats1)
            self.out.fillBranch("w_btag_SF_sys_down_lfstats2".format(self.syst_suffix),btagSF_down_lfstats2)

       #We will remove jets later so better count them now
        num_jets = len(good_jets)

        #Set up for the Higgs candicates!
        Higgsbb_cand_p4 = ROOT.TLorentzVector()
        HiggsZZ_cand_p4 = ROOT.TLorentzVector()
        HiggsZjet_cand_p4 = ROOT.TLorentzVector()
        HiggsZlep_cand_p4 = ROOT.TLorentzVector()
        HH_cand_p4 = ROOT.TLorentzVector()
        Higgs_cand_0 = ROOT.TLorentzVector()
        Zjet_cand_0 = ROOT.TLorentzVector()
        Zlep_cand_0 = ROOT.TLorentzVector()
        Higgs_cand_1 = ROOT.TLorentzVector()
        Higgsbb_candidate = []
        HiggsZjet_candidate = []
        HiggsZlep_candidate = []
        if len(good_leptons) == 2:
            HiggsZlep_cand_p4 = good_leptons[0].p4() + good_leptons[1].p4()
            HiggsZlep_candidate = [good_leptons[0],good_leptons[1]]
            self.out.fillBranch("Zlep_cand_pt{}".format(self.syst_suffix), HiggsZlep_cand_p4.Pt())
            self.out.fillBranch("Zlep_cand_eta{}".format(self.syst_suffix), HiggsZlep_cand_p4.Eta())
            self.out.fillBranch("Zlep_cand_phi{}".format(self.syst_suffix), HiggsZlep_cand_p4.Phi())
            self.out.fillBranch("Zlep_cand_mass{}".format(self.syst_suffix), HiggsZlep_cand_p4.M())

            HiggsZlep_candidate_sorted = sorted(HiggsZlep_candidate, key=lambda x: x.pt, reverse=True)
            try:
                _lead_lep_pt  = HiggsZlep_candidate_sorted[0].pt
                _lead_lep_eta = HiggsZlep_candidate_sorted[0].eta
                _lead_lep_phi = HiggsZlep_candidate_sorted[0].phi
                _lead_lep_iso = HiggsZlep_candidate_sorted[0].pfRelIso03_all
                _lead_lep_iso_chg = HiggsZlep_candidate_sorted[0].pfRelIso03_chg
            except:
                _lead_lep_pt  = -99.0
                _lead_lep_eta = -99.0
                _lead_lep_phi = -99.0
                _lead_lep_iso = -99.0
                _lead_lep_iso_chg = -99.0
            try:
                _trail_lep_pt  = HiggsZlep_candidate_sorted[1].pt
                _trail_lep_eta = HiggsZlep_candidate_sorted[1].eta
                _trail_lep_phi = HiggsZlep_candidate_sorted[1].phi
                _trail_lep_iso = HiggsZlep_candidate_sorted[1].pfRelIso03_all
                _trail_lep_iso_chg = HiggsZlep_candidate_sorted[1].pfRelIso03_chg
            except:
                _trail_lep_pt  = -99.0
                _trail_lep_eta = -99.0
                _trail_lep_phi = -99.0
                _trail_lep_iso = -99.0
                _trail_lep_iso_chg = -99.0

            self.out.fillBranch("leading_lep_pt{}".format(self.syst_suffix), _lead_lep_pt)
            self.out.fillBranch("leading_lep_eta{}".format(self.syst_suffix), _lead_lep_eta)
            self.out.fillBranch("leading_lep_phi{}".format(self.syst_suffix), _lead_lep_phi)
            self.out.fillBranch("leading_lep_iso{}".format(self.syst_suffix), _lead_lep_iso)
            self.out.fillBranch("leading_lep_iso_chg{}".format(self.syst_suffix), _lead_lep_iso_chg)
            self.out.fillBranch("trailing_lep_pt{}".format(self.syst_suffix), _trail_lep_pt)
            self.out.fillBranch("trailing_lep_eta{}".format(self.syst_suffix), _trail_lep_eta)
            self.out.fillBranch("trailing_lep_phi{}".format(self.syst_suffix), _trail_lep_phi)
            self.out.fillBranch("trailing_lep_iso{}".format(self.syst_suffix), _trail_lep_iso)
            self.out.fillBranch("trailing_lep_iso_chg{}".format(self.syst_suffix), _trail_lep_iso_chg)

        if num_jets >= 4:
            #Construct a Higgs boson candidate from b-tagged jets. We take the pair with mass closes to the Higgs mass
            if len(good_bjets) >= 2:
                Higgsbb_cand_p4 = good_bjets[0].p4() + good_bjets[1].p4() 
                Higgsbb_candidate = good_bjets[:2]
                #  for Hpair in itertools.combinations(good_bjets, 2):
                #     Higgs_cand_0 = Hpair[0].p4() + Hpair[1].p4()
                
                #     if abs(Higgs_cand_0.M()-self.Hmass) < abs(Higgsbb_cand_p4.M()-self.Hmass) or Higgsbb_cand_p4.M()==0.0:
                #             Higgsbb_cand_p4 = Higgs_cand_0
                #             Higgsbb_candidate = Hpair
                #We also look at the case where there are less than 2 b-tagged jets. Form a temp collection of the b-tagged jet and the other jets. 
            elif len(good_bjets) == 1:
                for jet in good_jets:
                    Higgs_cand_0 = good_bjets[0].p4() + jet.p4()
                    if jet.p4() == good_bjets[0].p4():
                        continue
                    if abs(Higgs_cand_0.M()-self.Hmass) < abs(Higgsbb_cand_p4.M()-self.Hmass) or Higgsbb_cand_p4.M()==0.0:
                        Higgsbb_cand_p4 = Higgs_cand_0
                        Higgsbb_candidate = [good_bjets[0],jet]
                        #The final case where we sadly have no b-tagged jets:(
            else:
                for Hpair in itertools.combinations(good_jets, 2):
                    Higgs_cand_0 = Hpair[0].p4() + Hpair[1].p4()
                
                    if abs(Higgs_cand_0.M()-self.Hmass) < abs(Higgsbb_cand_p4.M()-self.Hmass) or Higgsbb_cand_p4.M()==0.0:
                        Higgsbb_cand_p4 = Higgs_cand_0
                        Higgsbb_candidate = Hpair
                        #now we remove the jets that we already used from the collection so we dont use them twice
            if len(good_jets) >= 2:
                for jet in good_jets:
                    if jet.eta in [Higgsbb_candidate[0].eta,Higgsbb_candidate[1].eta] and jet.phi in [Higgsbb_candidate[0].phi,Higgsbb_candidate[1].phi]:
                        good_jets.remove(jet)
                    
            #Construct a Higgs boson candidate from jets and the lepton Z boson. We take the pair with mass closest to the Higgs boson mass
            if len(good_jets) >= 2 and len(good_leptons) >= 2: 
                for Zjetpair in itertools.combinations(good_jets, 2):#The possible Jet combinations for the Z boson
                    Zjet_cand_0 = Zjetpair[0].p4() + Zjetpair[1].p4()
                
                    for Zleppair in itertools.combinations(good_leptons, 2):#The possible lepton combinations for the Z boson
                        Zlep_cand_0 = Zleppair[0].p4() + Zleppair[1].p4()
                    
                        Higgs_cand_1 = Zjet_cand_0 + Zlep_cand_0
                        if abs(Higgs_cand_1.M()-self.Hmass) < abs(HiggsZZ_cand_p4.M()-self.Hmass) or HiggsZZ_cand_p4.M()==0.0:
                            HiggsZZ_cand_p4 = Higgs_cand_1
                            HiggsZjet_cand_p4 = Zjet_cand_0
                            HiggsZlep_cand_p4 = Zlep_cand_0
                            HiggsZjet_candidate = Zjetpair
                            HiggsZlep_candidate = Zleppair

            self.out.fillBranch("Higgsbb_cand_pt{}".format(self.syst_suffix), Higgsbb_cand_p4.Pt())
            self.out.fillBranch("Higgsbb_cand_eta{}".format(self.syst_suffix), Higgsbb_cand_p4.Eta())
            self.out.fillBranch("Higgsbb_cand_phi{}".format(self.syst_suffix), Higgsbb_cand_p4.Phi())
            self.out.fillBranch("Higgsbb_cand_mass{}".format(self.syst_suffix), Higgsbb_cand_p4.M())
            
            self.out.fillBranch("HiggsZZ_cand_pt{}".format(self.syst_suffix), HiggsZZ_cand_p4.Pt())
            self.out.fillBranch("HiggsZZ_cand_eta{}".format(self.syst_suffix), HiggsZZ_cand_p4.Eta())
            self.out.fillBranch("HiggsZZ_cand_phi{}".format(self.syst_suffix), HiggsZZ_cand_p4.Phi())
            self.out.fillBranch("HiggsZZ_cand_mass{}".format(self.syst_suffix), HiggsZZ_cand_p4.M())
            
            self.out.fillBranch("Zlep_cand_pt{}".format(self.syst_suffix), HiggsZlep_cand_p4.Pt())
            self.out.fillBranch("Zlep_cand_eta{}".format(self.syst_suffix), HiggsZlep_cand_p4.Eta())
            self.out.fillBranch("Zlep_cand_phi{}".format(self.syst_suffix), HiggsZlep_cand_p4.Phi())
            self.out.fillBranch("Zlep_cand_mass{}".format(self.syst_suffix), HiggsZlep_cand_p4.M())
            
            self.out.fillBranch("Zjet_cand_pt{}".format(self.syst_suffix), HiggsZjet_cand_p4.Pt())
            self.out.fillBranch("Zjet_cand_eta{}".format(self.syst_suffix), HiggsZjet_cand_p4.Eta())
            self.out.fillBranch("Zjet_cand_phi{}".format(self.syst_suffix), HiggsZjet_cand_p4.Phi())
            self.out.fillBranch("Zjet_cand_mass{}".format(self.syst_suffix), HiggsZjet_cand_p4.M())

            HH_cand_p4 = Higgsbb_cand_p4 + Higgsbb_cand_p4

            self.out.fillBranch("HH_cand_pt{}".format(self.syst_suffix), HH_cand_p4.Pt())
            self.out.fillBranch("HH_cand_eta{}".format(self.syst_suffix), HH_cand_p4.Eta())
            self.out.fillBranch("HH_cand_phi{}".format(self.syst_suffix), HH_cand_p4.Phi())
            self.out.fillBranch("HH_cand_mass{}".format(self.syst_suffix), HH_cand_p4.M())

            Higgsbb_candidate_sorted = sorted(Higgsbb_candidate, key=lambda x: x.pt, reverse=True)
            HiggsZlep_candidate_sorted = sorted(HiggsZlep_candidate, key=lambda x: x.pt, reverse=True)
            HiggsZjet_candidate_sorted = sorted(HiggsZjet_candidate, key=lambda x: x.pt, reverse=True)

#            Higgsbb_candidate.sort(key=lambda x: x.pt, reverse=True)
#            HiggsZlep_candidate.sort(key=lambda x: x.pt, reverse=True)
#            HiggsZjet_candidate.sort(key=lambda x: x.pt, reverse=True)

            #Lets look at the jets associated with the Higgs
            try:
                _lead_Hbb_pt  = Higgsbb_candidate_sorted[0].pt 
                _lead_Hbb_eta = Higgsbb_candidate_sorted[0].eta
                _lead_Hbb_phi = Higgsbb_candidate_sorted[0].phi
                _lead_Hbb_btag   = Higgsbb_candidate_sorted[0].btagDeepFlavB
            except:
                _lead_Hbb_pt  = -99.0
                _lead_Hbb_eta = -99.0
                _lead_Hbb_phi = -99.0
                _lead_Hbb_btag = -99.0
            try:
                _trail_Hbb_pt  = Higgsbb_candidate_sorted[1].pt 
                _trail_Hbb_eta = Higgsbb_candidate_sorted[1].eta 
                _trail_Hbb_phi = Higgsbb_candidate_sorted[1].phi
                _trail_Hbb_btag   = Higgsbb_candidate_sorted[1].btagDeepFlavB
            except:
                _trail_Hbb_pt  = -99.0
                _trail_Hbb_eta = -99.0
                _trail_Hbb_phi = -99.0
                _trail_Hbb_btag = -99.0

            self.out.fillBranch("leading_Hbb_pt{}".format(self.syst_suffix), _lead_Hbb_pt)
            self.out.fillBranch("leading_Hbb_eta{}".format(self.syst_suffix), _lead_Hbb_eta)
            self.out.fillBranch("leading_Hbb_phi{}".format(self.syst_suffix), _lead_Hbb_phi)
            self.out.fillBranch("leading_Hbb_btag{}".format(self.syst_suffix), _lead_Hbb_btag)            
            self.out.fillBranch("trailing_Hbb_pt{}".format(self.syst_suffix), _trail_Hbb_pt)
            self.out.fillBranch("trailing_Hbb_eta{}".format(self.syst_suffix), _trail_Hbb_eta)
            self.out.fillBranch("trailing_Hbb_phi{}".format(self.syst_suffix), _trail_Hbb_phi)
            self.out.fillBranch("trailing_Hbb_btag{}".format(self.syst_suffix), _trail_Hbb_btag)


            #Lets look at the leptons associated with the Z
            try:
                _lead_lep_pt  = HiggsZlep_candidate_sorted[0].pt
                _lead_lep_eta = HiggsZlep_candidate_sorted[0].eta 
                _lead_lep_phi = HiggsZlep_candidate_sorted[0].phi 
            except:
                _lead_lep_pt  = -99.0
                _lead_lep_eta = -99.0
                _lead_lep_phi = -99.0
            try:
                _trail_lep_pt  = HiggsZlep_candidate_sorted[1].pt
                _trail_lep_eta = HiggsZlep_candidate_sorted[1].eta 
                _trail_lep_phi = HiggsZlep_candidate_sorted[1].phi 
            except:
                _trail_lep_pt  = -99.0
                _trail_lep_eta = -99.0
                _trail_lep_phi = -99.0
            
            self.out.fillBranch("leading_lep_pt{}".format(self.syst_suffix), _lead_lep_pt)
            self.out.fillBranch("leading_lep_eta{}".format(self.syst_suffix), _lead_lep_eta)
            self.out.fillBranch("leading_lep_phi{}".format(self.syst_suffix), _lead_lep_phi)
            self.out.fillBranch("trailing_lep_pt{}".format(self.syst_suffix), _trail_lep_pt)
            self.out.fillBranch("trailing_lep_eta{}".format(self.syst_suffix), _trail_lep_eta)
            self.out.fillBranch("trailing_lep_phi{}".format(self.syst_suffix), _trail_lep_phi)
        
            #And the jets associated with the Z
            try:
                _lead_jet_pt  = HiggsZjet_candidate_sorted[0].pt
                _lead_jet_eta = HiggsZjet_candidate_sorted[0].eta 
                _lead_jet_phi = HiggsZjet_candidate_sorted[0].phi 
                _lead_jet_qgl = HiggsZjet_candidate_sorted[0].qgl
            except:
                _lead_jet_pt  = -99.0
                _lead_jet_eta = -99.0
                _lead_jet_phi = -99.0
                _lead_jet_qgl = -99.0
            
            try:
                _trail_jet_pt  = HiggsZjet_candidate_sorted[1].pt
                _trail_jet_eta = HiggsZjet_candidate_sorted[1].eta 
                _trail_jet_phi = HiggsZjet_candidate_sorted[1].phi 
                _trail_jet_qgl = HiggsZjet_candidate_sorted[1].qgl
            except:
                _trail_jet_pt  = -99.0
                _trail_jet_eta = -99.0
                _trail_jet_phi = -99.0
                _trail_jet_qgl = -99.0

            self.out.fillBranch("leading_jet_pt{}".format(self.syst_suffix), _lead_jet_pt)
            self.out.fillBranch("leading_jet_eta{}".format(self.syst_suffix), _lead_jet_eta)
            self.out.fillBranch("leading_jet_phi{}".format(self.syst_suffix), _lead_jet_phi)
            self.out.fillBranch("leading_jet_qgl{}".format(self.syst_suffix), _lead_jet_qgl)
            self.out.fillBranch("trailing_jet_pt{}".format(self.syst_suffix), _trail_jet_pt)
            self.out.fillBranch("trailing_jet_eta{}".format(self.syst_suffix), _trail_jet_eta)
            self.out.fillBranch("trailing_jet_phi{}".format(self.syst_suffix), _trail_jet_phi)
            self.out.fillBranch("trailing_jet_qgl{}".format(self.syst_suffix), _trail_jet_qgl)

        _dphi_j_met = tk.deltaPhi(good_jets[0], met.phi) if len(good_jets) else -99.0
        _lead_jet_pt = good_jets[0].pt if len(good_jets) else -99.0
        _lead_bjet_pt = good_bjets[0].pt if len(good_bjets) else -99.0


        self.out.fillBranch("jetHT{}".format(self.syst_suffix), jetHT)
        self.out.fillBranch("ngood_jets{}".format(self.syst_suffix), num_jets)
        self.out.fillBranch("ngood_jets_noHbb{}".format(self.syst_suffix), len(good_jets))
        self.out.fillBranch("ngood_bjets{}".format(self.syst_suffix), len(good_bjets))
        self.out.fillBranch("ngood_bjetsM{}".format(self.syst_suffix), len(good_bjetsM))
        self.out.fillBranch("ngood_bjetsT{}".format(self.syst_suffix), len(good_bjetsT))
        self.out.fillBranch("ngood_bjetsC{}".format(self.syst_suffix), len(good_bjetsC))
        self.out.fillBranch("lead_jet_pt{}".format(self.syst_suffix), _lead_jet_pt)
        self.out.fillBranch("lead_bjet_pt{}".format(self.syst_suffix), _lead_bjet_pt)
        self.out.fillBranch("delta_phi_j_met{}".format(self.syst_suffix), _dphi_j_met)


        # process taus
        had_taus = []
        for tau in taus:
            if tk.closest(tau, good_leptons)[1] < 0.4:
                continue
            # only hadronic tau decay
            if tau.decayMode != 5:
                continue
            if tau.pt > 18 and abs(tau.eta) <= 2.3:
                had_taus.append(tau)

        self.out.fillBranch("nhad_taus{}".format(self.syst_suffix), len(had_taus))
        self.out.fillBranch("lead_tau_pt{}".format(self.syst_suffix), had_taus[0].pt if len(had_taus) else 0)

        # Let remove the negative categories with no obvious meaning meaning
        # This will reduce the size of most of the bacground and data
        if (lep_category > 0 and event_category > 0 and num_jets >= 4  and len(good_leptons) == 2 and nextra_leptons == 0 and event_category < 10 and len(had_taus)==0): 
#            return True
            self.out.fillBranch("good_event{}".format(self.syst_suffix),1)
            return True
        else:
#            return False
            self.out.fillBranch("good_event{}".format(self.syst_suffix),0)
            self.out.fillBranch("met_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("met_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("ngood_leptons{}".format(self.syst_suffix), -99)
            self.out.fillBranch("nextra_leptons{}".format(self.syst_suffix), -99)
            self.out.fillBranch("lep_category{}".format(self.syst_suffix), -99) 
            self.out.fillBranch("event_category{}".format(self.syst_suffix), -99) 
            self.out.fillBranch("dR_l1l2{}".format(self.syst_suffix), -99)
            self.out.fillBranch("event_category_30{}".format(self.syst_suffix), -99)
            self.out.fillBranch("event_category_35{}".format(self.syst_suffix), -99)
            self.out.fillBranch("event_category_40{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("event_category_mini30{}".format(self.syst_suffix), -99)
            self.out.fillBranch("event_category_mini35{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("leading_Hbb_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_Hbb_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_Hbb_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_Hbb_btag{}".format(self.syst_suffix), -99)

            self.out.fillBranch("trailing_Hbb_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_Hbb_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_Hbb_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_Hbb_btag{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("leading_lep_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_lep_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_lep_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_lep_iso{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_lep_iso_chg{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_lep_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_lep_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_lep_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_lep_iso{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_lep_iso_chg{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_lep_flavor{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_lep_flavor{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("leading_jet_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_jet_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_jet_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leading_jet_qgl{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("trailing_jet_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_jet_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_jet_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("trailing_jet_qgl{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("met_filter{}".format(self.syst_suffix), -99)
            
            #Adding in the Higgs boson candidate variables
            self.out.fillBranch("Higgsbb_cand_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Higgsbb_cand_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Higgsbb_cand_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Higgsbb_cand_mass{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("HiggsZZ_cand_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("HiggsZZ_cand_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("HiggsZZ_cand_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("HiggsZZ_cand_mass{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("Zlep_cand_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Zlep_cand_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Zlep_cand_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Zlep_cand_mass{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("Zjet_cand_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Zjet_cand_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Zjet_cand_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("Zjet_cand_mass{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("HH_cand_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("HH_cand_eta{}".format(self.syst_suffix), -99)
            self.out.fillBranch("HH_cand_phi{}".format(self.syst_suffix), -99)
            self.out.fillBranch("HH_cand_mass{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("jetHT{}".format(self.syst_suffix), -99)
            self.out.fillBranch("ngood_jets{}".format(self.syst_suffix), -99)
            self.out.fillBranch("ngood_jets_noHbb{}".format(self.syst_suffix), -99)
            self.out.fillBranch("ngood_bjets{}".format(self.syst_suffix), -99)
            self.out.fillBranch("ngood_bjetsM{}".format(self.syst_suffix), -99)
            self.out.fillBranch("ngood_bjetsT{}".format(self.syst_suffix), -99)
            self.out.fillBranch("ngood_bjetsC{}".format(self.syst_suffix), -99)
            self.out.fillBranch("lead_jet_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("lead_bjet_pt{}".format(self.syst_suffix), -99)
            self.out.fillBranch("delta_phi_j_met{}".format(self.syst_suffix), -99)
            self.out.fillBranch("leptoncharge{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("nhad_taus{}".format(self.syst_suffix), -99)
            self.out.fillBranch("lead_tau_pt{}".format(self.syst_suffix), -99)
            
            self.out.fillBranch("w_muon_SF{}".format(self.syst_suffix), -99)
            self.out.fillBranch("w_muon_SFUp{}".format(self.syst_suffix), -99)
            self.out.fillBranch("w_muon_SFDown{}".format(self.syst_suffix), -99)
            self.out.fillBranch("w_electron_SF{}".format(self.syst_suffix), -99)
            self.out.fillBranch("w_electron_SFUp{}".format(self.syst_suffix), -99)
            self.out.fillBranch("w_electron_SFDown{}".format(self.syst_suffix), -99)
            self.out.fillBranch("w_btag_SF{}".format(self.syst_suffix),-99)
            return True
