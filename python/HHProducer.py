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


class HHProducer(Module):
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

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("met_pt{}".format(self.syst_suffix), "F")
        self.out.branch("met_phi{}".format(self.syst_suffix), "F")
        self.out.branch("ngood_leptons{}".format(self.syst_suffix), "I")
        self.out.branch("nextra_leptons{}".format(self.syst_suffix), "I")
        self.out.branch("lep_category{}".format(self.syst_suffix), "I")

        self.out.branch("leading_Hbb_pt{}".format(self.syst_suffix), "F")
        self.out.branch("leading_Hbb_eta{}".format(self.syst_suffix), "F")
        self.out.branch("leading_Hbb_phi{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_Hbb_pt{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_Hbb_eta{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_Hbb_phi{}".format(self.syst_suffix), "F")

        self.out.branch("leading_lep_pt{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_eta{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_phi{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_pt{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_eta{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_phi{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_flavor{}".format(self.syst_suffix), "I")

        self.out.branch("leading_jet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("leading_jet_eta{}".format(self.syst_suffix), "F")
        self.out.branch("leading_jet_phi{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_jet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_jet_eta{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_jet_phi{}".format(self.syst_suffix), "F")

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

        #DeltaR variables
        self.out.branch("dR_l1l2{}".format(self.syst_suffix), "F")
        self.out.branch("dR_l1j1{}".format(self.syst_suffix), "F")
        self.out.branch("dR_l1j2{}".format(self.syst_suffix), "F")
        self.out.branch("dR_l1b1{}".format(self.syst_suffix), "F")
        self.out.branch("dR_l1b2{}".format(self.syst_suffix), "F")
        self.out.branch("dR_l2j1{}".format(self.syst_suffix), "F")
        self.out.branch("dR_l2j2{}".format(self.syst_suffix), "F")
        self.out.branch("dR_l2b1{}".format(self.syst_suffix), "F")
        self.out.branch("dR_l2b2{}".format(self.syst_suffix), "F")
        self.out.branch("dR_j1j2{}".format(self.syst_suffix), "F")
        self.out.branch("dR_j1b1{}".format(self.syst_suffix), "F")
        self.out.branch("dR_j1b2{}".format(self.syst_suffix), "F")
        self.out.branch("dR_j2b1{}".format(self.syst_suffix), "F")
        self.out.branch("dR_j2b2{}".format(self.syst_suffix), "F")
        self.out.branch("dR_b1b2{}".format(self.syst_suffix), "F")

        #some other variables
        self.out.branch("delta_phi_ll{}".format(self.syst_suffix), "F")
        self.out.branch("delta_eta_ll{}".format(self.syst_suffix), "F")
        self.out.branch("delta_R_ll{}".format(self.syst_suffix), "F")
        self.out.branch("delta_phi_jj{}".format(self.syst_suffix), "F")
        self.out.branch("delta_eta_jj{}".format(self.syst_suffix), "F")
        self.out.branch("delta_R_jj{}".format(self.syst_suffix), "F")

        self.out.branch("ngood_jets{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_jets_noHbb{}".format(self.syst_suffix), "I")
        self.out.branch("ngood_bjets{}".format(self.syst_suffix), "I")
        self.out.branch("lead_jet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("lead_bjet_pt{}".format(self.syst_suffix), "F")
        self.out.branch("delta_phi_j_met{}".format(self.syst_suffix), "F")

        self.out.branch("nhad_taus{}".format(self.syst_suffix), "I")
        self.out.branch("lead_tau_pt{}".format(self.syst_suffix), "F")

        if self.isMC and len(self.syst_suffix)==0:
            self.out.branch("w_muon_SF{}".format(self.syst_suffix), "F")
            self.out.branch("w_muon_SFUp{}".format(self.syst_suffix), "F")
            self.out.branch("w_muon_SFDown{}".format(self.syst_suffix), "F")
            self.out.branch("w_electron_SF{}".format(self.syst_suffix), "F")
            self.out.branch("w_electron_SFUp{}".format(self.syst_suffix), "F")
            self.out.branch("w_electron_SFDown{}".format(self.syst_suffix), "F")


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def electron_id(self, electron, wp):
        pass_id = 0
        if (self.era == "2016" and wp == "80"):
            return electron.mvaSpring16GP_WP80
        elif (self.era == "2016" and wp == "90"):
            return electron.mvaSpring16GP_WP90

        elif (self.era == "2017" and wp == "80"):
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
        elif (self.era == "2017" and wp == "90"):
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
        elif (self.era == "2017" and wp == "WPL"):
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

        elif (self.era == "2018" and wp == "80"):
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
        elif (self.era == "2018" and wp == "90"):
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
        elif (self.era == "2018" and wp == "WPL"):
            try:
                pass_id = electron.mvaFall17V2Iso_WPL
            except:
                try:
                    pass_id = electron.mvaFall17V1Iso_WPL
                except:
                    try:
                        pass_id = electron.mvaFall17Iso_WPL
                    except ValueError:
                        print "[error] not mvaFall18 electron id found ... "

            return pass_id

    def btag_id(self, wp):
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
           and (flag.BadChargedCandidateFilter)
           and (flag.BadPFMuonFilter)
        )


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
        try:
            var_jet_pts = getattr(event,  "Jet_pt_{}".format(syst_var), None)
            if var_jet_pts:
                for i,jet in enumerate(jets):
                    jet.pt = var_jet_pts[i]
            else:
                print 'WARNING: jet pts with variation {}'
                'not available, using the nominal value'.format(syst_var)
        except:
            var_jet_pts = getattr(event,  "Jet_pt_nom", None)
            for i,jet in enumerate(jets):
                jet.pt = var_jet_pts[i]

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
        if self.isMC:
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

	muons.sort(key=lambda muon: muon.pt, reverse=True)
        electrons.sort(key=lambda el: el.pt, reverse=True)
        # Choose tight-quality e/mu for event categorization
        for idx,mu in enumerate(muons):
            isoLep   = mu.miniPFRelIso_all
            pass_ips = abs(mu.dxy) < 0.02 and abs(mu.dz) < 0.1
            pass_fid = abs(mu.eta) < 2.4 and mu.pt >= (20 if idx==0 else 10)
            pass_ids = mu.mediumId and isoLep <= 0.1
            if pass_fid and pass_ids and pass_ips:
                good_muons.append(mu)
        for idy,el in enumerate(electrons):
            id_CB = el.cutBased
            # changing to MVA based ID :
            if el.pt >= (25 if idy==0 else 20) and abs(el.eta) <= 2.5 and self.electron_id(el, "90"):
                good_electrons.append(el)

        # let sort the muons in pt
        good_muons.sort(key=lambda x: x.pt, reverse=True)
        good_electrons.sort(key=lambda x: x.pt, reverse=True)

        # Find any remaining e/mu that pass looser selection
        extra_leptons = []
        for mu in muons:
            isoLep   = mu.miniPFRelIso_all
            pass_ids = mu.softId and isoLep <= 0.4
            pass_fid = abs(mu.eta) < 2.4 and mu.pt >= 10
            if tk.closest(mu, good_muons)[1] < 0.01:
                continue
            if pass_fid and pass_ids:
                extra_leptons.append(mu)

        for el in electrons:
            pass_fid = abs(el.eta) < 2.5 and el.pt >= 10
            if tk.closest(el, good_electrons)[1] < 0.01:
                continue
            if pass_fid and self.electron_id(el, "WPL"):
                extra_leptons.append(el)

        good_leptons = good_electrons + good_muons
        good_leptons.sort(key=lambda x: x.pt, reverse=True)

        _lead_lep_pt = good_leptons[0].pt if len(good_leptons) else 0.0
        _lead_lep_eta = good_leptons[0].eta if len(good_leptons) else 0.0
        _trail_lep_pt = good_leptons[1].pt if len(good_leptons) >= 2 else 0.0
        _trail_lep_eta = good_leptons[1].eta if len(good_leptons) >= 2 else 0.0
        _lead_lep_pdgId = good_leptons[0].pdgId if len(good_leptons)  else 0.0
        _trail_lep_pdgId = good_leptons[1].pdgId if len(good_leptons) >= 2 else 0.0

        self.out.fillBranch("leading_lep_pt{}".format(self.syst_suffix), _lead_lep_pt)
        self.out.fillBranch("leading_lep_eta{}".format(self.syst_suffix), _lead_lep_eta)
        self.out.fillBranch("trailing_lep_pt{}".format(self.syst_suffix), _trail_lep_pt)
        self.out.fillBranch("trailing_lep_eta{}".format(self.syst_suffix), _trail_lep_eta)
	_leading_lep_flavor = 0
	if len(good_muons) and len(good_electrons):
		if good_muons[0].pt > good_electrons[0].pt: _leading_lep_flavor = 1
        self.out.fillBranch("leading_lep_flavor{}".format(self.syst_suffix),_leading_lep_flavor)

        ngood_leptons = len(good_leptons)
        nextra_leptons = len(extra_leptons)

        if False:
            print "number of leptons [all, good, extra]: ", ngood_leptons, " : ", nextra_leptons
            print "        CBId electrons : ", [e.cutBased for e in good_electrons]
            print "        WP90 electrons : ", [e.mvaFall17Iso_WP90 for e in good_electrons]
            print "             muons     : ", [e.tightId for e in good_muons]
            print "        lepton pts     : ", [e.pt for e in good_leptons]
        if ngood_leptons == 2:
            # constructing the signal region
            if (good_leptons[0].pdgId * good_leptons[1].pdgId) == -11*11:
                lep_category = 1 # EE category
            if (good_leptons[0].pdgId * good_leptons[1].pdgId) == -13*13:
                lep_category = 2 # MM category
            if (good_leptons[0].pdgId * good_leptons[1].pdgId) == -11*13:
                lep_category = 3 # EM category
            if (good_leptons[0].pdgId * good_leptons[1].pdgId) == 121:
                lep_category = 4 # EE SS category 
            if (good_leptons[0].pdgId * good_leptons[1].pdgId) == 169:
                lep_category = 5 # MM SS category
            if (good_leptons[0].pdgId * good_leptons[1].pdgId) == 11*13:
                lep_category = 6 # EM SS category
        else:
            lep_category = 7 #any other amount of leptons
            
        self.out.fillBranch("ngood_leptons{}".format(self.syst_suffix), ngood_leptons)
        self.out.fillBranch("nextra_leptons{}".format(self.syst_suffix), nextra_leptons)
        self.out.fillBranch("lep_category{}".format(self.syst_suffix), lep_category)

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
                if abs(good_leptons[0].pdgId) == 11:
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

        # process jet
        good_jets  = []
        good_bjets = []
        for jet in jets:
            if jet.btagDeepB > self.btag_id("loose"):
                jet.pt = jet.pt*jet.bRegCorr
            if jet.pt < 30.0 or  abs(jet.eta) > 2.4: # 4.7:
                continue
            if not jet.jetId:
                continue
            if tk.closest(jet, good_leptons)[1] < 0.4:
                continue
            good_jets.append(jet)
            # Count b-tag with medium WP DeepCSV
            # ref : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
            if jet.btagDeepB > self.btag_id("loose"):
                good_bjets.append(jet)        

        good_jets.sort(key=lambda jet: jet.pt, reverse=True)
        good_bjets.sort(key=lambda jet: jet.btagDeepB, reverse=True)
        #We will remove jets later so better count them now
	num_jets = len(good_jets)

        #Set up for the Higgs candicates!
        Higgsbb_cand_p4 = ROOT.TLorentzVector()
        HiggsZZ_cand_p4 = ROOT.TLorentzVector()
	HiggsZjet_cand_p4 = ROOT.TLorentzVector()
	HiggsZlep_cand_p4 = ROOT.TLorentzVector()
        Higgs_cand_0 = ROOT.TLorentzVector()
        Zjet_cand_0 = ROOT.TLorentzVector()
        Zlep_cand_0 = ROOT.TLorentzVector()
        Higgs_cand_1 = ROOT.TLorentzVector()
        Higgsbb_candidate = []
        HiggsZjet_candidate = []
        HiggsZlep_candidate = []

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
            #temp_jets = good_jets
            #temp_jets = []
            #temp_jets.append(good_bjets[0])
            #temp_jets.extend(good_jets)
            for jet in good_jets:
                Higgs_cand_0 = good_bjets[0].p4() + jet.p4()
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
                    # print('yoooooo')
                    good_jets.remove(jet)
        
        #Construct a Higgs boson candidate from jets and the lepton Z boson. We take the pair with mass closest to the Higgs boson mass
        if len(good_jets) >= 2 and len(good_leptons) >= 2: 
            for Zjetpair in itertools.combinations(good_jets, 2):#The possible Jet combinations for the Z boson
                Zjet_cand_0 = Zjetpair[0].p4() + Zjetpair[1].p4()

                for Zleppair in itertools.combinations(good_leptons, 2):#The possible lepton combinations for the Z boson
                    if lep_category > 2: continue#ensuring OSSF via the lep_category above
                    #if (Zleppair[0].pdgId * Zleppair[1].pdgId) not in [-169, -121]: continue#ensuring OSSF
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
 
        #Lets look at the jets associated with the Higgs
        try:
            _lead_Hbb_pt  = Higgsbb_candidate[0].pt 
            _lead_Hbb_eta = Higgsbb_candidate[0].eta
            _lead_Hbb_phi = Higgsbb_candidate[0].phi
        except:
            _lead_Hbb_pt  = -99.0
            _lead_Hbb_eta = -99.0
            _lead_Hbb_phi = -99.0
        try:
            _trail_Hbb_pt  = Higgsbb_candidate[1].pt 
            _trail_Hbb_eta = Higgsbb_candidate[1].eta 
            _trail_Hbb_phi = Higgsbb_candidate[1].phi 
        except:
            _trail_Hbb_pt  = -99.0
            _trail_Hbb_eta = -99.0
            _trail_Hbb_phi = -99.0

        self.out.fillBranch("leading_Hbb_pt{}".format(self.syst_suffix), _lead_Hbb_pt)
        self.out.fillBranch("leading_Hbb_eta{}".format(self.syst_suffix), _lead_Hbb_eta)
        self.out.fillBranch("leading_Hbb_phi{}".format(self.syst_suffix), _lead_Hbb_phi)
        self.out.fillBranch("trailing_Hbb_pt{}".format(self.syst_suffix), _trail_Hbb_pt)
        self.out.fillBranch("trailing_Hbb_eta{}".format(self.syst_suffix), _trail_Hbb_eta)
        self.out.fillBranch("trailing_Hbb_phi{}".format(self.syst_suffix), _trail_Hbb_phi)

        #Lets look at the leptons associated with the Z
        try:
            _lead_lep_pt  = HiggsZlep_candidate[0].pt
            _lead_lep_eta = HiggsZlep_candidate[0].eta 
            _lead_lep_phi = HiggsZlep_candidate[0].phi 
        except:
            _lead_lep_pt  = -99.0
            _lead_lep_eta = -99.0
            _lead_lep_phi = -99.0
        try:
            _trail_lep_pt  = HiggsZlep_candidate[1].pt
            _trail_lep_eta = HiggsZlep_candidate[1].eta 
            _trail_lep_phi = HiggsZlep_candidate[1].phi 
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
            _lead_jet_pt  = HiggsZjet_candidate[0].pt
            _lead_jet_eta = HiggsZjet_candidate[0].eta 
            _lead_jet_phi = HiggsZjet_candidate[0].phi 
        except:
            _lead_jet_pt  = -99.0
            _lead_jet_eta = -99.0
            _lead_jet_phi = -99.0
        try:
            _trail_jet_pt  = HiggsZjet_candidate[1].pt
            _trail_jet_eta = HiggsZjet_candidate[1].eta 
            _trail_jet_phi = HiggsZjet_candidate[1].phi 
        except:
            _trail_jet_pt  = -99.0
            _trail_jet_eta = -99.0
            _trail_jet_phi = -99.0

        self.out.fillBranch("leading_jet_pt{}".format(self.syst_suffix), _lead_jet_pt)
        self.out.fillBranch("leading_jet_eta{}".format(self.syst_suffix), _lead_jet_eta)
        self.out.fillBranch("leading_jet_phi{}".format(self.syst_suffix), _lead_jet_phi)
        self.out.fillBranch("trailing_jet_pt{}".format(self.syst_suffix), _trail_jet_pt)
        self.out.fillBranch("trailing_jet_eta{}".format(self.syst_suffix), _trail_jet_eta)
        self.out.fillBranch("trailing_jet_phi{}".format(self.syst_suffix), _trail_jet_phi)

        #Delta R (for first lepton)
        try:
            dR_l1l2 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi,)
        except: 
            dR_l1l2 = -99

        try:
            dR_l1j1 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi,)
        except:
            dR_l1j1 = -99

        try:
            dR_l1j2 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi,)
        except:
            dR_l1j2 = -99

        try:
            dR_l1b1 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi,)
        except:
            dR_l1b1 = -99

        try:
            dR_l1b2 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)
        except:
            dR_l1b2 = -99

       #Delta R (for second lepton)
        try:
            dR_l2j1 = tk.deltaR(HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi, HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi,)
        except:
            dR_l2j1 = -99

        try:
            dR_l2j2 = tk.deltaR(HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi, HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi,)
        except:
            dR_l2j2 = -99

        try:
            dR_l2b1 = tk.deltaR(HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi, Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi,)
        except:
            dR_l2b1 = -99

        try:
            dR_l2b2 = tk.deltaR(HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)
        except:
            dR_l2b2 = -99

        #Delta R (for first jet)
        try:
            dR_j1j2 = tk.deltaR(HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi, HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi,)
        except:
            dR_j1j2 = -99

        try:
            dR_j1b1 = tk.deltaR(HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi, Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi,)
        except:
            dR_j1b1 = -99

        try:
            dR_j1b2 = tk.deltaR(HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)
        except:
            dR_j1b2 = -99

        #Delta R (for second jet)
        try:
            dR_j2b1 = tk.deltaR(HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi, Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi,)
        except:
            dR_j2b1 = -99

        try:
            dR_j2b2 = tk.deltaR(HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)
        except:
            dR_j2b2 = -99

        #Delta R (b-tagged jets)
        try:  
            dR_b1b2 = tk.deltaR(Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)  
        except:  
            dR_b1b2 = -99


        self.out.fillBranch("dR_l1l2{}".format(self.syst_suffix), dR_l1l2)
        self.out.fillBranch("dR_l1j1{}".format(self.syst_suffix), dR_l1j1)
        self.out.fillBranch("dR_l1j2{}".format(self.syst_suffix), dR_l1j2)
        self.out.fillBranch("dR_l1b1{}".format(self.syst_suffix), dR_l1b1)
        self.out.fillBranch("dR_l1b2{}".format(self.syst_suffix), dR_l1b2)
        self.out.fillBranch("dR_l2j1{}".format(self.syst_suffix), dR_l2j1)
        self.out.fillBranch("dR_l2j2{}".format(self.syst_suffix), dR_l2j2)
        self.out.fillBranch("dR_l2b1{}".format(self.syst_suffix), dR_l2b1)
        self.out.fillBranch("dR_l2b2{}".format(self.syst_suffix), dR_l2b2)
        self.out.fillBranch("dR_j1j2{}".format(self.syst_suffix), dR_j1j2)
        self.out.fillBranch("dR_j1b1{}".format(self.syst_suffix), dR_j1b1)
        self.out.fillBranch("dR_j1b2{}".format(self.syst_suffix), dR_j1b2)
        self.out.fillBranch("dR_j2b1{}".format(self.syst_suffix), dR_j2b1)
        self.out.fillBranch("dR_j2b2{}".format(self.syst_suffix), dR_j2b2)
        self.out.fillBranch("dR_b1b2{}".format(self.syst_suffix), dR_b1b2)

        _dphi_j_met = tk.deltaPhi(good_jets[0], met.phi) if len(good_jets) else -99.0
        _lead_jet_pt = good_jets[0].pt if len(good_jets) else -99.0
        _lead_bjet_pt = good_bjets[0].pt if len(good_bjets) else -99.0

        self.out.fillBranch("ngood_jets{}".format(self.syst_suffix), num_jets)
        self.out.fillBranch("ngood_jets_noHbb{}".format(self.syst_suffix), len(good_jets))
        self.out.fillBranch("ngood_bjets{}".format(self.syst_suffix), len(good_bjets))
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
        # This will reduce the size of most of the background and data
	if (len(good_leptons) > 1 and len(good_jets) > 1):
            return True
        else:
            return False

