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


class ScaleFactorProducer(Module):
    def __init__(self, isMC, era, period, do_syst=False, syst_var='',):
        self.isMC = isMC
        self.period = period
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
        # 1 = dimuon channel, 2 = dielectron channel
        self.out.branch("event_category{}".format(self.syst_suffix), "I") 
        # used for ABCD method
        # 1 = A OS+ISO
        # 2 = B OS+non-ISO
        # 3 = C SS+ISO
        # 4 = D SS+non-ISO

        self.out.branch("met_filter{}".format(self.syst_suffix), "I")
        self.out.branch("leading_lep_pt{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_eta{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_phi{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_pt{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_eta{}".format(self.syst_suffix), "F")
        self.out.branch("trailing_lep_phi{}".format(self.syst_suffix), "F")
        self.out.branch("leading_lep_flavor{}".format(self.syst_suffix), "I")
        self.out.branch("trailing_lep_flavor{}".format(self.syst_suffix), "I")

        self.out.branch("DY_mass{}".format(self.syst_suffix),"F")
        self.out.branch("DY_pt{}".format(self.syst_suffix),"F")
        self.out.branch("electron_trigger{}".format(self.syst_suffix),"I")
        self.out.branch("muon_trigger{}".format(self.syst_suffix),"I")

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

    def electron_id_noiso(self, electron, wp):
        pass_id = 0
        if (self.era == "2016" and wp == "80"):
            return electron.mvaSpring16GP_WP80
        elif (self.era == "2016" and wp == "90"):
            return electron.mvaSpring16GP_WP90

        elif (self.era == "2017" and wp == "80"):
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
        elif (self.era == "2017" and wp == "90"):
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
        elif (self.era == "2017" and wp == "WPL"):
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

        elif (self.era == "2018" and wp == "80"):
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
        elif (self.era == "2018" and wp == "90"):
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
        elif (self.era == "2018" and wp == "WPL"):
            try:
                pass_id = electron.mvaFall17V2noIso_WPL
            except:
                try:
                    pass_id = electron.mvaFall17V1noIso_WPL
                except:
                    try:
                        pass_id = electron.mvaFall17noIso_WPL
                    except ValueError:
                        print "[error] not mvaFall18 electron id found ... "

            return pass_id


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


    def pass_electron_trigger(self, hlt):
        fired = 0
        if (self.era == "2016"):
            if self.isMC:
                fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele25_eta2p1_WPTight_Gsf or hlt.Ele27_eta2p1_WPLoose_Gsf or hlt.Ele27_WPTight_Gsf or hlt.Ele35_WPLoose_Gsf)
            else:
                if "Run2016H" == self.period:
                    fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele25_eta2p1_WPTight_Gsf or hlt.Ele27_eta2p1_WPLoose_Gsf or hlt.Ele27_WPTight_Gsf)
                else:
                    fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele25_eta2p1_WPTight_Gsf or hlt.Ele27_eta2p1_WPLoose_Gsf or hlt.Ele27_WPTight_Gsf or hlt.Ele35_WPLoose_Gsf)
        elif (self.era == "2017"):
            if self.isMC:
                fired = (hlt.Ele27_WPTight_Gsf or hlt.Ele32_WPTight_Gsf or hlt.Ele35_WPTight_Gsf or hlt.Ele40_WPTight_Gsf or hlt.Ele32_WPTight_Gsf_L1DoubleEG or hlt.Photon200 or hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or hlt.DiEle27_WPTightCaloOnly_L1DoubleEG or hlt.DoubleEle33_CaloIdL_MW or hlt.DoubleEle25_CaloIdL_MW or hlt.DoublePhoton70 or hlt.Ele115_CaloIdVT_GsfTrkIdT)
            else:
                if "Run2017B" == self.period:
                    fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele27_WPTight_Gsf or hlt.Ele35_WPTight_Gsf or hlt.Ele38_WPTight_Gsf or hlt.Ele40_WPTight_Gsf)
                elif "Run2017C" == self.period:
                    fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele27_WPTight_Gsf or hlt.Ele35_WPTight_Gsf or hlt.Ele38_WPTight_Gsf or hlt.Ele40_WPTight_Gsf)
                else:
                    fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or hlt.DiEle27_WPTightCaloOnly_L1DoubleEG or hlt.DoubleEle33_CaloIdL_MW or hlt.DoubleEle25_CaloIdL_MW or hlt.DoublePhoton70 or hlt.Ele115_CaloIdVT_GsfTrkIdT or hlt.Ele27_WPTight_Gsf or hlt.Ele32_WPTight_Gsf or hlt.Ele35_WPTight_Gsf or hlt.Ele32_WPTight_Gsf_L1DoubleEG or hlt.Photon200)
        elif (self.era == "2018"):
            if self.isMC:
                fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or hlt.DiEle27_WPTightCaloOnly_L1DoubleEG or hlt.DoubleEle33_CaloIdL_MW or hlt.DoubleEle25_CaloIdL_MW or hlt.DoubleEle27_CaloIdL_MW or hlt.DoublePhoton70 or hlt.Ele115_CaloIdVT_GsfTrkIdT or hlt.Ele27_WPTight_Gsf or hlt.Ele32_WPTight_Gsf or hlt.Ele35_WPTight_Gsf or hlt.Ele38_WPTight_Gsf or hlt.Ele40_WPTight_Gsf or hlt.Ele32_WPTight_Gsf_L1DoubleEG or hlt.Photon200)
            else:
                if "Run2018A" == self.period or "Run2018B" == self.period:
                    fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or hlt.DiEle27_WPTightCaloOnly_L1DoubleEG or hlt.DoubleEle33_CaloIdL_MW or hlt.DoubleEle25_CaloIdL_MW or hlt.DoubleEle27_CaloIdL_MW or hlt.DoublePhoton70 or hlt.Ele115_CaloIdVT_GsfTrkIdT or hlt.Ele27_WPTight_Gsf or hlt.Ele32_WPTight_Gsf or hlt.Ele35_WPTight_Gsf or hlt.Ele38_WPTight_Gsf or hlt.Ele40_WPTight_Gsf or hlt.Ele32_WPTight_Gsf_L1DoubleEG or hlt.Photon200 )
                else:
                    fired = (hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or hlt.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or hlt.DiEle27_WPTightCaloOnly_L1DoubleEG or hlt.DoubleEle33_CaloIdL_MW or hlt.DoubleEle25_CaloIdL_MW or hlt.DoubleEle27_CaloIdL_MW or hlt.DoublePhoton70 or hlt.Ele115_CaloIdVT_GsfTrkIdT or hlt.Ele27_WPTight_Gsf or hlt.Ele32_WPTight_Gsf or hlt.Ele35_WPTight_Gsf or hlt.Ele38_WPTight_Gsf or hlt.Ele40_WPTight_Gsf or hlt.Ele32_WPTight_Gsf_L1DoubleEG or hlt.Photon200)
        return fired

    def pass_muon_trigger(self, hlt):
        fired = 0
        if (self.era == "2016"):
            if self.isMC:
                fired = (hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or hlt.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ or hlt.IsoMu20 or hlt.IsoTkMu20 or hlt.IsoMu22 or hlt.IsoTkMu22 or hlt.IsoMu24 or hlt.IsoTkMu24)
            else:
                if "Run2016H" == self.period:
                    fired = (hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or hlt.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ or hlt.IsoMu20 or hlt.IsoTkMu20 or hlt.IsoMu22 or hlt.IsoTkMu22 or hlt.IsoMu24 or hlt.IsoTkMu24)
                else:
                    fired = (hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ or hlt.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ or hlt.IsoMu20 or hlt.IsoTkMu20 or hlt.IsoMu22 or hlt.IsoTkMu22 or hlt.IsoMu24 or hlt.IsoTkMu24)

        elif (self.era == "2017"):
            if self.isMC:
                fired = (hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 or hlt.Mu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL or hlt.IsoMu20 or hlt.IsoMu24 or hlt.IsoMu27 or hlt.IsoMu30 or hlt.Mu50)
            else:
                if "Run2017B" == self.period:
                    fired = (hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or hlt.Mu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL or hlt.IsoMu20 or hlt.IsoMu24 or hlt.IsoMu27)
                elif "Run2017C" == self.period:
                    fired = (hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or hlt.Mu8_TrkIsoVVL or hlt.Mu17_TrkIsoVVL or hlt.IsoMu20 or hlt.IsoMu24 or hlt.IsoMu27)
                else:
                    fired = (hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 or hlt.IsoMu24 or hlt.IsoMu27 or hlt.IsoMu30 or hlt.Mu50)

        elif (self.era == "2018"):
            if self.isMC:
                fired = (hlt.IsoMu24 or hlt.IsoMu27 or hlt.IsoMu30 or hlt.Mu50 or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8)
            else:
                if "Run2018A" == self.period or"Run2018B" == self.period:
                    fired = (hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8 or hlt.IsoMu24 or hlt.IsoMu27 or hlt.IsoMu30 or hlt.Mu50)
                else:
                    fired = (hlt.IsoMu24 or hlt.IsoMu27 or hlt.IsoMu30 or hlt.Mu50 or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8 or hlt.Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8)
        return fired

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
        hlt = Object(event, "HLT")

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
        lep_category = -1
        event_category = -1


        muons.sort(key=lambda muon: muon.pt, reverse=True)
        electrons.sort(key=lambda el: el.pt, reverse=True)
        # Choose loose/medium-quality e/mu for event categorization
        for idx,mu in enumerate(muons):
#            isoLep   = mu.pfRelIso04_all
            pass_ips = abs(mu.dxy) < 0.02 and abs(mu.dz) < 0.1
            pass_fid = abs(mu.eta) < 2.4 and mu.pt >= (20 if idx==0 else 10)
            pass_ids = mu.mediumId #and isoLep <= 0.25
            if pass_fid and pass_ids and pass_ips:
                good_muons.append(mu)
        for idy,el in enumerate(electrons):
            id_CB = el.cutBased
            # changing to MVA based ID :
            if el.pt >= (25 if idy==0 else 15) and abs(el.eta) <= 2.5 and self.electron_id_noiso(el, "90"):
                good_electrons.append(el)

        # let sort the muons in pt
        good_muons.sort(key=lambda x: x.pt, reverse=True)
        good_electrons.sort(key=lambda x: x.pt, reverse=True)

        # Find any remaining e/mu that pass looser selection
        extra_leptons = []
        for mu in muons:
            isoLep   = mu.pfRelIso04_all
            pass_ids = mu.looseId and isoLep <= 0.25
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
        _lead_lep_phi = good_leptons[0].phi if len(good_leptons) else 0.0
        _trail_lep_pt = good_leptons[1].pt if len(good_leptons) >= 2 else 0.0
        _trail_lep_eta = good_leptons[1].eta if len(good_leptons) >= 2 else 0.0
        _trail_lep_phi = good_leptons[1].phi if len(good_leptons) >= 2 else 0.0
        _lead_lep_pdgId = good_leptons[0].pdgId if len(good_leptons)  else 0.0
        _trail_lep_pdgId = good_leptons[1].pdgId if len(good_leptons) >= 2 else 0.0


	_leading_lep_flavor = 0
	if len(good_muons) and len(good_electrons):
		if good_muons[0].pt > good_electrons[0].pt: _leading_lep_flavor = 1

        self.out.fillBranch("leading_lep_pt{}".format(self.syst_suffix), _lead_lep_pt)
        self.out.fillBranch("leading_lep_eta{}".format(self.syst_suffix), _lead_lep_eta)
        self.out.fillBranch("leading_lep_phi{}".format(self.syst_suffix), _lead_lep_phi)    
        self.out.fillBranch("trailing_lep_pt{}".format(self.syst_suffix), _trail_lep_pt)
        self.out.fillBranch("trailing_lep_eta{}".format(self.syst_suffix), _trail_lep_eta)
        self.out.fillBranch("trailing_lep_phi{}".format(self.syst_suffix), _trail_lep_phi)
        self.out.fillBranch("leading_lep_flavor{}".format(self.syst_suffix), _leading_lep_flavor)

        ngood_leptons = len(good_leptons)
        nextra_leptons = len(extra_leptons)
        ngood_muons = len(good_muons)
        ngood_electrons = len(good_electrons)

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


        lep_category = 0
        event_category = 0
        if ngood_leptons < 2:
            lep_category = -1
            event_category = -1

        if ngood_muons == 2 and ngood_electrons == 0 and nextra_leptons == 0:
            lep_category = 1 #dimuon channel
            isoLep0   = good_muons[0].pfRelIso04_all
            isoLep1   = good_muons[1].pfRelIso04_all
            if (good_muons[0].pdgId * good_muons[1].pdgId == -13*13): # opposite sign: A or B region
                if (isoLep0 < 0.25 and isoLep1 < 0.25):
                    event_category = 1
                else:
                    event_category = 2
            else: #same sign: C or D region
                if (isoLep0 < 0.25 and isoLep1 < 0.25):
                    event_category = 3
                else:
                    event_category = 4
            if tk.deltaR(good_muons[0].eta, good_muons[0].phi,good_muons[1].eta, good_muons[1].phi,) < 0.3:
                event_category +=10 
        elif ngood_electrons == 2 and ngood_muons == 0 and nextra_leptons == 0:
            lep_category = 2 #dielectron channel
            isoLep0   = good_electrons[0].pfRelIso03_all
            isoLep1   = good_electrons[1].pfRelIso03_all
            if (good_electrons[0].pdgId * good_electrons[1].pdgId == -11*11): # opposite sign: A or B region
                if (isoLep0 < 0.15 and isoLep1 < 0.15):
                    event_category = 1
                else:
                    event_category = 2
            else: #same sign: C or D region
                if (isoLep0 < 0.15 and isoLep1 < 0.15):
                    event_category = 3
                else:
                    event_category = 4
            if tk.deltaR(good_electrons[0].eta, good_electrons[0].phi,good_electrons[1].eta, good_electrons[1].phi,) < 0.3:
                event_category +=10 


        self.out.fillBranch("ngood_leptons{}".format(self.syst_suffix), ngood_leptons)
        self.out.fillBranch("nextra_leptons{}".format(self.syst_suffix), nextra_leptons)
        self.out.fillBranch("lep_category{}".format(self.syst_suffix), lep_category)
        self.out.fillBranch("event_category{}".format(self.syst_suffix), event_category)

        DYcand =  ROOT.TLorentzVector()

        if (ngood_muons == 2 and ngood_electrons == 0 and nextra_leptons == 0):
            DYcand = good_muons[0].p4() + good_muons[1].p4()        
            self.out.fillBranch("DY_mass{}".format(self.syst_suffix),DYcand.M())
            self.out.fillBranch("DY_pt{}".format(self.syst_suffix),DYcand.Pt())

        if (ngood_electrons == 2 and ngood_muons == 0 and nextra_leptons == 0):
            DYcand = good_electrons[0].p4() + good_electrons[1].p4()
            self.out.fillBranch("DY_mass{}".format(self.syst_suffix),DYcand.M())
            self.out.fillBranch("DY_pt{}".format(self.syst_suffix),DYcand.Pt())

        # lepton trigger fired?
        ele_trigger = self.pass_electron_trigger(hlt)
        mu_trigger = self.pass_muon_trigger(hlt)

        self.out.fillBranch("electron_trigger{}".format(self.syst_suffix),ele_trigger)
        self.out.fillBranch("muon_trigger{}".format(self.syst_suffix),mu_trigger)


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

        # Let remove the negative categories with no obvious meaning meaning
        # This will reduce the size of most of the bacground and data
#        print "good leptons: ",len(good_leptons)
#        print "nextra_leptons: ",nextra_leptons
#        print "good mu: ",ngood_muons
#        print 'good el: ',ngood_electrons
        if (lep_category > 0 and event_category > 0  and len(good_leptons) == 2 and nextra_leptons == 0 and len(had_taus)==0): 
            return True
        else:
            return False
