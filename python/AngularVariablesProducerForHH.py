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

class AngularVariablesProducerForHH(Module):
    def __init__(self, do_syst=False, syst_var='',isMC=False):
        self.do_syst = do_syst
        self.syst_var = syst_var
        self.isMC = isMC
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

        #DeltaR variables
#        self.out.branch("dR_l1l2{}".format(self.syst_suffix), "F")
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

        # Helicity angles in the Collins-Sopper frame
        self.out.branch("cosThetaCS{}".format(self.syst_suffix),"F")
        self.out.branch("cosThetabHbb{}".format(self.syst_suffix),"F")
        self.out.branch("cosThetaZjjHzz{}".format(self.syst_suffix),"F")
        self.out.branch("cosThetaZllHzz{}".format(self.syst_suffix),"F")
        self.out.branch("phi1{}".format(self.syst_suffix),"F")
        self.out.branch("phi1_Zjj{}".format(self.syst_suffix),"F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getCosThetaStar_CS(self, h1, h2):
        """
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

        return math.cos(   CSaxis.Angle( h1.Vect().Unit() )    )
        """
        hh_lor = h1 + h2;
        hh = h1 + h2;

        h1_lor = h1;
        h_1 = h1

        h_1.Boost(-hh.BoostVector());

        return abs(h_1.CosTheta());

    def HelicityCosTheta(self, Booster, Boosted) :
        BoostVector = ROOT.TVector3( Booster.BoostVector() )
        Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() )
        return Boosted.CosTheta()

    def CosThetaAngles(self, Hj1, Hj2, Zj1, Zj2, ell1, ell2) :
        bb, zz, diHiggsCandidate = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()
        zz = Zj1 + Zj2 + ell1 + ell2
        bb = Hj1 + Hj2
        diHiggsCandidate = Hj1 + Hj2 + Zj1 + Zj2 + ell1 + ell2

        helicityThetas = []
        BoostedHgg, HHforBoost = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        HHforBoost.SetPtEtaPhiE(diHiggsCandidate.Pt(), diHiggsCandidate.Eta(), diHiggsCandidate.Phi(), diHiggsCandidate.Energy())
        BoostedHgg.SetPtEtaPhiE(zz.Pt(), zz.Eta(), zz.Phi(), zz.Energy())
        helicityThetas.append( self.HelicityCosTheta(HHforBoost, BoostedHgg) ) # CosThetaStar

        BoostedLeadingJet, HbbforBoost = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        HbbforBoost.SetPtEtaPhiE(bb.Pt(), bb.Eta(), bb.Phi(), bb.Energy())
        if (Hj1.Pt() >= Hj2.Pt()) : # is this leading jet ?
            BoostedLeadingJet.SetPtEtaPhiE(Hj1.Pt(), Hj1.Eta(), Hj1.Phi(), Hj1.Energy())
        else :
            BoostedLeadingJet.SetPtEtaPhiE(Hj2.Pt(), Hj2.Eta(), Hj2.Phi(), Hj2.Energy())
        helicityThetas.append( self.HelicityCosTheta(HbbforBoost, BoostedLeadingJet) ) # CosTheta_hbb


        BoostedZjj, BoostedZuu, BoostedLeadingZj, BoostedLeadingZu, HzzforBoost = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()
        HzzforBoost.SetPtEtaPhiE(zz.Pt(), zz.Eta(), zz.Phi(), zz.Energy())

        BoostedZjj.SetPtEtaPhiE((Zj1+Zj2).Pt(), (Zj1+Zj2).Eta(), (Zj1+Zj2).Phi(), (Zj1+Zj2).Energy())
        #if v == '': print 'before:', BoostedZjj.Pt(), BoostedZjj.Eta()
        helicityThetas.append( self.HelicityCosTheta(HzzforBoost, BoostedZjj) ) # CosTheta_zjj_hzz
        #if v == '': print 'after :', BoostedZjj.Pt(), BoostedZjj.Eta()

        BoostedZuu.SetPtEtaPhiE((ell1 + ell2).Pt(), (ell1 + ell2).Eta(), (ell1 + ell2).Phi(), (ell1 + ell2).Energy())
        helicityThetas.append( self.HelicityCosTheta(HzzforBoost, BoostedZuu) ) # CosTheta_zuu_hzz

        BoostedLeadingZj.SetPtEtaPhiE(Zj1.Pt(), Zj1.Eta(), Zj1.Phi(), Zj1.Energy())
        helicityThetas.append( self.HelicityCosTheta(HzzforBoost, BoostedLeadingZj) ) # CosTheta_zj1_hzz
        #if v == '': print 'CosTheta_zz:', helicityThetas[2], 'CosTheta_zz_test', helicityThetas[3]

        BoostedLeadingZu.SetPtEtaPhiE(ell1.Pt(), ell1.Eta(), ell1.Phi(), ell1.Energy())
        helicityThetas.append( self.HelicityCosTheta(HzzforBoost, BoostedLeadingZu) ) # CosTheta_zu1_hzz

        return helicityThetas

    def CosThetaAngles_ZZ(self,Zj1, Zj2, ell1, ell2) :
        zz = ROOT.TLorentzVector()
        zz = Zj1 + Zj2 + ell1 + ell2
        zjj = Zj1 + Zj2
        zellell = ell1 + ell2
        
        helicityThetas = []
        
        ZZforBoost, BoostedZuu, BoostedZjj = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()
        ZZforBoost.SetPtEtaPhiE(zz.Pt(), zz.Eta(), zz.Phi(), zz.Energy())
        
        BoostedZuu.SetPtEtaPhiE(zellell.Pt(), zellell.Eta(), zellell.Phi(), zellell.Energy())
        helicityThetas.append( HelicityCosTheta(ZZforBoost, BoostedZuu) ) # CosThetaStar_Zuu
        
        BoostedZjj.SetPtEtaPhiE(zjj.Pt(), zjj.Eta(), zjj.Phi(), zjj.Energy())
        helicityThetas.append( HelicityCosTheta(ZZforBoost, BoostedZjj) ) # CosThetaStar_Zjj
        
        return helicityThetas

    def CosThetaAngles_ZZ(self,Zj1, Zj2, ell1, ell2) :
        zz = ROOT.TLorentzVector()
        zz = Zj1 + Zj2 + ell1 + ell2
        zjj = Zj1 + Zj2
        zellell = ell1 + ell2

        helicityThetas = []

        ZZforBoost, BoostedZuu, BoostedZjj = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()
        ZZforBoost.SetPtEtaPhiE(zz.Pt(), zz.Eta(), zz.Phi(), zz.Energy())

        BoostedZuu.SetPtEtaPhiE(zellell.Pt(), zellell.Eta(), zellell.Phi(), zellell.Energy())
        helicityThetas.append( HelicityCosTheta(ZZforBoost, BoostedZuu) ) # CosThetaStar_Zuu

        BoostedZjj.SetPtEtaPhiE(zjj.Pt(), zjj.Eta(), zjj.Phi(), zjj.Energy())
        helicityThetas.append( HelicityCosTheta(ZZforBoost, BoostedZjj) ) # CosThetaStar_Zjj

        return helicityThetas

    def norm_planes_hi(self,partons, dihiggs) :
        boost_H = ROOT.TVector3( - dihiggs.BoostVector() )

        partons3v = []
        for i in range(len(partons)) :
            partons_i_boosted = ROOT.TLorentzVector()
            partons_i_boosted = partons[i]
            partons_i_boosted.Boost(boost_H)
            partons3v.append( partons_i_boosted.Vect().Unit() )

        vnorm = []
        R = ROOT.TRandom()
        for i in range(2) :
            rndm = R.Uniform(1)
            #if v == '': print 'i:', i, 'rndm', rndm
            if (rndm > 0.5) :
                vnorm.append( (partons3v[i*2].Cross(partons3v[i*2+1])).Unit() )
            else :
                vnorm.append( -1*(partons3v[i*2].Cross(partons3v[i*2+1])).Unit() )

        return vnorm

    def getPhi(self,Hj1, Hj2, Zj1, Zj2, ell1, ell2) :
        vPhi = []
        if (Hj1.Pt()<=0. or Hj2.Pt()<=0. or Zj1.Pt()<=0. or Zj2.Pt()<=0. or ell1.Pt()<=0. or ell2.Pt()<=0.) :
            vPhi.append(-3.5)
            vPhi.append(-3.5)
        else :
            zz, diHiggsCandidate = ROOT.TLorentzVector(), ROOT.TLorentzVector()
            zz = Zj1 + Zj2 + ell1 + ell2
            diHiggsCandidate = Hj1 + Hj2 + Zj1 + Zj2 + ell1 + ell2

            partons = []
            partons.append(Zj1 + Zj2) #leadingLepton
            partons.append(ell1 + ell2) #subleadingLepton
            partons.append(Hj1) #leadingJet
            partons.append(Hj2) #subleadingJet

            # Define hzz direction
            hzz = ROOT.TLorentzVector()
            hzz = Zj1 + Zj2 + ell1 + ell2

            boost_H = ROOT.TVector3( - diHiggsCandidate.BoostVector() )
            hzz.Boost(boost_H);

            hzz_vect = ROOT.TVector3( hzz.Vect().Unit() ) # hzz_vect has been boosted
            #hzz_vect_nobst = TVector3( zz.Vect().Unit() ) # hzz_vect_nobst has not been boosted

            #if v == '': print 'hzz_vect       :', hzz_vect.Pt(), hzz_vect.Eta(), hzz_vect.Phi()
            #if v == '': print 'hzz_vect_nobst :', hzz_vect_nobst.Pt(), hzz_vect_nobst.Eta(), hzz_vect_nobst.Phi()

            # Calculate the normal to Hzz and hbb decay plane
            vnorm = self.norm_planes_hi(partons, diHiggsCandidate)

            # Calculate Phi
            dsignhgg = hzz_vect.Dot(vnorm[1].Cross(vnorm[0]))/(abs(hzz_vect.Dot(vnorm[1].Cross(vnorm[0])))) # hzz_vect here is not important
            vPhi.append( dsignhgg*(-1) * math.acos(vnorm[0].Dot(vnorm[1])) )

            # Define z direction
            p1 = ROOT.TLorentzVector()
            p1.SetPxPyPzE(0, 0,  6500, 6500)
            z_vect = ROOT.TVector3( p1.Vect().Unit() )

            # Calcuate the normal to Hzz and z-direction plane
            zzprime = ROOT.TVector3((z_vect.Cross(hzz_vect)).Unit())
            #zzprime = TVector3((z_vect.Cross(hzz_vect_nobst)).Unit())

            # Calculate Phi1
            dsignhgg2 = hzz_vect.Dot(zzprime.Cross(vnorm[0]))/(abs(hzz_vect.Dot(zzprime.Cross(vnorm[0]))))  # hzz_vect here is not important
            if -1<=zzprime.Dot(vnorm[0]) and zzprime.Dot(vnorm[0])<=1:
                vPhi.append( dsignhgg2 * math.acos(zzprime.Dot(vnorm[0])) )
            else :
                vPhi.append(-3.5)

        return vPhi

    def getPhi_ZZ(self,Zj1, Zj2, ell1, ell2) :
        vPhi = []
        if (Zj1.Pt()<=0. or Zj2.Pt()<=0. or ell1.Pt()<=0. or ell2.Pt()<=0.) :
            vPhi.append(-3.5)
            vPhi.append(-3.5)
            vPhi.append(-3.5)
        else :
            zz = ROOT.TLorentzVector()
            zz = Zj1 + Zj2 + ell1 + ell2

            partons = []
            partons.append(ell1) #leadingLepton
            partons.append(ell2) #subleadingLepton
            partons.append(Zj1) #leadingJet
            partons.append(Zj2) #subleadingJet

            # Define Zuu direction
            zuu = ROOT.TLorentzVector()
            zuu = ell1 + ell2

            boost_ZZ = ROOT.TVector3( - zz.BoostVector() )
            zuu.Boost(boost_ZZ);
            zuu_vect = ROOT.TVector3( zuu.Vect().Unit() ) # zuu_vect has been boosted

            # Calculate the normal to Zuu and Zjj decay plane
            vnorm = self.norm_planes_hi(partons, zz)

            # Calculate Phi
                    #fixme Morse division by zero very small amount of the time
            if abs(zuu_vect.Dot(vnorm[1].Cross(vnorm[0])))>0:
                dsignhgg = zuu_vect.Dot(vnorm[1].Cross(vnorm[0]))/(abs(zuu_vect.Dot(vnorm[1].Cross(vnorm[0])))) # zuu_vect here is NOT important
            else:
                dsignhgg = 0
            if -1<=vnorm[0].Dot(vnorm[1]) and vnorm[0].Dot(vnorm[1])<=1:
                vPhi.append( dsignhgg*(-1) * math.acos(vnorm[0].Dot(vnorm[1])) )
            else:
                vPhi.append(-3.5)

            # Define z direction
            p1 = ROOT.TLorentzVector()
            p1.SetPxPyPzE(0, 0,  6500, 6500)
            z_vect = ROOT.TVector3( p1.Vect().Unit() )

            # Calcuate the normal to Zuu and z-direction plane
            zz1prime = ROOT.TVector3((z_vect.Cross(zuu_vect)).Unit())  # zuu_vect here IS important

            # Calculate Phi1_zuu
                    #fixme Morse division by zero very small amount of the time
            if abs(zuu_vect.Dot(zz1prime.Cross(vnorm[0])))>0:
                dsignhgg2 = zuu_vect.Dot(zz1prime.Cross(vnorm[0]))/(abs(zuu_vect.Dot(zz1prime.Cross(vnorm[0]))))  # zuu_vect here is NOT important
            else:
                dsignhgg2 = 0
            #print dsignhgg2,zz1prime.Dot(vnorm[0])
            vPhi.append( dsignhgg2 * math.acos(round(zz1prime.Dot(vnorm[0]),6)) )


            # Define Zjj direction
            zjj = ROOT.TLorentzVector()
            zjj = Zj1 + Zj2

            zjj.Boost(boost_ZZ);
            zjj_vect = ROOT.TVector3( zjj.Vect().Unit() ) # zjj_vect has been boosted

            # Calcuate the normal to Zjj and z-direction plane
            zz2prime = ROOT.TVector3((z_vect.Cross(zjj_vect)).Unit())  # zjj_vect here IS important

            # Calculate Phi1_zjj
            if abs(zjj_vect.Dot(zz2prime.Cross(vnorm[1])))>0:
                dsignhgg2 = zjj_vect.Dot(zz2prime.Cross(vnorm[1]))/(abs(zjj_vect.Dot(zz2prime.Cross(vnorm[1]))))  # zjj_vect here is NOT important
            else:
                dsignhgg2 = 0
            vPhi.append( dsignhgg2 * math.acos(zz2prime.Dot(vnorm[1])) )

        return vPhi    




    def analyze(self, event):

        
        # Z->ll
        _lead_lepton_pt = float(getattr(event,"leading_lep_pt{}".format(self.syst_suffix)))
        _lead_lepton_eta = float(getattr(event,"leading_lep_eta{}".format(self.syst_suffix)))
        _lead_lepton_phi = float(getattr(event,"leading_lep_phi{}".format(self.syst_suffix)))
        _trail_lepton_pt = float(getattr(event,"trailing_lep_pt{}".format(self.syst_suffix)))
        _trail_lepton_eta = float(getattr(event,"trailing_lep_eta{}".format(self.syst_suffix)))
        _trail_lepton_phi = float(getattr(event,"trailing_lep_phi{}".format(self.syst_suffix)))

        # H->bb
        _lead_Hjet_pt = float(getattr(event,"leading_Hbb_pt{}".format(self.syst_suffix)))
        _lead_Hjet_eta = float(getattr(event,"leading_Hbb_eta{}".format(self.syst_suffix)))
        _lead_Hjet_phi = float(getattr(event,"leading_Hbb_phi{}".format(self.syst_suffix)))
        _trail_Hjet_pt = float(getattr(event,"trailing_Hbb_pt{}".format(self.syst_suffix)))
        _trail_Hjet_eta = float(getattr(event,"trailing_Hbb_eta{}".format(self.syst_suffix)))
        _trail_Hjet_phi = float(getattr(event,"trailing_Hbb_phi{}".format(self.syst_suffix)))
        # Z->qq
        _lead_Zjet_pt = float(getattr(event,"leading_jet_pt{}".format(self.syst_suffix)))
        _lead_Zjet_eta = float(getattr(event,"leading_jet_eta{}".format(self.syst_suffix)))
        _lead_Zjet_phi = float(getattr(event,"leading_jet_phi{}".format(self.syst_suffix)))
        _trail_Zjet_pt = float(getattr(event,"trailing_jet_pt{}".format(self.syst_suffix)))
        _trail_Zjet_eta = float(getattr(event,"trailing_jet_eta{}".format(self.syst_suffix)))
        _trail_Zjet_phi = float(getattr(event,"trailing_jet_phi{}".format(self.syst_suffix)))
            
        electrons = list(Collection(event, "Electron"))
        muons = list(Collection(event, "Muon"))
        jets = list(Collection(event, "Jet"))
            
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

        #        print 'test: ',_lead_lepton_eta,_lead_lepton_phi

        for muon in muons:
            #           print 'mu: ',muon.eta,muon.phi
            if tk.deltaR( _lead_lepton_eta,_lead_lepton_phi, muon.eta, muon.phi,) < 0.0001:
                HiggsZlep_candidate.append(muon)
                HiggsZlep_candidate[-1].pt = _lead_lepton_pt
            if tk.deltaR( _trail_lepton_eta,_trail_lepton_phi, muon.eta, muon.phi,) < 0.0001:
                HiggsZlep_candidate.append(muon)
                HiggsZlep_candidate[-1].pt = _trail_lepton_pt

        for ele in electrons:
            if tk.deltaR( _lead_lepton_eta,_lead_lepton_phi, ele.eta, ele.phi,) < 0.0001:
                HiggsZlep_candidate.append(ele)
                HiggsZlep_candidate[-1].pt = _lead_lepton_pt
            if tk.deltaR( _trail_lepton_eta,_trail_lepton_phi, ele.eta, ele.phi,) < 0.0001:
                HiggsZlep_candidate.append(ele)
                HiggsZlep_candidate[-1].pt = _trail_lepton_pt
      
        for jet in jets:
            if tk.deltaR(_lead_Zjet_eta,_lead_Zjet_phi,jet.eta,jet.phi,) < 0.0001:
                HiggsZjet_candidate.append(jet)
                HiggsZjet_candidate[-1].pt = _lead_Zjet_pt
            if tk.deltaR(_trail_Zjet_eta,_trail_Zjet_phi,jet.eta,jet.phi,)< 0.0001:
                HiggsZjet_candidate.append(jet)
                HiggsZjet_candidate[-1].pt = _trail_Zjet_pt
            if tk.deltaR(_lead_Hjet_eta,_lead_Hjet_phi,jet.eta,jet.phi,)< 0.0001:
                Higgsbb_candidate.append(jet)
                Higgsbb_candidate[-1].pt = _lead_Hjet_pt
            if tk.deltaR(_trail_Hjet_eta,_trail_Hjet_phi,jet.eta,jet.phi,)< 0.0001:
                Higgsbb_candidate.append(jet)
                Higgsbb_candidate[-1].pt = _trail_Hjet_pt

        Higgs_cand_0 = Higgsbb_candidate[0].p4() + Higgsbb_candidate[1].p4()
        Zlep_cand_0 = HiggsZlep_candidate[0].p4() + HiggsZlep_candidate[1].p4()       
        Zjet_cand_0 = HiggsZjet_candidate[0].p4() + HiggsZjet_candidate[1].p4()
        Higgs_cand_1 = Zlep_cand_0 + Zjet_cand_0


        #        dR_l1l2 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi,)
        dR_l1j1 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi,)
        dR_l1j2 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi,)
        dR_l1b1 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi,)
        dR_l1b2 = tk.deltaR(HiggsZlep_candidate[0].eta, HiggsZlep_candidate[0].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)
        dR_l2j1 = tk.deltaR(HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi, HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi,)
        dR_l2j2 = tk.deltaR(HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi, HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi,)
        dR_l2b1 = tk.deltaR(HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi, Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi,)
        dR_l2b2 = tk.deltaR(HiggsZlep_candidate[1].eta, HiggsZlep_candidate[1].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)
        dR_j1j2 = tk.deltaR(HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi, HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi,)
        dR_j1b1 = tk.deltaR(HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi, Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi,)
        dR_j1b2 = tk.deltaR(HiggsZjet_candidate[0].eta, HiggsZjet_candidate[0].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)
        dR_j2b1 = tk.deltaR(HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi, Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi,)
        dR_j2b2 = tk.deltaR(HiggsZjet_candidate[1].eta, HiggsZjet_candidate[1].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)
        dR_b1b2 = tk.deltaR(Higgsbb_candidate[0].eta, Higgsbb_candidate[0].phi, Higgsbb_candidate[1].eta, Higgsbb_candidate[1].phi,)

        cosThetaStar_CS = self.getCosThetaStar_CS(Higgs_cand_0,Higgs_cand_1)
        cosTheta = self.CosThetaAngles(Higgsbb_candidate[0].p4(),Higgsbb_candidate[1].p4(),HiggsZjet_candidate[0].p4(),HiggsZjet_candidate[1].p4(),HiggsZlep_candidate[0].p4(),HiggsZlep_candidate[1].p4())
        Phi = self.getPhi(Higgsbb_candidate[0].p4(),Higgsbb_candidate[1].p4(),HiggsZjet_candidate[0].p4(),HiggsZjet_candidate[1].p4(),HiggsZlep_candidate[0].p4(),HiggsZlep_candidate[1].p4())
        PhiZZ = self.getPhi_ZZ(HiggsZjet_candidate[0].p4(),HiggsZjet_candidate[1].p4(),HiggsZlep_candidate[0].p4(),HiggsZlep_candidate[1].p4())

        #        self.out.fillBranch("dR_l1l2{}".format(self.syst_suffix), dR_l1l2)
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
        self.out.fillBranch("cosThetaCS{}".format(self.syst_suffix), cosThetaStar_CS)
        self.out.fillBranch("cosThetabHbb{}".format(self.syst_suffix),cosTheta[1])        
        self.out.fillBranch("cosThetaZjjHzz{}".format(self.syst_suffix),cosTheta[2])
        self.out.fillBranch("cosThetaZllHzz{}".format(self.syst_suffix),cosTheta[3])
        self.out.fillBranch("phi1{}".format(self.syst_suffix),Phi[-1])
        self.out.fillBranch("phi1_Zjj{}".format(self.syst_suffix),PhiZZ[-1])

        return True



