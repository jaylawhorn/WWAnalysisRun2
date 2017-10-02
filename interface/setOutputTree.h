#ifndef __setOutputTree__
#define __setOutputTree__

#include "TTree.h"
#include "TChain.h"
#include <vector>

class setOutputTree {

 public:

  TTree* fTree;
  std::vector<int>      LHEid;
  std::vector<float> LHEWeight;
  int run;
  int event;
  int nEvents;
  int nNegEvents;
  int lumi;
  int nPV;
  int issignal;
  int issignal_PuppiAK8;
  int issignal_AK4jetjet;
  int issignal_PuppiAK4jetjet;
  float wSampleWeight;
  float genWeight;
  float top1_NNLO_Weight;
  float top2_NNLO_Weight;
  float trig_eff_Weight;
  float id_eff_Weight;
  float gen_top1_pt;
  float gen_top2_pt;
  float totalEventWeight;
  float totalEventWeight_2;
  float totalEventWeight_3;
  float pu_Weight;
  float pu_Weight_up;
  float pu_Weight_down;
  float pfMET;
  float pfMET_jes_up;
  float pfMET_jes_dn;
  float pfMET_jer;
  float pfMET_jer_up;
  float pfMET_jer_dn;
  float pfMET_Phi;
  float pfMET_Corr; 
  float pfMET_Corr_phi; 
  float pfMET_Corr_Cov00; 
  float pfMET_Corr_Cov01; 
  float pfMET_Corr_Cov11; 
  float pfMET_Corr_jerup; 
  float pfMET_Corr_jerdn; 
  float pfMET_Corr_jenup; 
  float pfMET_Corr_jendn; 
  float pfMET_Corr_uncup; 
  float pfMET_Corr_uncdn; 
  float pfMET_Corr_jrsup; 
  float pfMET_Corr_jrsdn; 
  float pfMET_Corr_phijerup; 
  float pfMET_Corr_phijerdn; 
  float pfMET_Corr_phijenup; 
  float pfMET_Corr_phijendn; 
  float pfMET_Corr_phiuncup; 
  float pfMET_Corr_phiuncdn; 
  float pfMET_Corr_phijrsup; 
  float pfMET_Corr_phijrsdn; 
  float pfMETpuppi;
  float pfMETpuppi_jes_up;
  float pfMETpuppi_jes_dn;
  float pfMETpuppi_jer;
  float pfMETpuppi_jer_up;
  float pfMETpuppi_jer_dn;
  float pfMETpuppi_Phi;
  float pfMETpuppi_Cov00; 
  float pfMETpuppi_Cov01; 
  float pfMETpuppi_Cov11; 
  float pfMETpuppi_Corr; 
  float pfMETpuppi_Corr_phi; 
  float pfMETpuppi_Corr_Cov00; 
  float pfMETpuppi_Corr_Cov01; 
  float pfMETpuppi_Corr_Cov11; 
  float nu_pz_type0;
  float nu_pz_type2;
  float nu_pz_run2;
  float nu_pz_run2_oth;
  int nu_pz_run2_type;
  int nu_pz_isre;
  int type;
  float l_pt1;
  float l_eta1;
  float l_phi1;
  float l_e1;
  float l_charge1;
  float l_pt2;
  float l_eta2;
  float l_phi2;
  float l_e2;
  float l_charge2;
  float dilep_pt;
  float dilep_eta;
  float dilep_phi;
  float dilep_m;
  float ungroomed_AK8jet_pt;
  float ungroomed_AK8jet_pt_jes_up;
  float ungroomed_AK8jet_pt_jes_dn; 
  float ungroomed_AK8jet_pt_jer;
  float ungroomed_AK8jet_pt_jer_up;
  float ungroomed_AK8jet_pt_jer_dn;
  float ungroomed_AK8jet_eta;
  float ungroomed_AK8jet_phi;
  float ungroomed_AK8jet_e;
  float ungroomed_AK8jet_charge;
  float AK8jet_mass;
  float AK8jet_mass_pr;
  float AK8jet_mass_pr_jes_up;
  float AK8jet_mass_pr_jes_dn;
  float AK8jet_mass_pr_jer;
  float AK8jet_mass_pr_jer_up;
  float AK8jet_mass_pr_jer_dn;
  float AK8jet_mass_so;
  float AK8jet_pt_so;
  float AK8jet_mass_tr;
  float AK8jet_mass_fi;
  float AK8jet_tau2tau1;
  int   AK8_jetID_loose;
  float AK4_jetjet_pt;
  float AK4_jetjet_mass;
  float AK4_jetjet_deltaeta;
  float AK4_jetjet_deltaphi;
  float AK4_jetjet_deltar;
  float PuppiAK4_jetjet_pt;
  float PuppiAK4_jetjet_mass;
  float PuppiAK4_jetjet_deltaeta;
  float PuppiAK4_jetjet_deltaphi;
  float PuppiAK4_jetjet_deltar;
  float ttb_ungroomed_jet_pt;
  float ttb_ungroomed_jet_eta;
  float ttb_ungroomed_jet_phi;
  float ttb_ungroomed_jet_e;
  float ttb_jet_mass_pr;
  float ttb_jet_mass_so;
  float ttb_jet_pt_so;
  float ttb_jet_mass_tr;
  float ttb_jet_mass_fi;
  float ttb_jet_tau2tau1;
  float ttb_deltaeta_lak8jet;
  float ungroomed_PuppiAK8_jet_pt;
  float ungroomed_PuppiAK8_jet_pt_jes_up;
  float ungroomed_PuppiAK8_jet_pt_jes_dn; 
  float ungroomed_PuppiAK8_jet_pt_jer;
  float ungroomed_PuppiAK8_jet_pt_jer_up;
  float ungroomed_PuppiAK8_jet_pt_jer_dn;
  float ungroomed_PuppiAK8_jet_eta;
  float ungroomed_PuppiAK8_jet_phi;
  float ungroomed_PuppiAK8_jet_e;
  float ungroomed_PuppiAK8_jet_charge;
  float PuppiAK8_jet_mass;
  float PuppiAK8_jet_mass_pr;
  float PuppiAK8_jet_mass_pr_jes_up;
  float PuppiAK8_jet_mass_pr_jes_dn;
  float PuppiAK8_jet_mass_pr_jer;
  float PuppiAK8_jet_mass_pr_jer_up;
  float PuppiAK8_jet_mass_pr_jer_dn;
  float PuppiAK8_jet_mass_so;
  float PuppiAK8_jet_pt_so;
  float PuppiAK8_jet_mass_tr;
  float PuppiAK8_jet_mass_fi;
  float PuppiAK8_jet_tau2tau1;
  float PuppiAK8_jet_sj1_pt;
  float PuppiAK8_jet_sj1_eta;
  float PuppiAK8_jet_sj1_phi;
  float PuppiAK8_jet_sj1_m;
  float PuppiAK8_jet_sj1_q;
  float PuppiAK8_jet_sj2_pt;
  float PuppiAK8_jet_sj2_eta;
  float PuppiAK8_jet_sj2_phi;
  float PuppiAK8_jet_sj2_m;
  float PuppiAK8_jet_sj2_q;
  int  PuppiAK8_jetID_loose;
  float AK4_jet1_pt;
  float AK4_jet1_pt_jes_up;
  float AK4_jet1_pt_jes_dn;
  float AK4_jet1_pt_jer;
  float AK4_jet1_pt_jer_up;
  float AK4_jet1_pt_jer_dn;
  float AK4_jet1_eta;
  float AK4_jet1_phi;
  float AK4_jet1_e;
  float AK4_jet2_pt;
  float AK4_jet2_pt_jes_up;
  float AK4_jet2_pt_jes_dn;
  float AK4_jet2_pt_jer;
  float AK4_jet2_pt_jer_up;
  float AK4_jet2_pt_jer_dn;
  float AK4_jet2_eta;
  float AK4_jet2_phi;
  float AK4_jet2_e;
  float PuppiAK4_jet1_pt;
  float PuppiAK4_jet1_pt_jes_up;
  float PuppiAK4_jet1_pt_jes_dn;
  float PuppiAK4_jet1_pt_jer;
  float PuppiAK4_jet1_pt_jer_up;
  float PuppiAK4_jet1_pt_jer_dn;
  float PuppiAK4_jet1_eta;
  float PuppiAK4_jet1_phi;
  float PuppiAK4_jet1_e;
  float PuppiAK4_jet2_pt;
  float PuppiAK4_jet2_pt_jes_up;
  float PuppiAK4_jet2_pt_jes_dn;
  float PuppiAK4_jet2_pt_jer;
  float PuppiAK4_jet2_pt_jer_up;
  float PuppiAK4_jet2_pt_jer_dn;
  float PuppiAK4_jet2_eta;
  float PuppiAK4_jet2_phi;
  float PuppiAK4_jet2_e;
  int   isGen;
  float lep_pt_gen;
  float lep_eta_gen;
  float W_pt_gen;
  float W_pz_gen;
  float W_rap_gen;
  float nu_pz_gen;
  float nu_pt_gen;
  float nu_phi_gen;
  float nu_eta_gen;
  float hadW_pt_gen;
  float hadW_eta_gen;
  float hadW_phi_gen;
  float hadW_e_gen;
  float hadW_m_gen;
  float lepW_pt_gen;
  float lepW_eta_gen;
  float lepW_phi_gen;
  float lepW_e_gen;
  float lepW_m_gen;
  float WW_mass_gen;
  float WW_mT_gen;
  float WW_pT_gen;
  float AK8_pt_gen;
  float AK8_eta_gen;
  float AK8_phi_gen;
  float AK8_e_gen;
  float AK8_mass_gen;
  float AK8_pruned_mass_gen;
  float AK8_softdrop_mass_gen;
  float AK8_softdrop_pt_gen;
  float AK4_1_pt_gen;
  float AK4_1_eta_gen;
  float AK4_1_phi_gen;
  float AK4_1_e_gen;
  float AK4_1_mass_gen;
  float AK4_2_pt_gen;
  float AK4_2_eta_gen;
  float AK4_2_phi_gen;
  float AK4_2_e_gen;
  float AK4_2_mass_gen;
  float AK4_jj_DeltaEta_gen;
  float AK4_jj_mass_gen;
  float AK4_DR_GENRECO_11;
  float AK4_DR_GENRECO_12;
  float AK4_DR_GENRECO_21;
  float AK4_DR_GENRECO_22;
  float AK4Puppi_DR_GENRECO_11;
  float AK4Puppi_DR_GENRECO_12;
  float AK4Puppi_DR_GENRECO_21;
  float AK4Puppi_DR_GENRECO_22;
  float deltaR_Wjet_GenReco;
  float deltaR_lak8jet;
  float deltaphi_METak8jet;
  float deltaphi_Vak8jet;
  float deltaR_lPuppiak8jet;
  float deltaphi_METPuppiak8jet;
  float deltaphi_VPuppiak8jet;
  float deltaR_lak4jetjet;
  float deltaphi_METak4jetjet;
  float deltaphi_Vak4jetjet;
  float deltaR_lPuppiak4jetjet;
  float deltaphi_METPuppiak4jetjet;
  float deltaphi_VPuppiak4jetjet;
  float v_pt_type2;
  float v_pt_type0;
  float v_pt_run2;
  float v_eta_type2;
  float v_eta_type0;
  float v_eta_run2;
  float v_phi;
  float v_mt_type2;
  float v_mt_type0;
  float v_mt_run2;
  float v_mass_type2;
  float v_mass_type0;
  float v_mass_run2;
  float v_puppi_pt_type2;
  float v_puppi_pt_type0;
  float v_puppi_pt_run2;
  float v_puppi_eta_type2;
  float v_puppi_eta_type0;
  float v_puppi_eta_run2;
  float v_puppi_phi;
  float v_puppi_mt_type2;
  float v_puppi_mt_type0;
  float v_puppi_mt_run2;
  float v_puppi_mass_type2;
  float v_puppi_mass_type0;
  float v_puppi_mass_run2;
  float mass_lvj_type0;
  float mass_lvj_type0_met_jes_up;
  float mass_lvj_type0_met_jes_dn;
  float mass_lvj_type0_met_jer;
  float mass_lvj_type0_met_jer_up;
  float mass_lvj_type0_met_jer_dn;
  float mass_lvj_type0_met_PuppiAK8_jes_up;
  float mass_lvj_type0_met_PuppiAK8_jes_dn;
  float mass_lvj_type2;
  float mass_lvj_run2;
  float mass_lvj_type0_PuppiAK8;
  float mass_lvj_type2_PuppiAK8;
  float mass_lvj_run2_PuppiAK8;
  float mt_lvj_type0_PuppiAK8;
  float mt_lvj_type2_PuppiAK8;
  float mt_lvj_run2_PuppiAK8;
  float pt_lvj_type0_PuppiAK8;
  float pt_lvj_type2_PuppiAK8;
  float pt_lvj_run2_PuppiAK8;
  float eta_lvj_type0_PuppiAK8;
  float eta_lvj_type2_PuppiAK8;
  float eta_lvj_run2_PuppiAK8;
  float rapidity_lvj_type0_PuppiAK8;
  float rapidity_lvj_type2_PuppiAK8;
  float rapidity_lvj_run2_PuppiAK8;
  float phi_lvj_type0_PuppiAK8;
  float phi_lvj_type2_PuppiAK8;
  float phi_lvj_run2_PuppiAK8;
  float energy_lvj_type0_PuppiAK8;
  float energy_lvj_type2_PuppiAK8;
  float energy_lvj_run2_PuppiAK8;
  float mass_lvjj_type0_AK4;
  float mass_lvjj_type0_met_jes_up_AK4;
  float mass_lvjj_type0_met_jes_dn_AK4;
  float mass_lvjj_type2_AK4;
  float mass_lvjj_run2_AK4;
  float mass_lvjj_type0_PuppiAK4;
  float mass_lvjj_type0_met_jes_up_PuppiAK4;
  float mass_lvjj_type0_met_jes_dn_PuppiAK4;
  float mass_lvjj_type2_PuppiAK4;
  float mass_lvjj_run2_PuppiAK4;
  float mass_llj_PuppiAK8;
  float mass_leptonic_closerjet;
  float mass_ungroomedjet_closerjet;
  float AK8_closerjet_pt;
  float AK8_closerjet_eta;
  float AK8_closerjet_phi;
  float AK8_closerjet_e;
  int njets;
  int njetsPuppi;
  int nGoodAK8jets;
  int njets_unmerged;
  int njetsPuppi_unmerged;
  int nBTagJet_loose;
  int nBTagJet_medium;
  int nBTagJet_tight;
  int nBTagJetPuppi_loose;
  int nBTagJetPuppi_medium;
  int nBTagJetPuppi_tight;
  int nBTagJet_loose_unmerged;
  int nBTagJet_medium_unmerged;
  int nBTagJet_tight_unmerged;
  int nBTagJetPuppi_loose_unmerged;
  int nBTagJetPuppi_medium_unmerged;
  int nBTagJetPuppi_tight_unmerged;
  float vbf_maxpt_j1_pt;
  float vbf_maxpt_j1_pt_jes_up;
  float vbf_maxpt_j1_pt_jes_dn;
  float vbf_maxpt_j1_pt_jer;
  float vbf_maxpt_j1_pt_jer_up;
  float vbf_maxpt_j1_pt_jer_dn;
  float vbf_maxpt_j1_eta;
  float vbf_maxpt_j1_eta_jes_up;
  float vbf_maxpt_j1_eta_jes_dn;
  float vbf_maxpt_j1_eta_jer;
  float vbf_maxpt_j1_eta_jer_up;
  float vbf_maxpt_j1_eta_jer_dn;
  float vbf_maxpt_j1_phi;
  float vbf_maxpt_j1_e;
  float vbf_maxpt_j1_bDiscriminatorCSV;
  float vbf_maxpt_j2_pt;
  float vbf_maxpt_j2_pt_jes_up;
  float vbf_maxpt_j2_pt_jes_dn;
  float vbf_maxpt_j2_pt_jer;
  float vbf_maxpt_j2_pt_jer_up;
  float vbf_maxpt_j2_pt_jer_dn;
  float vbf_maxpt_j2_eta;
  float vbf_maxpt_j2_eta_jes_up;
  float vbf_maxpt_j2_eta_jes_dn;
  float vbf_maxpt_j2_eta_jer;
  float vbf_maxpt_j2_eta_jer_up;
  float vbf_maxpt_j2_eta_jer_dn;
  float vbf_maxpt_j2_phi;
  float vbf_maxpt_j2_e;
  float vbf_maxpt_j2_bDiscriminatorCSV;
  int   vbfPuppi_maxpt_j1_ID;
  int   vbfPuppi_maxpt_j2_ID;
  float vbf_maxpt_jj_pt;
  float vbf_maxpt_jj_eta;
  float vbf_maxpt_jj_phi;
  float vbf_maxpt_jj_m;
  float vbf_maxpt_jj_Deta;
  float vbf_maxpt_deltaR;
  float vbfPuppi_maxpt_j1_pt;
  float vbfPuppi_maxpt_j1_pt_jes_up;
  float vbfPuppi_maxpt_j1_pt_jes_dn;
  float vbfPuppi_maxpt_j1_pt_jer;
  float vbfPuppi_maxpt_j1_pt_jer_up;
  float vbfPuppi_maxpt_j1_pt_jer_dn;
  float vbfPuppi_maxpt_j1_eta;
  float vbfPuppi_maxpt_j1_eta_jes_up;
  float vbfPuppi_maxpt_j1_eta_jes_dn;
  float vbfPuppi_maxpt_j1_eta_jer;
  float vbfPuppi_maxpt_j1_eta_jer_up;
  float vbfPuppi_maxpt_j1_eta_jer_dn;
  float vbfPuppi_maxpt_j1_phi;
  float vbfPuppi_maxpt_j1_e;
  float vbfPuppi_maxpt_j1_bDiscriminatorCSV;
  float vbfPuppi_maxpt_j2_pt;
  float vbfPuppi_maxpt_j2_pt_jes_up;
  float vbfPuppi_maxpt_j2_pt_jes_dn;
  float vbfPuppi_maxpt_j2_pt_jer;
  float vbfPuppi_maxpt_j2_pt_jer_up;
  float vbfPuppi_maxpt_j2_pt_jer_dn;
  float vbfPuppi_maxpt_j2_eta;
  float vbfPuppi_maxpt_j2_eta_jes_up;
  float vbfPuppi_maxpt_j2_eta_jes_dn;
  float vbfPuppi_maxpt_j2_eta_jer;
  float vbfPuppi_maxpt_j2_eta_jer_up;
  float vbfPuppi_maxpt_j2_eta_jer_dn;
  float vbfPuppi_maxpt_j2_phi;
  float vbfPuppi_maxpt_j2_e;
  float vbfPuppi_maxpt_j2_bDiscriminatorCSV;
  float vbfPuppi_maxpt_jj_pt;
  float vbfPuppi_maxpt_jj_eta;
  float vbfPuppi_maxpt_jj_phi;
  float vbfPuppi_maxpt_jj_m;
  float vbfPuppi_maxpt_jj_Deta;
  float vbfPuppi_maxpt_deltaR;
  float jet2_pt;
  float jet2_eta;
  float jet2_phi;
  float jet2_e;
  float jet2_btag;
  float jet3_pt;
  float jet3_eta;
  float jet3_phi;
  float jet3_e;
  float jet3_btag;
  float deltaR_AK8_closestBtagJet;
  float deltaR_AK8_closestBtagJet_loose;
//  float deltaR_AK4;  
  float costheta1_type0;
  float costheta2_type0;
  float phi_type0;
  float phi1_type0;
  float costhetastar_type0;
  float VBSCentrality_type0;
  float costheta1_type2;
  float costheta2_type2;
  float phi_type2;
  float phi1_type2;
  float costhetastar_type2;
  float VBSCentrality_type2;
  float costheta1_run2;
  float costheta2_run2;
  float phi_run2;
  float phi1_run2;
  float costhetastar_run2;
  float VBSCentrality_run2;
  float LepWEta;
  float LepWRapidity;
  float HadWEta;
  float HadWRapidity;
  float WWEta;
  float WWRapidity;
  float ZeppenfeldWH;
  float RpT_type0;
  float ZeppenfeldWL_type0;
  float LeptonProjection_type0;
  float RpT_type2;
  float ZeppenfeldWL_type2;
  float LeptonProjection_type2;
  float RpT_run2;
  float ZeppenfeldWL_run2;
  float LeptonProjection_run2;
  float PtBalance_type0;
  float PtBalance_type2;
  float PtBalance_run2;
  float BosonCentrality_type0;
  float BosonCentrality_type2;
  float BosonCentrality_run2;
  float BosonCentralityPuppi_type0;
  float BosonCentralityPuppi_type2;
  float BosonCentralityPuppi_run2;

  setOutputTree(TTree* outputTree);
  //  setOutputTree(TTree *outputTree=0);
  //  setOutputTree(TFile *outputFile=0, std::string outputTreeName="WWTree");
  ~setOutputTree();

  void initializeVariables();
  
  void setBranches();

};

#endif
