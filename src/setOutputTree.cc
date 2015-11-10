#include "../interface/setOutputTree.h"

setOutputTree::setOutputTree(TTree* outTree){
  fTree = outTree;
  setBranches();
}

setOutputTree::~setOutputTree(){
  delete fTree;
}

void setOutputTree::initializeVariables()
{
  run=-999;
  event=-999;
  lumi=-999;
  njets=-999;
  nPV=-999;
  issignal=-999;
  wSampleWeight=-999;
  genWeight=1;
  totalEventWeight=-999;
  totalEventWeight_2=-999;
  totalEventWeight_3=-999;
  eff_and_pu_Weight=-999;
  eff_and_pu_Weight_2=-999;
  eff_and_pu_Weight_3=-999;
  pfMET=-999;
  pfMET_Phi=-999;
  nu_pz_type0=-999;
  nu_pz_type2=-999;
  nu_pz_run2=-999;
  nu_pz_run2_oth=-999;
  nu_pz_run2_type=-999;
  nu_pz_isre=1;
  l_pt=-999;
  l_eta=-999;
  l_phi=-999;
  l_e=-999;
  W_pt_gen=-999;
  W_pz_gen=-999;
  W_rap_gen=-999;
  genGravMass=-999;
  nu_pz_gen=-999;
  nu_pt_gen=-999;
  nu_phi_gen=-999;
  nu_eta_gen=-999;
  hadW_pt_gen=-999;
  hadW_eta_gen=-999;
  hadW_phi_gen=-999;
  hadW_e_gen=-999;
  hadW_m_gen=-999;
  lepW_pt_gen=-999;
  lepW_eta_gen=-999;
  lepW_phi_gen=-999;
  lepW_e_gen=-999;
  lepW_m_gen=-999;
  AK8_pt_gen=-999;
  AK8_eta_gen=-999;
  AK8_phi_gen=-999;
  AK8_e_gen=-999;
  AK8_pruned_mass_gen=-999;
  AK8_softdrop_mass_gen=-999;
  AK8_softdrop_pt_gen=-999;
  deltaR_lak8jet=-999;
  deltaphi_METak8jet=-999;
  deltaphi_Vak8jet=-999;
  v_pt=-999;
  v_eta=-999;
  v_phi=-999;
  v_mt=-999;
  mass_lvj_type0=-999;
  mass_lvj_type2=-999;
  mass_lvj_run2=-999;
  nBTagJet_loose=-999;
  nBTagJet_medium=-999;
  nBTagJet_tight=-999;
  ungroomed_jet_pt=-999;
  ungroomed_jet_eta=-999;
  ungroomed_jet_phi=-999;
  ungroomed_jet_e=-999;
  jet_mass_pr=-999;
  jet_mass_so=-999;
  jet_pt_so=-999;
  jet_mass_tr=-999;
  jet_mass_fi=-999;
  jet_tau2tau1=-999;
  ttb_ungroomed_jet_pt=-999;
  ttb_ungroomed_jet_eta=-999;
  ttb_ungroomed_jet_phi=-999;
  ttb_ungroomed_jet_e=-999;
  ttb_jet_mass_pr=-999;
  ttb_jet_mass_so=-999;
  ttb_jet_pt_so=-999;
  ttb_jet_mass_tr=-999;
  ttb_jet_mass_fi=-999;
  ttb_jet_tau2tau1=-999;
  ttb_deltaeta_lak8jet=-999;
  mass_leptonic_closerjet=-999;
  mass_ungroomedjet_closerjet=-999;
  AK8_closerjet_pt=-999;
  AK8_closerjet_eta=-999;
  AK8_closerjet_phi=-999;
  AK8_closerjet_e=-999;
  vbf_maxpt_j1_pt=-999;
  vbf_maxpt_j1_eta=-999;
  vbf_maxpt_j1_phi=-999;
  vbf_maxpt_j1_e=-999;
  vbf_maxpt_j1_bDiscriminatorCSV=-999;
  vbf_maxpt_j2_pt=-999;
  vbf_maxpt_j2_eta=-999;
  vbf_maxpt_j2_phi=-999;
  vbf_maxpt_j2_e=-999;
  vbf_maxpt_j2_bDiscriminatorCSV=-999;
  vbf_maxpt_jj_pt=-999;
  vbf_maxpt_jj_eta=-999;
  vbf_maxpt_jj_phi=-999;
  vbf_maxpt_jj_m=-999;
  jet2_pt=0;
  jet2_btag=0;
  jet3_pt=0;
  jet3_btag=0;
}

void setOutputTree::setBranches()
{
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("lumi",&lumi,"lumi/I");
  fTree->Branch("njets",&njets,"njets/I");
  fTree->Branch("nPV",&nPV,"nPV/I");
  fTree->Branch("issignal",&issignal,"issignal/I");
  fTree->Branch("wSampleWeight",&wSampleWeight,"wSampleWeight/F");
  fTree->Branch("genWeight",&genWeight,"genWeight/F");
  fTree->Branch("totalEventWeight",&totalEventWeight,"totalEventWeight/F");
  fTree->Branch("eff_and_pu_Weight",&eff_and_pu_Weight,"eff_and_pu_Weight/F");
  fTree->Branch("totalEventWeight_2",&totalEventWeight_2,"totalEventWeight_2/F");
  fTree->Branch("eff_and_pu_Weight_2",&eff_and_pu_Weight_2,"eff_and_pu_Weight_2/F");
  fTree->Branch("totalEventWeight_3",&totalEventWeight_3,"totalEventWeight_3/F");
  fTree->Branch("eff_and_pu_Weight_3",&eff_and_pu_Weight_3,"eff_and_pu_Weight_3/F");
  fTree->Branch("pfMET",&pfMET,"pfMET/F");
  fTree->Branch("pfMET_Phi",&pfMET_Phi,"pfMET_Phi/F");
  fTree->Branch("nu_pz_type0",&nu_pz_type0,"nu_pz_type0/F");
  fTree->Branch("nu_pz_type2",&nu_pz_type2,"nu_pz_type2/F");
  fTree->Branch("nu_pz_run2",&nu_pz_run2,"nu_pz_run2/F");
  fTree->Branch("nu_pz_run2_oth",&nu_pz_run2_oth,"nu_pz_run2_oth/F");
  fTree->Branch("nu_pz_run2_type",&nu_pz_run2_type,"nu_pz_run2_type/I");
  fTree->Branch("nu_pz_isre",&nu_pz_isre,"nu_pz_isre/I");
  fTree->Branch("l_pt",&l_pt,"l_pt/F");
  fTree->Branch("l_eta",&l_eta,"l_eta/F");
  fTree->Branch("l_phi",&l_phi,"l_phi/F");
  fTree->Branch("l_e",&l_e,"l_e/F");
  fTree->Branch("ungroomed_jet_pt",&ungroomed_jet_pt,"ungroomed_jet_pt/F");
  fTree->Branch("ungroomed_jet_eta",&ungroomed_jet_eta,"ungroomed_jet_eta/F");
  fTree->Branch("ungroomed_jet_phi",&ungroomed_jet_phi,"ungroomed_jet_phi/F");
  fTree->Branch("ungroomed_jet_e",&ungroomed_jet_e,"ungroomed_jet_e/F");
  fTree->Branch("jet_mass_pr",&jet_mass_pr,"jet_mass_pr");
  fTree->Branch("jet_mass_so",&jet_mass_so,"jet_mass_so");
  fTree->Branch("jet_pt_so",&jet_pt_so,"jet_pt_so");
  fTree->Branch("jet_mass_tr",&jet_mass_tr,"jet_mass_tr");
  fTree->Branch("jet_mass_fi",&jet_mass_fi,"jet_mass_fi");
  fTree->Branch("jet_tau2tau1",&jet_tau2tau1,"jet_tau2tau1");
  fTree->Branch("ttb_ungroomed_jet_pt",&ttb_ungroomed_jet_pt,"ttb_ungroomed_jet_pt/F");
  fTree->Branch("ttb_ungroomed_jet_eta",&ttb_ungroomed_jet_eta,"ttb_ungroomed_jet_eta/F");
  fTree->Branch("ttb_ungroomed_jet_phi",&ttb_ungroomed_jet_phi,"ttb_ungroomed_jet_phi/F");
  fTree->Branch("ttb_ungroomed_jet_e",&ttb_ungroomed_jet_e,"ttb_ungroomed_jet_e/F");
  fTree->Branch("ttb_jet_mass_pr",&ttb_jet_mass_pr,"ttb_jet_mass_pr");
  fTree->Branch("ttb_jet_mass_so",&ttb_jet_mass_so,"ttb_jet_mass_so");
  fTree->Branch("ttb_jet_pt_so",&ttb_jet_pt_so,"ttb_jet_pt_so");
  fTree->Branch("ttb_jet_mass_tr",&ttb_jet_mass_tr,"ttb_jet_mass_tr");
  fTree->Branch("ttb_jet_mass_fi",&ttb_jet_mass_fi,"ttb_jet_mass_fi");
  fTree->Branch("ttb_jet_tau2tau1",&ttb_jet_tau2tau1,"ttb_jet_tau2tau1");
  fTree->Branch("ttb_deltaeta_lak8jet",&ttb_deltaeta_lak8jet,"ttb_deltaeta_lak8jet/F");
  fTree->Branch("W_pt_gen",&W_pt_gen,"W_pt_gen");
  fTree->Branch("W_pz_gen",&W_pz_gen,"W_pz_gen");
  fTree->Branch("W_rap_gen",&W_rap_gen,"W_rap_gen");
  fTree->Branch("genGravMass",&genGravMass,"genGravMass");
  fTree->Branch("nu_pz_gen",&nu_pz_gen,"nu_pz_gen");
  fTree->Branch("nu_pt_gen",&nu_pt_gen,"nu_pt_gen");
  fTree->Branch("nu_phi_gen",&nu_phi_gen,"nu_phi_gen");
  fTree->Branch("nu_eta_gen",&nu_eta_gen,"nu_eta_gen");
  fTree->Branch("hadW_pt_gen",&hadW_pt_gen,"hadW_pt_gen");
  fTree->Branch("hadW_eta_gen",&hadW_eta_gen,"hadW_eta_gen");
  fTree->Branch("hadW_phi_gen",&hadW_phi_gen,"hadW_phi_gen");
  fTree->Branch("hadW_e_gen",&hadW_e_gen,"hadW_e_gen");
  fTree->Branch("hadW_m_gen",&hadW_m_gen,"hadW_m_gen");
  fTree->Branch("lepW_pt_gen",&lepW_pt_gen,"lepW_pt_gen");
  fTree->Branch("lepW_eta_gen",&lepW_eta_gen,"lepW_eta_gen");
  fTree->Branch("lepW_phi_gen",&lepW_phi_gen,"lepW_phi_gen");
  fTree->Branch("lepW_e_gen",&lepW_e_gen,"lepW_e_gen");
  fTree->Branch("lepW_m_gen",&lepW_m_gen,"lepW_m_gen");
  fTree->Branch("AK8_pt_gen",&AK8_pt_gen,"AK8_pt_gen");
  fTree->Branch("AK8_eta_gen",&AK8_eta_gen,"AK8_eta_gen");
  fTree->Branch("AK8_phi_gen",&AK8_phi_gen,"AK8_phi_gen");
  fTree->Branch("AK8_e_gen",&AK8_e_gen,"AK8_e_gen");
  fTree->Branch("AK8_pruned_mass_gen",&AK8_pruned_mass_gen,"AK8_pruned_mass_gen");
  fTree->Branch("AK8_softdrop_mass_gen",&AK8_softdrop_mass_gen,"AK8_softdrop_mass_gen");
  fTree->Branch("AK8_softdrop_pt_gen",&AK8_softdrop_pt_gen,"AK8_softdrop_pt_gen");
  fTree->Branch("deltaR_lak8jet",&deltaR_lak8jet,"deltaR_lak8jet/F");
  fTree->Branch("deltaphi_METak8jet",&deltaphi_METak8jet,"deltaphi_METak8jet/F");
  fTree->Branch("deltaphi_Vak8jet",&deltaphi_Vak8jet,"deltaphi_Vak8jet/F");
  fTree->Branch("v_pt",&v_pt,"v_pt/F");
  fTree->Branch("v_eta",&v_eta,"v_eta/F");
  fTree->Branch("v_phi",&v_phi,"v_phi/F");
  fTree->Branch("v_mt",&v_mt,"v_mt/F");
  fTree->Branch("mass_lvj_type0",&mass_lvj_type0,"mass_lvj_type0/F");
  fTree->Branch("mass_lvj_type2",&mass_lvj_type2,"mass_lvj_type2/F");
  fTree->Branch("mass_lvj_run2",&mass_lvj_run2,"mass_lvj_run2/F");
  fTree->Branch("nBTagJet_loose",&nBTagJet_loose,"nBTagJet_loose/I");
  fTree->Branch("nBTagJet_medium",&nBTagJet_medium,"nBTagJet_medium/I");
  fTree->Branch("nBTagJet_tight",&nBTagJet_tight,"nBTagJet_tight/I");
  fTree->Branch("mass_leptonic_closerjet",&mass_leptonic_closerjet,"mass_leptonic_closerjet/F");
  fTree->Branch("mass_ungroomedjet_closerjet",&mass_ungroomedjet_closerjet,"mass_ungroomedjet_closerjet/F");
  fTree->Branch("AK8_closerjet_pt",&AK8_closerjet_pt,"AK8_closerjet_pt/F");
  fTree->Branch("AK8_closerjet_eta",&AK8_closerjet_eta,"AK8_closerjet_eta/F");
  fTree->Branch("AK8_closerjet_phi",&AK8_closerjet_phi,"AK8_closerjet_phi/F");
  fTree->Branch("AK8_closerjet_e",&AK8_closerjet_e,"AK8_closerjet_e/F");
  fTree->Branch("vbf_maxpt_j1_pt",&vbf_maxpt_j1_pt,"vbf_maxpt_j1_pt/F");
  fTree->Branch("vbf_maxpt_j1_eta",&vbf_maxpt_j1_eta,"vbf_maxpt_j1_eta/F");
  fTree->Branch("vbf_maxpt_j1_phi",&vbf_maxpt_j1_phi,"vbf_maxpt_j1_phi/F");
  fTree->Branch("vbf_maxpt_j1_e",&vbf_maxpt_j1_e,"vbf_maxpt_j1_e/F");
  fTree->Branch("vbf_maxpt_j1_bDiscriminatorCSV",&vbf_maxpt_j1_bDiscriminatorCSV,"vbf_maxpt_j1_bDiscriminatorCSV/F");
  fTree->Branch("vbf_maxpt_j2_pt",&vbf_maxpt_j2_pt,"vbf_maxpt_j2_pt/F");
  fTree->Branch("vbf_maxpt_j2_eta",&vbf_maxpt_j2_eta,"vbf_maxpt_j2_eta/F");
  fTree->Branch("vbf_maxpt_j2_phi",&vbf_maxpt_j2_phi,"vbf_maxpt_j2_phi/F");
  fTree->Branch("vbf_maxpt_j2_e",&vbf_maxpt_j2_e,"vbf_maxpt_j2_e/F");
  fTree->Branch("vbf_maxpt_j2_bDiscriminatorCSV",&vbf_maxpt_j2_bDiscriminatorCSV,"vbf_maxpt_j2_bDiscriminatorCSV/F");
  fTree->Branch("vbf_maxpt_jj_pt",&vbf_maxpt_jj_pt,"vbf_maxpt_jj_pt/F");
  fTree->Branch("vbf_maxpt_jj_eta",&vbf_maxpt_jj_eta,"vbf_maxpt_jj_eta/F");
  fTree->Branch("vbf_maxpt_jj_phi",&vbf_maxpt_jj_phi,"vbf_maxpt_jj_phi/F");
  fTree->Branch("vbf_maxpt_jj_m",&vbf_maxpt_jj_m,"vbf_maxpt_jj_m/F");
  fTree->Branch("jet2_pt",&jet2_pt,"jet2_pt/F");
  fTree->Branch("jet2_btag",&jet2_btag,"jet2_btag/F");
  fTree->Branch("jet3_pt",&jet3_pt,"jet3_pt/F");
  fTree->Branch("jet3_btag",&jet3_btag,"jet3_btag/F");
}

