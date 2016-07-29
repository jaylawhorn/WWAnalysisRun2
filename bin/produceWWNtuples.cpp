#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"

#include "../interface/setInputTree.h"
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"
#include "../interface/readJSONFile.h"

using namespace std;

//*****PU WEIGHT***************

vector<double> generate_weights(TH1* data_npu_estimated, int isForSynch){
  // see SimGeneral/MixingModule/python/mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi.py; copy and paste from there:
  const double npu_probs[50] = { 0.000829312873542, 0.00124276120498,  0.00339329181587,  0.00408224735376, 0.00383036590008,
                                 0.00659159288946,  0.00816022734493,  0.00943640833116,  0.0137777376066,  0.017059392038,
                                 0.0213193035468,   0.0247343174676,   0.0280848773878,   0.0323308476564,  0.0370394341409,
                                 0.0456917721191,   0.0558762890594,   0.0576956187107,   0.0625325287017,  0.0591603758776,
                                 0.0656650815128,   0.0678329011676,   0.0625142146389,   0.0548068448797,  0.0503893295063,
                                 0.040209818868,    0.0374446988111,   0.0299661572042,   0.0272024759921,  0.0219328403791,
                                 0.0179586571619,   0.0142926728247,   0.00839941654725,  0.00522366397213, 0.00224457976761,
                                 0.000779274977993, 0.000197066585944, 7.16031761328e-05, 0.0,              0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.0, 0.0, 0.0, 0.0 };
  
  if (isForSynch==0) { //OFFICIAL RECIPE
    vector<double> result(50);
    double s = 0.0;
    for(int npu=0; npu<50; ++npu){
      double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
      result[npu] = npu_estimated / npu_probs[npu];
      s += npu_estimated;
    }
    // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    for(int npu=0; npu<50; ++npu){
      result[npu] /= s;
    }
    return result;
  }

  else { //THIS IS FOR THE SYNCH ONLY. THIS IS NOT THE OFFICIAL RECIPE!
    vector<double> result(50);
    for(int npu=0; npu<50; ++npu){
      if (data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu))==0.)
	result[npu] = 0.;
      else {
	double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
	result[npu] = npu_estimated;
      }
    }
    return result;
  }

}


//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  int isMC = atoi(argv[3]);
  std::string leptonName = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string numberOfEntries = argv[8];
  float genMass = atof(argv[9]);
  int applyTrigger = atoi(argv[10]);
  std::string jsonFileName = argv[11];
  int isLocal = atoi(argv[12]);
  
  float weight = std::atof(xSecWeight.c_str())/std::atof(numberOfEntries.c_str());
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }
 bool  verbose = 0;
  
  // std::string jsonFileName="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt";
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
  jsonMap = readJSONFile(jsonFileName);
  std::cout<<"JSON file: "<<jsonFileName<<std::endl;

  // define map with events                                                                                                                                     
  std::map<std::pair<int,std::pair<int,int> >,int> eventsMap;

  //applyTrigger=false;
  std::cout<<"apply trigger: "<<applyTrigger<<std::endl;

  TLorentzVector W,W_puppi,LEP;
  TLorentzVector Wjet1_AK4, Wjet2_AK4, TOT_Wjet, Wjets1, Wjets2, TOT_Wjet_Final;
  TLorentzVector NU0,NU1,NU2,NU0_puppi,NU1_puppi,NU2_puppi;
  TLorentzVector NU0_jes_up, NU0_jes_dn;
  TLorentzVector JET, JET_PuppiAK8, AK4;
  TLorentzVector JET_jes_up, JET_jes_dn, JET_PuppiAK8_jes_up, JET_PuppiAK8_jes_dn;
  TLorentzVector AK4_JET1,AK4_JET2;
  TLorentzVector AK4_JET1_jes_up, AK4_JET1_jes_dn;
  TLorentzVector AK4_JET2_jes_up, AK4_JET2_jes_dn;
  TLorentzVector PuppiAK4_JET1,PuppiAK4_JET2;
  TLorentzVector PuppiAK4_JET1_jes_up, PuppiAK4_JET1_jes_dn;
  TLorentzVector PuppiAK4_JET2_jes_up, PuppiAK4_JET2_jes_dn;
  TLorentzVector VBF1,VBF2,TOT;
  TLorentzVector ELE,MU;
  
  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> looseMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> looseEle;

  int ok=0, total=0;
  int runno=260627;
  int lumo=524;
  int evento=942309369;
  int count=0;
	
  setInputTree *ReducedTree = new setInputTree (inputTreeName.c_str());
  ReducedTree->Init();


  char command1[3000];
  //sprintf(command1, "xrd eoscms dirlist %s/%s/  | awk '{print \"root://eoscms.cern.ch/\"$5}' > listTemp_%s.txt", (inputFolder).c_str(), (inputFile).c_str(), outputFile.c_str());
  sprintf(command1, "ls /tmp/rasharma/ueos/user/r/rasharma/aQGC_Studies/FirstStepOutput/%s/*.root  | awk '{print $1}' > listTemp_%s.txt", (inputFile).c_str(), outputFile.c_str());
  std::cout<<command1<<std::endl;
  std::cout<<"Scale Factor for this sample = "<<weight<<endl;
  system(command1);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);

  int fileCounter=0;
  Long64_t totalEntries=0;
  
  if (isLocal==1) {
    ReducedTree->fChain->Add("ReducedSelection.root");
  }
  else {
    while (!rootList.eof())
    {
      char iRun_tW[700];
      rootList >> iRun_tW;
      ReducedTree->fChain->Add(iRun_tW);
      fileCounter++;
    }
  }
  
  std::cout<<"number of files found: "<<fileCounter-2<<std::endl;
  totalEntries=ReducedTree->fChain->GetEntries();
  std::cout<<"total entries: "<<totalEntries<<std::endl;
  
  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  
  //--------pile up file -----------------
  std::vector<double> weights_pu1; //these are made with our recipe
  TFile* pileupFile1 = TFile::Open("pileupDataRun2016B_71p3mb.root");
  TH1F* pileupHisto1 = (TH1F*)pileupFile1->Get("pileup");  
  weights_pu1 = generate_weights(pileupHisto1,0);
  pileupFile1->Close();
  
  std::vector<double> weights_pu2; //these are made with the official recipe
  TFile* pileupFile2 = TFile::Open("PUxSynch.root");  
  TH1F* pileupHisto2 = (TH1F*)pileupFile2->Get("puweights");
  weights_pu2 = generate_weights(pileupHisto2,1);
  pileupFile2->Close();


  //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: Starts -------------
  TFile* IDIsoEleA = TFile::Open("egammaEffi.txt_SF2D_ID_ISO_RunB_0pt5fb.root","READ");	// First Root file for ID&ISO for ele
  TH1F *hIDIsoEleA = (TH1F*)IDIsoEleA->Get("EGamma_SF2D");
  double IDIsoEleA_Lumi = 0.5;	// Lumi corresponding to above file

  TFile* IDIsoEleB = TFile::Open("egammaEffi.txt_SF2D_ID_ISO_RunB_5pt0fb.root","READ"); // Second Root file for ID&ISO for ele
  TH1F *hIDIsoEleB = (TH1F*)IDIsoEleB->Get("EGamma_SF2D");
  double IDIsoEleB_Lumi = 5.3;	// Lumi corresponding to above file

	// Variables to define: 
  double TempPt, TempEta;
  int binPt, binEta;
  double EleIDIsoScaleA, EleIDIsoScaleB;

  TFile* GSFCorrEle = TFile::Open("egammaEffi.txt_SF2D_GSF_tracking.root","READ");
  TH1F *hGSFCorrEle = (TH1F*)GSFCorrEle->Get("EGamma_SF2D");

  TFile* TriggerEffEle = TFile::Open("SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root","READ");
  TH1F *hTriggerEffEle = (TH1F*)TriggerEffEle->Get("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt_DATA");

  TFile* IDMu = TFile::Open("MuonID_Z_RunBCD_prompt80X_7p65.root","READ");
  TH1F *hIDMu = (TH1F*)IDMu->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio");

  TFile* IsoMu = TFile::Open("MuonIso_Z_RunBCD_prompt80X_7p65.root","READ");
  TH1F *hIsoMu = (TH1F*)IsoMu->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");

  TFile* TriggerEffMu = TFile::Open("SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root","READ");
  TH1F *hTriggerEffMu1 = (TH1F*)TriggerEffMu->Get("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt_DATA");
  TH1F *hTriggerEffMu2 = (TH1F*)TriggerEffMu->Get("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097/efficienciesDATA/abseta_pt_DATA");
  double MuonSF1, MuonSF2;
  double hTriggerEffMu1_Lumi = 0.622;	// Need to put actual val
  double hTriggerEffMu2_Lumi = 4.378;	// Need to put actual val
//IDIsoEle->Close();	GSFCorrEle->Close();	 IDMu->Close();		IsoMu->Close();		TriggerEffMu->Close();

  //---------------- Root Files for ID, ISO, Trigger, GSF correctiosn for ELE and MU both: ENDS -------------

  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  setOutputTree* WWTree = new setOutputTree(outTree);

  float top_NNLO_weight[2];

  std::ifstream badEventsFile;
  std::multimap<int,int> badEventsList;
  
  int run, lumi, evt;
  if (isMC==0) {
    if (strcmp(leptonName.c_str(),"el")==0)
      badEventsFile.open("SingleElectron_csc2015.txt");
    else
      badEventsFile.open("SingleMuon_csc2015.txt");      
    while(!badEventsFile.eof()) 
      {
        badEventsFile >> run >> lumi >> evt;
        badEventsList.insert(std::pair<int,int>(run,evt));
      }      
  }
  badEventsFile.close();

  int nNegEvents=0; 
  int nEvents=0;	nEvents=totalEntries;
  Long64_t jentry2=0;
  
  //---------start loop on events------------
  std::cout << "---------start loop on events------------" << std::endl;
  //for (Long64_t jentry=0; jentry<10000;jentry++,jentry2++)
  for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++,jentry2++)
  {
    //for (Long64_t jentry=531000; jentry<532000;jentry++,jentry2++) {
    
    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;
    ReducedTree->fChain->GetEntry(jentry);
    // if (Cut(ientry) < 0) continue;                                                                                                                           
    
    tightMuon.clear();
    tightEle.clear();
    looseMuon.clear();
    looseEle.clear();

    if (jentry2%1000 == 0) std::cout << "read entry: " << jentry2 <<"/"<<totalEntries<<std:: endl;
    
    //*********************************
    // JSON FILE AND DUPLIACTES IN DATA
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi  = ReducedTree->LumiBlockNum;
    
    /*    
    bool skipEvent = false;
    if( isMC==0 ) //apply json file
    {
      if(AcceptEventByRunAndLumiSection(ReducedTree->RunNum,ReducedTree->LumiBlockNum,jsonMap) == false) skipEvent = true;                           
      std::pair<int,Long64_t> eventLSandID(ReducedTree->LumiBlockNum,ReducedTree->EvtNum);            
      std::pair<int,std::pair<int,Long64_t> > eventRUNandLSandID(ReducedTree->RunNum,eventLSandID); 
      if( eventsMap[eventRUNandLSandID] == 1 ) skipEvent = true;                                             
      else eventsMap[eventRUNandLSandID] = 1;                                                          
    }
    if( skipEvent == true ) continue;    
    */    
    
    WWTree->initializeVariables(); //initialize all variables
    
    if (ReducedTree->passFilterHBHELooseRerun == 0) continue;
    if (ReducedTree->passFilterHBHEIsoRerun == 0) continue;
    if (ReducedTree->passFilterCSCHalo == 0) continue;
    if (ReducedTree->passFilterGoodVtx == 0) continue;
    if (ReducedTree->passFilterEEBadSC == 0) continue;
    
    WWTree->issignal = 0;
    WWTree->issignal_AK4jetjet=0;
    WWTree->isVBFJet = 0;
    WWTree->isVBFJet_NoB = 0;
    WWTree->isWJets = 0;
    WWTree->wSampleWeight = weight; //xsec/numberOfEntries
    WWTree->eff_and_pu_Weight = 1.; //temporary value
    WWTree->eff_and_pu_Weight_2 = 1.; //temporary value
    WWTree->eff_and_pu_Weight_3 = 1.; //temporary value
    WWTree->top1_NNLO_Weight = 1.;
    WWTree->top2_NNLO_Weight = 1.;
    WWTree->trig_eff_Weight = 1.;

    if (ReducedTree->genEventWeight>0)
      WWTree->genWeight=1.;
    else if (ReducedTree->genEventWeight<0) {
      WWTree->genWeight=-1.;
      nNegEvents++;
    }
    // std::cout<<ReducedTree->genEventWeight<<" "<<WWTree->genWeight<<std::endl;
    // WWTree->genWeight = ReducedTree->genEventWeight;
    
    //PILE-UP WEIGHT
    if (isMC==1) {
      if(ReducedTree->NVtx<int(weights_pu1.size())){
	WWTree->eff_and_pu_Weight = weights_pu1[ReducedTree->npT]; //official pu recipe
      }
      else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	std::cout<<"Warning! n_pu too big"<<std::endl;
	// throw logic_error("n_pu too big");
        WWTree->eff_and_pu_Weight = 0.;
      }
      
      if(ReducedTree->NVtx<int(weights_pu2.size())){
	WWTree->eff_and_pu_Weight_2 = weights_pu2[ReducedTree->NVtx]; //our pu recipe
      }
      else{ //should not happen as we have a weight for all simulated n_pu multiplicities!
	std::cout<<"Warning! n_pu too big"<<std::endl;
	// throw logic_error("n_pu too big");
        WWTree->eff_and_pu_Weight_2 = 0.;
      }    
      WWTree->eff_and_pu_Weight_3 = ReducedTree->PUWeight; //our pu recipe
    }
    
    //save event variables
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi  = ReducedTree->LumiBlockNum;
    
    WWTree->nPV = ReducedTree->NVtx;
    
    count=0;
    cutEff[0]++;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;
    
    
    /////////////////THE SELECTED LEPTON
    int nTightLepton=0;
    if (strcmp(leptonName.c_str(),"el")==0) { //electrons
      int passTrigger=0;
      float tempPt=0.;
      for (int i=0; i<ReducedTree->ElectronsNum; i++) {
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
	if (applyTrigger==1)
	  for (unsigned int t=0; t<ReducedTree->TriggerProducerTriggerNames->size(); t++)
	    if(//TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Ele27_eta2p1_WP75_Gsf_v") || 
	       TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Ele27_eta2p1_WPLoose_Gsf_v") 
	       //|| TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_Ele105_CaloIdVT_GsfTrkIdT_v") 
	       )
	      if (ReducedTree->TriggerProducerTriggerPass->at(t)==1) passTrigger=1; //trigger
	if (passTrigger==0 && applyTrigger==1) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
	if (ReducedTree->Electrons_isTight[i]==false) continue;       
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
        if (ReducedTree->ElectronsPt[i]<=30) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
        if (fabs(ReducedTree->ElectronsEta[i])>=2.1) continue;
	if (ReducedTree->ElectronsPt[i]<tempPt) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug ele: "<<i<<std::endl;
	ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
	tightEle.push_back(ELE);
	WWTree->l_pt  = ReducedTree->ElectronsPt[i];
	WWTree->l_eta = ReducedTree->ElectronsEta[i];
	WWTree->l_phi = ReducedTree->ElectronsPhi[i];	
	WWTree->l_e= ReducedTree->ElectronsE[i];	
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    } //end electrons
    else if (strcmp(leptonName.c_str(),"mu")==0){ //muons
      int passTrigger=0;
      float tempPt=0.;
      for (int i=0; i<ReducedTree->MuonsNum; i++) {
	if (applyTrigger==1)
	  for (unsigned int t=0; t<ReducedTree->TriggerProducerTriggerNames->size(); t++)
	    if(TString(ReducedTree->TriggerProducerTriggerNames->at(t)).Contains("HLT_IsoMu22_v"))
	      if (ReducedTree->TriggerProducerTriggerPass->at(t)==1) passTrigger=1; //trigger
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	if (passTrigger==0 && applyTrigger==1) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	if (ReducedTree->Muons_isTight[i]==false) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
        if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
        if (ReducedTree->MuonsPt[i]<25) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
        if (fabs(ReducedTree->MuonsEta[i])>=2.4) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
	tightMuon.push_back(MU);
	if (ReducedTree->MuonsPt[i]<tempPt) continue;
	if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug mu: "<<i<<std::endl;
	WWTree->l_pt  = ReducedTree->MuonsPt[i];
	WWTree->l_eta = ReducedTree->MuonsEta[i];
	WWTree->l_phi = ReducedTree->MuonsPhi[i];
	WWTree->l_e = ReducedTree->MuonsE[i];
	tempPt = WWTree->l_pt;
	nTightLepton++;
      }
    } //end muons
    if (nTightLepton==0) continue; //no leptons with required ID
    cutEff[1]++;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;
    
    LEP.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    
    //trigger SF for muon (from B2G-15-005, these are for the HLT_IsoMu20 trigger,
    //see https://indico.cern.ch/event/462268/contribution/9/attachments/1188638/1724574/2015.11.17_MuonPOG_SingleMuTrigEff_SF_KPLee_v2.pdf
    if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) {
    	//----------- Start Applying Trigger SF's for muons	--------------------
	TempPt = WWTree->l_pt;
	TempEta= fabs(WWTree->l_eta);
	if (TempPt > hTriggerEffMu1->GetYaxis()->GetXmax())		// Check if pt is not ouside upper limit; if so then to assign same SF as the last bin to high pt resacle tempPt to less then maximum range of pt defined in histo.
		TempPt = hTriggerEffMu1->GetYaxis()->GetXmax() - 1.0;	
	if (TempPt < hTriggerEffMu1->GetYaxis()->GetXmin())		// Check if pt is not ouside lower limit; if so then to assign same SF as the first bin to low pt resacle tempPt to greater then minimum range of pt defined in histo.
		TempPt = hTriggerEffMu1->GetYaxis()->GetXmin() + 1.0;
	binEta = hTriggerEffMu1->GetXaxis()->FindFixBin(TempEta);	// Get Eta bin 
	binPt  = hTriggerEffMu1->GetYaxis()->FindFixBin(TempPt);	// Get Pt Bin

	MuonSF1 = hTriggerEffMu1->GetBinContent( binEta, binPt );	// Get Scale factor corresponding to the pt and eta.
	if (verbose){
	cout<<"A: WWTree->l_pt = "<<WWTree->l_pt<<"\tTempPt = "<<TempPt<<"\tTempEta = "<<TempEta<<"\tMuonSF1 = "<<MuonSF1<<"\tMuonSF2 = "<<MuonSF2<<"\tWWTree->trig_eff_Weight = "<<WWTree->trig_eff_Weight<<endl;
}
	TempPt = WWTree->l_pt;
	TempEta= fabs(WWTree->l_eta);
	if (TempPt > hTriggerEffMu2->GetYaxis()->GetXmax())		// Check if pt is not ouside upper limit; if so then to assign same SF as the last bin to high pt resacle tempPt to less then maximum range of pt defined in histo.
		TempPt = hTriggerEffMu2->GetYaxis()->GetXmax() - 1.0;	
	if (TempPt < hTriggerEffMu2->GetYaxis()->GetXmin())		// Check if pt is not ouside lower limit; if so then to assign same SF as the first bin to low pt resacle tempPt to greater then minimum range of pt defined in histo.
		TempPt = hTriggerEffMu2->GetYaxis()->GetXmin() + 1.0;
	binEta = hTriggerEffMu2->GetXaxis()->FindFixBin(TempEta);	// Get Eta bin 
	binPt  = hTriggerEffMu2->GetYaxis()->FindFixBin(TempPt);	// Get Pt Bin

	MuonSF2 = hTriggerEffMu2->GetBinContent( binEta, binPt );	// Get Scale factor corresponding to the pt and eta.

	WWTree->trig_eff_Weight = (MuonSF1 * hTriggerEffMu1_Lumi + MuonSF2 * hTriggerEffMu2_Lumi )/ (hTriggerEffMu1_Lumi + hTriggerEffMu2_Lumi);

	if (verbose)
	cout<<"B: TempPt = "<<TempPt<<"\tTempEta = "<<TempEta<<"\tMuonSF1 = "<<MuonSF1<<"\tMuonSF2 = "<<MuonSF2<<"\tWWTree->trig_eff_Weight = "<<WWTree->trig_eff_Weight<<endl;
    }
    
    //trigger SF for electron (from B2G-15-005)
    if (strcmp(leptonName.c_str(),"el")==0 && isMC==1) {
    	//----------- Start Applying Trigger SF's for Electrons	--------------------
	TempPt = WWTree->l_pt;
	TempEta= WWTree->l_eta;
	if (TempPt > hTriggerEffEle->GetYaxis()->GetXmax())		// Check if pt is not ouside upper limit; if so then to assign same SF as the last bin to high pt resacle tempPt to less then maximum range of pt defined in histo.
		TempPt = hTriggerEffEle->GetYaxis()->GetXmax() - 1.0;	
	if (TempPt < hTriggerEffEle->GetYaxis()->GetXmin())		// Check if pt is not ouside lower limit; if so then to assign same SF as the first bin to low pt resacle tempPt to greater then minimum range of pt defined in histo.
		TempPt = hTriggerEffEle->GetYaxis()->GetXmin() + 1.0;
	binEta = hTriggerEffEle->GetXaxis()->FindFixBin(TempEta);	// Get Eta bin 
	binPt  = hTriggerEffEle->GetYaxis()->FindFixBin(TempPt);	// Get Pt Bin

	WWTree->trig_eff_Weight = hTriggerEffEle->GetBinContent( binEta, binPt );	// Get Scale factor corresponding to the pt and eta.
	if (verbose)
	cout<<"TempPt = "<<TempPt<<"\tTempEta = "<<TempEta<<"\tWWTree->trig_eff_Weight = "<<WWTree->trig_eff_Weight<<endl;
    }
    
    //ID & GSF efficiency SF for electrons (https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale)
    if (strcmp(leptonName.c_str(),"el")==0 && isMC==1) {
    	//----------- Start Applying ID, ISO SF's for electrons	--------------------
	TempPt = WWTree->l_pt;
	TempEta= WWTree->l_eta;
	if (TempPt > hIDIsoEleA->GetYaxis()->GetXmax())		// Check if pt is not ouside upper limit; if so then to assign same SF as the last bin to high pt resacle tempPt to less then maximum range of pt defined in histo.
		TempPt = hIDIsoEleA->GetYaxis()->GetXmax() - 1.0;	
	if (TempPt < hIDIsoEleA->GetYaxis()->GetXmin())		// Check if pt is not ouside lower limit; if so then to assign same SF as the first bin to low pt resacle tempPt to greater then minimum range of pt defined in histo.
		TempPt = hIDIsoEleA->GetYaxis()->GetXmin() + 1.0;
	binEta = hIDIsoEleA->GetXaxis()->FindFixBin(TempEta);	// Get Eta bin 
	binPt  = hIDIsoEleA->GetYaxis()->FindFixBin(TempPt);	// Get Pt Bin

	EleIDIsoScaleA = hIDIsoEleA->GetBinContent( binEta, binPt );	// Get Scale factor corresponding to the pt and eta.

	// Reset TempPt & TempEta to the pt of eta and pt of selected electron
	TempPt = WWTree->l_pt;
	TempEta= WWTree->l_eta;
	if (TempPt > hIDIsoEleB->GetYaxis()->GetXmax())
		TempPt = hIDIsoEleB->GetYaxis()->GetXmax() - 1.0;
	if (TempPt < hIDIsoEleB->GetYaxis()->GetXmin())
		TempPt = hIDIsoEleB->GetYaxis()->GetXmax() + 1.0;
	binEta = hIDIsoEleB->GetXaxis()->FindFixBin(TempEta);
	binPt  = hIDIsoEleB->GetYaxis()->FindFixBin(TempPt);

	EleIDIsoScaleB = hIDIsoEleB->GetBinContent( binEta, binPt );

	WWTree->id_eff_Weight = (EleIDIsoScaleA * IDIsoEleA_Lumi + EleIDIsoScaleB * IDIsoEleB_Lumi) / (IDIsoEleA_Lumi + IDIsoEleB_Lumi);
	//
    	//----------- END: Applying ID, ISO SF's for electrons	--------------------
		
    	//----------- Start Applying GSF/RECO SF's for electrons	--------------------
	// Reset TempPt & TempEta to the pt of eta and pt of selected electron
	TempPt = WWTree->l_pt;
	TempEta= WWTree->l_eta;
	if (TempPt > hGSFCorrEle->GetYaxis()->GetXmax())
		TempPt = hGSFCorrEle->GetYaxis()->GetXmax() - 1.0;
	if (TempPt < hGSFCorrEle->GetYaxis()->GetXmin())
		TempPt = hGSFCorrEle->GetYaxis()->GetXmax() + 1.0;
	binEta = hGSFCorrEle->GetXaxis()->FindFixBin(TempEta);
	binPt  = hGSFCorrEle->GetYaxis()->FindFixBin(TempPt);

	WWTree->id_eff_Weight = WWTree->id_eff_Weight*hGSFCorrEle->GetBinContent(binEta, binPt);
    	//----------- END: Applying GSF/RECO SF's for electrons	--------------------

    }
    
    //ID & ISO efficiency SF for muons (https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults)
    if (strcmp(leptonName.c_str(),"mu")==0 && isMC==1) {
    	//----------- Start Applying ID SF's for Muons	--------------------
	TempPt = WWTree->l_pt;
	TempEta= fabs(WWTree->l_eta);
	if (TempPt > hIDMu->GetYaxis()->GetXmax())		// Check if pt is not ouside upper limit; if so then to assign same SF as the last bin to high pt resacle tempPt to less then maximum range of pt defined in histo.
		TempPt = hIDMu->GetYaxis()->GetXmax() - 1.0;	
	if (TempPt < hIDMu->GetYaxis()->GetXmin())		// Check if pt is not ouside lower limit; if so then to assign same SF as the first bin to low pt resacle tempPt to greater then minimum range of pt defined in histo.
		TempPt = hIDMu->GetYaxis()->GetXmin() + 1.0;
	binEta = hIDMu->GetXaxis()->FindFixBin(TempEta);	// Get Eta bin 
	binPt  = hIDMu->GetYaxis()->FindFixBin(TempPt);	// Get Pt Bin

	WWTree->id_eff_Weight = hIDMu->GetBinContent(binEta, binPt);

    	//----------- Start Applying ISO SF's for Muons	--------------------
	TempPt = WWTree->l_pt;
	TempEta= fabs(WWTree->l_eta);
	if (TempPt > hIsoMu->GetYaxis()->GetXmax())		// Check if pt is not ouside upper limit; if so then to assign same SF as the last bin to high pt resacle tempPt to less then maximum range of pt defined in histo.
		TempPt = hIsoMu->GetYaxis()->GetXmax() - 1.0;	
	if (TempPt < hIsoMu->GetYaxis()->GetXmin())		// Check if pt is not ouside lower limit; if so then to assign same SF as the first bin to low pt resacle tempPt to greater then minimum range of pt defined in histo.
		TempPt = hIsoMu->GetYaxis()->GetXmin() + 1.0;
	binEta = hIsoMu->GetXaxis()->FindFixBin(TempEta);	// Get Eta bin 
	binPt  = hIsoMu->GetYaxis()->FindFixBin(TempPt);	// Get Pt Bin

	WWTree->id_eff_Weight = WWTree->id_eff_Weight * hIsoMu->GetBinContent(binEta, binPt);
    }
    
    //VETO ADDITIONAL LEPTONS
    int nLooseLepton=0;
    for (int i=0; i<ReducedTree->ElectronsNum; i++) {
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose el: "<<i<<std::endl; count++;
      if (ReducedTree->Electrons_isLoose[i]==false) continue; //NB: CHANGE TO VETO!!!
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose el: "<<i<<", pt: "<<ReducedTree->ElectronsPt[i]<<", eta: "<<ReducedTree->ElectronsEta[i]<<std::endl; count++;
      if (ReducedTree->ElectronsPt[i]<27) continue;       
      // if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<i<<std::endl; count++;
      // if (fabs(ReducedTree->ElectronsEta[i])>=2.5) continue;
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose el: "<<i<<std::endl; count++;
      ELE.SetPtEtaPhiE(ReducedTree->ElectronsPt[i],ReducedTree->ElectronsEta[i],ReducedTree->ElectronsPhi[i],ReducedTree->ElectronsE[i]);
      looseEle.push_back(ELE);      
      nLooseLepton++;
    }
    for (int i=0; i<ReducedTree->MuonsNum; i++) {
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      if (ReducedTree->Muons_isLoose[i]==false) continue;
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      if (fabs(ReducedTree->MuonsEta[i])>=2.4) continue;
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      if (ReducedTree->MuonsPt[i]<22) continue;
      if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug loose mu: "<<i<<std::endl; count++;
      MU.SetPtEtaPhiE(ReducedTree->MuonsPt[i],ReducedTree->MuonsEta[i],ReducedTree->MuonsPhi[i],ReducedTree->MuonsE[i]);
      looseMuon.push_back(MU);
      nLooseLepton++;
    }
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo)     std::cout<<"nlooseleptons: "<<nLooseLepton<<std::endl;
    if (nLooseLepton!=1) continue; //no additional leptons
    cutEff[2]++;
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;
    
    
    
    //////////////THE MET
    
    //preselection on met
    //if (ReducedTree->METPt < 30) continue;
    cutEff[3]++;
    // if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;
    
    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*

    float Wmass = 80.385;

    TLorentzVector W_Met, W_Met_jes_up, W_Met_jes_dn;
    
    W_Met.SetPxPyPzE(ReducedTree->METPt * TMath::Cos(ReducedTree->METPhi), ReducedTree->METPt * TMath::Sin(ReducedTree->METPhi), 0., sqrt(ReducedTree->METPt*ReducedTree->METPt));
    W_Met_jes_up.SetPxPyPzE(ReducedTree->METPtUp * TMath::Cos(ReducedTree->METPhiUp), ReducedTree->METPtUp * TMath::Sin(ReducedTree->METPhiUp), 0., sqrt(ReducedTree->METPtUp*ReducedTree->METPtUp));
    W_Met_jes_dn.SetPxPyPzE(ReducedTree->METPtDown * TMath::Cos(ReducedTree->METPhiDown), ReducedTree->METPtDown * TMath::Sin(ReducedTree->METPhiDown), 0., sqrt(ReducedTree->METPtDown*ReducedTree->METPtDown));
    
    if(LEP.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ;  }
    cutEff[4]++;
    // if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;
    
    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    METzCalculator NeutrinoPz_type0_jes_up;
    METzCalculator NeutrinoPz_type0_jes_dn;
    METzCalculator_Run2 NeutrinoPz_run2;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(LEP);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_type0_jes_up.SetMET(W_Met_jes_up);
    NeutrinoPz_type0_jes_up.SetLepton(LEP);
    NeutrinoPz_type0_jes_up.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_type0_jes_dn.SetMET(W_Met_jes_dn);
    NeutrinoPz_type0_jes_dn.SetLepton(LEP);
    NeutrinoPz_type0_jes_dn.SetLeptonType(leptonName.c_str());
    
    NeutrinoPz_run2.SetMET(W_Met);
    NeutrinoPz_run2.SetLepton(LEP);
    NeutrinoPz_run2.SetLeptonType(leptonName.c_str());
    
    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    //double pz2_type0 = NeutrinoPz_type0.getOther();  // Default one
    
    double pz1_run2 = NeutrinoPz_run2.Calculate();
    
    double pz1_type0_jes_up = NeutrinoPz_type0_jes_up.Calculate(); // Default one -> according to type0
    double pz1_type0_jes_dn = NeutrinoPz_type0_jes_dn.Calculate(); // Default one -> according to type0
    
    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));
    
    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi), nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi), nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }
    
    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(LEP);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());
    
    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    //double pz2_type2 = NeutrinoPz_type2.getOther();   // Default one
    
    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    // change the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2; 
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    
    if (NeutrinoPz_type2.IsComplex()) {// if this is a complex, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi), nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi), nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );
      
      if ( fabs((LEP+W_neutrino_1).M()-Wmass) < fabs((LEP+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }
    
    WWTree->pfMET = sqrt(ReducedTree->METPt*ReducedTree->METPt);
    WWTree->pfMET_jes_up = sqrt(ReducedTree->METPtUp*ReducedTree->METPtUp);
    WWTree->pfMET_jes_dn = sqrt(ReducedTree->METPtDown*ReducedTree->METPtDown);
    WWTree->pfMET_Phi = ReducedTree->METPhi;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;
    WWTree->nu_pz_run2 = pz1_run2;
    WWTree->nu_pz_isre = 1-NeutrinoPz_run2.IsComplex();
    WWTree->nu_pz_run2_oth = NeutrinoPz_run2.getOther();
    WWTree->nu_pz_run2_type = NeutrinoPz_run2.getType();
    
    
    /////////////////THE LEPTONIC W
    
    NU0.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type0,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type0*WWTree->nu_pz_type0));
    NU0_jes_up.SetPxPyPzE(ReducedTree->METPtUp*TMath::Cos(ReducedTree->METPhiUp),ReducedTree->METPtUp*TMath::Sin(ReducedTree->METPhiUp),pz1_type0_jes_up,TMath::Sqrt(WWTree->pfMET_jes_up*WWTree->pfMET_jes_up+pz1_type0_jes_up*pz1_type0_jes_up));
    NU0_jes_dn.SetPxPyPzE(ReducedTree->METPtDown*TMath::Cos(ReducedTree->METPhiDown),ReducedTree->METPtDown*TMath::Sin(ReducedTree->METPhiDown),pz1_type0_jes_dn,TMath::Sqrt(WWTree->pfMET_jes_dn*WWTree->pfMET_jes_dn+pz1_type0_jes_dn*pz1_type0_jes_dn));
    
    NU2.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_type2*WWTree->nu_pz_type2));
    NU1.SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_run2,TMath::Sqrt(WWTree->pfMET*WWTree->pfMET+WWTree->nu_pz_run2*WWTree->nu_pz_run2));
    
    W = LEP + NU2;
    
    WWTree->v_pt = W.Pt();
    WWTree->v_eta = W.Eta();
    WWTree->v_phi = W.Phi();
    WWTree->v_mt = TMath::Sqrt(2*LEP.Et()*NU2.Et()*(1-TMath::Cos(LEP.DeltaPhi(NU2))));


    
    ///////////THE JET - AK4
    //WWTree->njets=0;
    //WWTree->nVBFPairs=0;
    //int hadWpos = -1;
    //int ttb_jet_position=-1; //position of AK8 jet in ttbar-topology
    // if (ReducedTree->AK8JetsNum < 1 ) continue; 
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"debug: "<<count<<std::endl; count++;
    
    ////	Select all possible AK4 jet candidate
    //
    //float tempPt1 = 0.; int pos1 = -1;
    //float tempPt2 = 0.; int pos2 = -1;
    int nGoodAK4jets=0;
    WWTree->njets=0;
    std::vector<int> indexGoodAK4Jets;

    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
    {
      bool isCleanedJet = true;
      if (ReducedTree->Jets_PtCorr[i]<=30 || ReducedTree->JetsPt[i]<=20 || fabs(ReducedTree->JetsEta[i])>=2.4)  continue;
      if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
      if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.800) continue;
      
      //CLEANING FROM LEPTONS
      for (unsigned int j=0; j<tightEle.size(); j++) {
        if (deltaR(tightEle.at(j).Eta(), tightEle.at(j).Phi(),
                   ReducedTree->JetsEta[i], ReducedTree->JetsPhi[i]) < 0.3) {
          isCleanedJet = false;
        }
      }
      for (unsigned int j=0; j<tightMuon.size(); j++) {
        if (deltaR(tightMuon.at(j).Eta(), tightMuon.at(j).Phi(),
                   ReducedTree->JetsEta[i], ReducedTree->JetsPhi[i]) < 0.3) {
          isCleanedJet = false;
        }
      }
      
      if (isCleanedJet==false) continue;
      WWTree->njets++;
      nGoodAK4jets++;
      indexGoodAK4Jets.push_back(i); //save index of the "good" jets candidate
      }
      if (indexGoodAK4Jets.size()<4)  continue;
      cutEff[5]++;

    float DeltaEta = 0.;
    int nVBF1=-1, nVBF2=-1; //position of the two vbf jets

    int nGoodAK4VBFjets = 0;
    WWTree->nVBFPairs=0;

    //================================ Selection of VBF jets: BEGIN =====================================
    for(int i=0; i<indexGoodAK4Jets.size()-1;i++)
    {
        for(int j=i+1; j<indexGoodAK4Jets.size();j++)
        {
            VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodAK4Jets.at(i)],ReducedTree->JetsEta[indexGoodAK4Jets.at(i)],ReducedTree->JetsPhi[indexGoodAK4Jets.at(i)],ReducedTree->Jets_ECorr[indexGoodAK4Jets.at(i)]);
            VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodAK4Jets.at(j)],ReducedTree->JetsEta[indexGoodAK4Jets.at(j)],ReducedTree->JetsPhi[indexGoodAK4Jets.at(j)],ReducedTree->Jets_ECorr[indexGoodAK4Jets.at(j)]);

            if (DeltaEta > abs(VBF1.Eta()-VBF2.Eta()) || VBF1.Eta()*VBF2.Eta() > 0 || (VBF1+VBF2).M()<500) continue;
            if (abs(VBF1.Eta()-VBF2.Eta())<3.5) continue;

            DeltaEta = abs(VBF1.Eta()-VBF2.Eta()); //take the jet pair with largest DeltaEta
            nVBF1 = indexGoodAK4Jets.at(i); //save position of the 1st vbf jet
            nVBF2 = indexGoodAK4Jets.at(j); //save position of the 2nd vbf jet

    	nGoodAK4VBFjets++;
	WWTree->nVBFPairs++;
        }
    }

    if (nGoodAK4VBFjets > 0) 
	{
	WWTree->isVBFJet=1;
	cutEff[6]++;
	}

    if (nGoodAK4VBFjets > 0 && ReducedTree->Jets_bDiscriminatorICSV[nVBF1] < 0.935 && ReducedTree->Jets_bDiscriminatorICSV[nVBF2]< 0.935){
	WWTree->isVBFJet_NoB=1;
    	cutEff[7]++;
	}
      if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
      {
        // nVBF1=0; nVBF2=1;
        
        VBF1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF1],ReducedTree->JetsEta[nVBF1],ReducedTree->JetsPhi[nVBF1],ReducedTree->Jets_ECorr[nVBF1]);
        VBF2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nVBF2],ReducedTree->JetsEta[nVBF2],ReducedTree->JetsPhi[nVBF2],ReducedTree->Jets_ECorr[nVBF2]);
        TOT = VBF1 + VBF2;
	
        WWTree->vbf_maxpt_j1_pt = ReducedTree->Jets_PtCorr[nVBF1];
        WWTree->vbf_maxpt_j1_eta = ReducedTree->JetsEta[nVBF1];
        WWTree->vbf_maxpt_j1_phi = ReducedTree->JetsPhi[nVBF1];
        WWTree->vbf_maxpt_j1_e = ReducedTree->Jets_ECorr[nVBF1];
        WWTree->vbf_maxpt_j1_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nVBF1];
        WWTree->vbf_maxpt_j2_pt = ReducedTree->Jets_PtCorr[nVBF2];
        WWTree->vbf_maxpt_j2_eta = ReducedTree->JetsEta[nVBF2];
        WWTree->vbf_maxpt_j2_phi = ReducedTree->JetsPhi[nVBF2];
        WWTree->vbf_maxpt_j2_e = ReducedTree->Jets_ECorr[nVBF2];
        WWTree->vbf_maxpt_j2_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminatorICSV[nVBF2];
        WWTree->vbf_maxpt_jj_pt = TOT.Pt();
        WWTree->vbf_maxpt_jj_eta = TOT.Eta();
        WWTree->vbf_maxpt_jj_phi = TOT.Phi();
        WWTree->vbf_maxpt_jj_m = TOT.M();	
      }

    //================================ Selection of W-jets: STARTS  =====================================
    int nGoodAK4Wjets = 0;
    int nWjets1 = -1, nWjets2 = -1 ;
    double DeltaMassWindow = 80.;
    #if 0
    for(int i=0; i<indexGoodAK4Jets.size()-1;i++)
    {
        for(int j=i+1; j<indexGoodAK4Jets.size();j++)
        {
    if(indexGoodAK4Jets.at(i) == nVBF1 || indexGoodAK4Jets.at(j) == nVBF2 ) continue;
    if(indexGoodAK4Jets.at(j) == nVBF1 || indexGoodAK4Jets.at(i) == nVBF2 ) continue;
//      if (fabs(ReducedTree->JetsEta[indexGoodAK4Jets.at(i)])>=2.0)  continue;
//      if (fabs(ReducedTree->JetsEta[indexGoodAK4Jets.at(j)])>=3.0)  continue;
    //coutWjets++;
            Wjet1_AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodAK4Jets.at(i)],ReducedTree->JetsEta[indexGoodAK4Jets.at(i)],ReducedTree->JetsPhi[indexGoodAK4Jets.at(i)],ReducedTree->Jets_ECorr[indexGoodAK4Jets.at(i)]);
            Wjet2_AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodAK4Jets.at(j)],ReducedTree->JetsEta[indexGoodAK4Jets.at(j)],ReducedTree->JetsPhi[indexGoodAK4Jets.at(j)],ReducedTree->Jets_ECorr[indexGoodAK4Jets.at(j)]);
            TOT_Wjet = Wjet1_AK4 + Wjet2_AK4 ;
            //cout<<"Found Before Check!!!!"<<endl;


            if (DeltaMassWindow < abs(TOT_Wjet.M() - 35.)) continue;
            //cout<<"(Wjet1_AK4+Wjet2_AK4).M()-80. = "<<(Wjet1_AK4+Wjet2_AK4).M()<<endl;

//          cout<< "====> " << Wjet1_AK4.Mag() << "\t" << Wjet2_AK4.Mag() << "\t" << TOT_Wjet.Mag() <<endl;

            DeltaMassWindow = abs(TOT_Wjet.M()-35.); //take the jet pair with largest DeltaEta
            nWjets1 = indexGoodAK4Jets.at(i); //save position of the 1st W-jet
            nWjets2 = indexGoodAK4Jets.at(j); //save position of the 2nd W-jet
            nGoodAK4Wjets++;
        }
    }
    #endif
    for(int i=0; i<indexGoodAK4Jets.size();i++)
    {
	if(indexGoodAK4Jets.at(0) == nVBF1 || indexGoodAK4Jets.at(0) == nVBF2 ) continue;
	nGoodAK4Wjets++;
	if (nGoodAK4Wjets == 1)
	{
            Wjet1_AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodAK4Jets.at(i)],ReducedTree->JetsEta[indexGoodAK4Jets.at(i)],ReducedTree->JetsPhi[indexGoodAK4Jets.at(i)],ReducedTree->Jets_ECorr[indexGoodAK4Jets.at(i)]);
            nWjets1 = indexGoodAK4Jets.at(i); //save position of the 1st Wjets
            nGoodAK4Wjets++;
	}
	if (nGoodAK4Wjets == 2)
	{
            Wjet2_AK4.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[indexGoodAK4Jets.at(i)],ReducedTree->JetsEta[indexGoodAK4Jets.at(i)],ReducedTree->JetsPhi[indexGoodAK4Jets.at(i)],ReducedTree->Jets_ECorr[indexGoodAK4Jets.at(i)]);
            nWjets2 = indexGoodAK4Jets.at(i); //save position of the 2nd W-jet
            nGoodAK4Wjets++;
	}
    }

    if (nGoodAK4Wjets == 0){
    WWTree->isWJets = 1;
    cutEff[8]++;}
    if ( nWjets1 != -1 && nWjets2 != -1 ) {
    //cout<<nWjets1<<"\t"<<nWjets2<<endl;
    Wjets1.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nWjets1],ReducedTree->JetsEta[nWjets1],ReducedTree->JetsPhi[nWjets1],ReducedTree->Jets_ECorr[nWjets1]);
    Wjets2.SetPtEtaPhiE(ReducedTree->Jets_PtCorr[nWjets2],ReducedTree->JetsEta[nWjets2],ReducedTree->JetsPhi[nWjets2],ReducedTree->Jets_ECorr[nWjets2]);
    TOT_Wjet_Final = Wjets1 + Wjets2 ;

    WWTree->AK4_jet1_pt  = ReducedTree->Jets_PtCorr[nWjets1];
    WWTree->AK4_jet1_eta = ReducedTree->JetsEta[nWjets1];
    WWTree->AK4_jet1_phi = ReducedTree->JetsPhi[nWjets1];
    WWTree->AK4_jet1_e   = ReducedTree->Jets_ECorr[nWjets1];
    WWTree->AK4_jet1_pt_jes_up = (ReducedTree->Jets_PtCorr[nWjets1]/ReducedTree->Jets_AK4correction[nWjets1])*ReducedTree->Jets_AK4correctionUp[nWjets1];
    WWTree->AK4_jet1_pt_jes_dn = (ReducedTree->Jets_PtCorr[nWjets1]/ReducedTree->Jets_AK4correction[nWjets1])*ReducedTree->Jets_AK4correctionDown[nWjets1];

    WWTree->AK4_jet2_pt  = ReducedTree->Jets_PtCorr[nWjets2];
    WWTree->AK4_jet2_eta = ReducedTree->JetsEta[nWjets2];
    WWTree->AK4_jet2_phi = ReducedTree->JetsPhi[nWjets2];
    WWTree->AK4_jet2_e   = ReducedTree->Jets_ECorr[nWjets2];
    WWTree->AK4_jet2_pt_jes_up = (ReducedTree->Jets_PtCorr[nWjets2]/ReducedTree->Jets_AK4correction[nWjets2])*ReducedTree->Jets_AK4correctionUp[nWjets2];
    WWTree->AK4_jet2_pt_jes_dn = (ReducedTree->Jets_PtCorr[nWjets2]/ReducedTree->Jets_AK4correction[nWjets2])*ReducedTree->Jets_AK4correctionDown[nWjets2];


    WWTree->AK4_jetjet_pt = (TOT_Wjet_Final).Pt();
    WWTree->AK4_jetjet_mass = (TOT_Wjet_Final).M();
    WWTree->AK4_jetjet_deltaeta = deltaEta(Wjets1.Eta(),Wjets2.Eta());
    WWTree->AK4_jetjet_deltaphi = deltaPhi(Wjets1.Phi(),Wjets2.Phi());
    WWTree->AK4_jetjet_deltar = deltaR(Wjets1.Eta(),Wjets1.Phi(),Wjets2.Eta(),Wjets2.Phi());



    if (WWTree->AK4_jet2_pt>0) {
      WWTree->deltaR_lak4jetjet = deltaR((Wjets1+Wjets2).Eta(),(Wjets1+Wjets2).Phi(),LEP.Eta(),LEP.Phi());
      WWTree->deltaphi_METak4jetjet = deltaPhi((Wjets1+Wjets2).Phi(),NU2.Phi());
      WWTree->deltaphi_Vak4jetjet = deltaPhi((Wjets1+Wjets2).Phi(),W.Phi());
      if (WWTree->deltaR_lak4jetjet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak4jetjet)>2.0 && fabs(WWTree->deltaphi_Vak4jetjet)>2.0)
        WWTree->issignal_AK4jetjet=1;
    }


    WWTree->mass_lvjj_type0_AK4 = (LEP + NU0 + Wjets1 + Wjets2).M();
    WWTree->mass_lvjj_type2_AK4 = (LEP + NU2 + Wjets1 + Wjets2).M();
    WWTree->mass_lvjj_run2_AK4  = (LEP + NU1 + Wjets1 + Wjets2).M();
    
   } 

/*     
    //////////////////FOUR-BODY INVARIANT MASS
    WWTree->mass_lvj_type0 = (LEP + NU0 + JET).M();
    WWTree->mass_lvj_type2 = (LEP + NU2 + JET).M();
    WWTree->mass_lvj_run2  = (LEP + NU1 + JET).M();
    WWTree->mass_lvj_type0_met_jes_up = (LEP + NU0_jes_up + JET_jes_up).M();
    WWTree->mass_lvj_type0_met_jes_dn = (LEP + NU0_jes_dn + JET_jes_dn).M();
    
    WWTree->mass_lvj_type0_PuppiAK8 = (LEP + NU0_puppi + JET_PuppiAK8).M();
    WWTree->mass_lvj_type2_PuppiAK8 = (LEP + NU2_puppi + JET_PuppiAK8).M();
    WWTree->mass_lvj_run2_PuppiAK8  = (LEP + NU1_puppi + JET_PuppiAK8).M();
    WWTree->mass_lvj_type0_met_PuppiAK8_jes_up = (LEP + NU0_jes_up + JET_PuppiAK8_jes_up).M();
    WWTree->mass_lvj_type0_met_PuppiAK8_jes_dn = (LEP + NU0_jes_dn + JET_PuppiAK8_jes_dn).M();
    
    WWTree->mass_lvjj_type0_AK4 = (LEP + NU0 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_type2_AK4 = (LEP + NU2 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_run2_AK4  = (LEP + NU1 + AK4_JET1 + AK4_JET2).M();
    WWTree->mass_lvjj_type0_met_jes_up_AK4 = (LEP + NU0_jes_up + AK4_JET1_jes_up + AK4_JET2_jes_up).M();
    WWTree->mass_lvjj_type0_met_jes_dn_AK4 = (LEP + NU0_jes_dn + AK4_JET1_jes_dn + AK4_JET2_jes_dn).M();
    
    WWTree->mass_lvjj_type0_PuppiAK4 = (LEP + NU0 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_type2_PuppiAK4 = (LEP + NU2 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_run2_PuppiAK4  = (LEP + NU1 + PuppiAK4_JET1 + PuppiAK4_JET2).M();
    WWTree->mass_lvjj_type0_met_jes_up_PuppiAK4 = (LEP + NU0_jes_up + PuppiAK4_JET1_jes_up + PuppiAK4_JET2_jes_up).M();
    WWTree->mass_lvjj_type0_met_jes_dn_PuppiAK4 = (LEP + NU0_jes_dn + PuppiAK4_JET1_jes_dn + PuppiAK4_JET2_jes_dn).M();
    

    */
    
    
    /////////////////MC Infos
    if (isMC==1)
      {
	TLorentzVector hadW, lepW, temp;
	int posWhad =-1, posWlep =-1, posGenJet=-1;
	//	std::cout<<"entry: "<<iEntry<<" "<<GenNuNum<<std::endl;
	double deltaPhiOld=100.;
	WWTree->genGravMass=100.;	

	bool foundLeptonicW = false;

	for (int i=0; i<ReducedTree->GenBosonNum; i++) {
	  if (ReducedTree->GenBoson_isBosonLeptonic[i]==1) {
	    lepW.SetPtEtaPhiE(ReducedTree->GenBosonPt[i],ReducedTree->GenBosonEta[i],ReducedTree->GenBosonPhi[i],ReducedTree->GenBosonE[i]);
	    foundLeptonicW = true;
	    posWlep = i;
	  }
	  if (foundLeptonicW) break;
	}
	
	for (int i=0; i<ReducedTree->GenBosonNum; i++) {

	    hadW.SetPtEtaPhiE(ReducedTree->GenBosonPt[i],ReducedTree->GenBosonEta[i],ReducedTree->GenBosonPhi[i],ReducedTree->GenBosonE[i]);

	    if (fabs((hadW+lepW).M()-genMass)< fabs(WWTree->genGravMass-genMass)) { //found the gen graviton
	      WWTree->genGravMass=(hadW+lepW).M();	
	      posWhad=i; //save positions of the hadronic W
	    }

	  }	

	if (posWhad!=-1 && posWlep!=-1) {

	  float oldDR=100.;
          
	  if(WWTree->event==127270357) std::cout<<"debug: "<<std::endl;

	  for (int i=0; i<ReducedTree->GenJetsAK8Num; i++) {
	      if (deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],
			 ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i])< oldDR ) 
		{
		  posGenJet=i;		
		  oldDR = deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],ReducedTree->GenJetsAK8Eta[i], ReducedTree->GenJetsAK8Phi[i]);
		  if(WWTree->event==127270357) std::cout<<"debug: had "<<ReducedTree->GenJetsAK8Phi[i]<<" "<<ReducedTree->GenBosonPhi[posWhad]<<" "<<oldDR<<std::endl;
		}
	  }



	  hadW.SetPtEtaPhiE(ReducedTree->GenBosonPt[posWhad],ReducedTree->GenBosonEta[posWhad],ReducedTree->GenBosonPhi[posWhad],ReducedTree->GenBosonE[posWhad]);
	  lepW.SetPtEtaPhiE(ReducedTree->GenBosonPt[posWlep],ReducedTree->GenBosonEta[posWlep],ReducedTree->GenBosonPhi[posWlep],ReducedTree->GenBosonE[posWlep]);

	  WWTree->W_pt_gen = ReducedTree->GenBosonPt[posWlep];
	  WWTree->W_pz_gen = lepW.Pz();
	  WWTree->W_rap_gen = lepW.Rapidity();
	  
	  WWTree->hadW_pt_gen = ReducedTree->GenBosonPt[posWhad];
	  WWTree->hadW_eta_gen = ReducedTree->GenBosonEta[posWhad];
	  WWTree->hadW_phi_gen = ReducedTree->GenBosonPhi[posWhad];
	  WWTree->hadW_e_gen = ReducedTree->GenBosonE[posWhad];
	  WWTree->hadW_m_gen = hadW.M();

	  WWTree->lepW_pt_gen = ReducedTree->GenBosonPt[posWlep];
	  WWTree->lepW_eta_gen = ReducedTree->GenBosonEta[posWlep];
	  WWTree->lepW_phi_gen = ReducedTree->GenBosonPhi[posWlep];
	  WWTree->lepW_e_gen = ReducedTree->GenBosonE[posWlep];
	  WWTree->lepW_m_gen = lepW.M();

	  WWTree->AK8_pt_gen = ReducedTree->GenJetsAK8Pt[posGenJet];
	  WWTree->AK8_eta_gen = ReducedTree->GenJetsAK8Eta[posGenJet];
	  WWTree->AK8_phi_gen = ReducedTree->GenJetsAK8Phi[posGenJet];
	  WWTree->AK8_e_gen = ReducedTree->GenJetsAK8E[posGenJet];
	  WWTree->AK8_pruned_mass_gen = ReducedTree->GenJetsAK8_prunedMass[posGenJet];
	  WWTree->AK8_softdrop_mass_gen = ReducedTree->GenJetsAK8_softdropMass[posGenJet];
	  WWTree->AK8_softdrop_pt_gen = ReducedTree->GenJetsAK8_softdropPt[posGenJet];

          if (deltaR(ReducedTree->GenBosonEta[posWhad], ReducedTree->GenBosonPhi[posWhad],
                     WWTree->ungroomed_jet_eta, WWTree->ungroomed_jet_phi)<0.1)     ok++;
          total++;
	}

	deltaPhiOld=100.;
       	for (int i=0; i<ReducedTree->GenNuNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenNuPhi[i],WWTree->v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;	  
	  temp.SetPtEtaPhiE(ReducedTree->GenNuPt[i],ReducedTree->GenNuEta[i],ReducedTree->GenNuPhi[i],ReducedTree->GenNuE[i]);
	  WWTree->nu_pz_gen=temp.Pz();	  
	  WWTree->nu_pt_gen=temp.Pt();	  
	  WWTree->nu_phi_gen=temp.Phi();	  
	  WWTree->nu_eta_gen=temp.Eta();
	  deltaPhiOld = deltaPhi;
	}		

	top_NNLO_weight[0] = 1.;
	top_NNLO_weight[1] = 1.;
        if ( TString(outputFile.c_str()).Contains("TTbar_powheg") && ReducedTree->GenTopNum > 1) { //look at here: http://arxiv.org/pdf/1511.00549.pdf
	  for (int i=0; i<2; i++) {
	    if (ReducedTree->GenTopPt[i]<60)                                           top_NNLO_weight[i] = 1./0.97;
	    else if (ReducedTree->GenTopPt[i]>=60 && ReducedTree->GenTopPt[i]<100)     top_NNLO_weight[i] = 1./0.995;
	    else if (ReducedTree->GenTopPt[i]>=100 && ReducedTree->GenTopPt[i]<150)    top_NNLO_weight[i] = 1.;
	    else if (ReducedTree->GenTopPt[i]>=150 && ReducedTree->GenTopPt[i]<200)    top_NNLO_weight[i] = 1./1.025;
	    else if (ReducedTree->GenTopPt[i]>=200 && ReducedTree->GenTopPt[i]<260)    top_NNLO_weight[i] = 1./1.045;
	    else if (ReducedTree->GenTopPt[i]>=260 && ReducedTree->GenTopPt[i]<320)    top_NNLO_weight[i] = 1./1.065;
	    else                                                                       top_NNLO_weight[i] = 1./1.095;
	  }
	}
	WWTree->top1_NNLO_Weight*=float(top_NNLO_weight[0]);
	WWTree->top2_NNLO_Weight*=float(top_NNLO_weight[1]);

	WWTree->gen_top1_pt = ReducedTree->GenTopPt[0];
	WWTree->gen_top2_pt = ReducedTree->GenTopPt[1];
      }
    
//    WWTree->totalEventWeight = WWTree->genWeight*WWTree->eff_and_pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;
    //WWTree->totalEventWeight = WWTree->genWeight*WWTree->eff_and_pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->id_eff_Weight*WWTree->trig_eff_Weight;
    WWTree->totalEventWeight = WWTree->genWeight*WWTree->eff_and_pu_Weight*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->id_eff_Weight;
    WWTree->totalEventWeight_2 = WWTree->genWeight*WWTree->eff_and_pu_Weight_2*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;
    WWTree->totalEventWeight_3 = WWTree->genWeight*WWTree->eff_and_pu_Weight_3*WWTree->top1_NNLO_Weight*WWTree->top2_NNLO_Weight*WWTree->trig_eff_Weight;

    
    bool isBadEvent=false;
    if (isMC==0) {
      std::multimap<int,int>::iterator it = badEventsList.begin(); //filter bad events
      for (it=badEventsList.begin(); it!=badEventsList.end(); ++it) {
	if (it->first == WWTree->run && it->second == WWTree->event)
	  isBadEvent = true;
      }
    }
    if (isBadEvent) continue;
    
    //    WWTree->wSampleWeight = std::atof(xSecWeight.c_str())/nEvents; //xsec/numberOfEntries
    WWTree->nEvents = nEvents-nNegEvents;
    WWTree->nNegEvents = nNegEvents;
    
    
    /////////////////FILL THE TREE
    if(WWTree->event==evento && WWTree->run==runno && WWTree->lumi==lumo) std::cout<<"fill: "<<count<<std::endl; count++;
    outTree->Fill();
  }
  std::cout << "---------end loop on events------------" << std::endl;
  std::cout << std::endl;
  
  std::cout << "----------------------" << std::endl;
  std::cout << " SUMMARY" << std::endl;
  std::cout << "----------------------" << std::endl;
  std::cout << std::endl;
  std::cout<<"MC matching: "<<(float)ok/(float)total<<std::endl;
  std::cout<<"negative events: "<<nNegEvents<<std::endl;
  std::cout << std::endl;
  std::cout<<"(0) all events:        	"<<cutEff[0]<<std::endl
	   <<"(1) tight lepton:      	"<<cutEff[1]<<std::endl
	   <<"(2) loose lepton veto: 	"<<cutEff[2]<<std::endl
	   <<"(3) MET:               	"<<cutEff[3]<<std::endl
	   <<"(4) negative lep-MET:  	"<<cutEff[4]<<std::endl
	   <<"(5) At Least 4 AK4 jets:  "<<cutEff[5]<<std::endl
	   <<"(6) VBF pair found:      	"<<cutEff[6]<<std::endl
	   <<"(7) B-tag on VBF pair:    "<<cutEff[7]<<std::endl
	   <<"(7) Found W-jets:		"<<cutEff[8]<<std::endl;
  
  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
