#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPythia6Decayer.h>
#include <TPaveStats.h>
#endif

//
// macro to compute tracking uncertainty as combination of TPC tracking uncertainty and ITS/TPC matching uncertainty
//
// based on ComputeAcceptance Macro in $ALICE_PHYSICS/PWGHF/vertexingHF/macros
//

enum EDDecay{kD0Kpi,kDplusKpipi,kDstarD0pi,kDsKKpi,kLcpKpi,kLcK0Sp,kDplusKKpi};
enum EFidY{kFixedY,kPtDepY};
enum EPtShape{kFlat,kFONLL8TeV,kFONLL8TeVfeeddown,kFONLL7TeV,kPythia7TeV,kFONLL5TeV,kFONLL13TeVprompt,kFONLL13TeVfeeddown,kPythia13TeVprompt,kPythia13TeVfeeddown};
enum ECentrality{kCentral,kSemiCentral};
enum EDataset{k18r,k18q};
enum EProtonOption{kHistogram,kFit};


// Configuration
Int_t fDDecay=kLcpKpi;
Double_t fPtMinDau=0.1;
Double_t fEtaMaxDau=0.9;
Int_t fOptionYFiducial=kPtDepY;
Double_t fYMaxFidAccCut=0.8;
Int_t fPtShape=kFONLL5TeV;
TString fDecayTableFileName="decaytable_acc.dat"; 
Int_t fDebugLevel=0;
Int_t totTrials=1000000;


Bool_t CountKpi(TClonesArray *array, Int_t nentries, Int_t &nPions, Int_t &nKaons, Int_t &nPionsInAcc, Int_t &nKaonsInAcc);
Bool_t IsInFiducialAcceptance(Double_t pt, Double_t y);
Bool_t CountPKpi(TClonesArray *array, Int_t nentries, Int_t &nPions, Int_t &nKaons, Int_t &nProtons, Int_t &nPionsInAcc, Int_t &nKaonsInAcc, Int_t &nProtonsInAcc, Int_t &idLcResChan);
Bool_t GetTrackingUnc(Int_t i, TClonesArray *array, Int_t nentries, Double_t &uncTotal, Double_t &uncTPCTotal, Double_t &uncMatchTotal, Int_t protonOption);


// Pt-shape histograms
TH1D* LoadFONLL13TeV_promptD0();
TH1D* LoadFONLL13TeV_promptDplus();
TH1D* LoadFONLL13TeV_promptDstar();
TH1D* LoadFONLL13TeV_feeddownD();
TH1D* LoadFONLL13TeV_feeddownDstar();
TH1D* LoadPYTHIA13TeV_promptD0();
TH1D* LoadPYTHIA13TeV_promptDplus();
TH1D* LoadPYTHIA13TeV_promptDstar();
TH1D* LoadPYTHIA13TeV_promptDs();
TH1D* LoadPYTHIA13TeV_feeddownD0();
TH1D* LoadPYTHIA13TeV_feeddownDplus();
TH1D* LoadPYTHIA13TeV_feeddownDstar();
TH1D* LoadPYTHIA13TeV_feeddownDs();

// --- pass 1 ---
Int_t nBins = 21;
Double_t ptBins[22] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,99};
Double_t uncPionTPC[22] = {0.21,1.4,6,6.6,5.2,3.8,2.8,2.6,2.8,2.6,1.9,2.4,0.97,2.2,2.3,2.4,2,2.3,2.9,  2,2.};
Double_t uncProtonTPC[22] = { 0.042, 0.69, 1.4, 2.8, 3.9, 3.8, 3.1, 4.2, 4.2, 4.5, 4.3, 4.7, 3.5, 3.9, 3.5, 3.6, 3.6, 3.2, 3.3,	 3.2, 3.};

//Double_t uncPionMatch[22] = {0.21,1.4,6,6.6,5.2,3.8,2.8,2.6,2.8,2.6,1.9,2.4,0.97,2.2,2.3,2.4,2,2.3,2.9,  2,2.};
//Double_t uncProtonMatch[22] = { 0.042, 0.69, 1.4, 2.8, 3.9, 3.8, 3.1, 4.2, 4.2, 4.5, 4.3, 4.7, 3.5, 3.9, 3.5, 3.6, 3.6, 3.2, 3.3,	 3.2, 3.};
//pT bins (GeV /c) 	18q piK 	18q proton 	18r piK 	18r proton 	18q piK 	18q proton 	18r piK 	18r proton
//0.5 - 1.0 	4% 	8% 	4% 	9%           	4% 	7% 	4% 	7%
//1.0 - 2.0 	5% 	4% 	5% 	5%           	4% 	6% 	4% 	5%
//2.0 - 3.0 	6% 	4% 	6% 	4%           	4% 	6% 	4% 	5%
//3.0 - 4.0 	5% 	4% 	4% 	4%           	4% 	5% 	2% 	5%
//4.0 - 5.0 	3% 	4% 	3% 	5%           	3% 	5% 	2% 	3%
//5.0 - 6.0 	2% 	4% 	2% 	3%           	2% 	3% 	2% 	2%
//6.0 - 7.0 	4% 	2% 	3% 	2%           	2% 	2% 	1% 	2%
//7.0 - 8.0 	4% 	2% 	2% 	1%           	1% 	2% 	1% 	2%
//8.0 - 9.0 	2% 	2% 	2% 	3%           	3% 	2% 	2% 	2%
//9.0 - 10.0 	3% 	3% 	2% 	2%           	2% 	1% 	2% 	2%
//10.0 - 11.0 	2% 	3% 	2% 	4%        	2% 	1% 	2% 	1%
//11.0 - 12.0 	2% 	2% 	2% 	3%        	2% 	1% 	2% 	1%
//12.0 - 13.0 	2% 	3% 	2% 	2%        	3% 	2% 	2% 	2%
//13.0 - 14.0 	2% 	3% 	2% 	2%        	2% 	1% 	2% 	1%
Int_t nBinsMatch = 14;
Double_t ptBinsMatch[15] = {0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
// --- pass 1 ---
//Double_t uncPionMatch_010_q[14] = {4,5,6,5,3,2,4,4,2,3,2,2,2,2};
//Double_t uncProtonMatch_010_q[14] = {8,4,4,4,4,4,2,2,2,3,3,2,3,3};
//Double_t uncPionMatch_010_r[14] = {4,5,6,4,3,2,3,2,2,2,2,2,2,2};
//Double_t uncProtonMatch_010_r[14] = {0,5,4,4,5,3,2,1,3,2,4,3,2,2};
//Double_t uncPionMatch_3050_q[14] = {4,4,4,4,3,2,2,1,3,2,2,2,3,2};
//Double_t uncProtonMatch_3050_q[14] = {7,6,6,5,5,3,2,2,2,1,1,1,2,1};
//Double_t uncPionMatch_3050_r[14] = {4,4,4,2,2,2,1,1,2,2,2,2,2,2};
//Double_t uncProtonMatch_3050_r[14] = {7,5,5,5,3,2,2,2,2,2,1,1,2,1};
//Double_t uncPionMatch[14] = {0}; //to fill
//Double_t uncProtonMatch[14] = {0}; //to fill

// --- pass 3 ---
Double_t uncPionMatch_010_q[14] = {3,5,6,5,4,2,3,3,3,3,2,2,1,2};
Double_t uncProtonMatch_010_q[14] = {7,4,6,7,5,3,4,2,2,2,1,1,1,1};

Double_t uncPionMatch_010_r[14] = {3,5,6,5,4,2,3,2,3,2,2,2,1,1};
Double_t uncProtonMatch_010_r[14] = {10,4,6,7,6,3,4,2,3,2,1,2,1,1};

Double_t uncPionMatch_3050_q[14] = {3,4,4,4,3,2,1,1,1,1,1,1,1,1};
Double_t uncProtonMatch_3050_q[14] = {10,4,6,5,4,2,2,1,1,1,1,1,1,1};

Double_t uncPionMatch_3050_r[14] = {3,4,4,4,3,1,2,1,2,1,1,1,1,1};
Double_t uncProtonMatch_3050_r[14] = {10,4,6,5,5,2,3,1,1,1,1,1,1,1};

Double_t uncPionMatch[14] = {0}; //to fill
Double_t uncProtonMatch[14] = {0}; //to fill



//Double_t uncPionTPC[22] = {0.21,1.4,6,6.6,5.2,3.8,2.8,2.6,2.8,2.6,1.9,2.4,0.97,4.2,7.3,2.4,0,4.3,5.9,  10,2.};
//Double_t uncProtonTPC[22] = { 0.042, 0.69, 1.4, 2.8, 3.9, 3.8, 3.1, 4.2, 5.2, 3.5, 3.3, 3.7, 1.5, 1.9, 2.5, 3.6, 3.6, 12, 5.3,	 6.2, 3.};
//Double_t uncPionMatch[22] = {0.21,1.4,6,6.6,5.2,3.8,2.8,2.6,2.8,2.6,1.9,2.4,0.97,4.2,7.3,2.4,0,4.3,5.9,  10,2.};
//Double_t uncProtonMatch[22] = { 0.042, 0.69, 1.4, 2.8, 3.9, 3.8, 3.1, 4.2, 5.2, 3.5, 3.3, 3.7, 1.5, 1.9, 2.5, 3.6, 3.6, 12, 5.3,	 6.2, 3.};
TH1D *hTPCProton = 0x0;
TH1D *hTPCPion = 0x0;
TH1D *hMatchProton = new TH1D("hMatchProton","proton Match unnc.",nBinsMatch,ptBinsMatch);
TH1D *hMatchPion = new TH1D("hMatchPion","pion Match unnc.",nBinsMatch,ptBinsMatch);

TH2D *hLcProtonPt = new TH2D("hLcProtonPt","Lc pT vs proton pT",60,0,30,60,0,30);
TH2D *hLcPionPt = new TH2D("hLcPionPt","Lc pT vs pion pT",60,0,30,60,0,30);
TH2D *hLcKaonPt = new TH2D("hLcKaonPt","Lc pT vs kaon pT",60,0,30,60,0,30);

Int_t nPtBinsLc = 7;
Double_t ptBinsLc[8] = {0,1,2,4,6,8,12,24};
Double_t uncBinsLc[501] = {0};
TH2D *hUncTrackingTotal = 0x0;
TH2D *hUncTrackingTotalBinned = 0x0;

void ComputeUncertainty(Int_t dataset = k18r, Int_t centrality = kSemiCentral, Int_t protonOption = kFit,  TString pathTPC = ""){

  for(Int_t i=0;i<501;i++) uncBinsLc[i] = Double_t(i)*0.1;
  TH2D *hUncTrackingTotal = new TH2D("hUncTrackingTotal","total uncertainty vs Lc pT",300,0,30,500,0,50);
  hUncTrackingTotal->GetYaxis()->SetTitle("total uncertainty (%)"); 
  hUncTrackingTotal->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 
  TH2D *hUncTrackingTotalBinned = new TH2D("hUncTrackingTotalBinned","total uncertainty vs Lc pT", nPtBinsLc, ptBinsLc,500,uncBinsLc);
  hUncTrackingTotalBinned->GetYaxis()->SetTitle("total uncertainty (%)"); 
  hUncTrackingTotalBinned->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 
  TH2D *hUncTrackingTPCTotal = new TH2D("hUncTrackingTPCTotal","total uncertainty vs Lc pT",300,0,30,500,0,50);
  hUncTrackingTPCTotal->GetYaxis()->SetTitle("total uncertainty (%)"); 
  hUncTrackingTPCTotal->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 
  TH2D *hUncTrackingTPCTotalBinned = new TH2D("hUncTrackingTPCTotalBinned","total uncertainty vs Lc pT", nPtBinsLc, ptBinsLc,500,uncBinsLc);
  hUncTrackingTPCTotalBinned->GetYaxis()->SetTitle("total uncertainty (%)"); 
  hUncTrackingTPCTotalBinned->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 
  TH2D *hUncTrackingMatchTotal = new TH2D("hUncTrackingMatchTotal","total uncertainty vs Lc pT",300,0,30,500,0,50);
  hUncTrackingMatchTotal->GetYaxis()->SetTitle("total uncertainty (%)"); 
  hUncTrackingMatchTotal->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 
  TH2D *hUncTrackingMatchTotalBinned = new TH2D("hUncTrackingMatchTotalBinned","total uncertainty vs Lc pT", nPtBinsLc, ptBinsLc,500,uncBinsLc);
  hUncTrackingMatchTotalBinned->GetYaxis()->SetTitle("total uncertainty (%)"); 
  hUncTrackingMatchTotalBinned->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 





  hLcProtonPt->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 
  hLcProtonPt->GetYaxis()->SetTitle("proton p_{T} (GeV/c)"); 
  hLcPionPt->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 
  hLcPionPt->GetYaxis()->SetTitle("pion p_{T} (GeV/c)"); 
  hLcKaonPt->GetXaxis()->SetTitle("#Lambda_{c} p_{T} (GeV/c)"); 
  hLcKaonPt->GetYaxis()->SetTitle("pion p_{T} (GeV/c)"); 

  TFile *fTPCpi = 0x0;
  TFile *fTPCpr = 0x0;
  // ITS-TPC matching
  if(dataset==k18r) {
    if(centrality==kCentral)  {
     fTPCpi = new TFile(Form("%s/TrackingTPCCutUnc_cent_18r_pionCuts_3sigmaTPC.root",pathTPC.Data()));
     fTPCpr = new TFile(Form("%s/TrackingTPCCutUnc_cent_18r_protonCuts_3sigmaTPC.root",pathTPC.Data()));
      for(Int_t i=0;i<14;i++) {
        uncPionMatch[i] = uncPionMatch_010_r[i];
        uncProtonMatch[i] = uncProtonMatch_010_r[i];
      }
    }
    if(centrality==kSemiCentral)  {
     fTPCpi = new TFile(Form("%s/TrackingTPCCutUnc_semicent_18r_pionCuts_3sigmaTPC.root",pathTPC.Data()));
     fTPCpr = new TFile(Form("%s/TrackingTPCCutUnc_semicent_18r_protonCuts_3sigmaTPC.root",pathTPC.Data()));
      for(Int_t i=0;i<14;i++) {
        uncPionMatch[i] = uncPionMatch_3050_r[i];
        uncProtonMatch[i] = uncProtonMatch_3050_r[i];
      }
    }
  }
  if(dataset==k18q) {
    if(centrality==kCentral)  {
      fTPCpi = new TFile(Form("%s/TrackingTPCCutUnc_cent_18q_pionCuts_3sigmaTPC.root",pathTPC.Data()));
      fTPCpr = new TFile(Form("%s/TrackingTPCCutUnc_cent_18q_protonCuts_3sigmaTPC.root",pathTPC.Data()));
      for(Int_t i=0;i<14;i++) {
        uncPionMatch[i] = uncPionMatch_010_q[i];
        uncProtonMatch[i] = uncProtonMatch_010_q[i];
      }
    }
    if(centrality==kSemiCentral)  {
     fTPCpi = new TFile(Form("%s/TrackingTPCCutUnc_semicent_18q_pionCuts_3sigmaTPC.root",pathTPC.Data()));
     fTPCpr = new TFile(Form("%s/TrackingTPCCutUnc_semicent_18q_protonCuts_3sigmaTPC.root",pathTPC.Data()));
      for(Int_t i=0;i<14;i++) {
        uncPionMatch[i] = uncPionMatch_3050_q[i];
        uncProtonMatch[i] = uncProtonMatch_3050_q[i];
      }
    }
  }
  hTPCProton = (TH1D*)fTPCpr->Get("hunc_proton");
  hTPCPion = (TH1D*)fTPCpi->Get("hunc_pion");

  for(Int_t i=0;i<hTPCProton->GetNbinsX();i++) {
    if(hTPCProton->GetBinContent(i+1)>30) hTPCProton->SetBinContent(i+1,0.);
  }
  hTPCProton->Fit("pol1"); 

  // Fill histogram
  for(Int_t i=0;i<nBins;i++) {
  //  hTPCProton->SetBinContent(i+1,uncProtonTPC[i]);
  //  hTPCPion->SetBinContent(i+1,uncPionTPC[i]);
    hMatchProton->SetBinContent(i+1,uncProtonMatch[i]);
    hMatchPion->SetBinContent(i+1,uncPionMatch[i]);
  }
  // main function
  
  gSystem->Load("liblhapdf.so");      // Parton density functions
  gSystem->Load("libEGPythia6.so");   // TGenerator interface
  gSystem->Load("libpythia6.so");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations

  TPythia6Decayer* pdec=TPythia6Decayer::Instance();
  if(fDecayTableFileName.CompareTo("")!=0){
    if(fDecayTableFileName.Contains("ALICE_PHYSICS")){
      gSystem->Exec(Form("cp %s .",fDecayTableFileName.Data()));
      fDecayTableFileName.ReplaceAll("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/","./");
    }
    pdec->SetDecayTableFile(fDecayTableFileName.Data());
    pdec->ReadDecayTable();
  }
  pdec->Init();

  Int_t pdgCode=0;
  Int_t nPionDau=-1;
  Int_t nProtonDau=-1;
  Int_t nKaonDau=-1;
  TString outFileName="Acceptance_Toy_";
  if(fDDecay==kD0Kpi){
    pdgCode=421;
    nPionDau=1;
    nKaonDau=1;
    nProtonDau=0;
    outFileName.Append("D0Kpi_");
  }else if(fDDecay==kDplusKpipi){
    pdgCode=411;
    nPionDau=2;
    nKaonDau=1;
    nProtonDau=0;
    outFileName.Append("DplusKpipi_");
  }else if(fDDecay==kDstarD0pi){
    pdgCode=413;
    nPionDau=2;
    nKaonDau=1;
    nProtonDau=0;
    outFileName.Append("DStarD0pi_");
  }else if(fDDecay==kDsKKpi){
    pdgCode=431;
    nPionDau=1;
    nKaonDau=2;
    nProtonDau=0;
    outFileName.Append("DsKKpi_");
  }else if(fDDecay==kLcpKpi){
    pdgCode=4122;
    nPionDau=1;
    nKaonDau=1;
    nProtonDau=1;
    outFileName.Append("LcpKpi_");
  }else if(fDDecay==kLcK0Sp){
    pdgCode=4122;
    nPionDau=2;
    nKaonDau=0;
    nProtonDau=1;
    outFileName.Append("LcK0Sp_");
  }else if(fDDecay==kDplusKKpi){
    pdgCode=411;
    nPionDau=1;
    nKaonDau=2;
    nProtonDau=0;
    outFileName.Append("DplusKKpi_");
  }else{
    printf("ERROR: Wrong decay selected\n");
    return;
  }
  if(fOptionYFiducial==kFixedY) outFileName.Append(Form("yfid%02d_",(Int_t)(fYMaxFidAccCut*10)));
  else outFileName.Append("yfidPtDep_");
  TString periodCent = "";
  if(dataset==k18r) {
    periodCent.Append("18r_");
    if(centrality==kCentral) {
      periodCent.Append("central_");
    }
    if(centrality==kSemiCentral) {
      periodCent.Append("semiCentral_");
    }
  }
  if(dataset==k18q) {
    periodCent.Append("18q_");
    if(centrality==kCentral) {
      periodCent.Append("central_");
    }
    if(centrality==kSemiCentral) {
      periodCent.Append("semiCentral_");
    }
  }
  if(protonOption==kFit)  periodCent.Append("protonFit_");
  outFileName.Append(periodCent);


  outFileName.Append(Form("etaDau%02d_",(Int_t)(fEtaMaxDau*10)));
  outFileName.Append(Form("ptDau%d_",(Int_t)(fPtMinDau*1000)));
  TDatabasePDG* db=TDatabasePDG::Instance();
  Float_t massD=db->GetParticle(pdgCode)->Mass();
  TClonesArray *array = new TClonesArray("TParticle",100);

  TH2D* hPtVsYGen=new TH2D("hPtVsYGen","",400,0.,40.,20.,-1.,1.);
  hPtVsYGen->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hPtVsYGen->GetYaxis()->SetTitle("y");


  TF1*  funcPt=0x0;
  TH1D* histPt=0x0;
  if(fPtShape==kFONLL8TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.518046,3.01138,3.38914,1.75899); // Prompt
    outFileName.Append("FONLL8ptshape.root");
  }else if(fPtShape==kFONLL8TeVfeeddown){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.398252, 3.9603, 3.915, 1.51853); // FeedDown
    outFileName.Append("FONLL8ptshapeFeedDown.root");
  }else if(fPtShape==kFONLL7TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.322643,2.96275,2.30301,2.5);
    outFileName.Append("FONLL7ptshape.root");
  }else if(fPtShape==kFONLL5TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.302879,2.9750,3.68139,1.68855);
    outFileName.Append("FONLL5ptshape.root");
  }else if(fPtShape==kPythia7TeV){
    funcPt=new TF1("fFONLL","[0]*x/TMath::Power((1+TMath::Power(x/[1],[3])),[2])",0.,40.);
    funcPt->SetParameters(0.322643,1.94635,1.40463,2.5);
    outFileName.Append("PYTHIA7ptshape.root");  
  }else if(fPtShape==kFONLL13TeVprompt){
    if(fDDecay==kDplusKpipi){
      histPt = LoadFONLL13TeV_promptDplus();
      outFileName.Append("promptDplus");
    }else if (fDDecay==kDstarD0pi){
      histPt = LoadFONLL13TeV_promptDstar();
      outFileName.Append("promptDstar");
    }else{
      histPt = LoadFONLL13TeV_promptD0();
      outFileName.Append("promptD0");
    }
    outFileName.Append("FONLL13ptshape.root");

  }else if (fPtShape==kFONLL13TeVfeeddown){
    if (fDDecay==kDstarD0pi){
      histPt = LoadFONLL13TeV_feeddownDstar();
      outFileName.Append("feeddownDstar");
    }else{
      histPt = LoadFONLL13TeV_feeddownD();
      outFileName.Append("feeddownD");
    }
    outFileName.Append("FONLL13ptshape.root");
  }else if(fPtShape==kPythia13TeVprompt){
    if(fDDecay==kDplusKpipi){
      histPt = LoadPYTHIA13TeV_promptDplus();
      outFileName.Append("promptDplus");
    }else if (fDDecay==kDstarD0pi){
      histPt = LoadPYTHIA13TeV_promptDstar();
      outFileName.Append("promptDstar");
    }else if (fDDecay==kDsKKpi){
      histPt = LoadPYTHIA13TeV_promptDs();
      outFileName.Append("promptDs");
    }else{
      histPt = LoadPYTHIA13TeV_promptD0();
      outFileName.Append("promptD0");
    }
    outFileName.Append("PYTHIA13ptshape.root");
  }else if (fPtShape==kPythia13TeVfeeddown){
    if(fDDecay==kDplusKpipi){
      histPt = LoadPYTHIA13TeV_feeddownDplus();
      outFileName.Append("feeddownDplus");
    }else if (fDDecay==kDstarD0pi){
      histPt = LoadPYTHIA13TeV_feeddownDstar();
      outFileName.Append("feeddownDstar");
    }else if (fDDecay==kDsKKpi){
      histPt = LoadPYTHIA13TeV_feeddownDs();
      outFileName.Append("feeddownDs");
    }else{
      histPt = LoadPYTHIA13TeV_feeddownD0();
      outFileName.Append("feeddownD0");
    }
    outFileName.Append("PYTHIA13ptshape.root");
  }else{
    funcPt=new TF1("fFlat","pol0",0.,40.);
    funcPt->SetParameter(0,1.);
    outFileName.Append("flatpt.root");
  }

  if (funcPt) funcPt->SetNpx(10000);

  TRandom3* gener=new TRandom3(0);
  TLorentzVector* vec=new TLorentzVector();


  for(Int_t itry=0; itry<totTrials; itry++){
    if(itry%10000==0) printf("Event %d\n",itry);
    Float_t ptD = funcPt ? funcPt->GetRandom() : histPt->GetRandom();
    Float_t phiD=gener->Rndm()*2*TMath::Pi();
    Float_t yD=gener->Rndm()*2.-1.; // flat in -1<y<1
    Float_t px=ptD*TMath::Cos(phiD);
    Float_t py=ptD*TMath::Sin(phiD);
    Float_t mt=TMath::Sqrt(massD*massD+ptD*ptD);
    Float_t pz=mt*TMath::SinH(yD);
    Float_t E=TMath::Sqrt(massD*massD+px*px+py*py+pz*pz);

    // TLorentzVector* vec=new TLorentzVector(px,py,pz,E);
    vec->SetPxPyPzE(px,py,pz,E);
    pdec->Decay(pdgCode,vec);
    array->Clear();
    Int_t nentries = pdec->ImportParticles(array);
    TParticle* dmes=(TParticle*)array->At(0);
    Int_t nDaughters=dmes->GetNDaughters();
    if(fDDecay==kD0Kpi && nDaughters!=2) continue;
    if(fDDecay==kLcK0Sp && nentries>6) continue;
    Int_t nPionsInAcc=0;
    Int_t nProtonsInAcc=0;
    Int_t nKaonsInAcc=0;
    Int_t nPions=0;
    Int_t nProtons=0;
    Int_t nKaons=0;
    Int_t idLcResChan=0; // non resonant by default;
    Bool_t isOk=CountPKpi(array,nentries,nPions,nKaons,nProtons,nPionsInAcc,nKaonsInAcc,nProtonsInAcc, idLcResChan);
    Double_t uncTotal = 0;
    Double_t uncTPCTotal = 0;
    Double_t uncMatchTotal = 0;
    Bool_t ok = GetTrackingUnc(itry,array,nentries, uncTotal,uncTPCTotal,uncMatchTotal,protonOption);
    if(!ok) continue;
    hUncTrackingTotal->Fill(ptD, uncTotal);
    hUncTrackingTotalBinned->Fill(ptD, uncTotal);
    hUncTrackingTPCTotal->Fill(ptD, uncTPCTotal);
    hUncTrackingTPCTotalBinned->Fill(ptD, uncTPCTotal);
    hUncTrackingMatchTotal->Fill(ptD, uncMatchTotal);
    hUncTrackingMatchTotalBinned->Fill(ptD, uncMatchTotal);

    //if(itry%10000==0) printf("Event %d\n",itry);

    if(isOk){
      if(nPions==nPionDau && nKaons==nKaonDau && nProtons==nProtonDau){
	hPtVsYGen->Fill(ptD,yD);
      }
    }
    // delete vec;
  }

  TCanvas *c1 = new TCanvas();
  c1->SetLogz(); 
  hUncTrackingTotalBinned->GetXaxis()->SetRangeUser(4,24);
  hUncTrackingTotalBinned->Draw("COLZ");

  TCanvas *c2 = new TCanvas();
  c2->SetLogz(); 
  hLcProtonPt->GetXaxis()->SetRangeUser(4,24);
  hLcProtonPt->Draw("COLZ");

  TCanvas *c3 = new TCanvas();
  c3->SetLogz(); 
  hLcPionPt->GetXaxis()->SetRangeUser(4,24);
  hLcPionPt->Draw("COLZ");

  TCanvas *c4 = new TCanvas();
  TH1D *hUncVsPt = (TH1D*)hUncTrackingTotal->ProfileX();
  TH1D *hUncVsPtBinned = (TH1D*)hUncTrackingTotalBinned->ProfileX();
  hUncVsPt->GetYaxis()->SetTitle("total uncertainty (%)");
  hUncVsPt->Draw(); 
  hUncVsPtBinned->SetLineColor(kRed);
  hUncVsPtBinned->Draw("SAME");

  TCanvas *c5 = new TCanvas();
  c3->SetLogz(); 
  hLcKaonPt->GetXaxis()->SetRangeUser(4,24);
  hLcKaonPt->Draw("COLZ");

  TH1D *hUncTPCVsPt = (TH1D*)hUncTrackingTPCTotal->ProfileX();
  TH1D *hUncTPCVsPtBinned = (TH1D*)hUncTrackingTPCTotalBinned->ProfileX();
  TH1D *hUncMatchVsPt = (TH1D*)hUncTrackingMatchTotal->ProfileX();
  TH1D *hUncMatchVsPtBinned = (TH1D*)hUncTrackingMatchTotalBinned->ProfileX();

  TFile* outfil=new TFile(outFileName.Data(),"recreate");
  hPtVsYGen->Write();
  hTPCPion->Write();
  hTPCProton->Write();
  hMatchPion->Write();
  hMatchProton->Write();
  hUncTrackingTotal->Write();
  hUncTrackingTotalBinned->Write();
  hUncTrackingTPCTotal->Write();
  hUncTrackingTPCTotalBinned->Write();
  hUncTrackingMatchTotal->Write();
  hUncTrackingMatchTotalBinned->Write();
  hUncVsPt->Write();
  hUncVsPtBinned->Write();
  hUncTPCVsPt->Write();
  hUncTPCVsPtBinned->Write();
  hUncMatchVsPt->Write();
  hUncMatchVsPtBinned->Write();
  hLcProtonPt->Write();
  hLcPionPt->Write();
  hLcKaonPt->Write();
  c4->Write();
  c1->SaveAs(Form("Uncertainty_2d_%s.eps",periodCent.Data()));
  c2->SaveAs(Form("Pt_Lc_proton_%s.eps",periodCent.Data()));
  c3->SaveAs(Form("Pt_Lc_pion_%s.eps",periodCent.Data()));
  c5->SaveAs(Form("Pt_Lc_kaon_%s.eps",periodCent.Data()));
  c4->SaveAs(Form("Uncertainty_%s.eps",periodCent.Data()));

  outfil->Close();

}

//___________________________________________________
Bool_t IsInFiducialAcceptance(Double_t pt, Double_t y){
  // check fiducial acceptance

  if(fOptionYFiducial==kFixedY){
    if(TMath::Abs(y) > fYMaxFidAccCut) return kFALSE;
    else return kTRUE;
  }

  if(pt > 5.) {
    if (TMath::Abs(y) > 0.8) return kFALSE;
  } else {
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
    if (y < minFiducialY || y > maxFiducialY) return kFALSE;
  }

  return kTRUE;
}

Bool_t GetTrackingUnc(Int_t i, TClonesArray *array, Int_t nentries, Double_t &uncTotal, Double_t &uncTPCTotal, Double_t &uncMatchTotal, Int_t protonOption) {

  TParticle* dmes=(TParticle*)array->At(0);
  Double_t sumPx=0;
  Double_t sumPy=0;
  Double_t sumPz=0;
  Double_t nProtons = 0;
  Double_t nK0s = 0;
  Double_t nPions = 0;
  Double_t nKaons = 0;
  //cout<<"n daughters = "<<nentries<<endl;
  Double_t ptLc = dmes->Pt();
  for(int j=0; j<nentries; j++){
    TParticle * o = (TParticle*)array->At(j);
    Int_t pdgdau=TMath::Abs(o->GetPdgCode());
    //cout<<"pdg daughter "<<j<<" = "<<pdgdau<<endl;
    if(pdgdau==2212) { // proton
      //cout<<"found protonn"<<endl;
      nProtons++;
      Double_t ptProton = o->Pt();
      if(ptProton<0.5) return kFALSE; // cut on proton pt
      Double_t uncMatch = hMatchProton->GetBinContent(hMatchProton->FindBin(ptProton));
      Double_t uncTPC = 0;
      if(protonOption==kHistogram) uncTPC = hTPCProton->GetBinContent(hTPCProton->FindBin(ptProton));
      else if(protonOption==kFit) uncTPC = hTPCProton->GetFunction("pol1")->Eval(ptProton);
      uncMatchTotal += uncMatch;
      uncTPCTotal += uncTPC;
      uncTotal += TMath::Sqrt(uncMatch*uncMatch + uncTPC*uncTPC);
      hLcProtonPt->Fill(ptLc,ptProton);
      if(i%10000==0) {
        cout<<"pt proton = "<<ptProton<<endl;
        cout<<"unc proton = "<<uncMatch<<"(match) + "<<uncTPC<<"(TPC)"<<endl;
      }
    }
    if(pdgdau==211) { // pion
      //cout<<"found protonn"<<endl;
      nPions++;
      Double_t ptPion= o->Pt();
      if(ptPion<0.5) return kFALSE; // cut on proton pt
      Double_t uncMatch = hMatchPion->GetBinContent(hMatchPion->FindBin(ptPion));
      Double_t uncTPC = hTPCPion->GetBinContent(hTPCPion->FindBin(ptPion));
      uncMatchTotal += uncMatch;
      uncTPCTotal += uncTPC;
      uncTotal += TMath::Sqrt(uncMatch*uncMatch + uncTPC*uncTPC);
      hLcPionPt->Fill(ptLc,ptPion);
      if(i%10000==0) {
        cout<<"pt pion = "<<ptPion<<endl;
        cout<<"unc pion = "<<uncMatch<<"(match) + "<<uncTPC<<"(TPC)"<<endl;
      }
    }
    if(pdgdau==321) { // kaon
      //cout<<"found protonn"<<endl;
      nKaons++;
      Double_t ptKaon= o->Pt();
      if(ptKaon<0.5) return kFALSE; // cut on proton pt
      Double_t uncMatch = hMatchPion->GetBinContent(hMatchPion->FindBin(ptKaon));
      Double_t uncTPC = hTPCPion->GetBinContent(hTPCPion->FindBin(ptKaon));
      uncMatchTotal += uncMatch;
      uncTPCTotal += uncTPC;
      uncTotal += TMath::Sqrt(uncMatch*uncMatch + uncTPC*uncTPC);
      hLcKaonPt->Fill(ptLc,ptKaon);
      if(i%10000==0) {
        cout<<"pt kaon = "<<ptKaon<<endl;
        cout<<"unc kaon = "<<uncMatch<<"(match) + "<<uncTPC<<"(TPC)"<<endl;
      }
    }


    if(pdgdau==310) { //K0s
      //cout<<"found K0s"<<endl;
      nK0s++;
      TParticle *daugh1 = (TParticle*)array->At(j+1);
      if(!daugh1) {
        cout<<"nno first daughter of K0s found"<<endl;
        return kFALSE;
        //continue;
      }
      //cout<<"K0s n daughters = "<<o->GetNDaughters()<<endl;
      TParticle *daugh2 = (TParticle*)array->At(j+2);
      if(!daugh2) {
        cout<<"no second daughter of K0s found"<<endl;
        return kFALSE;
        //continue;
      }
      Double_t pdgKdaugh1 = TMath::Abs(daugh1->GetPdgCode());
      Double_t pdgKdaugh2 = TMath::Abs(daugh2->GetPdgCode());
      if(pdgKdaugh1 != 211 || pdgKdaugh2 != 211) {
        cout<<"k0 not decaying to pions"<<endl;
        return kFALSE;
      }
      Double_t ptKdaugh1 = daugh1->Pt();
      Double_t ptKdaugh2 = daugh2->Pt();
      if(ptKdaugh1<0.5) return kFALSE; // cut on pion pt
      if(ptKdaugh2<0.5) return kFALSE; // cut on pion pt

      //Double_t uncMatch1 = hMatchPion->GetBinContent(hMatchPion->FindBin(ptKdaugh1));
      //Double_t uncMatch2 = hMatchPion->GetBinContent(hMatchPion->FindBin(ptKdaugh2));
      Double_t uncMatch1 = 0.;
      Double_t uncMatch2 = 0.;
      Double_t uncTPC1 = hTPCPion->GetBinContent(hTPCPion->FindBin(ptKdaugh1));
      Double_t uncTPC2 = hTPCPion->GetBinContent(hTPCPion->FindBin(ptKdaugh2));
      uncMatchTotal += uncMatch1 + uncMatch2;
      uncTPCTotal += uncTPC1 + uncTPC2;
      uncTotal += TMath::Sqrt(uncMatch1*uncMatch1 + uncTPC1*uncTPC1) + TMath::Sqrt(uncMatch2*uncMatch2 + uncTPC2*uncTPC2 );
      hLcPionPt->Fill(ptLc,ptKdaugh1);
      hLcPionPt->Fill(ptLc,ptKdaugh2);
      if(i%10000==0) {
        cout<<"pt pion 1 = "<<ptKdaugh1<<"\tpt pion 2 = "<<ptKdaugh2<<endl;
        cout<<"unc pion 1 = "<<uncMatch1<<"(match) + "<<uncTPC1<<"(TPC)"<<endl;
        cout<<"unc pion 2 = "<<uncMatch2<<"(match) + "<<uncTPC2<<"(TPC)"<<endl;
        cout<<"unc pion total = "<<TMath::Sqrt(uncMatch1*uncMatch1 + uncTPC1*uncTPC1 + uncMatch2*uncMatch2 + uncTPC2*uncTPC2 )<<endl;
        cout<<"unc pion total (linear) = "<<TMath::Sqrt(uncMatch1*uncMatch1 + uncTPC1*uncTPC1) + TMath::Sqrt(uncMatch2*uncMatch2 + uncTPC2*uncTPC2 )<<endl;
      }
    }
//    if(pdgdau==211) {//pion
//      //cout<<"found Pion"<<endl;
//      nPi++;
//      Double_t ptKdaugh1 = o->Pt();
//
//      Double_t uncMatch1 = hMatchPion->GetBinContent(hMatchPion->FindBin(ptKdaugh1));
//      Double_t uncTPC1= hTPCPion->GetBinContent(hTPCPion->FindBin(ptKdaugh1));
//      uncTotal += TMath::Sqrt(uncMatch1*uncMatch1 + uncTPC1*uncTPC1);
//    }
  }
  if((( nProtons!=1 || nK0s!=1 ) && (nProtons!=1 || nPions!=2)) && (nProtons!=1 || nPions!=1 || nKaons!=1) ) {
    cout<<" decay products not right"<<endl;
    return kFALSE;
  }
  //cout<<" unc total = "<<uncTotal<<endl;
  return kTRUE;
}

//___________________________________________________
Bool_t CountKpi(TClonesArray *array, Int_t nentries, Int_t &nPions, Int_t &nKaons, Int_t &nPionsInAcc, Int_t &nKaonsInAcc){
  // count K and pi in Acc

  TParticle* dmes=(TParticle*)array->At(0);
  Double_t sumPx=0;
  Double_t sumPy=0;
  Double_t sumPz=0;
  for(int j=0; j<nentries; j++){
    TParticle * o = (TParticle*)array->At(j);
    Int_t pdgdau=TMath::Abs(o->GetPdgCode());
    if(fDebugLevel>0) printf("%d ",pdgdau);
    if(pdgdau==130) {
      if(fDebugLevel>0) printf("K0 dacaying into K0L\n");
      return kFALSE;
    }
    Float_t ptdau=TMath::Sqrt(o->Px()*o->Px()+o->Py()*o->Py());      
    Float_t etadau=o->Eta();
    if(pdgdau==211){ 
      nPions++;
      sumPx+=o->Px();
      sumPy+=o->Py();
      sumPz+=o->Pz();
    }
    if(pdgdau==321){ 
      nKaons++;
      sumPx+=o->Px();
      sumPy+=o->Py();
      sumPz+=o->Pz();
    }
    if(TMath::Abs(etadau)<fEtaMaxDau && ptdau>fPtMinDau){
      if(pdgdau==211) nPionsInAcc++;
      if(pdgdau==321) nKaonsInAcc++;
    }
  }
  if(fDebugLevel>0) printf("\n");
  if(TMath::Abs(sumPx-dmes->Px())>0.001 ||
     TMath::Abs(sumPy-dmes->Py())>0.001 ||
     TMath::Abs(sumPz-dmes->Pz())>0.001){
    printf("Momentum conservation violation\n");
    return kFALSE;
  }
  return kTRUE;
}


//___________________________________________________
Bool_t CountPKpi(TClonesArray *array, Int_t nentries, Int_t &nPions, Int_t &nKaons, Int_t &nProtons, Int_t &nPionsInAcc, Int_t &nKaonsInAcc, Int_t &nProtonsInAcc, Int_t &idLcResChan){
  // count K and pi in Acc

  TParticle* dmes=(TParticle*)array->At(0);
  Double_t sumPx=0;
  Double_t sumPy=0;
  Double_t sumPz=0;
  
  for(int j=0; j<nentries; j++){
    TParticle * o = (TParticle*)array->At(j);
    Int_t pdgdau=TMath::Abs(o->GetPdgCode());
    if(fDebugLevel>0) printf("%d ",pdgdau);
    if(pdgdau==130) {
      if(fDebugLevel>0) printf("K0 dacaying into K0L\n");
      return kFALSE;
    }
    Float_t ptdau=TMath::Sqrt(o->Px()*o->Px()+o->Py()*o->Py());      
    Float_t etadau=o->Eta();
    if(pdgdau==211){ 
      nPions++;
      sumPx+=o->Px();
      sumPy+=o->Py();
      sumPz+=o->Pz();
    }
    if(pdgdau==321){ 
      nKaons++;
      sumPx+=o->Px();
      sumPy+=o->Py();
      sumPz+=o->Pz();
    }
    if(pdgdau==2212){ 
      nProtons++;
      sumPx+=o->Px();
      sumPy+=o->Py();
      sumPz+=o->Pz();
    }
    if(pdgdau==313) idLcResChan=2; //K*0 
    else if(pdgdau==2224) idLcResChan=3; //Delta++
    else if(pdgdau==3124) idLcResChan=1;  //Lambda1520
    if(TMath::Abs(etadau)<fEtaMaxDau && ptdau>fPtMinDau){
      if(pdgdau==211) nPionsInAcc++;
      if(pdgdau==321) nKaonsInAcc++;
      if(pdgdau==2212) nProtonsInAcc++;
    }
  }
  if(fDebugLevel>0) printf("\n");
  if(TMath::Abs(sumPx-dmes->Px())>0.001 ||
     TMath::Abs(sumPy-dmes->Py())>0.001 ||
     TMath::Abs(sumPz-dmes->Pz())>0.001){
    printf("Momentum conservation violation\n");
    return kFALSE;
  }
  return kTRUE;
}



//___________________________________________________
TH1D* LoadFONLL13TeV_promptD0()
{
  TH1D *hFONLL13 = new TH1D("hFONLL13TeV_D0", "", 80, 0., 40.);
  Float_t val[80] = {
    1.4686e+08, 3.9542e+08, 5.4901e+08, 5.2166e+08, 4.1083e+08, 2.9968e+08, 2.1299e+08, 1.5057e+08, 1.0701e+08, 7.6919e+07,
    5.6121e+07, 4.1546e+07, 3.1184e+07, 2.3715e+07, 1.8253e+07, 1.4206e+07, 1.1175e+07, 8.8774e+06, 7.1169e+06, 5.7544e+06,
    4.6899e+06, 3.8509e+06, 3.1841e+06, 2.6499e+06, 2.2189e+06, 1.8687e+06, 1.5823e+06, 1.3467e+06, 1.1516e+06, 9.8933e+05,
    8.5356e+05, 7.3942e+05, 6.4302e+05, 5.6122e+05, 4.9154e+05, 4.3193e+05, 3.8074e+05, 3.3662e+05, 2.9846e+05, 2.6535e+05,
    2.3653e+05, 2.1136e+05, 1.8932e+05, 1.6997e+05, 1.5293e+05, 1.3789e+05, 1.2458e+05, 1.1277e+05, 1.0227e+05, 9.2913e+04,
    8.4560e+04, 7.7087e+04, 7.0388e+04, 6.4371e+04, 5.8956e+04, 5.4074e+04, 4.9667e+04, 4.5680e+04, 4.2068e+04, 3.8791e+04,
    3.5813e+04, 3.3103e+04, 3.0633e+04, 2.8380e+04, 2.6321e+04, 2.4436e+04, 2.2710e+04, 2.1127e+04, 1.9673e+04, 1.8336e+04,
    1.7106e+04, 1.5972e+04, 1.4926e+04, 1.3961e+04, 1.3069e+04, 1.2243e+04, 1.1479e+04, 1.0771e+04, 1.0113e+04, 9.5031e+03
  };
  for (Int_t ibin=0; ibin<80; ++ibin) hFONLL13->SetBinContent(ibin+1, val[ibin]);

  return hFONLL13;
}



//___________________________________________________
TH1D* LoadFONLL13TeV_promptDplus()
{
  TH1D *hFONLL13 = new TH1D("hFONLL13TeV_Dplus", "", 80, 0., 40.);
  Float_t val[80] = {
    1.5242e+08, 3.9396e+08, 5.3767e+08, 5.1263e+08, 4.0706e+08, 2.9926e+08, 2.1414e+08, 1.5230e+08, 1.0879e+08, 7.8507e+07,
    5.7460e+07, 4.2651e+07, 3.2093e+07, 2.4462e+07, 1.8865e+07, 1.4708e+07, 1.1588e+07, 9.2185e+06, 7.3998e+06, 5.9900e+06,
    4.8871e+06, 4.0167e+06, 3.3240e+06, 2.7686e+06, 2.3200e+06, 1.9552e+06, 1.6566e+06, 1.4107e+06, 1.2070e+06, 1.0374e+06,
    8.9550e+05, 7.7610e+05, 6.7519e+05, 5.8953e+05, 5.1652e+05, 4.5404e+05, 4.0037e+05, 3.5409e+05, 3.1405e+05, 2.7929e+05,
    2.4902e+05, 2.2258e+05, 1.9943e+05, 1.7909e+05, 1.6117e+05, 1.4535e+05, 1.3135e+05, 1.1892e+05, 1.0787e+05, 9.8024e+04,
    8.9230e+04, 8.1359e+04, 7.4302e+04, 6.7963e+04, 6.2256e+04, 5.7112e+04, 5.2465e+04, 4.8261e+04, 4.4452e+04, 4.0996e+04,
    3.7854e+04, 3.4995e+04, 3.2390e+04, 3.0012e+04, 2.7838e+04, 2.5849e+04, 2.4027e+04, 2.2355e+04, 2.0820e+04, 1.9408e+04,
    1.8108e+04, 1.6910e+04, 1.5805e+04, 1.4785e+04, 1.3842e+04, 1.2969e+04, 1.2161e+04, 1.1411e+04, 1.0716e+04, 1.0070e+04
  };
  for (Int_t ibin=0; ibin<80; ++ibin) hFONLL13->SetBinContent(ibin+1, val[ibin]);

  return hFONLL13;
}



//___________________________________________________
TH1D* LoadFONLL13TeV_promptDstar()
{
  TH1D *hFONLL13 = new TH1D("hFONLL13TeV_Dstar", "", 80, 0., 40.);
  Float_t val[80] = {
    1.2433e+08, 3.4512e+08, 5.0662e+08, 5.1020e+08, 4.2016e+08, 3.1661e+08, 2.3064e+08, 1.6632e+08, 1.2007e+08, 8.7340e+07,
    6.4329e+07, 4.8000e+07, 3.6287e+07, 2.7772e+07, 2.1492e+07, 1.6808e+07, 1.3278e+07, 1.0589e+07, 8.5178e+06, 6.9084e+06,
    5.6462e+06, 4.6479e+06, 3.8519e+06, 3.2124e+06, 2.6951e+06, 2.2738e+06, 1.9285e+06, 1.6437e+06, 1.4076e+06, 1.2108e+06,
    1.0459e+06, 9.0711e+05, 7.8968e+05, 6.8993e+05, 6.0484e+05, 5.3197e+05, 4.6933e+05, 4.1529e+05, 3.6850e+05, 3.2786e+05,
    2.9246e+05, 2.6152e+05, 2.3441e+05, 2.1058e+05, 1.8958e+05, 1.7104e+05, 1.5461e+05, 1.4003e+05, 1.2706e+05, 1.1550e+05,
    1.0517e+05, 9.5920e+04, 8.7626e+04, 8.0172e+04, 7.3461e+04, 6.7409e+04, 6.1940e+04, 5.6992e+04, 5.2508e+04, 4.8437e+04,
    4.4736e+04, 4.1367e+04, 3.8296e+04, 3.5492e+04, 3.2929e+04, 3.0583e+04, 2.8433e+04, 2.6460e+04, 2.4648e+04, 2.2981e+04,
    2.1446e+04, 2.0032e+04, 1.8727e+04, 1.7522e+04, 1.6407e+04, 1.5376e+04, 1.4420e+04, 1.3535e+04, 1.2713e+04, 1.1950e+04
  };
  for (Int_t ibin=0; ibin<80; ++ibin) hFONLL13->SetBinContent(ibin+1, val[ibin]);

  return hFONLL13;
}



//___________________________________________________
TH1D* LoadFONLL13TeV_feeddownD()
{ // B->D with B.R=1
  TH1D *hFONLL13 = new TH1D("FONLL13TeV_feeddownD", "", 80, 0., 40.);
  Float_t val[80] = {
    1.0310e+07, 2.6790e+07, 3.4480e+07, 3.4430e+07, 3.0200e+07, 2.4740e+07, 1.9600e+07, 1.5300e+07, 1.1880e+07, 9.2260e+06,
    7.1950e+06, 5.6470e+06, 4.4650e+06, 3.5580e+06, 2.8570e+06, 2.3120e+06, 1.8860e+06, 1.5490e+06, 1.2810e+06, 1.0660e+06,
    8.9260e+05, 7.5170e+05, 6.3650e+05, 5.4160e+05, 4.6310e+05, 3.9770e+05, 3.4300e+05, 2.9700e+05, 2.5820e+05, 2.2520e+05,
    1.9710e+05, 1.7310e+05, 1.5240e+05, 1.3460e+05, 1.1920e+05, 1.0590e+05, 9.4290e+04, 8.4160e+04, 7.5290e+04, 6.7500e+04,
    6.0650e+04, 5.4610e+04, 4.9270e+04, 4.4530e+04, 4.0320e+04, 3.6580e+04, 3.3240e+04, 3.0250e+04, 2.7570e+04, 2.5170e+04,
    2.3020e+04, 2.1070e+04, 1.9320e+04, 1.7740e+04, 1.6310e+04, 1.5010e+04, 1.3830e+04, 1.2760e+04, 1.1790e+04, 1.0900e+04,
    1.0090e+04, 9.3520e+03, 8.6760e+03, 8.0560e+03, 7.4890e+03, 6.9670e+03, 6.4880e+03, 6.0470e+03, 5.6410e+03, 5.2670e+03,
    4.9220e+03, 4.6030e+03, 4.3080e+03, 4.0350e+03, 3.7830e+03, 3.5480e+03, 3.3310e+03, 3.1290e+03, 2.9410e+03, 2.7670e+03
  };
  for (Int_t ibin=0; ibin<80; ++ibin) hFONLL13->SetBinContent(ibin+1, val[ibin]);

  return hFONLL13;
}



//___________________________________________________
TH1D* LoadFONLL13TeV_feeddownDstar()
{ // B->D* with B.R=1
  TH1D *hFONLL13 = new TH1D("FONLL13TeV_feeddownDstar", "", 80, 0., 40.);
  Float_t val[80] = {
    9.5260e+06, 2.5070e+07, 3.2890e+07, 3.3540e+07, 3.0010e+07, 2.5000e+07, 2.0090e+07, 1.5860e+07, 1.2430e+07, 9.7280e+06,
    7.6360e+06, 6.0240e+06, 4.7820e+06, 3.8240e+06, 3.0800e+06, 2.5000e+06, 2.0430e+06, 1.6810e+06, 1.3930e+06, 1.1610e+06,
    9.7390e+05, 8.2120e+05, 6.9610e+05, 5.9300e+05, 5.0750e+05, 4.3620e+05, 3.7650e+05, 3.2620e+05, 2.8370e+05, 2.4760e+05,
    2.1680e+05, 1.9050e+05, 1.6780e+05, 1.4830e+05, 1.3140e+05, 1.1680e+05, 1.0400e+05, 9.2870e+04, 8.3100e+04, 7.4530e+04,
    6.6990e+04, 6.0330e+04, 5.4440e+04, 4.9220e+04, 4.4580e+04, 4.0440e+04, 3.6760e+04, 3.3460e+04, 3.0510e+04, 2.7860e+04,
    2.5470e+04, 2.3330e+04, 2.1390e+04, 1.9640e+04, 1.8060e+04, 1.6630e+04, 1.5320e+04, 1.4140e+04, 1.3060e+04, 1.2080e+04,
    1.1190e+04, 1.0370e+04, 9.6190e+03, 8.9330e+03, 8.3050e+03, 7.7270e+03, 7.1970e+03, 6.7080e+03, 6.2580e+03, 5.8440e+03,
    5.4610e+03, 5.1080e+03, 4.7810e+03, 4.4790e+03, 4.1980e+03, 3.9390e+03, 3.6970e+03, 3.4740e+03, 3.2650e+03, 3.0720e+03
  };
  for (Int_t ibin=0; ibin<80; ++ibin) hFONLL13->SetBinContent(ibin+1, val[ibin]);

  return hFONLL13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_promptD0()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13TeV_D0", "", 40, 0., 40.);
  Float_t val[40] = {
    2617743, 3763836, 2235903, 1140807, 587707, 317881, 180756, 107762, 67299,42882,
    28550, 19690, 13817, 9754, 7205, 5299, 4100, 3177, 2386,1822,
    1484, 1248, 933, 798, 639, 534, 435, 405, 304,260,
    210, 196, 161, 143, 139, 99, 89, 80, 66,63
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_promptDplus()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13TeV_Dplus", "", 40, 0., 40.);
  Float_t val[40] = {
    1192379, 1736283, 1047485, 540499, 281523, 152153, 87150, 52200, 32709, 20991,
    13986, 9611, 6832, 4705, 3484, 2650, 2001, 1518, 1197, 924,
    735, 576, 497, 402, 292, 261, 214, 173, 169, 150,
    106, 102, 95, 70, 67, 47, 53, 46, 35, 35
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_promptDstar()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13TeV_Dstar", "", 40, 0., 40.);
  Float_t val[40] = {
    922985, 1419749, 904587, 484407, 257862, 142386, 82177, 49839, 31103, 20422,
    13493, 9272, 6549, 4746, 3344, 2575, 1959, 1493, 1161, 935,
    672, 585, 456, 415, 328, 260, 199, 180, 161, 134,
    109, 97, 103, 72, 63, 56, 40, 46, 38, 27
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_promptDs()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13TeV_Ds", "", 40, 0., 40.);
  Float_t val[40] = {
    346381, 519143, 320488, 167944, 87484, 47325, 26932, 16376, 10058, 6527,
    4347, 3041, 2112, 1521, 1069, 849, 609, 468, 359, 275,
    226, 197, 152, 124, 111, 90, 63, 60, 52, 46,
    39, 35, 26, 29, 16, 20, 10, 13, 18, 9
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_feeddownD0()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13_feeddownD0", "", 40, 0., 40.);
  Float_t val[40] = {
    1814030, 2939561, 2059721, 1197596, 684753, 397435, 237432, 146442, 93029, 60222,
    40940, 28861, 19992, 14344, 10602, 7742, 5781, 4395, 3421, 2732,
    2216, 1763, 1395, 1208, 888, 735, 599, 509, 418, 374,
    302, 281, 242, 200, 174, 158, 133, 104, 98, 88
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_feeddownDplus()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13_feeddownDplus", "", 40, 0., 40.);
  Float_t val[40] = {
    779389, 1268000, 899979, 527442, 302782, 177610, 105677, 66601, 42008, 27618,
    18802, 13024, 8975, 6530, 4773, 3495, 2708, 2020, 1634, 1212,
    1017, 809, 665, 507, 389, 355, 285, 250, 218, 180,
    143, 137, 106, 85, 71, 64, 57, 43, 50, 47
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_feeddownDstar()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13_feeddownDstar", "", 40, 0., 40.);
  Float_t val[40] = {
    665735, 1144214, 856647, 523704, 309061, 183569, 111994, 70187, 45266, 29601,
    20079, 14042, 10338, 7184, 5234, 4002, 2935, 2301, 1788, 1349,
    1088, 895, 748, 591, 480, 378, 302, 263, 216, 182,
    172, 122, 133, 100, 92, 87, 65, 50, 42, 45
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}



//___________________________________________________
TH1D* LoadPYTHIA13TeV_feeddownDs()
{
  TH1D *hPYTHIA13 = new TH1D("hPYTHIA13_feeddownDs", "", 40, 0., 40.);
  Float_t val[40] = {
    377925, 686435, 518802, 318819, 189396, 113651, 68963, 43048, 27776, 18177,
    12356, 8643, 6100, 4512, 3170, 2479, 1835, 1424, 1086, 825,
    659, 490, 465, 361, 308, 246, 181, 170, 117, 111,
    112, 91, 69, 66, 54, 57, 39, 33, 31, 25
  };
  for (Int_t ibin=0; ibin<40; ++ibin) hPYTHIA13->SetBinContent(ibin+1, val[ibin]);

  return hPYTHIA13;
}
