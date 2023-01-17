#include <TString.h>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLatex.h>
#include "iostream"
#include "iomanip"

enum tagflags {
  kIsPionFromK0s       = BIT(0),
  kIsPionFromL         = BIT(1),
  kIsProtonFromL       = BIT(2),
  kIsElectronFromGamma = BIT(3),
  kIsKaonFromKinks     = BIT(4),
  kIsKaonFromTOF       = BIT(5),
  kIsKaonFromTPC       = BIT(6),
  kIsKaonFromHMPID      = BIT(7),
  kIsDeuteronFromTPCTOF = BIT(8),
  kIsTritonFromTPCTOF   = BIT(9),
  kIsHe3FromTPCTOF      = BIT(10),
  kPositiveTrack        = BIT(14),
  kNegativeTrack        = BIT(15),
  kIsKaonFromOmega      = BIT(16),
  kIsPionFromTOF        = BIT(17),
  kIsProtonFromTOF      = BIT(18),
};

enum centrality {
  kCentral,kSemiCentral,kpp,kpPb
};
enum dataset{
  k18r, k18q, kBoth,k17pq,k16qt
};
enum speciescuts{
  kPionCuts,
  kProtonCuts,
  kAllCuts
};

// 
// Macro which loops over data and MC track trees, and checks variation in 
// efficiencies when varying track selection. Variations defined by DPG
//
// Note - Cut2 corresponds to crossed rows cut but this is no longer used in
// Pb-Pb 2018 data since MC/data description is not good - the cut is still
// done in this macro but not used to determine final uncertainty
//


void GetTrackingSyst(Int_t cent = kSemiCentral, Int_t dataset = kBoth , Int_t speciescuts = kPionCuts, TString fname_data = "", TString fname_mc = "", Bool_t isNsigmaCut = kTRUE) {

  Bool_t isTest = kFALSE;
  TString dirname_data, dirname_mc;
  TString datasetcent = "";
  if(cent==kCentral) {
    datasetcent+="0-10%";
    if(dataset==k18r) {
      dirname_data = Form("PWGHF_D2H_SystNsigmaPID_PbPb_010_kCentral_kINT7");
      dirname_mc = Form("PWGHF_D2H_SystNsigmaPID_PbPb_010_kMB");
      datasetcent+=", LHC18r";
    }
    else if(dataset==k18q) {
      dirname_data = Form("PWGHF_D2H_SystNsigmaPID_PbPb_010_kCentral_kINT7");
      dirname_mc = Form("PWGHF_D2H_SystNsigmaPID_PbPb_010_kMB");
      datasetcent+=", LHC18q";
    }
    else if(dataset==kBoth) {
      dirname_data = Form("PWGHF_D2H_SystNsigmaPID_PbPb_010_kCentral_kINT7");
      dirname_mc = Form("PWGHF_D2H_SystNsigmaPID_PbPb_010_kMB");
      datasetcent+=", LHC18q+LHC18r";
    }
  }
  if(cent==kSemiCentral) {
    datasetcent+="30-50%";
    if(dataset==k18r) {
      dirname_data = Form("PWGHF_D2H_SystNsigmaPID_PbPb_3050_kSemiCentral_kINT7");
      dirname_mc = Form("PWGHF_D2H_SystNsigmaPID_PbPb_3050_kMB");
      datasetcent+=", LHC18r";
    }
    else if(dataset==k18q) {
      dirname_data = Form("PWGHF_D2H_SystNsigmaPID_PbPb_3050_kSemiCentral_kINT7");
      dirname_mc = Form("PWGHF_D2H_SystNsigmaPID_PbPb_3050_kMB");
      datasetcent+=", LHC18q";
    }
    else if(dataset==kBoth) {
      dirname_data = Form("PWGHF_D2H_SystNsigmaPID_PbPb_3050_kSemiCentral_kINT7");
      dirname_mc = Form("PWGHF_D2H_SystNsigmaPID_PbPb_3050_kMB");
      datasetcent+=", LHC18q+LHC18r";
    }

  }
 if(cent==kpp) {
     datasetcent+="pp 5 TeV";
     dirname_data = Form("PWGHF_D2H_SystNsigmaPID_ppMB_kINT7");
     dirname_mc = Form("PWGHF_D2H_SystNsigmaPID_ppMB_kINT7");
     datasetcent+=", LHC17pq";
   }
 if(cent==kpPb) {
     datasetcent+="pPb 5 TeV";
     dirname_data = Form("PWGHF_D2H_SystNsigmaPID_pPb_kINT7");
     dirname_mc = Form("PWGHF_D2H_SystNsigmaPID_pPb_kINT7");
     datasetcent+=", LHC16qt";
   }

  TFile *f_data = new TFile(fname_data.Data());
  TFile *f_mc = new TFile(fname_mc.Data());
  TDirectoryFile *dir_data = (TDirectoryFile*)f_data->Get(dirname_data.Data());
  TDirectoryFile *dir_mc = (TDirectoryFile*)f_mc->Get(dirname_mc.Data());
  if(!dir_data || !dir_mc) {
    cout<<"can't find dir nname"<<endl;
    f_data->ls();
    f_mc->ls();
    return;
  }

  TTree *t_data = (TTree*)dir_data->Get("fPIDtree"); 
  TTree *t_mc = (TTree*)dir_mc->Get("fPIDtree"); 


  t_data->Print();

  // branch names tree data
  short n_sigma_TPC_pi_data, n_sigma_TPC_pi_mc;
  short n_sigma_TPC_K_data, n_sigma_TPC_K_mc;
  short n_sigma_TPC_p_data, n_sigma_TPC_p_mc;
  short n_sigma_TOF_p_data, n_sigma_TOF_p_mc;
  short n_sigma_TOF_K_data, n_sigma_TOF_K_mc;
  short n_sigma_TOF_pi_data, n_sigma_TOF_pi_mc;
  unsigned short pT_data, pT_mc;
  UChar_t NclusterTPC_data, NclusterTPC_mc;
  UChar_t NclusterPIDTPC_data, NclusterPIDTPC_mc;
  UChar_t NcrossedRowsTPC_data, NcrossedRowsTPC_mc;
  UChar_t NFindableTPC_data, NFindableTPC_mc;
  unsigned int tag_data, tag_mc;

  t_data->SetBranchAddress("pT",&pT_data);
  t_data->SetBranchAddress("NclusterTPC",&NclusterTPC_data);
  t_data->SetBranchAddress("NclusterPIDTPC",&NclusterPIDTPC_data);
  t_data->SetBranchAddress("NcrossedRowsTPC",&NcrossedRowsTPC_data);
  t_data->SetBranchAddress("NFindableTPC",&NFindableTPC_data);
  t_data->SetBranchAddress("n_sigma_TPC_pi",&n_sigma_TPC_pi_data);
  t_data->SetBranchAddress("n_sigma_TPC_p",&n_sigma_TPC_p_data);
  t_data->SetBranchAddress("n_sigma_TPC_K",&n_sigma_TPC_K_data);
  t_data->SetBranchAddress("n_sigma_TOF_pi",&n_sigma_TOF_pi_data);
  t_data->SetBranchAddress("n_sigma_TOF_p",&n_sigma_TOF_p_data);
  t_data->SetBranchAddress("n_sigma_TOF_K",&n_sigma_TOF_K_data);
  t_data->SetBranchAddress("tag",&tag_data);

  t_mc->SetBranchAddress("pT",&pT_mc);
  t_mc->SetBranchAddress("NclusterTPC",&NclusterTPC_mc);
  t_mc->SetBranchAddress("NclusterPIDTPC",&NclusterPIDTPC_mc);
  t_mc->SetBranchAddress("NcrossedRowsTPC",&NcrossedRowsTPC_mc);
  t_mc->SetBranchAddress("NFindableTPC",&NFindableTPC_mc);
  t_mc->SetBranchAddress("n_sigma_TPC_pi",&n_sigma_TPC_pi_mc);
  t_mc->SetBranchAddress("n_sigma_TPC_p",&n_sigma_TPC_p_mc);
  t_mc->SetBranchAddress("n_sigma_TPC_K",&n_sigma_TPC_K_mc);
  t_mc->SetBranchAddress("n_sigma_TOF_pi",&n_sigma_TOF_pi_mc);
  t_mc->SetBranchAddress("n_sigma_TOF_p",&n_sigma_TOF_p_mc);
  t_mc->SetBranchAddress("n_sigma_TOF_K",&n_sigma_TOF_K_mc);
  t_mc->SetBranchAddress("tag",&tag_mc);

  Int_t nbins = 12;
  Double_t bins[13] =  {0.,0.5,1.,1.5,2.,2.5,3,3.5, 4,5,6,8,15};

  TH1D *h_proton_all_data = new TH1D("h_proton_all_data","proton sample, no TPC cut",nbins,bins);
  TH1D *h_proton_Cut1_data = new TH1D("h_proton_Cut1_data","proton sample, cut 1",nbins,bins);
  TH1D *h_proton_Cut2_data = new TH1D("h_proton_Cut2_data","proton sample, cut 2",nbins,bins);
  TH1D *h_proton_Cut3_data = new TH1D("h_proton_Cut3_data","proton sample, cut 3",nbins,bins);
  TH1D *h_proton_Cut4_data = new TH1D("h_proton_Cut4_data","proton sample, cut 4",nbins,bins);
  TH1D *h_proton_Cut5_data = new TH1D("h_proton_Cut5_data","proton sample, cut 5",nbins,bins);
  TH1D *h_proton_Cut6_data = new TH1D("h_proton_Cut6_data","proton sample, cut 6",nbins,bins);
  h_proton_all_data->Sumw2();
  h_proton_Cut1_data->Sumw2();
  h_proton_Cut2_data->Sumw2();
  h_proton_Cut3_data->Sumw2();
  h_proton_Cut4_data->Sumw2();
  h_proton_Cut5_data->Sumw2();
  h_proton_Cut6_data->Sumw2();

  TH1D *h_kaon_all_data = new TH1D("h_kaon_all_data","kaon sample, no TPC cut",nbins,bins);
  TH1D *h_kaon_Cut1_data = new TH1D("h_kaon_Cut1_data","kaon sample, cut 1",nbins,bins);
  TH1D *h_kaon_Cut2_data = new TH1D("h_kaon_Cut2_data","kaon sample, cut 2",nbins,bins);
  TH1D *h_kaon_Cut3_data = new TH1D("h_kaon_Cut3_data","kaon sample, cut 3",nbins,bins);
  TH1D *h_kaon_Cut4_data = new TH1D("h_kaon_Cut4_data","kaon sample, cut 4",nbins,bins);
  TH1D *h_kaon_Cut5_data = new TH1D("h_kaon_Cut5_data","kaon sample, cut 5",nbins,bins);
  TH1D *h_kaon_Cut6_data = new TH1D("h_kaon_Cut6_data","kaon sample, cut 6",nbins,bins);
  h_kaon_all_data->Sumw2();
  h_kaon_Cut1_data->Sumw2();
  h_kaon_Cut2_data->Sumw2();
  h_kaon_Cut3_data->Sumw2();
  h_kaon_Cut4_data->Sumw2();
  h_kaon_Cut5_data->Sumw2();
  h_kaon_Cut6_data->Sumw2();

  TH1D *h_pion_all_data = new TH1D("h_pion_all_data","pion sample, no TPC cut",nbins,bins);
  TH1D *h_pion_Cut1_data = new TH1D("h_pion_Cut1_data","pion sample, cut 1",nbins,bins);
  TH1D *h_pion_Cut2_data = new TH1D("h_pion_Cut2_data","pion sample, cut 2",nbins,bins);
  TH1D *h_pion_Cut3_data = new TH1D("h_pion_Cut3_data","pion sample, cut 3",nbins,bins);
  TH1D *h_pion_Cut4_data = new TH1D("h_pion_Cut4_data","pion sample, cut 4",nbins,bins);
  TH1D *h_pion_Cut5_data = new TH1D("h_pion_Cut5_data","pion sample, cut 5",nbins,bins);
  TH1D *h_pion_Cut6_data = new TH1D("h_pion_Cut6_data","pion sample, cut 6",nbins,bins);
  h_pion_all_data->Sumw2();
  h_pion_Cut1_data->Sumw2();
  h_pion_Cut2_data->Sumw2();
  h_pion_Cut3_data->Sumw2();
  h_pion_Cut4_data->Sumw2();
  h_pion_Cut5_data->Sumw2();
  h_pion_Cut6_data->Sumw2();

  TH1D *h_charged_all_data = new TH1D("h_charged_all_data","charged sample, no TPC cut",nbins,bins);
  TH1D *h_charged_Cut1_data = new TH1D("h_charged_Cut1_data","charged sample, cut 1",nbins,bins);
  TH1D *h_charged_Cut2_data = new TH1D("h_charged_Cut2_data","charged sample, cut 2",nbins,bins);
  TH1D *h_charged_Cut3_data = new TH1D("h_charged_Cut3_data","charged sample, cut 3",nbins,bins);
  TH1D *h_charged_Cut4_data = new TH1D("h_charged_Cut4_data","charged sample, cut 4",nbins,bins);
  TH1D *h_charged_Cut5_data = new TH1D("h_charged_Cut5_data","charged sample, cut 5",nbins,bins);
  TH1D *h_charged_Cut6_data = new TH1D("h_charged_Cut6_data","charged sample, cut 6",nbins,bins);
  h_charged_all_data->Sumw2();
  h_charged_Cut1_data->Sumw2();
  h_charged_Cut2_data->Sumw2();
  h_charged_Cut3_data->Sumw2();
  h_charged_Cut4_data->Sumw2();
  h_charged_Cut5_data->Sumw2();
  h_charged_Cut6_data->Sumw2();

  TH1D *h_proton_all_mc = new TH1D("h_proton_all_mc","proton sample, no TPC cut",nbins,bins);
  TH1D *h_proton_Cut1_mc = new TH1D("h_proton_Cut1_mc","proton sample, cut 1",nbins,bins);
  TH1D *h_proton_Cut2_mc = new TH1D("h_proton_Cut2_mc","proton sample, cut 2",nbins,bins);
  TH1D *h_proton_Cut3_mc = new TH1D("h_proton_Cut3_mc","proton sample, cut 3",nbins,bins);
  TH1D *h_proton_Cut4_mc = new TH1D("h_proton_Cut4_mc","proton sample, cut 4",nbins,bins);
  TH1D *h_proton_Cut5_mc = new TH1D("h_proton_Cut5_mc","proton sample, cut 5",nbins,bins);
  TH1D *h_proton_Cut6_mc = new TH1D("h_proton_Cut6_mc","proton sample, cut 6",nbins,bins);
  h_proton_all_mc->Sumw2();
  h_proton_Cut1_mc->Sumw2();
  h_proton_Cut2_mc->Sumw2();
  h_proton_Cut3_mc->Sumw2();
  h_proton_Cut4_mc->Sumw2();
  h_proton_Cut5_mc->Sumw2();
  h_proton_Cut6_mc->Sumw2();

  TH1D *h_kaon_all_mc = new TH1D("h_kaon_all_mc","kaon sample, no TPC cut",nbins,bins);
  TH1D *h_kaon_Cut1_mc = new TH1D("h_kaon_Cut1_mc","kaon sample, cut 1",nbins,bins);
  TH1D *h_kaon_Cut2_mc = new TH1D("h_kaon_Cut2_mc","kaon sample, cut 2",nbins,bins);
  TH1D *h_kaon_Cut3_mc = new TH1D("h_kaon_Cut3_mc","kaon sample, cut 3",nbins,bins);
  TH1D *h_kaon_Cut4_mc = new TH1D("h_kaon_Cut4_mc","kaon sample, cut 4",nbins,bins);
  TH1D *h_kaon_Cut5_mc = new TH1D("h_kaon_Cut5_mc","kaon sample, cut 5",nbins,bins);
  TH1D *h_kaon_Cut6_mc = new TH1D("h_kaon_Cut6_mc","kaon sample, cut 6",nbins,bins);
  h_kaon_all_mc->Sumw2();
  h_kaon_Cut1_mc->Sumw2();
  h_kaon_Cut2_mc->Sumw2();
  h_kaon_Cut3_mc->Sumw2();
  h_kaon_Cut4_mc->Sumw2();
  h_kaon_Cut5_mc->Sumw2();
  h_kaon_Cut6_mc->Sumw2();

  TH1D *h_pion_all_mc = new TH1D("h_pion_all_mc","pion sample, no TPC cut",nbins,bins);
  TH1D *h_pion_Cut1_mc = new TH1D("h_pion_Cut1_mc","pion sample, cut 1",nbins,bins);
  TH1D *h_pion_Cut2_mc = new TH1D("h_pion_Cut2_mc","pion sample, cut 2",nbins,bins);
  TH1D *h_pion_Cut3_mc = new TH1D("h_pion_Cut3_mc","pion sample, cut 3",nbins,bins);
  TH1D *h_pion_Cut4_mc = new TH1D("h_pion_Cut4_mc","pion sample, cut 4",nbins,bins);
  TH1D *h_pion_Cut5_mc = new TH1D("h_pion_Cut5_mc","pion sample, cut 5",nbins,bins);
  TH1D *h_pion_Cut6_mc = new TH1D("h_pion_Cut6_mc","pion sample, cut 6",nbins,bins);
  h_pion_all_mc->Sumw2();
  h_pion_Cut1_mc->Sumw2();
  h_pion_Cut2_mc->Sumw2();
  h_pion_Cut3_mc->Sumw2();
  h_pion_Cut4_mc->Sumw2();
  h_pion_Cut5_mc->Sumw2();
  h_pion_Cut6_mc->Sumw2();

  TH1D *h_charged_all_mc = new TH1D("h_charged_all_mc","charged sample, no TPC cut",nbins,bins);
  TH1D *h_charged_Cut1_mc = new TH1D("h_charged_Cut1_mc","charged sample, cut 1",nbins,bins);
  TH1D *h_charged_Cut2_mc = new TH1D("h_charged_Cut2_mc","charged sample, cut 2",nbins,bins);
  TH1D *h_charged_Cut3_mc = new TH1D("h_charged_Cut3_mc","charged sample, cut 3",nbins,bins);
  TH1D *h_charged_Cut4_mc = new TH1D("h_charged_Cut4_mc","charged sample, cut 4",nbins,bins);
  TH1D *h_charged_Cut5_mc = new TH1D("h_charged_Cut5_mc","charged sample, cut 5",nbins,bins);
  TH1D *h_charged_Cut6_mc = new TH1D("h_charged_Cut6_mc","charged sample, cut 6",nbins,bins);
  h_charged_all_mc->Sumw2();
  h_charged_Cut1_mc->Sumw2();
  h_charged_Cut2_mc->Sumw2();
  h_charged_Cut3_mc->Sumw2();
  h_charged_Cut4_mc->Sumw2();
  h_charged_Cut5_mc->Sumw2();
  h_charged_Cut6_mc->Sumw2();


  TH1D *hunc_proton = new TH1D("hunc_proton","ucertainty proton",nbins,bins);
  TH1D *hunc_pion = new TH1D("hunc_pion","ucertainty pion",nbins,bins);
  TH1D *hunc_kaon = new TH1D("hunc_kaon","ucertainty kaon",nbins,bins);
  TH1D *hunc_charged = new TH1D("hunc_charged","ucertainty charged",nbins,bins);

  //
  // cuts which should be tested ( https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsTrackSystematicUncertainty )
  //
  //  additional cut on number TPC crossed rows > 120-(5/pt); (standard cut =70, no pt dependent);
  //  number of TPC clusters >0.65 x number of TPC crossed rows;
  //  ratio of crossed rows over findable clusters in the TPC >0.9. (standard cut =0.8)
  //  additional cut on number of clusters with TPC dE/dx signal > 0.5 x number of TPC crossed rows;
  //


  double cutClusterPercentOfCrossedRows = 0.75;
  double cutCrossedOverFindable = 0.9;
  double cutClusterWithdEdxPercentOfCrossedRows= 0.5;

  short nSigmaID = 10;

  //
  // Data
  //

  // loop over candidates 
  for(int i=0;i<t_data->GetEntries();i++) {
    if(i%1000000 == 0) cout<<"-- track "<<i<<endl;
    if(isTest && i>10000000) break;

    t_data->GetEntry(i);
    // convert to better format
    double pT = double( pT_data) / 1000.;
    double NclusterTPC = double ( NclusterTPC_data) ;
    double NclusterPIDTPC = double ( NclusterPIDTPC_data) ;
    double NcrossedRowsTPC = double ( NcrossedRowsTPC_data) ;
    double NFindableTPC = double ( NFindableTPC_data) ;

    if(pT<0.00001) continue;
    if(NFindableTPC < 1) continue;
    double cutPtDepCrossedRows = 120. - (5./pT);
    // is proton
    if((tag_data & (kIsProtonFromL | kIsProtonFromTOF))
        && ((isNsigmaCut && TMath::Abs(n_sigma_TPC_p_data) < 300) || !isNsigmaCut)) { 
    //  cout<<"is proton pT = "<<pT<<"\t crossed rows = "<<NcrossedRowsTPC<<"\tfindable clusters = "<<NFindableTPC<<"\tclusters = "<<NclusterTPC<<endl;

      // default cut
      //if(NcrossedRowsTPC > 70. && NclusterTPC > 50. && NcrossedRowsTPC / NFindableTPC > 0.8) { // pass1: used n cluseters but now not used as not described by MC
        h_proton_all_data->Fill(pT);
        // additional cut variations
        if(NcrossedRowsTPC > cutPtDepCrossedRows)  // pt dependent crossed rows
          h_proton_Cut1_data->Fill(pT);
        if(NclusterTPC > cutClusterPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of crossed rows
          h_proton_Cut2_data->Fill(pT);
        if(NcrossedRowsTPC / NFindableTPC > cutCrossedOverFindable)  // crossed over findable
          h_proton_Cut3_data->Fill(pT);
        if(NclusterPIDTPC  > cutClusterWithdEdxPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of number of clusters with signal 
          h_proton_Cut4_data->Fill(pT);
        if(NclusterPIDTPC > 40)  // N clusters with TPC signal > 40
          h_proton_Cut5_data->Fill(pT);
        if(NclusterPIDTPC > 60)  // N clusters with TPC signal > 60
          h_proton_Cut6_data->Fill(pT);
    }
    // is pion
    if((tag_data & (kIsPionFromK0s | kIsPionFromL | kIsPionFromTOF))  
        && ((isNsigmaCut && TMath::Abs(n_sigma_TPC_pi_data) < 300) || !isNsigmaCut)) { 
     // cout<<"is pion pT = "<<pT<<"\t crossed rows = "<<NcrossedRowsTPC<<"\tfindable clusters = "<<NFindableTPC<<"\tclusters = "<<NclusterTPC<<endl;

      // default cut
        h_pion_all_data->Fill(pT);
        // additional cut variations
        if(NcrossedRowsTPC > cutPtDepCrossedRows)  // pt dependent crossed rows
          h_pion_Cut1_data->Fill(pT);
        if(NclusterTPC > cutClusterPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of crossed rows
          h_pion_Cut2_data->Fill(pT);
        if(NcrossedRowsTPC / NFindableTPC > cutCrossedOverFindable)  // crossed over findable
          h_pion_Cut3_data->Fill(pT);
        if(NclusterPIDTPC > cutClusterWithdEdxPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of number of clusters with signal 
          h_pion_Cut4_data->Fill(pT);
        if(NclusterPIDTPC > 40)  // N clusters with TPC signal > 40
          h_pion_Cut5_data->Fill(pT);
        if(NclusterPIDTPC > 60)  // N clusters with TPC signal > 60
          h_pion_Cut6_data->Fill(pT);
    }
    // is pion
    if(tag_data & kIsKaonFromTOF)  {
      // cout<<"is pion pT = "<<pT<<"\t crossed rows = "<<NcrossedRowsTPC<<"\tfindable clusters = "<<NFindableTPC<<"\tclusters = "<<NclusterTPC<<endl;

      // default cut
      h_kaon_all_data->Fill(pT);
      // additional cut variations
      if(NcrossedRowsTPC > cutPtDepCrossedRows)  // pt dependent crossed rows
        h_kaon_Cut1_data->Fill(pT);
      if(NclusterTPC > cutClusterPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of crossed rows
        h_kaon_Cut2_data->Fill(pT);
      if(NcrossedRowsTPC / NFindableTPC > cutCrossedOverFindable)  // crossed over findable
        h_kaon_Cut3_data->Fill(pT);
      if(NclusterPIDTPC > cutClusterWithdEdxPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of number of clusters with signal 
        h_kaon_Cut4_data->Fill(pT);
      if(NclusterPIDTPC > 40)  // N clusters with TPC signal > 40
        h_kaon_Cut5_data->Fill(pT);
      if(NclusterPIDTPC > 60)  // N clusters with TPC signal > 60
        h_kaon_Cut6_data->Fill(pT);
    }

    // -- charged (all particles) ---
    if(( n_sigma_TOF_K_data > -nSigmaID && n_sigma_TOF_K_data < nSigmaID)
        || ( n_sigma_TOF_p_data > -nSigmaID && n_sigma_TOF_p_data < nSigmaID)
        || ( n_sigma_TOF_pi_data > -nSigmaID && n_sigma_TOF_pi_data < nSigmaID)) {
    //if( 
    //    ( n_sigma_TOF_pi_data > -nSigmaID && n_sigma_TOF_pi_data < nSigmaID)
    //  ) {
      if(i%10000==0) {
//        cout<<"\nTOF n sigma K = "<<n_sigma_TPC_K_data<<endl;
//        cout<<"TOF n sigma p = "<<n_sigma_TPC_p_data<<endl;
//        cout<<"TOF n sigma pi = "<<n_sigma_TPC_pi_data<<endl;
      }
      // default cut
      h_charged_all_data->Fill(pT);
      // additional cut variations
      if(NcrossedRowsTPC > cutPtDepCrossedRows)  // pt dependent crossed rows
        h_charged_Cut1_data->Fill(pT);
      if(NclusterTPC > cutClusterPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of crossed rows
        h_charged_Cut2_data->Fill(pT);
      if(NcrossedRowsTPC / NFindableTPC > cutCrossedOverFindable)  // crossed over findable
        h_charged_Cut3_data->Fill(pT);
      if(NclusterPIDTPC > cutClusterWithdEdxPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of number of clusters with signal 
        h_charged_Cut4_data->Fill(pT);
      if(NclusterPIDTPC > 40)  // N clusters with TPC signal > 40
        h_charged_Cut5_data->Fill(pT);
      if(NclusterPIDTPC > 60)  // N clusters with TPC signal > 60
        h_charged_Cut6_data->Fill(pT);
    }

  }

  
  //
  // MC
  //

  // loop over candidates 
  for(int i=0;i<t_mc->GetEntries();i++) {
    if(i%100000 == 0) cout<<"-- track "<<i<<endl;
    //if(isTest && i>100000) break;

    t_mc->GetEntry(i);
    // convert to better format
    double pT = double( pT_mc) / 1000.;
    double NclusterTPC = double ( NclusterTPC_mc) ;
    double NclusterPIDTPC = double ( NclusterPIDTPC_mc) ;
    double NcrossedRowsTPC = double ( NcrossedRowsTPC_mc) ;
    double NFindableTPC = double ( NFindableTPC_mc) ;

    if(pT<0.00001) continue;
    if(NFindableTPC < 1) continue;
    double cutPtDepCrossedRows = 120. - (5./pT);
    // is proton
    if((tag_mc & (kIsProtonFromL | kIsProtonFromTOF))
        && ((isNsigmaCut && TMath::Abs(n_sigma_TPC_p_mc) < 300) || !isNsigmaCut)) { 
      //cout<<"is proton pT = "<<pT<<"\t crossed rows = "<<NcrossedRowsTPC<<"\tfindable clusters = "<<NFindableTPC<<"\tclusters = "<<NclusterTPC<<endl;

      // default cut
      //if(NcrossedRowsTPC > 70. && NclusterTPC > 50. && NcrossedRowsTPC / NFindableTPC > 0.8) {
        h_proton_all_mc->Fill(pT);
        // additional cut variations
        if(NcrossedRowsTPC > cutPtDepCrossedRows)  // pt dependent crossed rows
          h_proton_Cut1_mc->Fill(pT);
        if(NclusterTPC > cutClusterPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of crossed rows
          h_proton_Cut2_mc->Fill(pT);
        if(NcrossedRowsTPC / NFindableTPC > cutCrossedOverFindable)  // crossed over findable
          h_proton_Cut3_mc->Fill(pT);
        if(NclusterPIDTPC  > cutClusterWithdEdxPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of number of clusters with signal 
          h_proton_Cut4_mc->Fill(pT);
        if(NclusterPIDTPC > 40)  // N clusters with TPC signal > 40
          h_proton_Cut5_mc->Fill(pT);
        if(NclusterPIDTPC > 60)  // N clusters with TPC signal > 60
          h_proton_Cut6_mc->Fill(pT);
    }
    // is pion
    if((tag_mc & (kIsPionFromK0s | kIsPionFromL | kIsPionFromTOF)) 
        && ((isNsigmaCut && TMath::Abs(n_sigma_TPC_pi_mc) < 300) || !isNsigmaCut)) { 
      //cout<<"is pion pT = "<<pT<<"\t crossed rows = "<<NcrossedRowsTPC<<"\tfindable clusters = "<<NFindableTPC<<"\tclusters = "<<NclusterTPC<<endl;

        h_pion_all_mc->Fill(pT);
        // additional cut variations
        if(NcrossedRowsTPC > cutPtDepCrossedRows)  // pt dependent crossed rows
          h_pion_Cut1_mc->Fill(pT);
        if(NclusterTPC > cutClusterPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of crossed rows
          h_pion_Cut2_mc->Fill(pT);
        if(NcrossedRowsTPC / NFindableTPC > cutCrossedOverFindable)  // crossed over findable
          h_pion_Cut3_mc->Fill(pT);
        if(NclusterPIDTPC  > cutClusterWithdEdxPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of number of clusters with signal 
          h_pion_Cut4_mc->Fill(pT);
        if(NclusterPIDTPC > 40)  // N clusters with TPC signal > 40
          h_pion_Cut5_mc->Fill(pT);
        if(NclusterPIDTPC > 60)  // N clusters with TPC signal > 60
          h_pion_Cut6_mc->Fill(pT);
    }
    // is kaon
    if(tag_mc & kIsKaonFromTOF) {
      //cout<<"is pion pT = "<<pT<<"\t crossed rows = "<<NcrossedRowsTPC<<"\tfindable clusters = "<<NFindableTPC<<"\tclusters = "<<NclusterTPC<<endl;

        h_kaon_all_mc->Fill(pT);
        // additional cut variations
        if(NcrossedRowsTPC > cutPtDepCrossedRows)  // pt dependent crossed rows
          h_kaon_Cut1_mc->Fill(pT);
        if(NclusterTPC > cutClusterPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of crossed rows
          h_kaon_Cut2_mc->Fill(pT);
        if(NcrossedRowsTPC / NFindableTPC > cutCrossedOverFindable)  // crossed over findable
          h_kaon_Cut3_mc->Fill(pT);
        if(NclusterPIDTPC  > cutClusterWithdEdxPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of number of clusters with signal 
          h_kaon_Cut4_mc->Fill(pT);
        if(NclusterPIDTPC > 40)  // N clusters with TPC signal > 40
          h_kaon_Cut5_mc->Fill(pT);
        if(NclusterPIDTPC > 60)  // N clusters with TPC signal > 60
          h_kaon_Cut6_mc->Fill(pT);
    }

    // -- charged (all particles) ---
    if(( n_sigma_TOF_K_mc > -nSigmaID && n_sigma_TOF_K_mc < nSigmaID)
        || ( n_sigma_TOF_p_mc > -nSigmaID && n_sigma_TOF_p_mc < nSigmaID)
        || ( n_sigma_TOF_pi_mc > -nSigmaID && n_sigma_TOF_pi_mc < nSigmaID)) {
    //if(
    //    ( n_sigma_TOF_pi_mc > -nSigmaID && n_sigma_TOF_pi_mc < nSigmaID)
    //  ) {

      if(i%10000==0) {
//        cout<<"\nTOF n sigma K = "<<n_sigma_TOF_K_mc<<endl;
//        cout<<"TOF n sigma p = "<<n_sigma_TOF_p_mc<<endl;
//        cout<<"TOF n sigma pi = "<<n_sigma_TOF_pi_mc<<endl;
      }

      h_charged_all_mc->Fill(pT);
      // additional cut variations
      if(NcrossedRowsTPC > cutPtDepCrossedRows)  // pt dependent crossed rows
        h_charged_Cut1_mc->Fill(pT);
      if(NclusterTPC > cutClusterPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of crossed rows
        h_charged_Cut2_mc->Fill(pT);
      if(NcrossedRowsTPC / NFindableTPC > cutCrossedOverFindable)  // crossed over findable
        h_charged_Cut3_mc->Fill(pT);
      if(NclusterPIDTPC  > cutClusterWithdEdxPercentOfCrossedRows * NcrossedRowsTPC)  // percetage of number of clusters with signal 
        h_charged_Cut4_mc->Fill(pT);
      if(NclusterPIDTPC > 40)  // N clusters with TPC signal > 40
        h_charged_Cut5_mc->Fill(pT);
      if(NclusterPIDTPC > 60)  // N clusters with TPC signal > 60
        h_charged_Cut6_mc->Fill(pT);
    }

  }

  //
  // Draw
  //

  TLatex info; info.SetNDC(); info.SetTextFont(43); info.SetTextSize(17);

  TCanvas *cproton_data = new TCanvas("cproton_data");
  h_proton_Cut1_data->Divide(h_proton_Cut1_data,h_proton_all_data,1.,1.,"B");
  h_proton_Cut1_data->SetLineColor(kRed);
  h_proton_Cut1_data->SetStats(0);
  h_proton_Cut1_data->GetYaxis()->SetTitle("TPC cut variation / TPC default cut");
  h_proton_Cut1_data->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_proton_Cut1_data->SetTitle("");
  h_proton_Cut1_data->GetYaxis()->SetRangeUser(0.8,1.1);
  h_proton_Cut1_data->Draw();
  h_proton_Cut2_data->Divide(h_proton_Cut2_data,h_proton_all_data,1.,1.,"B");
  h_proton_Cut2_data->SetLineColor(kBlue);
  //h_proton_Cut2_data->Draw("SAME");
  h_proton_Cut3_data->Divide(h_proton_Cut3_data,h_proton_all_data,1.,1.,"B");
  h_proton_Cut3_data->SetLineColor(kGreen+2);
  h_proton_Cut3_data->Draw("SAME");
  h_proton_Cut4_data->Divide(h_proton_Cut4_data,h_proton_all_data,1.,1.,"B");
  h_proton_Cut4_data->SetLineColor(kMagenta);
  h_proton_Cut4_data->Draw("SAME");
  h_proton_Cut5_data->Divide(h_proton_Cut5_data,h_proton_all_data,1.,1.,"B");
  h_proton_Cut5_data->SetLineColor(kCyan);
  h_proton_Cut5_data->Draw("SAME");
  h_proton_Cut6_data->Divide(h_proton_Cut6_data,h_proton_all_data,1.,1.,"B");
  h_proton_Cut6_data->SetLineColor(kOrange);
  h_proton_Cut6_data->Draw("SAME");

  TLegend *leg_data = new TLegend(0.45,0.88,0.87,0.7);
  leg_data->AddEntry(h_proton_Cut1_data,"TPC crossed rows > 120-5/p_{T}","l"); 
  //leg_data->AddEntry(h_proton_Cut2_data,"TPC clusters > 0.65 TPC crossed rows","l"); 
  leg_data->AddEntry(h_proton_Cut3_data,"cr. rows / findable clusters > 0.9","l"); 
  leg_data->AddEntry(h_proton_Cut4_data,"TPC clusters with dE/dx > 0.5 TPC crossed rows ","l"); 
  leg_data->AddEntry(h_proton_Cut5_data,"TPC clusters with dE/dx > 40","l"); 
  leg_data->AddEntry(h_proton_Cut6_data,"TPC clusters with dE/dx > 60","l"); 
  leg_data->SetLineColor(kWhite);
  leg_data->Draw("SAME");

  info.DrawLatex(0.25,0.84,Form("Data %s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"protons"); 

  TCanvas *cpion_data = new TCanvas("cpion_data");
  h_pion_Cut1_data->Divide(h_pion_Cut1_data,h_pion_all_data,1.,1.,"B");
  h_pion_Cut1_data->SetLineColor(kRed);
  h_pion_Cut1_data->SetStats(0);
  h_pion_Cut1_data->GetYaxis()->SetTitle("TPC cut variation / TPC default cut");
  h_pion_Cut1_data->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_pion_Cut1_data->SetTitle("");
  h_pion_Cut1_data->GetYaxis()->SetRangeUser(0.8,1.1);
  h_pion_Cut1_data->Draw();
  h_pion_Cut2_data->Divide(h_pion_Cut2_data,h_pion_all_data,1.,1.,"B");
  h_pion_Cut2_data->SetLineColor(kBlue);
  //h_pion_Cut2_data->Draw("SAME");
  h_pion_Cut3_data->Divide(h_pion_Cut3_data,h_pion_all_data,1.,1.,"B");
  h_pion_Cut3_data->SetLineColor(kGreen+2);
  h_pion_Cut3_data->Draw("SAME");
  h_pion_Cut4_data->Divide(h_pion_Cut4_data,h_pion_all_data,1.,1.,"B");
  h_pion_Cut4_data->SetLineColor(kMagenta);
  h_pion_Cut4_data->Draw("SAME");
  h_pion_Cut5_data->Divide(h_pion_Cut5_data,h_pion_all_data,1.,1.,"B");
  h_pion_Cut5_data->SetLineColor(kCyan);
  h_pion_Cut5_data->Draw("SAME");
  h_pion_Cut6_data->Divide(h_pion_Cut6_data,h_pion_all_data,1.,1.,"B");
  h_pion_Cut6_data->SetLineColor(kOrange);
  h_pion_Cut6_data->Draw("SAME");
  leg_data->Draw("SAME");
  info.DrawLatex(0.25,0.84,Form("Data %s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"pions"); 

  TCanvas *ckaon_data = new TCanvas("ckaon_data");
  h_kaon_Cut1_data->Divide(h_kaon_Cut1_data,h_kaon_all_data,1.,1.,"B");
  h_kaon_Cut1_data->SetLineColor(kRed);
  h_kaon_Cut1_data->SetStats(0);
  h_kaon_Cut1_data->GetYaxis()->SetTitle("TPC cut variation / TPC default cut");
  h_kaon_Cut1_data->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_kaon_Cut1_data->SetTitle("");
  h_kaon_Cut1_data->GetYaxis()->SetRangeUser(0.8,1.1);
  h_kaon_Cut1_data->Draw();
  h_kaon_Cut2_data->Divide(h_kaon_Cut2_data,h_kaon_all_data,1.,1.,"B");
  h_kaon_Cut2_data->SetLineColor(kBlue);
  //h_kaon_Cut2_data->Draw("SAME");
  h_kaon_Cut3_data->Divide(h_kaon_Cut3_data,h_kaon_all_data,1.,1.,"B");
  h_kaon_Cut3_data->SetLineColor(kGreen+2);
  h_kaon_Cut3_data->Draw("SAME");
  h_kaon_Cut4_data->Divide(h_kaon_Cut4_data,h_kaon_all_data,1.,1.,"B");
  h_kaon_Cut4_data->SetLineColor(kMagenta);
  h_kaon_Cut4_data->Draw("SAME");
  h_kaon_Cut5_data->Divide(h_kaon_Cut5_data,h_kaon_all_data,1.,1.,"B");
  h_kaon_Cut5_data->SetLineColor(kCyan);
  h_kaon_Cut5_data->Draw("SAME");
  h_kaon_Cut6_data->Divide(h_kaon_Cut6_data,h_kaon_all_data,1.,1.,"B");
  h_kaon_Cut6_data->SetLineColor(kOrange);
  h_kaon_Cut6_data->Draw("SAME");
  leg_data->Draw("SAME");
  info.DrawLatex(0.25,0.84,Form("Data %s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"kaons"); 

  TCanvas *ccharged_data = new TCanvas("ccharged_data");
  h_charged_Cut1_data->Divide(h_charged_Cut1_data,h_charged_all_data,1.,1.,"B");
  h_charged_Cut1_data->SetLineColor(kRed);
  h_charged_Cut1_data->SetStats(0);
  h_charged_Cut1_data->GetYaxis()->SetTitle("TPC cut variation / TPC default cut");
  h_charged_Cut1_data->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_charged_Cut1_data->SetTitle("");
  h_charged_Cut1_data->GetYaxis()->SetRangeUser(0.8,1.1);
  h_charged_Cut1_data->Draw();
  h_charged_Cut2_data->Divide(h_charged_Cut2_data,h_charged_all_data,1.,1.,"B");
  h_charged_Cut2_data->SetLineColor(kBlue);
  //h_charged_Cut2_data->Draw("SAME");
  h_charged_Cut3_data->Divide(h_charged_Cut3_data,h_charged_all_data,1.,1.,"B");
  h_charged_Cut3_data->SetLineColor(kGreen+2);
  h_charged_Cut3_data->Draw("SAME");
  h_charged_Cut4_data->Divide(h_charged_Cut4_data,h_charged_all_data,1.,1.,"B");
  h_charged_Cut4_data->SetLineColor(kMagenta);
  h_charged_Cut4_data->Draw("SAME");
  h_charged_Cut5_data->Divide(h_charged_Cut5_data,h_charged_all_data,1.,1.,"B");
  h_charged_Cut5_data->SetLineColor(kCyan);
  h_charged_Cut5_data->Draw("SAME");
  h_charged_Cut6_data->Divide(h_charged_Cut6_data,h_charged_all_data,1.,1.,"B");
  h_charged_Cut6_data->SetLineColor(kOrange);
  h_charged_Cut6_data->Draw("SAME");
  leg_data->Draw("SAME");
  info.DrawLatex(0.25,0.84,Form("Data %s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"charged particles"); 


  TCanvas *cproton_mc = new TCanvas("cproton_mc");
  h_proton_Cut1_mc->Divide(h_proton_Cut1_mc,h_proton_all_mc,1.,1.,"B");
  h_proton_Cut1_mc->SetLineColor(kRed);
  h_proton_Cut1_mc->SetStats(0);
  h_proton_Cut1_mc->GetYaxis()->SetTitle("TPC cut variation / TPC default cut");
  h_proton_Cut1_mc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_proton_Cut1_mc->SetTitle("");
  h_proton_Cut1_mc->GetYaxis()->SetRangeUser(0.8,1.1);
  h_proton_Cut1_mc->Draw();
  h_proton_Cut2_mc->Divide(h_proton_Cut2_mc,h_proton_all_mc,1.,1.,"B");
  h_proton_Cut2_mc->SetLineColor(kBlue);
  //h_proton_Cut2_mc->Draw("SAME");
  h_proton_Cut3_mc->Divide(h_proton_Cut3_mc,h_proton_all_mc,1.,1.,"B");
  h_proton_Cut3_mc->SetLineColor(kGreen+2);
  h_proton_Cut3_mc->Draw("SAME");
  h_proton_Cut4_mc->Divide(h_proton_Cut4_mc,h_proton_all_mc,1.,1.,"B");
  h_proton_Cut4_mc->SetLineColor(kMagenta);
  h_proton_Cut4_mc->Draw("SAME");
  h_proton_Cut5_mc->Divide(h_proton_Cut5_mc,h_proton_all_mc,1.,1.,"B");
  h_proton_Cut5_mc->SetLineColor(kCyan);
  h_proton_Cut5_mc->Draw("SAME");
  h_proton_Cut6_mc->Divide(h_proton_Cut6_mc,h_proton_all_mc,1.,1.,"B");
  h_proton_Cut6_mc->SetLineColor(kOrange);
  h_proton_Cut6_mc->Draw("SAME");

  info.DrawLatex(0.25,0.84,Form("MC %s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"protons"); 

  TLegend *leg_mc = new TLegend(0.45,0.88,0.87,0.7);
  leg_mc->AddEntry(h_proton_Cut1_mc,"TPC crossed rows > 120-5/p_{T}","l"); 
  //leg_mc->AddEntry(h_proton_Cut2_mc,"TPC clusters > 0.65 TPC crossed rows","l"); 
  leg_mc->AddEntry(h_proton_Cut3_mc,"cr. rows / findable clusters > 0.9","l"); 
  leg_mc->AddEntry(h_proton_Cut4_mc,"TPC clusters with dE/dx > 0.5 TPC crossed rows ","l"); 
  leg_mc->AddEntry(h_proton_Cut5_mc,"TPC clusters with dE/dx > 40","l"); 
  leg_mc->AddEntry(h_proton_Cut6_mc,"TPC clusters with dE/dx > 60","l"); 
  leg_mc->SetLineColor(kWhite);
  leg_mc->Draw("SAME");

  TCanvas *cpion_mc = new TCanvas("cpion_mc");
  h_pion_Cut1_mc->Divide(h_pion_Cut1_mc,h_pion_all_mc,1.,1.,"B");
  h_pion_Cut1_mc->SetLineColor(kRed);
  h_pion_Cut1_mc->SetStats(0);
  h_pion_Cut1_mc->GetYaxis()->SetTitle("TPC cut variation / TPC default cut");
  h_pion_Cut1_mc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_pion_Cut1_mc->SetTitle("");
  h_pion_Cut1_mc->GetYaxis()->SetRangeUser(0.8,1.1);
  h_pion_Cut1_mc->Draw();
  h_pion_Cut2_mc->Divide(h_pion_Cut2_mc,h_pion_all_mc,1.,1.,"B");
  h_pion_Cut2_mc->SetLineColor(kBlue);
  //h_pion_Cut2_mc->Draw("SAME");
  h_pion_Cut3_mc->Divide(h_pion_Cut3_mc,h_pion_all_mc,1.,1.,"B");
  h_pion_Cut3_mc->SetLineColor(kGreen+2);
  h_pion_Cut3_mc->Draw("SAME");
  h_pion_Cut4_mc->Divide(h_pion_Cut4_mc,h_pion_all_mc,1.,1.,"B");
  h_pion_Cut4_mc->SetLineColor(kMagenta);
  h_pion_Cut4_mc->Draw("SAME");
  h_pion_Cut5_mc->Divide(h_pion_Cut5_mc,h_pion_all_mc,1.,1.,"B");
  h_pion_Cut5_mc->SetLineColor(kCyan);
  h_pion_Cut5_mc->Draw("SAME");
  h_pion_Cut6_mc->Divide(h_pion_Cut6_mc,h_pion_all_mc,1.,1.,"B");
  h_pion_Cut6_mc->SetLineColor(kOrange);
  h_pion_Cut6_mc->Draw("SAME");
  leg_mc->Draw("SAME");

  info.DrawLatex(0.25,0.84,Form("MC %s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"pions"); 

  TCanvas *ckaon_mc = new TCanvas("ckaon_mc");
  h_kaon_Cut1_mc->Divide(h_kaon_Cut1_mc,h_kaon_all_mc,1.,1.,"B");
  h_kaon_Cut1_mc->SetLineColor(kRed);
  h_kaon_Cut1_mc->SetStats(0);
  h_kaon_Cut1_mc->GetYaxis()->SetTitle("TPC cut variation / TPC default cut");
  h_kaon_Cut1_mc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_kaon_Cut1_mc->SetTitle("");
  h_kaon_Cut1_mc->GetYaxis()->SetRangeUser(0.8,1.1);
  h_kaon_Cut1_mc->Draw();
  h_kaon_Cut2_mc->Divide(h_kaon_Cut2_mc,h_kaon_all_mc,1.,1.,"B");
  h_kaon_Cut2_mc->SetLineColor(kBlue);
  //h_kaon_Cut2_mc->Draw("SAME");
  h_kaon_Cut3_mc->Divide(h_kaon_Cut3_mc,h_kaon_all_mc,1.,1.,"B");
  h_kaon_Cut3_mc->SetLineColor(kGreen+2);
  h_kaon_Cut3_mc->Draw("SAME");
  h_kaon_Cut4_mc->Divide(h_kaon_Cut4_mc,h_kaon_all_mc,1.,1.,"B");
  h_kaon_Cut4_mc->SetLineColor(kMagenta);
  h_kaon_Cut4_mc->Draw("SAME");
  h_kaon_Cut5_mc->Divide(h_kaon_Cut5_mc,h_kaon_all_mc,1.,1.,"B");
  h_kaon_Cut5_mc->SetLineColor(kCyan);
  h_kaon_Cut5_mc->Draw("SAME");
  h_kaon_Cut6_mc->Divide(h_kaon_Cut6_mc,h_kaon_all_mc,1.,1.,"B");
  h_kaon_Cut6_mc->SetLineColor(kOrange);
  h_kaon_Cut6_mc->Draw("SAME");
  leg_mc->Draw("SAME");

  info.DrawLatex(0.25,0.84,Form("MC %s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"kaons"); 

  TCanvas *ccharged_mc = new TCanvas("ccharged_mc");
  h_charged_Cut1_mc->Divide(h_charged_Cut1_mc,h_charged_all_mc,1.,1.,"B");
  h_charged_Cut1_mc->SetLineColor(kRed);
  h_charged_Cut1_mc->SetStats(0);
  h_charged_Cut1_mc->GetYaxis()->SetTitle("TPC cut variation / TPC default cut");
  h_charged_Cut1_mc->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_charged_Cut1_mc->SetTitle("");
  h_charged_Cut1_mc->GetYaxis()->SetRangeUser(0.8,1.1);
  h_charged_Cut1_mc->Draw();
  h_charged_Cut2_mc->Divide(h_charged_Cut2_mc,h_charged_all_mc,1.,1.,"B");
  h_charged_Cut2_mc->SetLineColor(kBlue);
  //h_charged_Cut2_mc->Draw("SAME");
  h_charged_Cut3_mc->Divide(h_charged_Cut3_mc,h_charged_all_mc,1.,1.,"B");
  h_charged_Cut3_mc->SetLineColor(kGreen+2);
  h_charged_Cut3_mc->Draw("SAME");
  h_charged_Cut4_mc->Divide(h_charged_Cut4_mc,h_charged_all_mc,1.,1.,"B");
  h_charged_Cut4_mc->SetLineColor(kMagenta);
  h_charged_Cut4_mc->Draw("SAME");
  h_charged_Cut5_mc->Divide(h_charged_Cut5_mc,h_charged_all_mc,1.,1.,"B");
  h_charged_Cut5_mc->SetLineColor(kCyan);
  h_charged_Cut5_mc->Draw("SAME");
  h_charged_Cut6_mc->Divide(h_charged_Cut6_mc,h_charged_all_mc,1.,1.,"B");
  h_charged_Cut6_mc->SetLineColor(kOrange);
  h_charged_Cut6_mc->Draw("SAME");
  leg_mc->Draw("SAME");

  info.DrawLatex(0.25,0.84,Form("MC %s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"charged particles"); 


  TString post = "";
  if(cent==kCentral) post += "_cent";
  if(cent==kSemiCentral) post += "_semicent";
  if(dataset==k18r) post += "_18r";
  if(dataset==k18q) post += "_18q";
  if(dataset==kBoth) post += "_18qr";
  if(dataset==k17pq) post += "_17pq";
  if(dataset==k16qt) post += "_17qt";
  if(speciescuts==kPionCuts) post += "_pionCuts";
  if(speciescuts==kProtonCuts) post += "_protonCuts";
  if(isNsigmaCut) post += "_3sigmaTPC";

  cproton_data->SaveAs(Form("ProtonCutVar_Data%s.png",post.Data()));
  cpion_data->SaveAs(Form("PionCutVar_Data%s.png",post.Data()));
  ckaon_data->SaveAs(Form("KaonCutVar_Data%s.png",post.Data()));
  ccharged_data->SaveAs(Form("ChargedParticleCutVar_Data%s.png",post.Data()));
  cproton_mc->SaveAs(Form("ProtonCutVar_MC%s.png",post.Data()));
  cpion_mc->SaveAs(Form("PionCutVar_MC%s.png",post.Data()));
  ckaon_mc->SaveAs(Form("KaonCutVar_MC%s.png",post.Data()));
  ccharged_mc->SaveAs(Form("ChargedParticleCutVar_MC%s.png",post.Data()));

  // relative variation for systematics
  //

  TH1D *h_proton_Cut1_data_over_mc = (TH1D*)h_proton_Cut1_data->Clone("h_proton_Cut1_data_over_mc");
  TH1D *h_proton_Cut2_data_over_mc = (TH1D*)h_proton_Cut2_data->Clone("h_proton_Cut2_data_over_mc");
  TH1D *h_proton_Cut3_data_over_mc = (TH1D*)h_proton_Cut3_data->Clone("h_proton_Cut3_data_over_mc");
  TH1D *h_proton_Cut4_data_over_mc = (TH1D*)h_proton_Cut4_data->Clone("h_proton_Cut4_data_over_mc");
  TH1D *h_proton_Cut5_data_over_mc = (TH1D*)h_proton_Cut5_data->Clone("h_proton_Cut5_data_over_mc");
  TH1D *h_proton_Cut6_data_over_mc = (TH1D*)h_proton_Cut6_data->Clone("h_proton_Cut6_data_over_mc");

  TLegend *leg_ratio = new TLegend(0.2,0.4,0.55,0.14);
  leg_ratio->AddEntry(h_proton_Cut1_data_over_mc,"TPC crossed rows > 120-5/p_{T}","l"); 
  //leg_ratio->AddEntry(h_proton_Cut2_data_over_mc,"TPC clusters > 0.65 TPC crossed rows","l"); 
  leg_ratio->AddEntry(h_proton_Cut3_data_over_mc,"cr. rows / findable clusters > 0.9","l"); 
  leg_ratio->AddEntry(h_proton_Cut4_data_over_mc,"TPC clusters with dE/dx > 0.5 TPC crossed rows ","l"); 
  leg_ratio->AddEntry(h_proton_Cut5_data_over_mc,"TPC clusters with dE/dx > 40","l"); 
  leg_ratio->AddEntry(h_proton_Cut6_data_over_mc,"TPC clusters with dE/dx > 60","l"); 
  leg_ratio->SetLineColor(kWhite);

  h_proton_Cut1_data_over_mc->Divide(h_proton_Cut1_mc);
  h_proton_Cut2_data_over_mc->Divide(h_proton_Cut2_mc);
  h_proton_Cut3_data_over_mc->Divide(h_proton_Cut3_mc);
  h_proton_Cut4_data_over_mc->Divide(h_proton_Cut4_mc);
  h_proton_Cut5_data_over_mc->Divide(h_proton_Cut5_mc);
  h_proton_Cut6_data_over_mc->Divide(h_proton_Cut6_mc);

  TCanvas *c5 = new TCanvas("cproton_data_over_mc");
  h_proton_Cut1_data_over_mc->GetYaxis()->SetTitle("Data / MC");
  h_proton_Cut1_data_over_mc->Draw();
  //h_proton_Cut2_data_over_mc->Draw("SAME");
  h_proton_Cut3_data_over_mc->Draw("SAME");
  h_proton_Cut4_data_over_mc->Draw("SAME");
  h_proton_Cut5_data_over_mc->Draw("SAME");
  h_proton_Cut6_data_over_mc->Draw("SAME");
  leg_ratio->Draw("SAME");

  info.DrawLatex(0.25,0.84,Form("%s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"protons"); 

  TH1D *h_pion_Cut1_data_over_mc = (TH1D*)h_pion_Cut1_data->Clone("h_pion_Cut1_data_over_mc");
  TH1D *h_pion_Cut2_data_over_mc = (TH1D*)h_pion_Cut2_data->Clone("h_pion_Cut2_data_over_mc");
  TH1D *h_pion_Cut3_data_over_mc = (TH1D*)h_pion_Cut3_data->Clone("h_pion_Cut3_data_over_mc");
  TH1D *h_pion_Cut4_data_over_mc = (TH1D*)h_pion_Cut4_data->Clone("h_pion_Cut4_data_over_mc");
  TH1D *h_pion_Cut5_data_over_mc = (TH1D*)h_pion_Cut5_data->Clone("h_pion_Cut5_data_over_mc");
  TH1D *h_pion_Cut6_data_over_mc = (TH1D*)h_pion_Cut6_data->Clone("h_pion_Cut6_data_over_mc");

  h_pion_Cut1_data_over_mc->Divide(h_pion_Cut1_mc);
  h_pion_Cut2_data_over_mc->Divide(h_pion_Cut2_mc);
  h_pion_Cut3_data_over_mc->Divide(h_pion_Cut3_mc);
  h_pion_Cut4_data_over_mc->Divide(h_pion_Cut4_mc);
  h_pion_Cut5_data_over_mc->Divide(h_pion_Cut5_mc);
  h_pion_Cut6_data_over_mc->Divide(h_pion_Cut6_mc);

  TCanvas *c6 = new TCanvas("cpion_data_over_mc");
  h_pion_Cut1_data_over_mc->GetYaxis()->SetTitle("Data / MC");
  h_pion_Cut1_data_over_mc->Draw();
  //h_pion_Cut2_data_over_mc->Draw("SAME");
  h_pion_Cut3_data_over_mc->Draw("SAME");
  h_pion_Cut4_data_over_mc->Draw("SAME");
  h_pion_Cut5_data_over_mc->Draw("SAME");
  h_pion_Cut6_data_over_mc->Draw("SAME");
  leg_ratio->Draw("SAME");
  info.DrawLatex(0.25,0.84,Form("%s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"pions"); 

  TH1D *h_kaon_Cut1_data_over_mc = (TH1D*)h_kaon_Cut1_data->Clone("h_kaon_Cut1_data_over_mc");
  TH1D *h_kaon_Cut2_data_over_mc = (TH1D*)h_kaon_Cut2_data->Clone("h_kaon_Cut2_data_over_mc");
  TH1D *h_kaon_Cut3_data_over_mc = (TH1D*)h_kaon_Cut3_data->Clone("h_kaon_Cut3_data_over_mc");
  TH1D *h_kaon_Cut4_data_over_mc = (TH1D*)h_kaon_Cut4_data->Clone("h_kaon_Cut4_data_over_mc");
  TH1D *h_kaon_Cut5_data_over_mc = (TH1D*)h_kaon_Cut5_data->Clone("h_kaon_Cut5_data_over_mc");
  TH1D *h_kaon_Cut6_data_over_mc = (TH1D*)h_kaon_Cut6_data->Clone("h_kaon_Cut6_data_over_mc");

  h_kaon_Cut1_data_over_mc->Divide(h_kaon_Cut1_mc);
  h_kaon_Cut2_data_over_mc->Divide(h_kaon_Cut2_mc);
  h_kaon_Cut3_data_over_mc->Divide(h_kaon_Cut3_mc);
  h_kaon_Cut4_data_over_mc->Divide(h_kaon_Cut4_mc);
  h_kaon_Cut5_data_over_mc->Divide(h_kaon_Cut5_mc);
  h_kaon_Cut6_data_over_mc->Divide(h_kaon_Cut6_mc);

//  TF1 *fit_Cut1 = new TF1("fit_Cut1","[0]",2.5,8);
//  fit_Cut1->SetLineColor(kRed);
//  h_kaon_Cut1_data_over_mc->Fit("fit_Cut1","R");
//  TF1 *fit_Cut2 = new TF1("fit_Cut2","[0]",2.5,8);
//  fit_Cut2->SetLineColor(kBlue);
//  h_kaon_Cut2_data_over_mc->Fit("fit_Cut2","R");
//  TF1 *fit_Cut3 = new TF1("fit_Cut3","[0]",2.5,8);
//  fit_Cut3->SetLineColor(kGreen+2);
//  h_kaon_Cut3_data_over_mc->Fit("fit_Cut3","R");
//  TF1 *fit_Cut4 = new TF1("fit_Cut4","[0]",2.5,8);
//  fit_Cut4->SetLineColor(kMagenta);
//  h_kaon_Cut4_data_over_mc->Fit("fit_Cut4","R");
//  TF1 *fit_Cut5 = new TF1("fit_Cut5","[0]",2.5,8);
//  fit_Cut5->SetLineColor(kCyan);
//  h_kaon_Cut5_data_over_mc->Fit("fit_Cut5","R");
//  TF1 *fit_Cut6 = new TF1("fit_Cut6","[0]",2.5,8);
//  fit_Cut6->SetLineColor(kOrange);
//  h_kaon_Cut6_data_over_mc->Fit("fit_Cut6","R");


  TCanvas *c7 = new TCanvas("ckaon_data_over_mc");
  h_kaon_Cut1_data_over_mc->GetYaxis()->SetTitle("Data / MC");
  h_kaon_Cut1_data_over_mc->Draw();
  //h_kaon_Cut2_data_over_mc->Draw("SAME");
  h_kaon_Cut3_data_over_mc->Draw("SAME");
  h_kaon_Cut4_data_over_mc->Draw("SAME");
  h_kaon_Cut5_data_over_mc->Draw("SAME");
  h_kaon_Cut6_data_over_mc->Draw("SAME");
  leg_ratio->Draw("SAME");
  info.DrawLatex(0.25,0.84,Form("%s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"kaons"); 

  TH1D *h_charged_Cut1_data_over_mc = (TH1D*)h_charged_Cut1_data->Clone("h_charged_Cut1_data_over_mc");
  TH1D *h_charged_Cut2_data_over_mc = (TH1D*)h_charged_Cut2_data->Clone("h_charged_Cut2_data_over_mc");
  TH1D *h_charged_Cut3_data_over_mc = (TH1D*)h_charged_Cut3_data->Clone("h_charged_Cut3_data_over_mc");
  TH1D *h_charged_Cut4_data_over_mc = (TH1D*)h_charged_Cut4_data->Clone("h_charged_Cut4_data_over_mc");
  TH1D *h_charged_Cut5_data_over_mc = (TH1D*)h_charged_Cut5_data->Clone("h_charged_Cut5_data_over_mc");
  TH1D *h_charged_Cut6_data_over_mc = (TH1D*)h_charged_Cut6_data->Clone("h_charged_Cut6_data_over_mc");

  h_charged_Cut1_data_over_mc->Divide(h_charged_Cut1_mc);
  h_charged_Cut2_data_over_mc->Divide(h_charged_Cut2_mc);
  h_charged_Cut3_data_over_mc->Divide(h_charged_Cut3_mc);
  h_charged_Cut4_data_over_mc->Divide(h_charged_Cut4_mc);
  h_charged_Cut5_data_over_mc->Divide(h_charged_Cut5_mc);
  h_charged_Cut6_data_over_mc->Divide(h_charged_Cut6_mc);

  TCanvas *c8 = new TCanvas("ccharged_data_over_mc");
  h_charged_Cut1_data_over_mc->GetYaxis()->SetTitle("Data / MC");
  h_charged_Cut1_data_over_mc->Draw();
  //h_charged_Cut2_data_over_mc->Draw("SAME");
  h_charged_Cut3_data_over_mc->Draw("SAME");
  h_charged_Cut4_data_over_mc->Draw("SAME");
  h_charged_Cut5_data_over_mc->Draw("SAME");
  h_charged_Cut6_data_over_mc->Draw("SAME");
  leg_ratio->Draw("SAME");
  info.DrawLatex(0.25,0.84,Form("%s",datasetcent.Data())); 
  info.DrawLatex(0.25,0.80,"charged particles"); 



  cout<<"-- pion sigle track uncertainty--"<<endl;
  for(int i=0;i<h_pion_Cut1_data_over_mc->GetNbinsX();i++) {
    Double_t diff1 = TMath::Abs(h_pion_Cut1_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff2 = 0;
    Double_t diff3 = TMath::Abs(h_pion_Cut3_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff4 = TMath::Abs(h_pion_Cut4_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff5 = TMath::Abs(h_pion_Cut5_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff6 = TMath::Abs(h_pion_Cut6_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diffAll[6] = {diff1, diff2, diff3, diff4,diff5,diff6};
    Double_t syst = TMath::MaxElement(6,diffAll);
    cout<<setprecision(2);
    cout<<h_pion_all_data->GetBinLowEdge(i+1)<<" - "<<h_pion_all_data->GetBinLowEdge(i+1)+h_pion_all_data->GetBinWidth(i+1);
    cout<<"\t "<<syst<<endl;
    hunc_pion->SetBinContent(i+1,syst);
  }
  cout<<"-- kaon sigle track uncertainty--"<<endl;
  for(int i=0;i<h_kaon_Cut1_data_over_mc->GetNbinsX();i++) {
    Double_t diff1 = TMath::Abs(h_kaon_Cut1_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff2 = 0;
    Double_t diff3 = TMath::Abs(h_kaon_Cut3_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff4 = TMath::Abs(h_kaon_Cut4_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff5 = TMath::Abs(h_kaon_Cut5_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff6 = TMath::Abs(h_kaon_Cut6_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diffAll[6] = {diff1, diff2, diff3, diff4,diff5,diff6};
    Double_t syst = TMath::MaxElement(6,diffAll);
    cout<<setprecision(2);
    cout<<h_kaon_all_data->GetBinLowEdge(i+1)<<" - "<<h_kaon_all_data->GetBinLowEdge(i+1)+h_kaon_all_data->GetBinWidth(i+1);
    cout<<"\t "<<syst<<endl;
    hunc_kaon->SetBinContent(i+1,syst);
  }

  cout<<"-- proton sigle track uncertainty--"<<endl;
  for(int i=0;i<h_proton_Cut1_data_over_mc->GetNbinsX();i++) {
    Double_t diff1 = TMath::Abs(h_proton_Cut1_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff2 = 0;
    Double_t diff3 = TMath::Abs(h_proton_Cut3_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff4 = TMath::Abs(h_proton_Cut4_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff5 = TMath::Abs(h_proton_Cut5_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff6 = TMath::Abs(h_proton_Cut6_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diffAll[6] = {diff1, diff2, diff3, diff4, diff5, diff6};
    Double_t syst = TMath::MaxElement(6,diffAll);
    cout<<setprecision(2);
    cout<<h_proton_all_data->GetBinLowEdge(i+1)<<" - "<<h_proton_all_data->GetBinLowEdge(i+1)+h_proton_all_data->GetBinWidth(i+1);
    cout<<"\t "<<syst<<endl;

    hunc_proton->SetBinContent(i+1,syst);
  }
  cout<<"-- charged particle single track uncertainty--"<<endl;
  for(int i=0;i<h_charged_Cut1_data_over_mc->GetNbinsX();i++) {
    Double_t diff1 = TMath::Abs(h_charged_Cut1_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff2 = 0;
    Double_t diff3 = TMath::Abs(h_charged_Cut3_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff4 = TMath::Abs(h_charged_Cut4_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff5 = TMath::Abs(h_charged_Cut5_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diff6 = TMath::Abs(h_charged_Cut6_data_over_mc->GetBinContent(i+1) - 1) * 100;
    Double_t diffAll[6] = {diff1, diff2, diff3, diff4,diff5,diff6};
    Double_t syst = TMath::MaxElement(6,diffAll);
    cout<<setprecision(2);
    cout<<h_charged_all_data->GetBinLowEdge(i+1)<<" - "<<h_charged_all_data->GetBinLowEdge(i+1)+h_charged_all_data->GetBinWidth(i+1);
    cout<<"\t "<<syst<<endl;
    hunc_charged->SetBinContent(i+1,syst);
  }




  TFile *fout = new TFile(Form("TrackingTPCCutUnc%s.root",post.Data()),"RECREATE");

  h_pion_all_data->Write();
  h_pion_Cut1_data->Write();
  h_pion_Cut2_data->Write();
  h_pion_Cut3_data->Write();
  h_pion_Cut4_data->Write();
  h_pion_Cut5_data->Write();
  h_pion_Cut6_data->Write();
  h_pion_all_mc->Write();
  h_pion_Cut1_mc->Write();
  h_pion_Cut2_mc->Write();
  h_pion_Cut3_mc->Write();
  h_pion_Cut4_mc->Write();
  h_pion_Cut5_mc->Write();
  h_pion_Cut6_mc->Write();
  h_pion_Cut1_data_over_mc->Write(); 
  h_pion_Cut2_data_over_mc->Write(); 
  h_pion_Cut3_data_over_mc->Write(); 
  h_pion_Cut4_data_over_mc->Write(); 
  h_pion_Cut5_data_over_mc->Write(); 
  h_pion_Cut6_data_over_mc->Write(); 

  h_proton_all_data->Write();
  h_proton_Cut1_data->Write();
  h_proton_Cut2_data->Write();
  h_proton_Cut3_data->Write();
  h_proton_Cut4_data->Write();
  h_proton_Cut5_data->Write();
  h_proton_Cut6_data->Write();
  h_proton_all_mc->Write();
  h_proton_Cut1_mc->Write();
  h_proton_Cut2_mc->Write();
  h_proton_Cut3_mc->Write();
  h_proton_Cut4_mc->Write();
  h_proton_Cut5_mc->Write();
  h_proton_Cut6_mc->Write();
  h_proton_Cut1_data_over_mc->Write(); 
  h_proton_Cut2_data_over_mc->Write(); 
  h_proton_Cut3_data_over_mc->Write(); 
  h_proton_Cut4_data_over_mc->Write(); 
  h_proton_Cut5_data_over_mc->Write(); 
  h_proton_Cut6_data_over_mc->Write(); 

  h_kaon_all_data->Write();
  h_kaon_Cut1_data->Write();
  h_kaon_Cut2_data->Write();
  h_kaon_Cut3_data->Write();
  h_kaon_Cut4_data->Write();
  h_kaon_Cut5_data->Write();
  h_kaon_Cut6_data->Write();
  h_kaon_all_mc->Write();
  h_kaon_Cut1_mc->Write();
  h_kaon_Cut2_mc->Write();
  h_kaon_Cut3_mc->Write();
  h_kaon_Cut4_mc->Write();
  h_kaon_Cut5_mc->Write();
  h_kaon_Cut6_mc->Write();
  h_kaon_Cut1_data_over_mc->Write(); 
  h_kaon_Cut2_data_over_mc->Write(); 
  h_kaon_Cut3_data_over_mc->Write(); 
  h_kaon_Cut4_data_over_mc->Write(); 
  h_kaon_Cut5_data_over_mc->Write(); 
  h_kaon_Cut6_data_over_mc->Write(); 

  h_charged_all_data->Write();
  h_charged_Cut1_data->Write();
  h_charged_Cut2_data->Write();
  h_charged_Cut3_data->Write();
  h_charged_Cut4_data->Write();
  h_charged_Cut5_data->Write();
  h_charged_Cut6_data->Write();
  h_charged_all_mc->Write();
  h_charged_Cut1_mc->Write();
  h_charged_Cut2_mc->Write();
  h_charged_Cut3_mc->Write();
  h_charged_Cut4_mc->Write();
  h_charged_Cut5_mc->Write();
  h_charged_Cut6_mc->Write();
  h_charged_Cut1_data_over_mc->Write(); 
  h_charged_Cut2_data_over_mc->Write(); 
  h_charged_Cut3_data_over_mc->Write(); 
  h_charged_Cut4_data_over_mc->Write(); 
  h_charged_Cut5_data_over_mc->Write(); 
  h_charged_Cut6_data_over_mc->Write(); 


  hunc_pion->Write();
  hunc_proton->Write();
  hunc_kaon->Write();
  hunc_charged->Write();
  cproton_data->Write();
  cpion_data->Write();
  ckaon_data->Write();
  ccharged_data->Write();
  cproton_mc->Write();
  cpion_mc->Write();
  ckaon_mc->Write();
  ccharged_mc->Write();
  c5->Write();
  c6->Write();
  c7->Write();
  c8->Write();


  c5->SaveAs(Form("ProtonCutVar_DataOverMC%s.png",post.Data()));
  c6->SaveAs(Form("PionCutVar_DataOverMC%s.png",post.Data()));
  c7->SaveAs(Form("KaonCutVar_DataOverMC%s.png",post.Data()));
  c8->SaveAs(Form("ChargedParticleCutVar_DataOverMC%s.png",post.Data()));

    }
