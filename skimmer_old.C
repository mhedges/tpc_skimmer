#define skimmer_cxx
#include "skimmer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

void skimmer::Loop(TString FileName, TString OutputName)
//void skimmer(TString FileName, TString OutputName)
{
//   In a ROOT session, you can do:
//      Root > .L skimmer.C
//      Root > skimmer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   //if (fChain == 0) return;
   // Load and read BASF2 compatible TPC ROOT file (produced by Igal on KEKCC)
   //TFile df("/Volumes/KEK_USB/beast/phase1/data/TPC/tpc3_th50_data_1463360400.root","READ");
   TFile df(FileName,"READ");
   // Get the TTree
   TTree *dtr = (TTree*)df.Get("tree");

   class skimmer* event; 

   //TFile *ofile = new TFile("tpc_data.root", "RECREATE");

   TFile *ofile = new TFile(OutputName, "RECREATE");
   TTree *tr = new TTree("tr","TPC Event Data");

   // Initialize the TTree
   event->Init(dtr);

   // Activate all TBranches
   dtr->SetBranchStatus("*",1);

   // Make new TTree with relevant data
   
   int npoints, time_range, tot_sum; 
   int getentry;

   //int nrows = 336;
   //int ncol = 80;
   //int kMaxHits= nrows * ncol;

   int row[26880];
   int col[26880];
   int tot[26880];
   int bcid[26880];

   Float_t par_fit[5];
   Float_t par_fit_err[5];
   Float_t chi2, t_length, theta, phi;

   Double_t tstamp;
   
   tr->Branch("npoints",&npoints,"npoints/I");
   tr->Branch("row",&row,"row[npoints]/I");
   tr->Branch("col",&col,"col[npoints]/I");
   tr->Branch("bcid",&bcid,"bcid[npoints]/I");
   tr->Branch("tot",&tot,"tot[npoints]/I");
   tr->Branch("tstamp",&tstamp,"tstamp/D");
   tr->Branch("tot_sum",&tot_sum,"tot_sum/I");
   tr->Branch("time_range",&time_range,"time_range/I");
   tr->Branch("chi2",&chi2,"chi2/F");
   tr->Branch("t_length",&t_length,"t_length/F");
   tr->Branch("theta",&theta,"theta/F");
   tr->Branch("phi",&phi,"phi/F");
   tr->Branch("par_fit",&par_fit,"par_fit[5]/F");
   tr->Branch("par_fit_err",&par_fit_err,"par_fit_err[5]/F");

   //Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nbytes = 0, nb = 0;
   
   int nentries = dtr->GetEntriesFast();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      //if (jentry %1000 == 0)cout << "Event Counter: " << jentry << endl;

      //Long64_t ientry = LoadTree(jentry);

      //if (ientry < 0) break;
      //nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	  getentry = dtr->GetEntry(jentry);
	  npoints = event->MicrotpcMetaHits_m_pixNb[0];
	  tstamp = event->MicrotpcMetaHits_m_ts_start[0][0];
	  tot_sum = event->MicrotpcRecoTracks_m_totsum[0];
	  time_range = event->MicrotpcRecoTracks_m_time_range[0];
	  chi2 = event->MicrotpcRecoTracks_m_chi2[0];
	  t_length = event->MicrotpcRecoTracks_m_trl[0];
	  theta = event->MicrotpcRecoTracks_m_theta[0];
	  phi = event->MicrotpcRecoTracks_m_phi[0];
	  
	  //par_fit=event->MicrotpcRecoTracks_m_parFit[jentry];
	  //par_fit_err=event->MicrotpcRecoTracks_m_parFit_err[jentry];

	  for (int npars=0; npars<6;npars++){
	    par_fit[npars]=event->MicrotpcRecoTracks_m_parFit[0][npars];
	    par_fit_err[npars]=event->MicrotpcRecoTracks_m_parFit_err[0][npars];
	  }

	  
	  //for (int ii=0; ii<2; ii++){
	  //tstamp = event->MicrotpcMetaHits_m_ts_start[ii][0];
	  //if (tstamp > 0){
	  //  std::cout << "Timestamp entry is " << tstamp << " in event " << jentry << " Element = " << ii << endl;
	  //  std::cin.ignore();
	  //  std::cin.get();
	  //}
	  //}
	  
	  for (int pixn=0; pixn<npoints;pixn++){
		row[pixn]=event->MicrotpcDataHits_m_row[pixn];
		col[pixn]=event->MicrotpcDataHits_m_column[pixn];
		bcid[pixn]=event->MicrotpcDataHits_m_BCID[pixn];
		tot[pixn]=event->MicrotpcDataHits_m_TOT[pixn];
		//std::cout << "At pixel number " << pixn << std::endl;
		//std::cout << "Column number is " << col[pixn] << std::endl;
	  }
	  tr->Fill();


	  //if (jentry == 2) {
	  //  cout << "Currently at event " << jentry << "\n";
	  //  for (int i=0;i<MicrotpcDataHits_; i++){
	  //    evt->Fill(MicrotpcDataHits_m_column[i],MicrotpcDataHits_m_row[i],MicrotpcDataHits_m_TOT[i]);
	  //  }
	  //}
   }
   //tr->Write();
   ofile->Write();
   ofile->Close();

/*
//create a Tree file tree4.root
  TFile f("tree4.root","RECREATE");

  // Create a ROOT Tree
  TTree t4("t4","A Tree with Events");

  // Create a pointer to an Event object
  Event *event = new Event();

  // Create two branches, split one.
  t4.Branch("event_split", &event,16000,99);
  t4.Branch("event_not_split", &event,16000,0);
  */


//   if (fChain == 0) return;
//
//   Long64_t nentries = fChain->GetEntriesFast();
//
//   Long64_t nbytes = 0, nb = 0;
//   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//      Long64_t ientry = LoadTree(jentry);
//      if (ientry < 0) break;
//      nb = fChain->GetEntry(jentry);   nbytes += nb;
//      // if (Cut(ientry) < 0) continue;
//   }
}

//int main(int argc, char * argv[])
int main()
//int skimmer(int argc, char * argv[])
{
  //cout << "First argv is: " << argv[0] << endl;
  //cout << "Analyzing data from: " << argv[1] << endl;
  skimmer ana;
  ana.Loop("a", "b");
  //ana.Loop(argv[2]);
  //ana.Loop(argv[2], argv[3]);
  //skimmer::Loop(char* argv[0], char* argv[1]);
  //skimmer((char*) argv[0], (char)* argv[1]);
  return 1;
}
