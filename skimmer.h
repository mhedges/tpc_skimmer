//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 11 08:52:38 2016 by ROOT version 5.34/34
// from TTree tree/tree
// found on file: tpc_data_1460214023.root
//////////////////////////////////////////////////////////

#ifndef skimmer_h
#define skimmer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>
#include <TClonesArray.h>

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxMicrotpcDataHits = 30000;
   const Int_t kMaxMicrotpcMetaHits = 1;
   const Int_t kMaxMicrotpcRecoTracks = 1;
   const Int_t kMaxm_elements = 1;

class skimmer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //Belle2::EventMetaData *EventMetaData;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UInt_t          m_event;
   Int_t           m_run;
   Int_t           m_subrun;
   Int_t           m_experiment;
   Int_t           m_production;
   ULong64_t       m_time;
   string          m_parentLfn;
   Double_t        m_generatedWeight;
   UInt_t          m_errorFlag;
   Int_t           MicrotpcDataHits_;
   UInt_t          MicrotpcDataHits_fUniqueID[kMaxMicrotpcDataHits];   //[MicrotpcDataHits_]
   UInt_t          MicrotpcDataHits_fBits[kMaxMicrotpcDataHits];   //[MicrotpcDataHits_]
   UShort_t        MicrotpcDataHits_m_backgroundTag[kMaxMicrotpcDataHits];   //[MicrotpcDataHits_]
   Int_t           MicrotpcDataHits_m_column[kMaxMicrotpcDataHits];   //[MicrotpcDataHits_]
   Int_t           MicrotpcDataHits_m_row[kMaxMicrotpcDataHits];   //[MicrotpcDataHits_]
   Int_t           MicrotpcDataHits_m_BCID[kMaxMicrotpcDataHits];   //[MicrotpcDataHits_]
   Int_t           MicrotpcDataHits_m_TOT[kMaxMicrotpcDataHits];   //[MicrotpcDataHits_]
   Int_t           MicrotpcDataHits_m_detNb[kMaxMicrotpcDataHits];   //[MicrotpcDataHits_]
   Int_t           MicrotpcMetaHits_;
   UInt_t          MicrotpcMetaHits_fUniqueID[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   UInt_t          MicrotpcMetaHits_fBits[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   UShort_t        MicrotpcMetaHits_m_backgroundTag[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Int_t           MicrotpcMetaHits_m_detNb[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Int_t           MicrotpcMetaHits_m_pixNb[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Int_t           MicrotpcMetaHits_m_ts_nb[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_ts_start[kMaxMicrotpcMetaHits][10];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_ts_stop[kMaxMicrotpcMetaHits][10];   //[MicrotpcMetaHits_]
   Float_t         MicrotpcMetaHits_m_Temperature[kMaxMicrotpcMetaHits][4];   //[MicrotpcMetaHits_]
   Float_t         MicrotpcMetaHits_m_Pressure[kMaxMicrotpcMetaHits][2];   //[MicrotpcMetaHits_]
   Float_t         MicrotpcMetaHits_m_Flow[kMaxMicrotpcMetaHits][2];   //[MicrotpcMetaHits_]
   Float_t         MicrotpcMetaHits_m_SetFlow[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Float_t         MicrotpcMetaHits_m_GetFlow[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_IHER[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_PHER[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_tHER[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_flagHER[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_ILER[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_PLER[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_tLER[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Double_t        MicrotpcMetaHits_m_flagLER[kMaxMicrotpcMetaHits];   //[MicrotpcMetaHits_]
   Int_t           MicrotpcRecoTracks_;
   UInt_t          MicrotpcRecoTracks_fUniqueID[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   UInt_t          MicrotpcRecoTracks_fBits[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   UShort_t        MicrotpcRecoTracks_m_backgroundTag[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Int_t           MicrotpcRecoTracks_m_detNb[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Int_t           MicrotpcRecoTracks_m_pixnb[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_chi2[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_theta[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_phi[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_esum[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Int_t           MicrotpcRecoTracks_m_totsum[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_trl[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Int_t           MicrotpcRecoTracks_m_time_range[kMaxMicrotpcRecoTracks];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_parFit[kMaxMicrotpcRecoTracks][5];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_parFit_err[kMaxMicrotpcRecoTracks][5];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_cov[kMaxMicrotpcRecoTracks][25];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_impact_x[kMaxMicrotpcRecoTracks][4];   //[MicrotpcRecoTracks_]
   Float_t         MicrotpcRecoTracks_m_impact_y[kMaxMicrotpcRecoTracks][4];   //[MicrotpcRecoTracks_]
   Int_t           MicrotpcRecoTracks_m_side[kMaxMicrotpcRecoTracks][16];   //[MicrotpcRecoTracks_]
   Int_t           MicrotpcRecoTracks_m_partID[kMaxMicrotpcRecoTracks][6];   //[MicrotpcRecoTracks_]
 //Belle2::RelationContainer *MicrotpcRecoTracksToMicrotpcDataHits;
   //UInt_t          fUniqueID;
   //UInt_t          fBits;
   Int_t           m_elements_;
   UInt_t          m_elements_fUniqueID[kMaxm_elements];   //[m_elements_]
   UInt_t          m_elements_fBits[kMaxm_elements];   //[m_elements_]
   UInt_t          m_elements_m_from[kMaxm_elements];   //[m_elements_]
   vector<unsigned int> m_elements_m_to[kMaxm_elements];
   vector<float>   m_elements_m_weights[kMaxm_elements];
   string          m_fromName;
   Int_t           m_fromDurability;
   string          m_toName;
   Int_t           m_toDurability;

   // List of branches
   TBranch        *b_EventMetaData_fUniqueID;   //!
   TBranch        *b_EventMetaData_fBits;   //!
   TBranch        *b_EventMetaData_m_event;   //!
   TBranch        *b_EventMetaData_m_run;   //!
   TBranch        *b_EventMetaData_m_subrun;   //!
   TBranch        *b_EventMetaData_m_experiment;   //!
   TBranch        *b_EventMetaData_m_production;   //!
   TBranch        *b_EventMetaData_m_time;   //!
   TBranch        *b_EventMetaData_m_parentLfn;   //!
   TBranch        *b_EventMetaData_m_generatedWeight;   //!
   TBranch        *b_EventMetaData_m_errorFlag;   //!
   TBranch        *b_MicrotpcDataHits_;   //!
   TBranch        *b_MicrotpcDataHits_fUniqueID;   //!
   TBranch        *b_MicrotpcDataHits_fBits;   //!
   TBranch        *b_MicrotpcDataHits_m_backgroundTag;   //!
   TBranch        *b_MicrotpcDataHits_m_column;   //!
   TBranch        *b_MicrotpcDataHits_m_row;   //!
   TBranch        *b_MicrotpcDataHits_m_BCID;   //!
   TBranch        *b_MicrotpcDataHits_m_TOT;   //!
   TBranch        *b_MicrotpcDataHits_m_detNb;   //!
   TBranch        *b_MicrotpcMetaHits_;   //!
   TBranch        *b_MicrotpcMetaHits_fUniqueID;   //!
   TBranch        *b_MicrotpcMetaHits_fBits;   //!
   TBranch        *b_MicrotpcMetaHits_m_backgroundTag;   //!
   TBranch        *b_MicrotpcMetaHits_m_detNb;   //!
   TBranch        *b_MicrotpcMetaHits_m_pixNb;   //!
   TBranch        *b_MicrotpcMetaHits_m_ts_nb;   //!
   TBranch        *b_MicrotpcMetaHits_m_ts_start;   //!
   TBranch        *b_MicrotpcMetaHits_m_ts_stop;   //!
   TBranch        *b_MicrotpcMetaHits_m_Temperature;   //!
   TBranch        *b_MicrotpcMetaHits_m_Pressure;   //!
   TBranch        *b_MicrotpcMetaHits_m_Flow;   //!
   TBranch        *b_MicrotpcMetaHits_m_SetFlow;   //!
   TBranch        *b_MicrotpcMetaHits_m_GetFlow;   //!
   TBranch        *b_MicrotpcMetaHits_m_IHER;   //!
   TBranch        *b_MicrotpcMetaHits_m_PHER;   //!
   TBranch        *b_MicrotpcMetaHits_m_tHER;   //!
   TBranch        *b_MicrotpcMetaHits_m_flagHER;   //!
   TBranch        *b_MicrotpcMetaHits_m_ILER;   //!
   TBranch        *b_MicrotpcMetaHits_m_PLER;   //!
   TBranch        *b_MicrotpcMetaHits_m_tLER;   //!
   TBranch        *b_MicrotpcMetaHits_m_flagLER;   //!
   TBranch        *b_MicrotpcRecoTracks_;   //!
   TBranch        *b_MicrotpcRecoTracks_fUniqueID;   //!
   TBranch        *b_MicrotpcRecoTracks_fBits;   //!
   TBranch        *b_MicrotpcRecoTracks_m_backgroundTag;   //!
   TBranch        *b_MicrotpcRecoTracks_m_detNb;   //!
   TBranch        *b_MicrotpcRecoTracks_m_pixnb;   //!
   TBranch        *b_MicrotpcRecoTracks_m_chi2;   //!
   TBranch        *b_MicrotpcRecoTracks_m_theta;   //!
   TBranch        *b_MicrotpcRecoTracks_m_phi;   //!
   TBranch        *b_MicrotpcRecoTracks_m_esum;   //!
   TBranch        *b_MicrotpcRecoTracks_m_totsum;   //!
   TBranch        *b_MicrotpcRecoTracks_m_trl;   //!
   TBranch        *b_MicrotpcRecoTracks_m_time_range;   //!
   TBranch        *b_MicrotpcRecoTracks_m_parFit;   //!
   TBranch        *b_MicrotpcRecoTracks_m_parFit_err;   //!
   TBranch        *b_MicrotpcRecoTracks_m_cov;   //!
   TBranch        *b_MicrotpcRecoTracks_m_impact_x;   //!
   TBranch        *b_MicrotpcRecoTracks_m_impact_y;   //!
   TBranch        *b_MicrotpcRecoTracks_m_side;   //!
   TBranch        *b_MicrotpcRecoTracks_m_partID;   //!
   TBranch        *b_MicrotpcRecoTracksToMicrotpcDataHits_fUniqueID;   //!
   TBranch        *b_MicrotpcRecoTracksToMicrotpcDataHits_fBits;   //!
   TBranch        *b_MicrotpcRecoTracksToMicrotpcDataHits_m_elements_;   //!
   TBranch        *b_m_elements_fUniqueID;   //!
   TBranch        *b_m_elements_fBits;   //!
   TBranch        *b_m_elements_m_from;   //!
   TBranch        *b_m_elements_m_to;   //!
   TBranch        *b_m_elements_m_weights;   //!
   TBranch        *b_MicrotpcRecoTracksToMicrotpcDataHits_m_fromName;   //!
   TBranch        *b_MicrotpcRecoTracksToMicrotpcDataHits_m_fromDurability;   //!
   TBranch        *b_MicrotpcRecoTracksToMicrotpcDataHits_m_toName;   //!
   TBranch        *b_MicrotpcRecoTracksToMicrotpcDataHits_m_toDurability;   //!

   skimmer(TTree *tree=0);
   virtual ~skimmer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString FileName, TString OutputName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef skimmer_cxx
skimmer::skimmer(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tpc_data_1460214023.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tpc_data_1460214023.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

skimmer::~skimmer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skimmer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skimmer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void skimmer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_EventMetaData_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_EventMetaData_fBits);
   fChain->SetBranchAddress("m_event", &m_event, &b_EventMetaData_m_event);
   fChain->SetBranchAddress("m_run", &m_run, &b_EventMetaData_m_run);
   fChain->SetBranchAddress("m_subrun", &m_subrun, &b_EventMetaData_m_subrun);
   fChain->SetBranchAddress("m_experiment", &m_experiment, &b_EventMetaData_m_experiment);
   fChain->SetBranchAddress("m_production", &m_production, &b_EventMetaData_m_production);
   fChain->SetBranchAddress("m_time", &m_time, &b_EventMetaData_m_time);
   fChain->SetBranchAddress("m_parentLfn", &m_parentLfn, &b_EventMetaData_m_parentLfn);
   fChain->SetBranchAddress("m_generatedWeight", &m_generatedWeight, &b_EventMetaData_m_generatedWeight);
   fChain->SetBranchAddress("m_errorFlag", &m_errorFlag, &b_EventMetaData_m_errorFlag);
   fChain->SetBranchAddress("MicrotpcDataHits", &MicrotpcDataHits_, &b_MicrotpcDataHits_);
   fChain->SetBranchAddress("MicrotpcDataHits.fUniqueID", MicrotpcDataHits_fUniqueID, &b_MicrotpcDataHits_fUniqueID);
   fChain->SetBranchAddress("MicrotpcDataHits.fBits", MicrotpcDataHits_fBits, &b_MicrotpcDataHits_fBits);
   fChain->SetBranchAddress("MicrotpcDataHits.m_backgroundTag", MicrotpcDataHits_m_backgroundTag, &b_MicrotpcDataHits_m_backgroundTag);
   fChain->SetBranchAddress("MicrotpcDataHits.m_column", MicrotpcDataHits_m_column, &b_MicrotpcDataHits_m_column);
   fChain->SetBranchAddress("MicrotpcDataHits.m_row", MicrotpcDataHits_m_row, &b_MicrotpcDataHits_m_row);
   fChain->SetBranchAddress("MicrotpcDataHits.m_BCID", MicrotpcDataHits_m_BCID, &b_MicrotpcDataHits_m_BCID);
   fChain->SetBranchAddress("MicrotpcDataHits.m_TOT", MicrotpcDataHits_m_TOT, &b_MicrotpcDataHits_m_TOT);
   fChain->SetBranchAddress("MicrotpcDataHits.m_detNb", MicrotpcDataHits_m_detNb, &b_MicrotpcDataHits_m_detNb);
   fChain->SetBranchAddress("MicrotpcMetaHits", &MicrotpcMetaHits_, &b_MicrotpcMetaHits_);
   fChain->SetBranchAddress("MicrotpcMetaHits.fUniqueID", MicrotpcMetaHits_fUniqueID, &b_MicrotpcMetaHits_fUniqueID);
   fChain->SetBranchAddress("MicrotpcMetaHits.fBits", MicrotpcMetaHits_fBits, &b_MicrotpcMetaHits_fBits);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_backgroundTag", MicrotpcMetaHits_m_backgroundTag, &b_MicrotpcMetaHits_m_backgroundTag);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_detNb", MicrotpcMetaHits_m_detNb, &b_MicrotpcMetaHits_m_detNb);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_pixNb", MicrotpcMetaHits_m_pixNb, &b_MicrotpcMetaHits_m_pixNb);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_ts_nb", MicrotpcMetaHits_m_ts_nb, &b_MicrotpcMetaHits_m_ts_nb);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_ts_start[10]", MicrotpcMetaHits_m_ts_start, &b_MicrotpcMetaHits_m_ts_start);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_ts_stop[10]", MicrotpcMetaHits_m_ts_stop, &b_MicrotpcMetaHits_m_ts_stop);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_Temperature[4]", MicrotpcMetaHits_m_Temperature, &b_MicrotpcMetaHits_m_Temperature);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_Pressure[2]", MicrotpcMetaHits_m_Pressure, &b_MicrotpcMetaHits_m_Pressure);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_Flow[2]", MicrotpcMetaHits_m_Flow, &b_MicrotpcMetaHits_m_Flow);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_SetFlow", MicrotpcMetaHits_m_SetFlow, &b_MicrotpcMetaHits_m_SetFlow);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_GetFlow", MicrotpcMetaHits_m_GetFlow, &b_MicrotpcMetaHits_m_GetFlow);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_IHER", MicrotpcMetaHits_m_IHER, &b_MicrotpcMetaHits_m_IHER);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_PHER", MicrotpcMetaHits_m_PHER, &b_MicrotpcMetaHits_m_PHER);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_tHER", MicrotpcMetaHits_m_tHER, &b_MicrotpcMetaHits_m_tHER);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_flagHER", MicrotpcMetaHits_m_flagHER, &b_MicrotpcMetaHits_m_flagHER);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_ILER", MicrotpcMetaHits_m_ILER, &b_MicrotpcMetaHits_m_ILER);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_PLER", MicrotpcMetaHits_m_PLER, &b_MicrotpcMetaHits_m_PLER);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_tLER", MicrotpcMetaHits_m_tLER, &b_MicrotpcMetaHits_m_tLER);
   fChain->SetBranchAddress("MicrotpcMetaHits.m_flagLER", MicrotpcMetaHits_m_flagLER, &b_MicrotpcMetaHits_m_flagLER);
   fChain->SetBranchAddress("MicrotpcRecoTracks", &MicrotpcRecoTracks_, &b_MicrotpcRecoTracks_);
   fChain->SetBranchAddress("MicrotpcRecoTracks.fUniqueID", MicrotpcRecoTracks_fUniqueID, &b_MicrotpcRecoTracks_fUniqueID);
   fChain->SetBranchAddress("MicrotpcRecoTracks.fBits", MicrotpcRecoTracks_fBits, &b_MicrotpcRecoTracks_fBits);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_backgroundTag", MicrotpcRecoTracks_m_backgroundTag, &b_MicrotpcRecoTracks_m_backgroundTag);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_detNb", MicrotpcRecoTracks_m_detNb, &b_MicrotpcRecoTracks_m_detNb);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_pixnb", MicrotpcRecoTracks_m_pixnb, &b_MicrotpcRecoTracks_m_pixnb);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_chi2", MicrotpcRecoTracks_m_chi2, &b_MicrotpcRecoTracks_m_chi2);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_theta", MicrotpcRecoTracks_m_theta, &b_MicrotpcRecoTracks_m_theta);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_phi", MicrotpcRecoTracks_m_phi, &b_MicrotpcRecoTracks_m_phi);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_esum", MicrotpcRecoTracks_m_esum, &b_MicrotpcRecoTracks_m_esum);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_totsum", MicrotpcRecoTracks_m_totsum, &b_MicrotpcRecoTracks_m_totsum);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_trl", MicrotpcRecoTracks_m_trl, &b_MicrotpcRecoTracks_m_trl);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_time_range", MicrotpcRecoTracks_m_time_range, &b_MicrotpcRecoTracks_m_time_range);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_parFit[5]", MicrotpcRecoTracks_m_parFit, &b_MicrotpcRecoTracks_m_parFit);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_parFit_err[5]", MicrotpcRecoTracks_m_parFit_err, &b_MicrotpcRecoTracks_m_parFit_err);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_cov[25]", MicrotpcRecoTracks_m_cov, &b_MicrotpcRecoTracks_m_cov);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_impact_x[4]", MicrotpcRecoTracks_m_impact_x, &b_MicrotpcRecoTracks_m_impact_x);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_impact_y[4]", MicrotpcRecoTracks_m_impact_y, &b_MicrotpcRecoTracks_m_impact_y);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_side[16]", MicrotpcRecoTracks_m_side, &b_MicrotpcRecoTracks_m_side);
   fChain->SetBranchAddress("MicrotpcRecoTracks.m_partID[6]", MicrotpcRecoTracks_m_partID, &b_MicrotpcRecoTracks_m_partID);
//    fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_MicrotpcRecoTracksToMicrotpcDataHits_fUniqueID);
//    fChain->SetBranchAddress("fBits", &fBits, &b_MicrotpcRecoTracksToMicrotpcDataHits_fBits);
   fChain->SetBranchAddress("m_elements", &m_elements_, &b_MicrotpcRecoTracksToMicrotpcDataHits_m_elements_);
   fChain->SetBranchAddress("m_elements.fUniqueID", &m_elements_fUniqueID, &b_m_elements_fUniqueID);
   fChain->SetBranchAddress("m_elements.fBits", &m_elements_fBits, &b_m_elements_fBits);
   fChain->SetBranchAddress("m_elements.m_from", &m_elements_m_from, &b_m_elements_m_from);
   fChain->SetBranchAddress("m_elements.m_to", &m_elements_m_to, &b_m_elements_m_to);
   fChain->SetBranchAddress("m_elements.m_weights", &m_elements_m_weights, &b_m_elements_m_weights);
   fChain->SetBranchAddress("m_fromName", &m_fromName, &b_MicrotpcRecoTracksToMicrotpcDataHits_m_fromName);
   fChain->SetBranchAddress("m_fromDurability", &m_fromDurability, &b_MicrotpcRecoTracksToMicrotpcDataHits_m_fromDurability);
   fChain->SetBranchAddress("m_toName", &m_toName, &b_MicrotpcRecoTracksToMicrotpcDataHits_m_toName);
   fChain->SetBranchAddress("m_toDurability", &m_toDurability, &b_MicrotpcRecoTracksToMicrotpcDataHits_m_toDurability);
   Notify();
}

Bool_t skimmer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skimmer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skimmer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef skimmer_cxx
