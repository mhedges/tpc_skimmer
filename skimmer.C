/******************************************************************************
 * Program: TPC Phase I data skimmer
 * Author: Michael Hedges
 * Purpose: Skim Igal's BASF2 output for relevant data to be stored in BEAST
 * global ntuples for analysis
******************************************************************************/

#define skimmer_cxx
#include "skimmer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <algorithm>
//#include <thread>
#include <iostream>
#include <algorithm>
#include <vector>
#include <sys/time.h>

//#include "constants.h"

//#include "cantProceed.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TMath.h"
#include <TVector3.h>
#include <TVirtualFitter.h>

#define DEBUG 0
#define NROWS 336
#define NCOLS 80
#define NPERFILE 1000
#define ARRSIZE 2400

float weight[3] = {1., 1., 1.};
int iTPC=0;

void skimmer::fitTrack(TGraph2D *m_gr) {
  /* Find long aspect of track */
  int x_max_index, y_max_index;
  int x_min_index, y_min_index;
  double *x_vals, *y_vals, *z_vals;
  int npoints;
  npoints = m_gr->GetN();
  x_vals = m_gr->GetX(); y_vals = m_gr->GetY(); z_vals = m_gr->GetZ();

  for (int ii=0; ii<npoints; ii++){
	if ( x_vals[ii] == m_gr->GetXmax() )
	  x_max_index = ii;
	else if ( x_vals[ii] == m_gr->GetXmin() )
	  x_min_index = ii;
	if ( y_vals[ii] == m_gr->GetYmax() )
	  y_max_index = ii;
	else if ( x_vals[ii] == m_gr->GetYmin() )
	  y_min_index = ii;
  }

  int p_min_idx = 0, p_max_idx = 0;
  if ( ( m_gr->GetXmax() - m_gr->GetXmin() ) >=
	  ( m_gr->GetYmax() - m_gr->GetYmin() ) ) {
    p_min_idx = x_min_index;
    p_max_idx = x_max_index;
 } else {
    p_min_idx = y_min_index;
    p_max_idx = y_max_index;
  }
  
  // start fit track
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter *min = TVirtualFitter::Fitter(0, 5);  // Fitting with theta and phi
  min -> SetObjectFit(m_gr);
  min -> SetFCN(SumDistance2);
  
  // MAKE QUIET
  double p1=-1;
  min->ExecuteCommand("SET PRINTOUT",&p1, 1);

  double arglist[6] = {-1, 0, 0, 0, 0, 0};

  TVector3 temp_vector3 (x_vals[p_max_idx]-x_vals[p_min_idx],
			 y_vals[p_max_idx]-y_vals[p_min_idx],
			 z_vals[p_max_idx]-z_vals[p_min_idx]);
  double init_theta = temp_vector3.Theta();
  double init_phi   = temp_vector3.Phi();

  double pStart[5] = {x_vals[p_min_idx], y_vals[p_min_idx], z_vals[p_min_idx],
	init_theta, init_phi};
  min -> SetParameter(0, "x0",    pStart[0], 0.01, 0, 0);
  min -> SetParameter(1, "y0",    pStart[1], 0.01, 0, 0);
  min -> SetParameter(2, "z0",    pStart[2], 0.01, 0, 0);
  min -> SetParameter(3, "theta", pStart[3], 0.0001, 0, 0);
  min -> SetParameter(4, "phi",   pStart[4], 0.0001, 0, 0);
  
  arglist[0] = 1000; // number of fucntion calls
  arglist[1] = 0.01; // tolerance
  min -> ExecuteCommand("MIGRAD", arglist, 2);
  
  for (int iPar = 0; iPar < 5; iPar++){
    par_fit[iPar] = min -> GetParameter(iPar);
    par_fit_err[iPar] = min -> GetParError(iPar);
  }
  
  double theta, phi;
  phi   = par_fit[4]*180./3.14159;
  theta = par_fit[3]*180./3.14159;
  
  getTrackInfo();
  //getPID(); // Set PID flags

  /* Second arg comes from defs in constants.h */

  delete min;
}

void skimmer::SumDistance2(int &, double *, double & sum, double * par, int ) {
  TGraph2D * m_gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
  //assert (m_gr != 0);
  double * px = m_gr->GetX();
  double * py = m_gr->GetY();
  double * pz = m_gr->GetZ();
  int np = m_gr->GetN();
  sum = 0;
  for (int i  = 0; i < np; ++i) {
    double d = distance2(px[i],py[i],pz[i],par);
    sum += d;
  }
}

double skimmer::distance2(double px,double py,double pz, double *p) {
   TVector3 xp(px,py,pz);
   TVector3 x0(p[0], p[1], p[2]);
   TVector3 u (TMath::Sin(p[3])*TMath::Cos(p[4]), TMath::Sin(p[3])*TMath::Sin(p[4]), TMath::Cos(p[3]));

   double coeff = u*(xp-x0);
   TVector3 n = xp - x0 - coeff * u;

   double dx = n.x();
   double dy = n.y();
   double dz = n.z();
   double d2_x = TMath::Power(dx/weight[0], 2);
   double d2_y = TMath::Power(dy/weight[1], 2);
   double d2_z = TMath::Power(dz/weight[2], 2);
   double d2 = d2_x + d2_y + d2_z;

   return d2;
}

void skimmer::getTrackInfo(){
  // Get track length
  TVector3 fit_position(par_fit[0], par_fit[1], par_fit[2]);
  TVector3 unit_direction(TMath::Sin(par_fit[3])*TMath::Cos(par_fit[4]), TMath::Sin(par_fit[3])*TMath::Sin(par_fit[4]), TMath::Cos(par_fit[3]));
  
  float t_const = -1e-10;
  float x_point = fit_position.x() + t_const * unit_direction.x();
  float y_point = fit_position.y() + t_const * unit_direction.y();
  float z_point = fit_position.z() + t_const * unit_direction.z();
  TVector3 initial_position(x_point, y_point, z_point);
  
  float min_distance = 1e+30;
  float max_distance = -1e+30;
  float distance = 0.0;

  for (unsigned int iPoint = 0; iPoint < npoints; iPoint++){
    TVector3 position (col[iPoint], row[iPoint], bcid[iPoint]);
    TVector3 position_vect = position - initial_position;
    distance = unit_direction * position_vect;
    if (distance < min_distance) min_distance = distance;
    if (distance > max_distance) max_distance = distance;	
  }
  
  t_length = max_distance - min_distance;  

  // Get impact parameters
  //                        x=0    x=end    y=0    y=end
  float end_positions[4] = {0.0, 250.*80., 0.0, 50.*336.};
  for (int iPos = 0; iPos < 4; iPos++){
    t_const = 0.0;
    if ( iPos < 2 ) t_const = (end_positions[iPos]-fit_position.X())/unit_direction.X();
    else t_const = (end_positions[iPos]-fit_position.Y())/unit_direction.Y();
    
    // it might need to change
    if ( iPos < 2 ) impact_pars[iPos] = y_point;
    else impact_pars[iPos] = x_point; 
  } 
}
//
void skimmer::Loop(TString FileName, TString OutputName)
{
   TFile df(FileName,"READ");
   //
   // Get the TTree
   TTree *dtr = (TTree*)df.Get("tree");

   //TFile *ofile = new TFile("tpc_data.root", "RECREATE");

   TFile *ofile = new TFile(OutputName, "RECREATE");
   TTree *tr = new TTree("tr","TPC Event Data");

   //// Initialize the TTree
   Init(dtr);

   //// Activate all TBranches
   dtr->SetBranchStatus("*",1);

   //// Make TGraph for fitting the event
   TGraph2D *m_gr;
   //
   //// Make new TTree with relevant data
   tr->Branch("npoints",&npoints,"npoints/I");
   tr->Branch("row",&row,"row[npoints]/I");
   tr->Branch("col",&col,"col[npoints]/I");
   tr->Branch("bcid",&bcid,"bcid[npoints]/I");
   tr->Branch("tot",&tot,"tot[npoints]/I");
   tr->Branch("tstamp",&tstamp,"tstamp/D");
   tr->Branch("tot_sum",&tot_sum,"tot_sum/I");
   tr->Branch("sum_e",&sum_e,"sum_e/F");
   tr->Branch("time_range",&time_range,"time_range/I");
   tr->Branch("chi2",&chi2,"chi2/F");
   tr->Branch("t_length",&t_length,"t_length/F");
   tr->Branch("theta",&theta,"theta/F");
   tr->Branch("phi",&phi,"phi/F");
   tr->Branch("par_fit",&par_fit,"par_fit[6]/F");
   tr->Branch("par_fit_err",&par_fit_err,"par_fit_err[6]/F");
   tr->Branch("hitside",&hitside,"hitside[4]/I");
   tr->Branch("impact_pars",&impact_pars,"impact_pars[4]/F");

   int nentries = dtr->GetEntriesFast();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if (jentry %10000 == 0) cout << "Event Counter: " << jentry << endl;

      Long64_t ientry = LoadTree(jentry);

	  getentry = dtr->GetEntry(jentry);

	  npoints = MicrotpcMetaHits_m_pixNb[0];
	  tstamp = MicrotpcMetaHits_m_ts_start[0][0];
	  tot_sum = MicrotpcRecoTracks_m_totsum[0];
	  time_range = MicrotpcRecoTracks_m_time_range[0];
	  chi2 = MicrotpcRecoTracks_m_chi2[0];
	  t_length = MicrotpcRecoTracks_m_trl[0];
	  theta = MicrotpcRecoTracks_m_theta[0];
	  phi = MicrotpcRecoTracks_m_phi[0];
      sum_e = MicrotpcRecoTracks_m_esum[0];
	  
	  for (int pixn=0; pixn<npoints;pixn++){
	    row[pixn]=MicrotpcDataHits_m_row[pixn];
	    col[pixn]=MicrotpcDataHits_m_column[pixn];
	    bcid[pixn]=MicrotpcDataHits_m_BCID[pixn];
        cout << col[pixn] << row[pixn] << bcid[pixn] << endl;
	    tot[pixn]=MicrotpcDataHits_m_TOT[pixn];
		//m_gr->SetPoint(pixn, x, y, z);
	  

	  //for (int nsides=0; nsides<4; nsides++){
	  //   hitside[nsides] = MicrotpcRecoTracks_m_side[0][nsides];
	  //}

	  //for (int npars=0; npars<6;npars++){
	  //  par_fit[npars]=MicrotpcRecoTracks_m_parFit[0][npars];
	  //  par_fit_err[npars]=MicrotpcRecoTracks_m_parFit_err[0][npars];
	  //}

	   
	  }
	  double x,y,z;
	  for (int pixn=0; pixn<npoints;pixn++){
		x = col[pixn]; y = row[pixn]; z = bcid[pixn];
        cout << x << y << z<< endl;
		m_gr->SetPoint(pixn, x, y, z);
        cout << "Does it make it this far?" << endl;
	  }
	  //
	  
	  //Call fitter 
	  fitTrack(m_gr);

	  //Fill tree
	  tr->Fill();

	  //Delete track
	  m_gr -> Delete();
   }
   tr->Write();
   ofile->Write();
   ofile->Close();
}


//int main()
int main(int argc, char * argv[])
{
  skimmer s;
  //s.Loop("tpc4_th50_data_1463360400.root","test.root");
  s.Loop(argv[1], argv[2]);
  return 1;
}
