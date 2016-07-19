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

#include "TPCfitter.h"
#include "constants.h"

#include "cantProceed.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TMath.h"
#include <TVector3.h>
#include <TVirtualFitter.h>

#include "data_distributor.h"

#define DEBUG 0
#define NROWS 336
#define NCOLS 80
#define NPERFILE 1000
#define ARRSIZE 2400


//data_distributor *m_distributor = data_distributor::instance();
#define m_distributor (data_distributor::instance())


data_distributor::event_display_t evtDisplay;

struct timeval timeMark;
TGraph2D *m_gr;
waveformPV *m_waveformPV;
trackInfo *m_track = new trackInfo;   // Single instance of struct that will get overwritten
fitInfo m_fit[NPERFILE];              // NPERFILE instances to be dumped to nTuples
//fitInfo m_fit_copy[NPERFILE];      // Copy of m_fit to pass to writeNtuple thread

float weight[3] = {1., 1., 1.};
int iTPC=0;


/* Subscribe to TPC datastream. The thread never returns from this call. */
void tpcDatastreamSubscribe(char* arg) {
  char *pname;
  if(!strcmp(arg,"TPCmaster")){
    pname=epicsStrDup("BEAST:TPCMASTER:PACKED_TRACK");
    iTPC=3;
  }
  else {
    pname=epicsStrDup("BEAST:TPCDAQSTER:PACKED_TRACK");
    iTPC=2;
  }
  
  printf("\nBeginning subscription to %s\n\n", pname);
  m_track->iTrack=0;

  m_waveformPV = (waveformPV*)callocMustSucceed(1, sizeof(waveformPV), "Calloc must succeed"); 
  SEVCHK(ca_context_create(ca_disable_preemptive_callback),"ca_context_create");
  SEVCHK(ca_add_exception_event(exceptionCallback,NULL),"ca_add_exception_event");
  
  SEVCHK(ca_create_channel(pname,connectionCallback, m_waveformPV,20,&m_waveformPV->mychid),"ca_create_channel");
  
  m_waveformPV->data = (unsigned int*)malloc(dbr_size_n(DBR_LONG, 2400));
  
  SEVCHK(ca_replace_access_rights_event(m_waveformPV->mychid,accessRightsCallback),"ca_replace_access_rights_event");
  printf("...thread will not return\n");
  SEVCHK(ca_create_subscription(DBR_LONG,2400,m_waveformPV->mychid,DBE_VALUE,eventCallback,m_waveformPV,&m_waveformPV->myevid),"ca_create_subscription");
  
  SEVCHK(ca_pend_event(0.0),"ca_pend_event");
}

/* Called when subscribed PV captures an event */
void eventCallback(struct event_handler_args eha)
{
#if DEBUG
  printf("eventCallback %i\n", m_track->iTrack);
#endif
  chid chid = eha.chid;
  
  if(eha.status!=ECA_NORMAL) {
    printChidInfo(chid,"eventCallback");
  } 
  else {
    memcpy(m_waveformPV->data, eha.dbr, dbr_size_n(eha.type, eha.count));
    // printf("npoints: %i\n", (unsigned int)eha.count);

    // printf("CALLING UNPACK\n");
    unpackEvent();
    // printf("CALLING FIT\n");
    fitTrack(); // ABSOLUTELY DO NOT START NEW THREAD HERE
    //  printf("DONE WITH FIT\n");
    if((m_track->iTrack%10)==0) {
      printf("\nTrack %i of %i processed\n\n",m_track->iTrack, NPERFILE);
    }
    if(m_track->iTrack%(NPERFILE-1)==0 && m_track->iTrack>0) {
      //  printf("CALLING NTUPLE\n");
      writeNtuples(); 
      m_track->iTrack = 0;
    } else m_track->iTrack++;
  }
}


// Unpack single event from datastream
void unpackEvent() {
  m_fit[m_track->iTrack].npoints   = std::min((unsigned int)m_waveformPV->data[0],(unsigned int)ARRSIZE);
  //m_fit[m_track->iTrack].iTPC      = m_waveformPV->data[1];
  m_fit[m_track->iTrack].iTPC = iTPC;
  m_fit[m_track->iTrack].timestamp = (m_waveformPV->data[2])/1000.+1455662200; 

#if DEBUG
  printf("iTPC: %d, ts: %f, npoints: %d\n", m_fit[m_track->iTrack].iTPC, m_fit[m_track->iTrack].timestamp, m_fit[m_track->iTrack].npoints);
  
#endif

  m_gr = new TGraph2D(); /* Deleted in fitTrack() */

  m_track->bcid_min=1000000;
  m_track->bcid_max=-1;
  m_track->x_min=1000000;
  m_track->x_max=-1;
  m_track->y_min=1000000;
  m_track->y_max=-1;

  m_fit[m_track->iTrack].sumTOT=0;
  memset(evtDisplay.data(), 0, EVTDISPLENGTH);
  
  unsigned char tot, col, bcid;
  unsigned short row;


  for(unsigned int iPoint=0; iPoint<m_fit[m_track->iTrack].npoints; iPoint++) {
    tot  = (unsigned char)((m_waveformPV->data[iPoint+3] & 0x0000000F) >> 0);  
    row  = (unsigned short)((m_waveformPV->data[iPoint+3] & 0X00001FF0) >> 4);
    col  = (unsigned char)((m_waveformPV->data[iPoint+3] & 0x001FE000) >> 13);
    bcid = (unsigned char)((m_waveformPV->data[iPoint+3] & 0x1FE00000) >> 21);

    m_fit[m_track->iTrack].tot[iPoint]  = tot;
    m_fit[m_track->iTrack].row[iPoint]  = row;
    m_fit[m_track->iTrack].col[iPoint]  = col;
    m_fit[m_track->iTrack].bcid[iPoint] = bcid;

    //  printf("making evt display %i\n", iPoint);

    m_track->x[iPoint] = (float)col*20000./80.;
    m_track->y[iPoint] = (float)row*16800./336.;
    m_track->z[iPoint] = (float)bcid*250.;
    m_fit[m_track->iTrack].sumTOT+= tot;
    m_gr->SetPoint(iPoint,m_track->x[iPoint],m_track->y[iPoint],m_track->z[iPoint]);
    evtDisplay[col-1+NCOLS*(row-1)] = tot; // the -1 is because row and col start at 1

    //printf("B\n");
    if(m_track->x[iPoint]>m_track->x_max) {
      m_track->x_max = m_track->x[iPoint];
      m_track->x_max_idx = iPoint;
    }
    else if(m_track->x[iPoint]<m_track->x_min) {
      m_track->x_min = m_track->x[iPoint];
      m_track->x_min_idx = iPoint;
    }
    if(m_track->y[iPoint]>m_track->y_max) {
      m_track->y_max = m_track->y[iPoint];
      m_track->y_max_idx = iPoint;
    }
    else if(m_track->y[iPoint]<m_track->y_min) {
      m_track->y_min = m_track->y[iPoint];
      m_track->y_min_idx = iPoint;
    }
    if(bcid>m_track->bcid_max) m_track->bcid_max = bcid;
    else if(bcid<m_track->bcid_min) m_track->bcid_max = bcid;
  }
}

// Fit one track 
void fitTrack() {
  /* Find long aspect of track */
  int p_min_idx = 0, p_max_idx = 0;
  if ( (m_track->x_max-m_track->x_min) >= (m_track->y_max-m_track->y_min) ) {
    p_min_idx = m_track->x_min_idx;
    p_max_idx = m_track->x_max_idx;
 } else {
    p_min_idx = m_track->y_min_idx;
    p_max_idx = m_track->y_max_idx;
  }
  
  // start fit track
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter *min = TVirtualFitter::Fitter(0, 5);  // Fitting with theta and phi
  min -> SetObjectFit(m_gr);
  min -> SetFCN(SumDistance2);
  
  // MAKE QUIET
  double p1=-1;
  min ->ExecuteCommand("SET PRINTOUT",&p1, 1);

  double arglist[6] = {-1, 0, 0, 0, 0, 0};
  /*min -> ExecuteCommand("SET PRINT", arglist, 1);
  min -> ExecuteCommand("SET NOWARNINGS", arglist, 0);
  */

  TVector3 temp_vector3 (m_track->x[p_max_idx]-m_track->x[p_min_idx],
			 m_track->y[p_max_idx]-m_track->y[p_min_idx],
			 m_track->z[p_max_idx]-m_track->z[p_min_idx]);
  double init_theta = temp_vector3.Theta();
  double init_phi   = temp_vector3.Phi();

  double pStart[5] = {m_track->x[p_min_idx], m_track->y[p_min_idx], m_track->z[p_min_idx], init_theta, init_phi};
  min -> SetParameter(0, "x0",    pStart[0], 0.01, 0, 0);
  min -> SetParameter(1, "y0",    pStart[1], 0.01, 0, 0);
  min -> SetParameter(2, "z0",    pStart[2], 0.01, 0, 0);
  min -> SetParameter(3, "theta", pStart[3], 0.0001, 0, 0);
  min -> SetParameter(4, "phi",   pStart[4], 0.0001, 0, 0);
  
  arglist[0] = 1000; // number of fucntion calls
  arglist[1] = 0.01; // tolerance
  min -> ExecuteCommand("MIGRAD", arglist, 2);
  
  for (int iPar = 0; iPar < 5; iPar++){
    m_fit[m_track->iTrack].pars[iPar] = min -> GetParameter(iPar);
    m_fit[m_track->iTrack].errs[iPar] = min -> GetParError(iPar);
  }
  
  m_fit[m_track->iTrack].phi   = m_fit[m_track->iTrack].pars[4]*180./3.14159;
  m_fit[m_track->iTrack].theta = m_fit[m_track->iTrack].pars[3]*180./3.14159;
  
  getTrackInfo();
  getPID(); // Set PID flags

  /* Second arg comes from defs in constants.h */
  /* Make sure to set NaNs (or something) for wrong PID channels */

  m_distributor->setData(m_fit[m_track->iTrack].iTPC, chan_sumTOT,       m_fit[m_track->iTrack].sumTOT);
  m_distributor->setData(m_fit[m_track->iTrack].iTPC, chan_phi,          m_fit[m_track->iTrack].phi);
  m_distributor->setData(m_fit[m_track->iTrack].iTPC, chan_theta,        m_fit[m_track->iTrack].theta);
  m_distributor->setData(m_fit[m_track->iTrack].iTPC, chan_alpha_flag,   m_fit[m_track->iTrack].alphaFlag);
  m_distributor->setData(m_fit[m_track->iTrack].iTPC, chan_neutron_flag, m_fit[m_track->iTrack].neutronFlag);
  m_distributor->setData(m_fit[m_track->iTrack].iTPC, chanEvtDispl,      evtDisplay);

  m_distributor->incrementCount(m_fit[m_track->iTrack].iTPC, chan_total_rate);
  if(m_fit[m_track->iTrack].alphaFlag) 
    m_distributor->incrementCount(m_fit[m_track->iTrack].iTPC, chan_alpha_rate);
  else if(m_fit[m_track->iTrack].neutronFlag) 
    m_distributor->incrementCount(m_fit[m_track->iTrack].iTPC, chan_neutron_rate);


#if DEBUG
    printf("TPC[%d]: %d %d %d %d\n", m_fit[m_track->iTrack].iTPC, evtDisplay[0], evtDisplay[1], evtDisplay[2], evtDisplay[3]);
#endif

  m_distributor->channelDone(m_fit[m_track->iTrack].iTPC); /* Trigger processing of I/O Intr records */

  m_gr -> Delete();
  delete min;
}

void getPID() {
  m_fit[m_track->iTrack].hitside = getHitside();
  m_fit[m_track->iTrack].alphaFlag   = false;
  m_fit[m_track->iTrack].neutronFlag = false;

  if(m_fit[m_track->iTrack].hitside==11 && m_fit[m_track->iTrack].sumTOT>100)     m_fit[m_track->iTrack].alphaFlag = true;
  else if(m_fit[m_track->iTrack].hitside==0 && m_fit[m_track->iTrack].sumTOT>100) m_fit[m_track->iTrack].neutronFlag = true;

}

unsigned short getHitside() {
  int cut_dim = 500; //um (it should be 250*integer)
  int set_edge_row_up = 335-cut_dim/250*5, set_edge_row_dw = cut_dim/250*5;
  int set_edge_col_up = 79-cut_dim/250, set_edge_col_dw = cut_dim/250;
  
  int c_row_up = 0, c_row_dw = 0, c_col_up = 0, c_col_dw = 0;
  int c_row_up_col_up = 0, c_row_up_col_dw = 0, c_row_dw_col_up = 0, c_row_dw_col_dw = 0;
  int c_inside = 0;
  
  int ncorners = 0;
  int check_edge = 0, nhits_edge = 0;
  int icol, irow;
  int nedges;

  for (int iPoint = 0; iPoint < m_fit[m_track->iTrack].npoints; iPoint++){
    check_edge = 0;
    icol = m_fit[m_track->iTrack].col[iPoint];
    irow = m_fit[m_track->iTrack].row[iPoint];
    if ( irow > set_edge_row_up ) {
      c_row_up = 1;
      check_edge = 1;
    } 
    if ( irow < set_edge_row_dw ) {
      c_row_dw = 1;
      check_edge = 1;
    } 
    if ( icol > set_edge_col_up ) {
      c_col_up = 1;
      check_edge = 1;
    }
    if ( icol < set_edge_col_dw ) {
      c_col_dw = 1;
      check_edge = 1;
    }
    
    if ( check_edge == 1)  nhits_edge += 1;
    
    if ( irow > set_edge_row_up && icol > set_edge_col_up ) 
      c_row_up_col_up = 1;
    if ( irow > set_edge_row_up && icol < set_edge_col_dw ) 
      c_row_up_col_dw = 1;
    if ( irow < set_edge_row_dw && icol > set_edge_col_up ) 
      c_row_dw_col_up = 1;
    if ( irow < set_edge_row_dw && icol < set_edge_col_dw ) 
      c_row_dw_col_dw = 1;
    if (! ( irow < set_edge_row_dw || irow > set_edge_row_up || icol < set_edge_col_dw || icol > set_edge_col_up ) )
      c_inside = 1;
    
  }
  nedges = c_row_up + c_row_dw + c_col_up + c_col_dw;
  ncorners = c_row_up_col_up + c_row_up_col_dw + c_row_dw_col_up + c_row_dw_col_dw;
  nedges = nedges - ncorners;
  
  if ( c_inside == 1 ){
    if ( c_row_up_col_up+c_row_dw_col_up > 1 || 
	 c_row_up_col_dw+c_row_dw_col_dw > 1 || 
	 c_row_dw_col_up+c_row_dw_col_dw > 1 || 
	 c_row_up_col_up+c_row_up_col_dw > 1 ) 
      nedges += 1;
  }
  if ( nedges >= 2 && c_inside == 1){
    nedges = 2;
  }

  return (unsigned short)(c_row_up*1000 + c_row_dw*100 + c_col_up*10 + c_col_dw);
}

void printChidInfo(chid chid, const char *message)
{
  printf("\n%s\n",message);
  printf("pv: %s  type(%d) nelements(%ld) host(%s)",
	 ca_name(chid),ca_field_type(chid),ca_element_count(chid),
	 ca_host_name(chid));
  printf(" read(%d) write(%d) state(%d)\n",
	 ca_read_access(chid),ca_write_access(chid),ca_state(chid));
}

void exceptionCallback(struct exception_handler_args args)
{
  chid	chid = args.chid;
  long	stat = args.stat; /* Channel access status code*/
  const char  *channel;
  static char *noname = "unknown";
  
  channel = (chid ? ca_name(chid) : noname);
  
  if(chid) printChidInfo(chid,"exceptionCallback");
  printf("exceptionCallback stat %s channel %s\n",
	 ca_message(stat),channel);
}

void connectionCallback(struct connection_handler_args args)
{
    chid	chid = args.chid;

    printChidInfo(chid,"connectionCallback");
}

void accessRightsCallback(struct access_rights_handler_args args)
{
    printChidInfo(args.chid,"accessRightsCallback");
}


void SumDistance2(int &, double *, double & sum, double * par, int ) {
  TGraph2D * m_gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
  assert (m_gr != 0);
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

double distance2(double px,double py,double pz, double *p) {
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

void writeNtuples() {
  //  memcpy(m_fit, m_fit_copy, sizeof(m_fit));
  
  //  std::thread write_thread(writeNtuples_thread, m_fit_copy);
  //  write_thread.detach();

  // Not calling a thread saves copying large blocks of unused memory (all the [MAXHITS] arrays)
  // MAXHITS*5*NPERFILE bytes = 1.2MB for 100 tracks (regardless of nhits)
  // For right now, I won't even save hit data
  
  TFile *outfile = new TFile(Form("/mnt/iscsi/data/NTP/TPC/TPConline_%F.root", m_fit[0].timestamp),"RECREATE");
  std::cout << "Writing ntuples: " << Form("/mnt/iscsi/data/NTP/TPC/TPConline_%F.root", m_fit[0].timestamp) << std::endl;
  TTree *tout = new TTree("tout","tout");

  fitInfo tmpFit; // Use for copying memory version

  //                             Address of the pointers to the members of the struct
  tout->Branch("npoints",        &(tmpFit.npoints),     "npoints/i");
  tout->Branch("ts",             &(tmpFit.timestamp),   "ts/D");
  tout->Branch("sumTOT",         &(tmpFit.sumTOT),      "sumTOT/i");
  tout->Branch("length",         &(tmpFit.length),      "length/F");
  //  tout->Branch("energy",         &(tmpFit.energy),      "energy/F");
  tout->Branch("pars",           &(tmpFit.pars),        "pars[5]/F");
  tout->Branch("errs",           &(tmpFit.errs),        "errs[5]/F");
  tout->Branch("impact_pars",    &(tmpFit.impact_pars), "impact_pars[4]/F");
  tout->Branch("theta",          &(tmpFit.theta),       "theta/F");
  tout->Branch("phi",            &(tmpFit.phi),         "phi/F");
  tout->Branch("alphaFlag",      &(tmpFit.alphaFlag),   "alphaFlag/O");       // Bool_t
  tout->Branch("neutronFlag",    &(tmpFit.neutronFlag), "neutronFlag/O");     // Bool_t
  tout->Branch("iTPC",           &(tmpFit.iTPC),        "iTPC/b");            // UChar_t
  tout->Branch("hitside",        &(tmpFit.hitside),     "hitside/s");         // UShort_t
  //  tout->Branch("col",            &(tmpFit.col),         "col[npoints]/b");    // UChar_t
  //  tout->Branch("row",            &(tmpFit.row),         "row[npoints/s");     // UShort_t
  //  tout->Branch("tot",            &(tmpFit.tot),         "tot[npoints]/b");    // UChar_t
  //  tout->Branch("bcid",           &(tmpFit.bcid),        "bcid[npoints]/b");   // UChar_t

  // This copies each instance of the struct in memory to the static location of tmpFit. 
  // Performance-wise that's not great, but I don't see a way around it. 
  for(int iFit=0; iFit<NPERFILE; ++iFit) { 
    tmpFit = m_fit[iFit];
    tout->Fill(); 
  }
  
  tout->Write();
  outfile->Close();
}

// void writeNtuples_thread(fitInfo* fits) {
//   TFile *outfile = new TFile(Form("/mnt/iscsi/DATA/NTP/TPC/TPConline_%F.root", fits[0].timestamp),"RECREATE");
//   std::cout << "Writing ntuples: " << Form("/mnt/iscsi/DATA/NTP/TPC/TPConline_%F.root", fits[0].timestamp) << std::endl;
//   TTree *tout = new TTree("tout","tout");

//   fitInfo tmpFit;

//   tout->Branch("npoints",        &(tmpFit.npoints),     "npoints/i");
//   tout->Branch("ts",             &(tmpFit.timestamp),   "ts/D");
//   tout->Branch("sumTOT",         &(tmpFit.sumTOT),      "sumTOT/i");
//   tout->Branch("length",         &(tmpFit.length),      "length/F");
//   tout->Branch("energy",         &(tmpFit.energy),      "energy/F");
//   tout->Branch("pars",           &(tmpFit.pars),        "pars[5]/F");
//   tout->Branch("errs",           &(tmpFit.errs),        "errs[5]/F");
//   tout->Branch("impact_pars[4]", &(tmpFit.impact_pars), "impact_pars[4]/F");
//   tout->Branch("theta",          &(tmpFit.theta),       "theta/F");
//   tout->Branch("phi",            &(tmpFit.phi),         "phi/F");
//   tout->Branch("alphaFlag",      &(tmpFit.alphaFlag),   "alphaFlag/O");       // Bool_t
//   tout->Branch("neutronFlag",    &(tmpFit.neutronFlag), "neutronFlag/O");     // Bool_t
//   tout->Branch("iTPC",           &(tmpFit.iTPC),        "iTPC/b");            // UChar_t

//   for(int iFit=0; iFit<NPERFILE; ++iFit) { 
//     tmpFit = fits[iFit];
//     tout->Fill(); 
//   }
  
//   tout->Write();
//   outfile->Close();
// }

void getTrackInfo() {
  // Get track length
  TVector3 fit_position(m_fit[m_track->iTrack].pars[0], m_fit[m_track->iTrack].pars[1], m_fit[m_track->iTrack].pars[2]);
  TVector3 unit_direction (TMath::Sin(m_fit[m_track->iTrack].pars[3])*TMath::Cos(m_fit[m_track->iTrack].pars[4]), TMath::Sin(m_fit[m_track->iTrack].pars[3])*TMath::Sin(m_fit[m_track->iTrack].pars[4]), TMath::Cos(m_fit[m_track->iTrack].pars[3]));
  
  float t_const = -1e-10;
  float x_point = fit_position.x() + t_const * unit_direction.x();
  float y_point = fit_position.y() + t_const * unit_direction.y();
  float z_point = fit_position.z() + t_const * unit_direction.z();
  TVector3 initial_position(x_point, y_point, z_point);
  
  float min_distance = 1e+30;
  float max_distance = -1e+30;
  float distance = 0.0;

  for (unsigned int iPoint = 0; iPoint < m_fit[m_track->iTrack].npoints; iPoint++){
    TVector3 position (m_track->x[iPoint], m_track->y[iPoint], m_track->z[iPoint]);
    TVector3 position_vect = position - initial_position;
    distance = unit_direction * position_vect;
    if (distance < min_distance) min_distance = distance;
    if (distance > max_distance) max_distance = distance;	
  }
  
  m_fit[m_track->iTrack].length = max_distance - min_distance;  


  // Get impact parameters
  //                        x=0    x=end    y=0    y=end
  float end_positions[4] = {0.0, 250.*80., 0.0, 50.*336.};
  for (int iPos = 0; iPos < 4; iPos++){
    t_const = 0.0;
    if ( iPos < 2 ) t_const = (end_positions[iPos]-fit_position.X())/unit_direction.X();
    else t_const = (end_positions[iPos]-fit_position.Y())/unit_direction.Y();
    
    // it might need to change
    if ( iPos < 2 ) m_fit[m_track->iTrack].impact_pars[iPos] = y_point;
    else m_fit[m_track->iTrack].impact_pars[iPos] = x_point; 
  } 
}
