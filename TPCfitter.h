#pragma once

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph2D.h"
#include <vector>

#define MAXHITS 2400

struct trackInfo {
  float x[MAXHITS];
  float y[MAXHITS];
  float z[MAXHITS];
  float x_min, x_max, y_min, y_max;
  int iTrack; 
  int bcid_min, bcid_max;
  int x_min_idx, y_min_idx, x_max_idx, y_max_idx;
};

struct fitInfo {
  unsigned char  col[MAXHITS]; // unsigned char (unsigned char )
  unsigned short row[MAXHITS]; //unsigned short (unsigned short)
  unsigned char  tot[MAXHITS]; //unsigned char  (unsigned char )
  unsigned char  bcid[MAXHITS]; //unsigned char  (unsigned char )
  double timestamp;
  unsigned int npoints;
  float pars[5];
  float errs[5];
  float length;
  //  float energy; 
  float impact_pars[4]; //x=0, x=end, y=0, y=end
  float theta; // radian
  float phi;
  unsigned int sumTOT;
  unsigned int evtNum; // Not doing anything with this right now
  bool alphaFlag;
  bool neutronFlag;
  unsigned short hitside;
  unsigned short iTPC;
};


void fitTrack();
void getTrackInfo();
void getPID();
unsigned short getHitside();
void SumDistance2(int &, double *, double & sum, double * par, int );
double distance2(double px,double py,double pz, double *p);
void writeNtuples();
//void writeNtuples_thread();
