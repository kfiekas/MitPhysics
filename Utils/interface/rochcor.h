#include <iostream>
#include <TChain.h>
#include <TClonesArray.h>
#include <TString.h>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>


#ifndef ROCHCOR
#define ROCHCOR

class rochcor {
 public:
  rochcor();
  rochcor(int seed);
  ~rochcor();
  
  void momcor_mc(TLorentzVector&, float, float, int);
  void momcor_data(TLorentzVector&, float, float, int);
  
  void musclefit_data(TLorentzVector& , TLorentzVector&);
  
  float zptcor(float);
  int etabin(float);
  int phibin(float);
  
 private:
  
  TRandom3 eran;
  TRandom3 sran;
  
  static const float netabin[9];  

  //---------------------------------------------------------------------------------------------
  
  static const float dcor_bfA[8][8];  
  static const float dcor_maA[8][8];
  static const float mcor_bfA[8][8];
  static const float mcor_maA[8][8];
  static const float dcor_bfAer[8][8];  
  static const float dcor_maAer[8][8];
  static const float mcor_bfAer[8][8];
  static const float mcor_maAer[8][8];

  static const float dcor_bfB[8][8];  
  static const float dcor_maB[8][8];
  static const float mcor_bfB[8][8];
  static const float mcor_maB[8][8];
  static const float dcor_bfBer[8][8];  
  static const float dcor_maBer[8][8];
  static const float mcor_bfBer[8][8];
  static const float mcor_maBer[8][8];

  //=======================================================================================================
  
  static const float dmavgA[8][8];  
  static const float dpavgA[8][8];  
  static const float mmavgA[8][8];  
  static const float mpavgA[8][8];

  static const float dmavgB[8][8];  
  static const float dpavgB[8][8];  
  static const float mmavgB[8][8];  
  static const float mpavgB[8][8];
  
  //===============================================================================================
  //parameters for Z pt correction
  static const int nptbins=84;
  static const float ptlow[85];    
  
  static const float zptscl[84];
  static const float zptscler[84];

  float mptsys_mc_dm[8][8];
  float mptsys_mc_da[8][8];
  float mptsys_da_dm[8][8];
  float mptsys_da_da[8][8];

};

#endif
