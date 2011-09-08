#ifndef PhotonFix_Defined_hh
#define PhotonFix_Defined_hh

//-------------------------------------------------------//
// Project:  PhotonFix
// Author:   Paul Dauncey (p.dauncey@imperial.ac.uk)
// Modified: 11/07/2011
// Admins:   Paul Dauncey (p.dauncey@imperial.ac.uk)
//           Matt Kenzie (matthew.william.kenzie@cern.ch)
//-------------------------------------------------------//

/*
  Does post-reco fixes to ECAL photon energy and estimates resolution.
  This can run outside of the usual CMS software framework but requires 
  access to a file 'EcalGaps.dat' which must be in the same directory as 
  that used to run.

  To run within CMSSW use PhotonFixCMS.h (which can access the geometry 
  directly - go to "RecoEcal/EgammaCoreTools/plugins/PhotonFixCMS.h"
  for details.

  Before instantiating any objects of PhotonFix, the constants must be
  initialised in the first event using
    PhotonFix::initialise("3_8");
  
  The string gives the reco version used. Valid strings are 
  "3_8", "3_11", "4_2" and "Nominal", where the latter gives no correction 
  to the energy and a nominal resolution value. There is also "4_2e" which 
  provides corrections for electrons which are reconstructed as photons (to
  aid with testing the performance of these corrections in data).

  Make objects using
    PhotonFix a(energy,eta,phi,r9);
  where energy is the photon energy, eta and phi are the ECAL
  cluster positions (NB from the Supercluster object, _not_ the
  Photon object, as the latter gives eta and phi directions,
  not positions), and r9 is the R9 value of the SC.

  Get the corrected energy using
    a.fixedEnergy();
  and the resolution using
    a.sigmaEnergy();

*/

#include <iostream>
#include <string>

class PhotonFix {
 public:
  //PhotonFix(double e, double eta, double phi, double r9);

  PhotonFix() : _onePi(acos(-1.0)), _twoPi(2.0*acos(-1.0)), _initialised(false) {}
  // Must be called before instantiating any PhotonFix objects
  bool initialise(const std::string &s="Nominal", const std::string &infile="PhotonFix.dat");
  bool initialised() ;
  
  void setup(double e, double eta, double phi, double r9);
  
  bool isbarrel() const { return !_be; }
  
  // Corrected energy and sigma
  double fixedEnergy() const;
  double sigmaEnergy() const;
  
  // Input values
  double rawEnergy() const;
  double eta() const;
  double phi() const;
  double r9() const;

  // Derived EB crystal, submodule and module relative coordinates
  double etaC() const;
  double etaS() const;
  double etaM() const;

  double phiC() const;
  double phiS() const;
  double phiM() const;

  // Derived EE zeta, crystal, subcrystal and D-module relative coordinates
  double xZ() const;
  double xC() const;
  double xS() const;
  double xM() const;

  double yZ() const;
  double yC() const;
  double yS() const;
  double yM() const;

  // Return arrays containing positions of ecal gaps
  void barrelCGap(unsigned i, unsigned j, unsigned k, double c);
  void barrelSGap(unsigned i, unsigned j, unsigned k, double c);
  void barrelMGap(unsigned i, unsigned j, unsigned k, double c);
  void endcapCrystal(unsigned i, unsigned j, bool c); 
  void endcapCGap(unsigned i, unsigned j, unsigned k, double c);
  void endcapSGap(unsigned i, unsigned j, unsigned k, double c);
  void endcapMGap(unsigned i, unsigned j, unsigned k, double c);
  
  void print() const;

  // Input and output the fit parameters
  void setParameters(unsigned be, unsigned hl, const double *p);
  void getParameters(unsigned be, unsigned hl, double *p);

  void dumpParameters(std::ostream &o);
  void printParameters(std::ostream &o);

  // Utility functions
  double GetaPhi(double f0, double f1) const;
  double asinh(double s) const;  
  void dumpGaps(std::ostream &o);
  
 private:


  // Used by above; do not call directly
  bool initialiseParameters(const std::string &s);
  bool initialiseGeometry(const std::string &s, const std::string &infile);

  
   
   
  // Utility functions
  double dPhi(double f0, double f1) const;
  double aPhi(double f0, double f1) const;

  double expCorrection(double a, const double *p) const;
  double gausCorrection(double a, const double *p) const;

  // Actual data for each instantiated object
  unsigned _be,_hl;
  double _e,_eta,_phi,_r9;
  double _aC,_aS,_aM,_bC,_bS,_bM;
  
  // Constants
  const double _onePi;
  const double _twoPi;
  
  // Initialisation flag
  bool _initialised;
  
  // Parameters for fixes
  double _meanScale[2][2][4];
  double _meanAT[2][2][4];
  double _meanAC[2][2][4];
  double _meanAS[2][2][4];
  double _meanAM[2][2][4];
  double _meanBT[2][2][4];
  double _meanBC[2][2][4];
  double _meanBS[2][2][4];
  double _meanBM[2][2][4];
  double _meanR9[2][2][4];
  
  // Parameters for resolution
  double _sigmaScale[2][2][4];
  double _sigmaAT[2][2][4];
  double _sigmaAC[2][2][4];
  double _sigmaAS[2][2][4];
  double _sigmaAM[2][2][4];
  double _sigmaBT[2][2][4];
  double _sigmaBC[2][2][4];
  double _sigmaBS[2][2][4];
  double _sigmaBM[2][2][4];
  double _sigmaR9[2][2][4];
  
  // EB gap positions
  double _barrelCGap[169][360][2];
  double _barrelSGap[33][180][2];
  double _barrelMGap[7][18][2];
  
  // EE crystal existence and gap positions
  bool   _endcapCrystal[100][100];
  double _endcapCGap[2][7080][2];
  double _endcapSGap[2][264][2];
  double _endcapMGap[2][1][2];

};

#endif
