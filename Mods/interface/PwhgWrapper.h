#ifndef PWHEGWRAPPER_H_
#define PWHEGWRAPPER_H_

/*
c     INPUT
c     mh : Higgs boson mass (used in the POWHEG BOX generation)
c     gh : Higgs boson width (used in the POWHEG BOX generation)
c     mt : top quark mass
c     BWflag : 0    if the sample to reweight was produced with fixed Higgs width
c              1    if the sample to reweight was produced with running Higgs 
c                   width (this is the default in the POWHEG BOX)
c     m : virtuality of the produced Higgs boson resonance
c     OUTPUT
c     w : the reweighting factor 
*/


class pwhegwrapper{
public: 
  pwhegwrapper(){};
  ~pwhegwrapper(){};
  double getweight(double mh,double gh,double mt,double m,int BWflag);
};


#endif
