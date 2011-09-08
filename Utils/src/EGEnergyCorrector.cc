// $Id: EGEnergyCorrector.cc,v 1.12 2011/08/03 17:15:44 bendavid Exp $

#include "MitPhysics/Utils/interface/EGEnergyCorrector.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


ClassImp(mithep::EGEnergyCorrector)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
EGEnergyCorrector::EGEnergyCorrector() :
fReadereb(0),
fReaderee(0),
fMethodname("BDTG method"),
fIsInitialized(kFALSE),
fIsMC(kFALSE)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
EGEnergyCorrector::~EGEnergyCorrector()
{
  if (fReadereb) delete fReadereb;
  if (fReaderebvariance) delete fReaderebvariance;  
  if (fReaderee) delete fReaderee;
  if (fReadereevariance) delete fReadereevariance;
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::Initialize(Bool_t ismc, TString phfixstring, TString phfixfile, TString ebweights, TString ebvarweights, TString eeweights, TString eevarweights) {
    fIsInitialized = kTRUE;
    fIsMC = ismc;
    fPhFix.initialise(std::string(phfixstring),std::string(phfixfile));

    if (fReadereb) delete fReadereb;
    if (fReaderebvariance) delete fReaderebvariance;  
    if (fReaderee) delete fReaderee;
    if (fReadereevariance) delete fReadereevariance;    
    
    fReadereb = new TMVA::Reader( "!Color:!Silent:Error" );    
    fReaderebvariance = new TMVA::Reader( "!Color:!Silent:Error" );    
    fReaderee = new TMVA::Reader( "!Color:!Silent:Error" );    
    fReadereevariance = new TMVA::Reader( "!Color:!Silent:Error" );  
    
    TMVA::Reader *readers[4];
    readers[0]  = fReadereb;
    readers[1]  = fReaderebvariance;    
    readers[2] = fReaderee;
    readers[3] = fReadereevariance;
    
    
    for (UInt_t i=0; i<4; ++i) {
      readers[i]->AddVariable("ph.scrawe",&scrawe);
      readers[i]->AddVariable("ph.r9",&r9);
      readers[i]->AddVariable("ph.sceta",&sceta);
      readers[i]->AddVariable("ph.scphi",&scphi);
      readers[i]->AddVariable("ph.e5x5/ph.scrawe",&e5x5norm);

   
      readers[i]->AddSpectator("ph.isbarrel",&fSpectator);
      
      if (i<2) {
        readers[i]->AddVariable("ph.etac",&etac);
        readers[i]->AddVariable("ph.etas",&etas);
        readers[i]->AddVariable("ph.etam",&etam);
        readers[i]->AddVariable("ph.phic",&phic);
        readers[i]->AddVariable("ph.phis",&phis);
        readers[i]->AddVariable("ph.phim",&phim);
        readers[i]->AddSpectator("ph.xz",&fSpectator);
        readers[i]->AddSpectator("ph.xc",&fSpectator);
        readers[i]->AddSpectator("ph.xs",&fSpectator);
        readers[i]->AddSpectator("ph.xm",&fSpectator);
        readers[i]->AddSpectator("ph.yz",&fSpectator);
        readers[i]->AddSpectator("ph.yc",&fSpectator);
        readers[i]->AddSpectator("ph.ys",&fSpectator);
        readers[i]->AddSpectator("ph.ym",&fSpectator); 
      }
      else {
        readers[i]->AddVariable("(1.0-(!ismc)*0.072)*ph.scpse/ph.scrawe",&scpsenorm);
        readers[i]->AddSpectator("ph.etac",&fSpectator);
        readers[i]->AddSpectator("ph.etas",&fSpectator); 
        readers[i]->AddSpectator("ph.etam",&fSpectator);
        readers[i]->AddSpectator("ph.phic",&fSpectator);
        readers[i]->AddSpectator("ph.phis",&fSpectator);
        readers[i]->AddSpectator("ph.phim",&fSpectator);
        readers[i]->AddVariable("ph.xz",&xz);
        readers[i]->AddVariable("ph.xc",&xc);
        readers[i]->AddVariable("ph.xs",&xs);
        readers[i]->AddVariable("ph.xm",&xm);
        readers[i]->AddVariable("ph.yz",&yz);
        readers[i]->AddVariable("ph.yc",&yc);
        readers[i]->AddVariable("ph.ys",&ys);
        readers[i]->AddVariable("ph.ym",&ym);        
      }
      
      readers[i]->AddVariable("ph.hovere",&hovere);
      readers[i]->AddVariable("ph.scetawidth",&scetawidth);
      readers[i]->AddVariable("ph.scphiwidth",&scphiwidth);
      readers[i]->AddVariable("ph.sigietaieta",&sigietaieta);
      
      readers[i]->AddSpectator("ph.e5x5",&fSpectator);
      readers[i]->AddSpectator("ph.scpse",&fSpectator);   
      readers[i]->AddSpectator("ph.emax",&fSpectator);
      readers[i]->AddSpectator("ph.etop",&fSpectator);
      readers[i]->AddSpectator("ph.ebottom",&fSpectator);
      readers[i]->AddSpectator("ph.eleft",&fSpectator);
      readers[i]->AddSpectator("ph.eright",&fSpectator);
      readers[i]->AddSpectator("ph.ebottom",&fSpectator);
      readers[i]->AddSpectator("ph.e2x5max",&fSpectator);
      readers[i]->AddSpectator("ph.e2x5top",&fSpectator);
      readers[i]->AddSpectator("ph.e2x5bottom",&fSpectator);
      readers[i]->AddSpectator("ph.e2x5left",&fSpectator);
      readers[i]->AddSpectator("ph.e2x5right",&fSpectator);
      readers[i]->AddSpectator("ph.e2x5bottom",&fSpectator); 
      readers[i]->AddSpectator("ph.convp",&fSpectator); 
      readers[i]->AddSpectator("ph.convpt",&fSpectator); 
      readers[i]->AddSpectator("ph.convleadpt",&fSpectator); 
      readers[i]->AddSpectator("ph.gene",&fSpectator);   
      readers[i]->AddSpectator("ph.ispromptgen",&fSpectator); 
      readers[i]->AddSpectator("ph.hasphoton",&fSpectator);    
      readers[i]->AddSpectator("ph.haselectron",&fSpectator);    
      readers[i]->AddSpectator("vtxZ",&fSpectator); 
      readers[i]->AddSpectator("evt",&fSpectator); 
      readers[i]->AddSpectator("ismc",&fSpectator); 
      readers[i]->AddSpectator("ph.elept",&fSpectator);     
    }
    
    fReadereb->BookMVA("BDTG method",ebweights);
    fReaderebvariance->BookMVA("BDTG method",ebvarweights);
    fReaderee->BookMVA("BDTG method",eeweights);
    fReadereevariance->BookMVA("BDTG method",eevarweights);
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::CorrectEnergyWithError(Photon *p) {
  
  std::pair<double,double> correction = CorrectedEnergyWithError(p);
  FourVectorM mom = p->Mom();
  double scale = correction.first/mom.E();
  p->SetMom(scale*mom.X(), scale*mom.Y(), scale*mom.Z(), scale*mom.E());
  p->SetEnergyErr(correction.second);
  p->SetEnergyErrSmeared(correction.second);
  
  //printf("initial pt = %5f, regression pt = %5f\n",mom.Pt(),p->Pt());
  
  return;

}

//--------------------------------------------------------------------------------------------------
std::pair<double,double> EGEnergyCorrector::CorrectedEnergyWithError(const Photon *p) {
  
  const SuperCluster *s = p->SCluster();
  const BasicCluster *b = s->Seed();  
  
  r9 = p->R9();
  e5x5norm = p->E55()/s->RawEnergy();
  hovere = p->HadOverEm();
  sigietaieta = p->CoviEtaiEta();  
  
  double scpse = s->PreshowerEnergy();
  if (!fIsMC) {
    scpse *= (1.0-0.072);
  }
  
  scrawe = s->RawEnergy();
  scpsenorm = scpse/s->RawEnergy();
  sceta = s->Eta();
  scphi = s->Phi();
  scetawidth = s->EtaWidth();
  scphiwidth = s->PhiWidth();
  
  Bool_t isbarrel = (s->AbsEta()<1.5);
  
  fPhFix.setup(p->E(),sceta,scphi,r9); 
  const Float_t dval = -99.;
  if (fPhFix.isbarrel()) {
    etac = fPhFix.etaC();
    etas = fPhFix.etaS();
    etam = fPhFix.etaM();
    phic = fPhFix.phiC();
    phis = fPhFix.phiS();
    phim = fPhFix.phiM();
    xz = dval;
    xc = dval;
    xs = dval;
    xm = dval;
    yz = dval;
    yc = dval;
    ys = dval;
    ym = dval;
  }
  else {
    etac = dval;
    etas = dval;
    etam = dval;
    phic = dval;
    phis = dval;
    phim = dval;
    xz = fPhFix.xZ();
    xc = fPhFix.xC();
    xs = fPhFix.xS();
    xm = fPhFix.xM();
    yz = fPhFix.yZ();
    yc = fPhFix.yC();
    ys = fPhFix.yS();
    ym = fPhFix.yM();
  }     
  
  const Double_t varscale = 1.253;
  Double_t den;
  TMVA::Reader *reader;
  TMVA::Reader *readervar;
  if (isbarrel) {
    den = scrawe;
    reader = fReadereb;
    readervar = fReaderebvariance;
  }
  else {
    den = scrawe + scpse;
    reader = fReaderee;
    readervar = fReadereevariance;
  }
  
  Double_t ecor = reader->EvaluateRegression(0,fMethodname)*den;
  Double_t ecorerr = readervar->EvaluateRegression(0,fMethodname)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}