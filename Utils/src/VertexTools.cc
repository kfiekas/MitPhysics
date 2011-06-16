// $Id: VertexTools.cc,v 1.1 2011/05/16 13:26:48 bendavid Exp $

#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TVector3.h>

ClassImp(mithep::VertexTools)

using namespace mithep;

VertexTools* VertexTools::meobject = NULL;

//--------------------------------------------------------------------------------------------------
VertexTools::VertexTools(const char* str)  
{
  // Constructor.
  relname = str;
  reader = new TMVA::Reader( "!Color:!Silent" );    
  reader->AddVariable( "var1", &tmvar1 );
  reader->AddVariable( "var2", &tmvar2 );
  reader->AddVariable( "var3", &tmvar3 );
  reader->AddVariable( "var4", &tmvar4 );
  reader->AddVariable( "var5", &tmvar5 );
  reader->AddVariable( "var6", &tmvar6 ); 
  reader->BookMVA( "BDTG method",relname + TString("/src/MitPhysics/data/TMVAClassification_BDTG.weights.xml").Data());
  reader->BookMVA( "BDTD method",relname + TString("/src/MitPhysics/data/TMVAClassification_BDTD.weights.xml" ).Data());
  //reader->BookMVA( "CFMlpANN method", "/home/maxi/cms/root/TMVAClassification_CFMlpANN.weights.xml" );
  reader->BookMVA( "MLP method", relname + TString("/src/MitPhysics/data/TMVAClassification_MLP.weights.xml").Data());
  reader->BookMVA( "MLPBFGS method",relname + TString("/src/MitPhysics/data/TMVAClassification_MLPBFGS.weights.xml" ).Data());
}

double VertexTools::NewMass(const Photon* ph1, const Photon* ph2, const BaseVertex* vert)
{
  ThreeVector drv1 = (ThreeVector(ph1->SCluster()->Point()) - vert->Position()).Unit();
  FourVector pho1c(drv1.X()*ph1->E(),drv1.Y()*ph1->E(),drv1.Z()*ph1->E(),ph1->E());
  ThreeVector drv2 = (ThreeVector(ph2->SCluster()->Point()) - vert->Position()).Unit();
  FourVector pho2c(drv2.X()*ph2->E(),drv2.Y()*ph2->E(),drv2.Z()*ph2->E(),ph2->E());

  FourVector diboso = pho1c+pho2c;
  return diboso.M();
}

//--------------------------------------------------------------------------------------------------
VertexZarray VertexTools::ExtractZarray(const VertexCol* vcol, float zmin, float zmax,  const BaseVertex  *fBeamSpot)
{
  VertexZarray zs;
  if(vcol == NULL) return zs;

  for(unsigned vv = 0; vv < vcol->GetEntries(); vv++){    
    const Vertex* vert = vcol->At(vv);
    double zpos = vert->Z();
    if(fBeamSpot != NULL) 
      zpos = zpos - fBeamSpot->Z();
    if(zpos > zmin && zpos > zmin)
      zs.push_back(zpos);
  }
  return zs;
}

VertexZarray VertexTools::ExtractZarray(float zmin, float zmax, float step)
{
  VertexZarray zs;
  for(float zpos = zmin; zpos < zmax+step; zpos = zpos+step){
    zs.push_back(zpos);
  }
  return zs;
}

const Vertex* VertexTools::BestVtx( const PFCandidateCol *fPFJets, const VertexCol *c, 
				    const BaseVertex  *fBeamSpot, FourVector diboso) {
  
  if (!c || !c->GetEntries()) return NULL;
  
  double bestprob = -100.;
  const Vertex* bestvert = NULL;
  for(unsigned vv = 0; vv < c->GetEntries(); vv++){
    
    const Vertex* vert = c->At(vv);
    double zpos = vert->Z();     
    double prob = Prob(fPFJets, zpos, fBeamSpot, diboso);
    if(prob > bestprob){
      bestprob = prob;
      bestvert = vert;
    }
  }
  return bestvert;
}

double VertexTools::BestVtx( const PFCandidateCol *fPFJets, VertexZarray zcol,
			     const BaseVertex  *fBeamSpot, FourVector diboso){
  
  double bestprob = -100.;
  double bestz = -100.;
  for(unsigned vv = 0; vv < zcol.size(); vv++){
    double zpos = zcol[vv];     
    double prob = Prob(fPFJets, zpos, fBeamSpot, diboso);
    if(prob > bestprob){
      bestprob = prob;
      bestz = zpos;
    }
  }
  return bestz;
  
}

double VertexTools::Prob(const PFCandidateCol *fPFJets, double zpos,  
			 const BaseVertex  *fBeamSpot, FourVector diboso){
  
  double bosophi = diboso.Phi();
  double bosopt = diboso.Pt();
  
  Vertex* vert = new Vertex(0,0,zpos);
  double ZVC = vert->Z()-fBeamSpot->Z();
  
  Double_t sinsum = 0.;
  Double_t cossum = 0.;
  Double_t sumpt = 0.;
  Double_t ntrks = 0.;
  Double_t ntplus = 0.;
  Double_t ortminus = 0.;
  Double_t ortplus = 0.;
  Double_t bdplus = 0.;
  Double_t bdminus = 0.;
  Double_t zmean = 0.;
  Double_t zmeansq = 0.;
  Double_t ww = 0.;
  
  for(unsigned pfj = 0; pfj < fPFJets->GetEntries(); pfj++){
    const PFCandidate* pfca = fPFJets->At(pfj);
    if(! pfca->HasTrk()) continue;
    if(!(pfca->PFType() == PFCandidate::eX ||
	 pfca->PFType() == PFCandidate::eHadron || 
	 pfca->PFType() == PFCandidate::eElectron || 
	 pfca->PFType() == PFCandidate::eMuon) ) continue;
    const Track *t = pfca->Trk();
    if(fabs(t->DzCorrected(*fBeamSpot) - ZVC ) > 0.2) continue;
    
    if(pfca->Pt()<0.3 || pfca->Pt()>200) continue;
    
    zmean = zmean + pfca->Pt()* t->DzCorrected(*fBeamSpot);
    zmeansq = zmeansq + pfca->Pt()*t->DzCorrected(*fBeamSpot)*t->DzCorrected(*fBeamSpot);
    ntrks++;
    sumpt  = sumpt + pfca->Pt()*pfca->Pt();
    ww = ww + pfca->Pt();
    
    sinsum = sinsum + pfca->Pt()*TMath::Sin(t->Phi());
    cossum = cossum + pfca->Pt()*TMath::Cos(t->Phi());
  }
  if(ntrks < 2 || !(sumpt > 0.) ) return 0;
  
  Double_t phim = TMath::ATan2(sinsum,cossum);
  zmean = zmean/ww;
  zmeansq = sqrt(zmeansq/ww - zmean*zmean);
    
  //--------------------
  Double_t bosoproj = 0.;
  Float_t ymean = 0.;
  Float_t ymsq = 0.;
  for(unsigned pfj = 0; pfj < fPFJets->GetEntries(); pfj++){
    const PFCandidate* pfca = fPFJets->At(pfj);
    if(! pfca->HasTrk()) continue;
    if(!(pfca->PFType() == PFCandidate::eX ||
	 pfca->PFType() == PFCandidate::eHadron ||
	 pfca->PFType() == PFCandidate::eElectron ||
	 pfca->PFType() == PFCandidate::eMuon) ) continue;
    
    const Track *t = pfca->Trk();
    if(fabs(t->DzCorrected(*fBeamSpot) - ZVC ) > 0.2) continue;
    //if(fabs(t->DzCorrected(*fBeamSpot) - ZVC ) > 3*zwidth) continue;
    if(pfca->Pt()<0.3 || pfca->Pt()>200) continue;
    
    //Float_t phid = phim - t->Phi();
    Float_t phid = bosophi+3.14 - t->Phi();
    ymean = ymean + pfca->Pt()*TMath::Sin(phid);
    ymsq = ymsq + pow(pfca->Pt()*TMath::Sin(phid),2);
    bosoproj = bosoproj + pfca->Pt()*TMath::Cos(bosophi-t->Phi());
    if(TMath::Sin(phid) > 0.) ortplus = ortplus + pfca->Pt()*TMath::Sin(phid);
    else                      ortminus = ortminus + pfca->Pt()*fabs(TMath::Sin(phid));
    if(TMath::Cos(phid) > 0.) bdplus = bdplus + pfca->Pt()*TMath::Cos(phid);
    else                      bdminus = bdminus + pfca->Pt()*fabs(TMath::Cos(phid));
    if(TMath::Cos(phid) > 0.) ntplus = ntplus + TMath::Cos(phid);
  }
  
  Float_t phimsq = sqrt(ymsq/ntrks - pow(ymean/ntrks,2))*180/3.14;
  Double_t A2n = (ortplus+ortminus) > 0.01 ? (bdplus+bdminus)/(ortplus+ortminus) : (bdplus+bdminus)/0.01;
  Double_t A1n = bdplus/bosopt;
  Double_t A3n = ntplus > 0 ? bdplus/ntplus : 0.;
  //-------------------
  Double_t angle = 180/3.14*TVector2::Phi_0_2pi(bosophi-phim);
  VertexTools* vtool = VertexTools::instance("");
  
  vtool->tmvar1 = ntrks;
  vtool->tmvar2 = sumpt;
  vtool->tmvar3 = A1n;
  vtool->tmvar4 = angle;
  vtool->tmvar5 = phimsq;
  vtool->tmvar6 = A3n;
  
  double p1 = vtool->reader->GetProba ( "BDTD method" );
  double p2 = vtool->reader->GetProba ( "BDTG method" );
  double p3 = vtool->reader->GetProba ( "MLP method");
  double p4 = vtool->reader->GetProba ( "MLPBFGS method");
  
  double retval = p1*p2*p3*p4;
  
  return retval;
}

double VertexTools::VertexWidth(const Vertex* vert,  const BaseVertex  *fBeamSpot){
  double width = 0.;
  double zmean = 0.;
  double zmeansq = 0.;
  double ww = 0.;
  if(vert == NULL) return width;
  for(unsigned i = 0; i < vert->NTracks(); i++){
    const Track *t = vert->Trk(i);
    if(t->Pt() < 0.3 ) continue;
    if(t->Pt() > 500.) continue;
    ww = ww + t->Pt();
    zmean = zmean + t->Pt()*t->DzCorrected(*fBeamSpot);
    zmeansq = zmeansq + t->Pt()*t->DzCorrected(*fBeamSpot)*t->DzCorrected(*fBeamSpot);
  }
  if( !(ww > 0.) ) return width;
  zmean = zmean/ww;

  width = sqrt(zmeansq/ww - zmean*zmean);
  return width;
  
}
