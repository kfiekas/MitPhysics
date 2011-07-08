// $Id: VertexTools.cc,v 1.4 2011/07/04 13:48:50 fabstoec Exp $

#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
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
  ThreeVector drv1 = (ThreeVector(ph1->CaloPos()) - vert->Position()).Unit();
  //ThreeVector drv1 = (ThreeVector(ph1->SCluster()->Point()) - vert->Position()).Unit();
  FourVector pho1c(drv1.X()*ph1->E(),drv1.Y()*ph1->E(),drv1.Z()*ph1->E(),ph1->E());
  ThreeVector drv2 = (ThreeVector(ph2->CaloPos()) - vert->Position()).Unit();
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
  VertexTools* vtool = VertexTools::instance("");
  
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

    std::vector<const Track*>::iterator itt;
    itt = find ((vtool->excluded).begin(), (vtool->excluded).end(), t);
    if(itt != (vtool->excluded).end()) continue;
    
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

    std::vector<const Track*>::iterator itt;
    itt = find ((vtool->excluded).begin(), (vtool->excluded).end(), t);
    if(itt != (vtool->excluded).end()) continue;
    
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
  // not used... commented out by Fabian
  //Double_t A2n = (ortplus+ortminus) > 0.01 ? (bdplus+bdminus)/(ortplus+ortminus) : (bdplus+bdminus)/0.01;
  Double_t A1n = bdplus/bosopt;
  Double_t A3n = ntplus > 0 ? bdplus/ntplus : 0.;
  //-------------------
  Double_t angle = 180/3.14*TVector2::Phi_0_2pi(bosophi-phim);
  
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

void VertexTools::BanThisTrack(const Track* track){
  VertexTools* vtool = VertexTools::instance("");
  (vtool->excluded).push_back(track);
}
void VertexTools::Reset(){
  VertexTools* vtool = VertexTools::instance("");
  (vtool->excluded).clear();
}


//------------------------------------------------------------------------------------
// The below tools are from the H->2photons EPS2011 BaseLine (common) Selection
const Vertex* VertexTools::findVtxBasicRanking(const Photon*           ph1, 
					       const Photon*           ph2, 
					       const BaseVertex*       bsp,
					       const VertexCol*        vtcs,
					       const DecayParticleCol* conv) {
  
  // check if all input is valid
  if( !ph1 || !ph2 || !bsp || !vtcs ) return NULL;
  // CAUTION: We allow for passing NULL for the Conversions, in that case only the simple Ranking is used.

  // here we will store the idx of the best Vtx
  unsigned int bestIdx = 0;

  // using asd much as possible 'Globe' naming schemes...
  int*    ptbal_rank  = new int   [vtcs->GetEntries()];
  int*    ptasym_rank = new int   [vtcs->GetEntries()];
  int*    total_rank  = new int   [vtcs->GetEntries()];
  double* ptbal       = new double[vtcs->GetEntries()];
  double* ptasym      = new double[vtcs->GetEntries()];
  
  unsigned int numVertices = vtcs->GetEntries();
  double       ptgg        = 0.;                     // stored for later in the conversion

  // loop over all the vertices...
  for(unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    
    const Vertex* tVtx = vtcs->At(iVtx);
    ptbal      [iVtx] = 0.0;
    ptasym     [iVtx] = 0.0;
    ptbal_rank [iVtx] = 1;
    ptasym_rank[iVtx] = 1;
    
    // compute the photon momenta with respect to this Vtx
    FourVectorM newMomFst = ph1->MomVtx(tVtx->Position());
    FourVectorM newMomSec = ph2->MomVtx(tVtx->Position());
    
    FourVectorM higgsMom = newMomFst+newMomSec; 

    double ph1Eta = newMomFst.Eta();
    double ph2Eta = newMomSec.Eta();

    double ph1Phi = newMomFst.Phi();
    double ph2Phi = newMomSec.Phi();
    
    // loop over all tracks and computew variables for ranking...
    FourVectorM totTrkMom(0,0,0,0);
    for(unsigned int iTrk = 0; iTrk < tVtx->NTracks(); ++iTrk) {
      const Track* tTrk = tVtx->Trk(iTrk);
      
      // compute distance between Trk and the Photons
      double tEta = tTrk->Eta();
      double tPhi = tTrk->Phi();
      double dEta1 = TMath::Abs(tEta-ph1Eta);
      double dEta2 = TMath::Abs(tEta-ph2Eta);
      double dPhi1 = TMath::Abs(tPhi-ph1Phi);
      double dPhi2 = TMath::Abs(tPhi-ph2Phi);
      if(dPhi1 > M_PI) dPhi1 = 2*M_PI - dPhi1;
      if(dPhi2 > M_PI) dPhi2 = 2*M_PI - dPhi2;

      double dR1 = TMath::Sqrt(dEta1*dEta1+dPhi1*dPhi1);
      double dR2 = TMath::Sqrt(dEta2*dEta2+dPhi2*dPhi2);
            
      if(dR1 < 0.05 || dR2 < 0.05) continue;
      totTrkMom = totTrkMom + tTrk->Mom4(0);
    }

    // compute the ranking variables for this Vtx
    double ptvtx   = totTrkMom.Pt();
    double pthiggs = higgsMom.Pt();
    ptbal [iVtx]  = (totTrkMom.Px()*(newMomFst.Px()+newMomSec.Px()));
    ptbal [iVtx] += (totTrkMom.Py()*(newMomFst.Py()+newMomSec.Py()));
    ptbal [iVtx]  = -ptbal[iVtx]/pthiggs;
    ptasym[iVtx]  = (ptvtx - pthiggs)/(ptvtx + pthiggs);
    
    // do the little ranking acrobatics...
    for(unsigned int cVtx =0; cVtx < iVtx; ++cVtx) {
      if(ptbal [iVtx] > ptbal [cVtx])
	ptbal_rank[cVtx]++;
      else
	ptbal_rank[iVtx]++;
      if(ptasym [iVtx] > ptasym [cVtx])
	ptasym_rank[cVtx]++;
      else
	ptasym_rank[iVtx]++;
    }
  }
  
  // loop again over all Vertcices (*sigh*), compute total score and final rank
  // CAUTION: Total rank starts at '0', so the best ranked Vtx has RANK 0
  for(unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    ptasym_rank [iVtx] = ptbal_rank [iVtx]*ptasym_rank [iVtx]*(iVtx+1);
    total_rank  [iVtx] = 0;
    for(unsigned int cVtx =0; cVtx < iVtx; ++cVtx) {
      if(ptasym_rank [iVtx] > ptasym_rank [cVtx]) total_rank[iVtx]++;
      else if(ptasym_rank [iVtx] == ptasym_rank [cVtx]) {                   // use 'ptbal' as the tie-breaker
	if(ptbal_rank [iVtx] > ptbal_rank [cVtx]) total_rank[iVtx]++;
	else total_rank[cVtx]++;
      }
      else total_rank[cVtx]++;
    }
  }
  
  // delete the auxiliary dynamic arrays
  delete[] ptbal_rank  ;
  delete[] ptasym_rank ;
  delete[] ptbal       ;
  delete[] ptasym      ;
  
  // find the best ranked Vertex so far....
  for(unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    if(total_rank[iVtx] == 0) bestIdx = iVtx;
  }
  
  // check if there's a conversion among the pre-selected photons
  // ...this will return NULL in case conv==NULL
  const DecayParticle* conv1 = PhotonTools::MatchedCiCConversion(ph1, conv, 0.1, 0.1, 0.1);
  const DecayParticle* conv2 = PhotonTools::MatchedCiCConversion(ph2, conv, 0.1, 0.1, 0.1);

  if( conv1 && ( conv1->Prob() < 0.0005) ) conv1 = NULL;
  if( conv2 && ( conv2->Prob() < 0.0005) ) conv2 = NULL;
  
  double zconv  = 0.;
  double dzconv = 0.;
  
  //--------------------------------------------------------------------
  // start doing the Conversion acrobatics... 'copied' from the Globe...
  if(conv1 || conv2) {
    if( !conv2 ){
      const mithep::ThreeVector caloPos1(ph1->CaloPos());
      zconv  = conv1->Z0EcalVtx(bsp->Position(), caloPos1);
      if( ph1->IsEB() ) {
	double rho = conv1->Position().Rho();
	if     ( rho < 15. ) dzconv = 0.06;
	else if( rho < 60. ) dzconv = 0.67;
	else                 dzconv = 2.04;
      } else {
	double z = conv1->Position().Z();
	if     ( TMath::Abs(z) < 50. )   dzconv = 0.18;
	else if( TMath::Abs(z) < 100.)   dzconv = 0.61;
	else                 dzconv = 0.99;
      }
    } else if( !conv1 ) {
      const mithep::ThreeVector caloPos2(ph2->CaloPos());
      zconv  = conv2->Z0EcalVtx(bsp->Position(), caloPos2);
      if( ph2->IsEB() ) {
	double rho = conv2->Position().Rho();
	if     ( rho < 15. ) dzconv = 0.06;
	else if( rho < 60. ) dzconv = 0.67;
	else                 dzconv = 2.04;
      } else {
	double z = conv2->Position().Z();
	if     ( TMath::Abs(z) < 50. )   dzconv = 0.18;
	else if( TMath::Abs(z) < 100.)   dzconv = 0.61;
	else                 dzconv = 0.99;
      }
    } else {
      const mithep::ThreeVector caloPos1(ph1->CaloPos());
      double z1  = conv1->Z0EcalVtx(bsp->Position(), caloPos1);
      double dz1 = 0.;
      if( ph1->IsEB() ) {
	double rho = conv1->Position().Rho();
	if     ( rho < 15. ) dz1 = 0.06;
	else if( rho < 60. ) dz1 = 0.67;
	else                 dz1 = 2.04;
      } else {
	double z = conv1->Position().Z();
	if     ( TMath::Abs(z) < 50. )   dz1 = 0.18;
	else if( TMath::Abs(z) < 100.)   dz1 = 0.61;
	else                 dz1 = 0.99;
      }
      const mithep::ThreeVector caloPos2(ph2->CaloPos());
      double z2  = conv2->Z0EcalVtx(bsp->Position(), caloPos2);
      double dz2 = 0.;
      if( ph2->IsEB() ) {
	double rho = conv2->Position().Rho();
	if     ( rho < 15. ) dz2 = 0.06;
	else if( rho < 60. ) dz2 = 0.67;
	else                 dz2 = 2.04;
      } else {
	double z = conv2->Position().Z();
	if     ( TMath::Abs(z) < 50. )   dz2 = 0.18;
	else if( TMath::Abs(z) < 100.)   dz2 = 0.61;
	else                 dz2 = 0.99;
      }
      zconv  = ( 1./(1./dz1/dz1 + 1./dz2/dz2 )*(z1/dz1/dz1 + z2/dz2/dz2) ) ;  // weighted average
      dzconv = TMath::Sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
    }

    // loop over all ranked Vertices and choose the closest to the Conversion one
    int maxVertices = ( ptgg > 30 ? 3 : 5);
    double minDz = -1.;    
    for(unsigned int iVtx =0; iVtx < numVertices; ++iVtx) {
      if(total_rank[iVtx] < maxVertices) {
	const Vertex* tVtx = vtcs->At(iVtx);
	double tDz = TMath::Abs(zconv - tVtx->Z());
	if( (minDz < 0. || tDz < minDz) && ( tDz < dzconv ) ) {	  
	  minDz = tDz;
	  bestIdx = iVtx;
	}
      }    
    }
  }
  // END of Conversion Acrobatics
  //--------------------------------------------------------------------

  delete[] total_rank  ;
  return vtcs->At(bestIdx);
}
