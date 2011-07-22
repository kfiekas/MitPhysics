// $Id: PUReweighting.cc,v 1.1 2011/06/01 12:36:59 mzanetti Exp $


#include "MitPhysics/Utils/interface/PUReweightingMulti.h"
#include <TAxis.h>

ClassImp(mithep::PUReweightingMulti)

using namespace mithep;


PUReweightingMulti::PUReweightingMulti(const THnSparse *sparse, int maxn, double minprob) :
fSparse(sparse),
fNsBins(0),
fNDim(0),
fMaxN(maxn+1),
fMinProb(minprob),
fWeights(0),
fIndices(0),
fCache(0)
{

  Setup();

}

PUReweightingMulti::~PUReweightingMulti() {
  
  if (fWeights) delete [] fWeights;
  if (fIndices) delete [] fIndices;
  if (fCache) delete [] fCache;
  if (fOffsets) delete [] fOffsets;
  if (fIndOffsets) delete [] fIndOffsets;
  if (fFactors) delete [] fFactors;

  
}

void PUReweightingMulti::Setup() {
  
  //cache all values needed to quickly computed probability density as sum over bins of the sparse histogram
  
  fNsBins = fSparse->GetNbins();
  fNDim = fSparse->GetNdimensions();
  
  fOffsets = new Int_t[fNDim];
  fIndOffsets = new Int_t[fNDim];
  fFactors = new Double_t[fNDim];
  
  fWeights = new Double_t[fNsBins];
  fIndices = new UShort_t[fNDim*fNsBins];

  int nbinstotal = 0;
  int *cacheoffset = new int[fNDim];
  for (int idim=0; idim<fNDim; ++idim) {
    cacheoffset[idim] = nbinstotal;
    nbinstotal += fSparse->GetAxis(idim)->GetNbins()+1;
    fIndOffsets[idim] = fNsBins*idim;
  }
  fNCacheBins = nbinstotal;
  fCache = new Double_t[fMaxN*fNCacheBins];
  for (int idim=0; idim<fNDim; ++idim) {
    int nbins = fSparse->GetAxis(idim)->GetNbins()+1;
    for (int ibin=0; ibin<nbins; ++ibin) {
      double lambda = (ibin>0) ? exp(fSparse->GetAxis(idim)->GetBinCenter(ibin)) : 0.0;
      for (int n=0; n<fMaxN; ++n) {
        int globalbin = cacheoffset[idim] + ibin;
        //fCache[fMaxN*globalbin + n] = TMath::Poisson(n,lambda);
        double val = TMath::Poisson(n,lambda);
        fCache[fNCacheBins*n + globalbin] = (val > fMinProb) ? val : 0.0;
      }
    }
  }  
  
  int *coords = new int[fNDim];
  double sumweights = 0;
  for (int i=0; i<fNsBins; ++i) {
    sumweights += fSparse->GetBinContent(i,coords);
  }

  for (int i=0; i<fNsBins; ++i) {
    double weight = fSparse->GetBinContent(i,coords);
    fWeights[i] = weight/sumweights;
    for (int idim=0; idim<fNDim; ++idim) {
      //fIndices[fNDim*i + idim] = cacheoffset[idim] + coords[idim];
      fIndices[fNsBins*idim + i] = cacheoffset[idim] + coords[idim];
    }
  }
  
  delete[] cacheoffset;
  delete[] coords;

}

double PUReweightingMulti::TargetPdf(int *npus, int ndim) const {
    
  assert(ndim<=fNDim);
  
  for (int idim=0; idim<ndim; ++idim) {
    if (npus[idim] >= fMaxN) return 0.0;
    fOffsets[idim] = fNCacheBins*npus[idim];
  }
  
  double totalweight = 0.0;
  for (int i=0; i<fNsBins; ++i) {
    bool nonzero = true;
    for (int idim=0; idim<ndim; ++idim) {
      fFactors[idim] = fCache[fOffsets[idim]+fIndices[fIndOffsets[idim]+i]];
      if ( !(fFactors[idim]>0.0) ) {
        nonzero = false;
        break;
      }
    }
    if (nonzero) {
      double weight = fWeights[i];
      for (int idim=0; idim<ndim; ++idim) {
        weight *= fFactors[idim];
      }
      totalweight += weight;
    }
  }
  
  return totalweight;
}

TH3D *PUReweightingMulti::Target3D() {
 assert(fNDim>=3);
  
 int np[3];
  
 TH3D *target3d = new TH3D("target3d","",fMaxN,-0.5,fMaxN-0.5,fMaxN,-0.5,fMaxN-0.5,fMaxN,-0.5,fMaxN-0.5);
 for (int i=0; i<fMaxN; ++i) {
   np[0] = i;
   for (int j=0; j<fMaxN; ++j) {
     np[1] = j;
     for (int k=0; k<fMaxN; ++k) {
       np[2] = k;
       target3d->Fill(i,j,k,TargetPdf(np,3));
     }
   }
 }
 return target3d;

}

TH3D *PUReweightingMulti::Weights3D(const TH3D *source3d) {
 assert(fNDim>=3);
  
 TH3D *source3dclone = (TH3D*)source3d->Clone();
 source3dclone->Scale(1.0/source3dclone->GetSumOfWeights());
  
 int np[3];
 TH3D *target3d = new TH3D("target3dtemp","",fMaxN,-0.5,fMaxN-0.5,fMaxN,-0.5,fMaxN-0.5,fMaxN,-0.5,fMaxN-0.5);
 for (int i=0; i<fMaxN; ++i) {
   np[0] = i;
   for (int j=0; j<fMaxN; ++j) {
     np[1] = j;
     for (int k=0; k<fMaxN; ++k) {
       np[2] = k;
       if (source3dclone->GetBinContent(source3dclone->FindFixBin(i,j,k))>0.0) {
         target3d->Fill(i,j,k,TargetPdf(np,3));
       }
       else {
         target3d->Fill(i,j,k,0.0);
       }
     }
   }
 }
 target3d->Scale(1.0/target3d->GetSumOfWeights());
 
 TH3D *weights3d = new TH3D(*target3d / *source3dclone);
 delete source3dclone;
 return weights3d;

}