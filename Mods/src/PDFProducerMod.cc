// $Id: PDFProducerMod.cc,v 1.2 2009/08/11 10:50:38 loizides Exp $

#include "MitPhysics/Mods/interface/PDFProducerMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include "LHAPDF/LHAPDF.h"

using namespace mithep;

ClassImp(mithep::PDFProducerMod)

//--------------------------------------------------------------------------------------------------
PDFProducerMod::PDFProducerMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fMCEvInfoName(Names::gkMCEvtInfoBrn),
  fPDFName("cteq65.LHgrid"),
  fMCEventInfo(0)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
void PDFProducerMod::Process()
{
  // Process entries of the tree.

  LoadBranch(fMCEvInfoName);

  Double_t Q    = fMCEventInfo->Scale();
  Int_t    id1  = fMCEventInfo->Id1();
  Double_t x1   = fMCEventInfo->X1();
  Double_t pdf1 = fMCEventInfo->Pdf1();
  Int_t    id2  = fMCEventInfo->Id2();
  Double_t x2   = fMCEventInfo->X2();
  Double_t pdf2 = fMCEventInfo->Pdf2();
  
  UInt_t nmembers = LHAPDF::numberPDF() + 1;

  // Array to be filled
  FArrDouble *PDFArr = new FArrDouble(nmembers);

  if (fPrintDebug) 
    cout << "Start loop over PDF members:" << endl;

  for (UInt_t i=0; i<nmembers; ++i) {
    LHAPDF::usePDFMember(i);
    Double_t newpdf1 = LHAPDF::xfx(x1, Q, id1)/x1;
    Double_t newpdf2 = LHAPDF::xfx(x2, Q, id2)/x2;
    Double_t TheWeight = newpdf1/pdf1*newpdf2/pdf2;
    
    if (fPrintDebug) {
      cout << i << " --> " << pdf1 << " "   << pdf2 << " | "
                           << x1   << " "   << x2   << " | "
			   << id1  << " "   << id2  << " | "
			   << Q    << " : " <<  TheWeight << endl;
    }
    if (GetFillHist()) {
      hDPDFHisto[0]->Fill(TheWeight);
    }
    PDFArr->Add(TheWeight);
  }

  AddObjThisEvt(PDFArr, "PDFWeights");
}

//--------------------------------------------------------------------------------------------------
void PDFProducerMod::SlaveBegin()
{
  // Setup LHAPDF and book branch and histograms if wanted.

  gSystem->Setenv("LHAPATH","");
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(fPDFName.Data());
  LHAPDF::getDescription();

  ReqBranch(fMCEvInfoName, fMCEventInfo);

  if (GetFillHist()) {
    char sb[1024];
    sprintf(sb,"hDPDFHisto_%d", 0);  
    hDPDFHisto[0] = new TH1D(sb,sb,1,0,1); 
    AddOutput(hDPDFHisto[0]);
  }
}
