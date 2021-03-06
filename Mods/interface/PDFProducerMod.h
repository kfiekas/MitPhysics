//--------------------------------------------------------------------------------------------------
// $Id: PDFProducerMod.h,v 1.2 2010/03/13 20:50:01 ceballos Exp $
//
// PDFProducerMod
//
// Computes PDFs from LHAPDF based on input from MCEventInfo.
//
// Authors: G.Gomez-Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PDFPRODUCERMOD_H
#define MITPHYSICS_MODS_PDFPRODUCERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

class TH1D;
class TH2D;

namespace mithep 
{
  class MCEventInfo;
  class PDFProducerMod : public BaseMod
  {
    public:
      PDFProducerMod(const char *name="PDFProducerMod", 
                  const char *title="PDF producer module");

      void               SetMCEventInfoName(const char *s)   { fMCEvInfoName = s; }
      void               SetPDFName(const char *s)           { fPDFName      = s; }
      void               SetPrintDebug(Bool_t b)             { fPrintDebug   = b; }
      void               SetRunPDF(Bool_t b)                 { fRunPDF       = b; }
      void               SetIsData(Bool_t b)                 { fIsData       = b; }	 

    protected:
      void               Process();
      void               SlaveBegin();

      Bool_t             fPrintDebug;       //=true then print debug info
      TString            fMCEvInfoName;     //event info branch name
      TString            fPDFName;          //PDF name
      Bool_t             fRunPDF;           //=true run PDFs
      Bool_t             fIsData;           //=true then it does nothing (def=0)
      const MCEventInfo *fMCEventInfo;      //!event info branch pointer
      TH1D              *hDPDFHisto[10];    //!output histograms

    ClassDef(PDFProducerMod, 1) // Module to produce PDFs
  };
}
#endif
