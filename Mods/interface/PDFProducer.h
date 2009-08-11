//--------------------------------------------------------------------------------------------------
// $Id: PDFProducer.h,v 1.1 2009/08/10 20:27:47 ceballos Exp $
//
// PDFProducer
//
// Computes PDFs from LHAPDF based on input from MCEventInfo.
//
// Authors: G.Gomez-Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PDFProducer_H
#define MITPHYSICS_MODS_PDFProducer_H

#include "MitAna/TreeMod/interface/BaseMod.h" 

class TH1D;
class TH2D;

namespace mithep 
{
  class MCEventInfo;
  class PDFProducer : public BaseMod
  {
    public:
      PDFProducer(const char *name="PDFProducerMod", 
                  const char *title="PDF producer module");

      void               SetMCEventInfoName(const char *s)   { fMCEvInfoName = s; }
      void               SetPDFName(const char *s)           { fPDFName      = s; }
      void               SetPrintDebug(Bool_t b)             { fPrintDebug   = b; }

    protected:
      void               Process();
      void               SlaveBegin();

      Bool_t             fPrintDebug;       //=true then print debug info
      TString            fMCEvInfoName;     //event info branch name
      TString            fPDFName;          //PDF name
      const MCEventInfo *fMCEventInfo;      //!event info branch pointer
      TH1D              *hDPDFHisto[10];    //!output histograms

    ClassDef(PDFProducer, 1) // Module to produce PDFs
  };
}
#endif
