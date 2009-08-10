//--------------------------------------------------------------------------------------------------
// $Id: PDFProducer.h,v 1.2 2009/04/07 15:37:09 loizides Exp $
//
// PDFProducer
//
// Computes PDFs
//
// Authors: G. Gomez-Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PDFProducer_H
#define MITPHYSICS_MODS_PDFProducer_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCEventInfo.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class PDFProducer : public BaseMod
  {
    public:
      PDFProducer(const char *name="PDFProducer", 
                   const char *title="PDF information module");
      ~PDFProducer();

      void           SetPrintDebug(bool b)             { fPrintDebug = b;    }
      void           SetMCEventInfoName(const char *s) { fMCEvInfoName  = s; }
      void           SetPDFName(const char *s)         { fPDFName = s;       }

    protected:
      void               Process();
      void               SlaveBegin();

      Bool_t             fPrintDebug;       //=true then print debug info
      TString            fMCEvInfoName;     //event info branch name
      const MCEventInfo *fMCEventInfo;      //!event info branch pointer
      TString            fPDFName;          //PDF name
      TH1D              *hDPDFHisto[10];    //!output histograms

    ClassDef(PDFProducer, 1) // Module to produce PDFs
  };
}
#endif
