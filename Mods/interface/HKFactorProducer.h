//------------------------------------------------------------------------------
// $Id: HKFactorProducer.h,v 1.21 2009/04/03 22:42:59 ceballos Exp $
//
// HKFactorProducer
//
// Produces K Factors from LO to NNLO
//
//
// Authors: G. Gomez-Ceballos
//------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_HKFACTORPRODUCER_H
#define MITPHYSICS_MODS_HKFACTORPRODUCER_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitPhysics/Mods/interface/HWWKFactorList.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class HKFactorProducer : public BaseMod
  {
    public:
      HKFactorProducer(const char *name="HKFactorProducer", 
                   const char *title="KFactor information module");
      ~HKFactorProducer();

      void     SetDebug(Bool_t b)	       { fDebug 	   = b; }
      void     SetFillHist(Bool_t b)	       { fFillHist 	   = b; }
      void     SetProcessID(Int_t d)	       { fProcessID	   = d; }
      void     SetInputFilename(const char *s) { fInputFilename    = s; }

      void     SetMCBosonsName(const char *s)  { fMCBosonsName     = s; }
      void     SetMCPhotonsName(const char *s) { fMCEventInfoName  = s; }

    protected:
      void     Process();
      void     SlaveBegin();

      Bool_t fDebug;
      Bool_t fFillHist;
      HWWKfactorList* fPt_histo;
      Int_t  fProcessID;
      std::string fInputFilename;

      TString               fMCBosonsName;
      TString               fMCEventInfoName;
      const MCEventInfo    *fMCEventInfo;

      TH1D                 *hDHKFactor[10];

    ClassDef(HKFactorProducer, 1) // Module to gather generator information
  };
}
#endif
