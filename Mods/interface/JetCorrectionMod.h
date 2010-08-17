//--------------------------------------------------------------------------------------------------
// $Id: JetCorrectionMod.h,v 1.3 2010/05/10 15:15:43 bendavid Exp $
//
// JetCorrectionMod
//
// This module applies jet energy corrections on the fly at analysis level, using directly the
// FWLite oriented classes from CMSSW and the jet correction txt files from the release and/or
// checked out tags.  At the moment this is hardcoded to apply or re-apply L2+L3 corrections,
// but it is intended that this will become more configurable in the future.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETCORRECTIONMOD_H
#define MITPHYSICS_MODS_JETCORRECTIONMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/Jet.h"
#include "MitAna/DataCont/interface/Types.h"
 
class FactorizedJetCorrector;
namespace mithep 
{
  class JetCorrectionMod : public BaseMod
  {
    public:
      JetCorrectionMod(const char *name="JetCorrectionMod", 
               const char *title="Jet correction module");
     ~JetCorrectionMod();

      const char       *GetInputName()                 const      { return fJetsName;               }   
      const char       *GetCorrectedName()             const      { return GetCorrectedJetsName();  }     
      const char       *GetCorrectedJetsName()         const      { return fCorrectedJetsName;      }     
      const char       *GetOutputName()                const      { return GetCorrectedJetsName();  }
      void              AddCorrectionFromRelease(const std::string &path);
      void              AddCorrectionFromFile(const std::string &file);    
      void              SetCorrectedJetsName(const char *name)    { fCorrectedJetsName = name;      }     
      void              SetCorrectedName(const char *name)        { SetCorrectedJetsName(name);     }    
      void              SetInputName(const char *name)            { fJetsName = name;               }  
      void              SetOutputName(const char *name)           { SetCorrectedJetsName(name);     }          

    protected:
      void              SlaveBegin();
      void              Process();

      TString           fJetsName;              //name of jet collection (input)
      TString           fCorrectedJetsName;     //name of good jets collection (output)
      //std::vector<JetCorrectorParameters> fCorrectionParameters; //list of corrections files (full paths)
      std::vector<std::string> fCorrectionFiles; //list of jet correction files
      FactorizedJetCorrector *fJetCorrector;      //!CMSSW/FWLite jet corrections module
      BitMask8          fEnabledCorrectionMask; //bitmask of enabled corrections
      std::vector<Jet::ECorr> fEnabledCorrections; //vector of enabled corrections

      ClassDef(JetCorrectionMod, 1) // Jet identification module
  };
}
#endif
