//--------------------------------------------------------------------------------------------------
// $Id: ElectronCorrectionMod.h,v 1.7 2010/03/12 13:51:26 bendavid Exp $
//
// ElectronCorrectionMod
//
// Takes in a collection of electrons, *Copies* each electron (so can't compare pointers before and
// afterward) into a new collection called ModNames::gkCorrectedElectronsName.
//
// The corrections come from a text file with four columns: etaMin, etaMax, scale, resolution, where
// the correction "scale" and "resolution" apply between etaMin and etaMax. There can be as many eta
// bins as desired. There *must* be a newline at the end of *every* line in the file.
//
//
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_ELECTRONCORRECTIONMOD_H
#define MITPHYSICS_MODS_ELECTRONCORRECTIONMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
class TH1D;
class TNtuple;
class TRandom3;

namespace mithep 
{
  class ElectronCorrectionMod : public BaseMod
  {
    public:
      ElectronCorrectionMod(const char *name="ElectronCorrectionMod", 
                     const char *title="Correcting the Electrons");

    protected:
      void                     Begin();
      void                     Process();
      void                     SlaveBegin();
      void                     SlaveTerminate();
      void                     Terminate();

      void                     SetElectronsFromBranch(Bool_t b) { fElectronsFromBranch = b; }
      void                     SetInputName(TString name)   { fInElectronName = name; }
      void                     SetOutputName(TString name)   { fCorrectedElectronsName = name; }
      const char               *GetOutputName() const { return fCorrectedElectronsName; }
      void                     SetDoNotCorrect(Bool_t b) { fDoNotCorrect=b;}
      void                     SetScaleFactorType(const char* type) { fCorrectionType=type; }
      void                     SetCorrectData(Bool_t b) { fCorrectData=b; }
      void                     SetSmearMC(Bool_t b) { fSmearMC=b; }
      void                     SetCorrectionFile(const char *file) { fCorrectionFileName=file;}
      void                     SetCorrectionValues(Electron *el);

      Bool_t                   fElectronsFromBranch;    //where to get input electrons
      Bool_t                   fDoNotCorrect;           //false: don't do anything, true: do the corrections
      Bool_t                   fIsMC;                   //is this MC?
      Bool_t                   fCorrectData;            //true: correct data, false: correct MC
      Bool_t                   fSmearMC;                //smear the MC when you correct it?
      TString                  fInElectronName;         //!name of input electron collection
      TString                  fCorrectionType;         // "FromJosh", etc
      TString                  fCorrectionFileName;     //File from which to get corrections
      TString                  fCorrectedElectronsName; //!
      const ElectronCol       *fElectrons;              //!Input electron branch
      Float_t                  fScale;                  //Transient scale value. fScale should be the 
                                                        //number by which you multiply *data*.
      Float_t                  fResolution;             //Transient resolution value
      TRandom3                *fRand;                   //
      TNtuple                 *fCorrections;            //The actual numbers

      ClassDef(ElectronCorrectionMod, 1) // Correct the electrons
  };
}
#endif

