//--------------------------------------------------------------------------------------------------
// $Id: LeptonPlusIsoTrackSelMod.h,v 1.1 2009/03/22 09:04:13 loizides Exp $
//
// LeptonPlusIsoTrackSelMod
// 
// This module selects events containing one lepton and one isolated track
// Both gsfTracks and tracker tracks are considered.
//
// Authors: S. Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_SELMODS_LEPTONPLUSISOTRACKSELMOD_H
#define MITPHYSICS_SELMODS_LEPTONPLUSISOTRACKSELMOD_H

#include "MitAna/TreeMod/interface/BaseSelMod.h" 
#include <TH1D.h>

namespace mithep 
{
  class LeptonPlusIsoTrackSelMod : public BaseSelMod
  {
    public:
      LeptonPlusIsoTrackSelMod(const char *name="LeptonPlusIsoTrackSelMod", 
                    const char *title="Lepton plus isolated track selection module");

      const char              *GetLeptonColName()       const { return fLeptonColName;            }
      const char              *GetTrackerTrackColName() const { return fTrackerTrackColName;      }
      const char              *GetGsfTrackColName()     const { return fGsfTrackColName;          }
      Double_t                 GetLeptonEtaMin()        const { return fLeptonEtaMin;             }
      Double_t                 GetLeptonEtaMax()        const { return fLeptonEtaMax;             }
      Double_t                 GetLeptonPtMin()         const { return fLeptonPtMin;              }
      Double_t                 GetLeptonPtMax()         const { return fLeptonPtMax;              }
      Double_t                 GetTrackEtaMin()         const { return fTrackEtaMin;              }
      Double_t                 GetTrackEtaMax()         const { return fTrackEtaMax;              }
      Double_t                 GetTrackPtMin()          const { return fTrackPtMin;               }
      Double_t                 GetTrackPtMax()          const { return fTrackPtMax;               }
      void                     SetLeptonColName(const char *n)       { fLeptonColName=n;          }
      void                     SetTrackerTrackColName(const char *n) { fTrackerTrackColName=n;    }
      void                     SetGsfTrackColName(const char *n)     { fGsfTrackColName=n;        }
      void                     SetLeptonEtaMin(Double_t e)           { fLeptonEtaMin = e;         }
      void                     SetLeptonEtaMax(Double_t e)           { fLeptonEtaMax = e;         }
      void                     SetLeptonPtMin(Double_t pt)           { fLeptonPtMin = pt;         }
      void                     SetLeptonPtMax(Double_t pt)           { fLeptonPtMax = pt;         }
      void                     SetTrackEtaMin(Double_t e)            { fTrackEtaMin = e;          }
      void                     SetTrackEtaMax(Double_t e)            { fTrackEtaMax = e;          }
      void                     SetTrackPtMin(Double_t pt)            { fTrackPtMin = pt;          }
      void                     SetTrackPtMax(Double_t pt)            { fTrackPtMax = pt;          }

    protected:
      void                     Process();
      void                     SlaveBegin();

      TString                  fLeptonColName;        //name of input lepton collection
      TString                  fTrackerTrackColName;  //name of input lepton collection
      TString                  fGsfTrackColName;      //name of input lepton collection
      Double_t                 fLeptonPtMin;          //minimum pt required   (def = 0 GeV)
      Double_t                 fLeptonPtMax;          //maximum pt required    (def = 5000 GeV)
      Double_t                 fLeptonEtaMin;         //minimum eta required   (def = -10)
      Double_t                 fLeptonEtaMax;         //maximum eta required   (def = +10) 
      Double_t                 fTrackPtMin;           //minimum pt required    (def = 0 GeV)
      Double_t                 fTrackPtMax;           //maximum pt required    (def = 5000 GeV)
      Double_t                 fTrackEtaMin;          //minimum eta required   (def = -10)
      Double_t                 fTrackEtaMax;          //maximum eta required   (def = +10) 
      const ParticleCol       *fLeptonCol;            //!pointer to collection 
      const TrackCol          *fTrackerTrackCol;      //!pointer to collection 
      const TrackCol          *fGsfTrackCol;          //!pointer to collection 
      TH1D                    *fNAccCounters;         //!acceptance histogram

      ClassDef(LeptonPlusIsoTrackSelMod,1) // Lepton plus isolated track selection module
  };
}
#endif
