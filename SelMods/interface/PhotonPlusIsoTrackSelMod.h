//--------------------------------------------------------------------------------------------------
// $Id $
//
// PhotonPlusIsoTrackSelMod
// 
// This module selects events containing one lepton and one isolated track
// Both gsfTracks and tracker tracks are considered.
//
// Authors: S. Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_SELMODS_PHOTONPLUSISOTRACKSELMOD_H
#define MITPHYSICS_SELMODS_PHOTONPLUSISOTRACKSELMOD_H

#include "MitAna/TreeMod/interface/BaseSelMod.h" 
#include "MitAna/DataTree/interface/PhotonFwd.h" 
#include "MitAna/DataTree/interface/TrackFwd.h" 
class TH1D;

namespace mithep 
{
  class PhotonPlusIsoTrackSelMod : public BaseSelMod
  {
    public:
      PhotonPlusIsoTrackSelMod(const char *name="PhotonPlusIsoTrackSelMod", 
                               const char *title="Photon plus isolated track selection module");

      const char              *GetPhotonColName()       const { return fPhotonColName;            }
      const char              *GetTrackerTrackColName() const { return fTrackerTrackColName;      }
      const char              *GetGsfTrackColName()     const { return fGsfTrackColName;          }
      Double_t                 GetPhotonEtaMin()        const { return fPhotonEtaMin;             }
      Double_t                 GetPhotonEtaMax()        const { return fPhotonEtaMax;             }
      Double_t                 GetPhotonPtMin()         const { return fPhotonPtMin;              }
      Double_t                 GetPhotonPtMax()         const { return fPhotonPtMax;              }
      Double_t                 GetTrackEtaMin()         const { return fTrackEtaMin;              }
      Double_t                 GetTrackEtaMax()         const { return fTrackEtaMax;              }
      Double_t                 GetTrackPtMin()          const { return fTrackPtMin;               }
      Double_t                 GetTrackPtMax()          const { return fTrackPtMax;               }
      void                     SetPhotonColName(const char *n)       { fPhotonColName=n;          }
      void                     SetTrackerTrackColName(const char *n) { fTrackerTrackColName=n;    }
      void                     SetGsfTrackColName(const char *n)     { fGsfTrackColName=n;        }
      void                     SetPhotonEtaMin(Double_t e)           { fPhotonEtaMin = e;         }
      void                     SetPhotonEtaMax(Double_t e)           { fPhotonEtaMax = e;         }
      void                     SetPhotonPtMin(Double_t pt)           { fPhotonPtMin = pt;         }
      void                     SetPhotonPtMax(Double_t pt)           { fPhotonPtMax = pt;         }
      void                     SetTrackEtaMin(Double_t e)            { fTrackEtaMin = e;          }
      void                     SetTrackEtaMax(Double_t e)            { fTrackEtaMax = e;          }
      void                     SetTrackPtMin(Double_t pt)            { fTrackPtMin = pt;          }
      void                     SetTrackPtMax(Double_t pt)            { fTrackPtMax = pt;          }

    protected:
      void                     Process();
      void                     SlaveBegin();

      TString                  fPhotonColName;        //name of input lepton collection
      TString                  fTrackerTrackColName;  //name of input lepton collection
      TString                  fGsfTrackColName;      //name of input lepton collection
      Double_t                 fPhotonPtMin;          //minimum pt required   (def = 0 GeV)
      Double_t                 fPhotonPtMax;          //maximum pt required    (def = 5000 GeV)
      Double_t                 fPhotonEtaMin;         //minimum eta required   (def = -10)
      Double_t                 fPhotonEtaMax;         //maximum eta required   (def = +10) 
      Double_t                 fTrackPtMin;           //minimum pt required    (def = 0 GeV)
      Double_t                 fTrackPtMax;           //maximum pt required    (def = 5000 GeV)
      Double_t                 fTrackEtaMin;          //minimum eta required   (def = -10)
      Double_t                 fTrackEtaMax;          //maximum eta required   (def = +10) 
      const PhotonCol         *fPhotonCol;            //!pointer to collection 
      const TrackCol          *fTrackerTrackCol;      //!pointer to collection 
      const TrackCol          *fGsfTrackCol;          //!pointer to collection 
      TH1D                    *fNAccCounters;         //!acceptance histogram

      ClassDef(PhotonPlusIsoTrackSelMod,1) // Generic selection module
  };
}
#endif
