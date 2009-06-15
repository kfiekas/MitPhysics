//--------------------------------------------------------------------------------------------------
// $Id: JetPlusIsoTrackSelMod.h,v 1.2 2009/03/23 22:17:04 loizides Exp $
//
// JetPlusIsoTrackSelMod
// 
// This module selects events containing one lepton and one isolated track
// Both gsfTracks and tracker tracks are considered.
//
// Authors: S. Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_SELMODS_JETPLUSISOTRACKSELMOD_H
#define MITPHYSICS_SELMODS_JETPLUSISOTRACKSELMOD_H

#include "MitAna/TreeMod/interface/BaseSelMod.h" 
#include "MitAna/DataTree/interface/JetFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
class TH1D;

namespace mithep 
{
  class JetPlusIsoTrackSelMod : public BaseSelMod
  {
    public:
      JetPlusIsoTrackSelMod(const char *name="JetPlusIsoTrackSelMod", 
                            const char *title="Jet plus isolated track selection module");

      const char              *GetJetColName()          const { return fJetColName;               }
      const char              *GetTrackerTrackColName() const { return fTrackerTrackColName;      }
      const char              *GetGsfTrackColName()     const { return fGsfTrackColName;          }
      Double_t                 GetJetEtaMin()           const { return fJetEtaMin;                }
      Double_t                 GetJetEtaMax()           const { return fJetEtaMax;                }
      Double_t                 GetJetPtMin()            const { return fJetPtMin;                 }
      Double_t                 GetJetPtMax()            const { return fJetPtMax;                 }
      Double_t                 GetTrackEtaMin()         const { return fTrackEtaMin;              }
      Double_t                 GetTrackEtaMax()         const { return fTrackEtaMax;              }
      Double_t                 GetTrackPtMin()          const { return fTrackPtMin;               }
      Double_t                 GetTrackPtMax()          const { return fTrackPtMax;               }
      void                     SetJetColName(const char *n)          { fJetColName=n;             }
      void                     SetTrackerTrackColName(const char *n) { fTrackerTrackColName=n;    }
      void                     SetGsfTrackColName(const char *n)     { fGsfTrackColName=n;        }
      void                     SetJetEtaMin(Double_t e)              { fJetEtaMin = e;            }
      void                     SetJetEtaMax(Double_t e)              { fJetEtaMax = e;            }
      void                     SetJetPtMin(Double_t pt)              { fJetPtMin = pt;            }
      void                     SetJetPtMax(Double_t pt)              { fJetPtMax = pt;            }
      void                     SetTrackEtaMin(Double_t e)            { fTrackEtaMin = e;          }
      void                     SetTrackEtaMax(Double_t e)            { fTrackEtaMax = e;          }
      void                     SetTrackPtMin(Double_t pt)            { fTrackPtMin = pt;          }
      void                     SetTrackPtMax(Double_t pt)            { fTrackPtMax = pt;          }

    protected:
      void                     Process();
      void                     SlaveBegin();

      TString                  fJetColName;           //name of input lepton collection
      TString                  fTrackerTrackColName;  //name of input lepton collection
      TString                  fGsfTrackColName;      //name of input lepton collection
      Double_t                 fJetPtMin;             //minimum pt required   (def = 0 GeV)
      Double_t                 fJetPtMax;             //maximum pt required    (def = 5000 GeV)
      Double_t                 fJetEtaMin;            //minimum eta required   (def = -10)
      Double_t                 fJetEtaMax;            //maximum eta required   (def = +10) 
      Double_t                 fTrackPtMin;           //minimum pt required    (def = 0 GeV)
      Double_t                 fTrackPtMax;           //maximum pt required    (def = 5000 GeV)
      Double_t                 fTrackEtaMin;          //minimum eta required   (def = -10)
      Double_t                 fTrackEtaMax;          //maximum eta required   (def = +10) 
      const JetCol            *fJetCol;               //!pointer to collection 
      const TrackCol          *fTrackerTrackCol;      //!pointer to collection 
      const TrackCol          *fGsfTrackCol;          //!pointer to collection 
      TH1D                    *fNAccCounters;         //!acceptance histogram

      ClassDef(JetPlusIsoTrackSelMod,1) // Jet plus isolated track selection module
  };
}
#endif
