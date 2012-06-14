#include "MitPhysics/Skim/interface/H4lLightFlavorSkim.h"

using namespace mithep;

ClassImp(mithep::H4lLightFlavorSkim)

H4lLightFlavorSkim::H4lLightFlavorSkim(const char *name, const char *title) :
   BaseH4lSkim(name,title)
{
}

//--------------------------------------------------------------------------------------------------
H4lLightFlavorSkim::~H4lLightFlavorSkim()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------      
void H4lLightFlavorSkim::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void H4lLightFlavorSkim::SlaveBegin()
{
  ReqBranch(fMuonName,       fMuons);
  ReqBranch(fElectronName,   fElectrons);
  ReqBranch(fPrimVtxName,    fPrimVerts);
  ReqBranch(fPfCandidateName,fPfCandidates);
  ReqBranch(fTracksName,     fTracks);
}

//--------------------------------------------------------------------------------------------------
void H4lLightFlavorSkim::SlaveTerminate()
{
  BaseH4lSkim::SlaveTerminate();
}

//--------------------------------------------------------------------------------------------------
void H4lLightFlavorSkim::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void H4lLightFlavorSkim::BeginRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lLightFlavorSkim::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lLightFlavorSkim::Process()
{
  fTotal++;
  // Load branches
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fPfCandidateName);

  // Set fBestPv
  SetBestPv();

  // Look for tight leptons
  vector<const Muon*>     tightMuons;
  vector<const Electron*> tightElectrons;

  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    if (!PassWwMuonSel(mu))
      continue;
    if (fabs(mu->Ip3dPVSignificance()) >= 3)
      continue;
    tightMuons.push_back(mu);
  }

  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);
    if (fabs(ele->Ip3dPVSignificance()) >= 3)
      continue;
    if (!PassElecTagSel(ele))
      continue;
    tightElectrons.push_back(ele);
  }

  // Look for an additional opposite-flavor lepton
  vector<const Muon*> looseMuons;
  vector<const Electron*> looseElectrons;

  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    if (!PassMuonPreselNoIp(mu))
      continue;
    looseMuons.push_back(mu);
  }
  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);
    if (!PassElecPreselNoIp(ele))
      continue;
    looseElectrons.push_back(ele);
  }
  
  if ((tightMuons.size()     > 0 && looseElectrons.size() > 0) ||
      (tightElectrons.size() > 0 && looseMuons.size()     > 0)   ) {
    fSelected++;
  }
  else
    SkipEvent();
}
