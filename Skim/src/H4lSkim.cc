#include "MitPhysics/Skim/interface/H4lSkim.h"

using namespace mithep;

ClassImp(mithep::H4lSkim)
//--------------------------------------------------------------------------------------------------
H4lSkim::H4lSkim(const char *name, const char *title):
BaseH4lSkim(name,title)
{
}

//--------------------------------------------------------------------------------------------------
H4lSkim::~H4lSkim()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------      
void H4lSkim::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::SlaveBegin()
{
  ReqBranch(fMuonName,       fMuons);
  ReqBranch(fElectronName,   fElectrons);
  ReqBranch(fPrimVtxName,    fPrimVerts);
  ReqBranch(fPfCandidateName,fPfCandidates);  
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::SlaveTerminate()
{
  BaseH4lSkim::SlaveTerminate();
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::BeginRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lSkim::Process()
{
  fTotal++;
  // Load branches
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fPfCandidateName);

  // Set fBestPv
  SetBestPv();

  int nLepPt5=0;
  int nLepPt10=0;

  // loop through muons
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i); 
    if (!mu->IsTrackerMuon())
      continue; 
    if (fabs(mu->Eta()) > 2.5)
      continue; 
    if (mu->Pt() >5)
      nLepPt5++; 
    if (mu->Pt() >10)
      nLepPt10++; 
  }

  // loop through electrons.
  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);  
    if (fabs(ele->Eta()) > 2.5)
      continue;
    if (ele->Pt() >5)
      nLepPt5++;
    if (ele->Pt() >10)
      nLepPt10++;
  }

  if ( !(nLepPt5>=4 && nLepPt10>=2)) SkipEvent();
  else fSelected++;
}
