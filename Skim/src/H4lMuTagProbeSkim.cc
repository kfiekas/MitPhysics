#include "MitPhysics/Skim/interface/H4lMuTagProbeSkim.h"

using namespace mithep;

ClassImp(mithep::H4lMuTagProbeSkim)

H4lMuTagProbeSkim::H4lMuTagProbeSkim(const char *name, const char *title):
  BaseH4lSkim(name,title)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
H4lMuTagProbeSkim::~H4lMuTagProbeSkim()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------
void H4lMuTagProbeSkim::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void H4lMuTagProbeSkim::SlaveBegin()
{

  ReqBranch(fMuonName,            fMuons);
  ReqBranch(fElectronName,        fElectrons);
  ReqBranch(fPrimVtxName,         fPrimVerts);
  ReqBranch(fPfCandidateName,     fPfCandidates);
  ReqBranch(fTracksName,          fTracks);
}

//--------------------------------------------------------------------------------------------------
void H4lMuTagProbeSkim::SlaveTerminate()
{
  BaseH4lSkim::SlaveTerminate();
}

//--------------------------------------------------------------------------------------------------
void H4lMuTagProbeSkim::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void H4lMuTagProbeSkim::BeginRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lMuTagProbeSkim::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lMuTagProbeSkim::Process()
{
  fTotal++;

  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fPfCandidateName);
  LoadBranch(fTracksName);

  SetBestPv();

  // Find a tag
  bool hasTag = false;
  int nProbes = 0;
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    if(mu->Pt() > 5)
      nProbes++;
    if(!muon2012CutBasedIDTight(mu))
      continue;
    hasTag = true;
  }

  for(UInt_t i=0; i<fTracks->GetEntries(); ++i) {
    const mithep::Track *track = fTracks->At(i);
    
    // Check that the track is not associated with a muon.                                                                                                                                                      
    Bool_t isMuon = kFALSE;
    for(UInt_t j=0; j<fMuons->GetEntries(); ++j) {
      if(track == (fMuons->At(j)->TrackerTrk())) isMuon = kTRUE;
    }
    if(isMuon)
      continue;
    if(track->Pt() < 20)
      continue;
    if(fabs(track->Eta()) > 2.4)
      continue;

    nProbes++;
  }


  if(hasTag && nProbes>1)
    fSelected++;
  else
    SkipEvent();
}
