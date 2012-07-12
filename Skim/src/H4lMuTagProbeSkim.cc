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

  pfNoPileUpflag.clear();
  UInt_t pfnopu_size = makePFnoPUArray();
  assert(pfnopu_size == fPfCandidates->GetEntries());

  SetBestPv();

  // Find a tag
  bool hasTag = false;
  int nProbes = 0;
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    if(mu->Pt() > 5 && fabs(mu->Eta()) < 2.4)
      nProbes++;
    if(mu->Pt() < 20)
      continue;
    if(fabs(mu->Eta()) > 2.4)
      continue;
    if(!muon2012CutBasedIDTightForTagProbe(mu))
      continue;
    hasTag = true;
  }

  if(hasTag && nProbes>1)
    fSelected++;
  else
    SkipEvent();
}
