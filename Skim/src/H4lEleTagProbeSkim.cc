#include "MitPhysics/Skim/interface/H4lEleTagProbeSkim.h"

using namespace mithep;

ClassImp(mithep::H4lEleTagProbeSkim)

H4lEleTagProbeSkim::H4lEleTagProbeSkim(const char *name, const char *title):
  BaseH4lSkim(name,title)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
H4lEleTagProbeSkim::~H4lEleTagProbeSkim()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------      
void H4lEleTagProbeSkim::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void H4lEleTagProbeSkim::SlaveBegin()
{
  ReqBranch(fMuonName,            fMuons);
  ReqBranch(fElectronName,        fElectrons);
  ReqBranch(fPrimVtxName,         fPrimVerts);
  ReqBranch(fPfCandidateName,     fPfCandidates);
  
}

//--------------------------------------------------------------------------------------------------
void H4lEleTagProbeSkim::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void H4lEleTagProbeSkim::Terminate()
{
  BaseH4lSkim::SlaveTerminate();
}

//--------------------------------------------------------------------------------------------------
void H4lEleTagProbeSkim::BeginRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lEleTagProbeSkim::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void H4lEleTagProbeSkim::Process()
{
  fTotal++;

  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fPfCandidateName);

  SetBestPv();

  // Find a tag
  bool hasTag = false;
  int nProbes = 0;
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);
    if(ele->SCluster()->Et() > 5)
      nProbes++;

    if(!electron2012CutBasedIDMedium(ele))
      continue;
    hasTag = true;
  }

  if(hasTag && nProbes>1)
    fSelected++;
  else
    SkipEvent();
}
