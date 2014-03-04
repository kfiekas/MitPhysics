#include <TMath.h> 
#include <TFile.h>
#include <TRandom3.h>
#include <TSystem.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitPhysics/Utils/interface/MVATools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/StableData.h"

ClassImp(mithep::MVATools)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MVATools::MVATools():
  fReaderEndcap(0),
  fReaderBarrel(0),
  fMVAType     (MVATools::kNone)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MVATools::InitializeMVA(MVATools::IdMVAType type) {
  
  fMVAType = type;

  if (fReaderEndcap)
    delete fReaderEndcap;  
  if (fReaderBarrel)
    delete fReaderBarrel; 

  // no MVA needed if none requested  
  if (type == kNone)
    return;
  
  fReaderEndcap = new TMVA::Reader( "!Color:!Silent:Error" );       
  fReaderBarrel = new TMVA::Reader( "!Color:!Silent:Error" );   

  TString BarrelWeights;
  TString EndcapWeights;

  varNames.resize(0);
  
  switch (type) {
    
  case k2011IdMVA_HZg:

    EndcapWeights = gSystem->Getenv("MIT_DATA")+TString("/PhotonId_lowPt_EE_BDTG.weights.xml");
    BarrelWeights = gSystem->Getenv("MIT_DATA")+TString("/PhotonId_lowPt_EB_BDTG.weights.xml");

    varNames.push_back("sigieie");
    varNames.push_back("covieip");
    varNames.push_back("s4r"    );
    varNames.push_back("r9"     );
    varNames.push_back("sigeta" );
    varNames.push_back("sigphi" );
    varNames.push_back("pfgios" );
    varNames.push_back("pfciso" );
    varNames.push_back("rho"    );
    varNames.push_back("sceta"  );
    varNames.push_back("rawe"   );

    mvaVars.resize(varNames.size());

    std::cout<<"  Adding stuff here.... "<<std::endl;

    for (unsigned int iV = 0; iV < mvaVars.size(); ++iV) {      
      mvaVarMapEB.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
      mvaVarMapEE.insert(  std::pair<std::string,unsigned int>(varNames[iV], iV) );
    }
    
    std::cout<<"  ... done "<<std::endl;

    break;

  case k2011IdMVA:

    EndcapWeights =  (gSystem->Getenv("MIT_DATA")+
		      TString("/TMVAClassificationPhotonID_Endcap_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
		      TString("weights.xml"));
    BarrelWeights =  (gSystem->Getenv("MIT_DATA")+
		      TString("/TMVAClassificationPhotonID_Barrel_PassPreSel_Variable_10_BDTnCuts2000_BDT.")+
		      TString("weights.xml"));

    // set up the variable names
    mvaVars.resize(12);
    varNames.push_back("HoE"        );
    varNames.push_back("covIEtaIEta");
    varNames.push_back("tIso1abs"   );
    varNames.push_back("tIso3abs"   );
    varNames.push_back("tIso2abs"   );
    varNames.push_back("R9"         );
    varNames.push_back("absIsoEcal" );
    varNames.push_back("absIsoHcal" );
    varNames.push_back("NVertexes"  );
    varNames.push_back("ScEta"      );
    varNames.push_back("EtaWidth"   );
    varNames.push_back("PhiWidth"   );

    for (unsigned int iV = 0; iV < mvaVars.size(); ++iV) {      
      mvaVarMapEB.insert(std::pair<std::string,unsigned int>(varNames[iV],iV));
      mvaVarMapEE.insert(std::pair<std::string,unsigned int>(varNames[iV],iV));
    }
    
    break;

  case k2012IdMVA:
  case k2012IdMVA_globe:
    
    EndcapWeights =      (gSystem->Getenv("MIT_DATA")+
			  TString("/2012ICHEP_PhotonID_Endcap_BDT.weights_PSCorr.xml"));
    BarrelWeights =      (gSystem->Getenv("MIT_DATA")+
			  TString("/2012ICHEP_PhotonID_Barrel_BDT.weights.xml"));

    mvaVars.resize(12);
    varNames.push_back("ph.r9"			      );
    varNames.push_back("ph.sigietaieta"		      );
    varNames.push_back("ph.scetawidth"		      );
    varNames.push_back("ph.scphiwidth"		      );
    varNames.push_back("ph.idmva_CoviEtaiPhi"	      );
    varNames.push_back("ph.idmva_s4ratio"	      );
    varNames.push_back("ph.idmva_GammaIso"	      );
    varNames.push_back("ph.idmva_ChargedIso_selvtx"   );
    varNames.push_back("ph.idmva_ChargedIso_worstvtx" );
    varNames.push_back("ph.sceta"		      );
    varNames.push_back("rho"			      );
    varNames.push_back("ph.idmva_PsEffWidthSigmaRR"   );

    for( unsigned int iV = 0; iV < mvaVars.size() - 1; ++iV) {
      mvaVarMapEB.insert(std::pair<std::string,unsigned int>(varNames[iV],iV));
      mvaVarMapEE.insert(std::pair<std::string,unsigned int>(varNames[iV],iV));
    }
    
    // pre-shower only used for Endcaps
    mvaVarMapEE.insert( std::pair<std::string,unsigned int> (varNames[mvaVars.size()-1], mvaVars.size()-1));

    break;

  case k2013FinalIdMVA:
    
    EndcapWeights =      (gSystem->Getenv("MIT_DATA")+
			  TString("/2013FinalPaper_PhotonID_Endcap_BDT_TrainRangePT15.weights.xml"));
    BarrelWeights =      (gSystem->Getenv("MIT_DATA")+
			  TString("/2013FinalPaper_PhotonID_Barrel_BDT_TrainRangePT15.weights.xml"));

    mvaVars.resize(13);
    varNames.push_back("ph.scrawe"		      );
    varNames.push_back("ph.r9"			      );
    varNames.push_back("ph.sigietaieta"		      );
    varNames.push_back("ph.scetawidth"		      );
    varNames.push_back("ph.scphiwidth"		      );
    varNames.push_back("ph.idmva_CoviEtaiPhi"	      );
    varNames.push_back("ph.idmva_s4ratio"	      );
    varNames.push_back("ph.idmva_GammaIso"	      );
    varNames.push_back("ph.idmva_ChargedIso_selvtx"   );
    varNames.push_back("ph.idmva_ChargedIso_worstvtx" );
    varNames.push_back("ph.sceta"		      );
    varNames.push_back("rho"			      );
    varNames.push_back("ph.idmva_PsEffWidthSigmaRR"   );

    for (unsigned int iV = 0; iV < mvaVars.size() - 1; ++iV) {
      mvaVarMapEB.insert(std::pair<std::string,unsigned int>(varNames[iV],iV));
      mvaVarMapEE.insert(std::pair<std::string,unsigned int>(varNames[iV],iV));
    }
    
    // pre-shower only used for Endcaps
    mvaVarMapEE.insert( std::pair<std::string,unsigned int> (varNames[mvaVars.size()-1], mvaVars.size()-1));

    break;
    
  default:
    // no variables... better never called..
    std::cerr<<" MVATools: Trying to initialize with unknown type."<<std::endl;
    break;
  }


  // looping over both maps and adding Vars to BDT readers
  for (unsigned int iV = 0; iV < varNames.size(); ++iV) {
    std::map<std::string,unsigned int>::const_iterator it = mvaVarMapEB.find(varNames[iV]);
    if (it != mvaVarMapEB.end() )
      fReaderBarrel -> AddVariable((it->first).c_str(), &(mvaVars[it->second]));
    it = mvaVarMapEE.find(varNames[iV]);
    if (it != mvaVarMapEE.end() )
      fReaderEndcap -> AddVariable((it->first).c_str(), &(mvaVars[it->second]));
  }

  fReaderEndcap->BookMVA("BDT method",EndcapWeights);
  fReaderBarrel->BookMVA("BDT method",BarrelWeights);
  
  assert(fReaderEndcap);
  assert(fReaderBarrel);
}


//--------------------------------------------------------------------------------------------------
Double_t MVATools::GetMVAbdtValue(const Photon* p, const Vertex* vtx, const TrackCol* trackCol,
				  const VertexCol* vtxCol, Double_t _tRho, const PFCandidateCol *fPFCands,
				  const ElectronCol* els, Bool_t applyElectronVeto)
{
  // if there's no reader, or the type is kNone, return the default values of -99.
  if ((!fReaderBarrel || !fReaderEndcap) || fMVAType == kNone)
    return -99.;

  // we compute the variable names... make sure no confilcts when adding new variables...
  
  // check if it's a Barrel or EE photon
  bool isBarrel = ( p->SCluster()->AbsEta() < 1.5 );

  std::map<std::string,unsigned int>* theVars = ( isBarrel ? &mvaVarMapEB : &mvaVarMapEE );  
  
  // loop over all the variables in the map... and keep count (to make sure we have filled all variables)
  unsigned int varCounter = 0;
  for (std::map<std::string,unsigned int>::const_iterator iV = theVars->begin(); iV != theVars->end(); ++iV) {
    
    TString theVarName  = TString(iV->first);  
    float* theVarValue  = &(mvaVars[iV->second]);           // pointer to the variable...
    
    if( 
       !theVarName.CompareTo("HoE") 
       ) {
      (*theVarValue) = p->HadOverEm();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("covIEtaIEta") || !theVarName.CompareTo("ph.sigietaieta") || !theVarName.CompareTo("sigieie") 
	       ) {
      (*theVarValue) = p->CoviEtaiEta();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("R9") || !theVarName.CompareTo("ph.r9") || !theVarName.CompareTo("r9")
	       ) {
      (*theVarValue) = p->R9();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ScEta") || !theVarName.CompareTo("ph.sceta") || !theVarName.CompareTo("sceta")
	       ) {
      (*theVarValue) = p->SCluster()->Eta();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("rho")
	       ) {
      (*theVarValue) = _tRho;
       varCounter++;
    } else if (
	       !theVarName.CompareTo("tIso1abs")
	       ) {      
      double _ecalIso3 = p->EcalRecHitIsoDr03();
      double _hcalIso4 = p->HcalTowerSumEtDr04();  
      double _trackIso1 = IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, NULL, NULL, (!applyElectronVeto ? els : NULL) );//Question Ming:whyfPV->At(0) instead of selected vertex using ranking method?
      (*theVarValue) = _ecalIso3+_hcalIso4+_trackIso1 - 0.17*_tRho;
      varCounter++;
    } else if (
	       !theVarName.CompareTo("tIso2abs")
	       ) {
      double _ecalIso4 = p->EcalRecHitIsoDr04();
      double _hcalIso4 = p->HcalTowerSumEtDr04();
      unsigned int wVtxInd = 0;
      double _trackIso2 = IsolationTools::CiCTrackIsolation(p,vtx, 0.4, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, &wVtxInd,vtxCol, (!applyElectronVeto ? els : NULL) );
      (*theVarValue) = _ecalIso4+_hcalIso4+_trackIso2 - 0.52*_tRho;
      varCounter++;
    } else if (
	       !theVarName.CompareTo("tIso3abs")
	       ) {
      (*theVarValue) = IsolationTools::CiCTrackIsolation(p,vtx, 0.3, 0.02, 0.0, 0.0, 0.1, 1.0, trackCol, NULL, NULL, (!applyElectronVeto ? els : NULL) );
      varCounter++;
    } else if (
	       !theVarName.CompareTo("absIsoEcal")
	       ) {
      double _ecalIso3 = p->EcalRecHitIsoDr03();
      (*theVarValue) =  (_ecalIso3-0.17*_tRho);
      varCounter++;
    } else if (
	       !theVarName.CompareTo("absIsoHcal")
	       ) {
      double _hcalIso4 = p->HcalTowerSumEtDr04();
      (*theVarValue) = (_hcalIso4-0.17*_tRho);
      varCounter++;
    } else if (
	       !theVarName.CompareTo("NVertexes")
	       ) {

      (*theVarValue) = vtxCol->GetEntries();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("EtaWidth") || !theVarName.CompareTo("ph.scetawidth") || !theVarName.CompareTo("sigeta")
	       ) {
      (*theVarValue) = p->EtaWidth();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("PhiWidth") || !theVarName.CompareTo("ph.scphiwidth") || !theVarName.CompareTo("sigphi")
	       ) {
      (*theVarValue) = p->PhiWidth();
      varCounter++;

    } else if (
	       !theVarName.CompareTo("ph.idmva_CoviEtaiPhi") || !theVarName.CompareTo("covieip")
	       ) {
      (*theVarValue) = p->SCluster()->Seed()->CoviEtaiPhi();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_s4ratio") || !theVarName.CompareTo("s4r") 
	       ) {
      (*theVarValue) = p->S4Ratio();
      //(*theVarValue) = p->SCluster()->Seed()->E2x2()/p->E55();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_GammaIso") || !theVarName.CompareTo("pfgiso") 
	       ) {
      (*theVarValue) = IsolationTools::PFGammaIsolation(p,0.3,0,fPFCands);
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_ChargedIso_selvtx") || !theVarName.CompareTo("pfciso") 
	       ) {
      (*theVarValue) = IsolationTools::PFChargedIsolation(p,vtx,0.3,0,fPFCands);
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_ChargedIso_worstvtx")
	       ) {
      unsigned int wVtxInd = 0;
      (*theVarValue) = IsolationTools::PFChargedIsolation(p,vtx,0.3,0,fPFCands,&wVtxInd,vtxCol); 
      varCounter++;
    } else if (
	       !theVarName.CompareTo("ph.idmva_PsEffWidthSigmaRR")
	       ) {
      (*theVarValue) = p->EffSigmaRR();
      varCounter++;
    } else if (
	       !theVarName.CompareTo("rawe")|| !theVarName.CompareTo("ph.scrawe") 
	       ) {
      (*theVarValue) = p->SCluster()->RawEnergy();
      varCounter++;
    } else {
      // a variable is not know... copmplain!
      std::cerr<<" ERROR: MVA Evaluation called with unknown variable name >"<<theVarName<<">."<<std::endl;
    }
  }
  
  // now all the variables should be filled... check!
  if( varCounter != theVars->size() )
    std::cerr<<" ERROR: MVA Evaludation called and not all variables are filled."<<std::endl;

  // we're ready to compute the MVA value
  TMVA::Reader* reader = NULL;
  if (isBarrel)
    reader = fReaderBarrel;
  else
    reader = fReaderEndcap;
  
  assert(reader); 
  
  return (reader->EvaluateMVA("BDT method"));
}
