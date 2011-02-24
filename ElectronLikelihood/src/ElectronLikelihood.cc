#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include <stdlib.h>
#include <iostream>


ElectronLikelihood::ElectronLikelihood (TDirectory *EBlt15dir, TDirectory *EElt15dir,
                                        TDirectory *EBgt15dir, TDirectory *EEgt15dir,
					LikelihoodSwitches eleIDSwitches,
					std::string signalWeightSplitting,
					std::string backgroundWeightSplitting,
					bool splitSignalPdfs,
					bool splitBackgroundPdfs) :
  _EBlt15lh (new LikelihoodPdfProduct ("electronID_EB_ptLt15_likelihood",0,0)) ,
  _EElt15lh (new LikelihoodPdfProduct ("electronID_EE_ptLt15_likelihood",1,0)) ,
  _EBgt15lh (new LikelihoodPdfProduct ("electronID_EB_ptGt15_likelihood",0,1)) ,
  _EEgt15lh (new LikelihoodPdfProduct ("electronID_EE_ptGt15_likelihood",1,1)) ,
  m_eleIDSwitches (eleIDSwitches) ,
  m_signalWeightSplitting (signalWeightSplitting), 
  m_backgroundWeightSplitting (backgroundWeightSplitting),
  m_splitSignalPdfs (splitSignalPdfs), 
  m_splitBackgroundPdfs (splitBackgroundPdfs)  
{
  Setup (EBlt15dir, EElt15dir,
         EBgt15dir, EEgt15dir,
	 signalWeightSplitting, backgroundWeightSplitting,
	 splitSignalPdfs, splitBackgroundPdfs) ;
}



// --------------------------------------------------------



ElectronLikelihood::~ElectronLikelihood () {
  delete _EBlt15lh ;
  delete _EElt15lh ;
  delete _EBgt15lh ;
  delete _EEgt15lh ;
}



// --------------------------------------------------------


void 
ElectronLikelihood::Setup (TDirectory *EBlt15dir, TDirectory *EElt15dir,
                           TDirectory *EBgt15dir, TDirectory *EEgt15dir,
			   std::string signalWeightSplitting,
			   std::string backgroundWeightSplitting,
			   bool splitSignalPdfs,
			   bool splitBackgroundPdfs) 
{

  // ECAL BARREL LIKELIHOOD - Pt < 15 GeV region
  _EBlt15lh->initFromFile (EBlt15dir) ;

  _EBlt15lh->addSpecies ("electrons", 1.0) ;
  _EBlt15lh->addSpecies ("hadrons",  1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EBlt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EBlt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EBlt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EBlt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EBlt15lh->setSplitFrac ("electrons", "class0", 1.0) ;
    _EBlt15lh->setSplitFrac ("electrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EBlt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EBlt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EBlt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EBlt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EBlt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EBlt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EBlt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EBlt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EBlt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EBlt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EBlt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EBlt15lh->setSplitFrac ("hadrons", "class0", 1.0) ;
    _EBlt15lh->setSplitFrac ("hadrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EBlt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EBlt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EBlt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EBlt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EBlt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EBlt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EBlt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;



  // ECAL BARREL LIKELIHOOD - Pt >= 15 GeV region
  _EBgt15lh->initFromFile (EBgt15dir) ;

  _EBgt15lh->addSpecies ("electrons", 1.0 ) ;  
  _EBgt15lh->addSpecies ("hadrons",1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EBgt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EBgt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EBgt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EBgt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EBgt15lh->setSplitFrac ("electrons", "class0", 1.0);
    _EBgt15lh->setSplitFrac ("electrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EBgt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EBgt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EBgt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EBgt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EBgt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EBgt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EBgt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EBgt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EBgt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EBgt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EBgt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EBgt15lh->setSplitFrac ("hadrons", "class0", 1.0);
    _EBgt15lh->setSplitFrac ("hadrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EBgt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EBgt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EBgt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EBgt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EBgt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EBgt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EBgt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;


  // ECAL ENDCAP LIKELIHOOD - Pt < 15 GeV
  _EElt15lh->initFromFile (EElt15dir) ;

  _EElt15lh->addSpecies ("electrons", 1.0 ) ;
  _EElt15lh->addSpecies ("hadrons",1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EElt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EElt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EElt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EElt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EElt15lh->setSplitFrac ("electrons", "class0", 1.0);
    _EElt15lh->setSplitFrac ("electrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now";
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EElt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EElt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EElt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EElt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EElt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EElt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EElt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EElt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EElt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EElt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EElt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EElt15lh->setSplitFrac ("hadrons", "class0", 1.0);
    _EElt15lh->setSplitFrac ("hadrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EElt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EElt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EElt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EElt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EElt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EElt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EElt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;

  // ECAL ENDCAP LIKELIHOOD - Pt >= 15 GeV
  _EEgt15lh->initFromFile (EEgt15dir) ;

  _EEgt15lh->addSpecies ("electrons", 1.0 ) ;
  _EEgt15lh->addSpecies ("hadrons",1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EEgt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EEgt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EEgt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EEgt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EEgt15lh->setSplitFrac ("electrons", "class0", 1.0) ;
    _EEgt15lh->setSplitFrac ("electrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EEgt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EEgt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EEgt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EEgt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EEgt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EEgt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EEgt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EEgt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EEgt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EEgt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EEgt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EEgt15lh->setSplitFrac ("hadrons", "class0", 1.0) ;
    _EEgt15lh->setSplitFrac ("hadrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EEgt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EEgt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EEgt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EEgt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EEgt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EEgt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EEgt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;

}



// --------------------------------------------------------




// --------------------------------------------------------



float 
ElectronLikelihood::result (const LikelihoodMeasurements electron) const 
{

  //=======================================================
  // used classification:
  // nbrem clusters = 0         =>  0
  // nbrem clusters >= 1        =>  1
  //=======================================================

  std::vector<float> measurements ;
  if(m_eleIDSwitches.m_useDeltaPhi) measurements.push_back( electron.deltaPhi );
  if(m_eleIDSwitches.m_useDeltaEta) measurements.push_back( electron.deltaEta );
  if(m_eleIDSwitches.m_useEoverP) measurements.push_back( electron.eSuperClusterOverP );
  if(m_eleIDSwitches.m_useHoverE) measurements.push_back( electron.hadronicOverEm );
  if(m_eleIDSwitches.m_useSigmaEtaEta) measurements.push_back( electron.sigmaIEtaIEta );
  if(m_eleIDSwitches.m_useSigmaPhiPhi) measurements.push_back( electron.sigmaIPhiIPhi );
  if(m_eleIDSwitches.m_useFBrem) measurements.push_back( electron.fBrem );

  // Split using only the 1 / >1 cluster
  int nBremClusters=electron.nBremClusters;
  int bitVal = (nBremClusters==0) ? 0 : 1 ;
  
  char className[20] ;
  if(m_signalWeightSplitting.compare("class")==0) {
    sprintf (className,"class%d",bitVal);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  int subdet = electron.subdet;
  float thisPt =  electron.pt;

  if (subdet==0 && thisPt<15.)
    return _EBlt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==0 && thisPt>=15.)
    return _EBgt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==1 && thisPt<15.)
    return _EElt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==1 && thisPt>=15.)
    return _EEgt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else return -999. ;
}

