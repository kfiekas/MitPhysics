//------------------------------------------------------------------------------
// $Id: HWWKFactorList.h,v 1.01 2009/04/03 22:42:59 ceballos Exp $
//
// HWWKFactorList
//
// brief Resolution Map
// Basically just a TH1D with text I/O
//
// Authors: G. Gomez-Ceballos
//------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_HWWKFactorList_H
#define MITPHYSICS_MODS_HWWKFactorList_H

#include <TH1.h>

class HWWKfactorList : public TH1D {

 public:

  /// default constructor
  HWWKfactorList() : TH1D() {}
  
  /// create a map from text file mapfile
  HWWKfactorList(const char* name, const char* mapfile);
  
  /// create an empty map and initialize it 
  HWWKfactorList(const char* name, 
		  unsigned nbinspt, double minpt, double maxpt, double value);
  
  /// create a map from a 1d histogram
  HWWKfactorList(const TH1D& h) : TH1D(h) {}
 

  /// read text file
  bool ReadMapFile(const char* mapfile);

  /// write text file
  /// is not const because mapFile_ will be updated
  bool WriteMapFile(const char* mapfile);

  
  const char* GetMapFile() const {return mapFile_.c_str();}

  /// print this map
  friend std::ostream& operator<<(std::ostream& out, 
				  const HWWKfactorList& rm);


  // get alternative Kfactor
  inline double GetAlterKfactor() { return alternativeK_;};
  inline double GetAlterNNLOKfactor() { return alternativeNNLOK_;};
  
 private:
  static const unsigned lineSize_;
  std::string		mapFile_;

  double		alternativeK_;
  double		alternativeNNLOK_;

};
#endif
