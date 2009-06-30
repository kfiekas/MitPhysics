// $Id: ParticleFwd.h,v 1.1 2009/06/15 15:00:13 loizides Exp $

#ifndef MITPHYSICS_FAKEMODS_FAKEEVENTHEADERFWD_H
#define MITPHYSICS_FAKEMODS_FAKEEVENTHEADERFWD_H

#include "MitAna/DataCont/interface/CollectionFwd.h"
#include "MitAna/DataCont/interface/ArrayFwd.h"
#include "MitAna/DataCont/interface/ObjArrayFwd.h"

namespace mithep {
  class FakeEventHeader;
  typedef Collection<FakeEventHeader>          FakeEventHeaderCol;
  typedef Array<FakeEventHeader>               FakeEventHeaderArr;
  typedef ObjArray<FakeEventHeader>            FakeEventHeaderOArr;
}
#endif
