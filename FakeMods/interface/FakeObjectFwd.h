// $Id: ParticleFwd.h,v 1.1 2009/06/15 15:00:13 loizides Exp $

#ifndef MITPHYSICS_FAKEMODS_FAKEOBJECTFWD_H
#define MITPHYSICS_FAKEMODS_FAKEOBJECTFWD_H

#include "MitAna/DataCont/interface/CollectionFwd.h"
#include "MitAna/DataCont/interface/ArrayFwd.h"
#include "MitAna/DataCont/interface/ObjArrayFwd.h"

namespace mithep {
  class FakeObject;
  typedef Collection<FakeObject>          FakeObjectCol;
  typedef Array<FakeObject>               FakeObjectArr;
  typedef ObjArray<FakeObject>            FakeObjectOArr;
}
#endif
