#ifndef RHO_UTILITIES_H
#define RHO_UTILITIES_H

namespace mithep{

  struct RhoUtilities {
    enum RhoType {
      NONE = 0,
      MIT_RHO_VORONOI_LOW_ETA,
      MIT_RHO_VORONOI_HIGH_ETA,
      MIT_RHO_RANDOM_LOW_ETA,
      MIT_RHO_RANDOM_HIGH_ETA
    };
  };
}

#endif
