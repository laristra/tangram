/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef INTERSECT_SPLIT_RnD_H
#define INTERSECT_SPLIT_RnD_H

#include <type_traits>

#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"

// Dimension agnostic version of split_r2d, split_r3d functionality
// If we ever have it for 1D, we have to put a nested conditional
//
// template<int D>
// using ClipRnD =
//  typename std::conditional<D==1,
//                            ClipR1D,
//                            typename std::conditional<D==2,
//                                                      ClipR2D,
//                                                      ClipR3D>::type>
//                           >::type;


namespace Tangram {

template<int D, class CoordSys = Wonton::DefaultCoordSys>
using SplitRnD = typename std::conditional<D==2, SplitR2D<CoordSys>, SplitR3D>::type;
template<int D, class CoordSys = Wonton::DefaultCoordSys>
using ClipRnD = typename std::conditional<D==2, ClipR2D<CoordSys>, ClipR3D>::type;

}

#endif
