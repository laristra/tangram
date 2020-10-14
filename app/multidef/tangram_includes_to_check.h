/*
This file is part of the Ristra Tangram project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

// All the include files we want to check to make sure they won't
// cause multiply defined symbols. If some function does cause multiply
// defined symbol errors, the options are to
//
// 0. Make sure you have include guards so that multiple inclusions of
// the include file in a source file is ok
// 1. Move it to a .cc file
// 2. Make it inline
// 3. Make it static
// 4. Enclose it in its own namespace
//
// Option 3 and 4 will cause each translation unit (compiled source file) to
// have it's own copy of the function


#ifndef TANGRAM_MULTIDEF_H_
#define TANGRAM_MULTIDEF_H_

#include "tangram/driver/CellMatPoly.h"
#include "tangram/driver/driver.h"
#include "tangram/driver/write_to_gmv.h"

#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/intersect/split_rNd.h"

#include "tangram/reconstruct/cutting_distance_solver.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/reconstruct/nested_dissections.h"
#include "tangram/reconstruct/SLIC.h"
#include "tangram/reconstruct/VOF.h"

#include "tangram/support/bfgs.h"
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"

#endif
