// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_AETHER_H_
#define INCLUDE_AETHER_H_

/**************************************************************
 * 
 * This is the main header file for the aether model. It includes
 * all of the header files for the model.
 *
 **************************************************************/

// The armadillo library is to allow the use of 3d cubes and other
// array types, with array math built in. This eliminates loops!
#include <armadillo>
using namespace arma;

// Contains all information about time in the code and wall time:
#include "times.h"

// not done
#include "report.h"

// not done
#include "inputs.h"

// Defines physical and conversion constants
#include "constants.h"

// Needs to be eliminated (but not in this branch!)
#include "sizes.h"

// not done
#include "time_conversion.h"
// not done
#include "file_input.h"
// not done
#include "indices.h"
// not done
#include "read_f107_file.h"
// not done
#include "planets.h"
// not done
#include "grid.h"

// Contains the neutral states and derived quantities
#include "neutrals.h"
// not done
#include "ions.h"

// not done
#include "bfield.h"

// Defines the Extreme Ultraviolet radiation above the atmosphere
#include "euv.h"
// not done
#include "calc_euv.h"
// not done
#include "chemistry.h"
// not done
#include "output.h"
// not done
#include "advance.h"

// not done
#include "solvers.h"
// not done
#include "transform.h"

#endif  // INCLUDE_AETHER_H_
