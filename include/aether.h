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

// Contains a reporting system for the model include verbose level and timing:
#include "report.h"

// not done
#include "inputs.h"

// Defines physical and conversion constants
#include "constants.h"

// Needs to be eliminated (but not in this branch!)
#include "sizes.h"

// These are functions that convert time between different systems.
#include "time_conversion.h"
// Functions that assist in the reading of files.
#include "file_input.h"
// A class for keeping track of indices (1d vectors w/time)
#include "indices.h"
// A set of routines to read in a number of different indices
#include "read_indices_files.h"
// A class for keeping track of all of the planetary characteristics
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
