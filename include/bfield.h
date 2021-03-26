// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_BFIELD_H_
#define INCLUDE_BFIELD_H_

/**************************************************************
 * 
 * bfield.h:
 *
 **************************************************************/

struct bfield_info_type {
  float b[3];
  float lon;
  float lat;
};

/**************************************************************
 * Returns the location (lon/lat) of the magnetic pole.
 *
 * @param IsNorth Set to 1 if you want the north magnetic pole. Set
 *                    to 0 if you want the south.
 * @param planet Need to pass Planet class, so code can get info 
 *               about b-field
 * @param input Need to pass Input class, so code can get info
 *              about how user has configured things
 * @param report Need to pass Report class, so reporting can occur  
 **************************************************************/

fvec get_magnetic_pole(int IsNorth,
		       Planets planet,
		       Inputs input,
		       Report &report);

bfield_info_type get_bfield(float lon,
                            float lat,
                            float alt,
                            Planets planet,
                            Inputs input,
                            Report &report);

bfield_info_type get_dipole(float lon,
                            float lat,
                            float alt,
                            Planets planet,
                            Inputs input,
                            Report &report);

#endif // INCLUDE_BFIELD_H_
