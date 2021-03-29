// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#ifndef INCLUDE_BFIELD_H_
#define INCLUDE_BFIELD_H_

/**************************************************************
 * 
 * bfield.h:
 *
 *   Functions that are associated with the magnetic field of the
 * specific planet. These characteristics are set in the Planets
 * class (through the planetary CSV file). Currently, you can set
 * the magnetic field type to be dipolar, which allows for a
 * tilted and offset dipole field.
 *
 **************************************************************/

/**************************************************************
 * bfield_info_type is a structure that allows the specification
 * of the three components of the magnetic field, as well as the
 * latitude and longitude in the magnetic coordinate system.
 **************************************************************/

struct bfield_info_type {
  float b[3];   ///< three components of the magnetic field (lon, lat, rad)
  float lon;    ///< magnetic longitude (in radians)
  float lat;    ///< magnetic (invariant) latitude (in radians)
};

/**************************************************************
 * Returns the location (lon/lat) of the magnetic pole.
 *
 * \param IsNorth Set to 1 if you want the north magnetic pole. Set
 *                    to 0 if you want the south.
 * \param planet Need to pass Planet class, so code can get info 
 *               about b-field
 * \param input Need to pass Input class, so code can get info
 *              about how user has configured things
 * \param report Need to pass Report class, so reporting can occur  
 **************************************************************/

fvec get_magnetic_pole(int IsNorth,
		       Planets planet,
		       Inputs input,
		       Report &report);

/**************************************************************
 * get_bfield is the general call for getting the b-field vector
 * and the latitude and longitude of the position in the magnetic
 * coordinate system that is being used.
 * 
 * @param lon geographic longitude in radians
 *
 * @param lat geographic latitude in radians
 *
 * @param alt geographic altitude in meters
 *
 * @param planet Need to pass Planet class, so code can get info 
 *               about b-field
 *
 * @param input Need to pass Input class, so code can get info
 *              about how user has configured things
 *
 * @param report Need to pass Report class, so reporting can occur  
 **************************************************************/

bfield_info_type get_bfield(float lon,
                            float lat,
                            float alt,
                            Planets planet,
                            Inputs input,
                            Report &report);

/**************************************************************
 * get_dipole is the specific call for getting the b-field vector
 * and the latitude and longitude of the position in the magnetic
 * coordinate system for the planet-specific dipole. The planet-specific
 * characteristics are carried around in the Planets class.
 * 
 * @param lon geographic longitude in radians
 *
 * @param lat geographic latitude in radians
 *
 * @param alt geographic altitude in meters
 *
 * @param planet Need to pass Planet class, so code can get info 
 *               about b-field
 *
 * @param input Need to pass Input class, so code can get info
 *              about how user has configured things
 *
 * @param report Need to pass Report class, so reporting can occur  
 **************************************************************/

bfield_info_type get_dipole(float lon,
                            float lat,
                            float alt,
                            Planets planet,
                            Inputs input,
                            Report &report);

#endif // INCLUDE_BFIELD_H_
