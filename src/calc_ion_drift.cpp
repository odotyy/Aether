// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>

#include "../include/aether.h"

// --------------------------------------------------------------------------
// Calculate the electric field from the potential
// --------------------------------------------------------------------------

void Ions::calc_efield(Grid grid, Report &report) {

  // efield = - grad(potential)
  efield_vcgc = calc_gradient_vector(potential_scgc, grid);
  for (int64_t iComp = 0; iComp < 3; iComp++)
    efield_vcgc[iComp] = -efield_vcgc[iComp];
    
  // Remove component along b-field (should be zero, anyways!)
  fcube edotb = dot_product(efield_vcgc, grid.bfield_unit_vcgc);
  for (int64_t iComp = 0; iComp < 3; iComp++) 
    efield_vcgc[iComp] =
      efield_vcgc[iComp] - edotb % grid.bfield_unit_vcgc[iComp];
}

// --------------------------------------------------------------------------
// Calculate the E x B drift from the electric field and magnetic field
// --------------------------------------------------------------------------

void Ions::calc_exb_drift(Grid grid, Report &report) {
  fcube bmag2 =
    (grid.bfield_mag_scgc) % (grid.bfield_mag_scgc);
  exb_vcgc = cross_product(efield_vcgc, grid.bfield_vcgc);
  for (int64_t iComp = 0; iComp < 3; iComp++)
    exb_vcgc[iComp] = exb_vcgc[iComp] / bmag2;
}

// --------------------------------------------------------------------------
// 
// --------------------------------------------------------------------------

std::vector<fcube> Ions::calc_ion_electron_pressure_gradient(int64_t iIon,
							     Grid grid,
							     Report &report) {
  std::vector<fcube> pressure_gradient_vcgc;
  fcube total_pressure_scgc;

  // Total Pressure =
  //     Ion Pressure + Electron Pressure =
  //     (Ni * Ti + Ne * Te) * k
  
  total_pressure_scgc =
    (species[iIon].density_scgc %
     species[iIon].temperature_scgc +
     density_scgc %
     electron_temperature_scgc) *
    cKB;

  pressure_gradient_vcgc = calc_gradient_vector(total_pressure_scgc, grid);

  return pressure_gradient_vcgc;
}
  
// --------------------------------------------------------------------------
// Calculate the ion drift
// --------------------------------------------------------------------------

void Ions::calc_ion_drift(Neutrals neutrals, Grid grid, Report &report) {

  int64_t nX = grid.get_nX();
  int64_t nY = grid.get_nY();
  int64_t nZ = grid.get_nZ();
  
  calc_efield(grid, report);

  // This is for the electron drift motion:
  calc_exb_drift(grid, report);
  
  std::vector<fcube> gravity_vcgc = make_cube_vector(nX, nY, nZ, 3);
  std::vector<fcube> wind_forcing = make_cube_vector(nX, nY, nZ, 3);
  std::vector<fcube> total_forcing = make_cube_vector(nX, nY, nZ, 3);
  
  // This is assuming that the 3rd dim is radial.
  // Want actual gravity for 3rd dim
  // (negative is due to gravity being positive for some reason):
  gravity_vcgc[2] = -grid.gravity_scgc;

  calc_ion_neutral_coll_freq(neutrals, report);

  int64_t iIon, iNeutral;

  std::vector<fcube> grad_Pi_plus_Pe;
  fcube rho, rho_nuin;
  
  for (iIon = 0; iIon < nIons; iIon++) {

    // Need mass density for the current ion species:
    rho = species[iIon].mass * species[iIon].density_scgc;

    // Get gradient in pressure:
    grad_Pi_plus_Pe = calc_ion_electron_pressure_gradient(iIon,
							  grid,
							  report);

    // Continue assumption about gravity...
    gravity_vcgc[2] = -grid.gravity_scgc % rho;

    // Neutral Wind Forcing:
    for (int64_t iComp = 0; iComp < 3; iComp++)
      wind_forcing[iComp].zeros();
    for (iNeutral = 0; iNeutral < nSpecies; iNeutral++) {
      rho_nuin = rho % species[iIon].nu_ion_neutral_vcgc[iNeutral];
      for (int64_t iComp = 0; iComp < 3; iComp++)
	wind_forcing[iComp] = wind_forcing[iComp] +
	  rho_nuin % neutrals.velocity_vcgc[iComp];
    }

    // Total Forcing (sum everything - this is A_s):
    for (int64_t iComp = 0; iComp < 3; iComp++) {
      total_forcing[iComp] =
	- grad_Pi_plus_Pe[iComp]
	+ gravity_vcgc[iComp]
	+ wind_forcing[iComp];
    }
    
  }  

}


