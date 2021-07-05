// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md


#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iostream>

#include "../include/aether.h"

// *will have to change units in code*


// -----------------------------------------------------------------------------
// Function to Calculate ionization rate for 1 Ebin for 1 alt profile
// -----------------------------------------------------------------------------
fvec calculate_fang(float eflux,  // in ergs/cm2/s
                    float avee,   // in keV
                    float Ebin,   // eV
                    fvec rhoH, 
                    std::vector<float> Ci, 
                    float dE,     // eV
                    fvec H, 
                    Report &report)
                    {

    // report 
    std::string function = "calc_fang";
    static int iFunction = -1;
    report.enter(function, iFunction);

    //std::cout<<"eflux: "<<eflux<<", avee: "<<avee<<std::endl;
    float de = 0.035;  // keV

    float char_e = avee / 2;  // keV
    long double Q0 = eflux * 6.242e11;   // eV/cm2/s 
    float E0 = char_e * 1000;  // eV

    long double a = Q0 / 2 / (E0 * E0 * E0);  // cm2/s/eV2
    long double flux_max = a * Ebin * exp(-Ebin / E0);  //  /cm2/s/eV
    long double QE = flux_max * Ebin * dE; //  eV/cm2/s

    //std::cout<<" Ebin: "<<Ebin<<std::endl;
    //std::cout<<" final flux for calculations: "<<QE<<std::endl;
    //std::cout<<" dE: "<<dE<<std::endl;
    //std::cout<<"QE: "<<QE<<std::endl;

    // !!!!! ebin and QE in keV !!!!!! /1000 or not like in pyhton code?
    long double E_keV = Ebin / 1000.0;
    long double Q_keV = QE / 1000.0;

    fvec yE = (2.0/E_keV) * pow( (rhoH / (6e-6)), 0.7);
    //std::cout<<" yE: "<<std::endl;
    //yE.raw_print();

    fvec fyE =
      (Ci[0] * pow(yE, Ci[1])) % exp((-Ci[2] * pow(yE, Ci[3]))) + 
      (Ci[4] * pow(yE, Ci[5])) % exp((-Ci[6] * pow(yE, Ci[7])));

    fvec fac = Q_keV / de / H; 

    fvec qtot = fyE % fac; 
    
    report.exit(function);
    return qtot; 
}



// -----------------------------------------------------------------------------
// Calculate aurora
// -----------------------------------------------------------------------------
void calc_aurora(Grid grid, 
		 Neutrals &neutrals, 
		 Ions &ions, 
		 Report &report) {

  std::string function = "calc_aurora";
  static int iFunction = -1;
  report.enter(function, iFunction);

  // some variables
  int64_t iSpecies;
  int64_t iAlt, iLon, iLat;
  int64_t nLons = grid.get_nLons();
  int64_t nLats = grid.get_nLats();
  int64_t nAlts = grid.get_nAlts();

  // DENSITY INTEGRAL CALULATION ( done in calc_neutral_derived.cpp line
  // 170 rho_alt_int_scgc species[iSpecies].rho_alt_int_scgc =
  // integral3d * species[iSpecies].mass;

  // SET UP PIJ VALUES 
  static mat Pij = { {1.25, 1.45903, -2.42e-1, 5.95e-2},
		     {2.24, -4.23e-7, 1.36e-2, 2.53e-3},
		     {1.42, 1.45e-1, 1.70e-2, 6.40e-4},
		     {0.248775, -1.51e-1, 6.31e-9, 1.24e-3},
		     {-0.465119, -1.05e-1, -8.96e-2, 1.22e-2},
		     {3.86e-1, 1.75e-3, -7.43e-4, 4.61e-4},
		     {-6.45e-1, 8.50e-4, -4.29e-2, -2.99e-3},
		     {9.49e-1, 1.97e-1, -2.51e-1, -2.07e-3} }; 

  static std::vector<std::vector<float>> CiArray;
  static bool IsFirstTime = 1;

  // ENERGY BINS AND DE (E in eV)
  static float min = 100;
  static float max = 1000000;
  static float Emin = log(min);
  static float Emax = log(max);
  static int nBins = 101;
  static std::vector<float> auroral_energies;
  static std::vector<float> auroral_energy_widths;
  std::vector<float> Ci;

  if (IsFirstTime) {
  
    float lnE;

    for (int64_t iBin = 0; iBin < nBins; iBin++) {
      float energy = exp(Emin + iBin*(Emax-Emin)/(nBins-1));
      auroral_energies.push_back(energy);  
    }
    auroral_energy_widths = calc_bin_widths(auroral_energies);
    
    for (int64_t iBin = 0; iBin < nBins; iBin++) {

      lnE = log(auroral_energies[iBin]/1000.0);  // eV -> keV

      // loop through Pij values to get vector of Ci values
      for(int i = 0; i < 8; i++){
	float tot = 0;
	for(int j = 0; j < 4; j++){
	  tot = tot +  Pij.at(i,j) * pow(lnE , j); 
	}
	Ci.push_back(exp(tot)); 
      }
      
      CiArray.push_back(Ci);
    }

    IsFirstTime = 0;
  }
    
  fvec rhoH1d;
  fcube scale_height;
  fvec ionization1d;
  fvec H;

  fvec rho_tube; 

  fvec weighted_sum;
  fvec O_density_tube;
  fvec N2_density_tube;
  fvec O2_density_tube;
  fvec O;
  fvec N2;
  fvec O2;

  rhoH1d.set_size(nAlts);
  ionization1d.set_size(nAlts);

  scale_height = cKB * neutrals.temperature_scgc /
    (neutrals.mean_major_mass_scgc % abs(grid.gravity_scgc));

  float eflux;
  float avee;

  // loop through each altitude and calculate ionization
  for (iLon = 0; iLon < nLons ; iLon++) {
    for (iLat = 0; iLat < nLats ; iLat++) {

      eflux = ions.eflux(iLon, iLat);
      avee = ions.avee(iLon, iLat);

      if (eflux > 0.01) {
      
	rhoH1d.zeros(); 
	for (iSpecies=0; iSpecies < nSpecies; iSpecies++) {
	  rho_tube = neutrals.species[iSpecies].rho_alt_int_scgc.tube(iLon, iLat);
	  rhoH1d = rhoH1d + rho_tube;
	}

	// auroral_energy_widths.at(iBin), ions.eflux(iLon,iLat), ions.avee(iLon,iLat),
	ionization1d.zeros();
	H = scale_height.tube(iLon,iLat);
	fvec temp;

	for(int iBin = 0; iBin < nBins; iBin++){
	  Ci = CiArray[iBin];
	  temp = calculate_fang(eflux, avee,
				auroral_energies.at(iBin),
				rhoH1d,
				Ci,
				auroral_energy_widths.at(iBin),
				H,
				report);
	  
	  ionization1d = ionization1d + temp;  
	}

	// FINAL IONIZATION 
        
	O_density_tube = neutrals.species[iO_].density_scgc.tube(iLon, iLat);
	N2_density_tube = neutrals.species[iN2_].density_scgc.tube(iLon, iLat);
	O2_density_tube = neutrals.species[iO2_].density_scgc.tube(iLon, iLat);
            
	weighted_sum =
	  1.00 * O2_density_tube + 
	  0.92 * N2_density_tube +
	  0.56 * O_density_tube;
            
	// O ionization
	O = neutrals.species[iO_].ionization_scgc.tube(iLon, iLat);
	neutrals.species[iO_].ionization_scgc.tube(iLon, iLat) =
	  O + (0.56 * ionization1d % O_density_tube / weighted_sum);

	O = ions.species[iOP_].ionization_scgc.tube(iLon, iLat);
	ions.species[iOP_].ionization_scgc.tube(iLon, iLat) =
	  O + (0.56 * ionization1d % O_density_tube / weighted_sum);
             
	// N2 ionization
	N2 = neutrals.species[iN2_].ionization_scgc.tube(iLon, iLat);
	neutrals.species[iN2_].ionization_scgc.tube(iLon, iLat) =
	  N2 + (0.92 * ionization1d % N2_density_tube / weighted_sum);

	N2 = ions.species[iN2P_].ionization_scgc.tube(iLon, iLat);
	ions.species[iN2P_].ionization_scgc.tube(iLon, iLat) =
	  N2 + (0.92 * ionization1d % N2_density_tube / weighted_sum);

	// O2 ionization
	O2 = neutrals.species[iO2_].ionization_scgc.tube(iLon, iLat);
	neutrals.species[iO2_].ionization_scgc.tube(iLon, iLat) =
	  O2 + (1.00 * ionization1d % N2_density_tube / weighted_sum);

	O2 = ions.species[iO2P_].ionization_scgc.tube(iLon, iLat);
	ions.species[iO2P_].ionization_scgc.tube(iLon, iLat) =
	  O2 + (1.00 * ionization1d % N2_density_tube / weighted_sum);	

      }
            
    }
  }
 
  report.exit(function);
}














