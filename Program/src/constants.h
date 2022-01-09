/*
Monte Carlo program to simulate a truncated Lennard-Jones or
monatomic water liquid either in bulk, in contact with a small
solute, or confined to a slit. 

For information on how this program works, please consult the
accompanying tutorials or Section 4.2 of the following thesis:

M. K. Coe, Hydrophobicity Across Length Scales: The Role of
Surface Criticality, Ph.D. Thesis, University of Bristol (2021)
Available at: https://research-information.bris.ac.uk/ws/portalfiles/portal/304220732/Thesis_Mary_Coe.pdf

Copyright (c) 2022 Mary Coe
Made available under the MIT License.

This file contains definitions of constants used within the 
program. AIn most cases, reduced units are used.
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants {
	
	// Used in all simulations
	constexpr long nrands {10000};	 // Number of random numbers (used by some generators)
	constexpr double pi {3.14159};	
	constexpr double avagadro {6.022};
	constexpr long max_particles {34000};	// Maximum number of particles allowed within the simulation box
	constexpr double delta {1e-6};	// Used within histogram reweighting

	// Constants which relate to the truncated Lennard-Jones fluid
	constexpr double LJ_epsilon {3.3747};  
	constexpr double LJ_rc {2.5};
	constexpr double LJ_sigma {1.0};
	
	// Constants which relate to the monatomic water fluid
	constexpr double mW_sigma {2.3925};
	constexpr double mW_A {7.049556277};
	constexpr double mW_B {0.6022245584};
	constexpr double mW_gamma {1.2};
	constexpr double mW_a {1.8};
	constexpr double mW_aa {3.24};
	constexpr double mW_lambda {23.15}; 
	constexpr double mW_MW {18.015};
	
}

#endif
