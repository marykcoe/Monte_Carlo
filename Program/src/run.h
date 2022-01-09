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

This file contains functions to run the simulation. 
*/

#ifndef RUN_H
#define RUN_H

#include "insert.h"
#include "delete.h"
#include "constants.h"
#include "particle.h"
#include "output.h"
#include <cstdlib>
#include <math.h>
#include <filesystem>

using namespace constants;

particle generate_new_particle(sim_box* sim_env_ptr);
void Equilibrate(sim_box* sim_env_ptr);
void Run(sim_box* sim_env_ptr);

#endif