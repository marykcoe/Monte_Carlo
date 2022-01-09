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

This file contains the functions to test the program. 
*/

#ifndef TEST_H
#define TEST_H

#include "simbox.h"
#include "run.h"
#include "energetics.h"
#include "constants.h"
#include <fstream>
#include <filesystem>

void metropolis_test(sim_box* sim_env_ptr);
void cell_list_test(sim_box* sim_env_ptr);
void critical_point_test(sim_box* sim_env_ptr);
void lj_external_potential_test(sim_box* sim_env_ptr);

#endif