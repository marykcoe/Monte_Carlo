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

This file contains functions to calculate energies within the
system.
*/

#ifndef ENERGETICS_H
#define ENERGETICS_H

#include "particle.h"
#include "simbox.h"
#include "cell.h"
#include "constants.h"
#include "measures.h"
#include "output.h"

#include <deque>
#include <cstdlib>
#include <cmath>

using namespace constants;

double LJPotential (int pn, sim_box* sim_env_ptr);
double mWPotential (int pn, sim_box* sim_env_ptr);
double LJEnergy (sim_box* sim_env_ptr);
double mWEnergy(sim_box* sim_env_ptr);
double Vext_LJPotential(int pn, sim_box* sim_env_ptr);
double Vext_LJ( sim_box* sim_env_ptr);

#endif