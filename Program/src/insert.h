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

This file contains functions relating to the insertion of particles.
*/

#ifndef INSERT_H
#define INSERT_H

#include "simbox.h"
#include "particle.h"
#include "cell.h"
#include "energetics.h"
#include <cmath>

void insert_particle(sim_box* sim_env_ptr, particle p, double energy, double Vext);
double insertion_prob(sim_box* sim_env_ptr, double dE);
void decide_insertion_mc(particle p, sim_box* sim_env_ptr);
void decide_insertion(particle p, sim_box* sim_env_ptr);

#endif