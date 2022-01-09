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

This file contains functions to delete particles within the 
simulation box.
*/

#ifndef DELETE_H
#define DELETE_H

#include "particle.h"
#include "simbox.h"
#include "cell.h"
#include "energetics.h"

void delete_molecule(sim_box* sim_env_ptr, int pn, double energy, double Vext);
void decide_deletion_mc(int pn, sim_box* sim_env_ptr);
void decide_deletion(int pn, sim_box* sim_env_ptr);

#endif