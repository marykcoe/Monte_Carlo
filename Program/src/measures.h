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

This file contains functions which perform measures on and
within the system.
*/

#ifndef MEASURES_H
#define MEASURES_H

#include "constants.h"
#include "cell.h"
#include "particle.h"
#include "simbox.h"
#include "output.h"

#include <fstream>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <math.h>
#include <filesystem>

using namespace constants;

double pair_dist(double xa, double xb, int cx);
double three_body_dist(double xa, double xb, int cx, int L);
double three_body_dist(double xa, double xb, int cx);
double Find_distance(particle pa, particle pb, int i, int j, int k);
double Find_distance(particle pa, particle pb, int i, int j, int k, int Lx, int Ly, int Lz, bool slit);
double Find_angle(double ar2, double br, double cr);
void ang_dist(sim_box sim_env, int iter);
void gr(sim_box* sim_env_ptr, int sweep);
void planar_distribution(sim_box* sim_env_ptr, int sweep);
void solute_rdf(sim_box* sim_env_ptr, int sweep);
void local_compressibility(double profile[], sim_box* sim_env_ptr);
void local_susceptibility(double profile[], sim_box* sim_env_ptr);

#endif