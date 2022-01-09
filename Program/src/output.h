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

This file contains functions which write data gathered within the
simulation or Input/Output information to file.
*/

#ifndef OUTPUT_H
#define OUTPUT_H

#include "cell.h"
#include "particle.h"
#include "simbox.h"
#include "constants.h"

#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <filesystem>

using namespace constants;

void output_example_input_file(std::filesystem::path file_path) ;
void output_particle_positions(sim_box* sim_env_ptr, std::string name);
void output_particle_positions_xyz(sim_box* sim_env_ptr, std::string name);
void output_data(bool tmmc, bool sub_volume, bool LJ_Vext,long sweep[], 
        double dens[], long nparticles[], long sv_particles[],double energy[], double acc_ratio[], 
        std::filesystem::path file_path, double weight[], double Vext[]);
void output_sim_data(sim_box* sim_env_ptr, double end_energy);
void output_average_distribution(std::string name, sim_box* sim_env_ptr);
void output_profile(std::string name, sim_box* sim_env_ptr, bool compressibility);
void output_local_measures(std::string name, sim_box* sim_env_ptr);

#endif
