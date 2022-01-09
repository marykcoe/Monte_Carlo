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

This file contains functions of the sim_box class. 
*/

#ifndef SIMBOX_H
#define SIMBOX_H

#include "constants.h"
#include "particle.h"
#include "cell.h"
#include "solute.h"

#include <iostream>
#include <string>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <vector>
#include <filesystem>

using namespace constants;

class sim_box {
	
	/*Class structure for the simulation box. */
	
	private: std::filesystem::path output_file_path, load_file_path;
			 std::string test_type;

			 bool load, slit, bsolute, mW, mc, tmmc, load_weights, sub_volume, vext, test, equilibrate;
			 bool compressibility, susceptibility, positions, output_individual_distributions;

			 long Lx=0, Ly=0, Lz=0, sweeps, equilibrium_sweeps, save_frequency, distribution_frequency;
			 long sub_volume_particles, num_particles, rand_check, sub_volume_dz, num_cells, seed;
			 long mca, mcr;

			 double Rs, epsilon, mu, temperature, minimum_density, maximum_density, initial_density;
			 double surface_interaction, volume, accessible_volume, *random_numbers=new double[nrands];

			 std::vector<double>  collection_totals, weights;
			 std::vector<std::vector<double>> collection, transition;
			 
			 double running_energy,  maximum_particles, minimum_particles;
			 double running_Vext, dmu, dT, avg_acc;

			 long total_sweeps, acc_entries;
			
			 int sv_particles, sv_max_particles, sv_min_particles, profile_bins;
			 double svolume;

	public: 	
			std::vector<particle> particle_list;
            std::vector<cell> cell_list;
			solute vsolute;
			std::vector<double> compressibility_profile, susceptibility_profile;
			std::vector<std::vector<double>> average_distribution;
			double sum_compressibility, sum_susceptibility, sum_distributions;
           	
	sim_box();


	// Functions which set initial parameters from Input or Output files.
	void set_output_file_path(std::string file_path);
	void set_load_file_path(std::string file_path);
	void set_load(bool iload);
	void set_positions(bool ipos);
	void set_Lx(int iLx);
	void set_Ly(int iLy);
	void set_Lz(int iLz);
	void set_slit(bool islit);
	void set_solute(bool isolute);
	void set_solute_radius(double iRs);
	void set_mW(bool ift);
	void set_temperature(double temp);
	void set_chemical_potential(double imu);
	void set_epsilon(double ieps);
	void set_minimum_density(double dens);
	void set_maximum_density(double dens);
	void set_initial_density(double dens);
	void set_sweeps(long isweeps);
	void set_equilibrate(bool ieq);
	void set_save_frequency(int freq);
	void set_distribution_frequency(int freq);
	void set_output_distributions(bool oup);
	void set_profile_bins(int npb);
	void set_compressibility(bool comp);
	void set_susceptibility(bool susc);
	void set_mc(bool imc);
	void set_load_weights(bool iload);
	void set_tmmc(bool itmmc);
	void set_sub_volume(bool isub_v);
	void set_sub_volume_dz(int idz);
	void set_vext(bool ivext);
	void set_surface_interaction_strength(double ieps);
	void set_test(bool itest);
	void set_test_type(std::string itest_type);
	void set_seed(long iseed);
	void set_mca(long imca);
	void set_mcr(long imcr);	
	void set_total_sweeps(long itsweeps);

	// Initialising simulation functions
    void setup(void);
    void place_particles(bool load, bool positions);

	// Load functions
	void load_particles();
	void load_previous_weights(bool load);
	void load_previous_collection(void);
	void load_previous_distributions(void);
	void read_distribution_file(std::filesystem::path fn, bool density, bool compressibility);

	// Printing functions
    void print_simulation_parameters(void);

	// Particle list functions
	void reorder_particles(int pn);
	
	// Transition Matrix sampling functions
	void update_collection(int N, int Na, double prob);
	void update_transition(void);
	void update_weights(long sweeps);
	void output_weights(std::string name);
	
	// Return functions
	int return_num_cells(void);
	double return_volume(void);
	bool return_load(void);
	bool return_test(void);
	std::string return_test_type(void);
	int return_Lx(void);
	int return_Ly(void);
	int return_Lz(void);
	long return_seed(void);
	long return_total_sweeps(void);
	long return_sweeps_eqb(void);
	long return_sweeps(void);
	int return_distribution_frequency(void);
	bool return_output_distributions(void);
	int return_profile_bins(void);
	int return_save_frequency(void);
	bool return_compressibility(void);
	bool return_susceptibility(void);
	double return_dmu(void);
	double return_dT(void);
	std::filesystem::path return_output_file_path(void);
	std::filesystem::path return_load_file_path(void);
			
	// Functions for random number generation.
	void check_rnd(void);
	double get_rnd(void);
			
	// Functions for move acceptance.
	void inc_mca(void);
	void inc_mcr(void);
	long return_mca(void);
	long return_mcr(void);
	void update_avg_acc(void);
	double return_avg_acc(void);
	long return_acceptance_entries(void);
			
	// Functions for  running energy.
	void inc_running_energy(double energy);
	void dec_running_energy(double energy);
	double return_running_energy(void);
	void set_running_energy(double energy);
			
	// Functions for the number of particles within the system.
	void inc_nparticles(void);
	void dec_nparticles(void);
	int return_nparticles(void);
	void set_nparticles(int particles);
	double return_maximum_particles(void);
	double return_minimum_particles(void);
	double return_minimum_density(void);
	double return_maximum_density(void);
	
	// Functions for system parameters.
	double return_mu(void);
	double return_epsilon(void);
	double return_temperature(void);
	bool return_mW(void);
	bool return_mc(void);
	bool return_tmmc(void);
	bool return_sub_volume(void);
	bool return_load_weights(void);
	bool return_equilibrate(void);

	// Functions for sub-volume sampling.
	int return_sub_volume_dz(void);
	int return_sv_particles(void);
	void inc_sv_particles(void);
	void dec_sv_particles(void);
	
	// Functions for weights information.
	double return_weight(int N);
	double return_prob(int Na, int Nb);
	
	// Functions for geometry.
	bool return_slit(void);
	bool return_solute(void);
	void add_solute(double Rs);
	bool check_valid_position(particle);

	// Functions for external potential energy.
	bool return_LJ_Vext(void);
	double return_epsilon_wall(void);
	void inc_running_Vext(double);
	void dec_running_Vext(double);
	double return_running_Vext(void);
	void set_running_Vext(double);

	// Debugging functions
	void output_cells(std::string);	
};
#endif
