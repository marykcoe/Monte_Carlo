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

#include "output.h"

void output_example_input_file(std::filesystem::path file_path) {
	
	/*
	Outputs an example of the Input file required for the simulation to the file path given
	in the command line when launching the program. 
	*/
	std::filesystem::path path = file_path/"Input_example";
	std::ofstream output(path);

	std::cout<<"Some essential arguments were missing from your input file."<<std::endl;
	std::cout<<"An example input file has been saved to the input directory."<<std::endl;
	
	output<<"# Input/Output Parameters:"<<std::endl;
	output<<"# These parameters set the input/output file directories. If LOAD is True"<<std::endl;
	output<<"# then the program will load the previous simulation from the LOAD_FILE_PATH directory."<<std::endl;
	output<<"# If LOAD is POSITIONS, then the program will load a previous configuration of particles"<<std::endl;
	output<<"# only. The simulation will then proceed as a new simulation."<<std::endl;
	output<<"# OUTPUT_FILE_PATH is the file_path to the directory to save the output of the simulation."<<std::endl;
	output<<std::endl;
	output<<"OUTPUT_FILE_PATH ./output/file/path/goes/here/"<<std::endl;
	output<<"LOAD TRUE/POSITIONS/FALSE"<<std::endl;
	output<<"LOAD_FILE_PATH ./input/file/path/goes/here/"<<std::endl;

	output<<std::endl;
	output<<"# Box parameters: "<<std::endl;
	output<<"# The dimensions of the box are governed by LX, LY and LZ, which each "<<std::endl;
	output<<"# represent the size of the box in X,Y and Z, in multiples of rc*sigma."<<std::endl;
	output<<"# For example, LX = 6 for a Lennard-Jones fluid is equivalent to the length"<<std::endl;
	output<<"# of the box in the X dimension being 6*2.5*sigma = 15sigma. "<<std::endl;
	output<<std::endl;
	output<<"LX INT"<<std::endl;
	output<<"LY INT"<<std::endl;
	output<<"LZ INT"<<std::endl;

	output<<std::endl;
	output<<"# These parameters indicate whether the box features two planar surfaces at "<<std::endl;
	output<<"# either end of the box in the z dimension, and is therefore a slit, or whether "<<std::endl;
	output<<"# there is a solute in the box. "<<std::endl;
	output<<"# SLIT takes the value TRUE/FALSE, and SOLUTE is either FALSE, or if a float, "<<std::endl;
	output<<"# represents the radius of the solute in multiples of sigma (e.g SOLUTE 3.0 = solute "<<std::endl;
	output<<"# of radius 3sigma centered on the middle of the box.) "<<std::endl;
	output<<std::endl;
	output<<"SLIT TRUE/FALSE"<<std::endl;
	output<<"SOLUTE FLOAT/FALSE"<<std::endl;

	output<<std::endl;
	output<<"# Fluid parameters: "<<std::endl;
	output<<"# This program supports a truncated Lennard-Jones potential, with truncation at"<<std::endl;
	output<<"# 2.5sigma, and monatomic water, which is a short-ranged water potential. Further"<<std::endl;
	output<<"# details on this potential are given in: "<<std::endl;
	output<<"# V. Molinero, and E.B. Moore, J. Phys. Chem. B 113. 4008-16 (2009)"<<std::endl;
	output<<"# TEMPERATURE represents the temperature of the system. "<<std::endl;
	output<<"# The chemical potential is in units of mu/k_bT. Example parameters for each fluid are:"<<std::endl;
	output<<std::endl;
	output<<"# LJ Points given in: N. B. Wilding. Phys. Rev. E. 52 (1) 602-11. (1995)."<<std::endl;
	output<<"# and: R. Evans, M.C. Stewart, and N.B. Wilding, J. Chem. Phys. 147. 044701 (2017)."<<std::endl;
	output<<"# LJ Critical Point: TEMPERATURE = 1.1876, CHEMICAL_POTENTIAL = -2.778"<<std::endl;
	output<<"# LJ 0.775Tc Coexistence: TEMPERATURE = 0.91954, CHEMICAL POTENTIAL = -3.86595"<<std::endl;
	output<<std::endl;
	output<<"# mW Points found in: Available Soon."<<std::endl;
	output<<"# mW Critical Point: TEMPERATURE = 917.6K, CHEMICAL_POTENTIAL = −5.110"<<std::endl;
	output<<"# mW Coexistence Points: "<<std::endl;
	output<<"# TEMPERATURE = 426K CHEMICAL_POTENTIAL: -10.66336"<<std::endl;
	output<<"# TEMPERATURE = 360K CHEMICAL_POTENTIAL: -12.75623"<<std::endl;
	output<<"# TEMPERATURE = 300K CHEMICAL_POTENTIAL: -15.53438"<<std::endl;
	output<<std::endl;
	output<<"FLUID_TYPE LJ/mW"<<std::endl;
	output<<"CHEMICAL_POTENTIAL FLOAT"<<std::endl;
	output<<"TEMPERATURE FLOAT"<<std::endl;
	output<<std::endl;
	output<<"# Densities are in terms of sigma^3 for LJ and gcm^-3 for mW."<<std::endl;
	output<<"INITIAL_DENSITY FLOAT"<<std::endl;
	output<<"MIN_DENSITY FLOAT"<<std::endl;
	output<<"MAX_DENSITY FLOAT"<<std::endl;
	output<<std::endl;

	output<<"# Sampling parameters: "<<std::endl;
	output<<"# SAVE_FRENQUENCY is the number of sweeps to wait before sampling the density and energy."<<std::endl;
	output<<"# DISTRIBUTION_FREQUENCY is the number of sweeps to wait before outputting the distribution"<<std::endl;
	output<<"# of the system. For a bulk fluid, this is the radial distribution function, for a slit"<<std::endl;
	output<<"# geometry, this is a distribution along the z dimension of the box, and for a solute"<<std::endl;
	output<<"# this is radial distribution as measured from the center of the solute."<<std::endl;
	output<<"# OUTPUT_INDIVIDUAL_DISTRIBUTIONS refers to whether to output each calculated distribution"<<std::endl;
	output<<"# to file (if TRUE) or whether to output only averge distributions (FALSE)."<<std::endl;
	output<<"# PROFILE_BINS refers to the number of histogram bins to use when measuring the particle "<<std::endl;
	output<<"# distributions within the system."<<std::endl;
	output<<"# COMPRESSIBILITY refers to whether to measure the local compressibility profile. It takes "<<std::endl;
	output<<"# the value TRUE/FALSE."<<std::endl;
	output<<"# SUSCEPTIBILITY refers to whether to measure the local thermal susceptibility profile."<<std::endl;
	output<<"# It takes the value TRUE/FALSE."<<std::endl;
	output<<std::endl;
	output<<"SWEEPS INT"<<std::endl;
	output<<"EQUILIBRATE TRUE/FALSE"<<std::endl;
	output<<"SAVE_FREQUENCY INT"<<std::endl;
	output<<"DISTRIBUTION_FREQUENCY INT"<<std::endl;
	output<<"OUTPUT_INDIVIDUAL_DISTRIBUTIONS TRUE/FALSE"<<std::endl;
	output<<"PROFILE_BINS INT"<<std::endl;
	output<<"COMPRESSIBILITY TRUE/FALSE"<<std::endl;
	output<<"SUSCEPTIBILITY TRUE/FALSE"<<std::endl;
			
	output<<std::endl<<"# MC turns on multicanonical sampling. TMMC uses the transition matrix method to"<<std::endl;
	output<<"# calculate the weights during the simulation. LOAD_WEIGHTS loads a previous configuration of"<<std::endl;
	output<<"# weights. This file must have the name 'weights_intial' and must be located in the same "<<std::endl;
	output<<"# directory as the input file. SUB_VOLUME uses sub-volume sampling in combination with TMMC to"<<std::endl;
	output<<"# predict weights function. It takes either FALSE or the number of cells in the Z dimension to"<<std::endl;
	output<<"# use for sampling. If using MC, either TMMC or LOAD_WEIGHTs must be chosen (or both)."<<std::endl;
	output<<"# For more information refer to accompanying documentation."<<std::endl<<std::endl;
	output<<"MC TRUE/FALSE"<<std::endl;
	output<<"LOAD_WEIGHTS TRUE/FALSE"<<std::endl;
	output<<"TMMC TRUE/FALSE"<<std::endl;
	output<<"SUB_VOLUME FLOAT/FALSE"<<std::endl;
	output<<std::endl;

	output<<"# External Potential Parameters: "<<std::endl;
	output<<"# These govern the potential exerted by the planar walls of the slit or solute on the fluid."<<std::endl;
	output<<"# The external potential supported in a long-ranged Lennard-Jones potential, which"<<std::endl;
	output<<"# is shifted such that the minimum occurs at the surface of the wall/solute. The"<<std::endl;
	output<<"# actual solute/wall is hard. The SURFACE_INTERACTION_STRENGTH parameter is a"<<std::endl;
	output<<"# multiple of the FLUID_INTERACTION_POTENTIAL."<<std::endl<<std::endl;

	output<<"LJ_VEXT TRUE/FALSE"<<std::endl;
	output<<"SURFACE_INTERACTION_STRENGTH FLOAT"<<std::endl;

	output<<std::endl;
	output<<"# Testing Parameters: "<<std::endl;
	output<<"# Several tests can be run to check the program is working. These are: "<<std::endl;
	output<<"# METROPOLIS: Runs a short simulation with non-interacting particles. The"<<std::endl;
	output<<"# chemical potental should then match the natural log of the average density. This tests"<<std::endl;
	output<<"# the underlying acceptance and rejection rules are obeyed. "<<std::endl;
	output<<"# CELL_LIST: Sets-up the cell list and exits the program. This is then compared against"<<std::endl;
	output<<"# a separate cell list within the accompanying Python tools to ensure it is working as"<<std::endl;
	output<<"# expected."<<std::endl;
	output<<"# CRITICAL_POINT: Runs a short simulation at the critical point of either fluid. Compares"<<std::endl;
	output<<"# the average density of the system to that expected (obtained from publications). The"<<std::endl;
	output<<"# distribution of the density and energy are also plotted within the accompanying tools"<<std::endl;
	output<<"# to ensure they are double peaked. This test ensures that the fluids are behaving as "<<std::endl;
	output<<"# as expected."<<std::endl;
	output<<"# LJ_EXTERNAL: Calculates the external potential felt by particles read-in from file. "<<std::endl;
	output<<"# This is then compared against separate calculations within the accompanying Python tools."<<std::endl;
	output<<"# This test ensures that the external potential is being calculated as expected."<<std::endl;
	output<<"# INT: If an integer is provided then this is used as the seed for the simulation.  This"<<std::endl;
	output<<"# helps with debugging why a simulation failed."<<std::endl;
	output<<"TEST FALSE/METROPOLIS/CELL_LIST/CRITICAL_POINT/LJ_EXTERNAL/INT"<<std::endl;
}

void output_particle_positions(sim_box* sim_env_ptr, std::string name) {

	/* 
	Writes to file the positions of all particles within the system. These
	are written in terms of the cell in which they reside and their 
	co-ordinates within the cell. This file is designed to be used when
	restarting a simulation. 
	*/

	std::ofstream positions; std::filesystem::path out_name;

	// Final positions are output to main output folder. All other position files
	// are saved to the 'particle_positions' folder. 
	if (name == "final") out_name = sim_env_ptr->return_output_file_path()/("positions_" + name);
	else out_name = sim_env_ptr->return_output_file_path()/"/particle_positions"/("positions_" + name);

	positions.open(out_name, std::ios::out);
	positions<<std::fixed<<std::setprecision(9);
	
	// Output cell and position within cell of each particle.
	for (int m=0;m<sim_env_ptr->return_nparticles();m++) {
        positions<<m<<" "<<sim_env_ptr->particle_list[m].x<<" "<<sim_env_ptr->particle_list[m].y<<" "
                        <<sim_env_ptr->particle_list[m].z<<" "<<sim_env_ptr->particle_list[m].cell_id<<std::endl;
    }
	positions.close();
}

void output_particle_positions_xyz(sim_box* sim_env_ptr, std::string name) {

	/* 
	Writes to file the real-space positions of all particles within the system.
	This file is designed to be used when visualising the system. 
	*/
	
	double x_pos, y_pos, z_pos, length_mult; cell mcell;
	std::ofstream positions; std::filesystem::path out_name;

	// Calculate real-space distance multipliers.
	if (sim_env_ptr->return_mW()) length_mult = mW_sigma*mW_a;
	else length_mult = LJ_rc*LJ_sigma;

	// Final positions are output to main output folder. All other position files
	// are saved to the 'particle_positions' folder. 
	if (name == "final") out_name = sim_env_ptr->return_output_file_path()/("positions_" + name + ".xyz");
	else out_name = sim_env_ptr->return_output_file_path()/"particle_positions"/("positions_" + name + ".xyz");
	positions.open(out_name, std::ios::out);

	// Preamble detailing total number of particles and fluid.
	positions<<sim_env_ptr->return_nparticles()<<std::endl;
	if (sim_env_ptr->return_mW()) positions<<"monatomic water"<<std::endl;
	else positions<<"lennard-jones fluid"<<std::endl;
	
	// Output coordinates of each particle.
	for (int m=0;m<sim_env_ptr->return_nparticles();m++) {
		mcell = sim_env_ptr->cell_list[sim_env_ptr->particle_list[m].cell_id];
		x_pos = (mcell.i+sim_env_ptr->particle_list[m].x)*length_mult;
		y_pos = (mcell.j+sim_env_ptr->particle_list[m].y)*length_mult;
		z_pos = (mcell.k+sim_env_ptr->particle_list[m].z)*length_mult;
		positions<<"O"<<" "<<std::setw(5)<<x_pos<<" "<<std::setw(5)<<y_pos<<" "<<std::setw(5)<<z_pos<<std::endl;
	}
	positions.close();
}

void output_data(bool weighted, bool sub_volume, bool LJ_Vext, long sweep[], double dens[], long nparticles[], long svparticles[], 
							double energy[], double acc_ratio[], std::filesystem::path file_path, double weight[], double Vext[]) 
	{
	
	/*
	Writes to file the density, number of particles, running energy, move acceptance ratio and,
	if applicable, the running external potential energy, bias weight and number of sub-volume
	particles.  
	*/

	std::ofstream data_out;
	data_out.open(file_path/"data",std::ios::app);
	data_out<<std::fixed<<std::setprecision(6);

	if (weighted) {
		if (sub_volume) {
			if (LJ_Vext) for(int i=0; i<11; i++) data_out<<sweep[i]<<" "<<dens[i]<<" "<<nparticles[i]<<" "<<svparticles[i]<<" "<<weight[i]<<" "<<energy[i]<<" "<<Vext[i]<<" "<<acc_ratio[i]<<std::endl; 
			else for(int i=0; i<11; i++) data_out<<sweep[i]<<" "<<dens[i]<<" "<<nparticles[i]<<" "<<svparticles[i]<<" "<<weight[i]<<" "<<energy[i]<<" "<<acc_ratio[i]<<std::endl;
		}
		else {
			if (LJ_Vext) for(int i=0; i<11; i++) data_out<<sweep[i]<<" "<<dens[i]<<" "<<nparticles[i]<<" "<<weight[i]<<" "<<energy[i]<<" "<<Vext[i]<<" "<<acc_ratio[i]<<std::endl; 
			else for(int i=0; i<11; i++) data_out<<sweep[i]<<" "<<dens[i]<<" "<<nparticles[i]<<" "<<weight[i]<<" "<<energy[i]<<" "<<acc_ratio[i]<<std::endl;
		}
	}
	else {
		if (LJ_Vext) for(int i=0; i<11; i++) data_out<<sweep[i]<<" "<<dens[i]<<" "<<nparticles[i]<<" "<<energy[i]<<" "<<Vext[i]<<" "<<acc_ratio[i]<<std::endl; 
		else for(int i=0; i<11; i++) data_out<<sweep[i]<<" "<<dens[i]<<" "<<nparticles[i]<<" "<<energy[i]<<" "<<acc_ratio[i]<<std::endl;
	}
	data_out.close();
}

void output_average_distribution(std::string name, sim_box* sim_env_ptr) {

	/*
	Writes to file the average spatial density distribution collected throughout
	the course of the simulation. 
	*/

	std::ofstream distribution;
	double volume = 0.0, length_mult=0.0, dens_mult = 1.0, dr=0.0;

	// Calculate the multiplier needed to convert reduced unit
	// distances to real-space unit distances. In the case of a
	// monatomic water fluid, also calculate the multiplier 
	// needed to convert the reduced density to real-unit
	// density.
	if (sim_env_ptr->return_mW()) { length_mult = mW_a*mW_sigma; dens_mult = (mW_MW/avagadro)*10;}
	else length_mult = LJ_rc*LJ_sigma;

	// Calculate the width of a density bin.
	dr = sim_env_ptr->average_distribution[1][0] - sim_env_ptr->average_distribution[0][0];

	// Calculate the volume of a bin in the slit geometry.
	// The shape of this bin is a slab. 
	if (sim_env_ptr->return_slit()) {
			volume = sim_env_ptr->return_Lx()*sim_env_ptr->return_Ly();
			volume *= length_mult*length_mult*dr;
	}

	// Loop over all bins and write to file the average density.
	// For a Lennard-Jones fluid, this density is in reduced units.
	// For a monatomic water fluid, this is in gcm^-3.
	distribution.open(sim_env_ptr->return_output_file_path()/"averaged_distributions"/("average_distribution_" + name),std::ios::out);
	distribution<<std::fixed<<std::setprecision(22);
	for (int i=0; i<sim_env_ptr->return_profile_bins(); i++) {
		
		// Calculate the volume of a bin if in the bulk or solute geometry.
		// The shape of this bin is a spherical shell.
		if (!sim_env_ptr->return_slit()) {
			volume = (sim_env_ptr->average_distribution[i][0] +(dr/2.))*(sim_env_ptr->average_distribution[i][0] + (dr/2.))*dr;
			volume *= 4.*pi;
		}
		distribution<<sim_env_ptr->average_distribution[i][0]<<" "<<dens_mult*sim_env_ptr->average_distribution[i][1]/(sim_env_ptr->sum_distributions*volume)<<std::endl;
	}
	distribution.close();

	// Write to file the 'raw' profile. This is simply the total
	// number of particles in each bin. Raw profiles can be used
	// to continue collecting samples when restarting a simulation.
	distribution.open(sim_env_ptr->return_output_file_path()/"averaged_distributions"/("average_distribution_" + name + "_raw"),std::ios::out);
	distribution<<std::fixed<<std::setprecision(22);
	for (int i=0; i<sim_env_ptr->return_profile_bins(); i++) distribution<<sim_env_ptr->average_distribution[i][0]<<" "<<sim_env_ptr->average_distribution[i][1]<<std::endl;
	distribution<<"Sum = "<<sim_env_ptr->sum_distributions<<std::endl;

	distribution.close();
}

void output_profile(std::filesystem::path name, sim_box* sim_env_ptr, bool compressibility) {

	/*
	Writes to file local compressibility and thermal susceptibility profiles 
	collected during the course of the simulation.
	*/

	double deriv  = 0.0, volume = 0.0, length_mult = 0.0, dr = 0.0,prev_volume = 0.0;
	double normalised_profile = 0.0, normalised_density = 0.0, dens_mult = 1.0;
	std::ofstream distribution;

	// Calculate the multiplier needed to convert reduced unit
	// distances to real-space unit distances. In the case of a
	// monatomic water fluid, also calculate the multiplier 
	// needed to convert the reduced density to real-unit
	// density.
	if (sim_env_ptr->return_mW()) {length_mult = mW_a*mW_sigma; dens_mult = (mW_MW/avagadro)*10;}
	else length_mult = LJ_rc*LJ_sigma;

	// Calculate the width of a density bin.
	dr = sim_env_ptr->average_distribution[1][0] - sim_env_ptr->average_distribution[0][0];

	// Calculate the volume of a bin in the slit geometry.
	// The shape of this bin is a slab. 
	if (sim_env_ptr->return_slit()) {
			volume = sim_env_ptr->return_Lx()*sim_env_ptr->return_Ly();
			volume *= length_mult*length_mult*dr;
	}

	// Write to file the averaged local compressibility or local thermal 
	// susceptibility profiles. 
	distribution.open(name,std::ios::out);
	distribution<<std::fixed<<std::setprecision(22);

	for (int i=0; i<sim_env_ptr->return_profile_bins(); i++) {

		// Calculate the average density profile at a higher chemical potential for use
		// when calculating the local compressibility.
		if (compressibility){
			normalised_profile = sim_env_ptr->compressibility_profile[i]/sim_env_ptr->sum_compressibility;
		}

		// Calculate the average density profile at a higher temperature for use
		// when calculating the local thermal susceptibility.
		else {
			normalised_profile = sim_env_ptr->susceptibility_profile[i]/sim_env_ptr->sum_susceptibility;
		}

		// Calculate the averaged density profile at the simulation temperature and chemical
		// potential.
		normalised_density = sim_env_ptr->average_distribution[i][1]/sim_env_ptr->sum_distributions;

		// Calculate the difference between the average local compressibility/thermal susceptibility
		// profile and the average density profile. This is the numerator of the numerical derivative
		// for calculating the local compressibility/thermal susceptibility.
		deriv =  normalised_profile - normalised_density;

		// Calculate the volume of a bin if in the bulk or solute geometry.
		// The shape of this bin is a spherical shell.
		if (!sim_env_ptr->return_slit()) {
			volume = (sim_env_ptr->average_distribution[i][0] +(dr/2.))*(sim_env_ptr->average_distribution[i][0] + (dr/2.))*dr;
            volume *= 4.*pi;
		}

		// Write to file the average local compressibility/thermal susceptibility. This requires
		// dividing by the difference in chemical potential/temperature and multiplying by various
		// factors.
		distribution<<sim_env_ptr->average_distribution[i][0]<<" "<<dens_mult*deriv/(volume*delta)<<std::endl;
	}

	distribution.close();

	// Write to file the 'raw' profile. This is the unaveraged density profile at the higher
	// chemical potential/temperature. The purpose of this file is to allow for the simulation
	// to continue adding to the average profile upon restarting.
	name += "_raw";
	distribution.open(name ,std::ios::out);
	distribution<<std::fixed<<std::setprecision(22);

	if (compressibility) {
		for (int i=0; i<sim_env_ptr->return_profile_bins(); i++) distribution<<sim_env_ptr->average_distribution[i][0]<<" "<<sim_env_ptr->compressibility_profile[i]<<std::endl;
	}
	else {
		for (int i=0; i<sim_env_ptr->return_profile_bins(); i++) distribution<<sim_env_ptr->average_distribution[i][0]<<" "<<sim_env_ptr->susceptibility_profile[i]<<std::endl;
	}

	if (compressibility) distribution<<"Sum = "<<sim_env_ptr->sum_compressibility<<std::endl;
	else distribution<<"Sum = "<<sim_env_ptr->sum_susceptibility<<std::endl;

	distribution.close();
}

void output_local_measures(std::string name, sim_box* sim_env_ptr) {

	/*
	Writes to file local measures (the local compressibility/thermal susceptibility) that
	the simulation is collecting. 
	*/

	std::filesystem::path file_name;

	// Output local compressibility
	if (sim_env_ptr->return_compressibility()) {
		file_name = sim_env_ptr->return_output_file_path()/"averaged_distributions"/("local_compressibility_" + name);
		output_profile(file_name, sim_env_ptr, true);
	}

	// Output local thermal susceptibility
	if (sim_env_ptr->return_susceptibility()) {
		file_name = sim_env_ptr->return_output_file_path()/"averaged_distributions"/("local_thermal_susceptibility_" + name);
		output_profile(file_name, sim_env_ptr, false);
	}
}

void output_sim_data(sim_box* sim_env_ptr, double end_energy){
		
	/* 
	Writes to file simulation parameters and information at the end of the simulation.
	This 'Output' file is required if restarting a simulation. 
	*/

	// Output file is written to same folder as Input file.
	std::ofstream output;
	output.open(sim_env_ptr->return_output_file_path()/"Output",std::ios::out);
	output<<std::fixed<<std::setprecision(7);

	// Information about the simulation box size and geometry.
	output<<"# Box parameters: "<<std::endl;
	output<<"LX "<<sim_env_ptr->return_Lx()<<std::endl;
	output<<"LY "<<sim_env_ptr->return_Ly()<<std::endl;
	output<<"LZ "<<sim_env_ptr->return_Lz()<<std::endl;
	if (sim_env_ptr->return_slit()) output<<"SLIT TRUE"<<std::endl;
	else output<<"SLIT FALSE"<<std::endl;
	if (sim_env_ptr->return_solute()) output<<"SOLUTE "<<sim_env_ptr->vsolute.Rs<<std::endl;
	else output<<"SOLUTE FALSE"<<std::endl;
	output<<std::endl;

	// Information about the fluid simulated.
	output<<"# Fluid parameters: "<<std::endl;
	if (sim_env_ptr->return_mW()) output<<"FLUID_TYPE mW"<<std::endl;
	else output<<"FLUID_TYPE LJ"<<std::endl;
	output<<"CHEMICAL_POTENTIAL "<<sim_env_ptr->return_mu()<<std::endl;
	output<<"TEMPERATURE "<<sim_env_ptr->return_temperature()<<std::endl;
	output<<"FLUID_INTERACTION_STRENGTH "<<sim_env_ptr->return_epsilon()<<std::endl;
	output<<"MIN_DENSITY "<<sim_env_ptr->return_minimum_density()<<std::endl;
	output<<"MAX_DENSITY "<<sim_env_ptr->return_maximum_density()<<std::endl;
	output<<std::endl;

	// Information about sampling that occured within the system. 
	output<<"# Sampling parameters: "<<std::endl;
	if (sim_env_ptr->return_equilibrate()) output<<"EQUILIBRATE TRUE"<<std::endl;
	else output<<"EQUILIBRATE FALSE"<<std::endl;
	output<<"SWEEPS "<<sim_env_ptr->return_sweeps()<<std::endl;
	output<<"TOTAL_SWEEPS "<<sim_env_ptr->return_total_sweeps() + sim_env_ptr->return_sweeps() + sim_env_ptr->return_sweeps_eqb()<<std::endl;
	output<<"SAVE_FREQUENCY "<<sim_env_ptr->return_save_frequency()<<std::endl;
	output<<"DISTRIBUTION_FREQUENCY "<<sim_env_ptr->return_distribution_frequency()<<std::endl;
	output<<"PROFILE_BINS "<<sim_env_ptr->return_profile_bins()<<std::endl;
	if (sim_env_ptr->return_compressibility()) output<<"COMPRESSIBILITY TRUE"<<std::endl;
	else output<<"COMPRESSIBILITY FALSE"<<std::endl;
	if (sim_env_ptr->return_susceptibility()) output<<"SUSCEPTIBILITY TRUE"<<std::endl;
	else output<<"SUSCEPTIBILITY FALSE"<<std::endl;
	if (sim_env_ptr->return_mc()) output<<"MC TRUE"<<std::endl;
	else output<<"MC FALSE"<<std::endl;
	if (sim_env_ptr->return_tmmc()) output<<"TMMC TRUE"<<std::endl;
	else output<<"TMMC FALSE"<<std::endl;
	if (sim_env_ptr->return_sub_volume()) output<<"SUB_VOLUME "<<sim_env_ptr->return_sub_volume_dz()<<std::endl;
	else output<<"SUB_VOLUME FALSE"<<std::endl;
	if (sim_env_ptr->return_load_weights()) output<<"LOAD_WEIGHTS TRUE"<<std::endl;
	else output<<"LOAD_WEIGHTS FALSE"<<std::endl;
	output<<"AVERAGE_MOVE_ACCEPTANCE "<<100.*sim_env_ptr->return_avg_acc()/sim_env_ptr->return_acceptance_entries()<<std::endl;
	output<<std::endl;

	// Information about any external potential supplied to the fluid.
	output<<"# External Potential Parameters: "<<std::endl;
	if (sim_env_ptr->return_LJ_Vext()) { 
		output<<"LJ_VEXT TRUE"<<std::endl;
		output<<"SURFACE_INTERACTION_STRENGTH "<<sim_env_ptr->return_epsilon_wall()<<std::endl;
	}
	else output<<"LJ_VEXT FALSE"<<std::endl;
	output<<std::endl;

	// Information about the seed for the random number generator.
	output<<"SEED "<<sim_env_ptr->return_seed()<<std::endl;
	output.close();		
}

