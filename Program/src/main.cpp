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

This file runs the program.
*/

#include "constants.h"
#include "cell.h"
#include "particle.h"
#include "simbox.h"
#include "run.h"
#include "energetics.h"
#include "output.h"
#include "input.h"
#include "test.h"

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <filesystem>

using namespace constants;

void read_in_parameters(std::filesystem::path file_path, std::string file_name, sim_box* sim_env_ptr, bool load) {

	/*
	Opens and parses Input and Output files. Parsed lines are then
	used to set the parameters of the simulation.
	*/

	std::string line, arg1, arg2;
	std::ifstream input;
	bool space;
	
	// Try to open file. If the file cannot be found and is an Input
	// file, output an example Input file and exit. If the file 
	// cannot be found and is an Output file, warn user and exit.
	input.exceptions (std::ifstream::failbit | std::ifstream::badbit);
	try {
		input.open(file_path/file_name);
	}
	catch (std::ifstream::failure e) {
   		std::cerr<<"Could not locate "<<file_name<<std::endl;
		if (!load) {
			std::cerr<<"Writing example Input file to "<<file_path<<std::endl;
			output_example_input_file(file_path);
		}
		exit(0);
  	}
	
	// Read each line of the file, separating values by white space.
	// Lines which begin with a '#' are taken to be comments and are
	// thus ignored.
	try {
		while (std::getline(input, line)) {
			if (line[0] != '#' && line.length() != 0) { // Ignore comments and blank lines.
				for (auto i:line) {
					if (i == ' ') space = true;
					else if (!space)  arg1 = arg1 +i;
					else arg2= arg2 + i; 
				}
				set_input_paramater(sim_env_ptr, arg1, arg2, load); // Use parsed line to set simulation parameters.
				arg1.erase(); arg2.erase(); space = false;
			}
		}
	}
	catch (std::ifstream::failure e) { // End of file.
		input.close();
	}
}

int main(int argc, char* argv[]) {
	
	sim_box sim_env;
	sim_box* sim_env_ptr;
	std::filesystem::path file_path;
	std::ifstream input;
	bool space=false;
	double energy, vext;
	
	// Output message upon launching program.
	std::cout<<"*******************************************************************************"<<std::endl;
	std::cout<<"Grand Canonical Monte Carlo Program for simulating truncated Lennard-Jones"<<std::endl;
	std::cout<<"or monatomic water liquid either in bulk, confined to a slit, or in contact"<<std::endl;
	std::cout<<"with a solute."<<std::endl<<std::endl;
	std::cout<<"Copyright (c) Mary Coe (2022)"<<std::endl;
	std::cout<<"Made available under the MIT License."<<std::endl;
	std::cout<<"*******************************************************************************"<<std::endl<<std::endl;

	// Initialise the simulation box
	sim_env = sim_box();  sim_env_ptr = &sim_env;

	// Load the Input file.
	file_path = argv[1]; // File path to Input file is supplied as command line argument.
	
	// Read in simulation parameters from the Input file.
	read_in_parameters(file_path, "Input" ,sim_env_ptr, false);

	// In the simulation is a test simulation, launch correct test.
	if (sim_env_ptr->return_test()) {
		if (sim_env_ptr->return_test_type() == "metropolis") metropolis_test(sim_env_ptr);
		else if (sim_env_ptr->return_test_type() == "critical_point") critical_point_test(sim_env_ptr);
		else if (sim_env_ptr->return_test_type() == "cell_list") cell_list_test(sim_env_ptr);
		else if (sim_env_ptr->return_test_type() == "lj_external") lj_external_potential_test(sim_env_ptr);
	}

	// If the simulation is a continuation of a previous simulation, 
	// read the Output file of the previous simulation and set the
	//  necessary parameters.
	if (sim_env_ptr->return_load()) read_in_parameters(sim_env_ptr->return_load_file_path(), "Output", sim_env_ptr, true);

	// Set up simulation box
	sim_env_ptr->setup();

	// Output overview of simulation parameters.
	sim_env_ptr->print_simulation_parameters();

	// Equilibrate system if required.
	if (sim_env_ptr->return_equilibrate()) Equilibrate(sim_env_ptr);

	// Prepare to run simulation. First work out the initial energy of 
	// the system, taking into accoun the fluid-fluid interaction and
	// any supplied external potential.
	if (sim_env_ptr->return_mW()) energy = mWEnergy(sim_env_ptr);
	else energy = LJEnergy(sim_env_ptr);
	sim_env_ptr->set_running_energy(energy);

	if (sim_env_ptr->return_LJ_Vext()) {
		vext = Vext_LJ(sim_env_ptr);
		sim_env_ptr->set_running_Vext(vext);
	}
	else sim_env_ptr->set_running_Vext(0.0);

	// Run simulation. 
	Run(sim_env_ptr);

	// When simulation is finished, output final positions of 
	// particles and any other measures. 
	std::cout<<"Saving data..."<<std::endl;
	output_particle_positions(sim_env_ptr, "final");
	output_particle_positions_xyz(sim_env_ptr, "final");
	output_average_distribution("final", sim_env_ptr);
	if (sim_env_ptr->return_compressibility() || sim_env_ptr->return_susceptibility()) output_local_measures("final", sim_env_ptr);
	std::cout<<"Simulation complete."<<std::endl;	
}
