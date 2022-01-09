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

This file contains functions to parse the input and output files
for the program.
*/

#include "input.h"

void set_input_paramater(sim_box* sim_env_ptr, std::string name, std::string value, bool load) {

    /* 
    Uses arguments from Input and Output files to set values within the simulation
    environment instance. These inputs control features of the simulation. 
    */

    bool ignore;

    // Set file paths.
    if (name == "OUTPUT_FILE_PATH") sim_env_ptr->set_output_file_path(value);
    else if (name == "LOAD_FILE_PATH") sim_env_ptr->set_load_file_path(value);

    // Set whether and how to load previous simulation.
    else if (name == "LOAD") {
        if (value == "TRUE") sim_env_ptr->set_load(true);
        else if (value=="POSITIONS") sim_env_ptr->set_positions(true);
        else if (value != "FALSE") std::cout<<"Invalid LOAD value. LOAD set to FALSE."<<std::endl;
    }

    // Set dimensions of simulation box.
    else if (name == "LX") sim_env_ptr->set_Lx(std::stoi(value));
    else if (name == "LY") sim_env_ptr->set_Ly(std::stoi(value));
    else if (name == "LZ") sim_env_ptr->set_Lz(std::stoi(value));

    // Set geometry of box.
    else if (name == "SLIT") {
        if (value == "TRUE") sim_env_ptr->set_slit(true);
        else if (value != "FALSE") std::cout<<"Invalid SLIT value. SLIT set to FALSE."<<std::endl;
    }
     else if (name == "SOLUTE") {
        if (value != "FALSE") {
            sim_env_ptr->set_solute(true);
            sim_env_ptr->set_solute_radius(std::stof(value));
        }
    }

    // Set fluid type - either truncated Lennard-Jones fluid (LJ) or 
    // monatomic water (mW).
    else if (name == "FLUID_TYPE") {
        if (value == "mW") sim_env_ptr->set_mW(true);
        else if (value == "LJ") sim_env_ptr->set_mW(false);
        else {
            std::cout<<"Did not recognise fluid type. Defaulting to LJ."<<std::endl;
            sim_env_ptr->set_mW(false);
        }
    }

    // Set temperature and chemical potential. Note these are in reduced units.
    else if (name == "TEMPERATURE" && !load) sim_env_ptr->set_temperature(std::stof(value));
    else if (name == "CHEMICAL_POTENTIAL" && !load) sim_env_ptr->set_chemical_potential(std::stof(value));

    // Set density parameters. Min and Max densities prevent program trying to explore
    // non-physical and therefore unimportant regions when using multicanonical
    // sampling.
    else if (name == "MIN_DENSITY") sim_env_ptr->set_minimum_density(std::stof(value));
    else if (name == "MAX_DENSITY") sim_env_ptr->set_maximum_density(std::stof(value));
    else if (name == "INITIAL_DENSITY") sim_env_ptr->set_initial_density(std::stof(value));

    // Set sampling parameters. 
    else if (name == "SWEEPS" and !load) sim_env_ptr->set_sweeps(std::stoi(value));
    else if (name == "SAVE_FREQUENCY") sim_env_ptr->set_save_frequency(std::stoi(value));
    else if (name == "DISTRIBUTION_FREQUENCY") sim_env_ptr->set_distribution_frequency(std::stoi(value));
    else if (name == "OUTPUT_INDIVIDUAL_DISTRIBUTIONS") {
        if (value == "TRUE") sim_env_ptr->set_output_distributions(true);
        else if (value == "FALSE") sim_env_ptr->set_output_distributions(false);
        else {
                std::cout<<"Invalid OUTPUT_INDIVIDUAL_DISTRIBUTIONS value. OUTPUT_INDIVIDUAL_DISTRIBUTIONS set to TRUE."<<std::endl;
                sim_env_ptr->set_output_distributions(true);
        }
    }
    else if (name == "PROFILE_BINS") sim_env_ptr->set_profile_bins(std::stoi(value));
    else if (name == "EQUILIBRATE") {
        if (value == "TRUE") sim_env_ptr->set_equilibrate(true);
        else if (value == "FALSE") sim_env_ptr->set_equilibrate(false);
        else {
            std::cout<<"Invalid EQUILIBRATE value. EQUILIBRATE set to TRUE."<<std::endl;
            sim_env_ptr->set_equilibrate(true);
        }
    }

    // Set whether to use multicanonical sampling, and if so how to find the weights.
    else if (name == "MC") {
        if (value == "TRUE") sim_env_ptr->set_mc(true);
        else if (value != "FALSE") std::cout<<"Invalid MC value. MC set to FALSE."<<std::endl;
    }
    else if (name == "LOAD_WEIGHTS") {
        if (value == "TRUE") sim_env_ptr->set_load_weights(true);
        else if (value != "FALSE") std::cout<<"Invalid LOAD_WEIGHTS value. LOAD_WEIGHTS set to FALSE."<<std::endl;
    }
    else if (name == "TMMC") {
        if (value == "TRUE") sim_env_ptr->set_tmmc(true);
        else if (value != "FALSE") std::cout<<"Invalid TMMC value. TMMC set to FALSE."<<std::endl;
    }
    else if (name == "SUB_VOLUME") {
        if (value == "TRUE") {
            sim_env_ptr->set_sub_volume(true);
            sim_env_ptr->set_sub_volume_dz(std::stoi(value));
        }
        else if (value != "FALSE") std::cerr<<"Invalid SUB_VOLUME value. SUB_VOLUME set to FALSE."<<std::endl;
    }

    // Set whether to use external potential and, if so, strength of interactions in reduced energy units.
    else if (name == "LJ_VEXT") {
        if (value == "TRUE") sim_env_ptr->set_vext(true);
        else if (value != "FALSE") std::cout<<"Invalid LJ_VEXT value. LJ_VEXT set to FALSE."<<std::endl;
    }
    else if (name == "SURFACE_INTERACTION_STRENGTH") sim_env_ptr->set_surface_interaction_strength(std::stof(value));

    // Set whether to measure the local compressibility and thermal susceptibility.
    else if (name == "COMPRESSIBILITY") { 
            if (value == "TRUE") sim_env_ptr->set_compressibility(true);
            else if (value != "FALSE") std::cout<<"Invalid COMPRESSIBILITY value. COMPRESSIBILITY set to FALSE."<<std::endl;
    }
    else if (name == "SUSCEPTIBILITY") { 
            if (value == "TRUE") sim_env_ptr->set_susceptibility(true);
            else if (value != "FALSE") std::cout<<"Invalid SUSCEPTIBILITY value. SUSCEPTIBILITY set to FALSE."<<std::endl;
    }

    // Set whether simulation is a test and if so which type.
    else if (name == "TEST") {
        if (value == "FALSE") sim_env_ptr->set_test(false);
        else if (value == "METROPOLIS") {
            sim_env_ptr->set_test(true);
            sim_env_ptr->set_test_type("metropolis");
        }
        else if (value == "CRITICAL_POINT") {
            sim_env_ptr->set_test(true);
            sim_env_ptr->set_test_type("critical_point");
        }
        else if (value == "CELL_LIST") {
            sim_env_ptr->set_test(true);
            sim_env_ptr->set_test_type("cell_list");
        }
        else if (value == "LJ_EXTERNAL") {
            sim_env_ptr->set_test(true);
            sim_env_ptr->set_test_type("lj_external");
        }
        else {
            sim_env_ptr->set_test(true);
            sim_env_ptr->set_test_type("seed");
            sim_env_ptr->set_seed(std::stoi(value));
        }
    }

    // These are parameters read in from an Output file when loading a simulation.
    else if (name == "MCA") ignore = true;
    else if (name == "MCR") ignore = true;
    else if (name == "AVERAGE_MOVE_ACCEPTANCE") ignore = true;
    else if (name == "TOTAL_SWEEPS") sim_env_ptr->set_total_sweeps(std::stol(value));
    else if (name == "SEED" and !load) sim_env_ptr->set_seed(std::stoi(value));
    else if (name == "FLUID_INTERACTION_STRENGTH") ignore = true;
    else {
        std::cout<<"Argument invalid. Ignoring argument: "<<name<<std::endl;
    }
}
