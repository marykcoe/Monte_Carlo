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

#include "test.h"

void metropolis_test(sim_box* sim_env_ptr) {

    /*
    Purpose: To test the underlying Monte Carlo algorithm is working as expected. This test can
    also be used to test the incorporation of a random number generator has been successful.

    Method: The simulation is run for 0 fluid interaction strength, making the fluid effectively
    an ideal gas. In this case, the average density of the system should obey the equation
    chemical potential = ln(average density). 
    */

   // Set box parameters.
    sim_env_ptr->set_Lx(5); sim_env_ptr->set_Ly(5); sim_env_ptr->set_Lz(5);
    sim_env_ptr->set_slit(false); sim_env_ptr->set_solute(false); 
   
   // Set fluid parameters.
    sim_env_ptr->set_mW(false); sim_env_ptr->set_chemical_potential(-3.0); sim_env_ptr->set_temperature(1.0);
    sim_env_ptr->set_maximum_density(1.0); sim_env_ptr->set_minimum_density(0.0);
    sim_env_ptr->set_initial_density(0.0);

    // Set sampling parameters.
    sim_env_ptr->set_save_frequency(10); sim_env_ptr->set_sweeps(100000); 
    sim_env_ptr->set_distribution_frequency(5000);
    sim_env_ptr->set_equilibrate(false); sim_env_ptr->set_mc(false);
    sim_env_ptr->set_compressibility(false); sim_env_ptr->set_susceptibility(false);

    // Set other parameters.
    sim_env_ptr->set_load(false); sim_env_ptr->set_positions(false); 

    // Setup simulation box.
    sim_env_ptr->setup();

    // Set fluid interaction to 0.
    sim_env_ptr->set_epsilon(0.0); 

    // Run simulation.
    std::cout<<"Running metropolis test..."<<std::endl;
    std::cout<<"This will take a few minutes."<<std::endl;
    Run(sim_env_ptr);

    exit(0);
}

void cell_list_test(sim_box* sim_env_ptr) {

    /*
    Purpose: To test that the cell list has been set up correctly upon start-up.

    Method: A simulation box is setup for a user defined Lx, Ly and Lz. The cell list
    is generated and saved to file. Information saved includes position of cell 
    within cell list, location of cell in real-space and neighbours of cell. The program
    then exits. The accompanying Python program also creates a cell list. That output from
    the program is then compared to this Python generated version. If they match, then the
    test is passed.
    */
   std::cout<<"Running cell list test..."<<std::endl;
   sim_env_ptr->setup();
}

void critical_point_test(sim_box* sim_env_ptr) {

    /*
    Purpose: Critical points can be a very sensitive measure of whether the fluid potential is correct.
    For the fluids studied within this program, they are also well known. Hence, running a simulation
    at the critical point parameters should return a probability distribution which exhibits two peaks,
    one for liquid and one for vapour, and should have an average density within a specified range.

    Method: Runs a simulation at the critical point parameters for a specified fluid type. Density and
    energy statsitics are saved as the simulatio progresses. At the end of the simultion, the average
    density is printed to screen, along with that expected. The accompanying Python program then
    plots the density and energy distribution. If the average density is similar to that expected, 
    and the density and energy distributions exhibit two peaks, then it is likely that the fluid is
    behaving as expected and the test is therefore passed. There are two things to note. Firstly, the
    user will have to determine that there is a double peaked distribution. Secondly, this distribution
    is not expected to be perfect and the average density calculated is not expected to be exactly 
    equal to the reported result. Reasons for this are twofold - (i) the simulation run is very short
    and cannot be expected to provide super reliable statistics and (ii) finite size effects are not
    taken into account (we quote the infinite system critical density).
    */

   double energy = 0.0;

   // Set box parameters.
    sim_env_ptr->set_slit(false); sim_env_ptr->set_solute(false); 
   
   // Set fluid parameters.
   // mW critical point parameters taken from (available soon, see thesis above).
   // LJ critical point parameters taken from N. B. Wilding. Phys. Rev. E. 52(1) 602-11. (1995).
    if (sim_env_ptr->return_mW()) {
            sim_env_ptr->set_Lx(5); sim_env_ptr->set_Ly(5); sim_env_ptr->set_Lz(5);
            sim_env_ptr->set_chemical_potential(-5.011); sim_env_ptr->set_temperature(917.6);
            sim_env_ptr->set_maximum_density(0.75);
    }
    else {
        sim_env_ptr->set_Lx(5); sim_env_ptr->set_Ly(5); sim_env_ptr->set_Lz(5);
        sim_env_ptr->set_chemical_potential(-2.778); sim_env_ptr->set_temperature(1.1876);
        sim_env_ptr->set_maximum_density(0.55); 
    }
    sim_env_ptr->set_minimum_density(0.0); sim_env_ptr->set_initial_density(0.0);

    // Set sampling parameters.
    sim_env_ptr->set_save_frequency(10); sim_env_ptr->set_sweeps(10000000); 
    sim_env_ptr->set_distribution_frequency(100000);
    sim_env_ptr->set_equilibrate(true); sim_env_ptr->set_mc(false);
    sim_env_ptr->set_compressibility(false); sim_env_ptr->set_susceptibility(false);

    // Set other parameters.
    sim_env_ptr->set_load(false); sim_env_ptr->set_positions(false); 

    // Set up simulation box.
    sim_env_ptr->setup();

    // Run simulation
    if (sim_env_ptr->return_mW()) std::cout<<"Running critical point test for monatomic water..."<<std::endl;
    else std::cout<<"Running critical point test for lennard-jones fluid..."<<std::endl;
    std::cout<<"This will take approximately 60 minutes..."<<std::endl;

    // Equilibrate system first.
    Equilibrate(sim_env_ptr);

    // Run simulation with sampling.
    if (sim_env_ptr->return_mW()) energy = mWEnergy(sim_env_ptr);
	else energy = LJEnergy(sim_env_ptr);
	sim_env_ptr->set_running_energy(energy);
    Run(sim_env_ptr);

    exit(0);
}

void lj_external_potential_test(sim_box* sim_env_ptr) {

    /*
    Purpose: To confirm that the LJ external potential is calculating correctly.

    Method: Predetermined particle positions are read in for a given geometry (solute or slit)
    and the external potential felt by each particle is calculated. This is then compared against a 
    set of precalculated (from python) external potentials. 
    */

    double ext=0.0, z=0.0, mult, dx=0.0, dy=0.0, dz=0.0;
    std::filesystem::path path = sim_env_ptr->return_output_file_path()/"external_potential";
    std::ofstream output(path);
    output<<std::fixed<<std::setprecision(7);

    std::cout<<"Running lj external potential test..."<<std::endl;

    // Set essential simulation parameters but note we will not run a simulation.
    sim_env_ptr->set_sweeps(1); sim_env_ptr->set_save_frequency(1); sim_env_ptr->set_distribution_frequency(1);
    sim_env_ptr->setup();

    // Calculate length multiplier.
    if (sim_env_ptr->return_mW()) mult = constants::mW_a;
    else mult = constants::LJ_rc;

    // Calculate the external potential felt by each particle within the system and
    // write the result to file.
    for (int p=0; p<sim_env_ptr->return_nparticles(); p++) {
        ext = Vext_LJPotential(p, sim_env_ptr);

        if (sim_env_ptr->return_slit()) {
            z = mult*(sim_env_ptr->cell_list[sim_env_ptr->particle_list[p].cell_id].k + sim_env_ptr->particle_list[p].z);
        }
        else {
             dx = (sim_env_ptr->cell_list[sim_env_ptr->particle_list[p].cell_id].i + sim_env_ptr->particle_list[p].x)-sim_env_ptr->vsolute.x;
             dy = (sim_env_ptr->cell_list[sim_env_ptr->particle_list[p].cell_id].j + sim_env_ptr->particle_list[p].y)-sim_env_ptr->vsolute.y;
             dz = (sim_env_ptr->cell_list[sim_env_ptr->particle_list[p].cell_id].k + sim_env_ptr->particle_list[p].z)-sim_env_ptr->vsolute.z;
             z = mult*sqrt(dx*dx + dy*dy + dz*dz);
        }
        output<<z<<" "<<ext<<std::endl;
    }
    output.close();
    exit(0);
}
