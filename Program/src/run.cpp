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

This file contains functions to run the simulation. 
*/

#include "run.h"

particle generate_new_particle(sim_box* sim_env_ptr) {

	/* Generates a new particle at random for attempting insertion. */

    int pcell; double x,y,z;
    
	// Pick a cell at random.
    pcell = floor(sim_env_ptr->get_rnd()*(sim_env_ptr->return_num_cells()));

	// If the random number generator picked 1.0, then the cell picked will not
	// exist (we'll be off the end of the cell list), so instead pick the final cell.
	if (pcell>(sim_env_ptr->return_num_cells()-1)) pcell = sim_env_ptr->return_num_cells()-1;

	// Pick a position within the cell at random. 
	x = sim_env_ptr->get_rnd(); y = sim_env_ptr->get_rnd(); z = sim_env_ptr->get_rnd();

    return particle(x, y, z, pcell);
}

void Equilibrate(sim_box* sim_env_ptr) {
	
	/* 
	Equilibrate the system to ensure we reach thermodynamic equilibrium
	prior to sampling. This essentially means running the simulation for
	a given number of sweeps without sampling properties of the system. 
	*/
	
	long np, sweeps, isweep, reset_sweeps = sim_env_ptr->return_sweeps_eqb()/100;
	particle temp;
	bool valid = true;
	double rnd;

	// Run the simulation for a specified number of sweeps.
	for (sweeps=0;sweeps<sim_env_ptr->return_sweeps_eqb();sweeps++) {
		
		// A sweep is defined as a move attempt per number of cells.
		for (isweep = 0; isweep<sim_env_ptr->return_num_cells(); isweep++){

			rnd = sim_env_ptr->get_rnd(); // Generate a random number.

			// If the random number is greater than 0.5, attempt to insert a particle but
			// only if the number of particles is less than a specified maximum number of particles
			// allowed within the system.
			if (rnd > 0.5 && sim_env_ptr->return_nparticles()<sim_env_ptr->return_maximum_particles()) {

				temp = generate_new_particle(sim_env_ptr); // Generate new particle at random.

				// If the system contains a solute, check that the generated particle exists outside
				// the solute. If it does or if the system does not contain a solute, then the move
				// is valid and should be attempted. If not, immediately reject the attempted 
				// insertion move.
				if (sim_env_ptr->return_solute()) valid = sim_env_ptr->check_valid_position(temp);
				if (valid) {
					if (sim_env_ptr->return_mc()) decide_insertion_mc(temp,sim_env_ptr);	
               	 	else decide_insertion(temp,sim_env_ptr);
				}
				else sim_env_ptr->inc_mcr();
			}

			// If the random number is less than 0.5, attempt to delete a particle but
			// only if the number of particles is greater than a specified minimum number of particles
			// allowed within the system.
			else if (rnd < 0.5 && sim_env_ptr->return_nparticles()>sim_env_ptr->return_minimum_particles()) {

				// If the number of particles within the system is 1, then the only particle
				// we can attempt to delete is the solo particle.
				if (sim_env_ptr->return_nparticles() == 1) np = 0;

				// If the number of particles within the system is greater than 1, then choose
				// a particle to delete at random.
				else {
					np = floor(sim_env_ptr->get_rnd()*(sim_env_ptr->return_nparticles()));

					// If the random number generator returns 1.0, then the particle chosen is beyond
					// the number of particles in the system, hence instead choose the final particle
					// in the particle list.
					if (np == (sim_env_ptr->return_nparticles())) np = sim_env_ptr->return_nparticles()-1;
				}
				
				// Attempt to delete the particle.
				if (sim_env_ptr->return_mc()) decide_deletion_mc(np, sim_env_ptr);
                else decide_deletion(np, sim_env_ptr);
			}

			// Else reject the attempted move.
			else sim_env_ptr->inc_mcr();
		}

		// To prevent the number of accepted and rejected move tallies exceeding
		// their maximum size, we reset every few sweeps.
		if (sweeps%reset_sweeps == 0) {
			sim_env_ptr->set_mca(0); sim_env_ptr->set_mcr(0);
		}
	}	

	std::cout<<"Equilibration Complete"<<std::endl;

	// At the end of the equilibration period, calculate weights if using 
	// transition matrix method.
	if (sim_env_ptr->return_tmmc()) {
			sim_env_ptr->update_transition(); 
			sim_env_ptr->update_weights(0);
			std::cout<<"Updating Weights"<<std::endl;
	}
}


void Run(sim_box* sim_env_ptr) {
	
	/* 
	Runs the simulation, outputting data at specified intervals.
	*/
	
	long pn=0, snum=0, sweeps, isweep, avg_entries = 0, osweeps;
	double avg_dens=0, end_energy, avg_energy=0, acc_ratio, volume = sim_env_ptr->return_volume();
	double error = 0.0, rnd;
	particle temp;
	long sparticles[11], ssweep[11], ssvparticles[11], save_profiles = sim_env_ptr->return_sweeps()/20;
	double senergy[11],sVext[11],srat[11],sweight[11],sdens[11];
	bool valid=true;

	// Calculate the accessible volume within the system, if the system contains a solute.
	if (sim_env_ptr->return_solute()) volume -= (4./3.)*pi*pow(sim_env_ptr->vsolute.Rs,3);

	// Run the simulation for a specified number of sweeps.
	for (sweeps=0;sweeps<sim_env_ptr->return_sweeps();sweeps++) {

		// A sweep is defined as a move attempt per number of cells.
		for (isweep = 0; isweep<sim_env_ptr->return_num_cells();isweep++){

				rnd = sim_env_ptr->get_rnd(); // Generate a random number.

				// If the random number is greater than 0.5, attempt to insert a particle but
				// only if the number of particles is less than a specified maximum number of particles
				// allowed within the system.
				if ((rnd > 0.5) && (sim_env_ptr->return_nparticles() < sim_env_ptr->return_maximum_particles())) {
					
					temp = generate_new_particle(sim_env_ptr); // Generate new particle at random.

					// If the system contains a solute, check that the generated particle exists outside
					// the solute. If it does or if the system does not contain a solute, then the move
					// is valid and should be attempted. If not, immediately reject the attempted 
					// insertion move.
					if (sim_env_ptr->return_solute()) valid = sim_env_ptr->check_valid_position(temp);
					if (valid) {
						if (sim_env_ptr->return_mc()) decide_insertion_mc(temp, sim_env_ptr);
						else decide_insertion(temp, sim_env_ptr);
					}
					else sim_env_ptr->inc_mcr();
				}

				// If the random number is less than 0.5, attempt to delete a particle but
				// only if the number of particles is greater than a specified minimum number of particles
				// allowed within the system.
				else if ((rnd < 0.5) && sim_env_ptr->return_nparticles() > sim_env_ptr->return_minimum_particles()) {
					
					// If the number of particles within the system is 1, then the only particle
					// we can attempt to delete is the solo particle.
					if (sim_env_ptr->return_nparticles() == 1) pn = 0;
					else {
						// Choose a particle to attempt to delete at random.
						pn = floor(sim_env_ptr->get_rnd()*(sim_env_ptr->return_nparticles()));

						// If the random number generator returns 1.0, then the particle chosen is beyond
						// the number of particles in the system, hence instead choose the final particle
						// in the particle list.
						if (pn == (sim_env_ptr->return_nparticles())) pn = sim_env_ptr->return_nparticles()-1;
					}

					// Attempt to delete the particle
					if (sim_env_ptr->return_mc()) decide_deletion_mc(pn, sim_env_ptr);
					else decide_deletion(pn,sim_env_ptr);
				}

			// Else reject the attempted move.
			else sim_env_ptr->inc_mcr();
		}

		// Every given number of sweeps, the profile will save information about the 
		// current state of the simulation and write this to file if applicable.
		if ((sweeps%sim_env_ptr->return_save_frequency() == 0) && sweeps>0) {
			
			// Save the current density of the system. This is saved in real unit
			// for a monatomic water fluid (gcm^-3) and in reduced units for a 
			// truncated Lennard-Jones fluid.
			if (sim_env_ptr->return_mW()) sdens[snum]=((double)sim_env_ptr->return_nparticles()/(volume*mW_sigma*mW_sigma*mW_sigma))*(mW_MW/avagadro)*10;
			else sdens[snum] = (double)sim_env_ptr->return_nparticles()/volume;
			
			// If using multicanonical sampling, save the weight of the 
			// current number of particles within the system. If using the
			// transition matrix method to calculate the weights, then update
			// the necessary arrays.
        	if (sim_env_ptr->return_mc()){ 

				// The weight depends on whether the system is using the sub-volume 
				// sampling method.
				if (sim_env_ptr->return_sub_volume()) {
					sweight[snum] = sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles());
					ssvparticles[snum] = sim_env_ptr->return_sv_particles();
				}
				else sweight[snum] = sim_env_ptr->return_weight(sim_env_ptr->return_nparticles());

				if (sim_env_ptr->return_tmmc()) {
					sim_env_ptr->update_transition(); 
					sim_env_ptr->update_weights(sweeps);
				}
			}
			
			// Add the current system energy to the average energy total for the system.
			avg_energy += (sim_env_ptr->return_running_energy());

			// Save the average density.
			avg_dens+=sdens[snum];

			// Calculate and save the current move acceptance ratio within the system. 
			srat[snum] = ((double) sim_env_ptr->return_mca()/((double)(sim_env_ptr->return_mca())+(double)(sim_env_ptr->return_mcr())))*100.0;
			
			// Calculate the number of sweeps which have occurred, taking into account the previous
			// simulation if applicable.
			if (sim_env_ptr->return_load()) ssweep[snum] = sweeps + sim_env_ptr->return_total_sweeps();
			else ssweep[snum] = sweeps;

			// Save the current number of particles.
			sparticles[snum] = sim_env_ptr->return_nparticles();

			// Save the current energy and external potential energy of the system.
			senergy[snum] = sim_env_ptr->return_running_energy();
			if (sim_env_ptr->return_LJ_Vext()) {
				avg_energy += sim_env_ptr->return_running_Vext();
				sVext[snum] = sim_env_ptr->return_running_Vext();
			}

			// Every 11 samples, output the data collected. 
			if (snum%10==0 and snum>0) {
				output_data(sim_env_ptr->return_tmmc(),  sim_env_ptr->return_sub_volume(), sim_env_ptr->return_LJ_Vext(),ssweep, sdens, 
							sparticles, ssvparticles, senergy, srat, sim_env_ptr->return_output_file_path(), sweight, sVext);
				snum = 0;
				
				// Update the average move acceptance ratio.
				sim_env_ptr->update_avg_acc();
			}
			else snum++;

			avg_entries++; // This is used when calculated the average energy/density at the end of the simulation.
		}	
	
		// Every given number of sweeps, information about the current spatial density profile
		// and if applicable, the local compressibility and thermal susceptibility, is saved
		// and written to file. 
		if(sweeps%sim_env_ptr->return_distribution_frequency()==0) {

			// Calculate the number of sweeps, taking into account previous simulations.
			// This will be used in the output file name of the distribution. 
			if (sim_env_ptr->return_load()) osweeps = sweeps + sim_env_ptr->return_total_sweeps();
			else osweeps = sweeps;

			// Calculate the appropriate spatial density distribution for the geometry
			// of the system. 
			if (sim_env_ptr->return_slit()) planar_distribution(sim_env_ptr, osweeps);
			else if (sim_env_ptr->return_solute()) solute_rdf(sim_env_ptr, osweeps);
			else gr(sim_env_ptr, osweeps);

			// Write current averaged distributions to file.
			output_average_distribution(std::to_string(osweeps), sim_env_ptr);
			if (sim_env_ptr->return_compressibility() || sim_env_ptr->return_susceptibility()) output_local_measures(std::to_string(osweeps), sim_env_ptr);
			
			// Write current particle positions to file. 
			output_particle_positions(sim_env_ptr, std::to_string(osweeps));
			output_particle_positions_xyz(sim_env_ptr, std::to_string(osweeps));
		}
	}

	// At the end of the simulation, various statistics about the simulation are 
	// printed to the screen.
	// Calculate the final energy of the system. 
	if (sim_env_ptr->return_mW()) end_energy = mWEnergy(sim_env_ptr);
	else end_energy = LJEnergy(sim_env_ptr);
	if (sim_env_ptr->return_LJ_Vext()) end_energy += Vext_LJ(sim_env_ptr);

	std::cout << std::fixed;
    std::cout << std::setprecision(5);
	std::cout<<"Simulation Complete."<<std::endl<<std::endl;
	std::cout<<"Move Acceptance: "<<100.*sim_env_ptr->return_avg_acc()/sim_env_ptr->return_acceptance_entries()<<"%"<<std::endl;
	std::cout<<"Average Density: "<<avg_dens/avg_entries<<std::endl;
	std::cout<<"Average Energy: "<<avg_energy/avg_entries<<std::endl;

	// A good test that the simulation behaved as expected is to compare
	// the running energy to the final energy of the system.
	std::cout<<"Running and End energy agreed to a tolerance of "<<abs(sim_env_ptr->return_running_energy()+sim_env_ptr->return_running_Vext()-end_energy)<<std::endl;
	
	// Write the Output file.
	output_sim_data(sim_env_ptr, end_energy); 

	// Write to file the final weights of the system if using multicanonical sampling.
	if (sim_env_ptr->return_mc()) sim_env_ptr->output_weights("final");

	// If the program was performing certain tests, print the final results
	// of the test to the screen.
	if (sim_env_ptr->return_test()) {
		if (sim_env_ptr->return_test_type() == "metropolis") {
			error = abs(sim_env_ptr->return_mu()-log(avg_dens/avg_entries));
			std::cout<<std::endl<<"-----------------------------------------"<<std::endl;
			std::cout<<"Metropolis Test Results"<<std::endl;
			std::cout<<"Chemical Potential = "<<sim_env_ptr->return_mu()<<std::endl;
			std::cout<<"ln(<density>) = "<<log(avg_dens/avg_entries)<<std::endl;
			std::cout<<"Agreement = "<<error<<std::endl;
			error /= sim_env_ptr->return_mu();
			std::cout<<"Relative Error = "<<abs(error)<<std::endl;
			std::cout<<"-----------------------------------------"<<std::endl;
			if (abs(error)<1e-3) std::cout<<"TEST PASSED"<<std::endl;
			else std::cout<<"TEST FAILED"<<std::endl;
			std::cout<<"-----------------------------------------"<<std::endl;	
		}
		else if (sim_env_ptr->return_test_type() == "critical_point") {

			if (sim_env_ptr->return_mW()) error = abs(0.10396- avg_dens/avg_entries);
			else error = abs(0.3197 - avg_dens/avg_entries);
			std::cout<<std::endl<<"-----------------------------------------"<<std::endl;
			std::cout<<"Critical Point Test Results"<<std::endl;
			if (sim_env_ptr->return_mW()) std::cout<<"Critical Density = "<<0.322;
			else std::cout<<"Critical Density = "<<0.3197;
			std::cout<<" Average Simulation Density = "<<avg_dens/avg_entries<<std::endl;
			std::cout<<"Agreement = "<<error<<std::endl;
			if (sim_env_ptr->return_mW()) error /= 0.10396; 
			else error /= 0.3197;
			std::cout<<"Relative Error = "<<abs(error)<<std::endl;
			std::cout<<"PLEASE NOW PLOT TO CONFIRM DOUBLE PEAKED DISTRIBUTION."<<std::endl;
			std::cout<<"-----------------------------------------"<<std::endl;
		}
	}
}
