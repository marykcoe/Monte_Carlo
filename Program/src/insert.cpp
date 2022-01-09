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

#include "insert.h"

void insert_particle(sim_box* sim_env_ptr, particle p, double energy, double Vext) {
	
	/*
	Inserts particle into system, updating particle lists and running energy,
	*/
	
	// Add particle to relevant cell.
	sim_env_ptr->cell_list[p.cell_id].occupancy[sim_env_ptr->cell_list[p.cell_id].num_occupants] = sim_env_ptr->return_nparticles();
	sim_env_ptr->cell_list[p.cell_id].num_occupants+=1;

	sim_env_ptr->inc_nparticles(); // Increase number of particles
	sim_env_ptr->inc_running_energy(energy); // Increase running energy

	// Increase running external potential.
	if (sim_env_ptr->return_solute() || sim_env_ptr->return_slit()) sim_env_ptr->inc_running_Vext(Vext); 

	// If using sub-volume sampling, increase the number of sub-volume particles.
	if (sim_env_ptr->return_sub_volume() && sim_env_ptr->cell_list[p.cell_id].sv) { sim_env_ptr->inc_sv_particles(); } 

}

double insertion_prob(sim_box* sim_env_ptr, double dE) {
	
	/*
	Calculates the grand canonical probability of insertion.
	*/

	double prob;
	prob = exp((sim_env_ptr->return_mu()-dE));
	prob*= (sim_env_ptr->return_volume()/(sim_env_ptr->return_nparticles()+1));
    return prob;
}

void decide_insertion_mc(particle p, sim_box* sim_env_ptr){
	
	/* Decides whether to insert a given particle when using biased sampling,
	 and if so proceeds with insertion. */
		
	double dE, prob, prob_unbiased, Vext=0.0;

	// Add particle to particle list - this is necessary for energetics functions.
	sim_env_ptr->particle_list[sim_env_ptr->return_nparticles()] = p;
    
	// Calculate energy inserted particle would contribute to system.
	if (sim_env_ptr->return_mW()) dE = mWPotential(sim_env_ptr->return_nparticles(), sim_env_ptr);
    else dE = LJPotential(sim_env_ptr->return_nparticles(), sim_env_ptr);

	if (sim_env_ptr->return_LJ_Vext()) Vext = Vext_LJPotential(sim_env_ptr->return_nparticles(),sim_env_ptr);

	prob = insertion_prob(sim_env_ptr, dE+Vext); // Calculate unbiased probability of insertion. 

	if (prob>1.) prob_unbiased = 1.; // If probability is greater than one, set to 1.
	else prob_unbiased = prob;

    // Multiply insertion probability by appropriate multicanonical bias.
	if (sim_env_ptr->return_tmmc() && sim_env_ptr->cell_list[p.cell_id].sv) {
		
		// Transition matrix bias with sub-volume sampling.
		if (sim_env_ptr->return_sub_volume()) {
			sim_env_ptr->update_collection(sim_env_ptr->return_sv_particles(),2,prob_unbiased);
			sim_env_ptr->update_collection(sim_env_ptr->return_sv_particles(),1,1.-prob_unbiased);
			prob*=exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles())-sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()+1));
		}

		// Transition matrix bias.
		else {
			sim_env_ptr->update_collection(sim_env_ptr->return_nparticles(),2,prob_unbiased);
			sim_env_ptr->update_collection(sim_env_ptr->return_nparticles(),1,1.-prob_unbiased);
			prob*=exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles())-sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()+1));
		}
	}
	// Pre-set bias.
	else  if (!sim_env_ptr->return_tmmc()) prob*=exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles())-sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()+1));

	// Decide whether to accept or reject the insertion.
	if (prob>sim_env_ptr->get_rnd()) { // Accept insertion.
		insert_particle(sim_env_ptr, p, dE, Vext);
		sim_env_ptr->inc_mca();
	}
	else { // Reject insertion.
		sim_env_ptr->inc_mcr();
	}
}

void decide_insertion(particle p, sim_box* sim_env_ptr){
	
	/*
	Decides whether to insert a given particle when using unbiased sampling,
	and if so proceeds with insertion.
	*/
		
	double dE, prob, prob_unbiased, Vext=0.0;

	// Add particle to particle list - this is necessary for energetics functions.
	sim_env_ptr->particle_list[sim_env_ptr->return_nparticles()] = p;

    // Calculate energy inserted particle would contribute to system.
	if (sim_env_ptr->return_mW()) dE = mWPotential(sim_env_ptr->return_nparticles(), sim_env_ptr);
    else dE = LJPotential(sim_env_ptr->return_nparticles(), sim_env_ptr);

	if (sim_env_ptr->return_LJ_Vext()) Vext = Vext_LJPotential(sim_env_ptr->return_nparticles(),sim_env_ptr);

	prob = insertion_prob(sim_env_ptr, dE+Vext); // Calculate unbiased probability of insertion. 
	
	// Decide whether to accept or reject the insertion.
	if (prob>sim_env_ptr->get_rnd()) { // Accept insertion.
		insert_particle(sim_env_ptr, p, dE, Vext);
		sim_env_ptr->inc_mca();
	}
	else { // Reject insertion. 
		sim_env_ptr->inc_mcr();
	}
}