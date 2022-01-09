/*
Monte Carlo program to simulate a truncated Lennard-Jones or
monatomic water liquid either in bulk, in contact with a small
solute, or confined to a slit. 

For information on how this program works, please consult the
accompanying tutorials or Section 4.2 of the following thesis:

M. K. Coe, Hydrophobicity Across Length Scales: The Role of
Surface Criticality, Ph.D. Thesis, University of Bristol (2021)
Available at: (https://research-information.bris.ac.uk/ws/portalfiles/portal/304220732/Thesis_Mary_Coe.pdf

Copyright (c) 2022 Mary Coe
Made available under the MIT License.

This file contains functions to delete particles within the 
simulation box.
*/

#include "delete.h"
#include <cmath>

void delete_particle(sim_box* sim_env_ptr, int pn, double energy, double Vext) {

	/* Deletes particle and reorganises cell and particle lists. */

		// Used within sub-volume sampling. If the particle is within a sub-volume
		// cell, decrease the number of sub-volume particles.
		if (sim_env_ptr->return_sub_volume() && sim_env_ptr->cell_list[sim_env_ptr->particle_list[pn].cell_id].sv) {
			sim_env_ptr->dec_sv_particles();
		} 

		sim_env_ptr->reorder_particles(pn);
		sim_env_ptr->dec_running_energy(energy); 
		if (sim_env_ptr->return_solute() || sim_env_ptr->return_slit()) sim_env_ptr->dec_running_Vext(Vext);
		sim_env_ptr->dec_nparticles(); // Decrease number of particles
}

inline double deletion_prob(sim_box* sim_env_ptr, double dE) {

	/* Calculate the probability of deleting the particle. */

	return ((sim_env_ptr->return_nparticles())/sim_env_ptr->return_volume())*exp((dE-sim_env_ptr->return_mu()));
}

void decide_deletion_mc(int pn, sim_box* sim_env_ptr) {

	/*
	Decides whether to delete a given particle when using biased sampling,
	and if so proceeds with deletion. 
	*/
	 
	double dE, prob, prob_unbiased, Vext = 0.0;

	// Calculate energy the particle contributes to the system.
	if (sim_env_ptr->return_mW()) dE = mWPotential(pn, sim_env_ptr);
	else  dE = LJPotential(pn, sim_env_ptr);

	// Calculate the external potential energy the particle contributes to the system.
	if (sim_env_ptr->return_LJ_Vext()) Vext = Vext_LJPotential(pn,sim_env_ptr);

	// Calculate the probability of deletion.
	prob = deletion_prob(sim_env_ptr, dE+Vext);

	// If the probability exceeds 1, set to 1.
	if (prob>1.) prob_unbiased = 1.;
	else prob_unbiased = prob;

	// Add bias to probability of deletion.
	// Bias for transition matrix method.
	if (sim_env_ptr->return_tmmc() && sim_env_ptr->cell_list[sim_env_ptr->particle_list[pn].cell_id].sv) {

		// Using sub-volume sampling.
		if (sim_env_ptr->return_sub_volume()) {
			sim_env_ptr->update_collection(sim_env_ptr->return_sv_particles(),0,prob_unbiased);
			sim_env_ptr->update_collection(sim_env_ptr->return_sv_particles(),1,1.-prob_unbiased);
			prob*=exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles())-sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()-1));
		}

		// Sampling entire system.
		else {
			sim_env_ptr->update_collection(sim_env_ptr->return_nparticles(),0,prob_unbiased);
			sim_env_ptr->update_collection(sim_env_ptr->return_nparticles(),1,1.-prob_unbiased);
			prob*=exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles())-sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()-1));
		}
	}
	// Bias for pre-set weights.
	else  if (!sim_env_ptr->return_tmmc()) prob*=exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles())-sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()-1));
	
	// Decide whether to accept the move and implement decision.
	if (prob>sim_env_ptr->get_rnd())  { 
		// Accept move.
		sim_env_ptr->inc_mca();
		delete_particle(sim_env_ptr, pn, dE, Vext);
	}
	else {
		// Reject move.
		sim_env_ptr->inc_mcr();
	}
}

void decide_deletion(int pn, sim_box* sim_env_ptr) {

	/*
	Decides whether to delete a given particle when using unbiased sampling,
	and if so proceeds with deletion.
	*/
	 
	double dE, prob, prob_unbiased, Vext = 0.0;

	// Calculate energy the particle contributes to the system.
	if (sim_env_ptr->return_mW()) dE = mWPotential(pn, sim_env_ptr);
	else dE = LJPotential(pn, sim_env_ptr);

	// Calculate the external potential energy the particle contributes to the system.
	if (sim_env_ptr->return_LJ_Vext()) Vext = Vext_LJPotential(pn,sim_env_ptr);

	// Calculate the probability of deletion.
	prob = deletion_prob(sim_env_ptr, dE+Vext);

	// Decide whether to accept the move and implement decision.
	if (prob>sim_env_ptr->get_rnd())  { 
		//Accept the move
		sim_env_ptr->inc_mca();
		delete_particle(sim_env_ptr, pn, dE, Vext);
	}
	else { 
		//Reject the move
		sim_env_ptr->inc_mcr();
	}
}

