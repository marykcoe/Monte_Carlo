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

This file contains functions to calculate energies within the
system.
*/

#include "energetics.h"

double LJPotential (int pn, sim_box* sim_env_ptr) {
	
	/* 
	 Calculates the energy that a single particle contributes to the 
	 system using the 6-12 truncated Lennard-Jones potential.
	*/	
	
	double LJ=0., r2, r6, r12;
	int i, j, k, pcell=sim_env_ptr->particle_list[pn].cell_id;
	cell mcell;
	particle p = sim_env_ptr->particle_list[pn];

	// Loops over all particles in cells neighbouring and including that
	// of the particle to find particles with which the particle interacts.
	// Calculates the distance between interacting particles, and uses this
	// to calculate the Lennard-Jones potential between the two.
	for(int n=0;n<sim_env_ptr->cell_list[pcell].num_neighbours;n++) {

		// Get cell
		mcell = sim_env_ptr->cell_list[sim_env_ptr->cell_list[pcell].neighbours[n][0]];

		// Loop over all particles within the cell
		for (int m=0;m<mcell.num_occupants;m++) {
			if(pn == mcell.occupancy[m]) continue; // Ignore the particle in question

			// Find the distance between the interacting particles
			r2 = Find_distance(p, sim_env_ptr->particle_list[mcell.occupancy[m]],
								sim_env_ptr->cell_list[pcell].neighbours[n][1],
								sim_env_ptr->cell_list[pcell].neighbours[n][2],
								sim_env_ptr->cell_list[pcell].neighbours[n][3]);
			
			// If the particles are within interacting distance, calculate 
			// their energy and add to the total energy the particle 
			// contributes to the system.
			if (r2<=1.) {
				r2=1./(r2*LJ_rc*LJ_rc);
				r6=r2*r2*r2; r12=r6*r6;
				LJ+=((LJ_sigma*r12)-(LJ_sigma*r6));
			}
		}
	} 
	return sim_env_ptr->return_epsilon()*LJ;
}

double mWPotential (int pn, sim_box* sim_env_ptr) {

	/* 
	Calculates the energy that a single monatomic water particle contributes to
	the system. The monatomic water potential is a modified Stillinger-Weber 
	potential - see V Molinero and E B Moore J. Phys. Chem. B 113 pp 4008-16 (2009).
	*/
	 
	 double rij2, r4, two_body=0, three_body=0, rik2, rjk2, costheta, neighbour_info[3][50], temp;
	 int ncounter=0, pcellid = sim_env_ptr->particle_list[pn].cell_id, ci, cj, ck, pnei;
	 cell ncell, pcell = sim_env_ptr->cell_list[sim_env_ptr->particle_list[pn].cell_id];
	 particle neiparticle, p = sim_env_ptr->particle_list[pn]; 
	 
	 // Loop over all particles in cells neighbouring and including the particle
	 // in question. Calculate the distance between each particle and the particle 
	 // in question. If the particles are within the interaction radius, add the
	 // two-body potential to the total two-body energy, and add the neighbouring
	 // particle to an array of neighbours.
	 for (int n=0;n<pcell.num_neighbours;n++) {
		ncell = sim_env_ptr->cell_list[pcell.neighbours[n][0]]; // neighbouring cell

		// Loop over all particles in the neighbouring cell
		for(int m=0;m<ncell.num_occupants;m++){
			if(pn == ncell.occupancy[m]) continue; // This is the particle in question so we ignore

			// Calculate the distance between the particles
			rij2 = Find_distance(p,sim_env_ptr->particle_list[ncell.occupancy[m]],pcell.neighbours[n][1],
																pcell.neighbours[n][2], pcell.neighbours[n][3]);
			
			// If particles are within interaction distance, add their two-body interaction potential
			// to the total two-body interaction energy.
			if (rij2<1.) {
				rij2*=mW_aa;
				r4=rij2*rij2;
				rij2 = sqrt(rij2);
				two_body+=(mW_B*(1./r4)-1.)*exp(1./(rij2-mW_a)); 
			
				// Add neighbouring particle information to neighbour info array
				neighbour_info[0][ncounter] = ncell.occupancy[m]; // position of neighbouring particle in particle list
				neighbour_info[1][ncounter] = n; // position of neighbouring particle cell in the cell neighbour list
				neighbour_info[2][ncounter] = rij2; // distance between particles
				ncounter++; // number of neighbouring particles
			}
		}
	}

	int uncounter = ncounter;

	// Loop over all particles in the neighbouring particle list and calculate the three-body contribution
	// to the energy.
	
	for(int m=0;m<ncounter;m++) {
		for(int n=m+1;n<ncounter;n++) {

			// Calculate the distance between cells of the two neighbouring particles.
			ci = pcell.neighbours[(int)neighbour_info[1][m]][1]-pcell.neighbours[(int)neighbour_info[1][n]][1];
			cj = pcell.neighbours[(int)neighbour_info[1][m]][2]-pcell.neighbours[(int)neighbour_info[1][n]][2];
			ck = pcell.neighbours[(int)neighbour_info[1][m]][3]-pcell.neighbours[(int)neighbour_info[1][n]][3]; 

			// Find the distance between the two neighbouring particles.
			rjk2 = Find_distance(sim_env_ptr->particle_list[(int)neighbour_info[0][n]],
					sim_env_ptr->particle_list[(int)neighbour_info[0][m]],ci,cj,ck, sim_env_ptr->return_Lx(),
					sim_env_ptr->return_Ly(), sim_env_ptr->return_Lz(), sim_env_ptr->return_slit());
			rjk2*=mW_aa;

			// Use distances to find the angle between the three particles.
			costheta = Find_angle(rjk2,neighbour_info[2][m],neighbour_info[2][n]);
			
			// Finally work out the three-body energy contributed by the particle due to its neighbouring
			// particles.
			three_body+=(costheta + (1./3.))*(costheta + (1./3.))*exp((mW_gamma/(neighbour_info[2][m]-mW_a))
											+(mW_gamma/(neighbour_info[2][n]-mW_a)));
		}
	}

	// Some particles which are not neighbours with the particle in question will make contributions
	// to the three-body potential due to interactions with intermediate particles - see the section
	// 4.2.3 of the thesis listed in the header for more details.

	// Loop over all neighbouring particles and find particles which fall within their interaction
	// range. 
	for(int wm =0;wm<ncounter;wm++) {
		neiparticle = sim_env_ptr->particle_list[(int)neighbour_info[0][wm]]; // neighbouring particle
		ncell = sim_env_ptr->cell_list[neiparticle.cell_id]; //neighbouring particle's cell
		
		//Find the particle in questions's cell in the neighbouring particles's cell's neighbour list
		for(int c=0;c<ncell.num_neighbours;c++) {
			if(ncell.neighbours[c][0] == p.cell_id) {pnei = c; break;}
		} 
		
		// Loop over all cells neighbouring the neighbouring particles's cell, including the 
		// cell of the neighbouring particle.
		for(int n=0;n<ncell.num_neighbours;n++) { 
			pcell = sim_env_ptr->cell_list[ncell.neighbours[n][0]]; // Neighbouring particle's neighbouring cell
			
			// Loop over all particles in the neighbouring particle's neighbouring cell
			for(int m=0;m<pcell.num_occupants;m++) {

				if(pn==pcell.occupancy[m]) continue; // ignore particle in question
				if(neighbour_info[0][wm] == pcell.occupancy[m]) continue; // ignore particles we've already found energy from
				
				// Otherwise calculate the distance between the neighbouring particle and its neighbour
				rjk2 = Find_distance(neiparticle, sim_env_ptr->particle_list[pcell.occupancy[m]], ncell.neighbours[n][1],
										ncell.neighbours[n][2],ncell.neighbours[n][3]);
				
				// If these two particles are within interacting range, calculate the three-body energy
				// contribution they make with the particle in question.
				if(rjk2<1.0) {			
					rjk2=sqrt(rjk2)*mW_a;
					ci = ncell.neighbours[n][1]-ncell.neighbours[pnei][1]; 
					cj = ncell.neighbours[n][2]-ncell.neighbours[pnei][2];
					ck = ncell.neighbours[n][3]-ncell.neighbours[pnei][3]; 
					rik2 = Find_distance(p,sim_env_ptr->particle_list[pcell.occupancy[m]],ci,cj,ck, 
						sim_env_ptr->return_Lx(), sim_env_ptr->return_Ly(), sim_env_ptr->return_Lz(),
						sim_env_ptr->return_slit());
					rik2*=mW_aa;
					costheta = Find_angle(rik2,rjk2,neighbour_info[2][wm]);
					
					three_body+=(costheta + (1./3.))*(costheta + (1./3.))*exp((mW_gamma/(rjk2-mW_a))+(mW_gamma/(neighbour_info[2][wm]-mW_a)));
				}
			}
		}
	}

	return sim_env_ptr->return_epsilon()*(mW_A*two_body + mW_lambda*three_body);
}

double LJEnergy (sim_box* sim_env_ptr) {
	
	/* 
	Calculates total energy of system due to particles interacting
	via a 6-12 truncated Lennard-Jones potential.
	*/
	double LJ=0;
	
	for(int m=0;m<sim_env_ptr->return_nparticles();m++) LJ+=LJPotential(m, sim_env_ptr);
	
	return LJ/2.;
}

double mWEnergy(sim_box* sim_env_ptr) {
	
	/* 
	Calculates total energy of system due to particles interacting via a monatomic
	water potential - see V Molinero and E B Moore J. Phys. Chem. B 113 pp 4008-16 (2009).
	*/
	
	double two_body=0, three_body=0, r2,r4, nnr2, costheta;
	std::deque<particle> neighbours;
	int  ci=0,cj=0,ck=0, ncell, mcell, neighbour_counter=0, pcell;
	cell rcell;
	particle p;
	std::deque<int> cell_neighbour_pos;
	std::deque<double>  ndistances;
	
	// Loop over every particle.
	for(int wm=0;wm<sim_env_ptr->return_nparticles();wm++) {
		p = sim_env_ptr->particle_list[wm]; // Find the particle in the particle list
		pcell = p.cell_id; // Find the particle's cell

		// Loop over all neighbouring cells of the particle's cell, as well as its own.
		for (int n=0;n<sim_env_ptr->cell_list[pcell].num_neighbours;n++) {
			rcell = sim_env_ptr->cell_list[sim_env_ptr->cell_list[pcell].neighbours[n][0]];

			// Loop over all particles within the cell.
			for(int m=0;m<rcell.num_occupants;m++){
				if(wm == rcell.occupancy[m]) continue; // Ignore the particle in question, p

				// Find the distance between the particle and its neighbour
				r2 = Find_distance(p,sim_env_ptr->particle_list[rcell.occupancy[m]],sim_env_ptr->cell_list[pcell].neighbours[n][1],
								sim_env_ptr->cell_list[pcell].neighbours[n][2],sim_env_ptr->cell_list[pcell].neighbours[n][3]);
				
				// If the particles are within interacting distance, work out the two-body energy contribution
				if (r2<1.0) {
					neighbours.push_back(sim_env_ptr->particle_list[rcell.occupancy[m]]);
					cell_neighbour_pos.push_back(n);
					ndistances.push_back(sqrt(r2)*mW_a);
					r2*=mW_aa;
					r4=r2*r2;
					two_body+=(mW_B*(1./r4)-1.)*exp(1./((sqrt(r2)-mW_a)));
				}
			}
		}
	
		// Loop over all combinations of neighbouring particles to calculate the
		// three energy contribution.
		for(int m=0;m<neighbours.size();m++) {
			for(int n=m+1;n<neighbours.size();n++) {

				// Find the distance between the two neighbouring particles.
				mcell = cell_neighbour_pos[m];
				ncell = cell_neighbour_pos[n];
			
				ci = sim_env_ptr->cell_list[pcell].neighbours[mcell][1]-sim_env_ptr->cell_list[pcell].neighbours[ncell][1];
				cj = sim_env_ptr->cell_list[pcell].neighbours[mcell][2]-sim_env_ptr->cell_list[pcell].neighbours[ncell][2];
				ck = sim_env_ptr->cell_list[pcell].neighbours[mcell][3]-sim_env_ptr->cell_list[pcell].neighbours[ncell][3]; 
			
				nnr2 = Find_distance(neighbours[n],neighbours[m],ci,cj,ck,sim_env_ptr->return_Lx(),
								sim_env_ptr->return_Ly(), sim_env_ptr->return_Lz(), sim_env_ptr->return_slit());
				nnr2*=mW_aa;

				// Find the angle and therefore three-body energy contribution.
				costheta = Find_angle(nnr2,ndistances[n],ndistances[m]);
				three_body+=(costheta + (1./3.))*(costheta + (1./3.))*exp((mW_gamma/(ndistances[m]-mW_a))+(mW_gamma/(ndistances[n]-mW_a)));
			}
		}
		while(!neighbours.empty()) {neighbours.pop_back();cell_neighbour_pos.pop_back();ndistances.pop_back();}
	}
	return sim_env_ptr->return_epsilon()*((0.5*mW_A*two_body)+(mW_lambda*three_body));
}

double Vext_LJPotential(int pn, sim_box* sim_env_ptr)  {

	/* 
	Returns the external potential felt by a particle due to a solute/slit.
	The interaction potential is taken to be a 6-12 Lennard-Jones potential,
	which has been integrated over degrees of freedom of homogeneous density
	and shifted such that the minimum occurs at the surface of the solute or
	at the walls of the slit - see thesis in header for more information. Each 
	wall of the slit is assumed to have the same interaction potential. In all
	cases, the diameter of solute/slit-wall particles is taken to each that of
	a fluid particle.
	*/

	particle p = sim_env_ptr->particle_list[pn];
	cell pcell = sim_env_ptr->cell_list[p.cell_id];
	double V = 0.0, z = 0.0, zmax = sim_env_ptr->return_Lz(), mult, shift = pow(0.4,1./6.);
	double dx=0.0, dy=0.0, dz=0.0;

	if (sim_env_ptr->return_mW()) mult = mW_a;
	else mult = LJ_rc;

	zmax*=mult;

	// External potential due to slit geometry. Each wall contributes the
	// same potential which differs only from origin point.
	if (sim_env_ptr->return_slit()) {
		z = mult*(pcell.k + p.z); // distance of particle from left wall
		V = (2./15.)*(1./pow(z+shift,9.)) - 1./pow(z+shift,3.); // Potential due to left wall
		V += (2./15.)*(1/pow(zmax+shift-z,9.)) - 1./pow(zmax-z + shift,3.); // Potential due to right wall
	}

	// Solute interaction potential is more complicated however has only
	// one origin - the centre of the solute particle.
	else {
		dx = sim_env_ptr->vsolute.x - (pcell.i + p.x);
		dy = sim_env_ptr->vsolute.y - (pcell.j + p.y);
		dz = sim_env_ptr->vsolute.z - (pcell.k + p.z);

		z = dx*dx + dy*dy + dz*dz; z*=mult*mult; z = sqrt(z); // radial distance of particle from solute centre
		z+= shift; // We use planar potential shift as shift in spherical geometry cannot be calculated analytically

		// Potential due to solute
		V = (2./15.)*(1./pow(z- sim_env_ptr->vsolute.Rs,9.) -1./pow(z  + sim_env_ptr->vsolute.Rs,9.));
		V += (1./pow(z  + sim_env_ptr->vsolute.Rs,3.) - 1./pow(z  - sim_env_ptr->vsolute.Rs,3.));
		V += (3./20.) *(1./z) * (1./pow(z  + sim_env_ptr->vsolute.Rs,8.) - 1./pow(z - sim_env_ptr->vsolute.Rs,8.));
		V += 1.5 * (1./z) *  (1./pow(z  - sim_env_ptr->vsolute.Rs,2.) - 1./pow(z + sim_env_ptr->vsolute.Rs,2.));
	}

	return sim_env_ptr->return_epsilon()*sim_env_ptr->return_epsilon_wall()*V;
}

double Vext_LJ(sim_box* sim_env_ptr) {
	
	/*
	Calculates external potential energy of system due to all particles.
	*/

	double V = 0.0;

	for(int p=0; p<sim_env_ptr->return_nparticles(); p++) V+=Vext_LJPotential(p,sim_env_ptr);

	return V;
}

