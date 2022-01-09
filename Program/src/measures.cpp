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

This file contains functions which perform measures on and
within the system.
*/

#include "measures.h"

double pair_dist(double xa, double xb, int cx) {

	/* 
	Finds the distance in one dimension between a pair of particles
	which are in neighbouring cells.
	*/

    double r;
    if (cx == 0) r = xb-xa; // Case when particle's cell displacement is 0
    else r = 1.+cx*(xb-xa); // Case when particle's cell displacement is not 0
    return r;
}

double three_body_dist(double xa, double xb, int cx, int L) {
	
	/* 
	Finds the distance in one dimension between a pair of particles 
	which may not be in neighbouring cells. Accounts for periodic 
	boundary conditions.
	*/

	int cxx=0;
	double r;

	// Calculates the distance if the particle's cell displacement
	// is greater than 0. 
	if (cx>0) {
		if (cx>1) cxx = cx-1; // Works out how many full cells between particles.
		r = 1. + (xb - xa); // Distance between particles if they were in neighbouring cells.
		if (cxx>=(int)ceil(L/2)-1) {
			if ((1-(xb - xa))<r) r = 1-(xb - xa); // Account for periodic boundary conditions.
		}
	}
	
	// Calculates the distance if the particle's cell displacement
	// is less than 0. 
	else if (cx<0) {
		if(cx<-1) cxx = -1*(cx+1); // Works out how many full cells between particles.
		r= 1-(xb - xa); // Distance between particles if they were in neighbouring cells.	
		if (cxx>=(int)ceil(L/2)-1) {
			if ((1+(xb - xa))<r) r = 1+(xb - xa); // Account for periodic boundary conditions.
		}
	}

	// Calculates the distance if the particle's cell displacement
	// is 0. 
	else r = xb - xa; 

	return (r + cxx); 
}

double three_body_dist(double xa, double xb, int cx) {

	/* 
	Finds the distance in one dimension between a pair of particles 
	which may not be in neighbouring cells. Does not Account for 
	periodic boundary conditions.
	*/

	int cxx=0;
	double r;

	// Calculates the distance if the particle's cell displacement
	// is greater than 0. 
	if (cx>0) {
		if (cx>1) cxx = cx-1; // Works out how many full cells between particles.
		r = 1. + (xb - xa); // Distance between particles if they were in neighbouring cells.
	}

	// Calculates the distance if the particle's cell displacement
	// is less than 0. 
	else if (cx<0) {
		if(cx<-1) cxx = -1*(cx+1); // Works out how many full cells between particles.
		r = 1-(xb - xa); // Distance between particles if they were in neighbouring cells.
	}

	// Calculates the distance if the particle's cell displacement
	// is equal to 0. 
	else r = xb - xa; 

	return (r + cxx);
}

double Find_distance(particle pa, particle pb, int i, int j, int k) {
	
	/*
	Calculates the distance between a pair of particles in neighbouring cells
	using Pythagoras rule. Returns squared distance.
	*/
	
	double r2,ri,rj,rk;
	
    ri = pair_dist(pa.x,pb.x,i);
    rj = pair_dist(pa.y,pb.y,j);
    rk = pair_dist(pa.z,pb.z,k);
	
	r2 = ri*ri+rj*rj+rk*rk;
	return r2;
}

double Find_distance(particle pa, particle pb, int i, int j, int k, int Lx, int Ly, int Lz, bool slit) {
	
	/*
	Calculates the distance between a pair of particles which may not be in neighbouring cells
	using Pythagoras rule. Returns squared distance.
	*/
	 
	double r2, ri,rj,rk;

	ri = three_body_dist(pa.x,pb.x,i,Lx);
	rj = three_body_dist(pa.y,pb.y,j,Ly);
	if (slit) rk = three_body_dist(pa.z,pb.z,k); // Do not account for periodic boundary conditions.
	else rk = three_body_dist(pa.z,pb.z,k,Lz); // Account for periodic boundary conditions.
	
	r2 = (ri*ri) + (rj*rj) + (rk*rk);
	
	return r2;
}

/* Calculates cosine of angle between three particles using cosine rule. */
double Find_angle(double ar2, double br, double cr) { return ((cr*cr + br*br - ar2)/(2*br*cr)); }

void gr(sim_box* sim_env_ptr, int sweep) {
	
	/* 
	Finds the radial distribution function (g(r)) for every particule within the system.
	The result is added to the average g(r) for the system and if applicable, is output
	to file.
	*/
	 
	cell mcell, ncell;
	int ci, cj, ck, m, mm, b;
	double Lx = sim_env_ptr->return_Lx(), Ly = sim_env_ptr->return_Ly(), Lz = sim_env_ptr->return_Lz();
	double r2, mult, dx, dy, dz;
	double x0, y0, z0, x, y, z;
	double hist[sim_env_ptr->return_profile_bins()];

	std::ofstream grout;

	// Calculate real-space distance multipliers.
	// This allows distribution to be output in real rather than
	// reduced units.
	if (sim_env_ptr->return_mW()) mult = mW_a*mW_sigma;
	else mult = LJ_rc*LJ_sigma;
	
	// Prepare empty array for histogram.
	for (m=0; m<sim_env_ptr->return_profile_bins(); m++) hist[m] = 0.0;

	// Calculate the distance between every particle within the system.
	// Add this to the correct histogram bin.
	for (m=0; m<sim_env_ptr->return_nparticles()-1; m++) {
		
		// Calculate location of origin particle in real space.
		mcell = sim_env_ptr->cell_list[sim_env_ptr->particle_list[m].cell_id];
		x0 = mult*(sim_env_ptr->particle_list[m].x + mcell.i);
		y0 = mult*(sim_env_ptr->particle_list[m].y + mcell.j);
		z0 = mult*(sim_env_ptr->particle_list[m].z + mcell.k);

		for (mm = m+1; mm<sim_env_ptr->return_nparticles(); mm++) {

			// Calculate location of comparison particle in real space.
			ncell = sim_env_ptr->cell_list[sim_env_ptr->particle_list[mm].cell_id];
			x = mult*(sim_env_ptr->particle_list[mm].x + ncell.i);
			y = mult*(sim_env_ptr->particle_list[mm].y + ncell.j);
			z = mult*(sim_env_ptr->particle_list[mm].z + ncell.k);
			 
			// Find distance between the two particles, accounting for periodic
			// boundary conditions.
			dx = abs(x-x0); if (dx>Lx/2.) dx = Lx-dx;
			dy = abs(y-y0); if (dy>Ly/2.) dy = Ly-dy; 
			dz = abs(z-z0); if (dz>Lz/2.) dz = Lz-dz;
			r2 = dx*dx + dy*dy + dz*dz;
			r2 = sqrt(r2);

			// Find histogram bin distance falls into, and increment tally of bin.
			for(b=0;b<sim_env_ptr->return_profile_bins()-1;b++){
				if (r2>=sim_env_ptr->average_distribution[b][0] && r2<sim_env_ptr->average_distribution[b+1][0]) {hist[b] +=2; break;}
				else if (r2>sim_env_ptr->average_distribution[sim_env_ptr->return_profile_bins()-1][0]) {hist[sim_env_ptr->return_profile_bins()-1] +=2; break;}
			}
		}
	}

	// Add tally of distances to running total for averaging later.
	for (b=0; b<sim_env_ptr->return_profile_bins(); b++) 	{
		if (sim_env_ptr->return_mc()) {
			if (sim_env_ptr->return_sub_volume()) sim_env_ptr->average_distribution[b][1] += exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()))*hist[b];
			else sim_env_ptr->average_distribution[b][1] += exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()))*hist[b];
		}
		else sim_env_ptr->average_distribution[b][1] += hist[b];
	}
	
	// If the simulation uses multicanonical bias, add the appropriate weight to the
	// sum of distributions. If the simulation is unbiased, add the number of particles,
	// and therefore number of distributions calculated, to the sum of distributions.
	if (sim_env_ptr->return_mc()) {
		if (sim_env_ptr->return_sub_volume()) sim_env_ptr->sum_distributions += sim_env_ptr->return_nparticles()*exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()));
		else sim_env_ptr->sum_distributions += sim_env_ptr->return_nparticles()*exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()));
	}
	else sim_env_ptr->sum_distributions += sim_env_ptr->return_nparticles();


	// Write the collected distributions to file if applicable.
	if (sim_env_ptr->return_output_distributions()) {
		grout.open(sim_env_ptr->return_output_file_path()/"distributions"/("gr_" + std::to_string(sweep)), std::ios::app);
		
		for (m=0; m< sim_env_ptr->return_profile_bins(); m++) grout<<sim_env_ptr->average_distribution[m][0]<<" "<<hist[m]<<std::endl;

		grout<<"Particles = "<<sim_env_ptr->return_nparticles()<<std::endl;
		if (sim_env_ptr->return_sub_volume()) grout<<"Sub-Volume Particles = "<<sim_env_ptr->return_sv_particles()<<std::endl;

		if (sim_env_ptr->return_mc()) {
			if (sim_env_ptr->return_sub_volume()) grout<<"Weight = "<<sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles())<<std::endl;
			else grout<<"Weight = "<<sim_env_ptr->return_weight(sim_env_ptr->return_nparticles())<<std::endl;
		}
		grout.close();
	}

    // Calculate local compressibility and thermal susceptibility if applicable.
    if (sim_env_ptr->return_compressibility()) local_compressibility(hist, sim_env_ptr);
    if (sim_env_ptr->return_susceptibility()) local_susceptibility(hist, sim_env_ptr);
}

void planar_distribution(sim_box* sim_env_ptr, int sweep){

	/* 
	Calculates the spatial distribution of particles along the z-dimension within a slit.
	The result is added to the average distribution for the system and if applicable, is
	output to file. These results can then be used to calculate the density profile of 
	the system.
	*/

	double maxz = sim_env_ptr->return_Lz(), mult, binsize, z;
	double hist[sim_env_ptr->return_profile_bins()];
	int count=0, b;
	std::ofstream distribution;

	// Calculate real-space distance multipliers.
	// This allows distribution to be output in real rather than
	// reduced units.
	if (sim_env_ptr->return_mW()) mult = mW_a*mW_sigma;
	else mult = LJ_rc*LJ_sigma;

	// Prepare empty array for histogram.
	for(b = 0; b<sim_env_ptr->return_profile_bins(); b++) hist[b] = 0.0;

	// For this distribution, the system is split into slabs in the x-y plane of
	// equal width. The z-position of each particle is calculated, and the slab
	// it falls within identified. The histogram tally of this slab is then increased.
	for (int p = 0; p < sim_env_ptr->return_nparticles(); p++) {

		// z-position in real-space of particle.
		z = sim_env_ptr->cell_list[sim_env_ptr->particle_list[p].cell_id].k + sim_env_ptr->particle_list[p].z;
		z *= mult;

		// Slab particle falls within is found here, and the histogram tally updated.
		for (b = 0; b < sim_env_ptr->return_profile_bins()-1; b++) {
			if (z>sim_env_ptr->average_distribution[b][0] and z<sim_env_ptr->average_distribution[b+1][0]) {hist[b]+=1; count +=1; break;}
			else if (z>sim_env_ptr->average_distribution[sim_env_ptr->return_profile_bins()-1][0]) {hist[sim_env_ptr->return_profile_bins()-1]+=1; count +=1; break;}
		}
	}

	// Add tally of distances to running total for averaging later.
	for (b=0; b<sim_env_ptr->return_profile_bins(); b++) 	{
		if (sim_env_ptr->return_mc()) {
			if (sim_env_ptr->return_sub_volume()) sim_env_ptr->average_distribution[b][1] += exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()))*hist[b];
			else sim_env_ptr->average_distribution[b][1] += exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()))*hist[b];
		}
		else sim_env_ptr->average_distribution[b][1] += hist[b];
	}

	// If the simulation uses multicanonical bias, add the appropriate weight to the
	// sum of distributions. If the simulation is unbiased, add one to the sum of distributions.
	if (sim_env_ptr->return_mc()) {
		if (sim_env_ptr->return_sub_volume()) sim_env_ptr->sum_distributions += exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()));
		else sim_env_ptr->sum_distributions += exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()));
	}
	else sim_env_ptr->sum_distributions += 1.;

	// Write the collected distribution to file if applicable.
	if (sim_env_ptr->return_output_distributions()) {

		distribution.open(sim_env_ptr->return_output_file_path()/"distributions"/("z_distribution_" + std::to_string(sweep)),std::ios::out);

		for(int b = 0; b<sim_env_ptr->return_profile_bins(); b++) distribution<<sim_env_ptr->average_distribution[b][0]<<" "<<hist[b]<<std::endl;

		// The total number of particles and multicanonical weight, if applicable are also
		// written to file, to allow for post-simulation histogram reweighting.
		distribution<<"Particles = "<<sim_env_ptr->return_nparticles()<<std::endl;
		if (sim_env_ptr->return_sub_volume()) distribution<<"Sub-Volume Particles = "<<sim_env_ptr->return_sv_particles()<<std::endl;

		if (sim_env_ptr->return_mc()) {
			if (sim_env_ptr->return_sub_volume()) distribution<<"Weight = "<<sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles())<<std::endl;
			else distribution<<"Weight = "<<sim_env_ptr->return_weight(sim_env_ptr->return_nparticles())<<std::endl;
		}
		distribution.close();
	}

	// Calculate local compressibility and thermal susceptibility if applicable.
	if (sim_env_ptr->return_compressibility()) local_compressibility(hist, sim_env_ptr);
	if (sim_env_ptr->return_susceptibility()) local_susceptibility(hist, sim_env_ptr);

}

void solute_rdf(sim_box* sim_env_ptr, int sweep) {

	/* 
	Calculates the spatial distribution of particles along the radial axis extending from 
	the centre of the solute. The result is added to the average distribution for the system
	and if applicable, is output to file. These results can then be used to calculate the 
	density profile of the system.
	*/

	double length_mult, hist[sim_env_ptr->return_profile_bins()], r;
	double dx,dy,dz;
	int  total, b;
	particle part;
	cell pcell;
	std::ofstream distribution;

	// Calculate real-space distance multipliers.
	// This allows distribution to be output in real rather than
	// reduced units.
	if (sim_env_ptr->return_mW()) length_mult = mW_a*mW_sigma;
	else length_mult = LJ_rc*LJ_sigma;

	// Prepare empty array for histogram.
	for(b = 0; b<sim_env_ptr->return_profile_bins(); b++) hist[b] = 0.0;

	// Loop over all particles and find their distance from the centre of
	// the solute. The system is split into spherical shells from the 
	// centre of the solute. Find the shell the particle falls within 
	// and add to the histogram tally.
	for (int p = 0; p<sim_env_ptr->return_nparticles(); p++) {

		// Particle and cell particle falls within.
		part = sim_env_ptr->particle_list[p];
		pcell = sim_env_ptr->cell_list[part.cell_id];

		// Calculate distance from centre of solute in real-space
		// using Pythagoras.
		dx = (pcell.i + part.x - sim_env_ptr->vsolute.x);
		dy = (pcell.j + part.y - sim_env_ptr->vsolute.y);
		dz = (pcell.k + part.z - sim_env_ptr->vsolute.z);
		r = dx*dx + dy*dy + dz*dz;
		r = length_mult*sqrt(r);

		// Find histogram shell particle falls within and add one to the tally.
		for (b = 0; b<sim_env_ptr->return_profile_bins()-1; b++) {
			if (r>sim_env_ptr->average_distribution[b][0] and r<sim_env_ptr->average_distribution[b+1][0]) {hist[b]+=1; break;}
			else if (r>sim_env_ptr->average_distribution[sim_env_ptr->return_profile_bins()-1][0]) {hist[sim_env_ptr->return_profile_bins()-1]+=1; break;}
		}
	}

	// Write the collected distribution to file if applicable.
	if (sim_env_ptr->return_output_distributions()) {

		distribution.open(sim_env_ptr->return_output_file_path()/"distributions"/("solute_rdf_" + std::to_string(sweep)),std::ios::out);
		
		for(b = 0; b<sim_env_ptr->return_profile_bins(); b++)  distribution<<sim_env_ptr->average_distribution[b][0]<<" "<<hist[b]<<std::endl;

		// The total number of particles and multicanonical weight, if applicable are also
		// written to file, to allow for post-simulation histogram reweighting.
		distribution<<"Particles = "<<sim_env_ptr->return_nparticles()<<std::endl;
		if (sim_env_ptr->return_sub_volume()) distribution<<"Sub-Volume Particles = "<<sim_env_ptr->return_sv_particles()<<std::endl;

		if (sim_env_ptr->return_mc()) {
			if (sim_env_ptr->return_sub_volume()) distribution<<"Weight = "<<sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles())<<std::endl;
			else distribution<<"Weight = "<<sim_env_ptr->return_weight(sim_env_ptr->return_nparticles())<<std::endl;
		}
		distribution.close();
	}

	// Add tally of distances to running total for averaging later.
	for(b=0; b<sim_env_ptr->return_profile_bins(); b++) {
		if (sim_env_ptr->return_mc()) {
			if (sim_env_ptr->return_sub_volume()) sim_env_ptr->average_distribution[b][1] += exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()))*hist[b];
			else sim_env_ptr->average_distribution[b][1] += exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()))*hist[b];
		}
		else sim_env_ptr->average_distribution[b][1] += hist[b];
	}
	if (sim_env_ptr->return_mc()) {
		if (sim_env_ptr->return_sub_volume()) sim_env_ptr->sum_distributions += exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()));
		else sim_env_ptr->sum_distributions += exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()));
	}
	else sim_env_ptr->sum_distributions += 1.;

	// Calculate local compressibility and thermal susceptibility if applicable.
	if (sim_env_ptr->return_compressibility()) local_compressibility(hist, sim_env_ptr);
	if (sim_env_ptr->return_susceptibility()) local_susceptibility(hist, sim_env_ptr);

}

void local_compressibility(double profile[], sim_box* sim_env_ptr) {

	/*
	Using histogram reweighting in the chemical potential, finds the appropriate spatial
	distribution of particles at a slightly higher chemical potential than that of the
	simulation. Post-simulation, this can be used in conjunction with the spatial
	distribution at the simulation chemical potential to find the local compressibility.
	See T. Eckert, N.C.X. Stuhlmuller, F Sammuller and M. Schidt, Phys Rev. Lett. 125, 268004 (2020)
	for definition of the local compressibility. For information on histogram reweighting,
	see the thesis listed in the header, or A.H. Ferrenberg and R.H. Swendsen, Phys. Rev. Lett.
	61 2635-8 (1988).
	*/

	// Calculate the histogram reweighting weight.
	double weight = exp(sim_env_ptr->return_nparticles()*sim_env_ptr->return_dmu());

	// If using biased sampling, multiply this weight by the bias weight.
	if (sim_env_ptr->return_mc()) {
		if (sim_env_ptr->return_sub_volume()) weight*= exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()));
		else weight *= exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()));
	}

	// Multiply each histogram bin by the histogram weight. 
	for(int i=0;i<sim_env_ptr->return_profile_bins();i++) {
		sim_env_ptr->compressibility_profile[i] += weight*profile[i];
	}

	// Add the appropriate weight to the sum of weights for averaging later.
	if (sim_env_ptr->return_slit() or sim_env_ptr->return_solute()) sim_env_ptr->sum_compressibility += weight;
    else sim_env_ptr->sum_compressibility += sim_env_ptr->return_nparticles()*weight;
}

void local_susceptibility(double profile[], sim_box* sim_env_ptr) {

	/*
	Using histogram reweighting in the temperature, finds the appropriate spatial
	distribution of particles at a slightly higher temperature than that of the
	simulation. Post-simulation, this can be used in conjunction with the spatial
	distribution at the simulation temperature to find the local thermal susceptibility.
	See T. Eckert, N.C.X. Stuhlmuller, F Sammuller and M. Schidt, Phys Rev. Lett. 125, 268004 (2020)
	for definition of the local thermal susceptibility For information on histogram reweighting,
	see the thesis listed in the header, or A.H. Ferrenberg and R.H. Swendsen, Phys. Rev. Lett.
	61 2635-8 (1988).
	*/

	// Calculate the histogram reweighting weight.
	double weight = exp((sim_env_ptr->return_running_energy()+sim_env_ptr->return_running_Vext())*sim_env_ptr->return_dT());

	// If using biased sampling, multiply this weight by the bias weight.
	if (sim_env_ptr->return_mc()) {
		if (sim_env_ptr->return_sub_volume()) weight*= exp(sim_env_ptr->return_weight(sim_env_ptr->return_sv_particles()));
		else weight *= exp(sim_env_ptr->return_weight(sim_env_ptr->return_nparticles()));
	}

	// Multiply each histogram bin by the histogram weight. 
	for(int i=0;i<sim_env_ptr->return_profile_bins();i++) {
		sim_env_ptr->susceptibility_profile[i] += weight*profile[i];
	}

	// Add the appropriate weight to the sum of weights for averaging later.
	if (sim_env_ptr->return_slit() or sim_env_ptr->return_solute()) sim_env_ptr->sum_susceptibility += weight;
	else sim_env_ptr->sum_susceptibility += sim_env_ptr->return_nparticles()*weight;
}
