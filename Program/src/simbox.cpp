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

#include "simbox.h"

sim_box::sim_box() {
	
	/*
	Constructor for sim_box.
	*/

	load = false; slit = false; bsolute = false, equilibrate = true;
	compressibility = false; susceptibility = false, positions = false; 
	mc = false; tmmc = false, load_weights = false, sub_volume = false;
	vext = false; avg_acc = 0.0; acc_entries = 0;

	mca = 0; mcr = 0; rand_check=0;
	minimum_density = 0.0; maximum_density = 1.0; initial_density = 0.0;
	dmu = 0.0; dT = 0.0; running_energy = 0.0; running_Vext = 0.0;
	profile_bins = 200; surface_interaction = 0.0;
}

// Functions to set initial parameters from Input or Output files.
void sim_box::set_output_file_path(std::string file_path) {output_file_path = file_path;}
void sim_box::set_load_file_path(std::string file_path) {load_file_path = file_path;}
void sim_box::set_load(bool iload) {load = iload;}
void sim_box::set_positions(bool ipos) {positions = ipos;}
void sim_box::set_Lx(int iLx) {Lx = iLx;}
void sim_box::set_Ly(int iLy) {Ly = iLy;}
void sim_box::set_Lz(int iLz) {Lz = iLz;}
void sim_box::set_slit(bool islit) {slit = islit;}
void sim_box::set_solute(bool isolute) {bsolute = isolute;}
void sim_box::set_solute_radius(double iRs) {Rs = iRs;}
void sim_box::set_mW(bool ift) {mW = ift;}
void sim_box::set_temperature(double temp) {temperature = temp;}
void sim_box::set_chemical_potential(double imu) {mu = imu;}
void sim_box::set_epsilon(double ieps) {epsilon = ieps;}
void sim_box::set_minimum_density(double dens) {minimum_density = dens;}
void sim_box::set_maximum_density(double dens) {maximum_density = dens;}
void sim_box::set_initial_density(double dens) {initial_density = dens;}
void sim_box::set_sweeps(long isweeps) {sweeps = isweeps;}
void sim_box::set_equilibrate(bool ieq) {equilibrate = ieq;}
void sim_box::set_save_frequency(int freq) {save_frequency = freq;}
void sim_box::set_distribution_frequency(int freq) {distribution_frequency = freq;}
void sim_box::set_output_distributions(bool oup) {output_individual_distributions = oup;}
void sim_box::set_profile_bins(int npb) {profile_bins = npb;}
void sim_box::set_compressibility(bool comp) {compressibility = comp;}
void sim_box::set_susceptibility(bool susc) {susceptibility = susc;}
void sim_box::set_mc(bool imc) {mc = imc;}
void sim_box::set_load_weights(bool iload) {load_weights = iload;}
void sim_box::set_tmmc(bool itmmc) {tmmc = itmmc;}
void sim_box::set_sub_volume(bool isub_v) {sub_volume = isub_v;} 
void sim_box::set_sub_volume_dz(int idz) {sub_volume_dz = idz;}
void sim_box::set_vext(bool ivext) {vext = ivext;}
void sim_box::set_surface_interaction_strength(double ieps) {surface_interaction = ieps;}
void sim_box::set_test(bool itest) {test = itest;}
void sim_box::set_test_type(std::string itest_type) {test_type = itest_type;}
void sim_box::set_seed(long iseed) {seed = iseed;}
void sim_box::set_mca(long imca) {mca = imca;}
void sim_box::set_mcr(long imcr) {mcr = imcr;}
void sim_box::set_total_sweeps(long itsweeps) {total_sweeps = itsweeps;}

void sim_box::print_simulation_parameters(void) {

	/*
	Prints simulation parameters to the screen.
	*/

	std::cout<<std::endl<<"=========================================================="<<std::endl;
	std::cout<<"Simulation Parameters:"<<std::endl;

	std::cout<<std::endl<<"Fluid Parameters: "<<std::endl;
	std::cout<<"Fluid Potential: ";
	if (mW) std::cout<<"Monatomic Water"<<std::endl;
	else std::cout<<"Truncated Lennard-Jones"<<std::endl;
	std::cout<<"Chemical Potential: "<<std::setprecision(8)<<mu<<" Temperature: "<<temperature<<" Fluid Interaction Stength: "<<std::setprecision(8)<<epsilon<<std::endl;
	std::cout<<"Initial Density: "<<initial_density<<std::endl;
	std::cout<<"Maximum Density: "<<maximum_density<<std::endl;
	std::cout<<"Minimum Density: "<<minimum_density<<std::endl;

	std::cout<<std::endl<<"Box Parameters: "<<std::endl;
	std::cout<<"Box Dimenstions : "<<Lx<<"x"<<Ly<<"x"<<Lz;
	if (mW) std::cout<<" (a*sigma)^3"<<std::endl;
	else std::cout<<" (rc*sigma)^3"<<std::endl;
	std::cout<<"Slit: ";
	if (slit) std::cout<<"True"<<std::endl;
	else std::cout<<"False"<<std::endl;
	std::cout<<"Solute: ";
	if (bsolute) std::cout<<"Rs = "<<vsolute.Rs<<" sigma"<<std::endl;
	else std::cout<<"False"<<std::endl;

	std::cout<<std::endl<<"Sampling Parameters: "<<std::endl;
	std::cout<<"Multicanonical: ";
	if (mc) std::cout<<"True"<<std::endl;
	else std::cout<<"False"<<std::endl;
	std::cout<<"Transition Matrix: ";
	if (tmmc) std::cout<<"True"<<std::endl;
	else std::cout<<"False"<<std::endl;
	std::cout<<"Sub-Volume: ";
	if (sub_volume) std::cout<<sub_volume_dz<<" sigma"<<std::endl;
	else std::cout<<"False"<<std::endl;
	std::cout<<"Load Weight Function: ";
	if (load_weights) std::cout<<"True"<<std::endl;
	else std::cout<<"False"<<std::endl;
	std::cout<<"Save Frequency: "<<save_frequency<<" sweeps"<<std::endl;
	std::cout<<"Output Distribution Frequency: "<<distribution_frequency<<" sweeps"<<std::endl;
	std::cout<<"Number of Bins for Profile: "<<profile_bins<<std::endl;
	std::cout<<"Sweeps: "<<sweeps<<std::endl;
	std::cout<<"Equilibrate: ";
	if (equilibrate) std::cout<<"True"<<std::endl;
	else std::cout<<"False"<<std::endl;
	std::cout<<"Compressibility: ";
	if (compressibility) {
		std::cout<<"True"<<std::endl;
		std::cout<<"dmu = "<<dmu<<std::endl;
	}
	else std::cout<<"False"<<std::endl;
	std::cout<<"Susceptibility: ";
	if (susceptibility) {
		std::cout<<"True"<<std::endl;
		std::cout<<"dT = "<<dT<<std::endl;
	}
	else std::cout<<"False"<<std::endl;

	std::cout<<std::endl<<"File Parameters: "<<std::endl;
	std::cout<<"Output File Path: "<<output_file_path<<std::endl;
	if (load) std::cout<<"Load: True"<<std::endl<<"Load File Path: "<<load_file_path<<std::endl;
	else std::cout<<"Load: False"<<std::endl;

	std::cout<<std::endl<<"External Potential Parameters: "<<std::endl;
	std::cout<<"External Potential: ";
	if(vext) std::cout<<"Lennard-Jones"<<std::endl;
	else std::cout<<"None"<<std::endl;
	std::cout<<"Surface Interaction Strength: "<<surface_interaction<<" * fluid interaction strength"<<std::endl;

	std::cout<<std::endl<<"Seed: "<<seed<<std::endl;
	std::cout<<"=========================================================="<<std::endl;
}

void sim_box::setup(void) {
	
	/*
	Sets-up the simulation box according to the simulation parameters
	supplied in the Input and Output files. This includes organising
	the cell list, placing initial particles, adding a solute/slit and
	more. Simultaneously checks that the parameters supplied for the
	simulation are reasonable.
	*/

	int sv_lower, sv_upper, sv, nnbrs, ci, cj ,ck, id=0;
	bool upper = false;
	double cut_off, binsize, sigma, maxr;
	std::ofstream out;
	int correct_neighb=0;

	// Check geometry of box is sensible.
	if (slit && bsolute) {
		std::cerr<<"Slit with solute is not supported. Please set"<<std::endl;
		std::cerr<<"one of these parameters to FALSE."<<std::endl;
		if (load) std::cerr<<"Please check Output file."<<std::endl;
		else std::cerr<<"Please check Input file."<<std::endl;
		exit(0);
	}

	// Check box dimensions are sensible and calculate number of cells.
	if (Lx<=0 || Ly <=0 || Lz<=0) {
		std::cerr<<"Box dimensions were set to "<<Lx<<"x"<<Ly<<"x"<<Lz<<std::endl;
		std::cerr<<"Box dimensions must be greater than 0."<<std::endl;
		if (load) std::cerr<<"Please check Output file."<<std::endl;
		else std::cerr<<"Please check Input file."<<std::endl;
		exit(0);
	}
	else num_cells = Lx*Ly*Lz; 

	// Calculate the reduced energy from the supplied temperature.
	// This is dependent on fluid type.
	if (mW) {
		cut_off = mW_a; sigma = mW_sigma;
		epsilon = 6.189 * 4184./(6.02214076e23);
		epsilon /=(temperature*1.38064852e-23);
	}
	else {
		cut_off = LJ_rc; sigma = LJ_sigma;
		epsilon = 4./temperature;
	}

	volume = num_cells*pow(cut_off,3); // Calculate volume of simulation box.

	// If sub-volume sampling is enabled, calculate the upper and lower 
	// bounds for the cells in the z-direction.
	sv_lower = 0; sv_upper = Lz-1;
	if (sub_volume) {
		while ((sv_upper-sv_lower) >= sub_volume_dz) {
			if (upper) {sv_upper--; upper = false;}
			else {sv_lower++; upper = true;}
		}
	}

	// Create the cell list. Loop over the number of cells specified
	// in each dimension and add a cell with the co-ordinates to the
	// cell list.
	for (int i=0;i<Lx;i++) {
		for (int j=0;j<Ly;j++) {
			for (int k=0;k<Lz;k++) {
				// If periodic boundary conditions are enabled on all
				// sides of the box, then the number of neighbours 
				// each cell with have is 27.
				if (!slit) nnbrs = 27;

				// If using the slit geometry, periodic boundary
				// conditions must be removed on the x-y planes
				// at z=0 and z=Lz. Cells with one side at either
				// of these boundaries will have only 18 neighbours.
				else {
					if (k>0 and k<(Lz-1)) nnbrs = 27;
					else nnbrs = 18;
				}

				// If using sub-volume sampling, specify whether the cell is within the
				// sampling region.
				if (k>=sv_lower && k<=sv_upper) cell_list.push_back(cell(i,j,k,nnbrs,1));
				else cell_list.push_back(cell(i,j,k,nnbrs,0)); 
				id++;
			}
		}
	}
	
	// Find the neighbours of each cell, accounting for periodic boundary conditions.
	for(int c =0; c<num_cells; c++) {
		nnbrs=0;

		// Loop over all cells.
		for(int cc=0; cc<num_cells; cc++) {

			// If the cell positions in the cell list match, add
			// the cell to its own neighbour list.
			if (cc==c) {
				cell_list[c].neighbours[nnbrs][0]=cc;
				cell_list[c].neighbours[nnbrs][1]=0;
				cell_list[c].neighbours[nnbrs][2]=0;
				cell_list[c].neighbours[nnbrs][3]=0;
				nnbrs++;
			}

			// Else calculate the displacement between the cells in terms of
			// number of cells. If this displacement is (-1,0,1), add the 
			// cell to the cell's neighbour list. 
			else {
				
				// Calculate distance between cells in terms of cells.
				ci = int(floor(cc/(Ly*Lz))) - int(floor(c/(Ly*Lz)));
				cj = ((int)floor(cc/Lz)%Ly) - ((int)floor(c/Lz)%Ly);
				ck = cc%Lz - c%Lz;
				
				// If the cell displacement is the entire cell length of
				// the system, set the cell displacement to (-1,1)
				// accordingly. This takes into account the periodic 
				// boundary conditions. 
				if(ci == (Lx-1)) ci =-1; 
				if(ci == -1*(Lx-1)) ci=1;
				if(cj == (Ly-1)) cj=-1;
				if(cj == -1*(Ly-1)) cj=1;
				if(ck == (Lz-1)) ck=-1;
				if(ck == -1*(Lz-1)) ck=1;
				
				// Check if the cell is a neighbour (has displacements in all directions
				// of (-1,0,1)). If so, add to the cell's neighbour list.
				if (ci>=-1 and ci<=1 and cj>=-1 and cj<=1 and ck>=-1 and ck<=1) { 

					if (!slit) {
						cell_list[c].neighbours[nnbrs][0]=cc;
						cell_list[c].neighbours[nnbrs][1]=ci;
						cell_list[c].neighbours[nnbrs][2]=cj;
						cell_list[c].neighbours[nnbrs][3]=ck;
						nnbrs++;
					}

					// If using the slit geometry, remove periodic boundary conditions
					// in the z-direction.
					else {
						if (cell_list[c].k == 0 and ck == -1) continue;
						else if (cell_list[c].k==(Lz-1) and ck == 1) continue;
						else {
							cell_list[c].neighbours[nnbrs][0]=cc;
							cell_list[c].neighbours[nnbrs][1]=ci;
							cell_list[c].neighbours[nnbrs][2]=cj;
							cell_list[c].neighbours[nnbrs][3]=ck;
							nnbrs++;
						}
					}
				}
			}
			// When all neighbours have been found, break the loop and move onto the
			// next cell.
			if (nnbrs == cell_list[c].num_neighbours) break;
		}
	}
	
	// Write the cell list to file for post-processing and debugging purposes.
	out.open(output_file_path/"cell_list", std::ios::out);

	for(int c=0;c<num_cells;c++) {

		out<<"Cell "<<c<<" "<<cell_list[c].i<<" "<<cell_list[c].j<<" "<<cell_list[c].k<<std::endl;
		out<<"Number of Neighbouring Cells: "<<cell_list[c].num_neighbours<<std::endl;
		out<<"Inclusion in multicanonical sampling (if applicable): "<<cell_list[c].sv<<std::endl;

		for (int cc=0;cc<cell_list[c].num_neighbours;cc++) {
			out<<cell_list[c].neighbours[cc][0]<<" "<<cell_list[c].neighbours[cc][1]<<" ";
			out<<cell_list[c].neighbours[cc][2]<<" "<<cell_list[c].neighbours[cc][3]<<" ";
			out<<cell_list[cell_list[c].neighbours[cc][0]].i<<" "<<cell_list[cell_list[c].neighbours[cc][0]].j<<" ";
			out<<cell_list[cell_list[c].neighbours[cc][0]].k<<std::endl;
		}
	}

	out.close();

	// If performing the cell list test, the program can now exit.
	if ((test) && (test_type =="cell_list")) exit(0);

	// Determine the seed for the random number generator. 
	// If performing the seed test, load the seed from the Input file 
	// and check that it is sensible. 
	if (test_type == "seed") {
		if (seed <=0 || seed>=2147483647) {
			std::cerr<<"Seed specified is not valid."<<std::endl;
			std::cerr<<"Seed must be integer between 0 and 2147483647."<<std::endl;
			if (load) std::cerr<<"Please check seed in Output file."<<std::endl;
			else std::cerr<<"Please check test seed in Input file."<<std::endl;
			exit(0);
		}
	}
	
	// If no seed is supplied, use the current system time. Save the seed to
	// file immediately for debugging purposes.
	else  {
		time_t *seed_time = new time_t;
		time(seed_time);
		seed = (long)*seed_time;
		out.open(output_file_path/"seed", std::ios::out);
		out<<seed<<std::endl;
		out.close();
	}

	/*
	---------------------------------------------------------
	Insert lines to intitialise random number generator here.
	---------------------------------------------------------
	*/

	// If using the solute geometry, add the solute to the simulation
	// box, checking that the radius supplied is sensible (e.g. larger
	// than 0 and smaller than the system size).
	if (bsolute) {
		if (Rs<=0.0) {
			std::cerr<<"Solute must have radius greater than 0."<<std::endl;
			if (load) std::cerr<<"Please check Output file."<<std::endl;
			else std::cerr<<"Please check Input file."<<std::endl;
			exit(0);
		}
		else if (Rs >= Lx*cut_off || Rs >= Ly*cut_off || Rs >= Lz*cut_off) {
			std::cerr<<"Solute must have radius smaller than box size."<<std::endl;
			if (load) std::cerr<<"Please check Output file."<<std::endl;
			else std::cerr<<"Please check Input file."<<std::endl;
			exit(0);
		}
		else {
			vsolute = solute(Rs, Lx, Ly, Lz);
			accessible_volume = volume - (4./3.)*pi*Rs*Rs*Rs; // Calculate volume accessible to the fluid.
		}
	}

	// If using a geometry other than solute, the entire volume of the system is
	// accessible to the fluid.
	else accessible_volume = volume; 

	// Check the sampling parameters supplied in the Input and Output files
	// are sensible.
	if (sweeps<=0) {
		std::cerr<<"SWEEPS must be greater than 0."<<std::endl;
		std::cerr<<"Please check Input file."<<std::endl;
		exit(0);
	} 
	else if (save_frequency <= 0) {
		std::cerr<<"SAVE_FREQUENCY must be greater than 0."<<std::endl;
		std::cerr<<"Please check Input file."<<std::endl;
		exit(0);
	}
	else if (distribution_frequency<=0) {
		std::cerr<<"DISTRIBUTION_FREQUENCY must be greater than 0."<<std::endl;
		std::cerr<<"Please check Input file."<<std::endl;
		exit(0);
	}

	// Set the number of sweeps to perform during equilibration. Note that this does
	// not mean that equilibration will occur. Equilibration occurs only if 
	// specified in the Input file. 
	equilibrium_sweeps = 0.2*sweeps;

	// If measuring the local compressibility and/or local thermal susceptibility,
	// set up the histogram and necessary constants.
	if (compressibility) {
		dmu = delta; sum_compressibility = 0.0;
		for (int i=0; i<profile_bins; i++) compressibility_profile.push_back(0.0);
	}
	if (susceptibility) {
		dT = (1.-temperature/(temperature+delta)); sum_susceptibility = 0.0;
		for (int i=0; i<profile_bins; i++) susceptibility_profile.push_back(0.0);
	}

	// Calculate the maximum distance between two particles within the system.
	// This is necessary for splitting the simulation box into distance bins
	// for measuring spatial distributions within the system. 
	if (slit) maxr = Lz;
	else  maxr = sqrt(Lx*Lx + Ly*Ly + Lz*Lz)/2.;

	// Calculate the width of each distance bin for spatial distributions
	// and then determine the boundary of each bin.  
	binsize = cut_off*sigma*maxr/profile_bins;
	for (int i=0; i<profile_bins;i++) {
		average_distribution.push_back({i*binsize,0.0});
	}

	// If using the transition matrix method to calculate the weights for
	// multicanonical sampling, set up the necessary matrices.
	if (tmmc) {
		for (int i=0;i<max_particles;i++) {
			collection.push_back({0.0,0.0,0.0}); transition.push_back({0.0,0.0,0.0});
			collection_totals.push_back(0.0); weights.push_back(0.0);
		}
	}

	// If a maximum and minimum density have been supplied, calculate the corresponding
	// maximum and minimum number of particles the system can have. If a maximum and
	// minimum density have not been supplied, set the minimum number of particles to be
	// 0, and the maximum number to be equivalent to a density of 1.01gcm^-3 for monatomic
	// water and 0.75 for a truncated Lennard-Jones fluid. Check that the maximum 
	// number of particles is greater than the minimum number of particles, and if not
	// set the maximum and minimum numbers to the default.
	if (mW) {
		maximum_particles = 0.1*maximum_density*(avagadro/mW_MW)*pow(mW_sigma,3)*accessible_volume;
		minimum_particles = 0.1*minimum_density*(avagadro/mW_MW)*pow(mW_sigma,3)*accessible_volume;
		if (initial_density > minimum_density) num_particles = 0.1*initial_density*(avagadro/mW_MW)*pow(mW_sigma,3)*accessible_volume;
		else num_particles = minimum_particles;
		std::cout<<"Maximum Particles: "<<maximum_particles<<std::endl;
		std::cout<<"Minimum particles: "<<minimum_particles<<std::endl;
		if (minimum_particles >= maximum_particles) {
			std::cerr<<"Error: Maximum density must be larger than minimum density."<<std::endl;
			std::cerr<<"Setting maximum and minimum particles to default parameters."<<std::endl;
			maximum_particles = 0.1*1.01*(avagadro/mW_MW)*pow(mW_sigma,3)*accessible_volume;
			minimum_particles = 0.0;	
		}
	}
	else {
		maximum_particles = maximum_density*accessible_volume; 
		minimum_particles = minimum_density*accessible_volume; 
		if (initial_density > minimum_density) num_particles = initial_density*accessible_volume;
		else num_particles = minimum_particles;
		if (minimum_particles >= maximum_particles) {
			std::cerr<<"Error: Maximum density must be larger than minimum density."<<std::endl;
			std::cerr<<"Setting maximum and minimum particles to default parameters."<<std::endl;
			maximum_particles = 0.75*accessible_volume; 
			minimum_particles = 0.0;	
		}
	}

	// If the calculated maximum number of particles in the system exceeds the maximum allowed
	// number of particles, alert the user and exit the program. This implies that the user
	// has selected a system size which is too big and may run into memory issues.
	if (maximum_particles > max_particles) {
		std::cerr<<"Error: Maximum density of system exceeds available memory in particle list."<<std::endl;
		std::cerr<<"Please reduce the maximum density, or move to a smaller system size."<<std::endl;
		exit(0);
	}

	// Setup the particle list by filling it with particles. 
	for (int i=0;i<maximum_particles;i++) particle_list.push_back(particle(9999.9,9999.9,9999.9,maximum_particles));

	// Now place the initial particles in the system.
	place_particles(load, positions);

	// If loading a previous simulation, load the previous distributions to continue
	// adding to them.
	if (load) { 
		load_previous_distributions();
		equilibrate = false;
	}
	
	// If using multicanonical sampling, load the intiial weights, if applicable.
	if (load_weights) {
		if (load || positions) {
			load_previous_weights(true);
			if (tmmc) load_previous_collection();
		}
		else load_previous_weights(false);
	}
}

void sim_box::place_particles(bool load, bool positions){

	/* 
	Places the initial particles in the system. If loading a previous simulation,
	or an initial configuration of particles has been supplied, then the particle
	positions are loaded from file. Otherwise, particles are placed within the 
	system at random.
	*/
	
	int m=0, pcell, c;
	double x,y,z;
	particle p;
	bool valid = true;

	// If an initial configuration is supplied, load the configuration.
	if (load||positions) load_particles();

	// Else place particles in the system at random based on the initial density
	// or minimum density supplied.
	else {
		while (m<num_particles) {
			
			// Generate a cell and position within cell for the particle at random.
			c = floor(get_rnd()*num_cells);
			x = get_rnd(); y = get_rnd(); z = get_rnd();
			p = particle(x,y,z, c);

			// If using the solute geometry, check particle is outside of the solute.
			if (bsolute) valid = check_valid_position(p);

			// If the particle is outside the solute, or if using the slit of bulk
			// geometry, add the particle to the particle list and update the 
			// relevant cell's occupancy list and number of occupants.
			if (valid) {
				particle_list[m] = p;
				cell_list[c].occupancy[cell_list[c].num_occupants]=m; 
				cell_list[c].num_occupants+=1;	
				m++;

				// If using sub volume sampling, and particle is within sub volume,
				// increase sub_volume particle count.
				if (sub_volume) {
					if (cell_list[c].sv) sub_volume_particles +=1;
				}
			}
		}
	}
}


void sim_box::load_particles() {
	
	/*
	Reads in the positions of particles from file and places these within
	the simulation box. Assumes the name of the file is 'positions_final'
	and that it is located at the specified 'load file path'.
	*/

	int id, cell_id, m=0, c;
	double x,y,z;
	std::ifstream positions; 

	positions.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	// Try to open the particle positions file. If the file cannot be found
	// alert the user and exit the program.
	try {
		positions.open(load_file_path/"positions_final");
	}
	catch (std::ifstream::failure e) {
   		std::cerr <<"Could not locate positions file at: "<<load_file_path/"positions_final"<<std::endl;
		std::cerr<<"Please check file exists in this location."<<std::endl;
		exit(0);
  	}
	
	// Read in the position of each particle from the file and add the particle
	// to the relevant lists.
	try {
		while(positions>>id>>x>>y>>z>>c) {
			particle_list[m] = particle(x,y,z,c);
			cell_list[c].occupancy[cell_list[c].num_occupants] = m;
			cell_list[c].num_occupants+=1;
			m++;
		}
	}
	catch (std::ifstream::failure e) {
		positions.close();
		std::cout<<"Completed reading "<<load_file_path/"positions_final"<<std::endl;
	}
	num_particles = m;
}

void sim_box::reorder_particles(int pn){

	/*
	The particle list is organised such that the number of particles in the system
	indicates the last position in the particle list that a particle exists. When a 
	particle is deleted from the system, the bottom particle in the particle list
	must therefore be moved to the deleted particle's position in the list. This
	function does this, and then updates the relevant cell occupancy lists to
	reflect the new position of the particle.
	*/

	int cell_id, to_shift, i, pcell = particle_list[pn].cell_id;
	particle p = particle_list[pn];

	// If there is only one particle in the list, remove the particle from
	// its cell.
	if(num_particles==1) {
		cell_list[pcell].occupancy[0] = max_particles;
		cell_list[pcell].num_occupants-=1;
	}

	else {
		cell_id = particle_list[num_particles-1].cell_id; //Find the cell of the last particle in the particle list.
			
		// Update the cell of the last particle to reflect its new position in the
		// particle list.
		for(i=0;i<cell_list[cell_id].num_occupants;i++) {
			if(cell_list[cell_id].occupancy[i]==num_particles-1) {
				cell_list[cell_id].occupancy[i]=pn; 
				break;
			}
		}
		particle_list[pn] = particle_list[num_particles-1]; //Move the last particle to the empty position.

		for (i=0;i<cell_list[pcell].num_occupants;i++) {
			// Find the position of the particle to delete in its cell occupancy list.
			// Store this position, as we'll move the last particle in the cell's 
			// occupancy list to this position. 
			if(cell_list[pcell].occupancy[i]==pn) {
				to_shift  = i; 
				break;
			}	
		}
		
		// Move the last particle in the deleted particle's cell occupancy list to the 
		// vacated position. 
		cell_list[pcell].occupancy[to_shift] = cell_list[pcell].occupancy[cell_list[pcell].num_occupants-1];
		
		// Remove the last particle in the deleted particles's cell occupancy list.
		cell_list[pcell].occupancy[cell_list[pcell].num_occupants-1] = max_particles;

		// Decrease the occupancy of the cell.
		cell_list[pcell].num_occupants-=1;
	}
}

void sim_box::update_collection(int N, int Na, double prob) {

	/*
	Updates the collection matrix if using the transition matrix method
	of finding weights for multicanonical sampling. See
	G.R. Smith and A.D. Bruce, J. Phys. A: Math Gen. 28 6623-43 (1995)
	or
	D.J. Ashton and N.B. Wilding, Molecular Physics 109 999-1007 (2010) 
	for details on transition matrix method.
	*/

	// N represents the current number of particles. Na represents
	// the proposed number of particles, e.g. Na = 0 -> N-1 particles,
	//  Na = 1 -> N particles and Na = 2 -> N + 1 particles.
	collection[N][Na]+=prob;  collection_totals[N]+=prob;
}

void sim_box::update_transition(void) {

	/*
	Updates the transition matrix if using the transition matrix method
	of finding weights for multicanonical sampling. See
	G.R. Smith and A.D. Bruce, J. Phys. A: Math Gen. 28 6623-43 (1995)
	or
	D.J. Ashton and N.B. Wilding, Molecular Physics 109 999-1007 (2010) 
	for details on transition matrix method.
	*/

	// i = 0 -> j-1 particles in system.
	// i = 1 -> j particles in system.
	// i = 2 -> j+1 particles in system.
	for(int j = 0;j<max_particles; j++) {
		for(int i=0; i<3; i++) {
			if (collection_totals[j]>0.0) transition[j][i] = (collection[j][i])/(collection_totals[j]);
		}
	}
}

void sim_box::update_weights(long sweeps) {

	/*
	Updates weights for biasing if using transition matrix method for finding
	weights for multicanonical biasing. See
	G.R. Smith and A.D. Bruce, J. Phys. A: Math Gen. 28 6623-43 (1995)
	or
	D.J. Ashton and N.B. Wilding, Molecular Physics 109 999-1007 (2010) 
	for details on transition matrix method.
	*/

	double new_weights[max_particles];
	bool valid_weights = true;

	// All weights are calculated relative to one another. We set the weight to
	// compare to to 1.0.
	new_weights[0] = 1.0;
	
	// If possible, calculate the new weights according to the transition matrix method.
	// Else, set the new weight to the value of the previous weight.
	for(int i=0;i<max_particles-1;i++) {
		if ((transition[i][2] > 0.0) && (transition[i+1][0] > 0.0)) new_weights[i+1] = new_weights[i]*(transition[i][2]/transition[i+1][0]);
		else new_weights[i+1] = new_weights[i];
	}

	// Try to calculate the log of the weight. If nan or inf is returned,
	// the new weights do not represent a valid set. In this case, use the
	// current weights in the system. If the new weights do represent a
	// valid set, set the weights of the system to the new weights.
	for (int i=0;i<max_particles; i++) {
		new_weights[i] = log(new_weights[i]);
		if (isnan(new_weights[i]) || isinf(new_weights[i])) valid_weights = false;
		if (valid_weights) for(int i=0; i<max_particles; i++) weights[i] = new_weights[i];
	}
}

void sim_box::load_previous_collection(void){

	/*
	Reads in the collection matrix found during a previous simulation from file
	and sets this collection matrix of the current simulation to equal this.
	*/

	double ca,cb,cc; int m=0;
	std::ifstream collect;

	collect.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	// Try to open the file. If the file is not found, alert the user and 
	// exit the program.
	try {
		collect.open(load_file_path/"collection_final");
	}
	catch (std::ifstream::failure e) {
   		std::cerr <<"Could not locate collection file at: "<<load_file_path/"collection_final"<<std::endl;
		std::cerr<<"Please check file exists in this location."<<std::endl;
		exit(0);
  	}
	
	// Read contents of file and set the collection matrix of the current
	// system to this.
	try {
		while(collect>>ca>>cb>>cc) {
			collection[m][0] = ca;
			collection[m][1] = cb;
			collection[m][2] = cc;
			collection_totals[m] = ca+cb+cc;
			m++;
		}
	}
	catch (std::ifstream::failure e) {
		collect.close();
		std::cout<<"Completed reading "<<load_file_path/"collection_final"<<std::endl;
	}
}

void sim_box::load_previous_weights(bool load){

	/*
	Reads in weights to use within the simulation from file. If load is true,
	these weights are read in from a previous simulation. The weights file
	in this case is expected to exist in the 'load file path' directory and
	to be names 'weights_final'. If load is false, these weights are read 
	from a file called 'weights_initial' which is assumed to exist in the 
	'output file path' directory supplied to the system. 
	*/

	std::ifstream fweight; 
	double w; int m=0;
	std::filesystem::path path;
	
	// Determine the location and file name of the weights file to be read.
	if (load) {
		path = load_file_path/"weights_final";
	}
	else {
		path = output_file_path/"weights_initial";
	}

	fweight.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	// Attempt to open file. If the file does not exist in the expected location,
	// alert the user and exit the program. 
	try {
		fweight.open(path);
	}
	catch (std::ifstream::failure e) {
   		std::cerr <<"Could not locate weights file at: "<<path<<std::endl;
		std::cerr<<"Please check file exists in this location."<<std::endl;
		exit(0);
  	}
	
	// Read contents of file, and set weights accordingly.
	try {
		while(fweight>>w) { weights[m] = w; m++;}
	}
	catch (std::ifstream::failure e) {
		fweight.close();
		std::cout<<"Completed reading "<<path<<std::endl;
	}
}

void sim_box::read_distribution_file(std::filesystem::path path, bool density, bool compressibility) {

	/*
	Reads in distributions collected during a previous simulation and updates the 
	distributions within the current simulation to reflect these.
	*/

	std::string line, arg1, arg2; std::ifstream fdist; 
	int m=0; bool space=false;

	fdist.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	// Try to open previous spatial distribution file. If the file could not
	// be found, alert the user and continue with the simulation, starting 
	// collecting distribution results from scratch. 
	try {
		fdist.open(path);
	}
	catch (std::ifstream::failure e) {
   		std::cerr <<"Could not locate average distribution file at: "<<path<<std::endl;
		std::cerr<<"The simulation will continue without loaded distribution."<<std::endl;
  	}
	
	// Read contents of file and set distributions accordingly.
	try {
		while(getline(fdist,line)) { 

			// Values in file are separated by white space.
			// First check for lines containing histogram values.
			if (line[0] != 'S' && line.length() != 0) {
				for (auto i:line) {
					if (i == ' ') space = true;
					else if (!space)  arg1 = arg1 +i;
					else arg2 = arg2+ i; 
				}

				// Set the relevant array element.
				if (density) {
					average_distribution[m][0] = std::stod(arg1); average_distribution[m][1] = std::stod(arg2);
				}
				else if (compressibility) compressibility_profile[m] = std::stod(arg2);
				else susceptibility_profile[m] = std::stod(arg2);

				// Reset variables.
			 	m++; arg1.erase(); arg2.erase(); space = false;
			}

			// Check for lines detailing the sums of the number of distributions.
			// This is used when averaging. 
			else if (line.length() !=0) {

				// Values are separated by white space.
				for (auto i:line) {
					if (i == ' ') {
						if (arg1.length() > 3) space = true;
					}
					else if (!space)  {arg1 = arg1 + i;}
					else arg2 = arg2 + i; 
				}

				// Set the relevant sum. 
				if (density) {sum_distributions = std::stod(arg2); std::cout<<"Distribution sum: "<<sum_distributions<<std::endl;}
				else if (compressibility) {sum_compressibility = std::stod(arg2); std::cout<<"Compressibility sum: "<<sum_compressibility<<std::endl;}
				else {sum_susceptibility = std::stod(arg2); std::cout<<"Distribution susceptibility: "<<sum_susceptibility<<std::endl;}
			}
		}
	}

	catch (std::ifstream::failure e) {

		fdist.close();

		// If the sum of distributions was not found, alert the user and exit the program.
		if (density && sum_distributions <= 0.0) {std::cout<<"Error, did not load sum of distributions."<<std::endl; exit(0);}
		if (compressibility && sum_compressibility <= 0.0) {std::cout<<"Error, did not load sum of compressibility."<<std::endl; exit(0);}
		if (!density && !compressibility && sum_susceptibility <= 0.0) {std::cout<<"Error, did not load sum of susceptibility."<<std::endl; exit(0);}
		std::cout<<"Completed reading "<<path<<std::endl;
	}
}
void sim_box::load_previous_distributions(void){

	/*
	Reads in relevant distributions from a previous simulation from file.
	*/

	std::filesystem::path path;

	// Particle spatial distribution.
	path = load_file_path/"averaged_distributions"/"average_distribution_final_raw";
	read_distribution_file(path, true, false);

	// Local compressibility distribution.
	if (compressibility) {
		path = load_file_path/"averaged_distributions"/"local_compressibility_final_raw";
		read_distribution_file(path, false, true);
	}

	// Local Thermal Susceptibility distribution.
	if (susceptibility) {
		path = load_file_path/"averaged_distributions"/"local_thermal_susceptibility_final_raw";
		read_distribution_file(path, false, false);
	}
}

void sim_box::output_weights(std::string name) {
	
	/*
	Writes the current collection, transition and weight matrices to file.
	*/

	int i, j;
	std::ofstream hist; std::filesystem::path out_name;

	// Write to file the current collection matrix.
	if (name == "final") out_name = output_file_path/("collection_" + name);
	else out_name = output_file_path/"collection_matrices"/("collection_" + name);
	hist.open(out_name, std::ios::out);
	hist<<std::fixed<<std::setprecision(10);
	for(i = 0;i<max_particles; i++) {
		for(j=0; j<3;j++) {
			hist<<collection[i][j]<<" ";
		}
		hist<<std::endl;
	}
	hist<<std::endl;
	hist.close();

	// Write to file the current transition matrix.
	if (name == "final") out_name = output_file_path/("transition_" + name);
	else out_name = output_file_path/"transitions_matrices"/("transition_" + name);
	hist.open(out_name, std::ios::out);
	hist<<std::fixed<<std::setprecision(10);

	for(i = 0;i<max_particles; i++) {
		for(j=0; j<3;j++) {
			hist<<transition[i][j]<<" ";
		}
		hist<<std::endl;
	}
	hist<<std::endl;
	hist.close();

	// Write to file the current weights.
	if (name == "final") out_name = output_file_path/("weights_" + name);
	else out_name = output_file_path/"weight_distributions"/("weights_" + name);
	hist.open(out_name, std::ios::out);
	hist<<std::fixed<<std::setprecision(10);
	for(i = 0;i<max_particles; i++) hist<<weights[i]<<" "<<std::endl;
	hist.close();
}

void sim_box::output_cells(std::string name) {

	/*
	Writes the cell list, along with each cells current occupying particles,
	to file. Can be used for debugging.
	*/

	std::ofstream cell_info;
	particle part;
	std::filesystem::path path;
	
	path = output_file_path/("cell_information_" + name);

	cell_info.open(path, std::ios::out);
	for (int c=0; c<num_cells; c++) {
		cell_info << "Cell " <<c<<" i j k "<< cell_list[c].i<<" "<<cell_list[c].j<<" "<<cell_list[c].k;
		cell_info<<" occupancy "<<cell_list[c].num_occupants<<" sample: "<<cell_list[c].sv<<std::endl;
		for (int p=0; p<cell_list[c].num_occupants; p++) {
			part = particle_list[cell_list[c].occupancy[p]];
			cell_info<<p<<" "<<cell_list[c].occupancy[p]<<" "<<part.cell_id<<" "<<part.x;
			cell_info<<" "<<part.y<<" "<<part.z<<std::endl;
		}
	}
	cell_info.close();
}

bool sim_box::check_valid_position(particle temp) {

	/*
	Check's that a particle exists outside of the solute.
	*/

	double dx, dy, dz, r;
	bool valid;
	cell pcell = cell_list[temp.cell_id];

	// Calculate radial distance of particle from centre of solute.
	dx = vsolute.x - (pcell.i + temp.x);
	dy = vsolute.y - (pcell.j + temp.y);
	dz = vsolute.z - (pcell.k + temp.z);
	r = dx*dx + dy*dy + dz* dz; r= sqrt(r);
	if (mW) r *= mW_a;
	else r *= LJ_rc;

	// Compare this radial distance to the radius of the solute.
	if (r> vsolute.Rs) valid = true;
	else valid = false;

	return valid;
}

// Functions to return system parameters.
double sim_box::return_volume(void) {return volume;}
bool sim_box::return_load(void) {return load;}
bool sim_box::return_test(void) {return test;}
std::string sim_box::return_test_type(void) {return test_type;}
int sim_box::return_Lx(void) {return Lx;}
int sim_box::return_Ly(void) {return Ly;}
int sim_box::return_Lz(void) {return Lz;}
int sim_box::return_num_cells(void) {return num_cells;}
long sim_box::return_seed(void) {return seed;}
long sim_box::return_total_sweeps(void) {return total_sweeps;}
long sim_box::return_sweeps_eqb(void) {return equilibrium_sweeps;}
long sim_box::return_sweeps(void) {return sweeps;}
int sim_box::return_save_frequency(void) {return save_frequency;}
int sim_box::return_distribution_frequency(void) {return distribution_frequency;}
bool sim_box::return_output_distributions(void) {return output_individual_distributions;}
int sim_box::return_profile_bins(void) {return profile_bins;}
bool sim_box::return_compressibility(void) {return compressibility;}
bool sim_box::return_susceptibility(void) {return susceptibility;}
double sim_box::return_dmu(void) {return dmu;}
double sim_box::return_dT(void) {return dT;}
std::filesystem::path sim_box::return_output_file_path(void) {return output_file_path;}
std::filesystem::path sim_box::return_load_file_path(void) {return load_file_path;}

// Functions used for random number generation. 
void sim_box::check_rnd(void) {

	/*
	This function can be used if the random number generator outputs
	a number of random numbers at a time. After a random number is used,
	the program should check that there are still numbers left. If there
	are not, the program should generate more random numbers from the
	generator.
	*/

	if (nrands-rand_check<2) {
		/*
		-------------------------------------------------------
		Insert lines to generate new random numbers array here.
		-------------------------------------------------------
		*/
		rand_check = 0; 
	}
}
double sim_box::get_rnd(void) {

	/*
	Returns a random number.
	*/

	double rnd;
	/* 
	-----------------------------------------------------------
	Insert lines to call random number generator here.
	Set random number to variable rnd.
	If random number generator returns arrays of random number,
	insert call to check_rnd function, to ensure the program
	does not run out of random numbers.
	-----------------------------------------------------------
	*/

	return rnd;
}
			
// Functions which provide details of the number of accepted
// and rejected moves.
void sim_box::inc_mca(void) {mca=mca+1;} 
void sim_box::inc_mcr(void) {mcr=mcr+1;}
long sim_box::return_mca(void) {return mca;}
long sim_box::return_mcr(void) {return mcr;}

void sim_box::update_avg_acc(void) {

	/*
	Updates the average move acceptance total and resets the 
	accepted and rejected move counters.
	*/
	avg_acc += (double)mca/(double)mcr;
	acc_entries++;
	mca = 0; mcr = 0;
}

double sim_box::return_avg_acc(void) {return avg_acc;}
long sim_box::return_acceptance_entries(void) {return acc_entries;}
			
// Functions which provide details of the running energy of the system.
void sim_box::inc_running_energy(double energy) {running_energy=running_energy+energy;}
void sim_box::dec_running_energy(double energy) {running_energy=running_energy-energy;}
double sim_box::return_running_energy(void) {return running_energy;}
void sim_box::set_running_energy(double energy) {running_energy = energy;}
			
// Functions which provide details of or alter the number
// of particles within the system.
void sim_box::inc_nparticles(void) {num_particles+=+1;}
void sim_box::dec_nparticles(void) {num_particles-=1;}
int sim_box::return_nparticles(void) {return num_particles;}
void sim_box::set_nparticles(int particles) {num_particles=particles;}
double sim_box::return_maximum_particles(void) {return maximum_particles;}
double sim_box::return_minimum_particles(void) {return minimum_particles;}
double sim_box::return_minimum_density(void) {return minimum_density;}
double sim_box::return_maximum_density(void) {return maximum_density;}

// Functions related to system parameters.
double sim_box::return_mu(void) {return mu;}
double sim_box::return_epsilon(void) {return epsilon;}
double sim_box::return_temperature(void) {return temperature;}
bool sim_box::return_mW(void) {return mW;}
bool sim_box::return_mc(void) {return mc;}
bool sim_box::return_tmmc(void) {return tmmc;}
bool sim_box::return_equilibrate(void) {return equilibrate;}
bool sim_box::return_sub_volume(void) {return sub_volume;}
int sim_box::return_sub_volume_dz(void) {return sub_volume_dz;}

// Functions related to sub-volume sampling.
int sim_box::return_sv_particles(void) {return sv_particles;}
void sim_box::inc_sv_particles(void) {sv_particles+=1;}
void sim_box::dec_sv_particles(void) {sv_particles-=1;}

// Functions related to information about weights within the system.
bool sim_box::return_load_weights(void) {return load_weights;}
double sim_box::return_weight(int N) {return weights[N];}
double sim_box::return_prob(int Na, int Nb) {return transition[Na][Nb];}

// Functions related to information about any external potential applied to
// the system.
bool sim_box::return_LJ_Vext(void) {return vext;}
double sim_box::return_epsilon_wall(void) {return surface_interaction;}
void sim_box::inc_running_Vext(double Vext) {running_Vext=running_Vext+Vext;}
void sim_box::dec_running_Vext(double Vext) {running_Vext=running_Vext-Vext;}
double sim_box::return_running_Vext(void) {return running_Vext;}
void sim_box::set_running_Vext(double Vext) {running_Vext = Vext;}

// Functions which provide information about the geometry of the system. 
bool sim_box::return_slit(void) {return slit;}
bool sim_box::return_solute(void) {return bsolute;}
