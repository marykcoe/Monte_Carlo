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

This file contains the particle class.
*/

#ifndef PARTICLE_H
#define PARTICLE_H

class particle {
	public: 		
			double x, y, z; // position of molecule within cell
			int cell_id; 		//index of occupied cell in cell_list
	
			particle() {};
			particle(double xin, double yin, double zin, int icell_in);
};

#endif