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

This file defines the cell class, which is used to form the 
cell list. 
*/

#ifndef CELL_H
#define CELL_H

#include "constants.h"

using namespace constants;

class cell {
	
	/*Class structure for a cell*/
	
	public: 
			int occupancy[20]; //ids of occupants of the cell
			int num_occupants; //number of occupants in cell
			int i,j,k; //position of cell in real space w.r.t other cells (imagine 3D grid)
			int neighbours[27][4]; //Gives information of neighbouring cells within cell list
			int num_neighbours; // Number of neighbouring cells
			int sv; // Indicates if cell is within sub-volume for TMMC sampling
			
			cell() {};
			cell(int iin, int jin, int kin, int innbrsm, int isv); //constructor for cell
};

#endif