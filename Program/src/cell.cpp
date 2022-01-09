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

#include "cell.h"

cell::cell(int iin, int jin, int kin, int innbrs, int isv) {
    
    /*
    Defines a cell within the cell list.
    */

    // Position within simulation box. This is measured from
    // the bottom front left corner of the box.
    i = iin; j=jin; k=kin; 

    // Number of particles within the cell.
    num_occupants = 0; 

    // Number of neighbours cell has (depends on periodic boundary
    // conditions).
    num_neighbours = innbrs;

    // The occupancy list contains the position of particles in the 
    // particle list which are within the cell.
    for(int m=0;m<20;m++) occupancy[m] = max_particles;

    // List of cells which neighbour this one.
    for(int m=0; m<27; m++) {
        neighbours[m][0] = 999999; neighbours[m][1] = 999999; 
        neighbours[m][2] = 999999; neighbours[m][3] = 999999; 
    }

    // Determines whether the cell is included within the sub-volume
    // sampling.
    sv = isv;
}