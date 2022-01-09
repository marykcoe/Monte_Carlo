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

This file contains the contructor of the solute class. 
*/


#include "solute.h"

solute::solute(double iRs, int Lx, int Ly, int Lz) {

    /* Sets the position and radius of the solute.*/

    x = Lx/2.; y = Ly/2.; z = Lz/2.;
    Rs = iRs;
}