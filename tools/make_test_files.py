#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python tools to accompany Monte Carlo program to simulate
a truncated Lennard-Jones or monatomic water liquid either
in bulk, in contact with a small solute, or confined to a slit.

For information on how this program works, please consult the
accompanying tutorials or Section 4.2 of the following thesis:

M. K. Coe, Hydrophobicity Across Length Scales: The Role of
Surface Criticality, Ph.D. Thesis, University of Bristol (2021)
Available at: https://research-information.bris.ac.uk/ws/portalfiles/portal/304220732/Thesis_Mary_Coe.pdf

Copyright (c) 2022 Mary Coe
Made available under the MIT License.

This file contains tools for making test files for the program.
"""

import numpy as np
import sys
import os

def make_cell_list(Lx, Ly, Lz, slit, lj_vext = True, sv = False,
                   dz = None):

    """
    Generates a cell list in exact style of the Monte Carlo program.
    If lj_vext is False, the list is written to file.

    Parameters
    ----------
    Lx : int
        Length of simulation box in x-direction in units of cell.
    Ly : int
        Length of simulation box in y-direction in units of cell.
    Lz : int
        Length of simulation box in z-direction in units of cell.
    slit : bool
        Determines whether the geometry of the system is a slit
        (slit = True) or a solute/bulk (slit = False)
    lj_vext : bool, optional
        Determines whether the test is for the external potential
        (lj_vext = True) or whether the test is for the cell list
        (lj_vext = False). The deafult is True.
    sv : bool, optional
        Determines whether sub-volume sampling is enabled. The
        default is False.
    dz : None/int, optional
        If sub-volume sampling is enabled, dz must be supplied.
        This represents the number of cells in the z-direction
        which form the sub-volume to sample.

    Returns
    --------
    cell_list : array
        Cell list in where the columns in order represent the
        x position, y position, z position, number of neighbours and,
        if using sub-volume sampling, whether the cell is included
        within the sampling.
    """

    # Calculate the total number of cells and make an
    # array of the correct size.
    num_cells = Lx*Ly*Lz
    cell_list = np.zeros((num_cells,5), dtype=int)
    cell_id = 0

    # If using sub-volume sampling, calculate the boundaries of
    # the cells which fall into sub-volume region.
    if sv:
        if dz is None:
            print('Error: If sub_volume sampling is enabled (sv=True) then')
            print('dz must be supplied.')
            sys.exit(0)

        sv_lower = 0; sv_upper = Lz-1; upper = False;
        while ((sv_upper-sv_lower) >= dz):
            if upper:
                sv_upper-=1; upper = False
            else:
                sv_lower+=1; upper = True

    # Make the cell list in the same way as the Monte Carlo
    # program.
    for i in range(Lx):
        for j in range(Ly):
            for k in range(Lz):
                if not slit:
                    nnbrs = 27;
                else:
                    if (k>0 and k<(Lz-1)):
                        nnbrs = 27;
                    else:
                        nnbrs = 18;

                cell_list[cell_id,0] = i
                cell_list[cell_id,1] = j
                cell_list[cell_id,2] = k
                cell_list[cell_id,3] = nnbrs

                if sv:
                    if (k>=sv_lower) and (k<=sv_upper):
                        cell_list[cell_id,4] = 1
                    else:
                        cell_list[cell_id,4] = 0
                else:
                    cell_list[cell_id,4] = 1

                cell_id+=1

    # If performing the cell_list test, find the neighbours of each
    # cell and output to file
    if not lj_vext:

        # Make an empty output file for the cell list. This is done
        # to clear any previous test files which may have been
        # written.
        with open(os.path.join('Tests','cell_list','test_cell_list'), 'w') as out:
            pass

        # Find the neighbouring cells within periodic boundary conditions
        # for each cell in the list.
        for c in range(num_cells):

            nnbrs=0; neighbours = np.zeros((int(cell_list[c,3]),4), dtype=int)
            for cc in range(num_cells):

                if cc==c:
                    neighbours[nnbrs,0]=cc
                    neighbours[nnbrs,1]=0
                    neighbours[nnbrs,2]=0
                    neighbours[nnbrs,3]=0
                    nnbrs+=1
                else:
                    # Calculate distance between cells in terms of cells.
                    ci = int(np.floor(cc/(Ly*Lz))) - int(np.floor(c/(Ly*Lz)))
                    cj = (int(np.floor(cc/Lz)%Ly)) - (int(np.floor(c/Lz)%Ly))
                    ck = cc%Lz - c%Lz

                    # Take into account the periodic boundary conditions
                    if ci == (Lx-1):
                        ci =-1
                    if ci == -1*(Lx-1):
                        ci=1
                    if cj == (Ly-1):
                        cj=-1
                    if cj == -1*(Ly-1):
                        cj=1;
                    if ck == (Lz-1):
                        ck=-1
                    if ck == -1*(Lz-1):
                        ck=1;

                    # Finally, if the neighbour is found, add it to the neighbour list.
                    if (ci>=-1) and (ci<=1) and (cj>=-1) and (cj<=1) and (ck>=-1) and (ck<=1):

                        # Slit geometries do not include periodic conditions in Z dimension.
                        if not slit:
                            neighbours[nnbrs,0]=cc;
                            neighbours[nnbrs,1]=ci
                            neighbours[nnbrs,2]=cj
                            neighbours[nnbrs,3]=ck
                            nnbrs+=1;
                        else:
                            if (cell_list[c,2] == 0) and (ck == -1):
                                continue
                            elif (cell_list[c,2]==(Lz-1)) and (ck == 1):
                                continue
                            else:
                                neighbours[nnbrs,0]=cc
                                neighbours[nnbrs,1]=ci
                                neighbours[nnbrs,2]=cj
                                neighbours[nnbrs,3]=ck
                                nnbrs+=1;

                if (nnbrs == cell_list[c,3]):
                    break

            # Once the cell list has been generated, write it to file.
            with open(os.path.join('Tests','cell_list','test_cell_list'), 'a') as out:
                out.write(f'Cell {c} {cell_list[c,0]} {cell_list[c,1]} {cell_list[c,2]}\n')
                out.write(f'Number of Neighbouring Cells: {cell_list[c,3]}\n')
                out.write(f'Inclusion in multicanonical sampling (if applicable): {cell_list[c,4]}\n')

                for n in range(nnbrs):
                    out.write(f'{neighbours[n,0]} {neighbours[n,1]} {neighbours[n,2]} {neighbours[n,3]} ')
                    out.write(f'{cell_list[int(neighbours[n,0]),0]} {cell_list[int(neighbours[n,0]),1]} ')
                    out.write(f'{cell_list[int(neighbours[n,0]),2]}\n')

    return cell_list

def lj_vext(z, slit = True, Lz = None, Rs = None):

    """
    Calculates the external potential felt by a particle located
    at position z. This external potential takes the form of a
    shifted Lennard-Jones potential with degrees of freedom in
    which the density if homogeneous integrated out.

    Parameters
    ----------
    z : float
        Position of the particle. If using the slit geometry, this is
        the z-position of the particle relative to the left wall.
        If using the solute geometry, this is the radial distance of
        the particle from the centre of the solute.
    slit : bool, optional
        Determines whether to use a slit geometry. The default is
        True.
    Lz : None/float, optional
        If using the slit geometry, the length of the box in the
        z-direction must be supplied. The default is None.
    Rs : None/float, optional
        If using the solute geometry, the radius of the solute must
        be supplied. The default is None.
    """

    # Check parameters supplied are sensible.
    if (slit) and (Lz is None):
        print('Error. If using planar geometry, must supply length of')
        print('simulation box (Lz).')
        sys.exit(0)
    elif not slit and (Rs is None):
        print('Error. If using solute geometry, must supply radius of')
        print('solute (Rs).')
        sys.exit(0)

    # Calculate the location of the minimum of the potential
    # in the case of a planar wall.
    zshift = np.power(0.4,1./6.)

    # Calculate the external potential.
    if slit:
        vext = ((2.0/15.0)*np.power(1./(z+zshift),9.0)-\
                np.power(1./(z+zshift),3.0))
        vext += ((2.0/15.0)*np.power(1./(Lz+zshift - z),9.0)-\
                 np.power(1./(Lz+zshift - z),3.0))
    else:
        rR_plus = 1./(z+zshift+Rs); rR_minus = 1./(z+zshift-Rs);
        vext = (2.0/15.0)*(np.power(rR_minus,9.0)-np.power(rR_plus,9.0))
        vext += (3.0/20.0)*(1./(z+zshift))*(np.power(rR_plus,8.0)-np.power(rR_minus,8.0))
        vext += (np.power(rR_plus,3.0)-np.power(rR_minus,3.0))
        vext += (3.0/2.0)*(1./(z+zshift))*(np.power(rR_minus,2.0)-np.power(rR_plus,2.0))

    return vext

def make_lj_vext_test_files(Lx, Ly, Lz,
                            particles = 30, slit = True, Rs = None,
                            mW = True):

    """
    Generates positions for particles in a Monte Carlo system
    with a given geometry. These are then written to file in
    the format required by the Monte Carlo program.

    Parameters
    ----------
    Lx : int
        Length of the simulation box in the x-direction in terms
        of cells.
    Ly : int
        Length of the simulation box in the y-direction in terms
        of cells.
    Lz : int
        Length of the simulation box in the z-direction in terms
        of cells.
    particles : int, optional
        The number of particles to generate. The default is 30.
    slit : bool, optional
        Determines whether to use a slit geometry. The default is
        True.
    Rs : None/float, optional
        Determines whether to use a solute geometry, and if so the
        radius of solute to use. If None, a solute geometry is not
        used. The default is None.
    mW : bool, optional
        Determines whether to use a monatomic water fluid. If False,
        a truncated Lennard-Jones fluid is used instead. The default
        is True.

    Returns
    -------
    positions : array(float)
        The positions of the randomly generated particles.
    """

    # If using the solute geometry, it will be necessary to
    # check that the particle generated does not fall within
    # the solute itself. To do that, it is necessary to
    # know the cell list and cut-off radius of interaction.
    if not slit:
        cell_list = make_cell_list(Lx, Ly, Lz, False, lj_vext=True)
        if mW:
            rc = 1.8
        else:
            rc = 2.5

    # Calculate total number of cells.
    num_cells = Lx*Ly*Lz

    # Set up particle positions array.
    positions = np.zeros((particles,5))

    for i in range(particles):

        # If using the slit geometry, a random particle position
        # can be generated by picking a cell and position within
        # the cell at random.
        if slit:
            positions[i,0] = i
            positions[i,1] = np.random.random_sample()
            positions[i,2] = np.random.random_sample()
            positions[i,3] = np.random.random_sample()
            positions[i,4] = np.random.random_integers(low=0, high=num_cells)

        # If using a solute geometry, the real-space radial position
        # of the particle with respect to the centre of the solute
        # must be found. This is then compared to the radius of the
        # solute to ensure the particle falls outside the solute.
        else:
            valid = False
            while not valid:
                # Generate a particle at a random position.
                x = np.random.random_sample();
                y = np.random.random_sample();
                z = np.random.random_sample()
                cell = np.random.random_integers(low=0, high=num_cells)

                # Calculate the radial distance of the particle.
                dx = (x+cell_list[int(cell),0]) - Lx/2.
                dy = (y+cell_list[int(cell),1]) - Ly/2.
                dz = (z+cell_list[int(cell),2]) - Lz/2.
                r = rc*np.sqrt(dx*dx + dy*dy + dz*dz)

                # If this is greater than the radius of the solute, use
                # the particle.
                if r > Rs:
                    valid = True
                    positions[i,0] = i
                    positions[i,1] = x; positions[i,2] = y; positions[i,3] = z
                    positions[i,4] = cell

    # Write the particle positions to file in the correct format
    # for the Monte Carlo program.
    with open(os.path.join('Tests','lj_external','positions_final'), 'w') as out:
        for i in range(particles):
            out.write(f'{int(positions[i,0])} {positions[i,1]} {positions[i,2]} {positions[i,3]} {int(positions[i,4])}\n')

    return positions



