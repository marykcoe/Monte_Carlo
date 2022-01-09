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

This file contains tools for running tests on the program.
"""
import os
import sys

import tools.run as run
import tools.make_test_files as mtf

import filecmp
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors

def metropolis():

    """
    Runs the Metropolis test on the program.

    Purpose: To test the underlying Monte Carlo algorithm
    is working as expected. This test can also be used to
    test the incorporation of a random number generator has
    been successful.

    Method: A simulation is run for 0 fluid interaction
    strength, making the fluid effectively an ideal gas.
    For an ideal gas, the average density of the system
    should obey the equation
            chemical potential = ln(average density).
    If this is obeyed during the simulation, then the test
    is passed. If not, the user should check the random
    number generator.
    """

    # Set up output file directories.
    ofp = os.path.join('Tests','metropolis')
    if not os.path.exists(ofp):
        os.makedirs(ofp)

    # Write Input file
    with open(os.path.join(ofp, 'Input'), 'w') as out:

        # Input/Output parameters
        out.write('# Input/Output Parameters\n')
        out.write(f'OUTPUT_FILE_PATH {ofp}\n')

        # Test Parameters
        out.write('\n# Test Parameters\n')
        out.write('TEST METROPOLIS\n')

    # Run test
    run.run([ofp])

def critical_point(mW):

    """
    Runs critical point test on program.

    Purpose: Critical points can be a very sensitive measure of whether
    the fluid potential is correct. For the fluids studied within this
    program, they are also well known. Hence, running a simulation
    at the critical point parameters should return a probability
    distribution which exhibits two peaks, one for liquid and one for
    vapour, and should have an average density within a specified range.

    Method: Runs a simulation at the critical point parameters for a
    specified fluid type. Density and energy statsitics are saved as the
    simulatio progresses. At the end of the simultion, the average
    density is printed to screen, along with that expected. The accompanying
    Python program then plots the density and energy distribution.
    If the average density is similar to that expected, and the density
    and energy distributions exhibit two peaks, then it is likely that
    the fluid is behaving as expected and the test is therefore passed.
    There are two things to note. Firstly, the user will have to determine
    that there is a double peaked distribution. Secondly, this distribution
    is not expected to be perfect and the average density calculated
    is not expected to be exactly equal to the reported result. Reasons for
    this are twofold - (i) the simulation run is very short and cannot
    be expected to provide super reliable statistics and (ii) finite size
    effects are not taken into account (we quote the infinite system
    critical density).

    Parameters
    ----------
    mW : bool
        If True, the critical point test is performed for the
        monatomic water fluid. If False, the test is performed
        for the truncated Lennard-Jones fluid.
    """

    # Set up output file directories.
    if mW:
        ofp = os.path.join('Tests','mW_critical_point')
    else:
        ofp = os.path.join('Tests','lj_critical_point')
    if not os.path.exists(ofp):
        os.makedirs(ofp)

    # Write Input file
    with open(os.path.join(ofp, 'Input'), 'w') as out:

        # Input/Output parameters
        out.write('# Input/Output Parameters\n')
        out.write(f'OUTPUT_FILE_PATH {ofp}\n')

        # Test Parameters
        out.write('\n# Test Parameters\n')
        out.write('TEST CRITICAL_POINT\n')

        # Fluid Parameters
        out.write('\n# Fluid Parameters\n')
        if mW:
            out.write('FLUID_TYPE mW\n')
        else:
            out.write('FLUID_TYPE LJ\n')

    # Run test.
    run.run([ofp])

    # Plot distributions.
    if mW:
        volume = 0.0
    else:
        volume = np.power(5.*2.5,3.)
    data = np.genfromtxt(ofp + 'data', usecols=(2,3))

    densy, densx = np.histogram(data[:,0], bins = int(np.amax(data[:,0])-np.amin(data[:,0])),
                    range = (int(np.amin(data[:,0])),int(np.amax(data[:,0]))),
                    density=True)

    eny,enx = np.histogram(data[:,1], bins = 200,
                        range = (np.amin(data[:,1]),np.amax(data[:,1])), density=True)

    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.plot(densx[:-1]/volume, densy[:], color='cornflowerblue', ls=':', marker='o',
          markeredgecolor='cornflowerblue', markerfacecolor='None', lw=2,
          markersize=3)
    if mW:
        ax1.set_xlabel(r'$\rho\;(gcm^{-3})$')
    else:
        ax1.set_xlabel(r'$\rho\;(\sigma^{-3})$');
    ax1.set_ylabel(r'$P(\rho)$')
    ax1.set_xlim(np.amin(data[:,0])/volume, np.amax(data[:,0])/volume)
    ax1.tick_params(direction='in', top='True', right='True', which='both')
    plt.savefig(os.path.join(ofp, 'density_distribution.pdf'))
    plt.show()

    fig = plt.figure(2)
    ax2 = fig.add_subplot(111)
    ax2.plot(enx[:-1], eny[:], color='royalblue', ls=':', marker='s',
          markeredgecolor='royalblue', markerfacecolor='None', lw=2,
          markersize=3)
    ax2.set_xlabel(r'$\beta E$'); ax2.set_ylabel(r'$P(\beta E)$')
    ax2.set_xlim(np.amin(data[:,1]), np.amax(data[:,1]))
    ax2.tick_params(direction='in', top='True', right='True', which='both')
    plt.savefig(os.path.join(ofp, 'energy_distribution.pdf'))
    plt.show()


def cell_list(Lx=5, Ly=5, Lz=5, slit=False,  sv = False, dz = None):

    """
    Performs cell list test. User may optionally provide parameters
    of simulation box or use defaults.

    Purpose: To test that the cell list has been set up correctly upon
    start-up.

    Method: A simulation box is setup for a user defined Lx, Ly and Lz.
    The cell list is generated and saved to file. Information saved
    includes position of cell within cell list, location of cell in
    real-space and neighbours of cell. The program then exits. The
    accompanying Python program also creates a cell list. That output
    from the program is then compared to this Python generated version.
    If they match, then the test is passed.

    Parameters
    ----------

    Lx : int, optional
        Number of cells in x-direction. The default is 5.
    Ly : int, optional
        Number of cells in y-direction. The default is 5.
    Lz : int, optional
        Number of cells in z-direction. The default is 5.
    slit : bool, optional
        Determines whether the geometry of the system
        should be a slit (True) or a solute/bulk (False).
        The default is False.
    sv : bool, optional
        Determines whether sub-volume sampling is enabled.
        The default is False.
    dz : int, optional
        If sub-volume sampling is enabled, dz must be supplied.
        This represents the number of cells in the z-direction
        which form the sub-volume to sample. The default is None.
    """

    # Set up output file directories.
    ofp = os.path.join('Tests','cell_list')
    if not os.path.exists(ofp):
        os.makedirs(ofp)

    # Output test cell list
    mtf.make_cell_list(Lx, Ly, Lz, slit,
                                lj_vext=False, sv = sv, dz = dz)

    # Write Input file
    with open(os.path.join(ofp, 'Input'), 'w') as out:

        # Input/Output parameters
        out.write('# Input/Output Parameters\n')
        out.write(f'OUTPUT_FILE_PATH {ofp}\n')

        # Box parameters
        out.write('\n# Box Parameters\n')
        out.write(f'LX {Lx}\n')
        out.write(f'LY {Ly}\n')
        out.write(f'LZ {Lz}\n')

        if slit:
            out.write('SLIT TRUE\n')

        # Sampling parameters
        if sv:
            out.write('\n# Sampling Parameters\n')
            out.write(f'SUB_VOLUME {dz}\n')

        # Test Parameters
        out.write('\n# Test Parameters\n')
        out.write('TEST CELL_LIST\n')

    # Run test
    run.run([ofp])

    # Compare cell_lists and print result of test.
    same = filecmp.cmp(os.path.join(ofp, 'cell_list'),
                       os.path.join(ofp, 'test_cell_list'))
    print('\n----------------------')
    if same:
        print('CELL LIST TEST PASSED')
    else:
        print('CELL LIST TEST FAILED')
    print('----------------------')

def lj_external_potential(Lx=5, Ly=5, Lz=10, slit=True, mW=True,
                          T=426, surface_fluid_attraction=0.1, Rs=None,
                          particles = 30):

    """
    Runs the Lennard-Jones external potential test.

    Purpose: To confirm that the LJ external potential is calculating
    correctly.

    Method: Predetermined particle positions are read in for a given
    geometry (solute or slit) and the external potential felt by each
    particle is calculated. This is then compared against a set of
    precalculated (from python) external potentials.

    Parameters
    ----------
    Lx : int, optional
        Number of cells in x-direction. The default is 5.
    Ly : int, optional
        Number of cells in y-direction. The default is 5.
    Lz : int, optional
        Number of cells in z-direction. The default is 10.
    slit : bool, optional
        Determines whether to use the slit geometry. The default
        is True.
    mW : bool, optional
        Determine whether to use a monatomix water fluid (True)
        or a truncated Lennard-Jones fluid (False). The default
        is True.
    T : float, optional
        The temperature of the system. The default is 426K for the
        monatomic water fluid.
    surface_fluid_interaction : float, optional
        The strength of the surface-fluid interaction in units of the
        strength of the fluid-fluid interaction.
    Rs : float, optional
        Radius of the solute if using the solute geometry. The
        default is None.
    particles : int, optional
        The number of particles to calculate the potential for.

    """

    # Set up the output file directories.
    ofp = os.path.join('Tests', 'lj_external')
    if not os.path.exists(ofp):
        os.makedirs(ofp)

    # The solute geometry is used is a radius of the
    # solute is supplied.
    if not slit and (Rs is None):
        print('Error. If using solute geometry, you must supply the ' +
              'radius of the solute, Rs.')
        sys.exit(0)

    # Set up fluid parameters.
    if mW:
        rc = 1.8; sigma = 2.3925
        epsilon = 6.189 * 4184./(6.02214076e23);
        epsilon /=(T*1.38064852e-23)
    else:
        rc = 2.5; sigma = 1.0
        epsilon = 4./T

    epsilon_wall = surface_fluid_attraction*epsilon

    # Generate random positions for particles and output to file.
    positions = mtf.make_lj_vext_test_files(Lx, Ly, Lz, epsilon, epsilon_wall,
                                slit = slit, particles = particles, mW=mW, Rs=Rs)

    mtf.make_cell_list(Lx,Ly,Lz,slit,lj_vext=True)

    # Write Input file
    with open(os.path.join(ofp, 'Input'), 'w') as out:

        # Input/Output parameters
        out.write('# Input/Output Parameters\n')
        out.write(f'OUTPUT_FILE_PATH {ofp}\n')
        out.write('LOAD POSITIONS\n')
        out.write(f'LOAD_FILE_PATH {ofp}\n')

        # Box parameters
        out.write('\n# Box Parameters\n')
        out.write(f'LX {Lx}\n')
        out.write(f'LY {Ly}\n')
        out.write(f'LZ {Lz}\n')

        if slit:
            out.write('SLIT TRUE\n')
        else:
            out.write(f'SOLUTE {Rs}\n')

        # Fluid parameters
        out.write('\n# Fluid Parameters\n')
        if mW:
            out.write('FLUID_TYPE mW\n')
        else:
            out.write('FLUID_TYPE LJ\n')
        out.write(f'TEMPERATURE {T}\n')

        # External potential
        out.write('\n# External Potential Parameters\n')
        out.write('LJ_VEXT TRUE\n')
        out.write(f'SURFACE_INTERACTION_STRENGTH {surface_fluid_attraction}\n')

        # Test Parameters
        out.write('\n# Test Parameters\n')
        out.write('TEST LJ_EXTERNAL\n')

    # Run test
    run.run([ofp])

    pyext = np.zeros((particles, 2))

    # For each particle, calculate the external potential felt.
    for p in range(particles):

        if slit:
            z = rc*(positions[p,3] + cell_list[int(positions[p,4]),2])
            pyext[p,0] = z; pyext[p,1] = epsilon_wall*mtf.lj_vext(z,Lz=Lz*rc)
        else:
            dx = (positions[p,1]+cell_list[int(positions[p,4]),0]) - Lx/2.
            dy = (positions[p,2]+cell_list[int(positions[p,4]),1]) - Ly/2.
            dz = (positions[p,3]+cell_list[int(positions[p,4]),2]) - Lz/2.
            r = rc*np.sqrt(dx*dx + dy*dy + dz*dz)
            pyext[p,0] = r
            pyext[p,1] = epsilon_wall*mtf.lj_vext(r, planar=False, Rs=Rs)

    # Calculate continous external potential for plotting.
    if slit:
        z = np.linspace(0.0,Lz*rc, num=100)
        potential = np.zeros((100,2))
        for p in range(100):
            potential[p,0] = z[p]
            potential[p,1] = epsilon_wall*mtf.lj_vext(z[p], Lz=Lz*rc)
    else:
        r = np.linspace(Rs,(Lz/2.)*rc, num=100)
        potential = np.zeros((100,2))
        for p in range(100):
            potential[p,1] = epsilon_wall*mtf.lj_vext(r[p],planar=False, Rs=Rs)
            potential[p,0] = r[p]

    # Write python calculated external potential for each particle
    # to file for comparison to C++ program.
    with open(os.path.join(ofp, 'test_external_potential'), 'w') as out:
        for i in range(particles):
            out.write(f'{pyext[i,0]:.7f} {pyext[i,1]:.7f}\n')

    # Read in external potential from C++ program and
    # calculate difference.
    cext = np.genfromtxt(os.path.join(ofp , 'external_potential'))
    difference = abs(pyext[:,1]-cext[:,1])

    # Plot results.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if slit:
        minus = 0
    else:
        minus = Rs
    ax.plot(potential[:,0]-minus, potential[:,1], color='rosybrown', ls=':')
    ax.plot(pyext[:,0]-minus, pyext[:,1], ls='None', marker='^', markeredgecolor='mediumvioletred',
         markerfacecolor='None', markersize=5, label='Test')
    ax.plot(cext[:,0]-minus, cext[:,1], ls='None', marker='o', markeredgecolor='cornflowerblue',
         markerfacecolor='None', markersize=5, label='Program')
    if slit:
        ax.set_xlabel(r'$z/\sigma$')
        ax.set_xlim(0.0,Lz*rc)
        ax.set_ylabel(r'$V_{ext}(z)$')
    else:
        ax.set_xlabel(r'$(r-R_s)/\sigma$')
        ax.set_xlim(0.0,(Lz/2.)*rc - Rs)
        ax.set_ylabel(r'$V_{ext}(r)$')
    ax.legend(frameon=False, ncol=2)
    ax.set_ylim(1.1*np.amin(potential[:,1]),0.05)
    ax.tick_params(which='both', direction='in', top=True, right=True)
    plt.savefig(os.path.join(ofp, 'external_potential.pdf'))

    # Calculate the maximum difference between the python
    # and C++ calculations.
    max_dev = np.amax(difference)

    # Use this to determine whether the test is passed.
    print('\n------------------------------------------')
    if abs(max_dev) < 1e-7:
        print('EXTERNAL POTENTIAL TEST PASSED')
    else:
        print('EXTERNAL POTENTIAL TEST FAILED')
    print('------------------------------------------')









