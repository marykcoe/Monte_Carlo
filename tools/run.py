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

This file contains tools for running simulation with the program.
"""

import os
import subprocess
import time
import sys

def run(inp):

    """
    Calls Monte Carlo program. Leaves 2 second delay between
    launching simulations to ensure different seeds. Records
    the time taken for each simulation to run and writes to
    file.

    Parameters
    -----------
    inp : list(string)
        File paths to Input file(s) for simulation(s).
    """

    # Start simualtions leaving a delay in between.
    processes = []; timings = [];
    for i in inp:
        timings.append(time.time())
        processes.append(subprocess.Popen(['./MC', i]))
        time.sleep(2)

    # Wait for simulations to finish
    for p in processes:
        p.wait()

    # Output timings
    if len(inp)>1:
        fn = inp[0][:-2] + 'timings'
    else:
        fn = inp[0] + 'timings'

    with open(fn, 'a') as out:
        for t in range(len(timings)):
            out.write(f'Simulation {t} completed in {(time.time()-timings[t])/3600.0} hours\n')

def write_input_file(ofp,Lx, Ly, Lz, fluid_type, chem_pot, temperature, sweeps,
        equilibrate = True, save_freq = 10, dist_freq = 200,
        compressibility = False, susceptibility = False, mc = False,
        tmmc=True, load_weights = False, sub_volume = None,
        initial_dens = 0.0, max_dens = None, min_dens = None,
        load_mode=0, ifp = None, slit = False, solute=None,
        lj_vext = False, epsilon_wall = 0.0, test = None,
        output_individual_dist = True, profile_bins=200):

    """
    Writes the Input file for a simulation.
    Creates the necessary output file directories.

    Parameters
    ----------
    ofp : string
        output file path
    Lx : int
        Number of cells in x-dimension of box
    Ly : int
        Number of cells in y-dimension of box
    Lz : int
        Number of cells in z-dimension of box
    fluid_type : string
        Type of fluid to simulate. Truncated Lennard-Jones (fluid_type = LJ)
        and monatomic water (fluid_type = mW) supported.
    chem_pot : float
        chemical potential
    temperature (: float
        temperature
    sweeps : int
        Number of Monte Carlo sweeps to perform. A sweep is
        defined as one move attempt per number of cells.
    load_mode : int, optional
        Determines whether to load a previous simulation or an
        initial particle positions file.
        Options are:
            0 = do not load any previous sim
            1 = load previous sim
            2 = load a previous configuration of particles and
                run a new sim
        The default is 0.
    ifp : string, optional
        File path to load a previous simulation or particle positions
        from. IThis is a required field is load = 1 or 2.
    slit : bool, optional
       Determines whether to use a slit geometry. The default is False.
    solute : None/float, optional
        Determines whether to use a solute geometry. If None, the solute
        geometry is not used.  If a float, then the float is the radius
        of the solute in units of fluid particle diameters. The default is
        None.
    initial_dens : float, optional
        Initial density of fluid. If not loading a previous simulation
        or configuration, this determines the number of particles to
        initially populate the system with. This should be given in units
        of gcm^-3 for a monatomic water fluid , and as a number density
        for a truncated Lennard-Jones fluid.
    max_dens : float
        Maximum fluid density allowed in system. This should be given
        in units of gcm^-3 for a monatomic water fluid, and as a number
        density for a truncated Lennard-Jones fluid.
    min_dens : float, optional
        Minumum fluid density allowed in system. This should be given
        in units of gcm^-3 for a monatomic water fluid, and as a number
        density for a truncated Lennard-Jones fluid.
    equilibrate (: bool, optional
        Determines whether to run the simulation for a short period
        before sampling properties.. This should be set to True when
        starting a new simulation or when running a previous simulation
        at a new state point. The number of equilibration sweeps is
        hard-coded to be 20% of the total number of sweeps.
        The default is True.
    save_freq : int, optional
        Number of sweeps to leave between sampling the density and
        energy of the system. The default is 10.
    dist_freq : int, optional
        Number of sweeps to leave between sampling density and
        fluctuation distributions. The default is 200.
    profile_bins : int, optional
        The number of histogram bins to use when calculating particle
        distributions. The default is 200.
    compressibility : bool, optional
        Determines whether to measure the local compressibility during the
        simulation. The default is False.
    susceptibility : bool, optional
        Determines whether to measure the local thermal susceptibility
        as the simulation progresses. The default is False.
    mc : bool, optional
        Determines whether to use multicanonical sampling. The default
        is False.
    tmmc : bool, optional
        Determines whether to use transition matrix methods to
        calculate the weights when using multicanonical sampling.
        The default is True.
    load_weights : bool, optional
        Determines whether to use a precalculated weight distribution
        when using multicanonical sampling. The weight function must
        be saved to the Input file directory. The default is False.
    sub_volume : None/int, optional
        Determines whether to use sub-volume sampling.If None,
        sub-volume sampling is set to False. If an int is supplied,
        then this will represent the number of cells to use for
        sub-volume sampling. This parameter is only valid when using
        multicanonical sampling. The default is None.
    lj_vext : bool, optional
        Determines whether to apply an external potential, originating
        from the walls of the slit or solute, to the fluid. The default
    epsilon_wall : float, optional
        The strength of the surface-fluid interaction, in units of
        fluid-fluid interaction strength, if using an external
        potential. The default is 0.0.
    test : None/string, int, optional
        Determines whether the simulation is a test simulation, and
        if so what type. If None is supplied, the simulation is not
        a test. If a string is supplied, then the string represents
        the name of the test. For these types of test, please refer
        to the test module. If an int is supplied, then the int
        acts as the seed for the random number generator.

    """


    # Set up output directories

    if not os.path.exists(ofp):
        os.makedirs(ofp)

    if not os.path.exists(os.path.join(ofp, 'distributions')):
        os.makedirs(os.path.join(ofp, 'distributions'))

    if not os.path.exists(os.path.join(ofp ,'averaged_distributions')):
        os.makedirs(os.path.join(ofp, 'averaged_distributions'))

    if not os.path.exists(os.path.join(ofp ,'particle_positions')):
        os.makedirs(os.path.join(ofp, 'particle_positions'))

    # If using the transition matrix method with multicanonical
    # sampling, create directories to store necessary matrices.
    if mc and tmmc:
        if not os.path.exists(os.path.join(ofp , 'collection_matrices')):
            os.makedirs(os.path.join(ofp, 'collection_matrices'))
        if not os.path.exists(os.path.join(ofp, 'transition_matrices')):
            os.makedirs(os.path.join(ofp, 'transition_matrices'))
        if not os.path.exists(os.path.join(ofp, 'weight_distributions')):
            os.makedirs(os.path.join(ofp, 'weight_distributions'))

    # Write Input file
    with open(os.path.join(ofp, 'Input'), 'w') as out:

        # Input/Output parameters
        out.write('# Input/Output Parameters\n')
        out.write(f'OUTPUT_FILE_PATH {ofp}\n')
        if load_mode == 0:
            out.write('LOAD FALSE\n')
        elif load_mode == 1:
            out.write('LOAD TRUE\n')
        elif load_mode == 2:
            out.write('LOAD POSITIONS\n')
        else:
            print('Invalid load mode. Load set to False.')
            out.write('LOAD FALSE\n')

        if (load_mode>0):
            if ifp is None:
                print('Error. Cannot load previous simulation without knowledge')
                print('of file path to load from.')
                sys.exit(0)
            else:
                out.write(f'LOAD_FILE_PATH {ifp}\n')

        # Simulation box parameters
        out.write('\n# Box Parameters\n')
        out.write(f'LX {Lx}\n')
        out.write(f'LY {Ly}\n')
        out.write(f'LZ {Lz}\n')

        if (slit):
            if solute is not None:
                print('Error. Cannot have a solute and slit geometry.')
                sys.exit(0)
            else:
                out.write('SLIT TRUE\n')
        else:
            out.write('SLIT FALSE\n')

        if (solute is not None):
            out.write(f'SOLUTE {solute}\n')
        else:
            out.write('SOLUTE FALSE\n')

        # Fluid parameters
        out.write('\n# Fluid Parameters\n')
        out.write(f'FLUID_TYPE {fluid_type}\n')
        out.write(f'CHEMICAL_POTENTIAL {chem_pot}\n')
        out.write(f'TEMPERATURE {temperature}\n')
        if (initial_dens is not None):
            out.write(f'INITIAL_DENSITY {initial_dens}\n')
        if (min_dens is not None):
            out.write(f'MIN_DENSITY {min_dens}\n')
        if (max_dens is not None):
            out.write(f'MAX_DENSITY {max_dens}\n')

        # Sampling parameters
        out.write('\n# Sampling parameters\n')
        out.write(f'SWEEPS {sweeps}\n')
        if equilibrate:
            out.write('EQUILIBRATE TRUE\n')
        else:
            out.write('EQUILIBRATE FALSE\n')
        out.write(f'SAVE_FREQUENCY {save_freq}\n')
        out.write(f'DISTRIBUTION_FREQUENCY {dist_freq}\n')
        if output_individual_dist:
            out.write('OUTPUT_INDIVIDUAL_DISTRIBUTIONS TRUE\n')
        else:
            out.write('OUTPUT_INDIVIDUAL_DISTRIBUTIONS FALSE\n')
        out.write(f'PROFILE_BINS {profile_bins}\n')
        if compressibility:
            out.write('COMPRESSIBILITY TRUE\n')
        else:
            out.write('COMPRESSIBILITY FALSE\n')
        if susceptibility:
            out.write('SUSCEPTIBILITY TRUE\n')
        else:
            out.write('SUSCEPTIBILITY FALSE\n')
        if mc:
            out.write('MC TRUE\n')
        if tmmc:
            out.write('TMMC TRUE\n')
        else:
            out.write('TMMC FALSE\n')
        if load_weights:
            out.write('LOAD_WEIGHTS TRUE\n')
        else:
            out.write('LOAD_WEIGHTS FALSE\n')
        if sub_volume is not None:
            out.write(f'SUB_VOLUME {sub_volume}\n')
        else:
            out.write('SUB_VOLUME FALSE\n')

        # External potential parameters
        if (solute is not None) or (slit is True):
            out.write('\n# External Potential Parameters\n')
            if lj_vext:
                out.write('LJ_VEXT TRUE\n')
                out.write(f'SURFACE_INTERACTION_STRENGTH {epsilon_wall}\n')
            else:
                out.write('LJ_VEXT FALSE\n')

        # Test parameters
        if test is not None:
            if isinstance(test, int):
                out.write('TEST {test}\n')
            else:
                print('Error. This function only supports test seeds. For other')
                print('tests please use the test function. Test will be set to False.')

def organise_inputs(nprocs, fp, params, ifp = None, individual = True, wfp = None):

    """
    Organisises directories and Input files when running multiple
    simulations simultaneously.

    Parameters
    ----------
    nprocs : int
        Number of simulations which will run simultaneously. Each
        simulation will require its own processor and hence this also
        reflects the number of processors required.
    fp : (string)
       Path to directory to save outputs of all simulations.
    params : (dict)
        Dictionary which stores all input parameters for simulation.
        Required keys are {Lx, Ly, Lz, fluid_type, chem_pot, temperature,
        sweeps}. For information on the meaning of each key, and further
        optional keys, please see the 'write_input_file' function.
    ifp : None/string, optional
        Path to directory containing the results of a previous
        set of simulations, if continuing a previous simulation, this
        The default is None.
    individual : bool, optional
        Determines whether simulations should be loaded from the same
        previous simulation or different ones. If True, simulations
        will be loaded from ifp + sim/. If False, simulations will be
        loaded from the ifp directory. The default is True.
    wfp : None/string, optional
        Path to directory containing weight function to be loaded,
        if applicable. The default is None.

    Returns
    -------
    ofps : list(string)
        List of file paths to Input files. Each of these acts as the
        command line argument for the Monte Carlo simulation, which
        directs the program to the folder containing the Input file.
	"""

    # Set up input and output file directories.
    ofps = []; ifps = [];
    if not os.path.exists(fp):
        os.makedirs(fp)

    # For every simulation to run, set up an output file
    # directory at the path fp/processor_number.
    for sim in range(nprocs):
        ofps.append(os.path.join(fp, str(sim)))
        if "load_mode" in params:
            if params["load_mode"] > 0:
                if individual:
                    ifps.append(os.path.join(ifp, str(sim)))


    # If loading a weight function, copy the weights file located in
    # the wfp directory to the output file directory of each
    # simulation.
    if "load_weights" in params:
        if params["load_weights"] is True:
            if wfp is None:
                print('Error. No file path to initial weight function was supplied.')
                sys.exit(0)
            else:
                for sim in range(nprocs):
                    command = 'cp ' + wfp  + ' ' + os.path.join(ofps[sim], 'weights_initial')
                    os.popen(command)

    # Write the Input file for each simulation and save to the
    # appropriate output file directory. This also creates any
    # necessary directories for the simulation output files.
    if "load_mode" in params:
        if params["load_mode"]>0:
            for ofp, ifp in zip(ofps,ifps):
                write_input_file(ofp, **params, ifp = ifp)
        else:
            for ofp in ofps:
                write_input_file(ofp, **params)
    else:
        for ofp in ofps:
            write_input_file(ofp, **params)

    return ofps

def write_restart_output_file(ifp, last, collect_params = True,
                              params = None):

    """
    Writes an Output file for a prematurely ended simulation to the
    file directory 'ifp'. Requires the Input, data and seed files from
    the prematurely ended simulation.

    Parameters
    ----------
    ifp : string
        Path to directory containing input and output files of the
        prematurely ended simulation.
    last : int
        Sweep number of last saved particle configuration of the
        prematurely ended simulation.
    collect_params : bool, optional
        If True, a dictionary of simulation parameters is collected
        and returned. The default is True.
    params : dict, optional
        Dictionary containing parameters of prematurely ended
        simulation.

    Returns
    -------
    params : dict, optional
        Dictionary of parameters of prematurely ended simulation.

    """

    # If collecting the simulation parameters, create a label
    # conversion dictionary and an empty dictionary to store
    # the parameters.
    if collect_params:
        params = {}
        label_conversion = {"LX": "Lx", "LY": "Ly", "LZ":"Lz","SLIT":"slit",
            "SOLUTE": "solute", "TMMC": "tmmc", "FLUID_TYPE":"fluid_type",
            "CHEMICAL_POTENTIAL":"chem_pot", "TEMPERATURE":"temperature",
            "SAVE_FREQUENCY":"save_freq", "DISTRIBUTION_FREQUENCY":"dist_freq",
            "MAX_DENSITY": "max_dens", "MIN_DENSITY": "min_dens",
            "COMPRESSIBILITY": "compressibility", "SUSCEPTIBILITY":"susceptibility",
            "MC": "mc", "LOAD_WEIGHTS": "load_weights","SUB_VOLUME": "sub_volume",
            "LOAD_MODE": "load_mode", "TEST": "test", "EQUILIBRATE": "equilibrate",
            "OUTPUT_INDIVIDUAL_DISTRIBUTIONS": "output_individual_dist",
            "INITIAL_DENSITY": "initial_dens", "LJ_VEXT": "lj_vext",
            "EPSILON_WALL": "epsilon_wall", "SURFACE_INTERACTION_STRENGTH": "epsilon_wall",
            "TEST": "test", "PROFILE_BINS": "profile_bins"}

    # Make an Output file in the ifp directory and write the
    # necessary parameters.
    with open(os.path.join(ifp, "Output"), 'w') as out:

        # Write parameters which appear in both the Input and
        # Output files
        with open(os.path.join(ifp ,"Input"), 'r') as inp:
            for l, line in enumerate(inp):
                split_line = line.split()
                if (len(split_line) > 0):
                    if split_line[0] == 'OUTPUT_FILE_PATH':
                        continue;
                    elif split_line[0] == "LOAD":
                        continue;
                    elif split_line[0] == "INPUT_FILE_PATH":
                        continue;
                    elif split_line[0] == "SWEEPS":
                        out.write(f'SWEEPS {last}\n')
                        out.write(f'TOTAL_SWEEPS {last}\n')
                        params["sweeps"] = float(split_line[-1]) -last
                    else:
                        if collect_params:
                            if split_line[0] != '#' and split_line[0] != 'LOAD_FILE_PATH':
                                label = label_conversion[split_line[0]]
                                if split_line[0] == "SOLUTE":
                                    if split_line[-1]=='FALSE':
                                        params[label] = None
                                    else:
                                        params[label] = float(split_line[-1])
                                elif split_line[-1] == 'TRUE':
                                    params[label] = True
                                elif split_line[-1] == 'FALSE':
                                    params[label] = False
                                elif split_line[-1] == 'mW' or split_line[-1] == 'LJ':
                                    params[label] = split_line[-1]
                                else:
                                    params[label] = float(split_line[-1])
                        if split_line[0] != "LOAD_FILE_PATH":
                            out.write(line)

    # If collecting parameters, check for missing required
    # entries and then return the collected dictionary.
    if collect_params:
        if "mc" not in params:
            params["mc"] = False
        return params

def restart_simulation(nprocs, ofp, ifp, last, collect_params = True,
                       params = None, sweeps = None):

    """
    Organises Input/Output files and file directories in order to
    restart a simulation which ended prematurely.

    Parameters
    ----------
    nprocs : int
        Number of simulations which within the same main file directory
        which ended prematurely. Each simulation will have required its
        own processor and hence this also reflects the number of
        processors which were used.
    ofp : string
        Path to directory to save outputs of the restarted simulation.
        It is recommended that this be different from the directory
        of the prematurely ended simulation.
    ifp : string
        Path to directory containing input and output files of the
        prematurely ended simulation.
    last : int
        Sweep number of last saved particle configuration.
    params : None/dict, optional
        Dictionary containing parameters of prematurely ended
        simulation. If None, parameters will be collected from Input
        file of the prematurely ended simulation. The default is None.
    sweeps : None/int, optional
        Number of sweeps to complete in restarted simulation. If None,
        the number of sweeps that were not completed in prematurely
        ended simulation will be used.

    Returns
    -------
    ofps : list(string)
        List of paths to directories containing Input files for the
        simulations to restart. Each of these acts as the command line
        argument for the Monte Carlo simulation, which directs the
        program to the folder containing the Input file.
	"""

    # Set up file directories for simulations to restart.
    ofps = []; ifps = [];
    if not os.path.exists(ofp):
        os.makedirs(ofp)



    # Loop through all simulations to restart. In each case,
    # write the necessary Input file as well as the Output file
    # of the corresponding simulation which ended prematurely.
    # Rename files in the prematurely ended simulation to reflect
    # that they were the final files.
    for sim in range(nprocs):

        # Names of the prematurely ended simulation folder
        # and the restarted simulation folder.
        sifp = os.path.join(ifp,str(sim))
        sofp = os.path.join(ofp,str(sim))

        # Collate names of the file directories of the prematurely
        # ended simulations and the simulations to restart.
        ofps.append(sofp)
        ifps.append(sifp)

        # Make the output folder for the restarted simulation.
        if not os.path.exists(sofp):
            os.makedirs(sofp)

        # Write the Output file for the prematurely ended simulation.
        # Collect simulation parameters at the same time if applicable.
        if collect_params:
            params = write_restart_output_file(sifp, last,
                                               collect_params=True)
            params["load_mode"] = 1
            params["ifp"] = sifp
        else:
            write_restart_output_file(sifp, last,
                                      collect_params=False,
                                      params = params)

        # Set the number of sweeps to complete for the restarted
        # simulation.
        if sweeps is not None:
            params["sweeps"] = sweeps

        # Restarted simulations are assumed to be equilibrated.
        params["equilibrate"]=False

        # Write the Input file for the restarted simulation.
        write_input_file(sofp, **params)

        # When loading a simulation, the Monte Carlo program looks
        # for certain files - normally ending with '_final'. Here we
        # copy the 'last' sweep files such that they now have the
        # appropriate ending to their file name.
        command = 'cp ' + sifp + os.path.join('particle_positions'',''positions_' \
                    + str(last)) + ' ' + os.path.join(sifp ,'positions_final')
        os.popen(command)

        if last > params["dist_freq"]:
            command = 'cp ' + sifp + os.path.join('averaged_distributions','average_distribution_' + str(last) + '_raw')
            command += ' ' + sifp + os.path.join("averaged_distributions","average_distribution_final_raw")
            os.popen(command)

        if params["compressibility"]:
            command = 'cp ' + os.path.join(sifp,'averaged_distributions','local_compressibility_' + str(last) +' _raw')
            command += ' ' + os.path.join(sifp,"averaged_distributions","local_compressibility_final_raw")
            os.popen(command)

        if params["susceptibility"]:
            command = 'cp ' + os.path.join(sifp, 'averaged_distributions','local_thermal_susceptibility_' + str(last) + '_raw')
            command += ' ' + os.path.join(sifp, 'averaged_distributions','local_thermal_susceptibility_final_raw')
            os.popen(command)

        if params["mc"]:
            if params["tmmc"]:
                command = 'cp ' + os.path.join(sifp, 'collection_matrices','collection_' + str(last))
                command += ' ' + os.path.join(sifp, 'collection_final')
                os.popen(command)

                command = 'cp ' + os.path.join(sifp, 'weight_distributions','weights_' + str(last))
                command += ' ' + os.path.join(sifp, 'weights_final')
                os.popen(command)
            else:
                command = 'cp ' + os.path.join(sifp, 'weights_initial')
                command += ' ' + os.path.join(sofp, 'weights_initial')
                os.popen(command)
                if params is None:
                    params["tmmc"] = False

        elif params is None:
            params["tmmc"] = False

    return ofps
