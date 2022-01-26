<h2>Grand Canonical Monte Carlo Program for Truncated Lennard-Jones Fluid and Monatomic Water</h2>
<h3>Mary Coe<br>January 2022</h3>

This is a grand canonical Monte Carlo program to simulate truncated Lennard-Jones fluids and a monatomic water model wither in bulk, in contact with a solute or confined to an infinite slit. This program was written as part of my PhD and is made available under the MIT license in the hope that it will be useful for other students in the future.

<h4>How does the program work and what has it been used for?</h4>
Details of exactly how the program works and the theory behind it can be found in Section 4.2 of my thesis, available <a href = https://research-information.bris.ac.uk/ws/portalfiles/portal/304220732/Thesis_Mary_Coe.pdf>here</a>. It currently supports (click the links for details)

<b>Fluids: </b> <a href=https://doi.org/10.1103/PhysRevE.52.602>Truncated Lennard-Jones fluid</a> and <a href = https://doi.org/10.1021/jp805227c>monatomic water</a>.

<b>Geometries: </b>Bulk fluid, small solutes and an infinite slit.

<b>External Potentials: </b>Solutes and the walls of infinite slits are impenetrable to fluid particles. There is an option for these to interact with fluid particles via an attractive Lennard-Jones tail.

<b>Sampling: </b>Insertion and deletion moves, <a href=https://doi.org/10.1103/PhysRevE.52.602>multicanonical biasing</a>, either with preset weights or with weights determined through the course of the simulation using the <a href = https://doi.org/10.1088/0305-4470/28/23/015>transition matrix method</a>.

<b>Measures: </b>Density of particles within system, density profiles (radial distribution function for bulk fluids, distribution between slit walls or distribution around solute), local compressibility and thermal susceptibility profiles (see tutorial 4, section 4.2.10 of <a href = https://research-information.bris.ac.uk/ws/portalfiles/portal/304220732/Thesis_Mary_Coe.pdf>my thesis</a> or <a href = https://doi.org/10.1103/PhysRevLett.125.268004>this paper</a>).

This program has been used for:

<a href = https://doi.org/10.1103/PhysRevLett.128.045501>M. K. Coe, R. Evans and N. B. Wilding, Density depletion and enhanced fluctuations in water near hydrophobic solutes: Identifying the underlying physics. Phys. Rev, Lett. <b>128</b> 045501</a>

M. K. Coe, R. Evans and N. B. Wilding, The coexistence curve and surface tension of a monatomic water model. (Submitted to J. Chem. Phys.)

<a href = https://research-information.bris.ac.uk/ws/portalfiles/portal/304220732/Thesis_Mary_Coe.pdf>M. K. Coe, "Hydrophobicity across length scales: The role of surface criticality", PhD thesis, (University of Bristol, 2021)</a>

<h4>How is it structured?</h4>
The Monte Carlo program is written in C++ for efficiency. This is compiled via a Makefile, so you shouldn't have to interact with the C++ code much (the exception to this is adding a random number generator, see tutorial 1 for details). In addition to the Monte Carlo code, there is also a folder called <i>tools</i>, which contains several Python modules. These modules contain functions to set up Input files, organise and run simulations in parallel, and process, analyse and visualise results. It is recommended that you use these functions, and write any scripts to run the Monte Carlo in Python. 

<h4>What will be available and when?</h4>
The source code for the Monte Carlo, the associated Makefile and a selection of the Python tools are already available within the repository. In addition, there will be several tutorials to help you understand what the program can be used for an how. These are

<b>Tutorial 1 - Getting Started:</b> (available now) This tutorial teaches you how to add a random number generator and test the program is working.

<b>Tutorial 2 - Bulk Fluids:</b>(available soon) This tutorial will teach you how to use the program to study bulk fluids, particularly at the critical point and at liquid-vapour coexistence. It will also introduce biased sampling techniques.

<b>Tutorial 3 - Introducing Solute and Surfaces:</b> (tbc) This tutorial will teach you how to add a solute or a slit to the system and measure density profiles.

<b>Tutorial 4 - Measures:</b> (tbc) This tutorial will explain some of the other measures available within the program and toolkit, and how to implement them.

<h4>Any known problems?</h4>
Current known issues are:

1. At the end of the simulation, the program calculates the average move acceptance ratio incorrectly. 


