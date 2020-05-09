Atom-circuit: Solver for ground state and dynamics of circuits made of atoms. Can solve ring-lead systems, Y-junctions and linear chains. Supports both bosons and fermions. Solver is based on iTensor (https://itensor.org/ version 2.1.1) with DMRG and TEBD in C++. It models atomic (and possibly also other condensed matter systems) as a collection of coupled quasi one-dimensional chains.
This tool has been used to reveal Andreev reflections in Y-junctions with bosonic atoms, as well as the absence of the Aharonov-Bohm effect in bosonic rings. Further reference here https://iopscience.iop.org/article/10.1088/2058-9565/ab2e61/meta or here https://arxiv.org/abs/1807.03616
Install iTensor v.2.1.1 (Caution: Does not work with iTensor v3). Then, copy the files as is with the folder structure into the iTensor folder. 
Go to the code folder, then run make. g++ compiler for C++ is recommended. Run the solver in the code folder as ./quenchdynamics inputfileTest, where inputfileTest can be replaced with any other inputfile available. Output is a .dat file, which includes information such as currents and density of the time evolution.
Further explanations are available in the inputfiles itself as well as in the quenchdynamics.cc file.

Goal is to simulate quench dynamics in ring-lead and Y-junctions. Supports bosons and fermions.
First, it calculates ground state via DMRG. You can add a gaussian density perturbation in the source of the leads using the potential variable.
Then, it calculates the dynamics. The potential is quenched at t=0 and the resulting dynamics is simulated. The density and current in time is outputed to file.
The output .dat file is given automatically a name to identify it with the inputfile settings.
The output is as follows:
Line1: time where snapshots where taken. In each line of length numberSnaps
Line2: current bond dimension
Line3: computational time
Line4: current 1 measured
Line5: current 2 measured
Line6: current 3 measured
Line7: Overlap with initial state
Line8: Not used
Line9 to Number of sites +9: density for each site



When using this code, make sure to cite the itensor library as well as the author of this code:

Haug, T., Dumke, R., Kwek, L.-C. & Amico, L. Andreev-Reflection and Aharonov–Bohm dynamics in atomtronic circuits. Quantum Sci. Technol. 4, 045001 (2019).
