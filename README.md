# Cellular Automata for Quantum Mechanics
To run the simulation you need a working installation of [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) 3.4.0 or newer. Parameters for the simulation e.g. the size of the Lattice, are received in compile time so you must change them by directly editing qlb.cpp in the src folder. Compilation using the optimization flag -O3 of g++ increases considerably the efficiency and so we advised users to compile using it.

The project inspired on the work of David A. Meyer on [Quantum Cellular Automatas](https://arxiv.org/abs/quant-ph/9604003), the modification to it proposed by Eduardo Ortega (unpublished) and is an implementation of the latter.
