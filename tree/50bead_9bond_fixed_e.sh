#!/bin/bash

mkdir -p "50bead_fixed_e_3.0"
cd "50bead_fixed_e_3.0"
../constructPmatrix ../50bead_9bond_fixed_e.dat ../50bead_9bond_fixed_e.mtx 3.0
../disconnect/src/disconnectionDPS
convert tree.ps tree.3.0.gif
cd ../

mkdir -p "50bead_fixed_e_10.0"
cd "50bead_fixed_e_10.0"
../constructPmatrix ../50bead_9bond_fixed_e.dat ../50bead_9bond_fixed_e.mtx 10.0
../disconnect/src/disconnectionDPS
convert tree.ps tree.10.0.gif
cd ../

mkdir -p "50bead_fixed_e_22.0"
cd "50bead_fixed_e_22.0"
../constructPmatrix ../50bead_9bond_fixed_e.dat ../50bead_9bond_fixed_e.mtx 22.0
../disconnect/src/disconnectionDPS
convert tree.ps tree.22.0.gif
cd ../
