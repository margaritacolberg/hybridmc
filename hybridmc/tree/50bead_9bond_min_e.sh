#!/bin/bash

mkdir -p "50bead_min_e_1.0"
cd "50bead_min_e_1.0"
../constructPmatrix ../50bead_9bond_min_e.dat ../50bead_9bond_min_e.mtx 1.0
../disconnect/src/disconnectionDPS
convert tree.ps tree.1.0.gif
cd ../
