#!/bin/bash

mkdir -p "crambin_fixed_e_3.0"
cd "crambin_fixed_e_3.0"
../constructPmatrix ../crambin_fixed_e.dat ../crambin_fixed_e.mtx 3.0
../disconnect/src/disconnectionDPS
convert tree.ps tree.3.0.gif
cd ../

mkdir -p "crambin_fixed_e_5.5"
cd "crambin_fixed_e_5.5"
../constructPmatrix ../crambin_fixed_e.dat ../crambin_fixed_e.mtx 5.5
../disconnect/src/disconnectionDPS
convert tree.ps tree.5.5.gif
cd ../

mkdir -p "crambin_fixed_e_7.5"
cd "crambin_fixed_e_7.5"
../constructPmatrix ../crambin_fixed_e.dat ../crambin_fixed_e.mtx 7.5
../disconnect/src/disconnectionDPS
convert tree.ps tree.7.5.gif
cd ../
