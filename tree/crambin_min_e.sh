#!/bin/bash

mkdir -p "crambin_min_e_1.0"
cd "crambin_min_e_1.0"
../constructPmatrix ../crambin_min_e.dat ../crambin_min_e.mtx 1.0
../disconnect/src/disconnectionDPS
convert tree.ps tree.1.0.gif
cd ../
