#!/bin/bash

# Change to the hmc directory and run make
cd hybridmc
make

# run setup.py
$PYTHON ../setup.py install