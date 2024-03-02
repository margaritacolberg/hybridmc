#!/bin/bash

# Change to the hmc directory and run make
cd hybridmc || exit
make release

# Enter back to top level directory and run setup.py
cd ../
$PYTHON ../setup.py install