#!/bin/bash

# Change to the hmc directory and run make
cd "$SRC_DIR"/hybridmc || exit
make

# run setup.py
cd ../ && $PYTHON ../setup.py install
