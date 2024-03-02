#!/bin/bash

# Change to the hmc directory and run make
cd "$SRC_DIR"/hybridmc || exit
make

# run setup.py
$PYTHON ../setup.py install