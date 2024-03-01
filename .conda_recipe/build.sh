#!/bin/bash

# Change to the hmc directory and run make
cd $RECIPE_DIR/../hybridmc
make

# Go back to the root directory and run setup.py
cd $SRC_DIR
$PYTHON setup.py install
