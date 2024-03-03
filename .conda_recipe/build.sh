#!/bin/bash

# Change to the hmc directory and run make
cd hybridmc || exit
make release
cp release/hybridmc $PREFIX/bin/hybridmc.exe

# Enter back to top level directory and run setup.py
cd ../
$PYTHON setup.py install