# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

# cmake version needed
cmake_minimum_required(VERSION 3.26.4)

# set params for cmale
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# define project
project(hardspheres LANGUAGES CXX)

# find needed packages
find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)
find_package(HDF5 REQUIRED COMPONENTS CXX)

# add required directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# setup git
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_COMMIT)
git_local_changes(GIT_CLEAN_OR_DIRTY)
set(VERSION "${GIT_COMMIT}-${GIT_CLEAN_OR_DIRTY}")

# define output
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# allow testing
enable_testing()

# go into src cmake text
add_subdirectory(src)
