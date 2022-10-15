# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.13.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(hybridmc LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS filesystem program_options unit_test_framework)
find_package(HDF5 REQUIRED COMPONENTS CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_COMMIT)
git_local_changes(GIT_CLEAN_OR_DIRTY)
set(VERSION "${GIT_COMMIT}-${GIT_CLEAN_OR_DIRTY}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

enable_testing()

add_subdirectory(src)
