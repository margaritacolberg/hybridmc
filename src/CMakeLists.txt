# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

#-----------------------------------------------------------------------------------------------------

# create common library with all dependant source files for main
add_library(hybridmc_common OBJECT hardspheres.cc crankshaft.cc entropy.cc config.cc snapshot.cc)

target_include_directories(hybridmc_common SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

# Assign common libraries to hybridmc_common
target_link_libraries(hybridmc_common PUBLIC ${Boost_LIBRARIES} ${HDF5_LIBRARIES})
#------------------------------------------------------------------------------------------------------

# create executable called hybridmc that will run main.cc
add_executable(hybridmc main.cc)
# link all the dependancies (common libraries -- hybridmc_common) to hybridmc
target_link_libraries(hybridmc PRIVATE hybridmc_common)

#--------------------------------------------------------------
#pybind11 method for making wang landau port:
pybind11_add_module(wang_landau ./main.cc)
target_link_libraries(wang_landau PRIVATE hybridmc_common)
target_compile_definitions(wang_landau
        PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
#------------------------------------------------------------
#stuff for testing
#add_executable(hybridmc_test hybridmc_test.cc)
#target_link_libraries(hybridmc_test PRIVATE hybridmc_common)
#add_test(NAME hybridmc_test COMMAND hybridmc_test)
#-------------------------------------------------------------

set_property(SOURCE main.cc APPEND PROPERTY COMPILE_DEFINITIONS VERSION="${VERSION}")
