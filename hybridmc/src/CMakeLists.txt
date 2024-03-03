# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

#-----------------------------------------------------------------------------------------------------

# create common library with all dependant source files for main
add_library(hybridmc_common OBJECT molecule.cc hardspheres.cc crankshaft.cc entropy.cc config.cc snapshot.cc main_helpers.cc)

target_include_directories(hybridmc_common SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

# Assign common libraries to hybridmc_common
target_link_libraries(hybridmc_common PUBLIC ${Boost_LIBRARIES} ${HDF5_LIBRARIES})
#------------------------------------------------------------------------------------------------------

# create executable called hybridmc that will run main.cc
add_executable(hybridmc main.cc)
# link all the dependancies (common libraries -- hybridmc_common) to hybridmc
target_link_libraries(hybridmc PRIVATE hybridmc_common)

#stuff for testing
#add_executable(hybridmc_test hybridmc_test.cc hybridmc_test.cc)
#target_link_libraries(hybridmc_test PRIVATE hybridmc_common)
#add_test(NAME hybridmc_test COMMAND hybridmc_test)
#-------------------------------------------------------------

set_property(SOURCE main.cc APPEND PROPERTY COMPILE_DEFINITIONS VERSION="${VERSION}")
