# Copyright (c) 2018-2022 Margarita Colberg
# SPDX-License-Identifier: BSD-3-Clause

add_library(hybridmc_common OBJECT hardspheres.cc crankshaft.cc entropy.cc config.cc snapshot.cc)
# include library headers as SYSTEM headers to suppress compiler warnings
# https://gcc.gnu.org/onlinedocs/cpp/System-Headers.html
target_include_directories(hybridmc_common SYSTEM PUBLIC ${Boost_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(hybridmc_common PUBLIC ${Boost_LIBRARIES} ${HDF5_LIBRARIES})

add_executable(hybridmc main.cc)
target_link_libraries(hybridmc PRIVATE hybridmc_common)

#add_executable(hybridmc_test hybridmc_test.cc)
#target_link_libraries(hybridmc_test PRIVATE hybridmc_common)
#add_test(NAME hybridmc_test COMMAND hybridmc_test)

set_property(SOURCE main.cc APPEND PROPERTY COMPILE_DEFINITIONS VERSION="${VERSION}")
