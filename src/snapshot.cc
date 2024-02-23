// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// snapshot.cc saves the most recent positions and entropies of the hard
// spheres and state of random number generator in the event that the
// simulation is terminated before it converges

#include "snapshot.h"
#include <H5Cpp.h>
#include <cerrno>
#include <cstring>
#include <unistd.h>

void write_ensemble(H5::H5Location &file, const System &sys) {

  std::cout << "  Called new write_ensemble with ensemble size "
        << sys.nextEnsemble.size() << std::endl;
  // Write ensemble of structures

  for (size_t i = 0; i < sys.nextEnsemble.size(); ++i) {
        const std::string dataset_name = "molecule_" + std::to_string(i);
        auto& molecule = sys.nextEnsemble[i];
        auto positions = molecule.getPositions();
        hsize_t dims[2] = {positions.size(), 3};

        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = file.createDataSet(dataset_name,
                                                  H5::PredType::NATIVE_DOUBLE, dataspace);

        // Allocate memory for all positions
        std::vector<double> allPositions(positions.size() * positions[0].size());
        size_t posIndex = 0;
        for (const auto& position : positions) {
            for (const auto& coord : position) {
                allPositions[posIndex++] = coord;
            }
        }

        dataset.write(allPositions.data(), H5::PredType::NATIVE_DOUBLE);
  }
}


void write_snapshot(const std::string snapshot_name,
                    const std::vector<Vec3> &pos,
                    const std::vector<double> &s_bias,
                    Random &mt,
                    UpdateConfig &update_config) {
  H5::H5File file(snapshot_name, H5F_ACC_TRUNC);

  {
    const hsize_t mem_dims[2] = {pos.size(), 3};
    H5::DataSpace mem_space(2, mem_dims);
    H5::DataSet dataset_pos =
        file.createDataSet("pos", H5::PredType::NATIVE_DOUBLE, mem_space);

    dataset_pos.write(&pos[0], H5::PredType::NATIVE_DOUBLE, mem_space);
  }

  {
    const hsize_t mem_dims[2] = {s_bias.size(), 1};
    H5::DataSpace mem_space(2, mem_dims);
    H5::DataSet dataset_s_bias =
        file.createDataSet("s_bias", H5::PredType::NATIVE_DOUBLE, mem_space);

    dataset_s_bias.write(&s_bias[0], H5::PredType::NATIVE_DOUBLE, mem_space);
  }

  // save most recent state of random number generator as string to be able
  // to resume from this exact state of the program later (rather than using
  // new seed to resume, which generates a state different than the
  // interrupted state)
  {
    std::stringstream ss;
    ss << mt;

    H5::StrType datatype_mt(0, H5T_VARIABLE);
    H5::DataSet dataset_mt = file.createDataSet("mt", datatype_mt, H5S_SCALAR);

    dataset_mt.write(ss.str(), datatype_mt);
  }

  // flush written data from internal hdf5 buffer to os
  file.flush(H5F_SCOPE_GLOBAL);

  // retrieve os file descriptor for hdf5 file
  void *fd;
  file.getVFDHandle(&fd);

  // flush written data from os buffer to disk
  if (0 != fsync(*reinterpret_cast<int *>(fd))) {
    throw std::runtime_error(std::strerror(errno));
  }

  // output most recent configuration
  H5::Attribute attr_config =
      file.createAttribute("config", H5::PredType::NATIVE_UINT64, H5S_SCALAR);
  attr_config.write(H5::PredType::NATIVE_UINT64, &update_config.config);
}

void read_snapshot(const std::string snapshot_name, std::vector<Vec3> &pos,
                   std::vector<double> &s_bias, Random &mt,
                   UpdateConfig &update_config) {
  H5::H5File file(snapshot_name, H5F_ACC_RDONLY);

  {
    const H5::DataSet dataset_pos = file.openDataSet("pos");
    const H5::DataSpace filespace = dataset_pos.getSpace();

    hsize_t filedims[2];
    filespace.getSimpleExtentDims(filedims);
    pos.resize(filedims[0]);

    const hsize_t mem_dims[2] = {pos.size(), 3};
    H5::DataSpace mem_space(2, mem_dims);

    dataset_pos.read(pos.data(), H5::PredType::NATIVE_DOUBLE, mem_space);
  }

  {
    const H5::DataSet dataset_s_bias = file.openDataSet("s_bias");
    const H5::DataSpace filespace = dataset_s_bias.getSpace();

    hsize_t filedims[2];
    filespace.getSimpleExtentDims(filedims);
    s_bias.resize(filedims[0]);

    const hsize_t mem_dims[2] = {s_bias.size(), 1};
    H5::DataSpace mem_space(2, mem_dims);

    dataset_s_bias.read(s_bias.data(), H5::PredType::NATIVE_DOUBLE, mem_space);
  }

  {
    const H5::DataSet dataset_mt = file.openDataSet("mt");
    const H5::StrType datatype_mt(0, H5T_VARIABLE);

    std::string state;
    dataset_mt.read(state, datatype_mt, H5S_SCALAR);

    std::stringstream ss;
    ss.str(state);

    ss.exceptions(std::stringstream::failbit);

    ss >> mt;
  }

  H5::Attribute attr_config = file.openAttribute("config");
  attr_config.read(H5::PredType::NATIVE_UINT64, &update_config.config);
}

void read_input(const std::string input_name, std::vector<Vec3> &pos) {
  H5::H5File file(input_name, H5F_ACC_RDONLY);
  const H5::Group group = file.openGroup("unique_config");

  const H5::DataSet dataset_pos = group.openDataSet("pos");
  const hsize_t mem_dims[2] = {pos.size(), 3};
  H5::DataSpace mem_space(2, mem_dims);

  dataset_pos.read(pos.data(), H5::PredType::NATIVE_DOUBLE, mem_space);
}
