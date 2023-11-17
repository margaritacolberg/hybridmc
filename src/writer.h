// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause
//
// writer.h writes the positions of the beads, their velocities,
// configurations, distances between the beads, and the number of times each
// configuration is encountered to an output file

#ifndef HYBRIDMC_WRITER_H
#define HYBRIDMC_WRITER_H

#include "vec3.h"
#include <set>

class PosWriter {
  H5::DataSet dataset;
  unsigned int step;

public:
  PosWriter(H5::H5Location &file, const std::string &name,
            const unsigned int nbeads) {
    const hsize_t mem_dims[2] = {nbeads, 3};
    const hsize_t file_dims[3] = {0, nbeads, 3};
    const hsize_t max_dims[3] = {H5S_UNLIMITED, nbeads, 3};
    const hsize_t chunk_dims[3] = {10, nbeads, 3};
    H5::DataSpace mem_space(2, mem_dims);
    H5::DataSpace file_space(3, file_dims, max_dims);

    H5::DSetCreatPropList prop;
    prop.setChunk(3, chunk_dims);
    prop.setDeflate(6);

    dataset =
        file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, file_space, prop);

    step = 0;
  }

  void append(const std::vector<Vec3> &pos) {
    const unsigned int nbeads = pos.size();

    const hsize_t mem_dims[2] = {nbeads, 3};
    const hsize_t file_dims[3] = {step + 1, nbeads, 3};
    dataset.extend(file_dims);
    H5::DataSpace mem_space(2, mem_dims);
    H5::DataSpace file_space = dataset.getSpace();

    const hsize_t count[3] = {1, nbeads, 3};
    const hsize_t start[3] = {step, 0, 0};
    file_space.selectHyperslab(H5S_SELECT_SET, count, start);
    dataset.write(&pos[0], H5::PredType::NATIVE_DOUBLE, mem_space, file_space);

    step++;
  }
};

class VelWriter {
  H5::DataSet dataset;
  unsigned int step;

public:
  VelWriter(H5::H5Location &file, const std::string &name,
            const unsigned int nbeads) {
    const hsize_t mem_dims[2] = {nbeads, 3};
    const hsize_t file_dims[3] = {0, nbeads, 3};
    const hsize_t max_dims[3] = {H5S_UNLIMITED, nbeads, 3};
    const hsize_t chunk_dims[3] = {10, nbeads, 3};
    H5::DataSpace mem_space(2, mem_dims);
    H5::DataSpace file_space(3, file_dims, max_dims);

    H5::DSetCreatPropList prop;
    prop.setChunk(3, chunk_dims);
    prop.setDeflate(6);

    dataset =
        file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, file_space, prop);

    step = 0;
  }

  void append(const std::vector<Vec3> &vel) {
    const unsigned int nbeads = vel.size();

    const hsize_t mem_dims[2] = {nbeads, 3};
    const hsize_t file_dims[3] = {step + 1, nbeads, 3};
    dataset.extend(file_dims);
    H5::DataSpace mem_space(2, mem_dims);
    H5::DataSpace file_space = dataset.getSpace();

    const hsize_t count[3] = {1, nbeads, 3};
    const hsize_t start[3] = {step, 0, 0};
    file_space.selectHyperslab(H5S_SELECT_SET, count, start);
    dataset.write(&vel[0], H5::PredType::NATIVE_DOUBLE, mem_space, file_space);

    step++;
  }
};

class ConfigWriter {
  H5::DataSet dataset;
  unsigned int step;

public:
  ConfigWriter(H5::H5Location &file, const std::string &name) {
    const hsize_t mem_dims[1] = {1};
    const hsize_t file_dims[1] = {0};
    const hsize_t max_dims[1] = {H5S_UNLIMITED};
    const hsize_t chunk_dims[1] = {1000};
    H5::DataSpace mem_space(1, mem_dims);
    H5::DataSpace file_space(1, file_dims, max_dims);

    H5::DSetCreatPropList prop;
    prop.setChunk(1, chunk_dims);
    prop.setDeflate(6);

    dataset =
        file.createDataSet(name, H5::PredType::NATIVE_UINT64, file_space, prop);

    step = 0;
  }

  void append(const Config &config) {
    const hsize_t mem_dims[1] = {1};
    const hsize_t file_dims[1] = {step + 1};
    dataset.extend(file_dims);
    H5::DataSpace mem_space(1, mem_dims);
    H5::DataSpace file_space = dataset.getSpace();

    const hsize_t count[1] = {1};
    const hsize_t start[1] = {step};
    file_space.selectHyperslab(H5S_SELECT_SET, count, start);
    dataset.write(&config, H5::PredType::NATIVE_UINT64, mem_space, file_space);

    step++;
  }
};

class ConfigCountWriter {
  H5::DataSet dataset;
  unsigned int step;

public:
  ConfigCountWriter(H5::H5Location &file, const std::string &name,
                    const unsigned int nstates) {
    const hsize_t mem_dims[1] = {nstates};
    const hsize_t file_dims[2] = {0, nstates};
    const hsize_t max_dims[2] = {H5S_UNLIMITED, nstates};
    const hsize_t chunk_dims[2] = {100, nstates};
    H5::DataSpace mem_space(1, mem_dims);
    H5::DataSpace file_space(2, file_dims, max_dims);

    H5::DSetCreatPropList prop;
    prop.setChunk(2, chunk_dims);
    prop.setDeflate(6);

    dataset =
        file.createDataSet(name, H5::PredType::NATIVE_UINT64, file_space, prop);

    step = 0;
  }

  void append(const std::vector<uint64_t> &config_count) {
    const unsigned int nconfig = config_count.size();

    const hsize_t mem_dims[1] = {nconfig};
    const hsize_t file_dims[2] = {step + 1, nconfig};
    dataset.extend(file_dims);
    H5::DataSpace mem_space(1, mem_dims);
    H5::DataSpace file_space = dataset.getSpace();

    const hsize_t count[2] = {1, nconfig};
    const hsize_t start[2] = {step, 0};
    file_space.selectHyperslab(H5S_SELECT_SET, count, start);
    dataset.write(config_count.data(), H5::PredType::NATIVE_UINT64, mem_space,
                  file_space);

    step++;
  }
};

class DistWriter {
  H5::DataSet dataset;
  unsigned int step;

  public:
      DistWriter() = default;
      DistWriter(H5::H5Location &file, const std::string &name, const unsigned int nbonds) {

    const hsize_t mem_dims[1] = {nbonds};
    const hsize_t file_dims[2] = {0, nbonds};
    const hsize_t max_dims[2] = {H5S_UNLIMITED, nbonds};
    const hsize_t chunk_dims[2] = {10000, nbonds};
    H5::DataSpace mem_space(1, mem_dims);
    H5::DataSpace file_space(2, file_dims, max_dims);

    H5::DSetCreatPropList prop;
    prop.setChunk(2, chunk_dims);
    prop.setDeflate(6);

    // create a dataset called dist with double values indicating distance
    dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, file_space, prop);

    step = 0;
  }

  void append(const std::vector<double> &dist) {
    const unsigned int nbonds = dist.size();

    const hsize_t mem_dims[1] = {nbonds};
    const hsize_t file_dims[2] = {step + 1, nbonds};
    dataset.extend(file_dims);
    H5::DataSpace mem_space(1, mem_dims);
    H5::DataSpace file_space = dataset.getSpace();

    const hsize_t count[2] = {1, nbonds};
    const hsize_t start[2] = {step, 0};
    file_space.selectHyperslab(H5S_SELECT_SET, count, start);
    dataset.write(&dist[0], H5::PredType::NATIVE_DOUBLE, mem_space, file_space);

    step++;
  }

  int get_size() {
    H5::DataSpace dataspace = dataset.getSpace();
      hsize_t dims_out[2];
      dataspace.getSimpleExtentDims(dims_out, nullptr);
      return dims_out[0];
  }

};

#endif
