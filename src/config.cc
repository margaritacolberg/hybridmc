// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#include "config.h"


// TODO: add the k index for the rc value
NonlocalBonds::NonlocalBonds(const Pairs &ij) {
  for (auto const& [i, j, k]: ij) {


// as discussed on slack #hybridmc
// type check
    if (i > j) {
        ij_.emplace_back(std::make_tuple(j, i, k * k) );
    } else {
        ij_.emplace_back(std::make_tuple(i, j, k * k) );
    }
}


unsigned int NonlocalBonds::get_nbonds() const { return ij_.size(); }

// write the transient or permanent bead indices to the output file
void NonlocalBonds::write_hdf5(H5::H5Object &obj,
                               const std::string &name) const {
  const hsize_t dims[2] = {ij_.size(), 2};
  H5::DataSpace space(2, dims);
  H5::Attribute attr =
      obj.createAttribute(name, H5::PredType::NATIVE_UINT, space);
  if (!ij_.empty()) {
    attr.write(H5::PredType::NATIVE_UINT, &ij_[0]);
  }
}

void from_json(const nlohmann::json &json, NonlocalBonds &nonlocal_bonds) {
  nonlocal_bonds = json.get<NonlocalBonds::Pairs>();
}

// write the configuration of the protein as an int to the output file
UpdateConfigWriter::UpdateConfigWriter(H5::H5File &file) {
  file_dims[0] = 0;
  const hsize_t max_dims[1] = {H5S_UNLIMITED};
  file_space.setExtentSimple(1, file_dims, max_dims);

  H5::DSetCreatPropList dcpl;
  const hsize_t chunk_dims[1] = {1024};
  dcpl.setChunk(1, chunk_dims);

  H5::Group group = file.createGroup("config");
  dataset_int =
      group.createDataSet("int", H5::PredType::NATIVE_UINT64, file_space, dcpl);
}

void UpdateConfigWriter::append() {
  const hsize_t count[1] = {config_int.size()};
  H5::DataSpace mem_space(1, count);

  const hsize_t start[1] = {file_dims[0]};
  file_dims[0] += count[0];
  file_space.setExtentSimple(1, file_dims);
  file_space.selectHyperslab(H5S_SELECT_SET, count, start);

  dataset_int.extend(file_dims);
  dataset_int.write(&config_int[0], H5::PredType::NATIVE_UINT64, mem_space,
                    file_space);
}
