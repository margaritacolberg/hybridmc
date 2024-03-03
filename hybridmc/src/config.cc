// Copyright (c) 2018-2022 Margarita Colberg
// SPDX-License-Identifier: BSD-3-Clause

#include "config.h"

// TODO: add the k index for the rc value
NonlocalBonds::NonlocalBonds(const Pairs &ij) {
    for (auto const &[i, j, k]: ij) {

// as discussed on slack #hybridmc
// type check
        if (i > j) {
            ij_.emplace_back(std::make_tuple(j, i, k * k));
        }

        else {
            ij_.emplace_back(std::make_tuple(i, j, k * k));
        }
    }
}

// obtain the number of bonds in the given list of bonds
unsigned int NonlocalBonds::get_nbonds() const {return ij_.size();}

// obtain the bond_index'th bond in the master bonds list
std::tuple<unsigned int, unsigned int, double> NonlocalBonds::getBond(unsigned int bond_index) const {
    //std::cout << " Bond pair " << std::get<0>(ij_[bond_index]) << std::endl;
    return ij_[bond_index];
}

double NonlocalBonds::getrc(unsigned int bond_index) const {
    return std::sqrt(std::get<2>(ij_[bond_index]));
}

std::tuple<Config, double> NonlocalBonds::get_bond_mask(unsigned int i, unsigned int j) const {
    assert(j > i);

    //std::tuple<int,int,double> searchObject(i,j, 1.5);
    // check that i and j match with i and j from the json file
    const auto it = std::find_if(ij_.begin(), ij_.end(),
                                 [i,j](auto& e) {
                                     return (std::get<0>(e) == i and std::get<1>(e) == j);
                                 });

    // 0 if no bonds formed. pair i and j not in one of nonlcoal bonds lists
    if (it == ij_.end())
        return std::make_tuple(0, 0);

    double rc2 = std::get<2>(*it);

    // check that the number of configurations do not exceed 64 bits
    assert(ij_.size() <= std::numeric_limits<Config>::digits);

    // convert 1 from 32 to 64 bits and determine position of new bond
    // in bond pattern (using left shift operator)

    // make first element of return tuple the config and second the rc2
    Config returnConfig = Config(1) << (it - ij_.begin());
    return std::make_tuple(returnConfig, rc2);
}

// print out bonds list
void NonlocalBonds::printBonds()
{
    for (auto & i : ij_){
        std::cout << " Bond pair " << std::get<0>(i) << "-" << std::get<1>(i)
                  << " has rc2 =  " << std::get<2>(i) << std::endl;
    }
}

// write the transient or permanent bead indices to the output file
void NonlocalBonds::write_hdf5(H5::H5Object &obj,
                               const std::string &name) const {

  int nb = ij_.size();

  BondWrapper *bw = new BondWrapper[nb];
  for (int k=0;k<nb;k++)
  {
        bw[k].i = std::get<0>(ij_[k]);
        bw[k].j = std::get<1>(ij_[k]);
        bw[k].rc = std::sqrt(std::get<2>(ij_[k]));
  }

  H5::CompType mtype1( sizeof(BondWrapper) );
  mtype1.insertMember("bond i", HOFFSET(BondWrapper, i), H5::PredType::NATIVE_UINT32);
  mtype1.insertMember("bond j", HOFFSET(BondWrapper, j), H5::PredType::NATIVE_UINT32);
  mtype1.insertMember("rc", HOFFSET(BondWrapper, rc), H5::PredType::NATIVE_DOUBLE);

  const hsize_t dims[1] = {ij_.size()};
  H5::DataSpace space(1,dims);

  H5::Attribute attribute= obj.createAttribute(name,mtype1,space);
  if (!ij_.empty())
  {
      attribute.write(mtype1,&bw[0]);
  }

  delete []bw;
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
