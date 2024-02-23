#include "molecule.h"

void writeMoleculesToHDF5(const char* filename, const std::vector<Molecule>& molecules) {
    hid_t file_id, group_id, dataspace_id, dataset_id;
    hsize_t dims[2];
    const char* group_name = "/ensemble";

    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Error creating HDF5 file." << std::endl;
        return;
    }

    group_id = H5Gcreate2(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) {
        std::cerr << "Error creating group for molecules." << std::endl;
        H5Fclose(file_id);
        return;
    }

    for (size_t i = 0; i < molecules.size(); ++i) {
        const char* dataset_name = ("molecule_" + std::to_string(i)).c_str();
        auto& molecule = molecules[i];
        auto positions = molecule.getPositions();
        dims[0] = positions.size();
        dims[1] = 3;

        dataspace_id = H5Screate_simple(2, dims, NULL);
        if (dataspace_id < 0) {
            std::cerr << "Error creating dataspace for molecule " << i << std::endl;
            H5Gclose(group_id);
            H5Fclose(file_id);
            return;
        }

        dataset_id = H5Dcreate2(group_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dataset_id < 0) {
            std::cerr << "Error creating dataset for molecule " << i << std::endl;
            H5Sclose(dataspace_id);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return;
        }

        // Allocate memory for all positions
        std::vector<double> allPositions(positions.size() * positions[0].size());
        size_t posIndex = 0;
        for (const auto& position : positions) {
            for (const auto& coord : position) {
                allPositions[posIndex++] = coord;
            }
        }

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, allPositions.data());

        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }

    H5Gclose(group_id);
    H5Fclose(file_id);
}


std::vector<Molecule> readMoleculesFromHDF5(hid_t &file_id) {
    std::vector<Molecule> molecules;
    hid_t group_id = H5Gopen2(file_id, "ensemble", H5P_DEFAULT);

    hsize_t num_objs;
    H5Gget_num_objs(group_id, &num_objs);
    for (hsize_t i = 0; i < num_objs; ++i) {
        std::string dataset_name = "molecule_" + std::to_string(i);
        hid_t dataset_id = H5Dopen2(group_id, dataset_name.c_str(), H5P_DEFAULT);

        hid_t dataspace_id = H5Dget_space(dataset_id);

        hsize_t dims[3];
        H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

         // Allocate memory for all positions
        std::vector<double> allPositions( 3 * dims[0]);


        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, allPositions.data());

        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);

        Molecule molecule(dims[0]);

        std::vector<std::vector<double> > positions(dims[0]);
        size_t pos_index = 0;
        for (int i = 0; i< (int)dims[0];i++)
        {
            for (int j=0;j<3;j++) positions[i].push_back( allPositions[pos_index++] ) ;
        }
        molecule.setPositions(positions);
        molecules.push_back(molecule);
    }

    H5Gclose(group_id);
    H5Fclose(file_id);

    return molecules;
}


std::vector<Molecule> readMoleculesFromHDF5ByName(const char *filename)
{
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    return readMoleculesFromHDF5(file_id);
}

void printEnsemble(std::vector<Molecule> &ens)
{
    int n = ens.size();
    std::cout << "  Printing out ensemble of size " << n << std::endl;
    for (int i=0;i<n;i++){
        ens[i].printPositions(i);
    }
}
