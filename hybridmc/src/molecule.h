#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <H5Cpp.h>
#include "vec3.h"

#ifndef MOLECULE_H_
#define MOLECULE_H_

class Molecule {
private:
    int n;  // number of atoms
    std::vector<std::vector<double>> positions;

public:
    Molecule() {}
    Molecule(int numAtoms) : n(numAtoms) {
        //std::cout << " Molecule constructor called with numAtoms = " << numAtoms << std::endl;
        positions.resize(n);
    }

    ~Molecule() { n=0; positions.clear(); }

    void setPositions(const std::vector<std::vector<double>>& pos) {
        positions = pos;
    }

    void setPositions(const std::vector<Vec3> &pos)
    {

        n = pos.size();
        positions.resize(n);

        //std::cout << " Molecule setPositions for Vec3 called with length of " << n << std::endl;
        for (int i=0;i<n;i++)
        {
            positions[i].clear();
            positions[i].push_back(pos[i].x);
            positions[i].push_back(pos[i].y);
            positions[i].push_back(pos[i].z);
        }
    }

    void initializePositions(std::mt19937 &gen) {

        std::uniform_real_distribution<double> dis(0.0, 1.0);

        for (int i = 0; i < n; ++i) {
            std::vector<double> position(3);
            for (int j = 0; j < 3; ++j) {
                //double randomPosition = dis(gen); // Generate random number between 0 and 1
                position[j] = dis(gen);
            }
            positions[i] = position;
        }
    }


    const std::vector<std::vector<double>>& getPositions() const {
        return positions;
    }

    std::vector<Vec3> getVec3Positions() const {
        std::vector<Vec3> pos(n);
        for (int i=0;i<n;i++){
            pos[i].x = positions[i][0];
            pos[i].y = positions[i][1];
            pos[i].z = positions[i][2];
        }
        return pos;
    }

    void printPositions(int i=0) const {
        std::cout << std::endl << " Molecule " << i << " has positions:" << std::endl;
        for (const auto& pos : positions) {
            for (const auto& coord : pos) {
                std::cout << coord << " ";
            }
            std::cout << std::endl;
        }
    }
};

void writeMoleculesToHDF5(const char* filename, const std::vector<Molecule>& molecules);
std::vector<Molecule> readMoleculesFromHDF5(hid_t &file_id);
std::vector<Molecule> readMoleculesFromHDF5ByName(const char *filename);
void printEnsemble(std::vector<Molecule> &ensemble);
#endif