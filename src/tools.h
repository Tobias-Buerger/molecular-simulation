#ifndef __TOOLS_H
#define __TOOLS_H

#include "atoms.h"
#include <string>
#include <fstream>

class CSVWriter {
public:
    CSVWriter();
    CSVWriter(const std::string& filename);
    ~CSVWriter();

    /*
    Open the csv file.
    */
    void open(const std::string& filename);

    /*
    Write one row to the csv file.
    */
    void write(int frame, double time, double kin, double pot, double tmp, double stress, double strain);
private:
    std::ofstream file;
};


/*
Create a cubic lattice.
*/
inline Atoms create_cubic_lattice(int size_x, int size_y, int size_z, double lattice_distance) {
    Atoms atoms(size_x * size_y * size_z);
    for (size_t x = 0; x < size_x; x++) {
        for (size_t y = 0; y < size_y; y++) {
            for (size_t z = 0; z < size_z; z++) {
                Eigen::Vector3d pos = {(double) x, (double) y, (double) z};
                atoms.positions.col(x * size_y * size_z + y * size_z + z) = pos * lattice_distance;
            }
        }
    }
    return atoms;
}

#endif  // __TOOLS_H