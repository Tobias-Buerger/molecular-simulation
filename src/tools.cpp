#include "tools.h"

#include <string>
#include <fstream>

CSVWriter::CSVWriter() : file{std::ofstream()} {}

CSVWriter::CSVWriter(const std::string& filename) : file{std::ofstream()} {
    open(filename);
}

CSVWriter::~CSVWriter() {
    file.close();
}

void CSVWriter::open(const std::string& filename) {
    file.open(filename);
    file << "frame,simulated_time,kinetic_energy,potential_energy,total_energy,temperature,stress,strain\n";
    file.flush();
}

void CSVWriter::write(int frame, double time, double kin, double pot,
                      double tmp, double stress, double strain) {
    file << frame << ","
         << time << ","
         << kin << ","
         << pot << ","
         << kin + pot << ","
         << tmp << ","
         << stress << ","
         << strain << "\n";
    file.flush();
}