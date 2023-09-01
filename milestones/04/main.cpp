#include "xyz.h"
#include "verlet.h"
#include "lj.h"
#include "tools.h"
#include <argparse/argparse.hpp>
#include <filesystem>
#include <iostream>
#include <cmath>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("04");

    program.add_argument("--mass").default_value(1.0).scan<'g', double>();
    program.add_argument("--total_time").default_value(100.0).scan<'g', double>();
    program.add_argument("--timestep").default_value(0.001).scan<'g', double>();
    program.add_argument("--output_interval").default_value<int>(100).scan<'i', int>();
    program.add_argument("--traj_file").help("trajectory file").default_value("traj.xyz");
    program.add_argument("--csv_file").help("csv file for run statistics").default_value("data.csv");
    program.add_argument("input_file").default_value("lj54.xyz");
    program.add_argument("-q", "--quiet").default_value(false).implicit_value(true);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    const double mass = program.get<double>("--mass");
    const double epsilon = 1.0;
    const double sigma = 1.0;
    const double total_time = program.get<double>("--total_time") * std::sqrt(mass * sigma * sigma / epsilon);
    const double timestep = program.get<double>("--timestep") * std::sqrt(mass * sigma * sigma / epsilon);
    size_t timesteps = total_time / timestep;
    const int output_interval = program.get<int>("--output_interval");
    const bool quiet = program.get<bool>("-q");
    const double nan = 0.0 / 0.0;

    fs::path filepath = argv[0];
    filepath = filepath.parent_path();
    filepath.append(program.get<std::string>("input_file"));

    auto [names, positions, velocities]{read_xyz_with_velocities(filepath)};
    Atoms atoms(names, positions, velocities);

    std::ofstream traj(program.get<std::string>("--traj_file"));
    CSVWriter csv(program.get<std::string>("--csv_file"));

    for (size_t i = 0; i < timesteps; i++) {
        if (i % output_interval == 0) {
            write_xyz(traj, atoms);
        }

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
        double pot = lj_all(atoms, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        double kin = atoms.kinetic_energy();

        if (i % output_interval == 0) {
            csv.write(i / output_interval, i * timestep, kin, pot, nan, nan, nan);
            if (!quiet) {
                std::cout << "kin: " << kin << " pot: " << pot << " tot: " << kin + pot << std::endl;
            }
        }
    }

    traj.close();

    return 0;
}
