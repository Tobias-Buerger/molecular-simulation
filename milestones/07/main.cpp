#include "ducastelle.h"
#include "verlet.h"
#include "xyz.h"
#include "tools.h"
#include "thermostat.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <argparse/argparse.hpp>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("07");

    program.add_argument("--cutoff").default_value(7.0).scan<'g', double>();
    program.add_argument("--mass").default_value(196.97).scan<'g', double>();
    program.add_argument("--total_time").default_value(60000.0).scan<'g', double>();
    program.add_argument("--initial_tmp").default_value(100.0).scan<'g', double>();
    program.add_argument("--relaxation_time").default_value<size_t>(100).scan<'u', size_t>();
    program.add_argument("--delta_q").default_value(1.5).scan<'g', double>();
    program.add_argument("--output_interval").default_value<int>(200).scan<'i', int>();
    program.add_argument("input_file").default_value("cluster_923.xyz");
    program.add_argument("--traj_file").help("trajectory file").default_value("traj.xyz");
    program.add_argument("--csv_file").help("csv file for run statistics").default_value("data.csv");
    program.add_argument("-q", "--quiet").default_value(false).implicit_value(true);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    const double mass = program.get<double>("--mass") * 103.6;
    const double total_time = program.get<double>("--total_time");
    const double timestep = 1.0;
    const size_t timesteps = total_time / timestep;
    const double initial_tmp = program.get<double>("--initial_tmp") * 1e-5;
    const double cutoff = program.get<double>("--cutoff");
    const double delta_q = program.get<double>("--delta_q");
    const size_t relaxation_time = program.get<size_t>("--relaxation_time");
    const int output_interval = program.get<int>("--output_interval");
    const bool quiet = program.get<bool>("-q");
    const double nan = 0.0 / 0.0;

    fs::path filepath = argv[0];
    filepath = filepath.parent_path();
    filepath.append(program.get<std::string>("input_file"));

    auto [names, positions]{read_xyz(filepath)};
    Atoms atoms(names, positions);
    atoms.set_mass(mass);
    NeighborList neighborlist(cutoff);

    std::ofstream traj(program.get<std::string>("--traj_file"));
    CSVWriter csv(program.get<std::string>("--csv_file"));
    double tmp_sum = 0.0;
    double alpha = 0.01;

    if (!quiet) {
        std::cout << "setting initial temperature..." << std::endl;
    }
    for (size_t i = 0; i < 10000; i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces,
                     timestep, mass);
        neighborlist.update(atoms);
        double pot = ducastelle(atoms, neighborlist, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        double tmp = atoms.current_temperature();
        if (i < 1000) {
            berendsen_thermostat(atoms, initial_tmp, tmp, timestep, 100);
        }
        if (!quiet && (i % 100 == 0)) {
            std::cout << (double)i / 100.0 << "%" << std::endl;
        }
    }

    for (size_t i = 0; i < timesteps; i++) {
        if (i % output_interval == 0) {
            write_xyz(traj, atoms);
        }

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces,
                     timestep, mass);
        neighborlist.update(atoms);
        double pot = ducastelle(atoms, neighborlist, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        double kin = atoms.kinetic_energy();
        double tmp = atoms.current_temperature() * 1e5;

        if ((i + relaxation_time) % (relaxation_time * 2) == 0) {
            // add delta_q kinetic energy
            atoms.velocities *= std::sqrt(1.0 + delta_q / kin);
        }

        if (i % (relaxation_time * 2) < relaxation_time) {
            // measure average temperature: use exponential smoothing
            tmp_sum += tmp;
        }

        if ((i + relaxation_time) % output_interval == 0) {
            double avg_tmp = tmp_sum / relaxation_time;
            tmp_sum = 0.0;
            csv.write(i / output_interval, i * timestep, kin, pot, avg_tmp, nan, nan);
            if (!quiet) {
                std::cout << "frame: " << i / output_interval << " kin: " << kin << " pot: " << pot
                        << " tot: " << kin + pot << " tmp: " << avg_tmp << std::endl;
            }
        }
    }

    traj.close();

    return 0;
}
