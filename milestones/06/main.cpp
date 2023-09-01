#include "lj.h"
#include "thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include "tools.h"
#include <argparse/argparse.hpp>
#include <cmath>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    argparse::ArgumentParser program("06");

    program.add_argument("--mass").default_value(1.0).scan<'g', double>();
    program.add_argument("--total_time").default_value(100.0).scan<'g', double>();
    program.add_argument("--timestep").default_value(0.001).scan<'g', double>();
    program.add_argument("--relaxation_time").help("relaxation time as a multiple of timesteps").default_value(500.0).scan<'g', double>();
    program.add_argument("--target_temperature").help("target temperature for the thermostat in Kelvin").default_value(0.5).scan<'g', double>();
    program.add_argument("--lattice_distance").help("lattice distance as a multiple of sigma (sigma = 1)").default_value(1.08).scan<'g', double>();
    program.add_argument("--lattice").help("number of atoms in x,y,z dimension").nargs(3).default_value(std::vector<int>{3, 3, 3}).scan<'i', int>();
    program.add_argument("--alpha").help("temperature smoothing factor 1=no smoothing 0=max").default_value(0.1).scan<'g', double>();
    program.add_argument("--cutoff").help("cutoff for neighbor list").default_value(3.0).scan<'g', double>();
    program.add_argument("--output_interval").default_value<int>(100).scan<'i', int>();
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

    const double mass = program.get<double>("--mass");
    const double epsilon = 1.0;
    const double sigma = 1.0;
    const double total_time = program.get<double>("--total_time") * std::sqrt(mass * sigma * sigma / epsilon);
    const double timestep = program.get<double>("--timestep") * std::sqrt(mass * sigma * sigma / epsilon);
    const size_t timesteps = total_time / timestep;
    const double relaxation_time = program.get<double>("--relaxation_time") * timestep;
    const double target_temperature = program.get<double>("--target_temperature") * 1e-5;
    const double lattice_distance = program.get<double>("--lattice_distance") * sigma;
    const double alpha = program.get<double>("--alpha");
    const double cutoff = program.get<double>("--cutoff");
    const int output_interval = program.get<int>("--output_interval");
    const bool quiet = program.get<bool>("-q");
    const double nan = 0.0 / 0.0;
    const std::vector<int> lattice = program.get<std::vector<int>>("--lattice");

    Atoms atoms = create_cubic_lattice(lattice[0], lattice[1], lattice[2], lattice_distance);

    std::ofstream traj(program.get<std::string>("--traj_file"));
    CSVWriter csv(program.get<std::string>("--csv_file"));

    double avg_tmp = 0.0;
    for (size_t i = 0; i < timesteps; i++) {
        if (i % output_interval == 0) {
            write_xyz(traj, atoms);
        }

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces,
                     timestep, mass);
        double pot = lj_cutoff(atoms, cutoff, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        double kin = atoms.kinetic_energy();
        double current_tmp = atoms.current_temperature();

        berendsen_thermostat(atoms, target_temperature, current_tmp, timestep,
                             relaxation_time);

        avg_tmp = alpha * (current_tmp * 1e5) + (1 - alpha) * avg_tmp;

        if (i % output_interval == 0) {
            csv.write(i / output_interval, i * timestep, kin, pot, avg_tmp, nan, nan);
            if (!quiet) {
                std::cout << "kin: " << kin << " pot: " << pot << " tot: " << kin + pot << " tmp: " << avg_tmp << std::endl;
            }
        }
    }

    traj.close();

    return 0;
}
