#include "mpi_support.h"
#include "domain.h"
#include "ducastelle.h"
#include "verlet.h"
#include "xyz.h"
#include "tools.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <argparse/argparse.hpp>

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    auto mpi_guard = MPI::init_guard(&argc, &argv);
    
    argparse::ArgumentParser program("08");

    program.add_argument("--cutoff").default_value(10.0).scan<'g', double>();
    program.add_argument("--mass").default_value(196.97).scan<'g', double>();
    program.add_argument("--total_time").default_value(10000.0).scan<'g', double>();
    program.add_argument("--extra_space").help("additional space for the domain border").default_value(4.0).scan<'g', double>();
    program.add_argument("--domains").help("number of domains in x,y,z dimension").nargs(3).default_value(std::vector<int>{1, 2, 2}).scan<'i', int>();
    program.add_argument("--output_interval").default_value<int>(100).scan<'i', int>();
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

    double mass = program.get<double>("--mass") * 103.6;
    double total_time = program.get<double>("--total_time");
    double extra_space = program.get<double>("--extra_space");
    double timestep = 1.0;
    size_t timesteps = total_time / timestep;
    double cutoff = program.get<double>("--cutoff");
    std::vector<int> domains = program.get<std::vector<int>>("--domains");
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

    std::ofstream traj;
    CSVWriter csv;

    if (MPI::comm_rank(MPI_COMM_WORLD) == 0) {
        traj.open(program.get<std::string>("--traj_file"));
        csv.open(program.get<std::string>("--csv_file"));
    }

    Eigen::Array3d min_box = atoms.positions.rowwise().minCoeff();
    Eigen::Array3d max_box = atoms.positions.rowwise().maxCoeff();
    double shiftx, shifty, shiftz;
    double boxx, boxy, boxz;
    shiftx = -(double)min_box(0) + extra_space;
    shifty = -(double)min_box(1) + extra_space;
    shiftz = -(double)min_box(2) + extra_space;
    auto shift = Eigen::Array3d{shiftx, shifty, shiftz};
    atoms.positions.colwise() += shift;
    boxx = (double)max_box(0) - (double)min_box(0) + extra_space * 2;
    boxy = (double)max_box(1) - (double)min_box(1) + extra_space * 2;
    boxz = (double)max_box(2) - (double)min_box(2) + extra_space * 2;

    Domain domain(MPI_COMM_WORLD,
        {boxx, boxy, boxz},
        {domains[0], domains[1], domains[2]},
        {0, 0, 0});
    domain.enable(atoms);
    for (size_t i = 0; i < timesteps; i++) {
        if (i % output_interval == 0) {
            domain.disable(atoms);
            if (MPI::comm_rank(MPI_COMM_WORLD) == 0) {
                write_xyz(traj, atoms);
            }
            domain.enable(atoms);
        }

        verlet_step1(atoms.positions, atoms.velocities, atoms.forces,
                     timestep, mass);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, cutoff * 2.0);
        neighborlist.update(atoms);
        double local_pot = ducastelle(atoms, neighborlist, domain.nb_local(), cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
        double local_kin = atoms.kinetic_energy(domain.nb_local());

        double pot = MPI::allreduce(local_pot, MPI_SUM, MPI_COMM_WORLD);
        double kin = MPI::allreduce(local_kin, MPI_SUM, MPI_COMM_WORLD);

        if (i % output_interval == 0) {
            csv.write(i / output_interval, i * timestep, kin, pot, nan, nan, nan);
            domain.disable(atoms);
            if (MPI::comm_rank(MPI_COMM_WORLD) == 0) {
                if (!quiet) {
                    std::cout << "kin: " << kin << " pot: " << pot
                            << " tot: " << kin + pot << std::endl;
                }
            }
            domain.enable(atoms);
        }
    }

    traj.close();

    return 0;
}
