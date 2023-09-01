#include "mpi_support.h"
#include "domain.h"
#include "ducastelle.h"
#include "verlet.h"
#include "xyz.h"
#include "thermostat.h"
#include "tools.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <argparse/argparse.hpp>

namespace fs = std::filesystem;

double calc_local_stress(const Domain& domain, const Atoms& atoms, int dim) {
    double left_domain_boundary{domain.coordinate(dim) * domain.domain_length(dim) /
                                domain.decomposition(dim)};
    auto left_mask{atoms.positions.row(dim) < left_domain_boundary};
    double stress = 0.0;
    for (size_t i = 0; i < atoms.forces.cols(); i++) {
        if(left_mask[i]) {
            stress += atoms.forces(dim, i);
        }
    }
    return stress;
}

int main(int argc, char *argv[]) {
    auto mpi_guard = MPI::init_guard(&argc, &argv);
    
    argparse::ArgumentParser program("09");

    program.add_argument("--cutoff").default_value(7.0).scan<'g', double>();
    program.add_argument("--mass").default_value(196.97).scan<'g', double>();
    program.add_argument("--total_strain").default_value(25.0).scan<'g', double>();
    program.add_argument("--temperature").help("target temperature in Kelvin").default_value(0.0).scan<'g', double>();
    program.add_argument("--relaxation_time").default_value(500.0).scan<'g', double>();
    program.add_argument("--length_increase").help("how much more length to add every 100 timesteps on z-axis").default_value(0.05).scan<'g', double>();
    program.add_argument("--extra_space").help("additional space for the domain border").default_value(4.0).scan<'g', double>();
    program.add_argument("--domains").help("number of domains in x,y,z dimension").nargs(3).default_value(std::vector<int>{1, 1, 8}).scan<'i', int>();
    program.add_argument("--output_interval").default_value<int>(100).scan<'i', int>();
    program.add_argument("input_file").default_value("whisker_small.xyz");
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
    const double total_strain = program.get<double>("--total_strain");
    const double extra_space = program.get<double>("--extra_space");
    const double timestep = 1.0;
    const double cutoff = program.get<double>("--cutoff");
    const double length_increase = program.get<double>("--length_increase");
    const double target_temperature = program.get<double>("--temperature") * 1e-5;
    const double relaxation_time = program.get<double>("--relaxation_time");
    const size_t timesteps =  total_strain / length_increase * 100;
    const int output_interval = program.get<int>("--output_interval");
    const bool quiet = program.get<bool>("-q");
    const double nan = 0.0 / 0.0;
    std::vector<int> domains = program.get<std::vector<int>>("--domains");

    fs::path filepath = argv[0];
    filepath = filepath.parent_path();
    filepath.append(program.get<std::string>("input_file"));

    auto [names, positions]{read_xyz(filepath)};
    Atoms atoms(names, positions);
    atoms.set_mass(mass);
    NeighborList neighborlist(cutoff);

    std::cout << "atoms loaded" << std::endl;

    std::ofstream traj;
    CSVWriter csv;
    if (MPI::comm_rank(MPI_COMM_WORLD) == 0) {
        traj.open(program.get<std::string>("--traj_file"));
        csv.open(program.get<std::string>("--csv_file"));
    }

    Eigen::Array3d min_box = atoms.positions.rowwise().minCoeff();
    Eigen::Array3d max_box = atoms.positions.rowwise().maxCoeff();
    double shiftx, shifty, shiftz;
    double boxx = (double)max_box(0) - (double)min_box(0);
    double boxy = (double)max_box(1) - (double)min_box(1);
    double boxz = (double)max_box(2) - (double)min_box(2);
    shiftx = -(double)min_box(0) + extra_space;
    shifty = -(double)min_box(1) + extra_space;
    shiftz = 0.0;
    auto shift = Eigen::Array3d{shiftx, shifty, shiftz};
    atoms.positions.colwise() += shift;
    double domainx = boxx + extra_space * 2;
    double domainy = boxy + extra_space * 2;
    double domainz = (double)max_box(2) + (double)min_box(2);

    Domain domain(MPI_COMM_WORLD,
        {domainx, domainy, domainz},
        {domains[0], domains[1], domains[2]},
        {0, 0, 1});
    domain.enable(atoms);

    if (!quiet && MPI::comm_rank(MPI_COMM_WORLD) == 0) {
        std::cout << "setting initial temperature..." << std::endl;
    }
    for (size_t i = 0; i < 10000; i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces,
                     timestep, mass);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, cutoff * 2.0);
        neighborlist.update(atoms);
        double local_pot = ducastelle(atoms, neighborlist, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        double local_tmp = atoms.current_temperature(domain.nb_local());
        double tmp = MPI::allreduce(local_tmp, MPI_SUM, MPI_COMM_WORLD) / domain.size();

        if (i < 1000) {
            berendsen_thermostat(atoms, target_temperature, tmp, timestep, 100);
        }
        if (!quiet && (i % 100 == 0) && MPI::comm_rank(MPI_COMM_WORLD) == 0) {
            std::cout << (double)i / 100.0 << "%" << std::endl;
        }
    }

    double avg_stress = 0.0 / 0.0;

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
        double local_stress = calc_local_stress(domain, atoms, 2);
        verlet_step2(atoms.velocities, atoms.forces, timestep, mass);

        double local_kin = atoms.kinetic_energy(domain.nb_local());
        double local_tmp = atoms.current_temperature(domain.nb_local()) * 1e5;

        double pot = MPI::allreduce(local_pot, MPI_SUM, MPI_COMM_WORLD);
        double kin = MPI::allreduce(local_kin, MPI_SUM, MPI_COMM_WORLD);
        double tmp = MPI::allreduce(local_tmp, MPI_SUM, MPI_COMM_WORLD) / domain.size();
        double stress = MPI::allreduce(local_stress, MPI_SUM, MPI_COMM_WORLD) / (boxx * boxy * domains[2]);
        avg_stress += stress;

        berendsen_thermostat(atoms, target_temperature, tmp * 1e-5, timestep,
                             relaxation_time);

        if (i % output_interval == 0) {
            avg_stress /= output_interval;
            domain.disable(atoms);
            if (MPI::comm_rank(MPI_COMM_WORLD) == 0) {
                csv.write(i / output_interval, i * timestep, kin, pot, tmp, avg_stress, (i / 100) * length_increase);
                if (!quiet) {
                    std::cout << "frame: " << (i / 100) << " kin: " << kin << " pot: " << pot
                            << " tot: " << kin + pot << " tmp: " << tmp << " stress:" << avg_stress << std::endl;
                }
            }
            domain.enable(atoms);
            avg_stress = 0.0;
        }
        if (i % 100 == 0) {
            auto new_size = Eigen::Array3d{domain.domain_length(0), domain.domain_length(1), domain.domain_length(2) + length_increase};
            domain.scale(atoms, new_size);
        }
    }

    traj.close();

    return 0;
}
