#include "lj.h"
#include "neighbors.h"
#include <cmath>
#include <Eigen/Dense>

inline constexpr double w(double r, double epsilon, double sigma) {
    return 4 * epsilon *
           (std::pow(sigma / r, 12) - std::pow(sigma / r, 6));
}

inline constexpr double dr_w(double r, double epsilon, double sigma) {
    return 4 * epsilon *
           (6 * std::pow(1.0 / r, 7) * std::pow(sigma, 6) -
            12 * std::pow(1.0 / r, 13) * std::pow(sigma, 12));
}

double lj_all(Atoms &atoms, double epsilon, double sigma) {
    double potential = 0.0;
    for (size_t k = 0; k < atoms.nb_atoms(); k++) {
        Eigen::Vector3d force;
        force.setZero();
        for (size_t i = 0; i < atoms.nb_atoms(); i++) {
            if (i == k) {
                continue;
            }
            Eigen::Vector3d r_ik_vec = atoms.positions.col(i) - atoms.positions.col(k);
            double r_ik = r_ik_vec.norm();
            force += dr_w(r_ik, epsilon, sigma) * r_ik_vec.normalized();
            potential += w(r_ik, epsilon, sigma);
        }
        atoms.forces.col(k) = force;
    }
    return potential / 2.0;
}


double lj_cutoff(Atoms &atoms, double cutoff, double epsilon, double sigma) {
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);
    double energy_shift = - w(cutoff, epsilon, sigma);
    double potential = 0.0;
    atoms.forces.setZero();
    for (auto[k, i]: neighbor_list) {
        Eigen::Vector3d force = atoms.forces.col(k);
        Eigen::Vector3d r_ik_vec = atoms.positions.col(i) - atoms.positions.col(k);
        double r_ik = r_ik_vec.norm();
        force += dr_w(r_ik, epsilon, sigma) * r_ik_vec / r_ik;
        potential += w(r_ik, epsilon, sigma) + energy_shift;
        atoms.forces.col(k) = force;
    }
    return potential / 2.0;
}