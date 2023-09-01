#ifndef __ATOMS_H
#define __ATOMS_H

#include <cmath>
#include <algorithm>
#include "types.h"
 
class Atoms { 
public: 
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;
    Names_t names;

    // Atoms(const Atoms& a): positions{a.positions}, velocities{a.velocities}, forces{a.forces}, masses{a.masses}, names{a.names}, atom_mass{a.atom_mass} {}

    Atoms(size_t nb_atoms) : 
            positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms}, names(nb_atoms) {
        positions.setZero();
        velocities.setZero(); 
        forces.setZero();
        masses.fill(atom_mass);
        std::fill(names.begin(), names.end(), "H");
    }
 
    Atoms(const Positions_t &p) : 
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()}, names(p.cols()) { 
        velocities.setZero(); 
        forces.setZero();
        masses.fill(atom_mass);
        std::fill(names.begin(), names.end(), "H");
    }

    Atoms(const Names_t &names, Positions_t &p) : 
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()}, names{names} { 
        velocities.setZero(); 
        forces.setZero();
        masses.fill(atom_mass);
    }
 
    Atoms(const Positions_t &p, const Velocities_t &v) : 
            positions{p}, velocities{v}, forces{3, p.cols()}, masses{p.cols()}, names(p.cols()) { 
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.fill(atom_mass);
        std::fill(names.begin(), names.end(), "H");
    }

    Atoms(const Names_t &names, const Positions_t &p, const Velocities_t &v) : 
            positions{p}, velocities{v}, forces{3, p.cols()}, masses{p.cols()}, names{names} { 
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.fill(atom_mass);
    }

    /*
    Return the number of atoms
    */
    size_t nb_atoms() const { 
        return positions.cols(); 
    }

    /*
    Resize an Eigen::Array3Xd while keeping the information
    */
    void resize(Eigen::Array3Xd &array, int cols) {
        size_t copy_size = std::min<size_t>(array.cols(), cols);
        Eigen::Array3Xd copy = Eigen::Array3Xd{array};
        array.resize(3, cols);
        array(Eigen::all, Eigen::seq(0, copy_size - 1)) = copy(Eigen::all, Eigen::seq(0, copy_size - 1));
    }

    /*
    Resize the datastructure
    */
    void resize(size_t max_size) {
        resize(positions, max_size);
        resize(velocities, max_size);
        resize(forces, max_size);
        masses.resize(max_size);
        masses.fill(atom_mass);
    }

    /*
    Set the mass of all atoms
    */
    void set_mass(double mass) {
        atom_mass = mass;
        masses.fill(mass);
    }

    /*
    Calculates the kinectic energy in eV
    */
    double kinetic_energy() const {
        return atom_mass * velocities.cwiseAbs2().sum() / 2.0;
    }

    /*
    Calculates the kinectic energy in eV from local domain
    */
    double kinetic_energy(int nb_local) const {
        return atom_mass * velocities(Eigen::all, Eigen::seq(0, nb_local - 1)).cwiseAbs2().sum() / 2.0;
    }

    /*
    Return temperature in Kelvin * 1e-5
    */
    double current_temperature() const {
        return (2.0 / 3.0) * (kinetic_energy() / (nb_atoms() * BOLTZMANN_CONSTANT));
    }

    /*
    Return temperature in Kelvin * 1e-5
    */
    double current_temperature(int nb_local) const {
        return (2.0 / 3.0) * (kinetic_energy(nb_local) / (nb_local * BOLTZMANN_CONSTANT));
    }

    private:
    double atom_mass = 1.0;
};

#endif  // __ATOMS_H
