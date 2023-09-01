#include "verlet.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep, double mass) {
    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;
    x += vx * timestep;
    y += vy * timestep;
    z += vz * timestep;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep, double mass) {
    vx += 0.5 * fx * timestep / mass;
    vy += 0.5 * fy * timestep / mass;
    vz += 0.5 * fz * timestep / mass;
}

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep, double mass) {
    velocities += 0.5 * forces * timestep / mass;
    positions += velocities * timestep;
}

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double timestep, double mass) {
    velocities += 0.5 * forces * timestep / mass;
}

void verlet_step1(Atoms &atoms, double timestep, double mass) {
    verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, mass);
}

void verlet_step2(Atoms &atoms, double timestep, double mass) {
    verlet_step2(atoms.velocities, atoms.forces, timestep, mass);
}