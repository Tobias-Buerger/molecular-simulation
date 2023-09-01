#include "verlet.h"
#include "types.h"
#include "thermostat.h"
#include "lj.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>

// Test if velocities go near temperature 1
TEST(BerendsenTest, testTargetTemperature) {
    double target_temperature = 1.0;
    int nb_atoms = 10;
    Atoms atoms(nb_atoms);
    atoms.positions.setRandom();
    atoms.velocities.setRandom();

    for (size_t i = 0; i < 10'000; i++)
    {
        berendsen_thermostat(atoms, target_temperature, atoms.current_temperature(), 0.1, 1.0);
    }
    
    EXPECT_NEAR(atoms.current_temperature(), target_temperature, 1e-5);
}

// Test simulation with potential
TEST(BerendsenTest, testSimulation) {
    double target_temperature = 100.0;
    int nb_atoms = 10;
    Atoms atoms(nb_atoms);
    atoms.positions.setRandom();
    atoms.velocities.setRandom();

    double mass = 1.0;
    double epsilon = 1.0;
    double sigma = 1.0;
    double total_time = 10.0 * std::sqrt(mass * sigma * sigma / epsilon);
    double delta_time = 0.001 * std::sqrt(mass * sigma * sigma / epsilon);
    size_t timesteps = total_time / delta_time;
    double relaxation_time = 1000.0 * delta_time;

    for (size_t i = 0; i < timesteps; i++) {
        verlet_step1(atoms, delta_time, mass);
        double pot = lj_all(atoms, epsilon, sigma);
        verlet_step2(atoms, delta_time, mass);
        
        berendsen_thermostat(atoms, target_temperature, atoms.current_temperature(), delta_time, relaxation_time);
    }
    
    EXPECT_NEAR(atoms.current_temperature(), target_temperature, 1.0);
}