#include "verlet.h"
#include "types.h"
#include <gtest/gtest.h>
#include <Eigen/Dense>

// Test if a bouncing ball simulation conserves energy with doubles
TEST(VerletDoubleTest, EnergyConservingBall) {
    double x, y, z, vx, vy, vz, fx, fy, fz;
    x = z = vx = vz = fx = fz = 0;
    y = 5;
    fy = -10;
    vy = 0;
    double timestep = 1e-3;
    for(int step = 0; step < 10000; step++) {
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);
        // no force update
        verlet_step2(vx, vy, vz, fx, fy, fz, timestep);
        if(y <= 0)
            vy = -vy;
    }

    double energy_before = 1 * 10 * 5;
    double energy_after = 1 * 10 * y + 0.5 * 1 * vy * vy;

    ASSERT_NEAR(energy_after, energy_before, 1e-3);
}

// Test if a bouncing balls simulation conserves energy with Eigen
TEST(VerletEigenTest, EnergyConservingBall) {
    const size_t ball_count = 50;
    Positions_t positions(3, ball_count);
    Velocities_t velocities(3, ball_count);
    Forces_t forces(3, ball_count);
    double timestep = 1e-3;
    positions.fill(0);
    velocities.fill(0);
    forces.fill(0);
    forces.row(1).fill(-10);
    positions.row(1).fill(5);
    double pos0 = 0;
    for (auto &&pos : positions.row(0)) {
        pos = pos0++;
    }
    
    for(size_t step = 0; step < 10000; step++) {
        verlet_step1(positions, velocities, forces, timestep);
        // no force update
        verlet_step2(velocities, forces, timestep);
        for (size_t i = 0; i < ball_count; i++) {
            if (positions(1, i) <= 0) {
                velocities(1, i) *= -1;
            }
        }
    }

    double energy_before_sum = 1 * 10 * 5 * ball_count;
    auto energy_after{1 * 10 * positions.row(1) + 0.5 * 1 * velocities.row(1).pow(2)};
    double energy_after_sum = energy_after.sum();

    ASSERT_NEAR(energy_after_sum, energy_before_sum, 1e-3);
}
