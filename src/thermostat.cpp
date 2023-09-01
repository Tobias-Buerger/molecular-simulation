#include "thermostat.h"


void berendsen_thermostat(Atoms &atoms, double target_temperature, double current_temperature, double timestep,
                          double relaxation_time) {
    double lambda = std::sqrt(1.0 + (target_temperature / current_temperature - 1.0)
                                    * (timestep / relaxation_time));
    atoms.velocities *= lambda;
}