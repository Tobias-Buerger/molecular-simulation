#ifndef __THERMOSTAT_H
#define __THERMOSTAT_H

#include "types.h"
#include "atoms.h"


/*
Run Berendsen thermostat for a specific timestep
*/
void berendsen_thermostat(Atoms &atoms, double target_temperature, double current_temperature, double timestep,
                          double relaxation_time);


#endif  // __THERMOSTAT_H