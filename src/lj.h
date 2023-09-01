#ifndef __LJ_H
#define __LJ_H

#include "atoms.h"

/*
Calculate the Lennard-Jones potential between all atoms.
*/
double lj_all(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);


/*
Calculate the Lennard-Jones potential between atoms in cutoff range.
*/
double lj_cutoff(Atoms &atoms, double cutoff, double epsilon = 1.0, double sigma = 1.0);

#endif  // __LJ_H