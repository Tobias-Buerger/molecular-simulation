#ifndef __TYPES_H
#define __TYPES_H

#include <Eigen/Dense>
#include <vector>

/*
For more robust computation the constant is shifted by 1e5
*/
const double BOLTZMANN_CONSTANT = 8.617333262;

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Masses_t = Eigen::ArrayXd;
using Names_t = std::vector<std::string>;

#endif  // __TYPES_H