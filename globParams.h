#include <string>

#ifndef SIM_GLOBALCONSTS_H
#define SIM_GLOBALCONSTS_H

// Total number of years to integrate.
const int Nyear = 1;
// Total number of years over which to start.
const int Nyearstart = -1;
// Starting year.
const int Ystart = 1993;
// Diffusion constant.
const float D = 100.0;
// File with grid cell information.
const std::string cells = "cells_test.csv";

#endif //SIM_GLOBALCONSTS_H