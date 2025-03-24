#ifndef SIM_GLOBALCONSTS_H
#define SIM_GLOBALCONSTS_H

#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <cmath>

// Total number of years to integrate.
#define NYEAR 1
// Total number of years over which to start.
#define NYEARSTART 1
#define NPART 100
// Starting year.
//const int Ystart = 1993;
// Diffusion constant.
#define D 100.0
// File with grid cell information.
//#define CELLS "cells_test.csv"
#define NLON 368
#define NLAT 200
#define LONMIN -1.428988325
#define LONMAX 0.172351264
#define LONRES 0.004363323
#define VELRES 0.25
#define LATMIN 20.125
#define LATMAX 69.875

#define CELLLONMIN -5563.64787134093
#define CELLLATMIN 4590.6176765106

#define A 6378.140
#define B 6356.7500
#define E2 0.0066960376522835485
#define E1 0.0016796375940427949

extern std::mt19937_64 rng;

extern float mu_lat(float lat);

#endif //SIM_GLOBALCONSTS_H