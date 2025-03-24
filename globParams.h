#ifndef SIM_GLOBALCONSTS_H
#define SIM_GLOBALCONSTS_H

#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

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
#define LONMIN -81.875
#define LONMAX 9.875
#define VELRES 0.25
#define LATMIN 20.125
#define LATMAX 69.875

#define CELLLONMIN -5563.64787134093
#define CELLLATMIN 4590.6176765106

#endif //SIM_GLOBALCONSTS_H