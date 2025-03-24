#include "globParams.h"

std::mt19937_64 rng;

float mu_lat(float lat){

	return(log(std::abs(1.0/cos(lat) + tan(lat))));

};