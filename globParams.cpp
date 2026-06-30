#include "globParams.h"

float mu_lat(float lat){

	return(log(std::abs(1.0/cos(lat) + tan(lat))));

};

float lat_mu(float mu){

	return(M_PI_2-2.0*atan(exp(-mu)));

};

float sgn(float val){

	return((0.0 < val) - (val < 0.0));

};