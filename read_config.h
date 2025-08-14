#ifndef CONFIG
#define CONFIG

/*
	structure that contains the configuration-file parameters
	function that reads the parameters from the configuration file
*/

#include "globParams.h"

struct config_params{

	float x0;
	float y0;
	float r;
	std::string w;

};

// defined in "read_config.cpp"
struct config_params read_config();

#endif