#include "read_config.h"

struct config_params read_config(){

	// open the configuration file
	std::ifstream in(strname(CONFIGPATH));

	// variables for the parameter name
	std::string param;

	std::size_t pos;

	// define the structure for the configuration parameters
	struct config_params myparams;

	// read the lines in the configuration file
	while(!in.eof()){

		// read the first word on the line
		in >> param;

		// the horizontal initial position (degrees)
		if(param == "x0"){
			std::string x0;
			in >> x0;
			myparams.x0 = std::stof(x0,&pos);
		// the vertical initial position (degrees)
		} else if(param == "y0"){
			std::string y0;
			in >> y0;
			myparams.y0 = std::stof(y0,&pos);
		} else if(param == "srf"){
			std::string r;
			in >> r;
			myparams.r = sqrt(std::stof(r,&pos)/M_PI);
		}

	}

	return(myparams); 

}