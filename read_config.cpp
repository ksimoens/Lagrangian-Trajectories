#include "read_config.h"

struct config_params read_config(){

	// open the configuration file
	std::ifstream in(strname(CONFIGPATH));

	// variables for the file line
	std::string line;

	std::size_t pos;

	// define the structure for the configuration parameters
	struct config_params myparams;

	// read the lines in the configuration file
	//while(!in.eof()){
	while(getline( in, line )){

		std::stringstream ss(line);
		std::string s;
		std::string params[2];
		int i=0;
		while (getline(ss, s, ' ')) {
		// store token string in the vector
        	params[i] = s;
        	i++;
    	}

		// the horizontal initial position (degrees)
		if(params[0] == "x0"){
			myparams.x0 = std::stof(params[1],&pos);
		// the vertical initial position (degrees)
		} else if(params[0] == "y0"){
			myparams.y0 = std::stof(params[1],&pos);
		} else if(params[0] == "srf"){
			myparams.r = sqrt(std::stof(params[1],&pos)/M_PI);
		} else if(params[0] == "write"){
			myparams.w = params[1];
			std::replace(myparams.w.begin(),myparams.w.end(),'%',' ');
		} else if(params[0] == "velocity"){
			myparams.v = params[1];
			std::replace(myparams.v.begin(),myparams.v.end(),'%',' ');
		}

	}

	return(myparams); 

}