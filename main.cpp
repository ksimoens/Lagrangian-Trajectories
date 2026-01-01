#include "globParams.h"
#include "Grid.h"
#include "read_config.h"

int main(){

	auto t_start = std::chrono::high_resolution_clock::now();

	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
	rng.seed(1);

	struct config_params myparams = read_config();
	#ifndef NETWORK
		Grid grid = Grid(myparams.x0,myparams.y0,myparams.r,myparams.v);
	#else
		Grid grid = Grid(myparams.x0,myparams.y0,myparams.r,myparams.v,myparams.net);
	#endif
	/*
	auto t_end = std::chrono::high_resolution_clock::now();
	double dt_init = std::chrono::duration<double, std::milli>(t_end-t_start).count();
	
	t_start = std::chrono::high_resolution_clock::now();
	grid.do_simulation();
	t_end = std::chrono::high_resolution_clock::now();
	double dt_sim = std::chrono::duration<double, std::milli>(t_end-t_start).count();

	std::ifstream file(("output/"+myparams.w+".nc").c_str());
	if(file.good()){
		std::remove(("output/"+myparams.w+".nc").c_str());
	}

	grid.write_simulation(myparams.w,dt_init,dt_sim);
	*/
	return 0;

}