#include "globParams.h"
#include "Grid.h"
#include "read_config.h"

int main(){

	auto t_start = std::chrono::high_resolution_clock::now();

	struct config_params myparams = read_config();
	#ifdef CIRCULAR
		Grid grid = Grid(myparams.x0,myparams.y0,myparams.r,myparams.v);
	#elif NETWORK
		Grid grid = Grid(myparams.x0,myparams.y0,myparams.v,myparams.net);
	#elif LYAPUNOV
		Grid grid = Grid(myparams.v);
	#endif
	
	auto t_end = std::chrono::high_resolution_clock::now();
	double dt_init = std::chrono::duration<double, std::milli>(t_end-t_start).count();
	
	t_start = std::chrono::high_resolution_clock::now();
	grid.do_simulation();
	t_end = std::chrono::high_resolution_clock::now();
	double dt_sim = std::chrono::duration<double, std::milli>(t_end-t_start).count();

	std::ifstream file((myparams.w+".nc").c_str());
	if(file.good()){
		std::remove((myparams.w+".nc").c_str());
	}
	
	grid.write_simulation(myparams.w,dt_init,dt_sim);
	
	return 0;

}