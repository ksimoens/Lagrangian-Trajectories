#include "globParams.h"
#include "Grid.h"

int main(){

	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
	rng.seed(ss);

	Grid grid;

	for(int i=0;i<365*NYEAR;i++){
		grid.timestep(0);
		std::cout << grid.get_particles()[0].getPos().getX() << std::endl;
	}
	
	return 0;

}