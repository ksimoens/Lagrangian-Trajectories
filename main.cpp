#include "globParams.h"
#include "Grid.h"

int main(){

	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
	rng.seed(1);

	Grid grid;
	for(int i=0;i<365*NYEAR-1;i++){
		grid.timestep(i);
		std::cout << i << std::endl;
	}

	return 0;

}