#include "Grid.h"
#include "globParams.h"

Grid::Grid(){
	vels = new Vec[NLON*NLAT*365*(NYEAR)]();
	//velslice = new Vec[NLON*NLAT*2]();
	particles = new Particle[NPART]();
	mus = new float[NLAT]();

	fill_vels();
	initial_particles();
	get_mus();
}


void Grid::fill_vels(){

	std::uniform_real_distribution<float> unif(0, 0.2);

	#pragma omp parallel for num_threads(4)
	//std::cout << omp_get_num_threads() << std::endl;
    for(int i=0;i < NLON*NLAT*365*(NYEAR);i++){
    	//std::cout << omp_get_num_threads() << std::endl;
		vels[i] = Vec(unif(rng)-0.1,unif(rng)-0.1);
	}

	//std::cout << vels[60+NLON*(84+NLAT*127)].getX() << std::endl;

}

void Grid::initial_particles(){

	std::uniform_real_distribution<float> unif(0, 1);

	#pragma omp parallel for
	for(int i=0;i < NPART;i++){
		particles[i] = Particle(CELLLONMIN+100*unif(rng),CELLLATMIN+100*unif(rng));
	}

	std::cout << "particles" << std::endl;

}

void Grid::get_mus(){

	for(int i=0;i < NLAT;i++){
		mus[i] = mu_lat((LATMIN+i*VELRES)/180.0*M_PI);
	}

}

/*
void Grid::get_time_slice(int t){

	#pragma omp parallel for collapse (3)
	for(int k=0;k < 2;k++){
		for(int j=0;j < NLAT;j++){
			for(int i=0;i < NLON;i++){
				velslice[i+NLON*(j+NLAT*k)] = vels[i+NLON*(j+NLAT*(k+t))];
			}
		}
	}

}

void Grid::timestep(int t){

	get_time_slice(t);

	#pragma omp parallel for
	for(int i=0;i<NPART;i++){
		particles[i].RK_move(velslice,mus,t+1);
	}

}
*/

void Grid::do_simulation(){

	std::cout << "parallel" << std::endl;
	#pragma omp parallel for num_threads(4)
	for(int i=0;i<NPART;i++){
		//std::cout << omp_get_num_threads() << std::endl;
		std::cout << i << std::endl;
		particles[i].make_trajectory(vels,mus);
	}

}