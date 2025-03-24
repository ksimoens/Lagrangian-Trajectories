#include "Grid.h"
#include "globParams.h"

Grid::Grid(){
	vels = new Vec[NLON*NLAT*365*(NYEAR+NYEARSTART)]();
	particles = new Particle[NPART]();
	mus = new float[NLAT]();

	fill_vels();
	initial_particles();
	get_mus();
}


void Grid::fill_vels(){

	std::uniform_real_distribution<float> unif(0, 1);

    for(int i=0;i < NLON*NLAT*365*(NYEAR+NYEARSTART);i++){
		*(vels + i) = Vec(unif(rng),unif(rng));
	}

}

void Grid::initial_particles(){

	std::uniform_real_distribution<float> unif(0, 1);

	for(int i=0;i < NPART;i++){
		*(particles+i) = Particle(CELLLONMIN+100*unif(rng),CELLLATMIN+100*unif(rng));
	}

}

void Grid::get_mus(){

	for(int i=0;i < NLAT;i++){
		mus[i] = mu_lat((LATMIN+i*VELRES)/180.0*M_PI);
	}

}

Vec* Grid::get_time_slice(int t){

	Vec* slice = new Vec[NLON*NLAT];

	for(int i=0;i < NLON;i++){
		for(int j=0;j < NLAT;j++){
			slice[i+NLON*j] = vels[i+NLON*(j+NLAT*t)];
		}
	}

	return slice;

}
