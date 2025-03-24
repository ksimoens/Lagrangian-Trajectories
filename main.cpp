#include "globParams.h"
#include "Vec.h"
#include "Particle.h"


int main(){

	std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    std::uniform_real_distribution<float> unif(0, 1);

    int Ntot = NLON*NLAT*(NYEARSTART+NYEAR)*365;

	Vec* velgrid = new Vec[Ntot];

	for(int i=0;i < Ntot;i++){
		*(velgrid + i) = Vec(unif(rng),unif(rng));
	}

	Particle* particles = new Particle[NPART];

	for(int i=0;i < NPART;i++){
		particles[i] = Particle(1.0,1.0);
	}

	delete[] velgrid;
	delete[] particles;

	return 0;

}