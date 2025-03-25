#ifndef GRID_H
#define GRID_H

#include "globParams.h"
#include "Vec.h"
#include "Particle.h"

class Grid{

	private:
		Vec* vels;
		Particle* particles;
		float* mus;
		
		void fill_vels();
		void initial_particles();
		void get_mus();
		Vec* get_time_slice(int t);

	public:
		Grid();
		~Grid(){delete[] vels; delete[] particles; delete[] mus;};

		void timestep(int t);
		Particle* get_particles(){return particles;};

};

#endif