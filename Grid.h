#ifndef GRID_H
#define GRID_H

#include "globParams.h"
#include "Vec.h"
#include "Particle.h"

class Grid{

	private:
		Vec* vels;
		//Vec* velslice;
		Particle* particles;
		float* mus;
		Vec pos0;
		float radius;
		
		void fill_vels();
		void initial_particles();
		size_t calc_ndays(int current_year);
		void fill_vels_year(int year);
		void get_mus();
		//void get_time_slice(int t);

	public:
		Grid(float x0,float y0);
		#ifdef CIRCULAR
			Grid(float x0,float y0,float r);
		#endif
		~Grid(){delete[] vels; delete[] mus;vels=0; mus=0;};

		//void timestep(int t);
		Particle* get_particles(){return particles;};
		void do_simulation();
		void set_pos0(float x0,float y0){pos0=Vec(x0,y0);};
		#ifdef CIRCULAR
			void set_radius(float r){radius=r;};
		#endif

};

#endif