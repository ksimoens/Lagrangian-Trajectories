#ifndef PARTICLE_H
#define PARTICLE_H

#include "globParams.h"
#include "Vec.h"

class Particle{

	private:
		Vec pos{};
		Vec* path_pos{};
		Vec* path_vel{};
		int starttime{};

		void trans_pos();
		float fun_lon(float x0,float lat);
		float fun_lat(float mu);
		float get_mu(float y0);
		int get_lon_index(Vec pos0);
		int get_lat_index(Vec pos0, float* mus);
		Vec interpol(Vec pos0,Vec* velgrid,float* mus,int t);
		Vec interpol(Vec pos0,Vec* velgrid,float* mus,int k,int t);
		float lat_mu(float mu);

	public:
		Particle();
		Particle(float x0, float y0, int t0);
		~Particle(){delete[] path_vel;path_vel=0;delete[] path_pos;path_pos=0;};

		Vec getPos(){return pos;};
		Vec* getPathPos(){return path_pos;};
		Vec* getPathVel(){return path_vel;};
		void setPos(Vec pos0){pos = pos0;};

		void RK_move(Vec* velgrid, float* mus, int t);

		void make_trajectory(Vec* velgrid, float* mus);

		#ifdef STOREVEL
			void get_initial_pos(Vec pos0,float r1,float r2,float r0,int t0,Vec* velgrid,float* mus);
		#endif

};

#endif