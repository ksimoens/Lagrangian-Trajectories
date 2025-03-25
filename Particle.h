#ifndef PARTICLE_H
#define PARTICLE_H

#include "globParams.h"
#include "Vec.h"

class Particle{

	private:
		Vec pos{};

		void trans_pos();
		float fun_lon(float x0,float lat);
		float fun_lat(float mu);
		float get_mu(float y0);
		int get_lon_index(Vec pos0);
		int get_lat_index(Vec pos0, float* mus);
		Vec interpol(Vec pos0,Vec* velslice,float* mus);
		Vec interpol(Vec pos0,Vec* velslice,float* mus,int t);
		float lat_mu(float mu);

	public:
		Particle();
		Particle(float x0, float y0);

		Vec getPos(){return pos;};
		void setPos(Vec pos0){pos = pos0;};

		void RK_move(Vec* velgrid, float* mus);

};

#endif