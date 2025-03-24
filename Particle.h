#ifndef PARTICLE_H
#define PARTICLE_H

#include "globParams.h"
#include "Vec.h"

class Particle{

	private:
		Vec pos{};
		Vec vel{};

		float get_vel(Vec* velgrid);
		void trans_pos();
		float fun_lon(float x0,float lat);
		float fun_lat(float mu);
		float get_mu(float y0);

	public:
		Particle();
		Particle(float x0, float y0);

		Vec getPos(){return pos;};
		Vec getVel(){return vel;};

		void RK_move(Vec* velgrid);

};

#endif