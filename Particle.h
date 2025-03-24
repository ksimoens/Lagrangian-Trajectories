#ifndef PARTICLE_H
#define PARTICLE_H

#include "globParams.h"
#include "Vec.h"

class Particle{

	private:
		Vec pos{};
		Vec vel{};

		float get_vel(Vec* velgrid);

	public:
		Particle();
		Particle(float x0, float y0);

		void RK_move(Vec* velgrid);

};

#endif