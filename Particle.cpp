#include "Particle.h"
#include "globParams.h"

Particle::Particle(){
	pos = Vec(0.0,0.0);
	vel = Vec(0.0,0.0);
}

Particle::Particle(float x0,float y0){
	pos = Vec(x0,y0);
	vel = Vec(0.0,0.0);
}

//Particle::RK_move(Vec* velgrid){



//}