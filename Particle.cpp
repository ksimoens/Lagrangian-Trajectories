#include "Particle.h"
#include "globParams.h"

Particle::Particle(){
	pos = Vec(0.0,0.0);
	vel = Vec(0.0,0.0);
}

Particle::Particle(float x0,float y0){
	pos = Vec(x0,y0);
	trans_pos();
	vel = Vec(0.0,0.0);
}

float Particle::fun_lon(float x0, float lat){

	return( x0*sqrt(1.0-E2*pow(sin(lat),2))/A/cos(lat) );

}

float Particle::fun_lat(float mu){

	return( mu + (3.0*E1/2.0-27.0*pow(E1,3)/32.0)*sin(2.0*mu) + (21.0*pow(E1,2)/16.0-55.0*pow(E1,4)/32.0)*sin(4.0*mu) + (151.0*pow(E1,3)/96.0)*sin(6.0*mu) + (1097.0*pow(E1,4)/512.0)*sin(8.0*mu) );

}

float Particle::get_mu(float y0){

	return( y0/A/(1.0-E2/4.0 - 3.0*pow(E2,2)/64.0 - 5.0*pow(E2,3)/256.0) );

}

void Particle::trans_pos(){

	float mu = get_mu(pos.getY());
	float lat = fun_lat(mu);

	pos.setX(fun_lon(pos.getX(),lat));
	pos.setY(mu_lat(lat));

}

//Particle::RK_move(Vec* velgrid){



//}