#include "Particle.h"
#include "globParams.h"

Particle::Particle(){
	pos = Vec(0.0,0.0);
	path = new Vec[NYEAR*365];
	path[0] = pos;
}

Particle::Particle(float x0,float y0){
	pos = Vec(x0,y0);
	trans_pos();
	path = new Vec[NYEAR*365];
	path[0] = pos;
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

int Particle::get_lon_index(Vec pos0){

	int i = 0;

	if(pos0.getX() < LONMIN || pos0.getX() > LONMAX){
		i = -1;
	} else{
		i = floor((pos0.getX()-LONMIN)/LONRES);
	}

	return(i);

}

int Particle::get_lat_index(Vec pos0, float* mus){

	int i = 0;

	if(pos0.getY() < MUMIN || pos0.getY() > MUMAX){
		i = -1;
	} else{
		for(int j=0;j < NLAT;j++){
			if(pos0.getY() < mus[j]){
				break;
			}
			i = j;
		}
	}

	return(i);

}

Vec Particle::interpol(Vec pos0,Vec* velgrid,float* mus,int k,int t){
	
	int i = get_lon_index(pos0);
	int j = get_lat_index(pos0,mus);
	Vec intervel;

	if(i == -1 || j == -1){
		intervel.setX(-999.0);
		intervel.setY(-999.0);
		return(intervel);
	}

	Vec edges[4];
	float x1 = LONMIN+LONRES*i;
	float x2 = x1+LONRES;
	float y1 = mus[j];
	float y2 = mus[j+1];

	edges[0] = velgrid[i+NLON*(j+(t+k)*NLAT)];
	edges[1] = velgrid[i+NLON*(j+1+(t+k)*NLAT)];
	edges[2] = velgrid[(i+1)+NLON*(j+1+(t+k)*NLAT)];
	edges[3] = velgrid[(i+1)+NLON*(j+(t+k)*NLAT)];

	intervel = 1.0/LONRES/(y2-y1)*
				(edges[0]*(x2-pos0.getX())*(y2-pos0.getY()) + 
					edges[1]*(x2-pos0.getX())*(pos0.getY()-y1) +
					edges[2]*(pos0.getX()-x1)*(pos0.getY()-y1) +
					edges[3]*(pos0.getX()-x1)*(y2-pos0.getY()));

	return(intervel);

}

Vec Particle::interpol(Vec pos0,Vec* velgrid,float* mus,int t){

	int i = get_lon_index(pos0);
	int j = get_lat_index(pos0,mus);
	Vec intervel;

	if(i == -1 || j == -1){
		intervel.setX(-999.0);
		intervel.setY(-999.0);
		return(intervel);
	}

	Vec edges[4];
	float x1 = LONMIN+LONRES*i;
	float x2 = x1+LONRES;
	float y1 = mus[j];
	float y2 = mus[j+1];

	edges[0] = (velgrid[i+NLON*(j+NLAT*t)] +
					velgrid[i+NLON*(j+NLAT*(t+1))])/2;
	edges[1] = (velgrid[i+NLON*((j+1)+NLAT*t)] + 
					velgrid[i+NLON*(j+1+NLAT*(t+1))])/2;
	edges[2] = (velgrid[(i+1)+NLON*(j+1+NLAT*t)] +
					velgrid[(i+1)+NLON*(j+1+NLAT*(t+1))])/2;
	edges[3] = (velgrid[(i+1)+NLON*(j+NLAT*t)] + 
					velgrid[(i+1)+NLON*(j+NLAT*(t+1))])/2;

	intervel = 1.0/LONRES/(y2-y1)*
				(edges[0]*(x2-pos0.getX())*(y2-pos0.getY()) + 
					edges[1]*(x2-pos0.getX())*(pos0.getY()-y1) +
					edges[2]*(pos0.getX()-x1)*(pos0.getY()-y1) +
					edges[3]*(pos0.getX()-x1)*(y2-pos0.getY()));

	return(intervel);

}

float Particle::lat_mu(float mu){

	return(M_PI/2.0-2.0*atan(exp(-mu)));

}



void Particle::RK_move(Vec* velgrid,float* mus,int t){

	float K = sqrt(2*D);
	std::normal_distribution<float> norm(0.0,sqrt(DT));
	Vec dW;
	dW.setX(norm(rng));
	dW.setY(norm(rng));

	Vec v1 = interpol(pos,velgrid,mus,0,t);
	float num_v1 = 1.0/R/cos(lat_mu(pos.getY()));
	Vec p2 = pos + (DT/2.0*v1 + K/2.0*dW)*num_v1;

	Vec v2 = interpol(p2,velgrid,mus,t);
	float num_v2 = 1.0/R/cos(lat_mu(p2.getY()));
	Vec p3 = pos + (DT/2.0*v2 + K/2.0*dW)*num_v2;

	Vec v3 = interpol(p3,velgrid,mus,t);
	float num_v3 = 1.0/R/cos(lat_mu(p3.getY()));
	Vec p4 = pos + (DT*v3 + K/2.0*dW)*num_v3;

	Vec v4 = interpol(p4,velgrid,mus,1,t);
	float num_v4 = 1.0/R/cos(lat_mu(p4.getY()));

	pos += DT/6.0*(v1*num_v1 + 2.0*v2*num_v2 + 2.0*v3*num_v3 + v4*num_v4) +
			K*dW/6.0*(num_v1 + 2.0*num_v2 + 2.0*num_v3 + num_v4);	

	if(pos.getX() < -10.0){
		pos.setX(-999.0);
		pos.setY(-999.0);
	}

	path[t] = pos;

} 

void Particle::make_trajectory(Vec* velgrid,float* mus){

	for(int t=0;t<NYEAR*365-1;t++){

		RK_move(velgrid,mus,t);

	}

}