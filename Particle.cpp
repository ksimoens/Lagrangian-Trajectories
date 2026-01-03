#include "Particle.h"
#include "globParams.h"

Particle::Particle(){
	this->pos = Vec(0.0,0.0);
	this->starttime = 0;
	#ifdef STOREPOS
		this->path_pos = new Vec[NYEAR*365+1];
		this->path_pos[0] = this->pos;
	#else
		this->path_pos = 0;
	#endif

	#ifdef STOREVEL
		this->path_vel = new Vec[NYEAR*365+1];
	#else 
		this->path_vel = 0;
	#endif
}

Particle::Particle(float x0,float y0,int t0){
	this->pos = Vec(x0,y0);
	this->starttime = t0;
	#ifndef NETWORK
	trans_pos();
	#endif
	#ifdef STOREPOS
		this->path_pos = new Vec[NYEAR*365+1];
		this->path_pos[0].setX(this->pos.getX());
		this->path_pos[0].setY(this->pos.getY());
	#else
		this->path_pos = 0;
	#endif

	#ifdef STOREVEL
		this->path_vel = new Vec[NYEAR*365+1];
	#else 
		this->path_vel = 0;
	#endif
}

void Particle::get_initial_pos(Vec pos0,float r1,float r2,float r0,int t0,Vec* vels,float* mus){

	#ifdef CIRCULAR
		float r_rand = r0*sqrt(r1);
		float h_rand = 2.0*M_PI*r2;
		this->pos.setX(pos0.getX()+r_rand*cos(h_rand));
		this->pos.setY(pos0.getY()+r_rand*sin(h_rand));
		trans_pos();
		this->starttime = t0;
	#endif

	#ifdef STOREVEL
	/*	Vec vel_inter = interpol(this->pos,vels,mus,0,this->starttime);
		this->path_vel[0].setX(vel_inter.getX());
		this->path_vel[0].setY(vel_inter.getY());*/
	#endif

	#ifdef STOREPOS
		this->path_pos[0].setX(this->pos.getX());
		this->path_pos[0].setY(this->pos.getY());
	#endif

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

	float mu = get_mu(this->pos.getY());
	float lat = fun_lat(mu);

	this->pos.setX(fun_lon(this->pos.getX(),lat));
	this->pos.setY(mu_lat(lat));

}

#ifdef NETWORK
void Particle::xy_to_lonmu(){

	Vec newpos = Vec(0.0,0.0);

	if(abs(this->pos.getY()) >= M_PI_4){
		newpos.setX( (this->pos.getX() - (abs(this->pos.getY())-M_PI_4)/(abs(this->pos.getY())-M_PI_2)*(std::fmod(this->pos.getX(),M_PI_2)-M_PI_4) - M_PI ) );
		newpos.setY( mu_lat(M_PI_2 - acos((1.0-1.0/3.0*pow(2.0-4.0*abs(this->pos.getY()/M_PI),2))*sgn(this->pos.getY()))) );
	} else{
		newpos.setX(this->pos.getX() - M_PI);
		newpos.setY( mu_lat(M_PI_2 - acos(8.0/3.0/M_PI*this->pos.getY())) );
	}

	this->pos.setX(newpos.getX());
	this->pos.setY(newpos.getY());

}
#endif

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
				i = j-1;
				break;
			}
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

	return(M_PI_2-2.0*atan(exp(-mu)));

}



void Particle::RK_move(Vec* velgrid,float* mus,int t){

	float K = sqrt(2*D);
	std::normal_distribution<float> norm(0.0,sqrt(DT));
	Vec dW;
	dW.setX(norm(rng));
	dW.setY(norm(rng));

	Vec v1 = interpol(this->pos,velgrid,mus,0,t);
	#ifdef STOREVEL
		path_vel[t-this->starttime] = v1;
	#endif
	float num_v1 = 1.0/R/cos(lat_mu(this->pos.getY()));
	Vec p2 = this->pos + (DT/2.0*v1 + K/2.0*dW)*num_v1;

	Vec v2 = interpol(p2,velgrid,mus,t);
	float num_v2 = 1.0/R/cos(lat_mu(p2.getY()));
	Vec p3 = this->pos + (DT/2.0*v2 + K/2.0*dW)*num_v2;

	Vec v3 = interpol(p3,velgrid,mus,t);
	float num_v3 = 1.0/R/cos(lat_mu(p3.getY()));
	Vec p4 = this->pos + (DT*v3 + K/2.0*dW)*num_v3;

	Vec v4 = interpol(p4,velgrid,mus,1,t);
	float num_v4 = 1.0/R/cos(lat_mu(p4.getY()));

	this->pos += DT/6.0*(v1*num_v1 + 2.0*v2*num_v2 + 2.0*v3*num_v3 + v4*num_v4) +
			K*dW/6.0*(num_v1 + 2.0*num_v2 + 2.0*num_v3 + num_v4);	

	if(v1.getX() < -100.0 || v2.getX() < -100.0 || v3.getX() < -100.0 || v4.getX() < -100.0){
		this->pos.setX(-999.0);
		this->pos.setY(-999.0);
	}
	
	#ifdef STOREPOS
		this->path_pos[t-this->starttime+1] = this->pos;
	#endif 

} 

#ifdef NETWORK
float Particle::lonmu_to_x(){

	return(cos(this->pos.getX()+M_PI)*cos(lat_mu(this->pos.getY())));

}

float Particle::lonmu_to_y(){

	return(sin(this->pos.getX()+M_PI)*cos(lat_mu(this->pos.getY())));

}

float Particle::lonmu_to_z(){

	return(sin(lat_mu(this->pos.getY())));

}

void Particle::get_hlp_coord(int* basehp,int* x,int* y,float vx,float vy,float vz){

	float phi, phi_t, sector, xx, yy;
	int offset, column;
	float TwTh = 2.0/3.0;
	float root3 = sqrt(3.0);

	phi = atan2(vy,vx);
	if(phi < 0){
		phi += 2.0*M_PI;
	}
	phi_t = fmod(phi,M_PI_2);

	// definitely polar 
	if((vz >= TwTh) || (vz <= -TwTh)){

		bool north;
		float coz;
		float kx,ky;

		// north
		if(vz > 0){
			north = true;
		// south
		}else{
			north = false;
			vz *= -1.0;
		}

		coz = std::hypot(vx,vy);

		kx = (coz / sqrt(1.0 + vz)) * root3 * fabs(NSIDE * (2.0 * phi_t - M_PI) / M_PI);
		ky = (coz / sqrt(1.0 + vz)) * root3 * NSIDE * 2.0 * phi_t / M_PI;

		if (north) {
            xx = NSIDE - kx;
            yy = NSIDE - ky;
        } else {
            xx = ky;
            yy = kx;
        }

        *x = std::min(NSIDE-1, (int)floor(xx));
        *y = std::min(NSIDE-1, (int)floor(yy));

        sector = (phi - phi_t) / (M_PI_2);
        offset = (int)round(sector);
        offset = ((offset % 4) + 4) % 4;
        column = offset;

        if(north){
            *basehp = column;
        } else {
            *basehp = 8 + column;
        }

    // polar or equatorial
	} else{

		float zunits, phiunits, u1, u2;

		zunits = (vz + TwTh) / (4.0 / 3.0);
        phiunits = phi_t / M_PI_2;

		u1 = zunits + phiunits;
		u2 = zunits - phiunits + 1.0;

		// x is the northeast direction, y is the northwest.
		xx = u1 * NSIDE;
		yy = u2 * NSIDE;

		sector = (phi - phi_t) / (M_PI_2);
        offset = (int)round(sector);
        offset = ((offset % 4) + 4) % 4;

        // we're looking at a square in z,phi space with an X dividing it.
        // we want to know which section we're in.
        // xx ranges from 0 in the bottom-left to 2Nside in the top-right.
        // yy ranges from 0 in the bottom-right to 2Nside in the top-left.
        // (of the phi,z unit box)
        if (xx >= NSIDE) {
            xx -= NSIDE;
            // north polar
            if (yy >= NSIDE) {
                yy -= NSIDE;
                *basehp = offset;
            // right equatorial
            } else {
                *basehp = ((offset + 1) % 4) + 4;
            }
        } else {
        	// left equatorial
            if (yy >= NSIDE) {
                yy -= NSIDE;
                *basehp = offset + 4;
            // south polar
            } else {
                *basehp = 8 + offset;
            }
        }

		*x = std::max(0, std::min(NSIDE-1, (int)floor(xx)));
        *y = std::max(0, std::min(NSIDE-1, (int)floor(yy)));

	}

}

int Particle::get_pixel_id(hlp_coord hlp){

	int frow, F1, v, ring, index;

	frow = hlp.basehp / 4;
	F1 = frow + 2;
	v = hlp.x + hlp.y;
	ring = F1*NSIDE - v - 1;

	// north polar
	if (ring <= NSIDE) {

        index = (NSIDE - 1 - hlp.y);
        // offset from the other big healpixes
        index += ((hlp.basehp % 4) * ring);
        // offset from the other rings
        index += (int)ring*(ring-1)*2;

       // south polar
    } else if (ring >= (int)3*NSIDE) {

        // Here I first flip everything so that we label the pixels
        // at zero starting in the southeast corner, increasing to the
        // west and north, then subtract that from the total number of
        // healpixels.
        int ri = (int)4*NSIDE - ring;
        index = (ri-1) - hlp.x;
        // big healpixes
        index += ((3-(hlp.basehp % 4)) * ri);
        // other rings
        index += (int)ri*(ri-1)*2;
        // flip!
        index = 12*(int)NSIDE*NSIDE - 1 - index;

    // equatorial
    } else {
        int s, F2, h;
        s = (ring - NSIDE) % 2;
        F2 = 2*((int)hlp.basehp % 4) - (frow % 2) + 1;
        h = hlp.x - hlp.y;

        index = (F2 * NSIDE + h + s) / 2;
        // offset from the north polar region:
        index += (int)NSIDE * (NSIDE - 1) * 2;
        // offset within the equatorial region:
        index += (int)NSIDE * 4 * (ring - NSIDE);
        // handle healpix #4 wrap-around
        if ((hlp.basehp == 4) && (hlp.y > hlp.x)){
        	index += (4 * NSIDE - 1);
        }
    }

    return(index);

}

int Particle::get_network_id(std::set<int> IDvec){

	float vx,vy,vz;
	vx = lonmu_to_x();
	vy = lonmu_to_y();
	vz = lonmu_to_z();

	hlp_coord hlp;
	get_hlp_coord(&hlp.basehp,&hlp.x,&hlp.y,vx,vy,vz);

	int id;
	id = get_pixel_id(hlp);

	return(std::distance(IDvec.begin(),IDvec.find(id)));

}

void Particle::update_network(int t,std::set<int> IDvec,int* network,int Nstart,int i,int j){

	int i_id = get_network_id(IDvec);
	if(i_id != NCELL){
		network[i_id+NCELL*(j+i*NPART)] = std::min(network[i_id+NCELL*(j+i*NPART)],t);
	}

}

void Particle::make_trajectory(Vec* velgrid,float* mus,std::set<int> IDvec,int* network,int Nstart,int i,int j){

	for(int t=0;t<NYEAR*365;t++){
		RK_move(velgrid,mus,t+this->starttime);
		if( (this->pos.getX() > NETLONMIN) & (this->pos.getX() < NETLONMAX) & (this->pos.getY() > NETLATMIN) & (this->pos.getY() < NETLATMAX)){
			update_network(t+1,IDvec,network,Nstart,i,j);
		}
	}

}

#else
void Particle::make_trajectory(Vec* velgrid,float* mus){

	for(int t=0;t<NYEAR*365;t++){
		RK_move(velgrid,mus,t+this->starttime);

	}

	#ifdef STOREVEL
		Vec last_vel = interpol(this->pos,velgrid,mus,0,NYEAR*365+this->starttime);
		this->path_vel[NYEAR*365].setX(last_vel.getX());
		this->path_vel[NYEAR*365].setY(last_vel.getY());
	#endif

}
#endif