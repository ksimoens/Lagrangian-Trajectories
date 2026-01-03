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

		#ifdef NETWORK
			struct hlp_coord{
				int basehp=-1;
				int x=-1;
				int y=-1;
			};
			float xy_to_lon();
			void update_network(int t,std::set<int> IDvec,int* network,int Nstart,int i,int j);
			int get_network_id(std::set<int> IDvec);
			float lonmu_to_x();
			float lonmu_to_y();
			float lonmu_to_z();
			void get_hlp_coord(int* basehp,int* x,int* y,float vx,float vy,float vz);
			int get_pixel_id(hlp_coord hlp);
		#endif

	public:
		Particle();
		Particle(float x0, float y0, int t0);
		~Particle(){delete[] path_vel;path_vel=0;delete[] path_pos;path_pos=0;};

		Vec getPos(){return pos;};
		Vec* getPathPos(){return path_pos;};
		Vec* getPathVel(){return path_vel;};
		void setPos(Vec pos0){pos = pos0;};
		int get_starttime(){return this->starttime;};
		void set_starttime(int t0){this->starttime=t0;};
		float lat_mu(float mu);

		void RK_move(Vec* velgrid, float* mus, int t);

		#ifdef NETWORK
			void make_trajectory(Vec* velgrid, float* mus, std::set<int> IDvec, int* network, int Nstart, int i, int j);
		#else
			void make_trajectory(Vec* velgrid, float* mus);
		#endif

		void get_initial_pos(Vec pos0,float r1,float r2,float r0,int t0,Vec* velgrid,float* mus);

		#ifdef NETWORK
			void xy_to_lonmu();
		#endif

};

#endif