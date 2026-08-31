#ifndef PARTICLE_H
#define PARTICLE_H

#include "globParams.h"
#include "Vec.h"

class Particle{

	private:
		Vec pos{};
		#if defined(STOREPOS) || defined(LYAPUNOV)
			Vec* path_pos{};
		#endif
		#ifdef STOREVEL
			Vec* path_vel{};
		#endif
		int* velmask{};
		float* vecnum{};
		Vec* vecintervel{};
		Vec* posintermed{};
		int starttime{};
		#ifdef BROWNIAN
			Vec pos0{};
			float* distances{};
		#endif
		#ifdef LYAPSST
			float* path_SST{};
		#endif
		#ifdef LYAPCIRC
			int tau{};
			float radius{};
		#endif
		#ifdef LYAPRATIO
			int tau{};
			int tau_SST{};
			float SSTdiff0{};
			float radius{};
		#endif
		#ifdef LYAPINFINITY
			float radius{};
		#endif
		#ifdef TRACERPATH
			float* path_SST{};
			float* path_pos{};
		#endif
		#ifdef DISTRIBUTION
			int ntarget{};
			Particle* targets{};
			#ifdef DISTVEL
				float sumvel;
				float* arrivals{};
			#else
				int* arrivals{};
			#endif
		#endif

		void trans_pos();
		float fun_lon(float x0,float lat);
		float fun_lat(float mu);
		float get_mu(float y0);
		#if defined(SST) || defined(BROWNIAN) || defined(LYAPCIRC) || defined(LYAPRATIO) || defined(LYAPINFINITY) || defined(TRACERPATH) || defined(DISTRIBUTION)
			float fun_x(float lon,float lat);
			float fun_y(float lat);
		#endif
		int get_lon_index(Vec pos0,float lonmin,float lonmax,float lonres);
		int get_lat_index(float lat,float latmin,float latmax,float latres);
		Vec interpol(Vec pos0,float lat,Vec* velgrid,int t);
		Vec interpol(Vec pos0,float lat,Vec* velgrid,int k,int t);
		#ifndef LYAPCIRC
			#ifndef LYAPRATIO
				#ifndef TRACERPATH
					void RK_move(Vec* velgrid, int t,Vec dW);
				#endif
			#endif
		#endif
		void update_pos(float K,Vec dW);
		#if defined(BROWNIAN) || defined(DISTRIBUTION)
			float haversine(Vec pos1);
		#endif

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

		void free_memory();

	public:
		Particle();
		//Particle(float x0, float y0, int t0);
		/*#if defined(STOREPOS) || defined(LYAPUNOV)
		~Particle(){delete[] path_vel;path_vel=0;delete[] path_pos;path_pos=0;
					delete[] velmask;velmask=0;delete[] vecnum;vecnum=0;
					delete[] vecintervel;vecintervel=0;
					delete[] posintermed;posintermed=0;
					delete[] distances;distances=0;
					delete[] path_SST;path_SST=0;};*/
		~Particle(){free_memory();};

		Vec getPos(){return pos;};
		#if defined(STOREPOS) || defined(LYAPUNOV)
			Vec* getPathPos(){return path_pos;};
		#endif
		#ifdef STOREVEL
			Vec* getPathVel(){return path_vel;};
		#endif
		#ifdef LYAPSST
			float* getPathSST(){return path_SST;};
		#endif
		#if defined(LYAPCIRC) || defined(LYAPRATIO)
			int getTau(){return this->tau;};
			void setTau(int t){this->tau = t;};
			float getRadius(){return this->radius;};
		#endif
		#if defined(LYAPRATIO)
			int getTauSST(){return this->tau_SST;};
			void setTauSST(int t){this->tau_SST = t;};
			void setSSTdiff0(float dSST0){this->SSTdiff0 = dSST0;};
			float getSSTdiff0(){return(this->SSTdiff0);};
		#endif
		#ifdef BROWNIAN
			float* getDistances(){return this->distances;};
		#endif
		#ifdef LYAPINFINITY
			float getRadius(){return this->radius;};
		#endif
		#ifdef TRACERPATH
			void setSSTdiff(float dSST0,int t){this->path_SST[t] = dSST0;};
			void setDistance(Vec pos1,int t);
			float* getPathSST(){return path_SST;};
			float* getPathPos(){return path_pos;};
		#endif
		#ifdef DISTRIBUTION
			void set_targets(Particle* vec_target,int ntarget);
			Particle* get_targets(){return this->targets;};
			#ifdef DISTVEL
				float* get_arrivals(){return this->arrivals;};
			#else
				int* get_arrivals(){return this->arrivals;};
			#endif
		#endif
		void setPos(Vec pos0){pos = pos0;};
		int get_starttime(){return this->starttime;};
		void set_starttime(int t0){this->starttime=t0;};
		
		#if defined(SST) || defined(LYAPSST) || defined(LYAPRATIO) || defined(LYAPINFINITY) || defined(TRACERPATH)
			float interpol(float* SSTgrid,int t);
		#endif

		#ifdef NETWORK
			void make_trajectory(Vec* velgrid, std::set<int> IDvec, int* network, int Nstart, int i, int j,std::mt19937_64 &rng);
		#endif
		#if defined(CIRCULAR)
			void make_trajectory(Vec* velgrid,std::mt19937_64 &rng);
		#endif
		#if defined(LYAPUNOV) || defined(SST) || defined(LYAPINFINITY)
			void make_trajectory(Vec* velgrid,std::mt19937_64 &rng,int Ntime);
		#endif
		#ifdef BROWNIAN
			void make_trajectory(Vec* velgrid,std::mt19937_64 &rng,int* outtimes);
		#endif
		#ifdef LYAPSST
			void make_trajectory(Vec* velgrid,float* SSTs,std::mt19937_64 &rng,int Ntime);
		#endif
		#if defined(DISTRIBUTION)
			void make_trajectory(Vec* velgrid,std::mt19937_64 &rng,float r);
		#endif
		#if defined(LYAPCIRC) || defined(LYAPRATIO) || defined(TRACERPATH)
			void RK_move(Vec* velgrid, int t,Vec dW);
		#endif

		#if defined(LYAPCIRC) || defined(SST) || defined(LYAPRATIO) || defined(LYAPINFINITY) || defined(TRACERPATH)
			float haversine(Vec pos1);
		#endif

		void get_initial_pos(Vec pos0,float r1,float r2,float r0,int t0);

		#ifdef NETWORK
			void xy_to_lonmu();
		#endif

};

#endif