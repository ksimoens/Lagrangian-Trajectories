#ifndef GRID_H
#define GRID_H

#include "globParams.h"
#include "Vec.h"
#include "Particle.h"

class Grid{

	private:
		Vec* vels;
		//Vec* velslice;
		Particle* particles;
		Vec pos0;
		float radius;
		int* network;
		int Nstart;
		std::set<int> IDvec;
		float* SSTbeg;
		float* SSTend;
		float* SSTs;
		int* outtimes;
		float* vecR;
		int npart;
		
		void fill_vels(std::string veldir);
		void initial_particles();
		#ifdef DAY
			size_t calc_ndays(int current_year);
			void fill_vels_year(int year,std::string veldir);
		#elif HOUR
			size_t calc_nhours(int current_month,int current_year);
			void fill_vels_month(int year,int month,std::string veldir);
		#endif
		#if defined(LYAPSST) || defined(LYAPRATIO)
			size_t calc_ndays(int current_month,int current_year);
		#endif
		//void get_mus(std::string veldir);
		//void get_time_slice(int t);

		#ifdef NETWORK
			void get_cell_ids(std::string netdir);
			void initial_network();
		#endif

		#if defined(LYAPUNOV) || defined(SST)
			float euclidean(Vec pos0,Vec pos1);
		#endif

		#ifdef SST
			void fill_SSTs(std::string SSTbegdir,std::string SSTenddir);
		#endif

		#if defined(LYAPSST) || defined(LYAPRATIO)
			float get_SSTdiff(float SST0,float SST1);
		#endif
		#if defined(LYAPSST) || defined(LYAPRATIO)
			void fill_SSTs(std::string SSTbegdir);
			void fill_SSTs_month(int year,int month,std::string SSTdir);
		#endif

	public:
		#ifdef CIRCULAR
			Grid(float x0,float y0,float r,std::string veldir);
		#endif
		#ifdef NETWORK
			Grid(float x0,float y0,std::string veldir,std::string netdir);
		#endif
		#if defined(LYAPUNOV) || defined(LYAPCIRC)
			Grid(std::string veldir);
		#endif
		#ifdef SST
			Grid(float r,std::string veldir,std::string SSTbegdir,std::string SSTenddir);
		#endif
		#ifdef BROWNIAN
			Grid(float r,std::string veldir);
		#endif
		#if defined(LYAPSST) || defined(LYAPRATIO)
			Grid(std::string veldir,std::string SSTbegdir);
		#endif
		~Grid(){delete[] vels; vels = 0; 
				delete[] network; network = 0; 
				delete[] particles; particles = 0; 
				delete[] SSTbeg; SSTbeg = 0; 
				delete[] SSTend; SSTend = 0; 
				delete[] SSTs; SSTs = 0;
				delete[] vecR; vecR = 0;};

		//void timestep(int t);
		Particle* get_particles(){return particles;};
		void do_simulation();
		void set_pos0(float x0,float y0){pos0=Vec(x0,y0);};
		#ifdef CIRCULAR
			void set_radius(float r){radius=r;};
		#endif
		void write_simulation(std::string w,double dt_init,double dt_sim);

};

#endif