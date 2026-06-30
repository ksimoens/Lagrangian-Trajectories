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
		
		void fill_vels(std::string veldir);
		void initial_particles();
		#ifdef DAY
			size_t calc_ndays(int current_year);
			void fill_vels_year(int year,std::string veldir);
		#elif HOUR
			size_t calc_nhours(int current_month);
			void fill_vels_month(int month,std::string veldir);
		#endif
		//void get_mus(std::string veldir);
		//void get_time_slice(int t);

		#ifdef NETWORK
			void get_cell_ids(std::string netdir);
			void initial_network();
		#endif

		#ifdef LYAPUNOV
			float haversine(Vec pos0,Vec pos1);
			float euclidean(Vec pos0,Vec pos1);
		#endif

	public:
		#ifdef CIRCULAR
			Grid(float x0,float y0,float r,std::string veldir);
		#elif NETWORK
			Grid(float x0,float y0,std::string veldir,std::string netdir);
		#elif LYAPUNOV
			Grid(std::string veldir);
		#endif
		~Grid(){delete[] vels; delete[] network; delete[] particles; particles = 0; vels=0; network=0;};

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