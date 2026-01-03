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
		float* mus;
		Vec pos0;
		float radius;
		int* network;
		int Nstart;
		std::set<int> IDvec;
		
		void fill_vels(std::string veldir);
		void initial_particles();
		size_t calc_ndays(int current_year);
		void fill_vels_year(int year,std::string veldir);
		void get_mus(std::string veldir);
		//void get_time_slice(int t);

		#ifdef NETWORK
			void get_cell_ids(std::string netdir);
			void initial_network();
		#endif

	public:
		Grid(float x0,float y0,std::string veldir);
		#ifdef CIRCULAR
			Grid(float x0,float y0,float r,std::string veldir);
		#endif
		#ifdef NETWORK
			Grid(float x0,float y0,float r,std::string veldir,std::string netdir);
		#endif
		~Grid(){delete[] vels; delete[] mus; delete[] network; vels=0; mus=0; network=0;};

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