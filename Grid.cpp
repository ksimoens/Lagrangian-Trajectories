#include "Grid.h"
#include "Particle.h"
#include "globParams.h"

Grid::Grid(float x0,float y0){

	this->vels = new Vec[NLON*NLAT*calc_ndays(NYEAR+NYEARSTART+YSTART)]();
	//velslice = new Vec[NLON*NLAT*2]();
	this->particles = new Particle[NPART*(calc_ndays(NYEARSTART+YSTART)/DTSTART)];
	this->mus = new float[NLAT]();
	this->pos0 = Vec(x0,y0);
	this->radius = 0.0;

	//fill_vels();
	initial_particles();
	get_mus();
}

#ifdef CIRCULAR
Grid::Grid(float x0,float y0,float r){

	this->vels = new Vec[NLON*NLAT*calc_ndays(NYEAR+NYEARSTART+YSTART)]();
	//velslice = new Vec[NLON*NLAT*2]();
	this->particles = new Particle[NPART*(calc_ndays(NYEARSTART+YSTART)/DTSTART)]();
	this->mus = new float[NLAT]();
	this->pos0 = Vec(x0,y0);
	this->radius = r;

	fill_vels();
	get_mus();
	initial_particles();

}
#endif

size_t Grid::calc_ndays(int current_year){

	size_t nday = 0;
	for(int year=YSTART;year<current_year;year++){

		if((year%4==0) & (year!=2000)){
			nday += 366;
		} else{
			nday += 365;
		}

	}

	return(nday);

}

void Grid::fill_vels(){

	for(int i=0;i<NYEAR+NYEARSTART;i++){
		fill_vels_year(i);
	}

}

void Grid::fill_vels_year(int year){

	size_t nday;
	size_t nday_before = calc_ndays(year+YSTART);

	if(((year+YSTART)%4==0) & ((year+YSTART)!=2000)){
		nday = 366;
	} else{
		nday = 365;
	}

	float grid_time_x[NLAT][NLON];
	float grid_time_y[NLAT][NLON];

	std::string ystr = std::to_string(year+YSTART);
	netCDF::NcFile dataFile("../../network/curr_vel_transf/vel_"+std::to_string(year+YSTART)+".nc", netCDF::NcFile::read);

	netCDF::NcVar velxVar;
	velxVar = dataFile.getVar("u");
	netCDF::NcVar velyVar;
	velyVar = dataFile.getVar("v");

	std::vector<size_t> startp,countp;
	startp.push_back(0);
	startp.push_back(0);
	startp.push_back(0);
	countp.push_back(1);
	countp.push_back(NLAT);
	countp.push_back(NLON);

	for (size_t day = 0; day < nday; day++){
		// Read the data one record at a time.
		startp[0]=day;
		velxVar.getVar(startp,countp,grid_time_x);
		velyVar.getVar(startp,countp,grid_time_y);
     
		for(int ilon=0;ilon<NLON;ilon++){
			for(int ilat=0;ilat<NLAT;ilat++){
				this->vels[ilon+NLON*(ilat+NLAT*(day+nday_before))] = 
					Vec(grid_time_x[ilat][ilon],grid_time_y[ilat][ilon]);
			}
		}
     
	}


}

void Grid::initial_particles(){

	std::uniform_real_distribution<float> unif(0, 1);
	#pragma omp parallel for
	for(size_t i=0;i<(calc_ndays(NYEARSTART+YSTART)/DTSTART);i++){
		for(int j=0;j<NPART;j++){
			float r1 = unif(rng);
			float r2 = unif(rng);
			this->particles[i*NPART+j].get_initial_pos(pos0,r1,r2,radius,i,vels,mus);
		}
	}

}

void Grid::get_mus(){

	float mus_in[NLAT];

	netCDF::NcFile dataFile("../../network/curr_vel_transf/vel_1993.nc", netCDF::NcFile::read);

	netCDF::NcVar muVar;
	muVar = dataFile.getVar("mu");
	muVar.getVar(mus_in);

	for(int i=0;i<NLAT;i++){
		this->mus[i] = mus_in[i];
	}

}

void Grid::do_simulation(){

	#pragma omp parallel for num_threads(4)
	for(size_t i=0;i<(calc_ndays(NYEARSTART+YSTART)/DTSTART);i++){
		for(int j=0;j<NPART;j++){
			this->particles[j+NPART*i].make_trajectory(vels,mus);
		}
	}

}