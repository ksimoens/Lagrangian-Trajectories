#include "Grid.h"
#include "globParams.h"

Grid::Grid(){

	vels = new Vec[NLON*NLAT*calc_ndays(NYEAR+NYEARSTART+YSTART)]();
	//velslice = new Vec[NLON*NLAT*2]();
	particles = new Particle[NPART]();
	mus = new float[NLAT]();

	fill_vels();
	initial_particles();
	get_mus();
}

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
				vels[ilon+NLON*(ilat+NLAT*(day+nday_before))] = 
					Vec(grid_time_x[ilat][ilon],grid_time_y[ilat][ilon]);
			}
		}
     
	}


}

void Grid::initial_particles(){

	std::uniform_real_distribution<float> unif(0, 1);

	#pragma omp parallel for
	for(int i=0;i < NPART;i++){
		particles[i] = Particle(CELLLONMIN+100*unif(rng),CELLLATMIN+100*unif(rng));
	}

	//std::cout << "particles" << std::endl;

}

void Grid::get_mus(){

	for(int i=0;i < NLAT;i++){
		mus[i] = mu_lat((LATMIN+i*VELRES)/180.0*M_PI);
	}

}

/*
void Grid::get_time_slice(int t){

	#pragma omp parallel for collapse (3)
	for(int k=0;k < 2;k++){
		for(int j=0;j < NLAT;j++){
			for(int i=0;i < NLON;i++){
				velslice[i+NLON*(j+NLAT*k)] = vels[i+NLON*(j+NLAT*(k+t))];
			}
		}
	}

}

void Grid::timestep(int t){

	get_time_slice(t);

	#pragma omp parallel for
	for(int i=0;i<NPART;i++){
		particles[i].RK_move(velslice,mus,t+1);
	}

}
*/

void Grid::do_simulation(){

	//std::cout << "parallel" << std::endl;
	#pragma omp parallel for num_threads(4)
	for(int i=0;i<NPART;i++){
		//std::cout << omp_get_num_threads() << std::endl;
		//std::cout << i << std::endl;
		particles[i].make_trajectory(vels,mus);
	}

}