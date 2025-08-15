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
			this->particles[i*NPART+j].get_initial_pos(pos0,r1,r2,radius,i*DTSTART,vels,mus);
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

void Grid::write_simulation(std::string w,double dt_init,double dt_sim){

	auto t_start = std::chrono::high_resolution_clock::now();
	netCDF::NcFile data("output/"+w+".nc", netCDF::NcFile::replace);

	data.putAtt("title","Lagrangian simulation Northern Atlantic Ocean");
	time_t timestamp;
	time(&timestamp);
	data.putAtt("clock time",ctime(&timestamp));
	Particle p = Particle(this->pos0.getX(),this->pos0.getY(),0);
	data.putAtt("central starting position","("+
					std::to_string(p.getPos().getX()/M_PI*180)+" ; "+
					std::to_string(p.lat_mu(p.getPos().getY())/M_PI*180)+") degrees");
	#ifdef CIRCULAR
		data.putAtt("initial condition","circular");
		data.putAtt("radius",std::to_string(this->radius)+" km");
	#endif


	#if defined(STOREPOS) || defined(STOREVEL)
		data.putAtt("simulation type","full trajectory");

		netCDF::NcDim startDim = data.addDim("start", calc_ndays(NYEARSTART+YSTART)/DTSTART);
		netCDF::NcDim partDim = data.addDim("particles", NPART);
		netCDF::NcDim yearsDim = data.addDim("years", NYEAR);
		netCDF::NcDim daysDim = data.addDim("days", 365);

		std::vector<netCDF::NcDim> dimVector;
		dimVector.push_back(startDim);
		dimVector.push_back(yearsDim);
		dimVector.push_back(daysDim);
		dimVector.push_back(partDim);

		std::vector<netCDF::NcDim> dimVector_start;
		dimVector_start.push_back(startDim);
		std::vector<netCDF::NcDim> dimVector_years;
		dimVector_years.push_back(yearsDim);
		std::vector<netCDF::NcDim> dimVector_days;
		dimVector_days.push_back(daysDim);
		std::vector<netCDF::NcDim> dimVector_init;
		dimVector_init.push_back(startDim);
		dimVector_init.push_back(partDim);

		netCDF::NcVar yearsVar = data.addVar("years", netCDF::ncInt, dimVector_years);
		yearsVar.putAtt("units", "years");
		int vec_years[NYEAR];
		for(int i=0;i<NYEAR;i++){
			vec_years[i] = i;
		}
		std::vector<size_t> startp_years,countp_years;
		startp_years.push_back(0);
		countp_years.push_back(NYEAR);
		yearsVar.putVar(startp_years,countp_years,vec_years);

		netCDF::NcVar daysVar = data.addVar("days", netCDF::ncInt, dimVector_days);
		daysVar.putAtt("units", "days");
		int vec_days[365];
		for(int i=0;i<365;i++){
			vec_days[i] = i;
		}
		std::vector<size_t> startp_days,countp_days;
		startp_days.push_back(0);
		countp_days.push_back(365);
		daysVar.putVar(startp_days,countp_days,vec_days);

		netCDF::NcVar initxVar = data.addVar("lon_0", netCDF::ncFloat, dimVector_init);
		initxVar.putAtt("units", "radians");
		netCDF::NcVar inityVar = data.addVar("mu_0", netCDF::ncFloat, dimVector_init);
		inityVar.putAtt("units", "radians");
		int vec_initx[calc_ndays(NYEARSTART+YSTART)/DTSTART][NPART];
		int vec_inity[calc_ndays(NYEARSTART+YSTART)/DTSTART][NPART];
		for(size_t i=0;i<calc_ndays(NYEARSTART+YSTART)/DTSTART;i++){
			for(int j=0;j<NPART;j++){
				vec_initx[i][j] = particles[i*NPART+j].getPathPos()[0].getX();
				vec_inity[i][j] = particles[i*NPART+j].getPathPos()[0].getY();
			}
		}
		std::vector<size_t> startp_init,countp_init;
		startp_init.push_back(0);
		startp_init.push_back(0);
		countp_init.push_back(calc_ndays(NYEARSTART+YSTART)/DTSTART);
		countp_init.push_back(NPART);
		initxVar.putVar(startp_init,countp_init,vec_initx);
		inityVar.putVar(startp_init,countp_init,vec_inity);

		std::vector<size_t> startp,countp;
		startp.push_back(0);
		startp.push_back(0);
		startp.push_back(0);
		startp.push_back(0);
		countp.push_back(1);
		countp.push_back(1);
		countp.push_back(365);
		countp.push_back(NPART);

		std::vector<size_t> startp_start,countp_start;
		startp_start.push_back(0);
		countp_start.push_back(calc_ndays(NYEARSTART+YSTART)/DTSTART);

		std::string store_str = "";
		#ifdef STOREPOS
			store_str += "position ";

			netCDF::NcVar posxVar = data.addVar("lon", netCDF::ncFloat, dimVector);
			posxVar.putAtt("units", "radians");

			float vec_part_posx[365][NPART];
			netCDF::NcVar posyVar = data.addVar("mu", netCDF::ncFloat, dimVector);
			posyVar.putAtt("units", "radians");
			float vec_part_posy[365][NPART];
		#endif

		#ifdef STOREVEL
			store_str += " velocity";

			netCDF::NcVar velxVar = data.addVar("u", netCDF::ncFloat, dimVector);
			velxVar.putAtt("units", "m/s");
			float vec_part_velx[365][NPART];
			netCDF::NcVar velyVar = data.addVar("v", netCDF::ncFloat, dimVector);
			velyVar.putAtt("units", "m/s");
			float vec_part_vely[365][NPART];
		#endif

		data.putAtt("stored data",store_str);

		netCDF::NcVar startVar = data.addVar("start", netCDF::ncInt, dimVector_start);
		startVar.putAtt("units", "days");
		int vec_start[calc_ndays(NYEARSTART+YSTART)/DTSTART];

		//#pragma omp parallel for num_threads(4)
		for(size_t i=0;i < (calc_ndays(NYEARSTART+YSTART)/DTSTART);i++){
			vec_start[i] = i*DTSTART;
			startp[0]=i;

			for(int l=0;l<NYEAR;l++){
				startp[1]=l;

				for(int j=0;j<365;j++){

					for(int k=0;k < NPART;k++){

						#ifdef STOREPOS
							vec_part_posx[j][k] = particles[i*NPART+k].getPathPos()[l*365+j+1].getX();
							vec_part_posy[j][k] = particles[i*NPART+k].getPathPos()[l*365+j+1].getY();
						#endif
						#ifdef STOREVEL
							vec_part_velx[j][k] = particles[i*NPART+k].getPathVel()[l*365+j+1].getX();
							vec_part_vely[j][k] = particles[i*NPART+k].getPathVel()[l*365+j+1].getY();
						#endif
					}				
				}

				#ifdef STOREPOS
					posxVar.putVar(startp,countp,vec_part_posx);
					posyVar.putVar(startp,countp,vec_part_posy);
				#endif
				#ifdef STOREVEL
					velxVar.putVar(startp,countp,vec_part_velx);
					velyVar.putVar(startp,countp,vec_part_vely);
				#endif
			}
		}

		startVar.putVar(startp_start,countp_start,vec_start);

		data.putAtt("length of trajectories",std::to_string(NYEAR)+" years");
		std::string start_str;
		if(NYEARSTART == 1){
			start_str = " year";
		} else{
			start_str = " years";
		}
		data.putAtt("length of initialisation",strname(NYEARSTART)+start_str);
		data.putAtt("time between releases",std::to_string(DTSTART)+" days");
		data.putAtt("number of particles per release",strname(NPART));
		data.putAtt("starting date","01-01-"+std::to_string(YSTART));
		data.putAtt("diffusion constant",std::to_string(D)+" m^2/s");
		auto t_end = std::chrono::high_resolution_clock::now();
		double dt_writ = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		data.putAtt("initialisation wall time",std::to_string(dt_init/1000)+" seconds");
		data.putAtt("simulation wall time",std::to_string(dt_sim/1000)+" seconds");
		data.putAtt("output wall time",std::to_string(dt_writ/1000)+" seconds");

	#endif

}