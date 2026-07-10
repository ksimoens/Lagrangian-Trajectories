#include "Grid.h"
#include "Particle.h"
#include "globParams.h"

#ifdef CIRCULAR
Grid::Grid(float x0,float y0,float r,std::string veldir){

	#ifdef DAY
		this->Nstart = calc_ndays(NYEARSTART+YSTART)/DTSTART;
		this->vels = new Vec[NLON*NLAT*calc_ndays(NYEAR+NYEARSTART+YSTART)]();
	#elif HOUR
		this->Nstart = 1;
		this->vels = new Vec[NLON*NLAT*calc_nhours(MSTART,YSTART)]();
	#endif
	//velslice = new Vec[NLON*NLAT*2]();
	this->particles = new Particle[NPART*Nstart]();
	this->pos0 = Vec(x0,y0);
	this->radius = r;
	this->network = 0;
	this->SSTs = 0;

	fill_vels(veldir);
	initial_particles();

}
#endif

#ifdef NETWORK
Grid::Grid(float x0,float y0,std::string veldir,std::string netdir){

	#ifdef DAY
		this->Nstart = calc_ndays(NYEARSTART+YSTART)/DTSTART;
		this->vels = new Vec[NLON*NLAT*calc_ndays(NYEAR+NYEARSTART+YSTART)]();
	#elif HOUR
		this->Nstart = 1;
		this->vels = new Vec[NLON*NLAT*calc_nhours(MSTART,YSTART)]();
	#endif
	//velslice = new Vec[NLON*NLAT*2]();
	this->particles = new Particle[NPART*Nstart]();
	this->pos0 = Vec(x0,y0);
	this->SSTs = 0;

	this->network = new int[NPART*this->Nstart*NCELL]();
	get_cell_ids(netdir);
	initial_network();

	fill_vels(veldir);
	initial_particles();

}
#endif

#ifdef LYAPUNOV
Grid::Grid(std::string veldir){

	#ifdef DAY
		this->Nstart = calc_ndays(NYEARSTART+YSTART)/DTSTART;
		this->vels = new Vec[NLON*NLAT*calc_ndays(NYEAR+NYEARSTART+YSTART)]();
	#elif HOUR
		this->Nstart = 1;
		this->vels = new Vec[NLON*NLAT*calc_nhours(MSTART,YSTART)]();
	#endif
	
	int nlon = (int)((LYAPLONMAX-LYAPLONMIN)/LYAPLONRES);
	int nlat = (int)((LYAPLATMAX-LYAPLATMIN)/LYAPLATRES);
	this->particles = new Particle[(nlon-2)*(nlat-2)]();
	this->network = 0;
	this->SSTs = 0;

	fill_vels(veldir);
	initial_particles();

}
#endif

#ifdef SST
Grid::Grid(float r,std::string veldir,std::string SSTdir){

	#ifdef DAY
		this->Nstart = calc_ndays(NYEARSTART+YSTART)/DTSTART;
		this->vels = new Vec[NLON*NLAT*calc_ndays(NYEAR+NYEARSTART+YSTART)]();
	#elif HOUR
		this->Nstart = 1;
		this->vels = new Vec[NLON*NLAT*calc_nhours(MSTART,YSTART)]();
	#endif
	
	int nlon = (int)((SSTLONMAX-SSTLONMIN)/SSTLONRES);
	int nlat = (int)((SSTLATMAX-SSTLATMIN)/SSTLATRES);
	this->particles = new Particle[(nlon-2)*(nlat-2)*NPART]();
	this->radius = r;
	this->network = 0;
	this->SSTs = new float[SSTGRIDNLON*SSTGRIDNLAT];

	fill_vels(veldir);
	initial_particles();
	fill_SSTs(SSTdir);

}
#endif

#ifdef DAY
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
#elif HOUR
size_t Grid::calc_nhours(int current_month,int current_year){

	if(current_month < 1){
		current_month = 12;
		current_year--;
	}

	size_t nhour = 0;
	
	int nmonth = 0;
	int month_i = MSTART;
	int year_i = YSTART;

	while((month_i != current_month) || (year_i != current_year)){
		nmonth++;
		month_i--;
		if(month_i < 1){
			month_i = 12;
			year_i--;
		}
	}

	nmonth = NMONTH - nmonth;

	month_i = current_month;
	year_i = current_year;

	for(int i=nmonth-1;i >= 0;i--){
		if((month_i == 1) || (month_i == 3) || (month_i == 5) ||
			(month_i == 7) || (month_i == 8) || (month_i == 10) ||
			(month_i == 12)){
			nhour += 24*31;
		} else if((month_i != 2)){
			if((year_i % 4 == 0) & (year_i != 2000)){
				nhour += 24*29;
			} else{
				nhour += 24*28;
			}
		} else{
			nhour += 24*30;
		}
		month_i--;
		if(month_i < 1){
			month_i = 12;
			year_i--;
		}
	}

	return(nhour);

}
#endif

void Grid::fill_vels(std::string veldir){

	#ifdef DAY
		for(int i=0;i<NYEAR+NYEARSTART;i++){
			fill_vels_year(i,veldir);
		}
	#elif HOUR
		
		int vec_month[NMONTH];
		int vec_year[NMONTH];
		int year_i = YSTART;
		int month_i = MSTART;

		vec_month[NMONTH-1] = month_i;
		vec_year[NMONTH-1] = year_i;

		for(int i=NMONTH-2;i>=0;i--){
			month_i -= 1;
			if(month_i < 1){
				month_i = 12;
				year_i -= 1;
			}
			vec_month[i] = month_i;
			vec_year[i] = year_i;
		}

		for(int i=0;i < NMONTH;i++){
			fill_vels_month(vec_year[i],vec_month[i],veldir);
		}

	#endif

}

#ifdef DAY
void Grid::fill_vels_year(int year,std::string veldir){

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
	netCDF::NcFile dataFile(veldir+"/vel_"+std::to_string(year+YSTART)+".nc", netCDF::NcFile::read);
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
#elif HOUR
void Grid::fill_vels_month(int year,int month,std::string veldir){

	size_t nhour;
	size_t nhour_before = calc_nhours(month-1,year);

	if((month == 1) || (month == 3) || (month == 5) ||
		(month == 7) || (month == 8) || (month == 10) ||
		(month == 12)){
		nhour = 24*31;
	} else if((month == 2)){
		if((year % 4 == 0) & (year != 2000)){
			nhour = 24*29;
		} else{
			nhour = 24*28;
		}
	} else{
		nhour = 24*30;
	}

	float grid_time_x[NLAT][NLON];
	float grid_time_y[NLAT][NLON];
	std::string mstr = std::to_string(month);
	size_t n_zero = 2;
	mstr = std::string(n_zero - std::min(n_zero, mstr.length()), '0') + mstr;
	
	netCDF::NcFile dataFile(veldir+"/vel_"+std::to_string(year)+"_"+mstr+".nc", netCDF::NcFile::read);
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

	for (size_t hour = 0; hour < nhour; hour++){
		// Read the data one record at a time.
		startp[0]=hour;
		velxVar.getVar(startp,countp,grid_time_x);
		velyVar.getVar(startp,countp,grid_time_y);
     
		for(int ilon=0;ilon<NLON;ilon++){
			for(int ilat=0;ilat<NLAT;ilat++){
				this->vels[ilon+NLON*(ilat+NLAT*(hour+nhour_before))] = 
					Vec(grid_time_x[ilat][ilon],grid_time_y[ilat][ilon]);
			}
		}
     
	}

}
#endif

void Grid::initial_particles(){

	#ifdef CIRCULAR

	#pragma omp parallel 
	{
		std::random_device rd;
		std::seed_seq seed{static_cast<int>(rd()),omp_get_thread_num()};
		std::mt19937_64 rng(seed);
		std::uniform_real_distribution<float> unif(0, 1);

		#pragma omp for
		for(size_t i=0;i<this->Nstart;i++){
			for(int j=0;j<NPART;j++){
				float r1 = unif(rng);
				float r2 = unif(rng);
				this->particles[i*NPART+j].get_initial_pos(pos0,r1,r2,radius,i*DTSTART);
			}
		}
	}
	#endif

	#ifdef NETWORK
	int Nsub = sqrt(NPART);
	float ds = M_PI_2/NSIDE;
	for(int i=0;i<this->Nstart;i++){
		int j = 0;

			for(int m=1;m <= Nsub;m++){
				float y_k = this->pos0.getY() - ds/2.0 + m*ds/Nsub/2.0;
				for(int n=1;n <= m;n++){
					float x_k = this->pos0.getX() - ((m-1.0)/2.0)*ds/Nsub + (n-1)*ds/Nsub;
					this->particles[i*NPART+j].setPos(Vec(x_k,y_k));
					this->particles[i*NPART+j].xy_to_lonmu();
					this->particles[i*NPART+j].set_starttime(i*DTSTART);
					j++;
				}
			}

			for(int m=1;m <= (Nsub-1);m++){
				float y_k = this->pos0.getY() + ds/2.0 - m*ds/Nsub/2.0;
				for(int n=1;n <= m;n++){
					float x_k = this->pos0.getX() - ((m-1.0)/2.0)*ds/Nsub + (n-1)*ds/Nsub;
					this->particles[i*NPART+j].setPos(Vec(x_k,y_k));
					this->particles[i*NPART+j].xy_to_lonmu();
					this->particles[i*NPART+j].set_starttime(i*DTSTART);
					j++;
				}
			}

	}
	#endif

	#ifdef LYAPUNOV

	int nlon = (int)((LYAPLONMAX-LYAPLONMIN)/LYAPLONRES);
	int nlat = (int)((LYAPLATMAX-LYAPLATMIN)/LYAPLATRES);

	//#pragma omp parallel 
	{
		//#pragma omp for
		int k = 0;
		for(int ilat=1;ilat<(nlat-1);ilat++){
			for(int ilon=1;ilon<(nlon-1);ilon++){
				this->particles[k].get_initial_pos(Vec(LYAPLONMIN+ilon*LYAPLONRES,mu_lat(LYAPLATMIN+ilat*LYAPLATRES)),0.0,0.0,0.0,0);
				k++;
			}
		}
	}

	#endif

	#ifdef SST

	int nlon = (int)((SSTLONMAX-SSTLONMIN)/SSTLONRES);
	int nlat = (int)((SSTLATMAX-SSTLATMIN)/SSTLATRES);

	#pragma omp parallel 
	{

		std::random_device rd;
		std::seed_seq seed{static_cast<int>(rd()),omp_get_thread_num()};
		std::mt19937_64 rng(seed);
		std::uniform_real_distribution<float> unif(0, 1);

		#pragma omp for
		for(int ilat=1;ilat<(nlat-1);ilat++){
			for(int ilon=1;ilon<(nlon-1);ilon++){
				for(int j=0;j<NPART;j++){
					float r1 = unif(rng);
					float r2 = unif(rng);
					this->particles[j+NPART*(ilon-1+(nlon-2)*(ilat-1))].get_initial_pos(Vec(SSTLONMIN+ilon*SSTLONRES,SSTLATMIN+ilat*SSTLATRES),r1,r2,this->radius,0);
				}
			}
		}
	}

	#endif

}

/*void Grid::get_mus(std::string veldir){

	float mus_in[NLAT];

	netCDF::NcFile dataFile(veldir+"/vel_1993.nc", netCDF::NcFile::read);

	netCDF::NcVar muVar;
	muVar = dataFile.getVar("mu");
	muVar.getVar(mus_in);

	for(int i=0;i<NLAT;i++){
		this->mus[i] = mus_in[i];
	}

}*/

#ifdef NETWORK

void Grid::get_cell_ids(std::string netdir){

	std::fstream fin;

	fin.open(netdir+".csv",std::ios::in);

	std::string temp,line,word;

	std::getline(fin,line);

	for(int i = 0;i<NCELL;i++){

		std::getline(fin,line);
		std::stringstream s(line);
		std::getline(s,word,',');
		this->IDvec.insert(std::stoi(word));

	}

}

void Grid::initial_network(){

	for(int i=0;i < NPART*this->Nstart*NCELL;i++){

		this->network[i] = 999999;

	}

}

#endif

#ifdef SST

void Grid::fill_SSTs(std::string SSTdir){

	float grid_SST[SSTGRIDNLAT][SSTGRIDNLON];
	netCDF::NcFile dataFile(SSTdir+".nc", netCDF::NcFile::read);
	netCDF::NcVar SSTVar;
	SSTVar = dataFile.getVar("SST");

	std::vector<size_t> startp,countp;
	startp.push_back(0);
	startp.push_back(0);
	countp.push_back(SSTGRIDNLAT);
	countp.push_back(SSTGRIDNLON);

	SSTVar.getVar(startp,countp,grid_SST);
     
	for(int ilon=0;ilon<SSTGRIDNLON;ilon++){
		for(int ilat=0;ilat<SSTGRIDNLAT;ilat++){
			this->SSTs[ilon+SSTGRIDNLON*ilat] = grid_SST[ilat][ilon];
		}
	}
     

}

#endif

void Grid::do_simulation(){

	#pragma omp parallel
	{

		std::random_device rd;
		std::seed_seq seed{static_cast<int>(rd()),2};//omp_get_thread_num()};
		std::mt19937_64 rng(seed);

		#if defined(CIRCULAR) || defined(NETWORK)
			#pragma omp for
			for(int i=0;i<this->Nstart;i++){
				std::cout << i << std::endl;
				for(int j=0;j<NPART;j++){
					#ifdef NETWORK
						this->particles[j+NPART*i].make_trajectory(this->vels,this->IDvec,this->network,this->Nstart,i,j,rng);
					#else
						this->particles[j+NPART*i].make_trajectory(this->vels,rng);
					#endif
				}
			}
		#endif
		#ifdef LYAPUNOV

			int nlon = (int)((LYAPLONMAX-LYAPLONMIN)/LYAPLONRES);
			int nlat = (int)((LYAPLATMAX-LYAPLATMIN)/LYAPLATRES);
			#pragma omp for
			for(int j=0;j<((nlon-2)*(nlat-2));j++){
			//for(int j=0;j<1000;j++){
					//std::cout << j << std::endl;
					this->particles[j].make_trajectory(this->vels,rng);
			}

		#endif
		#ifdef SST
			int nlon = (int)((SSTLONMAX-SSTLONMIN)/SSTLONRES);
			int nlat = (int)((SSTLATMAX-SSTLATMIN)/SSTLATRES);
			#pragma omp for
			for(int j=0;j<((nlon-2)*(nlat-2)*NPART);j++){
			//for(int j=0;j<1000;j++){
					//std::cout << j << std::endl;
					this->particles[j].make_trajectory(this->vels,rng);
			}

		#endif
	}
}

#if defined(LYAPUNOV) || defined(SST)
float Grid::haversine(Vec pos0,Vec pos1){

	int mask0 = (pos0.getX() < -100.0) ? 1 : 0;
	int mask1 = (pos1.getX() < -100.0) ? 1 : 0;

	float dlon = pos0.getX()-pos1.getX();
	float lat0 = lat_mu(pos0.getY());
	float lat1 = lat_mu(pos1.getY());
	float dlat = lat1-lat0;

	float a = sin(dlat/2)*sin(dlat/2)+cos(lat0)*cos(lat1)*sin(dlon/2)*sin(dlon/2);

	return(2.0*atan2(sqrt(a),sqrt(1-a))/M_PI*180.0*(1-mask0)*(1-mask1)+
			(-999.0)*((1-mask0)*mask1+(1-mask1)*mask0+mask1*mask0));

}

float Grid::euclidean(Vec pos0,Vec pos1){

	int mask0 = (pos0.getX() < -100.0) ? 1 : 0;
	int mask1 = (pos1.getX() < -100.0) ? 1 : 0;

	float lat0 = lat_mu(pos0.getY());
	float lat1 = lat_mu(pos1.getY());

	float d = sqrt(pow(pos0.getX()-pos1.getX(),2)+pow(lat0-lat1,2));

	return(d/M_PI*180.0*(1-mask0)*(1-mask1)+
			(-999.0)*((1-mask0)*mask1+(1-mask1)*mask0+mask1*mask0));

}
#endif

#ifdef CIRCULAR
void Grid::write_simulation(std::string w,double dt_init,double dt_sim){

	netCDF::NcFile data(w+".nc", netCDF::NcFile::replace);

	data.putAtt("title","Lagrangian simulation Northern Atlantic Ocean");
	time_t timestamp;
	time(&timestamp);
	data.putAtt("clock time",ctime(&timestamp));

	Particle p = Particle(this->pos0.getX(),this->pos0.getY(),0);
	data.putAtt("central starting position","("+
					std::to_string(p.getPos().getX()/M_PI*180)+" ; "+
					std::to_string(lat_mu(p.getPos().getY())/M_PI*180)+") degrees");
	#ifdef CIRCULAR
		data.putAtt("initial condition","circular");
		data.putAtt("radius",std::to_string(this->radius)+" km");
	#endif

	#if defined(STOREPOS) || defined(STOREVEL)
		auto t_start = std::chrono::high_resolution_clock::now();
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
		float vec_initx[calc_ndays(NYEARSTART+YSTART)/DTSTART][NPART];
		float vec_inity[calc_ndays(NYEARSTART+YSTART)/DTSTART][NPART];
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

#endif

#ifdef LYAPUNOV

void Grid::write_simulation(std::string w,double dt_init,double dt_sim){

	int nlon = (int)((LYAPLONMAX-LYAPLONMIN)/LYAPLONRES);
	int nlat = (int)((LYAPLATMAX-LYAPLATMIN)/LYAPLATRES);

	netCDF::NcFile data(w+".nc", netCDF::NcFile::replace);

	data.putAtt("title","Lyapunov Exponents Northern Atlantic Ocean");
	time_t timestamp;
	time(&timestamp);
	data.putAtt("clock time",ctime(&timestamp));

	auto t_start = std::chrono::high_resolution_clock::now();
	data.putAtt("simulation type","Lyapunov");

	netCDF::NcDim lonDim = data.addDim("lon", (nlon-4));
	netCDF::NcDim latDim = data.addDim("lat", (nlat-4));

	std::vector<netCDF::NcDim> dimVector;
	dimVector.push_back(latDim);
	dimVector.push_back(lonDim);

	std::vector<netCDF::NcDim> dimVector_lon;
	dimVector_lon.push_back(lonDim);
	std::vector<netCDF::NcDim> dimVector_lat;
	dimVector_lat.push_back(latDim);

	netCDF::NcVar lonVar = data.addVar("lon", netCDF::ncFloat, dimVector_lon);
	lonVar.putAtt("units", "degrees");
	netCDF::NcVar latVar = data.addVar("lat", netCDF::ncFloat, dimVector_lat);
	latVar.putAtt("units", "degrees");
	float vec_lon[(nlon-4)];
	float vec_lat[(nlat-4)];
	
	for(int j=2;j<(nlon-2);j++){
		vec_lon[j-2] = (LYAPLONMIN+j*LYAPLONRES)/M_PI*180;
	}
	for(int j=2;j<(nlat-2);j++){
		vec_lat[j-2] = (LYAPLATMIN+j*LYAPLATRES)/M_PI*180;
	}

	std::vector<size_t> startp_lon,countp_lon;
	startp_lon.push_back(0);
	countp_lon.push_back((nlon-4));
	std::vector<size_t> startp_lat,countp_lat;
	startp_lat.push_back(0);
	countp_lat.push_back((nlat-4));
	lonVar.putVar(startp_lon,countp_lon,vec_lon);
	latVar.putVar(startp_lat,countp_lat,vec_lat);

	std::vector<size_t> startp,countp;
	startp.push_back(0);
	startp.push_back(0);
	countp.push_back(nlat-4);
	countp.push_back(nlon-4);

	netCDF::NcVar distVar = data.addVar("tau", netCDF::ncInt, dimVector);
	distVar.putAtt("units", "hours");
	int mat_tau[(nlat-4)][(nlon-4)];

	#pragma omp parallel for
	for(int ilat=2;ilat<(nlat-2);ilat++){
		for(int ilon=2;ilon<(nlon-2);ilon++){

				float vec_dist[4];
				int mask = 0;
				float mean_dist = 0;
				int c = 0;

				mat_tau[ilat-2][ilon-2] = -999;

				for(int t=NMONTH*28*24;t >= 0;t--){

					vec_dist[0] = euclidean(this->particles[(ilon-1)+(nlon-2)*(ilat-1)].getPathPos()[t],
											this->particles[(ilon-2)+(nlon-2)*(ilat-1)].getPathPos()[t]);
					vec_dist[1] = euclidean(this->particles[(ilon-1)+(nlon-2)*(ilat-1)].getPathPos()[t],
											this->particles[(ilon-1)+(nlon-2)*(ilat)].getPathPos()[t]);
					vec_dist[2] = euclidean(this->particles[(ilon-1)+(nlon-2)*(ilat-1)].getPathPos()[t],
											this->particles[(ilon)+(nlon-2)*(ilat-1)].getPathPos()[t]);
					vec_dist[3] = euclidean(this->particles[(ilon-1)+(nlon-2)*(ilat-1)].getPathPos()[t],
											this->particles[(ilon-1)+(nlon-2)*(ilat-2)].getPathPos()[t]);

					for(int x=0;x<4;x++){
						mask = (vec_dist[x] < -100.0) ? 1 : 0;
						mean_dist += vec_dist[x]*(1-mask);
						c += (1-mask);
					}
					
					mask = (c == 0) ? 1 : 0;
					c = (c == 0) ? 1 : c;

					if(mask == 1){
						mat_tau[ilat-2][ilon-2] = -999;
						mean_dist = 0.0;
						c = 0;
						break;
					}else if(mean_dist/c > DEND){
						mat_tau[ilat-2][ilon-2] = 28*24*NMONTH-t;
						mean_dist = 0.0;
						c = 0;
						break;
					}

					mean_dist = 0.0;
					c = 0;

				}

		}
	}	

				
	distVar.putVar(startp,countp,mat_tau);
				
	data.putAtt("length of trajectories",std::to_string(NMONTH)+" months");
	data.putAtt("number of particles per release",std::to_string((nlon-2)*(nlat-2)));
	data.putAtt("starting date","31-03-2026");
	data.putAtt("diffusion constant",std::to_string(D)+" m^2/s");
	data.putAtt("final distance",std::to_string(DEND)+" degrees");
	auto t_end = std::chrono::high_resolution_clock::now();
	double dt_writ = std::chrono::duration<double, std::milli>(t_end-t_start).count();
	data.putAtt("initialisation wall time",std::to_string(dt_init/1000)+" seconds");
	data.putAtt("simulation wall time",std::to_string(dt_sim/1000)+" seconds");
	data.putAtt("output wall time",std::to_string(dt_writ/1000)+" seconds");

}

#endif

#ifdef NETWORK
void Grid::write_simulation(std::string w,double dt_init,double dt_sim){

	auto t_start = std::chrono::high_resolution_clock::now();
	netCDF::NcFile data(w+".nc", netCDF::NcFile::replace);

	data.putAtt("title","Lagrangian network simulation Northern Atlantic Ocean");
	time_t timestamp;
	time(&timestamp);
	data.putAtt("clock time",ctime(&timestamp));
	Particle p = Particle(this->pos0.getX(),this->pos0.getY(),0);
	p.xy_to_lonmu();
	//data.putAtt("cell ID",)
	data.putAtt("cell central position","("+
					std::to_string(p.getPos().getX()/M_PI*180)+" ; "+
					std::to_string(lat_mu(p.getPos().getY())/M_PI*180)+") degrees");
	
	data.putAtt("initial condition","healpix");

	data.putAtt("simulation type","network");

	netCDF::NcDim startDim = data.addDim("start", this->Nstart);
	netCDF::NcDim partDim = data.addDim("particles", NPART);
	netCDF::NcDim cellDim = data.addDim("cells", NCELL);

	std::vector<netCDF::NcDim> dimVector;
	dimVector.push_back(cellDim);
	dimVector.push_back(startDim);
	dimVector.push_back(partDim);

	std::vector<netCDF::NcDim> dimVector_start;
	dimVector_start.push_back(startDim);

	std::vector<size_t> startp,countp;
	startp.push_back(0);
	startp.push_back(0);
	startp.push_back(0);
	countp.push_back(1);
	countp.push_back(this->Nstart);
	countp.push_back(NPART);

	std::vector<size_t> startp_start,countp_start;
	startp_start.push_back(0);
	countp_start.push_back(this->Nstart);
	netCDF::NcVar startVar = data.addVar("start", netCDF::ncInt, dimVector_start);
	startVar.putAtt("units", "days");
	int vec_start[this->Nstart];
	for(int i=0;i < this->Nstart;i++){
		vec_start[i] = i*DTSTART;
	}
	startVar.putVar(startp_start,countp_start,vec_start);

	netCDF::NcVar tminVar = data.addVar("tmin", netCDF::ncInt, dimVector);
	tminVar.putAtt("units", "days");
	int vec_part_tmin[this->Nstart][NPART];

	//#pragma omp parallel for num_threads(4)
	for(size_t c=0;c < NCELL;c++){
		
		startp[0]=c;

			for(int s=0;s<this->Nstart;s++){

				for(size_t k=0;k < NPART;k++){

					vec_part_tmin[s][k] = this->network[c+NCELL*(k+s*NPART)];
					
				}				
			}

		
		tminVar.putVar(startp,countp,vec_part_tmin);
					
	}

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


}
#endif

#ifdef SST

void Grid::write_simulation(std::string w,double dt_init,double dt_sim){

	int nlon = (int)((SSTLONMAX-SSTLONMIN)/SSTLONRES);
	int nlat = (int)((SSTLATMAX-SSTLATMIN)/SSTLATRES);

	netCDF::NcFile data(w+".nc", netCDF::NcFile::replace);

	data.putAtt("title","SST origins Northern Atlantic Ocean");
	time_t timestamp;
	time(&timestamp);
	data.putAtt("clock time",ctime(&timestamp));

	auto t_start = std::chrono::high_resolution_clock::now();
	data.putAtt("simulation type","SST");

	netCDF::NcDim lonDim = data.addDim("lon", (nlon-2));
	netCDF::NcDim latDim = data.addDim("lat", (nlat-2));

	std::vector<netCDF::NcDim> dimVector;
	dimVector.push_back(latDim);
	dimVector.push_back(lonDim);

	std::vector<netCDF::NcDim> dimVector_lon;
	dimVector_lon.push_back(lonDim);
	std::vector<netCDF::NcDim> dimVector_lat;
	dimVector_lat.push_back(latDim);

	netCDF::NcVar lonVar = data.addVar("lon", netCDF::ncFloat, dimVector_lon);
	lonVar.putAtt("units", "degrees");
	netCDF::NcVar latVar = data.addVar("lat", netCDF::ncFloat, dimVector_lat);
	latVar.putAtt("units", "degrees");
	float vec_lon[(nlon-2)];
	float vec_lat[(nlat-2)];
	
	for(int j=1;j<(nlon-1);j++){
		vec_lon[j-1] = (SSTLONMIN+j*SSTLONRES)/M_PI*180;
	}
	for(int j=1;j<(nlat-1);j++){
		vec_lat[j-1] = (SSTLATMIN+j*SSTLATRES)/M_PI*180;
	}

	std::vector<size_t> startp_lon,countp_lon;
	startp_lon.push_back(0);
	countp_lon.push_back((nlon-2));
	std::vector<size_t> startp_lat,countp_lat;
	startp_lat.push_back(0);
	countp_lat.push_back((nlat-2));
	lonVar.putVar(startp_lon,countp_lon,vec_lon);
	latVar.putVar(startp_lat,countp_lat,vec_lat);

	std::vector<size_t> startp,countp;
	startp.push_back(0);
	startp.push_back(0);
	countp.push_back(nlat-2);
	countp.push_back(nlon-2);

	netCDF::NcVar meandistVar = data.addVar("meandistance", netCDF::ncFloat, dimVector);
	meandistVar.putAtt("units", "kilometres");
	netCDF::NcVar vardistVar = data.addVar("vardistance", netCDF::ncFloat, dimVector);
	vardistVar.putAtt("units", "kilometres squared");
	netCDF::NcVar meanSSTVar = data.addVar("meanSST", netCDF::ncFloat, dimVector);
	meanSSTVar.putAtt("units", "degrees celcius");
	netCDF::NcVar meanSSTdiffVar = data.addVar("meanSSTdiff", netCDF::ncFloat, dimVector);
	meanSSTVar.putAtt("units", "degrees celcius");
	netCDF::NcVar varSSTVar = data.addVar("varSST", netCDF::ncFloat, dimVector);
	varSSTVar.putAtt("units", "degrees celcius squared");
	float mat_meandist[(nlat-2)][(nlon-2)];
	float mat_vardist[(nlat-2)][(nlon-2)];
	float mat_meanSST[(nlat-2)][(nlon-2)];
	float mat_varSST[(nlat-2)][(nlon-2)];
	float mat_meanSSTdiff[(nlat-2)][(nlon-2)];

	#pragma omp parallel for
	for(int ilat=1;ilat<(nlat-1);ilat++){
		for(int ilon=1;ilon<(nlon-1);ilon++){

				Particle p0 = Particle();
				p0.get_initial_pos(Vec(SSTLONMIN+ilon*SSTLONRES,SSTLATMIN+ilat*SSTLATRES),0,0,0,0);
				float temp_0 = p0.interpol(this->SSTs);

				float s_dist = 0.0;
				float s2_dist = 0.0;
				float s_sst = 0.0;
				float s_shift_sst = 0.0;
				float s2_sst = 0.0;
				float sdiff_sst = 0.0;
				int mask_dist = 0;
				int c_dist = 0;
				int mask_sst = 0;
				int c_sst = 0;
				float shift_sst = this->particles[0+NPART*(ilon-1+(nlon-2)*(ilat-1))].interpol(this->SSTs);
				Vec pos0 = Vec(SSTLONMIN+ilon*SSTLONRES,mu_lat(SSTLATMIN+ilat*SSTLATRES));
				float dist_j,temp_j;

				for(int j=0;j<NPART;j++){

					dist_j = haversine(pos0,this->particles[j+NPART*(ilon-1+(nlon-2)*(ilat-1))].getPos());
					mask_dist = (dist_j < -100.0) ? 1 : 0;
					s_dist += (1-mask_dist)*dist_j;
					s2_dist += (1-mask_dist)*pow(dist_j,2);
					c_dist += (1-mask_dist); 

					temp_j = this->particles[j+NPART*(ilon-1+(nlon-2)*(ilat-1))].interpol(this->SSTs);
					mask_sst = (temp_j < -100.0) ? 1 : 0;
					s_sst += (1-mask_sst)*temp_j;
					s_shift_sst += (1-mask_sst)*(temp_j-shift_sst);
					s2_sst += (1-mask_sst)*pow(temp_j-shift_sst,2);
					sdiff_sst += (1-mask_sst)*(temp_0-temp_j);
					c_sst += (1-mask_sst); 

				}

				mask_dist = (c_dist < (int)(NPART/10)) ? 1 : 0;
				c_dist = (mask_dist == 1) ? 1 : c_dist;
				mat_meandist[ilat-1][ilon-1] = (1-mask_dist)*s_dist/c_dist/180.0*M_PI*R/1000.0 + (-999.0)*mask_dist;
				mat_vardist[ilat-1][ilon-1] = (1-mask_dist)*(s2_dist/c_dist-pow(s_dist/c_dist,2))/180.0/180.0*M_PI*M_PI*R*R/1000.0/1000.0 + (-999.0)*mask_dist;

				mask_sst = (c_sst < (int)(NPART/10)) ? 1 : 0;
				c_sst = (mask_sst == 1) ? 1 : c_sst;
				mat_meanSST[ilat-1][ilon-1] = (1-mask_sst)*s_sst/c_sst + (-999.0)*mask_sst;
				mat_varSST[ilat-1][ilon-1] = (1-mask_sst)*(s2_sst/c_sst-pow(s_shift_sst/c_sst,2)) + (-999.0)*mask_sst;
				mat_meanSSTdiff[ilat-1][ilon-1] = (1-mask_sst)*sdiff_sst/c_sst + (-999.0)*mask_sst;

		}
	}	

				
	meandistVar.putVar(startp,countp,mat_meandist);
	vardistVar.putVar(startp,countp,mat_vardist);
	meanSSTVar.putVar(startp,countp,mat_meanSST);
	varSSTVar.putVar(startp,countp,mat_varSST);
	meanSSTdiffVar.putVar(startp,countp,mat_meanSSTdiff);
				
	data.putAtt("length of trajectories",std::to_string(NMONTH)+" months");
	data.putAtt("number of grid points",std::to_string((nlon-2)*(nlat-2)));
	data.putAtt("number of particles per release",std::to_string(NPART));
	data.putAtt("starting date","31-03-2026");
	data.putAtt("diffusion constant",std::to_string(D)+" m^2/s");
	auto t_end = std::chrono::high_resolution_clock::now();
	double dt_writ = std::chrono::duration<double, std::milli>(t_end-t_start).count();
	data.putAtt("initialisation wall time",std::to_string(dt_init/1000)+" seconds");
	data.putAtt("simulation wall time",std::to_string(dt_sim/1000)+" seconds");
	data.putAtt("output wall time",std::to_string(dt_writ/1000)+" seconds");

}

#endif