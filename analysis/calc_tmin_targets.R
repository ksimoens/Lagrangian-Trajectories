library(tidyverse)
library(ncdf4)
library(sf)
sf_use_s2(FALSE)
library(sp)

lat_mu <- function(mu){

	return(pi/2 - 2*atan(exp(-mu)))

}

make_targets_dist <- function(dist,R,pos0,sarg,crs_atl){

	N <- floor(2*pi/1.25/(2*R/dist))
	
	coord0 <- (pos0 %>% st_transform(crs=crs_atl) %>%
				st_coordinates() %>% as.numeric())/1e3
	
	sf_targ <- matrix(ncol=2,nrow=N) %>% as.data.frame() %>%
				rename(x=V1,y=V2)

	for(i in 0:(N-1)){
		sf_targ$x[i+1] <- coord0[1] + dist*cos(i*2*pi/N)
		sf_targ$y[i+1] <- coord0[2] + dist*sin(i*2*pi/N)
	}

	sf_targ <- (1e3*sf_targ) %>% st_as_sf(crs=crs_atl,coords=c("x","y")) %>%
					st_buffer(R*1e3)

	Atot <- st_area(sf_targ[1,]) %>% as.numeric()

	i_inter <- which(as.numeric(st_intersects(sf_targ,sarg,sparse=FALSE)) != 0)

	sf_targ <- sf_targ %>% slice(i_inter)

	i_in <- which((Atot - as.numeric(st_area(st_intersection(sf_targ,sarg)))) < 1)

	sf_targ <- sf_targ %>% slice(i_in)

	return(sf_targ)

}

make_targets <- function(dists,R,pos0,sarg,crs_atl){

	sf_targ <- make_targets_dist(dists[1],R,pos0,sarg,crs_atl) %>%
				mutate(distance=dists[1])
	for(i in 2:length(dists)){
		sf_targ <- rbind(sf_targ,
						make_targets_dist(dists[i],R,pos0,sarg,crs_atl)	%>%
							mutate(distance=dists[i])
						)
	}

	return(sf_targ)

}

crs_atl <- CRS("+proj=sinu +lon_0=0 +a=6378140.0 +b=6356750.0 +units=m +no_defs")

sarg <- st_read("sargasso_shapefile/WC_13_EBSA.shp") %>%
			st_geometry() %>% st_convex_hull() %>% 
			st_transform(crs=crs_atl) %>%
			st_buffer(-200000) 

dists <- seq(500,1000,100)
R <- sqrt(100*100/pi)

dir <- "/media/kobe/Windows/spectrum/output/"

list_sim <- list.files(dir)

for(sim_i in 1:length(list_sim)){

	print(sim_i)

nc_sim <- nc_open(paste0(dir,list_sim[sim_i]))

lon_mat <- nc_sim %>% ncvar_get("lon")
lon_mat[lon_mat==-999] <- NA
mu_mat <- nc_sim %>% ncvar_get("mu")
mu_mat[mu_mat==-999] <- NA
#u_mat <- nc_sim %>% ncvar_get("u")
#u_mat[abs(u_mat) > 5] <- NA
#v_mat <- nc_sim %>% ncvar_get("v")
#v_mat[abs(v_mat) > 5] <- NA

atts <- nc_sim %>% ncatt_get(0)
pos0 <- atts$'central starting position'
pos0 <- pos0 %>% substr(2,nchar(.)-9)
pos0 <- data.frame(x=unlist(strsplit(pos0,split=";"))[1] %>% as.numeric(),
					y=unlist(strsplit(pos0,split=";"))[2] %>% as.numeric()) %>%
			st_as_sf(crs=4326,coords=c("x","y"))

sf_targ <- make_targets(dists,R,pos0,sarg,crs_atl)

Nsim <- dim(lon_mat)[4]
Npart <- dim(lon_mat)[1]
Nyear <- dim(lon_mat)[3]
Ntarg <- nrow(sf_targ)

nc_sim %>% nc_close()

mat_dists <- array(NA, c(Npart,Ntarg,Nsim))

for(sim in 1:Nsim){
	print(sim)

year <- 1

sf_pos <- lon_mat[,,year,sim] %>% as.data.frame() %>% 
			mutate(ID=1:nrow(.)) %>%
			pivot_longer(1:(ncol(.)-1),names_to="time",values_to="lon") %>%
			mutate(time=sub(".","",time) %>% as.numeric()) %>%
			mutate(time=time+(year-1)*365) %>%
			mutate(lon=lon/pi*180) %>%
			mutate(lat=mu_mat[,,year,sim] %>% as.data.frame() %>% 
				mutate(ID=(sim-1)*Npart+1:nrow(.)) %>%
				pivot_longer(1:(ncol(.)-1),names_to="time",values_to="lat") %>%
				mutate(lat=lat_mu(lat)/pi*180) %>% pull(lat))

for(year in 2:Nyear){

sf_pos <- rbind(sf_pos,
			lon_mat[,,year,sim] %>% as.data.frame() %>% 
			mutate(ID=1:nrow(.)) %>%
			pivot_longer(1:(ncol(.)-1),names_to="time",values_to="lon") %>%
			mutate(time=sub(".","",time) %>% as.numeric()) %>%
			mutate(time=time+(year-1)*365) %>%
			mutate(lon=lon/pi*180) %>%
			mutate(lat=mu_mat[,,year,sim] %>% as.data.frame() %>% 
				mutate(ID=(sim-1)*Npart+1:nrow(.)) %>%
				pivot_longer(1:(ncol(.)-1),names_to="time",values_to="lat") %>%
				mutate(lat=lat_mu(lat)/pi*180) %>% pull(lat))
			)

}

sf_pos <- sf_pos %>% filter(!is.na(lon)) %>% 
			st_as_sf(crs=4326,coords=c("lon","lat"))

inter <- st_intersects(sf_targ %>% st_transform(crs=4326),
						sf_pos,sparse=TRUE) %>% as.data.frame() %>%
						mutate(time=sf_pos$time[col.id],
								ID=sf_pos$ID[col.id]) %>%
						group_by(row.id,ID) %>%
						summarise(tmin=min(time)) %>%
						pivot_wider(names_from=row.id,values_from=tmin)

if(nrow(inter) != Npart){
	ID_miss <- which(!(1:1000 %in% inter$ID))
	mat_miss <- matrix(ncol=ncol(inter),nrow=length(ID_miss))
	mat_miss[,1] <- ID_miss
	df_miss <- mat_miss %>% as.data.frame() 
	names(df_miss) <- names(inter)

	inter <- rbind(inter,df_miss) %>% arrange(ID) 
			
}

inter <- inter %>% select(-c(ID)) %>% as.matrix()

mat_dists[,,sim] <- inter

}

dir_targ <- "/media/kobe/Windows/spectrum/target/"
partdim <- ncdim_def("part","paricle ID",1:Npart) 
startdim <- ncdim_def("start","starting date",1:Nsim) 
targdim <- ncdim_def("target","target ID",1:nrow(sf_targ))

tminvar <- ncvar_def("tmin","days",list(partdim,targdim,startdim))
lonvar <- ncvar_def("target longitude","degrees east",list(targdim))
latvar <- ncvar_def("target latitude","degrees north",list(targdim))
distvar <- ncvar_def("target distance","kilometres",list(targdim))

ncout <- nc_create(paste0(dir_targ,list_sim[sim_i]),
					list(tminvar,lonvar,latvar,distvar),
					force_v4=TRUE)

df_targ <- cbind(sf_targ %>% st_drop_geometry(), 
				sf_targ %>% st_centroid() %>% st_transform(crs=4326) %>%
					st_coordinates() %>% as.data.frame())

ncout %>% ncvar_put(.,tminvar,mat_dists)
ncout %>% ncvar_put(.,lonvar,df_targ %>% pull(X))
ncout %>% ncvar_put(.,latvar,df_targ %>% pull(Y))
ncout %>% ncvar_put(.,distvar,df_targ %>% pull(distance))

for(i in 1:13){

	ncout %>% ncatt_put(.,0,names(atts)[i],atts[[i]])

}

ncout %>% nc_close()

}