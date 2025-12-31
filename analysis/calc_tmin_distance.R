library(tidyverse)
library(ncdf4)
library(sf)
sf_use_s2(FALSE)
library(sp)

lat_mu <- function(mu){

	return(pi/2 - 2*atan(exp(-mu)))

}

crs_atl <- CRS("+proj=sinu +lon_0=0 +a=6378140.0 +b=6356750.0 +units=m +no_defs")

sarg <- st_read("sargasso_shapefile/WC_13_EBSA.shp") %>%
			st_geometry() %>% st_convex_hull() %>% 
			st_transform(crs=crs_atl) %>%
			st_buffer(-200000) %>% st_transform(crs=4326)

dists <- seq(500,1000,50)

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

pos0 <- nc_sim %>% ncatt_get(0)
pos0 <- pos0$'central starting position'
pos0 <- pos0 %>% substr(2,nchar(.)-9)
pos0 <- data.frame(x=unlist(strsplit(pos0,split=";"))[1] %>% as.numeric(),
					y=unlist(strsplit(pos0,split=";"))[2] %>% as.numeric()) %>%
			st_as_sf(crs=4326,coords=c("x","y"))

Nsim <- dim(lon_mat)[4]
Npart <- dim(lon_mat)[1]
Nyear <- dim(lon_mat)[3]

nc_sim %>% nc_close()

mat_dists <- array(NA, c(Npart,length(dists),Nsim))

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

#sf_pos <- sf_pos %>% slice(which((st_intersects(sf_pos,sarg,sparse=FALSE) %>% as.numeric()) == 1))

sf_pos <- sf_pos %>% mutate(distance=(st_distance(.,pos0) %>% as.numeric())/1e3)

for(d in 1:length(dists)){

df_dist <- sf_pos %>% st_drop_geometry() %>% filter(distance > dists[d]) %>%
			group_by(ID) %>% summarise(tmin=min(time))

if(nrow(df_dist) != Npart){
df_dist <- rbind(df_dist,
		data.frame(ID=which(!(1:Npart %in% (df_dist %>% pull(ID)))),tmin=NA)) %>%
		arrange(ID)
}

mat_dists[,d,sim] <- df_dist %>% pull(tmin)

}

}

dir_dist <- "/media/kobe/Windows/spectrum/distance/"
partdim <- ncdim_def("part","paricle ID",1:Npart) 
startdim <- ncdim_def("start","starting date",1:Nsim) 
distdim <- ncdim_def("distance","distance (km)",dists)

tminvar <- ncvar_def("tmin","days",list(partdim,distdim,startdim))

ncout <- nc_create(paste0(dir_dist,list_sim[sim_i]),list(tminvar),force_v4=TRUE)

ncout %>% ncvar_put(.,tminvar,mat_dists)

for(i in 1:13){

	ncout %>% ncatt_put(.,0,names(atts)[i],atts[[i]])

}

ncout %>% nc_close()
}