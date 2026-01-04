library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(ncdf4)

xy_to_lon <- function(x,y){

	lon <- rep(NA,length(x))
	lon[abs(y) < pi/4] <- x[abs(y) < pi/4]
	lon[abs(y) >= pi/4] <- x[abs(y) >= pi/4] - (abs(y[abs(y) >= pi/4])-pi/4)/(abs(y[abs(y) >= pi/4])-pi/2)*(x[abs(y) >= pi/4]%%(pi/2)-pi/4)

	return((lon-pi)/pi*180)

}

xy_to_lat <- function(y){

	lat <- rep(NA,length(y))
	lat[abs(y) < pi/4] <- 8/3/pi*y[abs(y) < pi/4]
	lat[abs(y) >= pi/4] <- (1-1/3*(2-4*abs(y[abs(y) >= pi/4])/pi)^2)*sign(y[abs(y) >= pi/4])

	return((pi/2-acos(lat))/pi*180)

}

get_hor_lines <- function(df_line,sf_line,vec_lon,vec_lat){


lon_line <- vec_lon[vec_lon > min(df_line$lon)-0.5 & 
								vec_lon < max(df_line$lon)+0.5]
lat_line <- vec_lat[vec_lat > min(df_line$lat)-0.5 & 
								vec_lat < max(df_line$lat)+0.5]

df_vel_line <- expand.grid(lon_line,lat_line) %>%
				rename(lon=1,lat=2)

sf_vel_line <- df_vel_line %>% st_as_sf(crs=4326,coords=c("lon","lat"))

hor_lines <- data.frame(x=c(lon_line[1],lon_line[2]),
						y=c(lat_line[1],lat_line[1])) %>%
				mutate(ID=1) %>% st_as_sf(crs=4326,coords=c("x","y")) %>%
				group_by(ID) %>% summarise(do_union=FALSE) %>%
				st_cast("LINESTRING")

k <- 1
for(ilat in 1:length(lat_line)){

	for(ilon in 1:(length(lon_line)-1)){

		hor_lines <- rbind(hor_lines,
				data.frame(x=c(lon_line[ilon],lon_line[ilon+1]),
						y=c(lat_line[ilat],lat_line[ilat])) %>%
				mutate(ID=k) %>% st_as_sf(crs=4326,coords=c("x","y")) %>%
				group_by(ID) %>% summarise(do_union=FALSE) %>%
				st_cast("LINESTRING")
				)
		k <- k + 1

	}

}

hor_lines <- hor_lines %>% slice(-1) %>%
				slice(which(st_intersects(.,sf_line,sparse=FALSE)))

	return(hor_lines)

}

get_ver_lines <- function(df_line,sf_line,vec_lon,vec_lat){


lon_line <- vec_lon[vec_lon > min(df_line$lon)-0.5 & 
								vec_lon < max(df_line$lon)+0.5]
lat_line <- vec_lat[vec_lat > min(df_line$lat)-0.5 & 
								vec_lat < max(df_line$lat)+0.5]

df_vel_line <- expand.grid(lon_line,lat_line) %>%
				rename(lon=1,lat=2)

sf_vel_line <- df_vel_line %>% st_as_sf(crs=4326,coords=c("lon","lat"))

ver_lines <- data.frame(x=c(lon_line[1],lon_line[1]),
						y=c(lat_line[1],lat_line[2])) %>%
				mutate(ID=1) %>% st_as_sf(crs=4326,coords=c("x","y")) %>%
				group_by(ID) %>% summarise(do_union=FALSE) %>%
				st_cast("LINESTRING")

k <- 1
for(ilon in 1:length(lon_line)){

	for(ilat in 1:(length(lat_line)-1)){

		ver_lines <- rbind(ver_lines,
				data.frame(x=c(lon_line[ilon],lon_line[ilon]),
						y=c(lat_line[ilat],lat_line[ilat+1])) %>%
				mutate(ID=k) %>% st_as_sf(crs=4326,coords=c("x","y")) %>%
				group_by(ID) %>% summarise(do_union=FALSE) %>%
				st_cast("LINESTRING")
				)
		k <- k + 1

	}

}

ver_lines <- ver_lines %>% slice(-1) %>%
				slice(which(st_intersects(.,sf_line,sparse=FALSE)))

	return(ver_lines)

}

get_perp_velocity_hor <- function(coords_i,u_mat,v_mat,vec_lon,vec_lat,angle){

	i_lat <- which(vec_lat == coords_i[2])
i_lon <- max(which(vec_lon < coords_i[1]))

u_inter <- (u_mat[i_lon,i_lat]*(vec_lon[i_lon+1]-coords_i[1]) + 
			u_mat[i_lon+1,i_lat]*(coords_i[1]-vec_lon[i_lon]))/0.25

v_inter <- (v_mat[i_lon,i_lat]*(vec_lon[i_lon+1]-coords_i[1]) + 
			v_mat[i_lon+1,i_lat]*(coords_i[1]-vec_lon[i_lon]))/0.25

vel_perp <- u_inter*cos(pi/2-angle) + v_inter*cos(pi-angle)

return(vel_perp)

}

get_perp_velocity_ver <- function(coords_i,u_mat,v_mat,vec_lon,vec_lat,angle){

	i_lon <- which(vec_lon == coords_i[1])
i_lat <- max(which(vec_lat < coords_i[2]))

u_inter <- (u_mat[i_lon,i_lat]*(vec_lat[i_lat+1]-coords_i[2]) + 
			u_mat[i_lon,i_lat+1]*(coords_i[2]-vec_lat[i_lat]))/0.25

v_inter <- (v_mat[i_lon,i_lat]*(vec_lat[i_lat+1]-coords_i[2]) + 
			v_mat[i_lon,i_lat+1]*(coords_i[2]-vec_lat[i_lat]))/0.25

vel_perp <- u_inter*cos(pi/2-angle) + v_inter*cos(pi-angle)

return(vel_perp)

}

get_hor_perp_velocity <- function(df_line,sf_line,vec_lon,vec_lat,u_mat,v_mat,angle){

	hor_lines <- get_hor_lines(df_line,sf_line,vec_lon,vec_lat)

pts_inter <- st_intersection(hor_lines,sf_line) %>% st_coordinates() %>%
				as.data.frame()

vec_vel_perp <- c()

for(i in 1:nrow(pts_inter)){

coords_i <- pts_inter %>% slice(i) %>% as.numeric()
vec_vel_perp <- c(vec_vel_perp,get_perp_velocity_hor(coords_i,u_mat,v_mat,vec_lon,vec_lat,angle))

}

return(vec_vel_perp)

}

get_ver_perp_velocity <- function(df_line,sf_line,vec_lon,vec_lat,u_mat,v_mat,angle){

	ver_lines <- get_ver_lines(df_line,sf_line,vec_lon,vec_lat)

pts_inter <- st_intersection(ver_lines,sf_line) %>% st_coordinates() %>%
				as.data.frame()

vec_vel_perp <- c()

for(i in 1:nrow(pts_inter)){

coords_i <- pts_inter %>% slice(i) %>% as.numeric()
vec_vel_perp <- c(vec_vel_perp,get_perp_velocity_ver(coords_i,u_mat,v_mat,vec_lon,vec_lat,angle))

}

return(vec_vel_perp)

}

get_angle <- function(df_line){

	dx <- df_line$lon[2] - df_line$lon[1]
	dy <- df_line$lat[2] - df_line$lat[1] 

	if(dx < 0){
		return(pi+atan(dy/dx))
	} else{
		return(atan(dy/dx))
	}

}

get_line_flow <- function(df_line,vec_lon,vec_lat,u_mat,v_mat){

	sf_line <- df_line %>% st_as_sf(crs=4326,coords=c("lon","lat")) %>%
			mutate(ID=1) %>% group_by(ID) %>% summarise(do_union=FALSE) %>%
			st_cast("LINESTRING")

angle <- get_angle(df_line)

vec_vel_perp <- c(get_hor_perp_velocity(df_line,sf_line,vec_lon,vec_lat,u_mat,v_mat,angle),
					 get_ver_perp_velocity(df_line,sf_line,vec_lon,vec_lat,u_mat,v_mat,angle))

return( abs(as.numeric(st_distance(df_line %>% slice(1) %>% st_as_sf(coords=c("lon","lat"),crs=4326),
					df_line %>% slice(2) %>% st_as_sf(coords=c("lon","lat"),crs=4326)))*mean(vec_vel_perp,na.rm=TRUE)))

}

Nside <- 64
ds <- pi/Nside/2

df_grid <- read.csv("../../healpix/test_grid.csv")
sf_grid <- st_read("../../healpix/test_grid/test_grid.shp")

nc_vel <- nc_open("vel_average.nc")

vec_lon <- nc_vel %>% ncvar_get("lon")
vec_lat <- nc_vel %>% ncvar_get("lat")
u_mat <- nc_vel %>% ncvar_get("u")
v_mat <- nc_vel %>% ncvar_get("v")

nc_vel %>% nc_close()

df_grid <- df_grid %>% mutate(flow=NA)

for(i in 1:nrow(df_grid)){
	print(i)
	cell <- df_grid %>% slice(i)
	sf_cell <- data.frame(x=c(cell$x-ds/2,cell$x,cell$x+ds/2,cell$x),
							y=c(cell$y,cell$y+ds/2,cell$y,cell$y-ds/2)) %>%
				mutate(lon=xy_to_lon(x,y),lat=xy_to_lat(y))

	flow_i <- get_line_flow(sf_cell %>% slice(1:2),vec_lon,vec_lat,u_mat,v_mat)
	flow_i <- flow_i + get_line_flow(sf_cell %>% slice(2:3),vec_lon,vec_lat,u_mat,v_mat)
	flow_i <- flow_i + get_line_flow(sf_cell %>% slice(3:4),vec_lon,vec_lat,u_mat,v_mat)
	flow_i <- flow_i + get_line_flow(sf_cell %>% slice(c(4,1)),vec_lon,vec_lat,u_mat,v_mat)

	df_grid$flow[i] <- flow_i

}

sf_grid <- sf_grid %>% mutate(flow=df_grid %>% pull(flow))

sf_land <- st_read("../../../land/ne_50m_land/ne_50m_land.shp") %>%
			st_crop(data.frame(x=c(-80,-80,-30,-30),y=c(25,70,70,25)) %>%
						st_as_sf(crs=4326,coords=c("x","y")))

p <- ggplot() +
		geom_sf(data=sf_land) +
		geom_sf(data=sf_grid,aes(col=flow/1e4,fill=flow/1e4)) +
		scale_x_continuous(expand=c(0,0),limits=c(-80,-25)) +
		scale_y_continuous(expand=c(0,0),limits=c(27.5,50)) +
		scale_colour_viridis_c(name="transport (Iud)",option="magma",na.value="transparent") +
		scale_fill_viridis_c(name="transport (Iud)",option="magma",na.value="transparent") +
		theme_bw() +
		theme(axis.title=element_blank(),legend.position="top")

p %>% ggsave("map_flow.png",.,device="png",width=15,height=10,units="cm")