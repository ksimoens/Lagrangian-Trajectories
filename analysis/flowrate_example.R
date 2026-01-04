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

Nside <- 64
ds <- pi/Nside/2

df_grid <- read.csv("../../healpix/test_grid.csv")
sf_grid <- st_read("../../healpix/test_grid/test_grid.shp")

nc_vel <- nc_open("/media/kobe/shared/spectrum/vel_raw/vel_1993.nc")

vec_lon <- nc_vel %>% ncvar_get("longitude")
vec_lat <- nc_vel %>% ncvar_get("latitude")

nc_vel %>% nc_close()

cell <- df_grid %>% slice(1)
sf_cell <- data.frame(x=c(cell$x-ds/2,cell$x,cell$x+ds/2,cell$x),
						y=c(cell$y,cell$y+ds/2,cell$y,cell$y-ds/2)) %>%
			mutate(lon=xy_to_lon(x,y),lat=xy_to_lat(y))

df_line <- sf_cell %>% slice(1:2)
sf_line <- df_line %>% st_as_sf(crs=4326,coords=c("lon","lat")) %>%
			mutate(ID=1) %>% group_by(ID) %>% summarise(do_union=FALSE) %>%
			st_cast("LINESTRING")

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

hor_total <- hor_lines
pts_inter <- st_intersection(hor_lines,sf_line)
pts_veloc <- st_intersection(sf_vel_line,hor_lines)

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

ver_total <- ver_lines
pts_inter <- rbind(pts_inter,st_intersection(ver_lines,sf_line))
pts_veloc <- rbind(pts_veloc,st_intersection(sf_vel_line,ver_lines))

# top right

df_line <- sf_cell %>% slice(2:3)
sf_line <- df_line %>% st_as_sf(crs=4326,coords=c("lon","lat")) %>%
			mutate(ID=1) %>% group_by(ID) %>% summarise(do_union=FALSE) %>%
			st_cast("LINESTRING")

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

hor_total <- rbind(hor_total,hor_lines)
pts_inter <- rbind(pts_inter,st_intersection(hor_lines,sf_line))
pts_veloc <- rbind(pts_veloc,st_intersection(sf_vel_line,hor_lines))

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

ver_total <- rbind(ver_total,ver_lines)
pts_inter <- rbind(pts_inter,st_intersection(ver_lines,sf_line))
pts_veloc <- rbind(pts_veloc,st_intersection(sf_vel_line,ver_lines))

# bottom right

df_line <- sf_cell %>% slice(3:4)
sf_line <- df_line %>% st_as_sf(crs=4326,coords=c("lon","lat")) %>%
			mutate(ID=1) %>% group_by(ID) %>% summarise(do_union=FALSE) %>%
			st_cast("LINESTRING")

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

hor_total <- rbind(hor_total,hor_lines)
pts_inter <- rbind(pts_inter,st_intersection(hor_lines,sf_line))
pts_veloc <- rbind(pts_veloc,st_intersection(sf_vel_line,hor_lines))

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

ver_total <- rbind(ver_total,ver_lines)
pts_inter <- rbind(pts_inter,st_intersection(ver_lines,sf_line))
pts_veloc <- rbind(pts_veloc,st_intersection(sf_vel_line,ver_lines))

# bottom left

df_line <- sf_cell %>% slice(c(1,4))
sf_line <- df_line %>% st_as_sf(crs=4326,coords=c("lon","lat")) %>%
			mutate(ID=1) %>% group_by(ID) %>% summarise(do_union=FALSE) %>%
			st_cast("LINESTRING")

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

hor_total <- rbind(hor_total,hor_lines)
pts_inter <- rbind(pts_inter,st_intersection(hor_lines,sf_line))
pts_veloc <- rbind(pts_veloc,st_intersection(sf_vel_line,hor_lines))

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

ver_total <- rbind(ver_total,ver_lines)
pts_inter <- rbind(pts_inter,st_intersection(ver_lines,sf_line))
pts_veloc <- rbind(pts_veloc,st_intersection(sf_vel_line,ver_lines))

p <- ggplot() + geom_sf(data=sf_grid %>% slice(1),col="red",linewidth=1) +
		geom_sf(data=hor_total,col="black") +
		geom_sf(data=ver_total,col="blue") +
		geom_sf(data=pts_veloc,col="grey") +
		geom_sf(data=pts_inter,col="green") +
		theme_bw()

p %>% ggsave("flowrate_example.png",.,device="png",width=10,height=10,units="cm")
