library(ncdf4)

list_nc <- list.files("/media/kobe/shared/spectrum/vel_raw/",full.names=TRUE)

vec_lon <- list_nc[1] %>% nc_open(.) %>% ncvar_get("longitude")
vec_lat <- list_nc[1] %>% nc_open(.) %>% ncvar_get("latitude")
Nlon <- vec_lon %>% length()
Nlat <- vec_lat %>% length()

u_mat <- array(rep(NA,Nlon*Nlat*length(list_nc)),c(Nlon,Nlat,length(list_nc)))
v_mat <- array(rep(NA,Nlon*Nlat*length(list_nc)),c(Nlon,Nlat,length(list_nc)))

for(i in 1:length(list_nc)){

	print(i)

	nc_vel <- nc_open(list_nc[i])

	u_mat[,,i] <- nc_vel %>% ncvar_get("uo") %>% 
					apply(.,c(1,2),mean,na.rm=TRUE)
	v_mat[,,i] <- nc_vel %>% ncvar_get("vo") %>% 
					apply(.,c(1,2),mean,na.rm=TRUE)

	nc_vel %>% nc_close()

}

u_mat <- u_mat %>% apply(.,c(1,2),mean,na.rm=TRUE)
v_mat <- v_mat %>% apply(.,c(1,2),mean,na.rm=TRUE)

londim <- ncdim_def("lon","degrees_east",vec_lon)
latdim <- ncdim_def("lat","degrees_north",vec_lat)

u_var <- ncvar_def("u","m/s",list(londim,latdim))
v_var <- ncvar_def("v","m/s",list(londim,latdim))

nc_out <- nc_create("vel_average.nc",list(u_var,v_var))
ncvar_put(nc_out,u_var,u_mat)
ncvar_put(nc_out,v_var,v_mat)

nc_out %>% nc_close()