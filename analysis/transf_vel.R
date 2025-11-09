library(tidyverse)
library(ncdf4)

mu_lat <- function(lat){

	return(log(abs(1/cos(lat)+tan(lat))))

}

for(year in 1994:2024){

	print(year)

nc_vel <- nc_open(paste0("/media/kobe/Shared Partition/spectrum/vel_raw/vel_",year,".nc"))

mat_u <- nc_vel %>% ncvar_get("uo")
mat_u[is.na(mat_u)] <- -999
mat_v <- nc_vel %>% ncvar_get("vo")
mat_v[is.na(mat_v)] <- -999

lon <- nc_vel %>% ncvar_get("longitude")
lon <- lon/180*pi
lat <- nc_vel %>% ncvar_get("latitude")
mu <- mu_lat(lat/180*pi)

nc_close(nc_vel)

londim <- ncdim_def("lon","radians east",lon)
mudim <- ncdim_def("mu","transformed radians north",mu)
timedim <- ncdim_def("time","days",1:dim(mat_u)[3])
uvar <- ncvar_def("u","m/s",list(londim,mudim,timedim))
vvar <- ncvar_def("v","m/s",list(londim,mudim,timedim))
nc_new <- nc_create(paste0("/media/kobe/Shared Partition/spectrum/vel_transf/vel_",year,".nc"),list(uvar,vvar))
nc_new %>% ncvar_put(uvar,mat_u)
nc_new %>% ncvar_put(vvar,mat_v)

nc_new %>% nc_close()
}