library(tidyverse)
library(ncdf4)

nc_SST <- nc_open("/media/kobe/shared/spectrum/SST/SST_2026_03_31.nc")

mat_temp <- nc_SST %>% ncvar_get("to")
mat_temp[is.na(mat_temp)] <- -999

vec_lon <- nc_SST %>% ncvar_get("longitude")
vec_lat <- nc_SST %>% ncvar_get("latitude")

nc_SST %>% nc_close()

londim <- ncdim_def("lon","radians east",vec_lon)
latdim <- ncdim_def("lat","radians north",vec_lat)
SSTvar <- ncvar_def("SST","°C",list(londim,latdim),prec="float")
nc_new <- nc_create("/media/kobe/shared/spectrum/SST/SST_2026_03_31_transf.nc",list(SSTvar),force_v4=TRUE)
nc_new %>% ncvar_put(SSTvar,mat_temp)

nc_new %>% nc_close()