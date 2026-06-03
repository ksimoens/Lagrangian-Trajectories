library(tidyverse)
library(sf)
sf_use_s2(FALSE)

if(F){
d <- read.csv("test_slope.csv") %>% filter(R2 != 0)

p <- d %>% ggplot() + geom_point(aes(x=time,y=R2),col="red") +
			#geom_point(aes(x=time,y=coef*1e4),col="blue") +
			scale_x_continuous(name="time") +
			theme_bw()

p %>% ggsave("test_slope.png",.,device="png",width=10,height=6.67,units="cm")
}

d <- read.csv("out_slope/network_slopes_11498.csv")
d[d==0] <- NA

sf_grid <- st_read("../../healpix/test_grid/test_grid.shp")
sf_land <- st_read("../../../land/ne_50m_land/ne_50m_land.shp") %>%
			st_crop(data.frame(x=c(-80,-80,-30,-30),y=c(25,50,50,25)) %>%
						st_as_sf(crs=4326,coords=c("x","y")))

sf_grid <- sf_grid %>% mutate(R2=d %>% pull(R2),slope=d %>% pull(coef))

p <- ggplot() + geom_sf(data=sf_land) +
					geom_sf(data=sf_grid[342,],col="green",fill="green") +
					geom_sf(data=sf_grid,aes(col=slope*365,fill=slope*365)) +
					scale_colour_viridis_c(name="1/year",option="magma",na.value="transparent") +
					scale_fill_viridis_c(name="1/year",option="magma",na.value="transparent") +
					theme_bw() +
					scale_x_continuous(expand=c(0,0),limits=c(NA,-25)) +
					scale_y_continuous(expand=c(0,0)) 

p %>% ggsave("map_R2.png",.,device="png",width=10,height=5,units="cm")



