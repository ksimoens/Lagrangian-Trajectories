library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(ggpubr)
library(scales)

if(F){
vec_ID <- read.csv("../../healpix/test_grid.csv") %>% pull(ID)

list_files <- list.files("out_slope")

df_netw <- matrix(ncol=10,nrow=0) %>% as.data.frame() %>%
			rename(from=1,to=2,R2=3,slope=4,tmin_min=5,tmin_max=6,
					tmin_med=7,tmin_mn=8,tmin_low=9,tmin_hig=10)

for(i in 1:length(list_files)){

	print(i)

	from_i <- unlist(strsplit(list_files[i],split="_"))[3]
	from_i <- unlist(strsplit(from_i,split="[.]"))[1] %>% as.numeric()
	ID_i <- which(vec_ID == from_i)

	df_netw <- rbind(df_netw,
				read.csv(paste0("out_slope/",list_files[i])) %>%
					mutate(from=ID_i,to=1:nrow(.)) %>%
					filter(R2 != 0) %>%
					select(from,to,R2,slope,tmin_min,tmin_max,tmin_med,
							tmin_mn,tmin_low,tmin_hig))


}

df_netw %>% write.csv("network_test_grid.csv",row.names=FALSE)
}

df_netw <- read.csv("network_test_grid.csv") %>% 
			left_join(read.csv("../../plots/network/test_grid/test_grid_flow.csv") %>%
						mutate(from=1:nrow(.)) %>% select(-ID),by="from") %>%
			mutate(rate=slope*flow)

th_R2 <- 0.99

df_netw <- df_netw %>% filter(R2 > th_R2)

df_sum <- df_netw %>% group_by(from) %>% 
			summarise(out_degree=length(to),out_strength=sum(rate)) %>%
			rename(ID=from) %>%
			left_join(df_netw %>% group_by(to) %>%
						summarise(in_degree=length(from),in_strength=sum(rate)) %>%
						rename(ID=to),by="ID")

sf_grid <- st_read("../../healpix/test_grid/test_grid.shp")
sf_grid <- sf_grid %>% mutate(in_degree=df_sum %>% pull(in_degree),
								in_strength=df_sum %>% pull(in_strength),
								out_degree=df_sum %>% pull(out_degree),
								out_strength=df_sum %>% pull(out_strength))

sf_land <- st_read("../../../land/ne_50m_land/ne_50m_land.shp") %>%
			st_crop(data.frame(x=c(-80,-80,-30,-30),y=c(25,50,50,25)) %>%
						st_as_sf(crs=4326,coords=c("x","y")))

p_in_dg <- ggplot() + geom_sf(data=sf_grid,aes(col=in_degree/nrow(sf_grid),
								fill=in_degree/nrow(sf_grid))) +
			geom_sf(data=sf_land) +
			scale_colour_viridis_c(name="in degree",option="magma") +
			scale_fill_viridis_c(name="in degree",option="magma") +
			scale_x_continuous(expand=c(0,0),limits=c(NA,-25)) +
			scale_y_continuous(expand=c(0,0)) +
			theme_bw() +
			theme(panel.background=element_rect(fill="light blue"))

p_in_str <- ggplot() + geom_sf(data=sf_grid,aes(col=in_strength/max(in_strength),
								fill=in_strength/max(in_strength))) +
			geom_sf(data=sf_land) +
			scale_colour_viridis_c(name="in strength",option="magma") +
			scale_fill_viridis_c(name="in strength",option="magma") +
			scale_x_continuous(expand=c(0,0),limits=c(NA,-25)) +
			scale_y_continuous(expand=c(0,0)) +
			theme_bw() +
			theme(panel.background=element_rect(fill="light blue"))

p_out_dg <- ggplot() + geom_sf(data=sf_grid,aes(col=out_degree/nrow(sf_grid),
								fill=out_degree/nrow(sf_grid))) +
			geom_sf(data=sf_land) +
			scale_colour_viridis_c(name="out degree",option="magma") +
			scale_fill_viridis_c(name="out degree",option="magma") +
			scale_x_continuous(expand=c(0,0),limits=c(NA,-25)) +
			scale_y_continuous(expand=c(0,0)) +
			theme_bw() +
			theme(panel.background=element_rect(fill="light blue"))

p_out_str <- ggplot() + geom_sf(data=sf_grid,aes(col=out_strength/max(out_strength),
								fill=out_strength/max(out_strength))) +
			geom_sf(data=sf_land) +
			scale_colour_viridis_c(name="out strength",option="magma") +
			scale_fill_viridis_c(name="out strength",option="magma") +
			scale_x_continuous(expand=c(0,0),limits=c(NA,-25)) +
			scale_y_continuous(expand=c(0,0)) +
			theme_bw() +
			theme(panel.background=element_rect(fill="light blue"))

p_src_snk <- ggplot() + geom_sf(data=sf_grid,aes(col=out_strength-in_strength,
								fill=out_strength-in_strength)) +
			geom_sf(data=sf_land) +
			scale_colour_gradient2(name="out-in strength",low=muted("blue"),
									high=muted("red"),mid="grey",midpoint=0) +
			scale_fill_gradient2(name="out-in strength",low=muted("blue"),
									high=muted("red"),mid="grey",midpoint=0) +
			scale_x_continuous(expand=c(0,0),limits=c(NA,-25)) +
			scale_y_continuous(expand=c(0,0)) +
			theme_bw() +
			theme(panel.background=element_rect(fill="light blue"))

q <- ggarrange(p_in_dg,p_in_str,p_out_dg,p_out_str,p_src_snk,ncol=2,nrow=3)
q %>% ggsave("network_test_grid.png",.,device="png",width=24,height=6.67*3,units="cm",bg="white")

