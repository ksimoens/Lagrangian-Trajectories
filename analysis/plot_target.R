library(tidyverse)
library(ncdf4)
library(ggpubr)
library(sf)
sf_use_s2(FALSE)
library(sp)

expdist <- function(t,T){

	return(exp(-(t-2)/T)/(T*(1-exp(-(5-2)/T))))

}

maxLL <- function(T,times){

	return(-sum(log(expdist(times,T))))

}

make_hist <- function(vec){

	vec <- vec[!is.na(vec)]
	vec <- vec[vec != 0]
	vec <- vec[vec != Inf]
	vec <- vec[vec != -Inf]

	breaks <- seq(min(vec),max(vec),length=15)
	hist <- hist(vec,breaks=breaks,plot=FALSE)

	#breaks_log2 <- c(hist_log$breaks[2:length(hist_log$breaks)],NA)
	#mids <- sqrt(hist_log$breaks*breaks_log2)

	return(data.frame(counts=hist$mids,density=hist$density) %>%
			filter(density != 0))

}

crs_atl <- CRS("+proj=sinu +lon_0=0 +a=6378140.0 +b=6356750.0 +units=m +no_defs")

sarg <- st_read("sargasso_shapefile/WC_13_EBSA.shp") %>%
			st_geometry() %>% st_convex_hull() %>% 
			st_transform(crs=crs_atl) %>%
			st_buffer(-200000)

sf_init <- (read.csv("../initial/sargasso_points_pos.csv")*1e3) %>%
			st_as_sf(crs=crs_atl,coords=c("x","y"))

R <- sqrt(100*100/pi)*1e3

if(T){

df_coef <- matrix(ncol=5,nrow=0) %>% as.data.frame() %>%
			rename(a_rms=V1,b_rms=V2,a_D=V3,b_D=V4,simID=V5)

distdir <- "/media/kobe/Windows/spectrum/target/"
plotdir <- "/media/kobe/Windows/spectrum/target_plot/"

list_dist <- list.files(distdir)

for(sim_i in 6:length(list_dist)){

	print(sim_i)

	df_hist <- matrix(ncol=5,nrow=0) %>% as.data.frame() %>%
			rename(counts=V1,density=V2,distance=V3,lon=V4,lat=V5)

df_velD <- matrix(ncol=6,nrow=0) %>% as.data.frame() %>%
			rename(distance=V1,lon=V2,lat=V3,T=V4,D=V5,velrms=V6)

df_fit <- matrix(ncol=5,nrow=0) %>% as.data.frame() %>%
			rename(x=V1,y=V2,distance=V3,lon=V4,lat=V5)

nc_dist <- nc_open(paste0(distdir,list_dist[sim_i]))

dists <- nc_dist %>% ncvar_get("target distance")
lons <- nc_dist %>% ncvar_get("target longitude")
lats <- nc_dist %>% ncvar_get("target latitude")

velrms_mat <- nc_dist %>% ncvar_get("velrms")
tmin_mat <- nc_dist %>% ncvar_get("tmin")

p_list <- list()

for(i in 1:length(dists)){

	times <- as.numeric(tmin_mat[,i,])
	times <- times[!is.na(times)]/365
	times <- times[times > 2]

	df_hist_sub <- times %>% make_hist() %>% 
						mutate(distance=dists[i],lon=lons[i],lat=lats[i],
								simID=sim_i)
	df_hist <- rbind(df_hist,df_hist_sub)

	T_i <- optim(1,maxLL,method="Nelder-Mead",times=times)$par

	df_fit_sub <- data.frame(x=seq(min(times),max(times),length=100)) %>%
					mutate(y=expdist(x,T_i)) %>%
					mutate(distance=dists[i],lon=lons[i],lat=lats[i],
							simID=sim_i)
	df_fit <- rbind(df_fit,df_fit_sub)

	p <- ggplot() + 
			geom_point(data=df_hist_sub,aes(x=counts,y=density)) +
			geom_line(data=df_fit_sub,aes(x=x,y=y)) +
			scale_x_continuous(name=bquote(t[1]~(years))) +
			scale_y_continuous(name=bquote(P(t[1])~(years^{-1})),trans="log10") +
			labs(title=sprintf("%i km; T = %.2f years",dists[i],T_i)) +
			theme_bw()

	p_list <- append(p_list,list(p))


	df_velD <- rbind(df_velD,
				data.frame(velrms=as.numeric(velrms_mat[,i,])) %>%
					mutate(T=T_i) %>%
					mutate(D=dists[i]^2/T_i*1e6/365/24/60/60) %>%
					mutate(N=length(times)) %>%
					mutate(distance=dists[i],lon=lons[i],lat=lats[i]))

}

df_velD <- df_velD %>% filter(!is.na(velrms) & velrms < 5)

lm_rms <- lm(data=df_velD,log10(velrms)~log10(distance))

coef_rms <- summary(lm_rms)$coefficients[1:2]

lm_D <- lm(data=df_velD %>% select(distance,D,T) %>% unique() %>%
					filter(T < 100),
					log10(D)~log10(distance))

coef_D <- summary(lm_D)$coefficients[1:2]

df_coef <- rbind(df_coef,
			data.frame(a_rms=coef_rms[2],b_rms=coef_rms[1],
						a_D=coef_D[2],b_D=coef_D[1],simID=sim_i))

df_fit_i <- data.frame(x=10^seq(log10(475),log10(1025),length=100)) %>%
				mutate(y=10^(coef_D[1] + coef_D[2]*log10(x)))

p <- ggplot() +
		geom_boxplot(data=df_velD %>% select(lon,lat,D,T,distance) %>% unique() %>%
							filter(T < 100),
						aes(x=distance,y=D,group=distance)) +
		geom_line(data=df_fit_i,aes(x=x,y=y),col="red",linewidth=1) +
		scale_x_continuous(name="distance (km)",trans="log10") +
		scale_y_continuous(name=bquote(D~(m^2/s)),trans="log10") +
		theme_bw() +
		labs(title=sprintf("slope: %.2f",coef_D[2]))

p_list <- append(p_list,list(p))

df_fit_i <- data.frame(x=10^seq(log10(475),log10(1025),length=100)) %>%
				mutate(y=10^(coef_rms[1] + coef_rms[2]*log10(x)))

p <- ggplot() +
		geom_boxplot(data=df_velD,aes(x=distance,y=velrms,group=distance)) +
		geom_line(data=df_fit_i,aes(x=x,y=y),col="red",linewidth=1) +
		scale_x_continuous(name="distance (km)",trans="log10") +
		scale_y_continuous(name=bquote(sqrt("<"~v[tot]^2~">")~(m/s)),trans="log10") +
		theme_bw() +
		labs(title=sprintf("slope: %.2f",coef_rms[2]))

p_list <- append(p_list,list(p))

p <- ggplot() + 
		geom_line(data=df_hist,aes(x=counts,y=density,col=factor(distance),group=interaction(lon,lat))) +
		scale_x_continuous(name=bquote(t[1]~(years))) +
		scale_y_continuous(name=bquote(P(t[1])~(years^{-1})),trans="log10") +
		scale_colour_viridis_d(name="distance (km)",option="magma") +
		theme_bw() +
		theme(panel.background=element_rect(fill="lightgrey"),
				legend.position="top",legend.title=element_blank()) +
		guides(col = guide_legend(nrow = 1))

p_list <- append(p_list,list(p))

p <- ggplot() + geom_sf(data=sarg,col="red",fill=NA,linewidth=1) +
					geom_sf(data=sf_init[sim_i,]) +
					geom_sf(data=df_velD %>% select(lon,lat,D,T) %>% unique() %>%
									filter(T < 100) %>%
									st_as_sf(crs=4326,coords=c("lon","lat")) %>%
									st_transform(crs=crs_atl) %>% st_buffer(R),
									aes(fill=D,col=D)) +
					scale_fill_viridis_c(name=bquote(D~(m^2/s)),trans="log10",option="magma") +
					scale_colour_viridis_c(name=bquote(D~(m^2/s)),trans="log10",option="magma") +
					theme_bw() +
					theme(panel.background=element_rect(fill="lightblue"),
							legend.position="top")  

p_list <- append(p_list,list(p))

p <- ggplot() + geom_sf(data=sarg,col="red",fill=NA,linewidth=1) +
					geom_sf(data=sf_init[sim_i,]) +
					geom_sf(data=df_velD %>% select(lon,lat,velrms) %>% group_by(lon,lat) %>%
									summarise(velrms=mean(velrms)) %>%
									st_as_sf(crs=4326,coords=c("lon","lat")) %>%
									st_transform(crs=crs_atl) %>% st_buffer(R),
									aes(fill=velrms,col=velrms)) +
					scale_fill_viridis_c(name=bquote(bar(sqrt("<"~v[tot]^2~">"))~(m/s)),option="magma") +
					scale_colour_viridis_c(name=bquote(bar(sqrt("<"~v[tot]^2~">"))~(m/s)),option="magma") +
					theme_bw() +
					theme(panel.background=element_rect(fill="lightblue"),
							legend.position="top")  

p_list <- append(p_list,list(p))

q <- ggarrange(plotlist=p_list,ncol=5,nrow=ceiling(length(p_list)/5))
q %>% ggsave(paste0(plotdir,substr(list_dist[sim_i],1,nchar(list_dist[sim_i])-3),".pdf"),.,
				device=cairo_pdf,width=5*10,height=ceiling(length(p_list)/5)*6.67,units="cm",
				limitsize=FALSE)

}
}

p <- df_coef %>% filter(a_D > 0 & a_rms > 0) %>% ggplot() +
		geom_point(aes(x=a_D,y=a_rms)) +
		scale_x_continuous(name="slope D") +
		scale_y_continuous(name=bquote(slope~sqrt("<"~v[tot]^2~">"))) +
		theme_bw()

p %>% ggsave("compare_slopes_target.png",.,device="png",width=11,height=10,units="cm")

sf_init <- sf_init %>% mutate(simID=1:nrow(.)) %>%
						left_join(df_coef %>% select(a_D,a_rms,simID),by="simID")
#sf_init$a_D[1:5] <- c(2.76,0.27,1.46,0.21,1.98)
#sf_init$a_rms[1:5] <- c(0.01,-0.21,-0.03,-0.13,0.18)

p <- ggplot() + geom_sf(data=sarg,col="red",fill=NA,linewidth=1) +
					geom_sf(data=sf_init,aes(col=a_D),size=2) +
					scale_colour_viridis_c(name="slope D",option="magma") +
					theme_bw() +
					theme(panel.background=element_rect(fill="lightblue"),
							legend.position="top") 

p %>% ggsave("map_slopeD_target.png",.,device="png",width=15,height=10,units="cm")

p <- ggplot() + geom_sf(data=sarg,col="red",fill=NA,linewidth=1) +
					geom_sf(data=sf_init,aes(col=a_rms),size=2) +
					scale_colour_viridis_c(name=bquote(slope~sqrt("<"~v[tot]^2~">")),option="magma") +
					theme_bw() +
					theme(panel.background=element_rect(fill="lightblue"),
							legend.position="top") 

p %>% ggsave("map_slopeVelrms_target.png",.,device="png",width=15,height=10,units="cm")
