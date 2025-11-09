library(tidyverse)
library(ncdf4)
library(ggpubr)

expdist <- function(t,T){

	return(exp(-(t-1)/T)/(T*(1-exp(-(5-1)/T))))

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

if(F){
df_hist <- matrix(ncol=4,nrow=0) %>% as.data.frame() %>%
			rename(counts=V1,density=V2,distance=V3,simID=V4)

df_coef <- matrix(ncol=5,nrow=0) %>% as.data.frame() %>%
			rename(a_rms=V1,b_rms=V2,a_D=V3,b_D=V4,simID=V5)

df_velD <- matrix(ncol=4,nrow=0) %>% as.data.frame() %>%
			rename(distance=V1,D=V2,velrms=V3,simID=V4)

df_fit <- matrix(ncol=4,nrow=0) %>% as.data.frame() %>%
			rename(x=V1,y=V2,distance=V3,simID=V4)

distdir <- "/media/kobe/Windows/spectrum/distance/"

list_dist <- list.files(distdir)

for(sim_i in 1:length(list_dist)){

	print(sim_i)

nc_dist <- nc_open(paste0(distdir,list_dist[sim_i]))

dists <- nc_dist %>% ncvar_get("distance")

velrms_mat <- nc_dist %>% ncvar_get("velrms")
tmin_mat <- nc_dist %>% ncvar_get("tmin")

df_dist <- matrix(ncol=4,nrow=0) %>% as.data.frame() %>%
			rename(velrms=V1,D=V2,N=V3,distance=V4)

for(i in 1:length(dists)){

	times <- as.numeric(tmin_mat[,i,])
	times <- times[!is.na(times)]/365
	times <- times[times > 1]

	df_hist <- rbind(df_hist,times %>% make_hist() %>% 
						mutate(distance=dists[i],simID=sim_i))

	T_i <- optim(1,maxLL,method="Nelder-Mead",times=times)$par

	df_fit <- rbind(df_fit,
				data.frame(x=seq(min(times),max(times),length=100)) %>%
					mutate(y=expdist(x,T_i)) %>%
					mutate(distance=dists[i],simID=sim_i))

	df_dist <- rbind(df_dist,
				data.frame(velrms=as.numeric(velrms_mat[,i,])) %>%
					mutate(D=dists[i]^2/T_i*1e6/365/24/60/60) %>%
					mutate(N=length(times)) %>%
					mutate(distance=dists[i]))

}

df_dist <- df_dist %>% filter(!is.na(velrms))

lm_rms <- lm(data=df_dist,log10(velrms)~log10(distance))

coef_rms <- summary(lm_rms)$coefficients[1:2]

lm_D <- lm(data=df_dist %>% select(distance,D) %>% unique(),
					log10(D)~log10(distance))

coef_D <- summary(lm_D)$coefficients[1:2]

df_coef <- rbind(df_coef,
			data.frame(a_rms=coef_rms[2],b_rms=coef_rms[1],
						a_D=coef_D[2],b_D=coef_D[1],simID=sim_i))

df_velD <- rbind(df_velD,df_dist %>% group_by(distance) %>%
					summarise(D=mean(D),velrms=10^mean(log10(velrms))) %>%
					mutate(simID=sim_i))

}
}

p <- df_hist %>% ggplot() + 
		geom_line(aes(x=counts,y=density,col=factor(distance),group=interaction(simID,distance))) +
		scale_x_continuous(name=bquote(t[1]~(years))) +
		scale_y_continuous(name=bquote(P(t[1])~(years^{-~1})),trans="log10") +
		scale_colour_viridis_d(name="distance (km)",option="magma") +
		theme_bw()

p %>% ggsave("hist_distance.png",.,device="png",width=15,height=10,units="cm")

p <- df_velD %>% ggplot() + 
					geom_boxplot(aes(x=distance,y=velrms,group=distance)) +
					scale_x_continuous(name="distance (km)") +
					scale_y_continuous(name=bquote(sqrt("<"~v[tot]^2~">")~(m/s)),trans="log10") +
					theme_bw()

p %>% ggsave("summary_velrms_distance.png",.,device="png",width=15,height=10,units="cm")

p <- df_velD %>% ggplot() + 
					geom_boxplot(aes(x=distance,y=D,group=distance)) +
					scale_x_continuous(name="distance (km)") +
					scale_y_continuous(name=bquote(D~(m^2/s)),trans="log10") +
					theme_bw()

p %>% ggsave("summary_D_distance.png",.,device="png",width=15,height=10,units="cm")

phist_list <- list()
pD_list <- list()
pvel_list <- list()
IDs <- df_hist %>% pull(simID) %>% unique()
for(id in IDs){

	p <- ggplot() + 
		geom_point(data=df_hist %>% filter(simID==id),aes(x=counts,y=density,col=factor(distance))) +
		geom_line(data=df_fit %>% filter(simID==id),aes(x=x,y=y,col=factor(distance))) +
		scale_x_continuous(name=bquote(t[1]~(years))) +
		scale_y_continuous(name=bquote(P(t[1])~(years^{-1})),trans="log10") +
		scale_colour_viridis_d(name="distance (km)",option="magma") +
		theme_bw() +
		theme(panel.background=element_rect(fill="lightgrey")) +
		guides(col = guide_legend(nrow = 1))

	phist_list <- append(phist_list,list(p))

	df_fit_i <- data.frame(x=10^seq(log10(475),log10(1025),length=100)) %>%
				mutate(y=10^((df_coef %>% filter(simID==id) %>% pull(b_D)) +
							(df_coef %>% filter(simID==id) %>% pull(a_D))*log10(x)))

	p <- ggplot() +
		geom_point(data=df_velD %>% filter(simID==id),aes(x=distance,y=D)) +
		geom_line(data=df_fit_i,aes(x=x,y=y)) +
		scale_x_continuous(name="distance (km)",trans="log10") +
		scale_y_continuous(name=bquote(D~(m^2/s)),trans="log10") +
		theme_bw() +
		labs(title=sprintf("slope: %.2f",df_coef %>% filter(simID==id) %>% pull(a_D)))

	pD_list <- append(pD_list,list(p)) 

	df_fit_i <- data.frame(x=10^seq(log10(475),log10(1025),length=100)) %>%
				mutate(y=10^((df_coef %>% filter(simID==id) %>% pull(b_rms)) +
							(df_coef %>% filter(simID==id) %>% pull(a_rms))*log10(x)))

	p <- ggplot() +
		geom_point(data=df_velD %>% filter(simID==id),aes(x=distance,y=velrms)) +
		geom_line(data=df_fit_i,aes(x=x,y=y)) +
		scale_x_continuous(name="distance (km)",trans="log10") +
		scale_y_continuous(name=bquote(sqrt("<"~v[tot]^2~">")~(m/s)),trans="log10") +
		theme_bw() +
		labs(title=sprintf("slope: %.2f",df_coef %>% filter(simID==id) %>% pull(a_rms)))

	pvel_list <- append(pvel_list,list(p))

}

q <- ggarrange(plotlist=phist_list,ncol=5,nrow=ceiling(length(p_list)/5),
				common.legend=TRUE,legend="top")
q %>% ggsave("fit_hist.pdf",.,device=cairo_pdf,width=5*10,height=ceiling(length(p_list)/5)*6.67,
				units="cm",bg="white")

q <- ggarrange(plotlist=pD_list,ncol=5,nrow=ceiling(length(pD_list)/5))
q %>% ggsave("summary_D_distance.pdf",.,device=cairo_pdf,width=5*10,height=ceiling(length(p_list)/5)*6.67,
				units="cm",bg="white")

q <- ggarrange(plotlist=pvel_list,ncol=5,nrow=ceiling(length(pvel_list)/5))
q %>% ggsave("summary_velrms_distance.pdf",.,device=cairo_pdf,width=5*10,height=ceiling(length(p_list)/5)*6.67,
				units="cm",bg="white")

p <- df_coef %>% filter(a_D > 0 & a_rms > 0) %>% ggplot() +
		geom_point(aes(x=a_D,y=a_rms)) +
		scale_x_continuous(name="slope D") +
		scale_y_continuous(name=bquote(slope~sqrt("<"~v[tot]^2~">"))) +
		theme_bw()

p %>% ggsave("compare_slopes.png",.,device="png",width=11,height=10,units="cm")