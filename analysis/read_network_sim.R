library(tidyverse)
library(ncdf4)
library(sf)
sf_use_s2(FALSE)
library(ggpubr)
library(e1071)

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

sf_network <- st_read("../../healpix/test_grid/test_grid.shp")
df_network <- read.csv("../../healpix/test_grid.csv")
sf_land <- st_read("../../../land/ne_50m_land/ne_50m_land.shp") %>%
			st_crop(data.frame(x=c(-80,-80,-30,-30),y=c(25,70,70,25)) %>%
						st_as_sf(crs=4326,coords=c("x","y")))

nc_sim <- nc_open("/media/kobe/Windows/spectrum/network/sim_network_test_7647.nc")

str_pos0 <- nc_sim %>% ncatt_get(0)
str_pos0 <- str_pos0$'cell central position'
str_pos0 <- unlist(strsplit(substr(str_pos0,2,nchar(str_pos0)-9),split=" ; ")) %>%
				as.numeric()

str_pos0 <- c(-58.790321,43.406856)

df_network <- df_network %>% mutate(diff=abs(lon-str_pos0[1])+abs(lat-str_pos0[2]))
cell_id <- which(df_network$diff==min(df_network$diff))

mat_tmin <- nc_sim %>% ncvar_get("tmin")
mat_tmin[mat_tmin == 999999] <- NA
vec_tstart <- nc_sim %>% ncvar_get("start")

nc_sim %>% nc_close()

if(F){
plist_tmin <- list()
plist_arr <- list()

for(i in 1:dim(mat_tmin)[3]){

	print(i)

	if(i != cell_id){

		mat_sub <- mat_tmin[,,i]
		vec_tmin <- c(mat_sub)
		vec_tmin <- vec_tmin[!is.na(vec_tmin)]
		
		if(length(vec_tmin) > 1000){

			df_hist <- (vec_tmin/365) %>% make_hist()

			p <- df_hist %>% ggplot() + 
					geom_segment(aes(x=counts,y=min(density)/10,
								yend=density)) +
					geom_point(aes(x=counts,y=density),col="red") +
					scale_x_continuous(name=bquote(t[1]~(years))) +
					scale_y_continuous(name=bquote(P(t[1])~(years^{-1})),trans="log10") +
					theme_bw() +
					labs(title=paste0("ID ",i))

			plist_tmin <- append(plist_tmin,list(p))

			for(j in 1:dim(mat_sub)[2]){
				mat_sub[,j] <- mat_sub[,j] + vec_tstart[j]
			}

			vec_tarr <- c(mat_sub)
			vec_tarr <- vec_tarr[!is.na(vec_tarr)]

			df_cumsum <- vec_tarr %>% table() %>% as.data.frame() %>%
							rename(tarr=1,counts=2) %>%
							mutate(counts=counts/dim(mat_sub)[1]/dim(mat_sub)[2]) %>%
							mutate(tarr=tarr %>% as.character() %>% as.numeric())

			p <- df_cumsum %>% ggplot() +
					geom_point(aes(x=tarr/365,y=counts)) +
					scale_x_continuous(name="time (years)") +
					scale_y_continuous(name="fraction arrived",trans="log10") +
					theme_bw() +
					labs(title=paste0("ID ",i))

			plist_arr <- append(plist_arr,list(p))

		}else{
			plist_tmin <- append(plist_tmin,list(NULL))
			plist_arr <- append(plist_arr,list(NULL))
		}

	} else{

		vec_sum <- ((!is.na(mat_tmin)) %>% colSums(.,dims=1) %>% 
					colSums())/121/1024
		
		vec_sum[cell_id] <- NA
		vec_sum[vec_sum == 0] <- NA

		sf_network <- sf_network %>% mutate(arr=vec_sum)

		p <- ggplot() +
			geom_sf(data=sf_land) +
			geom_sf(data=sf_network,aes(col=arr,fill=arr)) +
			scale_x_continuous(expand=c(0,0),limits=c(-80,-25)) +
			scale_y_continuous(expand=c(0,0),limits=c(27.5,50)) +
			scale_colour_viridis_c(name="arrived",option="magma",na.value="transparent") +
			scale_fill_viridis_c(name="arrived",option="magma",na.value="transparent") +
			theme_bw() +
			theme(axis.title=element_blank(),legend.position="top")

		plist_tmin <- append(plist_tmin,list(p))
		plist_arr <- append(plist_arr,list(p))
	}

}

q <- ggarrange(plotlist=plist_tmin,ncol=20,nrow=ceiling(length(plist_tmin)/20))
q %>% ggsave("tmin_test.pdf",.,device=cairo_pdf,width=20*10,
				height=ceiling(length(plist_tmin)/20)*6.67,units="cm",limitsize=FALSE)

q <- ggarrange(plotlist=plist_arr,ncol=20,nrow=ceiling(length(plist_arr)/20))
q %>% ggsave("arrived_test.pdf",.,device=cairo_pdf,width=20*10,
				height=ceiling(length(plist_arr)/20)*6.67,units="cm",limitsize=FALSE)

}

df_sum <- data.frame(ID=rep(NA,dim(mat_tmin)[3]),
						Narr=rep(NA,dim(mat_tmin)[3]),
						Nrate=rep(NA,dim(mat_tmin)[3]),
						Tmin=rep(NA,dim(mat_tmin)[3]),
						Tmax=rep(NA,dim(mat_tmin)[3]),
						Tmed=rep(NA,dim(mat_tmin)[3]),
						Tmean=rep(NA,dim(mat_tmin)[3]),
						Tvar=rep(NA,dim(mat_tmin)[3]),
						Tskew=rep(NA,dim(mat_tmin)[3]),
						Tkurt=rep(NA,dim(mat_tmin)[3]),
						Tlow=rep(NA,dim(mat_tmin)[3]),
						Thig=rep(NA,dim(mat_tmin)[3])
						)

for(i in 1:dim(mat_tmin)[3]){

	print(i)
	df_sum$ID[i] <- i

	if(i != cell_id){

		mat_sub <- mat_tmin[,,i]
		vec_tmin <- c(mat_sub)
		vec_tmin <- vec_tmin[!is.na(vec_tmin)]

		df_sum$Tmin[i] <- min(vec_tmin)
		df_sum$Tmax[i] <- max(vec_tmin)
		df_sum$Tmed[i] <- median(vec_tmin)
		df_sum$Tmean[i] <- mean(vec_tmin)
		df_sum$Tvar[i] <- var(vec_tmin)
		df_sum$Tskew[i] <- skewness(vec_tmin)
		df_sum$Tkurt[i] <- kurtosis(vec_tmin)
		df_sum$Tlow[i] <- quantile(vec_tmin,0.25)
		df_sum$Thig[i] <- quantile(vec_tmin,0.75)

		for(j in 1:dim(mat_sub)[2]){
				mat_sub[,j] <- mat_sub[,j] + vec_tstart[j]
		}

		vec_tarr <- c(mat_sub)
		vec_tarr <- vec_tarr[!is.na(vec_tarr)]

		df_cumsum <- vec_tarr %>% table() %>% as.data.frame() %>%
						rename(tarr=1,counts=2) %>%
						mutate(counts=counts/dim(mat_sub)[1]/dim(mat_sub)[2]) %>%
						mutate(tarr=tarr %>% as.character() %>% as.numeric())

		dc <- log10(max(df_cumsum$counts)) - log10(min(df_cumsum$counts))
		df_cumsum <- df_cumsum %>% 
			filter(counts > 10^(log10(max(counts))-0.3*dc))

		df_sum$Nrate[i] <- 10^mean(log10(df_cumsum$counts))
		df_sum$Narr[i] <- length(vec_tarr)/dim(mat_sub)[1]/dim(mat_sub)[2]

	}

}