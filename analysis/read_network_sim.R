library(tidyverse)
library(ncdf4)
library(sf)
sf_use_s2(FALSE)
library(ggpubr)

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

nc_sim <- nc_open("/media/kobe/Windows/spectrum/network/test_sim.nc")

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
							mutate(Narr=cumsum(counts)/dim(mat_sub)[1]/dim(mat_sub)[2]) %>%
							mutate(tarr=tarr %>% as.character() %>% as.numeric())

			p <- df_cumsum %>% ggplot() +
					geom_line(aes(x=tarr/365,y=Narr)) +
					scale_x_continuous(name="time (years)") +
					scale_y_continuous(name="fraction arrived") +
					theme_bw() +
					labs(title=paste0("ID ",i))

			plist_arr <- append(plist_arr,list(p))

		}else{
			plist_tmin <- append(plist_tmin,list(NULL))
			plist_arr <- append(plist_arr,list(NULL))
		}

	} else{
		plist_tmin <- append(plist_tmin,list(NULL))
		plist_arr <- append(plist_arr,list(NULL))
	}

}

q <- ggarrange(plotlist=plist_tmin,ncol=20,nrow=ceiling(length(plist_tmin)/20))
q %>% ggsave("tmin_test.pdf",.,device=cairo_pdf,width=20*10,
				height=ceiling(length(plist_tmin)/20)*6.67,units="cm",limitsize=FALSE)

q <- ggarrange(plotlist=plist_arr,ncol=20,nrow=ceiling(length(plist_arr)/20))
q %>% ggsave("arrived_test.pdf",.,device=cairo_pdf,width=20*10,
				height=ceiling(length(plist_arr)/20)*6.67,units="cm",limitsize=FALSE)