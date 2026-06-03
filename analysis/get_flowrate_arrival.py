import numpy as np
import netCDF4 as ncdf
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

def read_arrival_time(nc_file):

	nc_data = ncdf.Dataset(nc_file,'r')
	mat_tmin = nc_data["tmin"][:,:,:]

	Ntot = mat_tmin.shape[1]*mat_tmin.shape[2]

	vec_start = nc_data["start"][:]

	nc_data.close()

	for i in range(len(vec_start)):

		mat_tmin[:,i,:] += vec_start[i]

	return(mat_tmin)

def get_tmin_cumul(mat_tmin,index):

	Ntot = mat_tmin.shape[1]*mat_tmin.shape[2]
	vec_tmin = mat_tmin[index,:,:].flatten()
	vec_tmin = vec_tmin[vec_tmin < 999999]

	df_c = pd.Series(vec_tmin).value_counts().reset_index().rename(columns={"index":"tmin"}).sort_values(by="tmin")
	df_c["cumcount"] = np.cumsum(df_c["count"])

	return((df_c["tmin"].values,df_c["cumcount"].values/Ntot))

def plot_deriv(mat_tmin):

	fig_cum,axs_cum = plt.subplots(int(np.ceil(mat_tmin.shape[0]/20)),20,
				figsize=(6*20,int(np.ceil(mat_tmin.shape[0]/20))*4))
	fig_der,axs_der = plt.subplots(int(np.ceil(mat_tmin.shape[0]/20)),20,
				figsize=(6*20,int(np.ceil(mat_tmin.shape[0]/20))*4))

	k = 0
	for i in range(int(np.ceil(mat_tmin.shape[0]/20))):

		for j in range(20):

			print(k)

			if(k==385):
				break

			tup_tarr = get_tmin_cumul(mat_tmin,k)

			spl = sp.interpolate.splrep(tup_tarr[0],tup_tarr[1],s=0.01)

			inter_cumul = sp.interpolate.splev(np.arange(1,np.max(tup_tarr[0])+1),spl)
			deriv_cumul = sp.interpolate.splev(np.arange(1,np.max(tup_tarr[0])+1),spl,der=1)

			n = 365
			deriv_cumsum = np.cumsum(deriv_cumul)
			vec_deriv_cumsum = deriv_cumsum[n:]-deriv_cumsum[:-n]
			vec_deriv_cumsum = vec_deriv_cumsum[n - 1:]/n
			vec_time_cumsum = np.arange(n,int(len(deriv_cumul)))[n-1:]

			axs_cum[i,j].plot(tup_tarr[0]/365,tup_tarr[1],linewidth=4,color="black")
			axs_cum[i,j].plot(np.arange(1,np.max(tup_tarr[0])+1)/365,inter_cumul,linewidth=2,color="red")
			axs_cum[i,j].set_xlabel("time (years)")
			axs_cum[i,j].set_ylabel("fraction arrived") 
			axs_cum[i,j].grid(True)
			axs_cum[i,j].set_title("ID: "+str(k))

			axs_der[i,j].plot(np.arange(1,np.max(tup_tarr[0])+1)/365,deriv_cumul*365,linewidth=4,color="black")
			axs_der[i,j].plot(vec_time_cumsum/365,vec_deriv_cumsum*365,color="red",linewidth=2)
			axs_der[i,j].set_xlabel("time (years)")
			axs_der[i,j].set_ylabel("derivative fraction arrived (1/year)") 
			axs_der[i,j].grid(True)
			axs_der[i,j].set_title("ID: "+str(k))

			k += 1

		if(k==385):
			break
	
	fig_cum.subplots_adjust(bottom=0.01,top=0.99,left=0.01,right=0.99,hspace=0.3)
	fig_cum.savefig("cumul.pdf",format="pdf")

	fig_der.subplots_adjust(bottom=0.01,top=0.99,left=0.01,right=0.99,hspace=0.3)
	fig_der.savefig("cumul_derivative.pdf",format="pdf")

	'''plt.figure(figsize=(6, 4))

	plt.plot(np.arange(1,np.max(tup_tarr[0])+1)/365,deriv_cumul*365)
	plt.xlabel("time (years)")
	plt.ylabel("derivative fraction arrived (1/year)") 
	plt.grid(True)
	plt.savefig("cumul_deriv.pdf",format="pdf")'''


def main():

	nc_file = "/media/kobe/Windows/spectrum/network/sim_network_test_7647.nc"

	mat_tmin = read_arrival_time(nc_file)

	tup_tarr = get_tmin_cumul(mat_tmin,1)

	plot_deriv(mat_tmin)



if __name__ == '__main__':
	main()