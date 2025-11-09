using NCDatasets,Statistics,DataStructures

veldir = "/media/kobe/Windows/spectrum/output/"
#distdir = "/media/kobe/Windows/spectrum/distance/"
distdir = "/media/kobe/Windows/spectrum/target/"

list_sim = readdir(veldir)

for sim_i in 3:length(list_sim)

	println(sim_i)

nc_sim = NCDataset(veldir*list_sim[sim_i])

u_mat = nc_sim["u"][:,:,:,:]
v_mat = nc_sim["v"][:,:,:,:]

close(nc_sim)

u_mat_new = zeros(Float32,(size(u_mat,1),size(u_mat,2)*size(u_mat,3),size(u_mat,4)))
v_mat_new = zeros(Float32,(size(u_mat,1),size(u_mat,2)*size(u_mat,3),size(u_mat,4)))

for i in 1:size(u_mat,3)

	u_mat_new[:,((i-1)*size(u_mat,2)+1):((i)*size(u_mat,2)),:] =
		u_mat[:,:,i,:]
	v_mat_new[:,((i-1)*size(v_mat,2)+1):((i)*size(v_mat,2)),:] =
		v_mat[:,:,i,:]

end

nc_dist = NCDataset(distdir*list_sim[sim_i],"a")

tmin_mat = nc_dist["tmin"][:,:,:]
velrms_mat = zeros(Float32,size(tmin_mat))

for i in 1:size(tmin_mat,1)
	for j in 1:size(tmin_mat,2)
		for k in 1:size(tmin_mat,3)

			if(!isnan(tmin_mat[i,j,k]))

				velrms_mat[i,j,k] = sqrt(mean(u_mat_new[i,1:round(Int,tmin_mat[i,j,k]),k].^2+
												v_mat_new[i,1:round(Int,tmin_mat[i,j,k]),k].^2))

			else

				velrms_mat[i,j,k] = NaN

			end

		end
	end
end

#velrms_var = defVar(nc_dist,"velrms",Float32,("part","distance","start"),
velrms_var = defVar(nc_dist,"velrms",Float32,("part","target","start"),
	 			attrib = OrderedDict("units" => "m/s"))
velrms_var[:,:,:] = velrms_mat

close(nc_dist)

end