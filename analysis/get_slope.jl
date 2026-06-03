using NCDatasets, DataStructures, DataFrames, GLM, CSV, Statistics

function get_slope(mat_tmin::Matrix{Int32},vec_start::Vector{Int32})

	mat_tmin[mat_tmin .== 999999] .= -999999

	vec_tmin = vec(mat_tmin)
	vec_tmin = vec_tmin[vec_tmin .> 0]

	tmin_min = findmin(vec_tmin)[1]*1.0
	tmin_max = findmax(vec_tmin)[1]*1.0
	tmin_med = median(vec_tmin)*1.0
	tmin_mn  = mean(vec_tmin)*1.0
	tmin_low = quantile(vec_tmin,0.025)*1.0
	tmin_hig = quantile(vec_tmin,0.975)*1.0

	for i in 1:length(vec_start)

		mat_tmin[:,i] = mat_tmin[:,i].+vec_start[i]

	end

	vec_tmin = vec(mat_tmin)

	vec_tmin = vec_tmin[vec_tmin .>  0]
	c = counter(vec_tmin)
	dict_time = sort(collect(c), by = x->x[1])

	vec_time = first.(dict_time)
	time_max = findmax(vec_time)[1]
	vec_cumsum = last.(dict_time)*1.0
	vec_cumsum = cumsum(vec_cumsum)*(1.0/size(mat_tmin,1)/size(mat_tmin,2))

	vec_R2 = zeros(Float64,ceil(Int,(time_max-365)/10))
	vec_coef = zeros(Float64,ceil(Int,(time_max-365)/10))

	k = 1
	for i in 1:10:(time_max-365)

		#println(k)

		df = DataFrame(time=vec_time[(vec_time .>= i) .& (vec_time .<= i+364)],
						cumsum=vec_cumsum[(vec_time .>= i) .& (vec_time .<= i+364)])

		if(size(df,1) >= 100)
			ols = lm(@formula(cumsum  ~ time),df)
			
			vec_R2[k] = r2(ols)
			vec_coef[k] = coef(ols)[2]
		end

		k += 1

	end

	i_max = findmax(vec_R2)

	return([i_max[1],vec_coef[i_max[2]],tmin_min,tmin_max,tmin_med,tmin_mn,
				tmin_low,tmin_hig])


end

function get_slopes(mat_tmin::Array{Int32,3},vec_start::Vector{Int32},i_cell::Int,ID::Int)

	mat_out = zeros(Float64,(size(mat_tmin,3),8))
	#vec_R2 = zeros(Float64,size(mat_tmin,3))
	#vec_coef = zeros(Float64,size(mat_tmin,3))

	for i in 1:size(mat_tmin,3)

		#println(i)

		if(i != i_cell)
			mat_out[i,:] = get_slope(mat_tmin[:,:,i],vec_start)
		end

	end

	df_out = DataFrame(mat_out,["R2","slope","tmin_min","tmin_max","tmin_med",
								"tmin_mn","tmin_low","tmin_hig"])
	df_out[!,"tmin_min"] = round.(Int,df_out[!,"tmin_min"])
	df_out[!,"tmin_max"] = round.(Int,df_out[!,"tmin_max"])

	CSV.write("out_slope/network_slopes_"*string(ID)*".csv",df_out)

end

function do_simulation(sim_dir::String,sim_file::String,vec_ID::Vector{Int})

	ID = split(sim_file,"_")[4]
	ID = parse(Int,split(ID,".")[1])
	i_ID = findall(x->x==ID,vec_ID)[1]

	nc_sim = NCDataset(sim_dir*sim_file)

	mat_tmin = nc_sim["tmin"][:,:,:]
	vec_start = nc_sim["start"][:]

	close(nc_sim)

	get_slopes(mat_tmin,vec_start,i_ID,ID)

end

function main()

	#sim_dir = "/media/kobe/Windows/spectrum/network/"
	sim_dir = "output/"

	if(!isdir("out_slope"))
		mkdir("out_slope")
	end

	sim_files = readdir(sim_dir)

	#vec_ID = CSV.read("../../healpix/test_grid.csv",DataFrame)[:,"ID"]
	vec_ID = CSV.read("test_grid.csv",DataFrame)[:,"ID"]

	for i in 1:length(sim_files)
		println(i)
		do_simulation(sim_dir,sim_files[i],vec_ID)
	end

end

main()
