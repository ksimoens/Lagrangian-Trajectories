import copernicusmarine

for year in range(1998,2025):

  print(year)

  copernicusmarine.subset(
    dataset_id="cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m",
    variables=["uo", "vo"],
    minimum_longitude=-90,
    maximum_longitude=15,
    minimum_latitude=10,
    maximum_latitude=70,
    start_datetime=str(year)+"-01-01T00:00:00",
    end_datetime=str(year)+"-12-31T23:59:59",
    minimum_depth=15,
    maximum_depth=15,
    output_filename = "/media/kobe/Shared Partition/spectrum/vel_raw/vel_"+str(year)+".nc"
  )