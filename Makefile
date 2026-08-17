#Total number of years to integrate.
OPT += -DNYEAR=1
#Total number of months to integrate.
#OPT += -DNMONTH=1
#Total number of years over which to start.
OPT += -DNYEARSTART=1
#Number of days between starting moments.
OPT += -DDTSTART=15
#Number of particles.
OPT += -DNPART=1000
#Starting year.
OPT += -DYSTART=1993
#Starting month.
OPT += -DMSTART=3
#Diffusion constant.
OPT += -DD=0.0
#Timestep (seconds)
OPT += -DDT=86400.0
#Shape of velocity grid.
OPT += -DNLON=420
OPT += -DNLAT=240
#Extent of velocity grid.
OPT += -DLONMIN=-1.568614665
OPT += -DLONMAX=0.259617726
OPT += -DLONRES=0.004363323
OPT += -DLATMIN=0.1767146
OPT += -DLATMAX=1.219549
OPT += -DLATRES=0.004363323
OPT += -DDAY
#OPT += -DHOUR
#Particle properties.
#OPT += -DSTOREPOS
#OPT += -DSTOREVEL
#OPT += -DCIRCULAR
#OPT += -DNETWORK
#OPT += -DLYAPUNOV
#OPT += -DSST
#OPT += -DBROWNIAN
#OPT += -DLYAPSST
#OPT += -DLYAPCIRC
#OPT += -DLYAPRATIO
#OPT += -DLYAPINFINITY
#OPT += -DTRACERPATH
OPT += -DDISTRIBUTION
#Network properties.
#OPT += -DNCELL=385
#OPT += -DNSIDE=64
#OPT += -DNETLONMIN=-1.155371
#OPT += -DNETLONMAX=-0.4154256
#OPT += -DNETLATMIN=0.4493061
#OPT += -DNETLATMAX=0.9618006
#Lyapunov properties.
#OPT += -DDEND=1.0
#Brownian properties.
#OPT += -DNTIMES=11
#Lyapunov SST properties.
#OPT += -DSCALEFACTOR=1.4142
#Lyapunov circle properties.
#OPT += -DNCIRC=3
OPT += -DRADIUS=500
#Tracer properties
#OPT += -DTRACER_CHL
#OPT += -DTRACER_SST
#Output grid (Lyapunov / SST / Brownian / LyapSST)
#OPT += -DOUTLONMIN=-1.568614665
OPT += -DOUTLONMIN=-0.95993108859
#OPT += -DOUTLONMAX=0.259617726
OPT += -DOUTLONMAX=-0.61086523819
#OPT += -DOUTLONRES=0.004363323
OPT += -DOUTLONRES=0.06981317007
#OPT += -DOUTLATMIN=0.1767146
OPT += -DOUTLATMIN=0.47996554429
#OPT += -DOUTLATMAX=1.219549
OPT += -DOUTLATMAX=0.82903139469
#OPT += -DOUTLATRES=0.004363323
OPT += -DOUTLATRES=0.06981317007
#Shape of SST data grid
OPT += -DSSTGRIDNLON=1261
OPT += -DSSTGRIDNLAT=720
#Extent of SST data grid
OPT += -DSSTGRIDLONMIN=-1.570796
OPT += -DSSTGRIDLONMAX=0.2617994
OPT += -DSSTGRIDLONRES=0.001454485
OPT += -DSSTGRIDLATMIN=0.1759873
OPT += -DSSTGRIDLATMAX=1.22173
OPT += -DSSTGRIDLATRES=0.001454485

#Variables for transformation.
OPT += -DA=6378.140
OPT += -DB=6356.7500
OPT += -DE2=0.0066960376522835485
OPT += -DE1=0.0016796375940427949
OPT += -DR=6367445.0

#Configuration file.
OPT += -DCONFIGPATH=config

export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

OBJS = main.o globParams.o Grid.o Particle.o Vec.o read_config.o
CFILES = main.cpp globParams.cpp Grid.cpp Particle.cpp Vec.cpp read_config.cpp
OPT_OPENMP = -fopenmp
EXEC = main
CXX = g++
OPTIMIZE = -O3 -Wall -g
INCL = Makefile globParams.h Grid.h Particle.h Vec.h read_config.h
NCDFINCL = -lnetcdf-cxx4 -lnetcdf
LIBS = -L/usr/local/lib

OPTIONS = $(OPTIMIZE) $(OPT) $(OPT_OPENMP)
CXXFLAGS = $(OPTIONS)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LIBS) $(NCDFINCL) -o $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC)