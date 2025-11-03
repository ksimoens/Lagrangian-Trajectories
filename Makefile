#Total number of years to integrate.
OPT += -DNYEAR=5
#Total number of years over which to start.
OPT += -DNYEARSTART=1
#Number of days between starting moments.
OPT += -DDTSTART=20
#Number of particles.
OPT += -DNPART=1000
#Starting year.
OPT += -DYSTART=1993
#Diffusion constant.
OPT += -DD=100.0
#Timestep (seconds)
OPT += -DDT=86400.0
#Shape of velocity grid.
OPT += -DNLON=420
OPT += -DNLAT=240
#Extent of velocity grid.
OPT += -DLONMIN=-1.568614665
OPT += -DLONMAX=0.259617726
OPT += -DLONRES=0.004363323
OPT += -DMUMIN=0.1776416
OPT += -DMUMAX=1.7290554
#Particle properties.
OPT += -DSTOREPOS
OPT += -DSTOREVEL
OPT += -DCIRCULAR

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