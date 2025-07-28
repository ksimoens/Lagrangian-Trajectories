#Total number of years to integrate.
OPT += -DNYEAR=25
#Total number of years over which to start.
OPT += -DNYEARSTART=1
#Number of particles.
OPT += -DNPART=100
#Starting year.
OPT += -DYSTART=1993
#Diffusion constant.
OPT += -DD=100.0
#Timestep (seconds)
OPT += -DDT=86400.0
#Shape of velocity grid.
OPT += -DNLON=368
OPT += -DNLAT=200
#Extent of velocity grid.
OPT += -DLONMIN=-1.428988325
OPT += -DLONMAX=0.172351264
OPT += -DLONRES=0.004363323
OPT += -DVELRES=0.25
OPT += -DLATMIN=20.125
OPT += -DLATMAX=69.875
OPT += -DMUMIN=0.3587011
OPT += -DMUMAX=1.7290554

#Extent of starting location.
OPT += -DCELLLONMIN=-5563.64787134093
OPT += -DCELLLATMIN=4590.6176765106

#Variables for transformation
OPT += -DA=6378.140
OPT += -DB=6356.7500
OPT += -DE2=0.0066960376522835485
OPT += -DE1=0.0016796375940427949
OPT += -DR=6367445.0

OBJS = main.o globParams.o Grid.o Particle.o Vec.o
CFILES = main.cpp globParams.cpp Grid.cpp Particle.cpp Vec.cpp
OPT_OPENMP = -fopenmp
EXEC = main
CXX = g++
OPTIMIZE = -O3 -Wall -g
INCL = Makefile globParams.h Grid.h Particle.h Vec.h

OPTIONS = $(OPTIMIZE) $(OPT) $(OPT_OPENMP)
CXXFLAGS = $(OPTIONS)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC)