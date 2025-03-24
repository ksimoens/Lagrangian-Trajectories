output: main.o Particle.o Vec.o Grid.o globParams.o
	g++ main.o Particle.o Vec.o Grid.o globParams.o -o simulation

Vec.o: Vec.cpp
	g++ -c Vec.cpp

Particle.o: Particle.cpp
	g++ -c Particle.cpp

Grid.o: Grid.cpp
	g++ -c Grid.cpp

globParams.o: globParams.cpp
	g++ -c globParams.cpp

main.o: main.cpp
	g++ -c main.cpp

clean:
	rm *.o simulation
