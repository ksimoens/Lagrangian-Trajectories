output: main.o Particle.o Vec.o
	g++ main.o Particle.o Vec.o -o simulation

Vec.o: Vec.cpp
	g++ -c Vec.cpp

Particle.o: Particle.cpp
	g++ -c Particle.cpp

main.o: main.cpp
	g++ -c main.cpp

clean:
	rm *.o simulation
