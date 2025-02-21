output: main.o 
	g++ main.o -o simulation

main.o: main.cpp
	g++ -c main.cpp

clean:
	rm *.o simulation
