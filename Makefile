all: simu

simu: simu.cc Makefile
	g++ -std=c++17 -O3 -Wall simu.cc -o simu

clean:
	rm -f simu
