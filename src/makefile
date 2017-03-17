SOURCES=$(wildcard *.cpp)
CFLAG=-std=c++11
GMPFLAG=-lgmpxx -lgmp

build:
	g++ $(SOURCES) -o main $(CFLAG) $(GMPFLAG)

clean:
	rm main
	rm *.gch
	rm *.o
