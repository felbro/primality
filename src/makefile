CC=g++
DBG=-Wall -pedantic
OPTI=-Ofast -march=native
CFLAG=-std=c++11 $(OPTI) $(DBG)
GMPFLAG=-lgmpxx -lgmp

main.out: main.cpp prime_lib.o
	$(CC) main.cpp prime_lib.o -static -o main.out $(CFLAG) $(GMPFLAG)

prime_lib.o: prime_lib.cpp prime_lib.h prime_1000.h
	$(CC) -c prime_lib.cpp $(CFLAG)

clean:
	rm -f *.out *.o *.gch *.s
