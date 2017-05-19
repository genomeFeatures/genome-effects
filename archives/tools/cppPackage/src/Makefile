
SHELL = /bin/sh

CC= g++

CPP_CORE = ./lib
CCFLAGS = -L/ubuntu/usr/lib -lgsl -lgslcblas -L${CPP_CORE} -lcpp_core

nnmf: ./lib/libcpp_core.a main.o nnmf.o
	${CC}  -w -O  main.o nnmf.o  ${CCFLAGS} -o nnmf


main.o: main.cpp
	${CC} -c -w -static -O main.cpp

nnmf.o: nnmf.cpp
	${CC} -c -w -static -O nnmf.cpp

clean:
	rm *.o
