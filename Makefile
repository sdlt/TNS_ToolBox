# Makefile

CC = gcc
CF = -O3 -fpic -Wall -pedantic -ffast-math -std=gnu11
CL = -lm -lgsl -lgslcblas -lfftw3 -fopenmp

objects1 = fftlog_PT biasing halofit TNSCorrections RealPower
objects2 = fftwlog TNS

all: rp tns libtns

$(objects1): %: %.c
	$(CC) $(CF) -c $< $(CL)

$(objects2): %: %.c
	$(CC) $(CF) -c $< $(CL)

rp: $(objects1)
	$(CC) $(CF) fftlog_PT.o biasing.o halofit.o TNSCorrections.o RealPower.o -o RealPower $(CL)

tns: $(objects2)
	$(CC) $(CF) fftwlog.o TNS.o -o TNS $(CL)

libtns: $(objects2)
	$(CC) $(CF) -shared fftwlog.o TNS.o -o libTNS.so $(CL)

clean:
	rm -f *~ *.o

purge:
	rm -f RealPower TNS libTNS.so *~ *.o
