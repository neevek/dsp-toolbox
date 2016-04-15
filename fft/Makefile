CC=gcc
CPP=g++
CFLAGS=-c -Wall
CPPFLAGS=-c -Wall -std=c++11
LIBS=-lpng
BIN=fftgen

all: ${BIN}

fftgen: fftgen.o png_image.o kiss_fftr.o kiss_fft.o
	${CPP} fftgen.o png_image.o kiss_fftr.o kiss_fft.o ${LIBS} -o ${BIN}

fftgen.o: fftgen.cc
	${CPP} ${CPPFLAGS} -o fftgen.o fftgen.cc
png_image.o: ../util/png_image.cc
	${CPP} ${CPPFLAGS} -o png_image.o ../util/png_image.cc
kiss_fftr.o: kiss_fftr.c
	${CC} ${CFLAGS} -o kiss_fftr.o kiss_fftr.c
kiss_fft.o: kiss_fft.c
	${CC} ${CFLAGS} -o kiss_fft.o kiss_fft.c

clean: 
	rm -f *.o fftgen 
