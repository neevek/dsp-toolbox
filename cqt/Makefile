CPP=g++
CPPFLAGS=-c -Wall -std=c++11
LIBS=-lpng
BIN=cqtgen

all: ${BIN}

cqtgen: cqtgen.o png_image.o ConstantQTransform.o FourierTransform.o
	${CPP} cqtgen.o png_image.o ConstantQTransform.o FourierTransform.o ${LIBS} -o ${BIN}

cqtgen.o: cqtgen.cc
	${CPP} ${CPPFLAGS} -o cqtgen.o cqtgen.cc
png_image.o: ../util/png_image.cc
	${CPP} ${CPPFLAGS} -o png_image.o ../util/png_image.cc
ConstantQTransform.o: ConstantQTransform.cxx
	${CPP} ${CPPFLAGS} -o ConstantQTransform.o ConstantQTransform.cxx
FourierTransform.o: FourierTransform.cxx
	${CPP} ${CPPFLAGS} -o FourierTransform.o FourierTransform.cxx

clean: 
	rm -f *.o cqtgen 
