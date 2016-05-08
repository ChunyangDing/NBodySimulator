CC = g++
PRODUCTFLAGS=-Wall -pedantic
LINKFLAGS=-lfftw3 -lm
DEBUGFLAGS=-g
TARGET=main

#Use this option to debug code.
#Make sure the code compiles without warning under this option before using 'make final'.
production: main.cxx
	${CC} main.cxx ${DEBUGFLAGS} ${LINKFLAGS} ${PRODUCTFLAGS} -o ${TARGET}

#Use this to compile the final product.
final: main.cxx
	${CC} main.cxx -o ${TARGET}
#Use this to delete compiled executable file.
clean:
	rm -rf ${TARGET}
