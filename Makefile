CC = g++
PRODUCTFLAGS=-Wall -pedantic
DEBUGFLAGS=-g
TARGET=main

production: main.cxx
	${CC} main.cxx ${PRODUCTFLAGS} -o ${TARGET}
clean:
	rm -rf ${TARGET}
