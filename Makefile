CC = g++
PRODUCTFLAGS=-Wall -pedantic
LINKFLAGS=-lfftw3 -lm
DEBUGFLAGS=-g
TARGET=pm

#Use this option to debug code.
#Make sure the code compiles without warning under this option before using 'make final'.
production: pm.cxx
	${CC} pm.cxx ${DEBUGFLAGS} ${LINKFLAGS} ${PRODUCTFLAGS} -o ${TARGET}


#Use this option to debug code further.
debug: pm.cxx
	valgrind -v --leak-check=yes --track-origins=yes ./pm


#Use this to compile the final product.
final: pm.cxx
	${CC} pm.cxx -o ${TARGET}

#Use this to delete compiled executable file.
clean:
	rm -rf ${TARGET}