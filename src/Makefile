CC=g++                                 #Name of compiler
FLAGS=-std=c++20 -w                    #Compiler flags
TARGET=SUPAASP 				           #Name of executable
OBJECTS=main.o calc_routines.o files.o #List of object files
LIBS=                                  #Libraries, not needed

#First target in Makefile is default
${TARGET}:${OBJECTS} #Make target from objects
	@echo "Linking..."
	${CC} ${FLAGS} ${OBJECTS} ${LIBS} -o ${TARGET}
	@make clean

main.o : main.cpp
	${CC} ${FLAGS} ${LIBS} -c main.cpp

calc_routines.o : calc_routines.cpp calc_routines.h
	${CC} ${FLAGS} ${LIBS} -c calc_routines.cpp

files.o : files.cpp files.h
	${CC} ${FLAGS} ${LIBS} -c files.cpp

clean: #No targets just run shell command to rm object files and emacs backups
	@rm -f *.o *~

cleantarget: #Delete the exectuable
	@rm -f ${TARGET}