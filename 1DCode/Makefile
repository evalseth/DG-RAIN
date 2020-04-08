# Files
EXEC := 1DSimulation
SRC  := $(wildcard *.c)
OBJ  := $(patsubst %.c, %.o, $(SRC))

# Options
CC 	:= /usr/local/bin/gcc-9 ##gcc
CFLAGS	:= -g -Wall -std=c99 -O2 #-pg
INCLUDE := -I/usr/local/include -I/usr/local/include/einspline -I/usr/local/include/gsl -I$$GSL_INC -I$$EINSPLINE_DIR/include/
LDFLAGS := -L/usr/local/lib #-L$$GSL_LIB -L$$EINSPLINE_DIR/lib/ 
LDLIBS 	:= -lm -lgsl -lgslcblas -leinspline #-pg


# Rules
$(EXEC) : $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o : %.c
	$(CC) $(CFLAGS) $(DEBUG) $(WD) $(INCLUDE) -lstdc++ -c $<

main.o: globals.h SimulationSteps.h
intialize.o: globals.h
store_mesh.o: mathfunctions.h globals.h 
time_evolution.o: globals.h mathfunctions.h oneTimeStep.h
computeL.o: mathfunctions.h globals.h
minmod.o: globals.h mathfunctions.h  
boundary_conditions.o: globals.h
numericalFlux1D.o: mathfunctions.h globals.h
1DInnerProducts.o: mathfunctions.h globals.h
constitutive_equations.o: mathfunctions.h
wetDry.o: globals.h

# clean directive 'make clean'
.PHONY: clean neat echo
clean: neat
	$(RM) $(OBJ) $(EXEC)
neat:
	$(RM) *~.*~
