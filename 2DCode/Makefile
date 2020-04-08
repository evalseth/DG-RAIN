# Files
EXEC := Simulation
SRC  := $(wildcard *.c)
OBJ  := $(patsubst %.c, %.o, $(SRC))

# Options
CC 	:= /usr/local/bin/gcc-9 ##gcc 
CFLAGS	:= -g -Wall -std=c99 -O2 #-DDEBUG 
INCLUDE := -I/usr/local/include -I/usr/local/include/gsl  -I$$GSL_INC 
LDFLAGS := -L/usr/local/lib #-L$$GSL_LIB
LDLIBS 	:= -lm -lgsl -lgslcblas  


# Rules
$(EXEC) : $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o : %.c
	$(CC) $(CFLAGS) $(DEBUG) $(WD) $(INCLUDE) -lstdc++ -c $<

main.o: SimulationSteps.h
extract_mesh_attributes.o: MeshAttributes.h mathfunctions.h
intialize.o: MeshAttributes.h
time_evolution.o: mathfunctions.h oneTimeStep.h
compute2DL.o: mathfunctions.h MeshAttributes.h
SlopeLimiter2D.o: mathfunctions.h MeshAttributes.h
boundary_conditions.o: MeshAttributes.h
extractXAxisEdges.o: MeshAttributes.h

# clean directive 'make clean'
.PHONY: clean neat echo
clean: neat
	$(RM) $(OBJ) $(EXEC)
neat:
	$(RM) *~.*~
