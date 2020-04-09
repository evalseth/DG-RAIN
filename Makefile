# Files
EXEC := KinField 
SRC  := $(wildcard *.c)
OBJ  := $(patsubst %.c, %.o, $(SRC))

# Options
CC 	:= /usr/local/bin/gcc-9 
CFLAGS	:= -g -Wall -std=c99 
INCLUDE := -I/usr/local/include -I$$GSL_INC -I$$PYTHON_INC
LDFLAGS := -L/usr/local/lib -L$$GSL_LIB -L$$PYTHON_LIB
LDLIBS 	:= -lgsl -lgslcblas -lpython3.7 #-lm


# Rules
$(EXEC) : $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o : %.c
	$(CC) $(CFLAGS) $(DEBUG) $(WD) $(INCLUDE) -lstdc++ -c $<

main.o: Watershed.h SimulationSteps.h Globals.h
process_watershed.o: Globals.h
computeL.c: Watershed.h Globals.h Math_Functions.h Constitutive_Equations.h Numerical_Flux.h Infiltration.h
Infiltration.c: Infiltration.h

# clean directive 'make clean'
.PHONY: clean neat echo
clean: neat
	$(RM) $(OBJ) $(EXEC)
neat:
	$(RM) *~.*~
