#Files
EXEC := KinVenant
SRC	 := $(wildcard *.c)
OBJ  := $(patsubst %.c, %.o, $(SRC))

# Options
CC 			:= gcc
CFLAGS 	:= -g -Wall -std=c99
INCLUDE := -I$$GSL_INC -I/workspace/prapti/EINSPLINE/include/
LDFLAGS := -L$$GSL_LIB -L/workspace/prapti/EINSPLINE/lib/
LDLIBS 	:= -lm -lgsl -lgslcblas -leinspline

# Rules
$(EXEC) : $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

main.o: Channels_KinRoutes.h	Simulation_Steps.h
basis_function_routines.o: Channels_KinRoutes.h Globals.h
computeL.o: Channels_KinRoutes.h Math_Functions.h Globals.h Constitutive_Equations.h Numerical_Flux.h
constitutive_equations.o: Math_Functions.h
initialize.o: Channels_KinRoutes.h Globals.h Math_Functions.h
minmod.o: Math_Functions.h Constitutive_Equations.h Globals.h Channels_KinRoutes.h
numerical_flux.o: Math_Functions.h Constitutive_Equations.h Channels_KinRoutes.h
time_evolution.o: Channels_KinRoutes.h Globals.h Math_Functions.h One_Time_Step.h Constitutive_Equations.h


# clean directive 'make clean'
.PHONY: clean neat echo
clean: neat
	$(RM) $(OBJ) $(EXEC)
neat:
	$(RM) *~.*~
