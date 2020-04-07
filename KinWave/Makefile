#Files
EXEC := KinWave
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

# clean directive 'make clean'
.PHONY: clean neat echo
clean: neat
	$(RM) $(OBJ) $(EXEC)
neat:
	$(RM) *~.*~
