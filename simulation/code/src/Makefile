## Modify only this: used external libraries (nr or gnuplot)
_EXTOBJ := gnuplot_i.c nrutil.c 
PROGRAM_NAME := nonlocal.exe
###################################

## External libraries
# Numerical Recipes
NR_LIBS=../nr/recipes
NR_HEAD=../nr/include

# Gnuplot
GNDIR=../gnuplot_c

# Include
IDIR=include

## Compilation options
CC=gcc
CFLAGS=-g -Wall -Wstrict-prototypes  -I$(NR_HEAD)  -I$(GNDIR) -I$(IDIR)

# Included libraries
COMPLIBS=-lm -lgsl -lgslcblas 
COMPLIBS_PARALLEL=-lm -lgsl -lgslcblas -fopenmp

## Object directory
OBDIR=obj

# Objects
OBJ := $(patsubst %.c,$(OBDIR)/%.o,$(wildcard *.c))
GOBJ := $(subst $(GNDIR),$(OBDIR),$(patsubst %.c,%.o,$(wildcard $(GNDIR)/*.c)))
NOBJ := $(subst $(NR_LIBS),$(OBDIR),$(patsubst %.c,%.o,$(wildcard $(NR_LIBS)/*.c)))
EXTOBJ := $(patsubst %.c,$(OBDIR)/%.o,$(_EXTOBJ))

$(PROGRAM_NAME): $(OBJ) $(GOBJ) $(NOBJ)
	gcc -o $@ $(OBJ) $(EXTOBJ)  $(COMPLIBS) $(CFLAGS)

$(OBDIR)/%.o: %.c 
	$(CC) -c -o $@ $< $(COMPLIBS) $(CFLAGS)

$(OBDIR)/%.o: $(GNDIR)/%.c
	$(CC) -c -o $@ $< $(COMPLIBS) $(CFLAGS)

$(OBDIR)/%.o: $(NR_LIBS)/%.c
	$(CC) -c -o $@ $< $(COMPLIBS) $(CFLAGS)

clean:
	rm obj/* $(PROGRAM_NAME) *~
cleandata: 
	rm -r results

