CC=gcc
CFLAGS=-g -Wall -Wstrict-prototypes

COMPLIBS=-lm -lgsl -lgslcblas 
COMPLIBS_PARALLEL=-lm -lgsl -lgslcblas -fopenmp
DEPS=qif.h utils.h winfree.h theta.h intg.h gnuplot_i.o theory.h

LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
OBJ=qif_64.o prtyp_64.o utils_64.o winfree_64.o theta_64.o intg_64.o gnuplot_i_64.o theory_64.o
OBJP=qif_64p.o prtyp_64p.o utils_64p.o winfree_64p.o theta_64p.o intg_64p.o gnuplot_i_64p.o theory_64p.o
%_64.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(COMPLIBS)
%_64p.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(COMPLIBS_PARALLEL) 
prtype_64.exe: $(OBJ)
	$(CC) -o $@ $^ $(COMPLIBS) 
	rm -r current && mkdir current && cp *.exe *.conf *.R *.p ./current
parallel: $(OBJP)
	$(CC) -o prtype_64p.exe $^ $(COMPLIBS_PARALLEL)
	rm -r current && mkdir current && cp *.conf *.R *.p ./current && mv *.exe ./current
	make clean

clean:
	rm *_64* 
else
OBJ=qif_32.o prtyp_32.o utils_32.o winfree_32.o theta_32.o intg_32.o gnuplot_i_32.o theory_32.o
OBJP=qif_32p.o prtyp_32p.o utils_32p.o winfree_32p.o theta_32p.o intg_32p.o gnuplot_i_32.o theory_64p,o
%_32.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(COMPLIBS)
%_32p.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(COMPLIBS_PARALLEL)
prtype_32.exe: $(OBJ)
	$(CC) -o $@ $^ $(COMPLIBS)
	rm -r current && mkdir current && cp *.exe *.conf *.R *.p ./current
parallel: $(OBJP)
	$(CC) -o prtype_32p.exe $^ $(COMPLIBS_PARALLEL)
	rm -r current && mkdir current && cp *.exe *.conf *.R *.p ./current
	make clean
clean:
	rm *_32*
endif

cleandata: 
	rm *.dat *.txt

