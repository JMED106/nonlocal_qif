CC=gcc
CFLAGS=-g -Wall -Wstrict-prototypes

COMPLIBS=-lm -lgsl -lgslcblas 
COMPLIBS_PARALLEL=-lm -lgsl -lgslcblas -fopenmp
DEPS=qif.h utils.h  intg.h gnuplot_i.o theory.h file.h

LBITS := $(shell getconf LONG_BIT)
ifeq ($(LBITS),64)
OBJ=qif_64.o main_64.o utils_64.o intg_64.o gnuplot_i_64.o theory_64.o file_64.o
OBJP=qif_64p.o main_64p.o utils_64p.o intg_64p.o gnuplot_i_64p.o theory_64p.o file_64p.o
%_64.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(COMPLIBS)
%_64p.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(COMPLIBS_PARALLEL) 
nonlocal_FR_64.exe: $(OBJ)
	$(CC) -o $@ $^ $(COMPLIBS) 
	rm -r current && mkdir current && cp *.exe *.conf *.R ./current
	make clean
parallel: $(OBJP)
	$(CC) -o nonlocal_FR_64p.exe $^ $(COMPLIBS_PARALLEL)
	rm -r current && mkdir current && cp *.conf *.R ./current && mv *.exe ./current
	make clean

clean:
	rm *_64* 
else
OBJ=qif_32.o main_32.o utils_32.o  intg_32.o gnuplot_i_32.o theory_32.o file_32.o
OBJP=qif_32p.o main_32p.o utils_32p.o  intg_32p.o gnuplot_i_32.o theory_64p,o file_32p.o
%_32.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(COMPLIBS)
%_32p.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(COMPLIBS_PARALLEL)
nonlocal_FR_32.exe: $(OBJ)
	$(CC) -o $@ $^ $(COMPLIBS)
	rm -r current && mkdir current && cp *.exe *.conf *.R  ./current
	make clean
parallel: $(OBJP)
	$(CC) -o nonlocal_FR_32p.exe $^ $(COMPLIBS_PARALLEL)
	rm -r current && mkdir current && cp *.exe *.conf *.R  ./current
	make clean
clean:
	rm *_32*
endif

cleandata: 
	rm *.dat *.txt

