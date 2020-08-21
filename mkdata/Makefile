CC=mpicc
CFLAGS=-g

FC=mpif90
FFLAGS=-g

DEPS=Makefile arralloc.h cioutils.h

SERIALOBJ=mkinfile.o arralloc.o

PARALLELCOBJ=parallelcmkinfile.o arralloc.o cioutils.o
PARALLELCOBJINDIVID=parallelcmkinfileindivid.o arralloc.o cioutils.o

PARALLELFOBJ=parallelfmkinfile.o fioutils.o

all: mkinfile parallelcmkinfile parallelfmkinfile parallelcmkinfileindivid
.PHONY:all

%.o: %.c $(DEPS) 
	$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.f90 Makefile
	$(FC) -c -o $@ $< $(FFLAGS)

mkinfile: $(SERIALOBJ) 
	$(CC) -o $@ $(SERIALOBJ) $(CFLAGS)

parallelcmkinfile: $(PARALLELCOBJ)
	$(CC) -o $@ $(PARALLELCOBJ) $(CFLAGS)

parallelcmkinfileindivid: $(PARALLELCOBJINDIVID)
	$(CC) -o $@ $(PARALLELCOBJINDIVID) $(CFLAGS)

parallelfmkinfile: $(PARALLELFOBJ)
	$(FC) -o $@ $(PARALLELFOBJ) $(FFLAGS)

clean:
	rm -fr *.o mkinfile parallelcmkinfile parallelcmkinfileindivid parallelfmkinfile
