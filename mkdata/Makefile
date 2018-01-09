CC=cc
CFLAGS=-g

FC=ftn
FFLAGS=-g

DEPS=Makefile arralloc.h cioutils.h

SERIALOBJ=mkinfile.o arralloc.o

PARALLELCOBJ=parallelcmkinfile.o arralloc.o cioutils.o

PARALLELFOBJ=parallelfmkinfile.o fioutils.o

all: mkinfile parallelcinfile parallelfinfile
.PHONY:all

%.o: %.c $(DEPS) 
	$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.f90 Makefile
	$(FC) -c -o $@ $< $(FFLAGS)

mkinfile: $(SERIALOBJ) 
	$(CC) -o $@ $(SERIALOBJ) $(CFLAGS)

parallelcmkinfile: $(PARALLELCOBJ)
	$(CC) -o $@ $(PARALLELCOBJ) $(CFLAGS)

parallelfmkinfile: $(PARALLELFOBJ)
	$(FC) -o $@ $(PARALLELFOBJ) $(FFLAGS)

clean:
	rm -fr *.o mkinfile parallelcmkinfile parallelfmkinfile
