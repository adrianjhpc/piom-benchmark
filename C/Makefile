CC=mpicc
CFLAGS=-g

DEPS=Makefile arralloc.h ioutils.h

OBJ=ioutils.o arralloc.o 
SSEND=mpiio_ssend
BCAST=mpiio_bcast
MASTER_VEC=mpiio_master_vector
MASTER_SUB=mpiio_master_subarray
READ_VEC=mpiio_read_vector
READ_SUB=mpiio_read_subarray
READALL_VEC=mpiio_readall_vector
READALL_SUB=mpiio_readall_subarray
SUB_READ_VEC=mpiio_sub_read_vector
SUB_READALL_VEC=mpiio_sub_readall_vector
SUB_READ_SUB=mpiio_sub_read_subarray
SUB_READALL_SUB=mpiio_sub_readall_subarray
INDIVID=mpiio_individual
INDIVID_FILES=individual_files

all: $(SSEND) $(BCAST) $(MASTER_VEC) $(MASTER_SUB) $(READ_VEC) $(READ_SUB) $(READALL_SUB) $(READALL_VEC) $(SUB_READ_VEC) $(SUB_READALL_VEC) $(SUB_READ_SUB) $(SUB_READALL_SUB) $(INDIVID) ${INDIVID_FILES}

.PHONY:all

%.o: %.c $(DEPS) 
	$(CC) -c -o $@ $< $(CFLAGS)

$(BCAST): $(OBJ) $(BCAST).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(SSEND): $(OBJ) $(SSEND).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(INDIVID): $(OBJ) $(INDIVID).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(INDIVID_FILES): $(OBJ) $(INDIVID_FILES).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(MASTER_VEC): $(OBJ) $(MASTER_VEC).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(MASTER_SUB): $(OBJ) $(MASTER_SUB).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(READ_VEC): $(OBJ) $(READ_VEC).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(READ_SUB): $(OBJ) $(READ_SUB).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(READALL_VEC): $(OBJ) $(READALL_VEC).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(READALL_SUB): $(OBJ) $(READALL_SUB).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(SUB_READ_VEC): $(OBJ) $(SUB_READ_VEC).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(SUB_READALL_VEC): $(OBJ) $(SUB_READALL_VEC).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(SUB_READ_SUB): $(OBJ) $(SUB_READ_SUB).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)

$(SUB_READALL_SUB): $(OBJ) $(SUB_READALL_SUB).o
	$(CC) -o $@ $@.o $(OBJ) $(CFLAGS)


clean:
	rm -fr *.o core* $(SSEND) $(BCAST) $(MASTER_VEC) $(MASTER_SUB) $(READ_VEC) $(READ_SUB) $(READALL_SUB) $(READALL_VEC) $(SUB_READ_VEC) $(SUB_READALL_VEC) $(SUB_READ_SUB) $(SUB_READALL_SUB) $(INDIVID) ${INDIVID_FILES}
