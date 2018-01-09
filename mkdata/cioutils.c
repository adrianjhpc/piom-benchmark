#include <stdio.h>
#include <stdlib.h>
#include "cioutils.h"
#include "mpi.h"


void checkandgetarguments(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *barrier, int size, int rank){

  if ( argc != 6 ){
    /* We print argv[0] assuming it is the program name */
    if(rank == 0){
      printf("usage: %s nx ny xprocs yprocs barrier\n", argv[0] );
      printf("This application expects you to provide the size of the input data set (nx*ny), the number of processes you want in each dimension (xproc and yproc), and whether to use a barrier before timing (barrier, 0 is no barrier, 1 uses a barrier)\n");
    }
    MPI_Finalize();
    exit(-1);
    
  }else{
    getargs(argc,argv,nx,ny,xprocs,yprocs,nxp,nyp,barrier,size,rank);
  }

}
void getargs(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *barrier, int size, int rank){

  *nx = atoi(argv[1]);
  *ny = atoi(argv[2]);
  *xprocs = atoi(argv[3]);
  *yprocs = atoi(argv[4]);
  *barrier = atoi(argv[5]);
  
  if((*xprocs)*(*yprocs) != size){
    if(rank == 0){
      printf("Number of processes being used (%d) does not match number requested in the two dimension (%d x %d)\n",size,*xprocs,*yprocs);
    }
    MPI_Finalize();
    exit(-1);
  }
  
  *nxp = (*nx)/(*xprocs);
  *nyp = (*ny)/(*yprocs);
  
  if((*nxp)*(*xprocs) != *nx){
    if(rank == 0){
      printf("Number of processes in x dimension (xprocs) does not exactly divide the x dimension of the input (nx).  Quitting\n");
    }
    MPI_Finalize();
    exit(-1);
  }
  
  if((*nyp)*(*yprocs) != *ny){
    if(rank == 0){
      printf("Number of processes in y dimension (yprocs) does not exactly divide the y dimension of the input (ny).  Quitting\n");
    }
    MPI_Finalize();
    exit(-1);
  }
  
  if(*barrier != 0 && *barrier != 1){
    if(rank == 0){
      printf("Barrier should be 0 or 1\n");
    }
    MPI_Finalize();
    exit(-1);
  }

  if(rank == 0){
    printf("Running on %d processes\n",size);
    printf("nx: %d ny: %d nxp: %d nyp: %d xprocs: %d yprocs: %d\n",*nx,*ny,*nxp,*nyp,*xprocs,*yprocs);
    printf("barrier: %d\n",*barrier);
  }

}


void createfilename(char *filename, char *basename, int nx, int ny)
{
  sprintf(filename, "%s%04dx%04d.dat", basename, nx, ny);
}

#define INITDATAVAL 0.5

void initarray(void *ptr, int nx, int ny)
{
  int i, j;

  float *data = (float *) ptr;

  for (i=0; i < nx*ny; i++)
    {
      data[i] = INITDATAVAL;
    }
}

#define NDIM 2

void initpgrid(void *ptr, int nxproc, int nyproc)
{
  MPI_Comm gridcomm;
  MPI_Comm comm = MPI_COMM_WORLD;

  int dims[NDIM];
  int periods[NDIM] = {0, 0};

  int reorder = 0;

  int *pcoords = (int *) ptr;

  int i;

  dims[0] = nxproc;
  dims[1] = nyproc;

  MPI_Cart_create(comm, NDIM, dims, periods, reorder, &gridcomm);

  for (i=0; i < nxproc*nyproc; i++)
    {
      MPI_Cart_coords(gridcomm, i, NDIM, pcoords+NDIM*i);
    }

  MPI_Comm_free(&gridcomm);

}
