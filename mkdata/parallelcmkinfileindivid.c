#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "cioutils.h"

#include "arralloc.h"

#define NDIM 2

int main(int argc, char **argv){
  /*
   *  pcoords stores the grid positions of each process
   */

  int **pcoords;
  
  /*
   *  x contains the local data only
   */
  float **x;

  int rank, size;
  int i, j;

  int nx, ny, nxp, nyp, xprocs, yprocs, barrier;

  char filename[MAXFILENAME];

  size_t size_written;
  //  int count, blocklength, stride, istart, jstart, floatsize;

  FILE *fh;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  checkandgetarguments(argc, argv, &nx, &ny, &xprocs, &yprocs, &nxp, &nyp, &barrier, size, rank);

  pcoords = (int **)arralloc(sizeof(int), 2, size, NDIM);
  x = (float **)arralloc(sizeof(float),2,nxp,nyp);

  /*
   *  Work out the coordinates of all the processes in the grid and
   *  print them out
   */

  initpgrid(&pcoords[0][0], xprocs, yprocs);

  if (rank == 0)
    {
      printf("Running on %d process(es) in a %d x %d grid\n",
	     size, xprocs, yprocs);
      printf("\n");
    }

  /*
   *  Initialise the arrays to a grey value
   */

  initarray(&x[0][0], rank, nxp, nyp);

  if(barrier){
    MPI_Barrier(comm);
  }


  /*
   *  Construct name of input file
   */

  createindividualfilename(filename, "cinput", nx, ny, nxp, nyp, rank);

  /*
   * Open file
   */
  fh = fopen(filename,"wb");
  if(fh == NULL){
    printf("Error opening file on rank %d\n", rank);
  }
  
  /*
   *  Write all the data for this process (ie nxp*nyp floats)
   */

  for(i=0; i<nxp; i++){
      size_written = fwrite(&x[i][0], sizeof(float), nyp, fh); 
      if(size_written != nyp){
    	printf("Write error on rank %d\n", rank);
      }
  }

  /*
   *  Close file
   */

  fclose(fh);

  MPI_Finalize();
  free(x);
  free(pcoords);

}
