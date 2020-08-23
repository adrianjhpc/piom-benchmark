#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

#include "ioutils.h"

#include "arralloc.h"

#define NDIM 2

int main(int argc, char **argv)
{
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

  int istart, jstart;
  
  double starttime, endtime, totaltime;

  char inputfilename[MAXFILENAME], outputfilename[MAXFILENAME];

  MPI_Comm comm = MPI_COMM_WORLD;
  
  long int offset;
  int datasize;
  
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

      for (i=0; i < size; i++)
	{
	  printf("Process %2d has grid coordinates (%2d, %2d)\n",
		 i, pcoords[i][0], pcoords[i][1]);
	}
      printf("\n");
    }

  /*
   *  Initialise the array to a grey value
   */
  initarray(&x[0][0], nxp, nyp);

  createfilename(inputfilename, argv[1], nx, ny, nxp, nyp, rank);

  if(barrier){
    MPI_Barrier(comm);
  }

  starttime = MPI_Wtime();
  
  datasize = sizeof(float);
  offset = 0;

  for(i=0; i<nxp; i++){
    iochunkread (inputfilename, &x[i][0], nyp, offset);
    offset = offset + nyp*datasize;
  }

  endtime = MPI_Wtime();

  int stop = 0;

  float initial_data = rank * 0.5;

  for(i=0; i<nxp; i++){
    for(j=0; j<nyp; j++){
      if(x[i][j] != initial_data){
        printf("%d error %d %d %lf %lf\n",rank,i,j,x[i][j],initial_data);
        stop = 1;
      }
      initial_data = initial_data + 1.0;
      if(stop == 1){
        break;
      }
    }
    if(stop == 1){
      break;
    }  
  }
#ifdef DEBUG
  createfilename(outputfilename, "coutput", nx, ny, nxp, nyp, rank);
  iowrite(filename, &x[0][0], nxp*nyp);
#endif

  totaltime = endtime - starttime;

  dotimings(totaltime, rank, size);

  MPI_Finalize();

  free(x);
  free(pcoords);
}
