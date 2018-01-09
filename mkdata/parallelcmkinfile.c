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

  /*
   *  Variables needed for MPI-IO
   */

  int count, blocklength, stride, istart, jstart, floatsize;
  MPI_Datatype my_mpi_vector;

  MPI_File fh;
  MPI_Offset disp;
  MPI_Status status;

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

      for (i=0; i < size; i++)
	{
	  printf("Process %2d has grid coordinates (%2d, %2d)\n",
		 i, pcoords[i][0], pcoords[i][1]);
	}
      printf("\n");
    }

  /*
   *  Initialise the arrays to a grey value
   */

  initarray(&x[0][0]  , nxp, nyp);

  if(barrier){
    MPI_Barrier(comm);
  }

  /*
   *  Define the nxp x nyp vector for this distribution
   *  Note that it is the same for every process. To ensure that each
   *  process write different data from the file, they use different disps
   */

  count = nxp;
  blocklength = nyp;
  stride = ny;

  MPI_Type_vector(count, blocklength, stride, MPI_FLOAT, &my_mpi_vector);

  /*  Commit it before use */

  MPI_Type_commit(&my_mpi_vector);

  /*
   *  Construct name of input file
   */

  createfilename(filename, "cinput", nx, ny);


  /*
   *  Open the file for write only and attach to file handle fh
   *  No IO hints are passed since MPI_INFO_NULL is specified
   */

  if (MPI_File_open(comm, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY,
		    MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    {
      printf("Open error on rank %d\n", rank);
    }

  /*
   *  Set view for this process using the vector type with an appropriate
   *  value if disp (computed in BYTES!)
   */

  istart = pcoords[rank][0]*nxp;
  jstart = pcoords[rank][1]*nyp;

  MPI_Type_size(MPI_FLOAT, &floatsize);

  disp = istart;
  disp = disp*(MPI_Offset)ny;
  disp = disp + (MPI_Offset)jstart;
  disp = disp*(MPI_Offset)floatsize;

  /*
   *  Set view for this process using appropriate datatype
   */

  if (MPI_File_set_view(fh, disp, MPI_FLOAT, my_mpi_vector, "native",
			MPI_INFO_NULL) != MPI_SUCCESS)
    {
      printf("View error on rank %d\n", rank);
    }

  /*
   *  Write all the data for this process (ie nxp*nyp floats)
   */

  if (MPI_File_write_all(fh, &x[0][0], nxp*nyp, MPI_FLOAT, &status) != MPI_SUCCESS)
    {
      printf("Write error on rank %d\n", rank);
    }

  /*
   *  Close file
   */

  if (MPI_File_close(&fh) != MPI_SUCCESS)
    {
      printf("Close error on rank %d\n", rank);
    }


  MPI_Finalize();
  free(x);
  free(pcoords);

}
