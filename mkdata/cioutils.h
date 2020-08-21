/*
 *  The maximum length of a file name
 */

#define MAXFILENAME 200

void initarray(void *ptr, int rank, int nx, int ny);
void initpgrid(void *ptr, int nxproc, int nyproc);

void checkandgetarguments(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *barrier, int size, int rank);
void getargs(int argc, char **argv, int *nx, int *ny, int *xprocs, int *yprocs, int *nxp, int *nyp, int *barrier, int size, int rank);

void createfilename(char *filename, char *basename, int nx, int ny);
void createindividualfilename(char *filename, char *basename, int nx, int ny, int nxp, int nyp, int rank);
