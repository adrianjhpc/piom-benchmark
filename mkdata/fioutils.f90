subroutine createfilename(filename, basename, nx, ny)

  implicit none

  character*(*) :: filename, basename

  integer :: nx, ny

  write(filename, fmt='(i10.10,''x'',i10.10,''.dat'')') nx, ny

  filename = basename//filename

end


subroutine initarray(data, nx, ny)

  implicit none

  real, parameter :: initdataval = 0.5

  integer :: nx, ny

  real data(nx*ny)

  integer :: i

  do i = 1, nx*ny
    data(i) = initdataval
  end do

end

subroutine initpgrid(pcoords, nxproc, nyproc)

  use mpi

  implicit none

  integer, parameter :: ndim = 2

  integer :: nxproc, nyproc
  integer, dimension(ndim, nxproc*nyproc) :: pcoords

  integer, dimension(ndim) :: dims
  logical, dimension(ndim) :: periods
  integer, dimension(ndim) :: local_pcoords

  integer :: i, j, ierr
  integer :: comm = MPI_COMM_WORLD
  integer :: gridcomm

  logical :: reorder

  periods = (/ .false., .false. /)
  reorder = .false.

  dims(1) = nxproc
  dims(2) = nyproc

  call MPI_CART_CREATE(comm, ndim, dims, periods, reorder, gridcomm, ierr)

  do i = 1, nxproc*nyproc
    ! Copy the pcoords to a local array to get around some gfortran interface checking errors 
    ! when passing slices of the global pcoords array
    do j = 1, ndim
      local_pcoords(j) = pcoords(j, i)
    end do
    call MPI_CART_COORDS(gridcomm, i-1, ndim, local_pcoords, ierr)
  end do

  call MPI_COMM_FREE(gridcomm, ierr)

end subroutine initpgrid


subroutine checkandgetarguments(nx, ny, xprocs, yprocs, nxp, nyp, barrier, mpi_size, rank)

  use mpi

  implicit none

  integer, intent(in) :: mpi_size, rank
  integer, intent(out) :: nx, ny, xprocs, yprocs, nxp, nyp, barrier
  integer :: i, numargs, ierr
  character(len=32) :: arg

  numargs = command_argument_count()
  if(numargs .ne. 5) then
     if(rank .eq. 0) then
        write(*,*) "usage: nx ny xprocs yprocs barrier"
        write(*,*) "This application expects you to provid the size of the input data set (nx*ny),"  
        write(*,*) "the number of processes you want in each dimension (xprocs and yprocs), "
        write(*,*) "and whether to use a barrier before timing (barrier, 0 is no barrier, 1 uses a barrier)"
     end if
     call MPI_FINALIZE(ierr)
     stop
  else
     call getargs(nx, ny, xprocs, yprocs, nxp, nyp, barrier, mpi_size, rank)
  end if

end subroutine checkandgetarguments

subroutine getargs(nx, ny, xprocs, yprocs, nxp, nyp, barrier, mpi_size, rank)

  use mpi

  implicit none

  integer, intent(in) :: mpi_size, rank
  integer, intent(out) :: nx, ny, xprocs, yprocs, nxp, nyp, barrier

  integer :: i, numargs, ierr
  character(len=32) :: arg
 
  call getarg(1,arg)
  arg = trim(arg)
  read(arg, '(I10)') nx
  call getarg(2,arg)
  arg = trim(arg)
  read(arg, '(I10)') ny
  call getarg(3,arg)
  arg = trim(arg)
  read(arg, '(I10)') xprocs
  call getarg(4,arg)
  arg = trim(arg)
  read(arg, '(I10)') yprocs
  call getarg(5,arg)
  arg = trim(arg)
  read(arg, '(I2)') barrier
  if(xprocs*yprocs .ne. mpi_size) then
     if(rank .eq. 0) then
        write(*,*) 'The specified xprocs and yprocs assignment does not match the total number of processes being used.'
        write(*,*) 'xprocs is ',xprocs,' yprocs is ',yprocs,' but total processes used is ',mpi_size,' which does not match xprocs * yprocs'
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if
  nxp = nx/xprocs
  if(xprocs*nxp .ne. nx) then
     if(rank .eq. 0) then
        write(*,*) 'nx does not exactly divide by xprocs.  Stopping the code!'
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if
  nyp = ny/yprocs
  if(yprocs*nyp .ne. ny) then
     if(rank .eq. 0) then
        write(*,*) 'ny does not exactly divide by yprocs. Stopping the code!'
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if
  if(barrier .ne. 0 .and. barrier .ne. 1) then
     if(rank .eq. 0) then
        write(*,*) 'barrier should be 1 (use a barrier prior to starting the run) or 0 (do not use a barrier)'
     end if
     call MPI_FINALIZE(ierr)
     stop
  end if
  if(rank .eq. 0) then
     write(*,*) 'Running on',mpi_size,'processes'
     write(*,*) 'nx:',nx,'ny:',ny,'nxp:',nxp,'nyp:',nyp,'xprocs:',xprocs,'yprocs:',yprocs
     write(*,*) 'barrier:',barrier
  end if
  
end subroutine getargs
