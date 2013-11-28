program heat
  use MainRoutines
  implicit none
  
  ! Block array for grid creation.
  type (BlockType), allocatable :: BlocksCollection(:)

  ! Proc array to contain blocks.
  type (Proc), allocatable :: Procs(:)

  ! Start MPI, find my id and total number of procs.
  call MPI_Init(ierror)
  call MPI_Comm_Rank(mpi_comm_world,MyID,ierror)
  call MPI_Comm_Size(mpi_comm_world,mpi_nprocs,ierror)

  ! Begin preprocessing.
  ! We use the first processor to write our initialization files.
  if (MyID == 0) then
    ! Allocate an initial array to hold all blocks and an additional
    ! array to store block information to go onto each proc.
    allocate(BlocksCollection(nBlocks))
    allocate(Procs(mpi_nprocs))

    ! First we create our blocks and pack them with nodes.
    call initialize_grid(BlocksCollection)

    ! Hand out blocks to processors.
    call distribute_blocks(BlocksCollection, Procs)

    ! We then write a connectivity file.
    call write_configuration_file(Procs)

    ! We then write a grid file and initial temperature file.
    call plotProcs(Procs)

    ! Deallocate our initialization array.
    deallocate(BlocksCollection, Procs)
  end if

  ! Hold all processors until we're done writing our initialization files.
  call MPI_Barrier(barrier, ierror)
  ! End preprocessing.

  ! Begin solver.
  ! We first read in the connectivity file.
  call read_configuration_file(Blocks)

  ! We then initialize the solver.
  call initialization

  ! Hold until we're all ready to start.
  call MPI_Barrier(barrier, ierror)

  ! Processor 0 times our iterations until convergence.
  if (MyID == 0) then
    call start_clock()
  end if

  write(*,*), 'Processor ', MyID, ' starting solver.'

  ! Each processor starts the solver.
  call solve(Blocks)
  call MPI_Barrier(barrier, ierror)

  if (MyID == 0) then
    call end_clock()
  end if
  ! End solver.

  ! Write some results to file/screen.
  !   call output(Blocks)

  ! Write final temperature distribution.
  call plot3D(Blocks)
  !
  ! Might as well be proper and cleanup before we leave.
  deallocate(Blocks)

  call MPI_Finalize(ierror)
end program heat
