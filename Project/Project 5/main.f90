PROGRAM heat
  USE MainRoutines
  IMPLICIT NONE

  ! Block array for grid creation.
  TYPE (BlockType), ALLOCATABLE :: BlocksCollection(:)

  ! Proc array to contain blocks.
  TYPE (Proc), POINTER :: Procs(:)

  ! Start MPI, find my id and total number of procs.
  CALL MPI_Init(ierror)
  CALL MPI_Comm_Rank(mpi_comm_world,MyID,ierror)
  CALL MPI_Comm_Size(mpi_comm_world,mpi_nprocs,ierror)

  ! Begin preprocessing.
  ! We use the first processor to write our initialization files.
  IF (MyID == 0) THEN
    ! Allocate an initial array to hold all blocks and an additional
    ! array to store block information to go onto each proc.
    ALLOCATE(BlocksCollection(nBlocks))
    ALLOCATE(Procs(mpi_nprocs))

    ! First we create our blocks and pack them with nodes.
    CALL initialize_grid(BlocksCollection)

    ! Hand out blocks to processors.
    CALL distribute_blocks(BlocksCollection, Procs)

    ! We then write a connectivity file.
    CALL write_configuration_file(Procs)

    ! We then write a grid file and initial temperature file.
    CALL plotProcs(Procs)

    ! Deallocate our initialization array.
    DEALLOCATE(BlocksCollection, Procs)
  END IF

  ! Hold all processors until we're done writing our initialization files.
  CALL MPI_Barrier(mpi_comm_world, ierror)
  ! End preprocessing.

  ! Begin solver.
  ! We first read in the connectivity file.
  CALL read_configuration_file(Blocks)

  ! We then initialize the solver.
  CALL initialization

  ! Hold until we're all ready to start.
  CALL MPI_Barrier(mpi_comm_world, ierror)
  WRITE(*,*), 'Processor ', MyID, ' starting solver.'

  ! Each processor starts the solver.
  CALL solve(Blocks)

  ! Write final temperature distribution.
  CALL plot3D(Blocks)

  ! Might as well be proper and cleanup before we leave.
  DEALLOCATE(Blocks)

  CALL MPI_Finalize(ierror)
END PROGRAM heat
