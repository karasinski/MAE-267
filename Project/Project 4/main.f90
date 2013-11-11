program heat
  use MainRoutines
  implicit none

  ! Block array for grid creation.
  type (BlockType), allocatable :: BlocksCollection(:)

  ! Proc array to contain blocks.
  type (Proc), allocatable :: Procs(:)

  ! Block array for solver.
  type (BlockType), allocatable :: Blocks(:)
  integer :: step = 0

  ! Set up our grid size, set number of blocks and processors.
  call SetGridSize(101)
  call SetNumberOfBlocks(10,10)
  call SetNumberOfProcs(4)

  ! Allocate an initial array to hold all blocks and an additional
  ! array to store block information to go onto each proc.
  allocate(BlocksCollection(nBlocks))
  allocate(Procs(nProcs))

  ! First we create our blocks and pack them with nodes.
  call initialize_grid(BlocksCollection)

  ! Hand out blocks to processors.
  call distribute_blocks(BlocksCollection, Procs, 2)

  ! We then write a connectivity file.
  call write_configuration_file(Procs)

  ! We then write a grid file and initial temperature file.
  call plotProcs(Procs)

  ! Deallocate our initialization array.
  deallocate(BlocksCollection, Procs)
!
!  ! We then initialize the solver.
!  allocate(Blocks(1:nBlocks))
!  call initialization(Blocks)
!
!  ! Time our iterations until convergence.
!  call start_clock()
!  call solve(Blocks, step)
!  call end_clock()
!
!  ! Write some results to file/screen.
!  call output(Blocks, step)
!
!  ! Write final temperature distribution.
!  call plot3D(Blocks, "f")
!
!  ! Might as well be proper and cleanup before we leave.
!  deallocate(Blocks)
end program heat
