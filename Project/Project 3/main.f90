program heat
  use MainRoutines
  implicit none

  ! Block array for grid creation.
  type (BlockType), allocatable :: BlocksCollection(:)

  ! Block array for solver.
  type (BlockType), allocatable :: Blocks(:)
  integer :: step = 0

  ! Set up our grid size and allocate our arrays for our grid points and grid cells.
  call SetGridSize(101)
  call SetNumberOfBlocks(10,10)
  allocate(BlocksCollection(nBlocks))

  ! First we create our blocks and pack them with nodes.
  call initialize_grid(BlocksCollection)

  ! We then write a connectivity file.
  call write_configuration_file(BlocksCollection)

  ! We then write a grid file and initial temperature file.
  call plot3D(BlocksCollection, "i")

  ! We then initialize the solver.
  allocate(Blocks(1:nBlocks))
  call initialization(Blocks)

  ! Time our iterations until convergence.
  call start_clock()
  call solve(Blocks, step)
  call end_clock()

  ! Write final temperature distribution.
  call output(Blocks, step)
  call plot3D(Blocks, "f")

  ! Might as well be proper and cleanup before we leave.
  deallocate(BlocksCollection)
end program heat
