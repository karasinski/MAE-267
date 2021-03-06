program heat
  use MainRoutines
  implicit none

  ! Block array for grid creation.
  type (BlockType), allocatable :: BlocksCollection(:)

  ! Block array for solver.
  type (BlockType), allocatable :: Blocks(:)

  ! Set up our grid size and allocate our arrays for our grid points and grid cells.
  allocate(BlocksCollection(nBlocks))

  ! First we create our blocks and pack them with nodes.
  call initialize_grid(BlocksCollection)

  ! We then write a connectivity file.
  call write_configuration_file(BlocksCollection)

  ! We then write a grid file and initial temperature file.
  call plot3D(BlocksCollection, "i")

  ! Deallocate our initialization array.
  deallocate(BlocksCollection)

  ! We then initialize the solver.
  allocate(Blocks(1:nBlocks))
  call initialization(Blocks)

  ! Time our iterations until convergence.
  call start_clock()
  call solve(Blocks)
  call end_clock()

  ! Write some results to file/screen.
  call output(Blocks)

  ! Write final temperature distribution.
  call plot3D(Blocks, "f")

  ! Might as well be proper and cleanup before we leave.
  deallocate(Blocks)
end program heat
