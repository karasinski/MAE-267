program heat
  use clock
  use MainRoutines
  use plot3D_module
  implicit none

!  type (GridPoint), allocatable :: Points(:,:)
!  type (GridCell),  allocatable :: Cells(:,:)
  type (BlockType), allocatable :: BlocksCollection(:)
  type (BlockType), allocatable :: Blocks(:)
  integer :: step = 0

  ! Set up our grid size and allocate our arrays for our grid points and grid cells.
  call SetGridSize(101)
  call SetNumberOfBlocks(1,2)
  allocate(BlocksCollection(nBlocks))

  ! First we create our blocks and pack them with nodes.
  call initialize_grid(BlocksCollection)

  ! We then write a connectivity file.
  call write_configuration_file(BlocksCollection)

  ! We then write a grid file and initial temperature file.
  call plot3D(BlocksCollection)

  ! We then initialize the solver.
  allocate(Blocks(1:nBlocks))
  call initialization(Blocks)

  ! Time our iterations until convergence.
  call start_clock()
  call solve(Blocks, step)
  call end_clock()

  ! Write final temperature distribution.
!  call output(Blocks, step)
  call plot3D(Blocks)

  ! Might as well be proper and cleanup before we leave.
  deallocate(BlocksCollection)
end program heat
