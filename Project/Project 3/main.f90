program heat
  use clock
  use MainRoutines
  use plot3D_module
  implicit none

  type (GridPoint), allocatable :: Points(:,:)
  type (GridCell),  allocatable :: Cells(:,:)
  type (BlockType), allocatable :: Blocks(:,:)
  type (BlockType), allocatable :: BlocksCollection(:,:)
  integer :: step = 0

  ! Set up our grid size and allocate our arrays for our grid points and grid cells.
  call SetGridSize(101)
  call SetNumberOfBlocks(5,4)
!  allocate(Points(1:IMAX, 1:JMAX))
!  allocate(Cells(1:IMAX-1, 1:JMAX-1))
!  allocate(Blocks(1:M, 1:N))
  allocate(BlocksCollection(1:M,1:N))

!  call initialization(Points, Cells)
  call create_blocks(BlocksCollection)
  call initialize_block_grid(BlocksCollection)
  call initialize_block_temp(BlocksCollection)
!  call make_blocks(Points, BlocksCollection)
!  call initialize_blocks(Blocks, Points, Cells)
!  call start_clock()
!  call solve(Blocks, step)
!  call end_clock()
!  call make_blocks(Points, BlocksCollection)
!  call output(Blocks, step)
  call plot3D(BlocksCollection)

  ! Might as well be proper and cleanup before we leave.
!  deallocate(Points, Cells, Blocks)
end program heat
