program heat
  use clock
  use MainRoutines
  use plot3D_module
  implicit none

  type (GridPoint), target, allocatable :: Points(:,:)
  type (GridCell),  target, allocatable :: Cells(:,:)
  type (GridPoint), allocatable :: Blocks(:,:,:,:)
  integer :: step = 0

  ! Set up our grid size and allocate our arrays for our grid points and grid cells.
  call SetGridSize(101)
  call SetNumberOfBlocks(10, 10)
  allocate(Points(1:IMAX, 1:JMAX))
  allocate(Cells(1:IMAX-1, 1:JMAX-1))
  allocate(Blocks(1:M, 1:N, 1:(1 + (IMAX - 1) / N), 1:(1 + (JMAX - 1) / M)))

  call initialization(Points, Cells)
  call make_blocks(Points, Blocks)
  !  call start_clock()
  !  call solve(Points, Cells, step)
  !  call end_clock()
  call output(Points, step)
  !  call make_blocks(Points, Blocks)
  call plot3D(Blocks)

  ! Might as well be proper and cleanup before we leave.
  deallocate(Points)
  deallocate(Cells)
end program heat
