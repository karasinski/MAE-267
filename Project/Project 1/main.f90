program heat
  use clock
  use MainRoutines
  implicit none

  type (GridPoint), target, allocatable :: Points(:,:)
  type (GridCell),  target, allocatable :: Cells(:,:)
  integer :: step = 0

  ! Set up our grid size and allocate our arrays for our grid points and grid cells.
  call SetGridSize(101)
  allocate(Points(1:IMAX, 1:JMAX))
  allocate(Cells(1:IMAX-1, 1:JMAX-1))

  call initialization(Points, Cells)
  call start_clock()
  call solve(Points, Cells, step)
  call end_clock()
  call output(Points, step)

  ! Might as well be proper and cleanup before we leave.
  deallocate(Points)
  deallocate(Cells)
end program heat
