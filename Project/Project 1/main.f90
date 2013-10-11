program heat
  use constants
  use GridPointModule
  use GridCellModule

  integer :: i, j, max_i = 0, max_j = 0, step = 1

  type (GridPoint), pointer :: Point
  type (GridPoint), target, allocatable :: Points(:,:)

  type (GridCell), pointer :: Cell
  type (GridCell), target, allocatable :: Cells(:,:)

  real*8, pointer :: Temperature(:,:), tempTemperature(:,:)
  real*8 :: maxDiff = 0., residual = 999. !Arbitrary large number for residual.

  ! Set up our grid size, grid points, grid cells and our arrays.
  call SetGridSize(101)
  allocate(Points(1:IMAX, 1:JMAX))
  allocate(Cells(1:IMAX-1, 1:JMAX-1))

  !  Initialize grid.
  do i = 1, IMAX
    do j = 1, JMAX
      call initialize_points(Points(i,j), i, j)
    end do
  end do

  !  Initialize Cells.
  do i = 1, IMAX-1
    do j = 1, JMAX-1
      call initialize_cells(Cells(i,j), Points, i, j)
    end do
  end do

  !  Set up Dirichlet condition.
  do i = 1, IMAX-1
    j = 1
    call set_temperature(Cells(i,j), abs(cos(pi * Points(i,j)%xp)) + 1.)
    j = IMAX-1
    call set_temperature(Cells(i,j), 5. * (sin(pi * Points(i,j)%xp) + 1.))
  end do

  do j = 1, JMAX-1
    i = 1
    call set_temperature(Cells(i,j), 3. * Points(i,j)%yp + 2.)
    i = IMAX-1
    call set_temperature(Cells(i,j), 3. * Points(i,j)%yp + 2.)
  end do
  !  End Dirichlet condition.

  !  Set some useful pointers.
  Temperature => Cells%T
  tempTemperature => Cells%tempT
  !  End set up.

  !  Begin main loop, stop if we hit our mark or after 100,000 iterations.
  do while (residual >= .00001 .and. step <= 100000)
    Temperature = tempTemperature

    write(*,*), 'step = ', step

    do i = 2, size(Cells,1) - 1
      do j = 2, size(Cells,2) - 1
        call update_temperature(Cells, i, j)
      end do
    end do

    residual = maxval(abs(Temperature - tempTemperature))
    write(*, *), "residual ", residual

    step = step + 1
  end do
  !  End main loop.

  ! Let's find the last cell to change temperature and write some output.
  open (unit = 1, file = "data.dat")
  do i = 1, size(Cells,1)
    do j = 1, size(Cells,2)
      Cell => Cells(i,j)

      if (abs(Cell%T - Cell%tempT) > maxDiff) then
        maxDiff = abs(Cell%T - Cell%tempT)
        max_i = i
        max_j = j
      end if

      write (1,'(I5, 5X, I3, 5X, I3, 5X, F10.8)'), step, Cell%i, Cell%j, Cell%T
    end do
  end do

  write (*,*), "Max residual = ", maxDiff, "At i ", max_i, ", j ", max_j
  write (*,*), "Max temperature = ", maxval(Temperature), " Low temperature = ", minval(Temperature)

  close(1)
  ! End output.

end program heat
