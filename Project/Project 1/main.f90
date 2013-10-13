program heat
  use constants
  use GridPointModule
  use GridCellModule
  use UpdateTemperature

  implicit none

  integer :: i, j, max_i = 0, max_j = 0, step = 1

  type (GridPoint), pointer :: Point
  type (GridPoint), target, allocatable :: Points(:,:)

  type (GridCell), pointer :: Cell
  type (GridCell), target, allocatable :: Cells(:,:)

  real*8, pointer :: Temperature(:,:), tempTemperature(:,:)
  real*8 :: maxDiff = 0., residual = 999. !Arbitrary large number for initial residual.

  ! Set up our grid size, grid points, grid cells and our arrays.
  call SetGridSize(101)
  allocate(Points(1:IMAX, 1:JMAX))
  allocate(Cells(1:IMAX-1, 1:JMAX-1))

  !  Initialize grid.
  do j = 1, JMAX
    do i = 1, IMAX
      call initialize_points(Points(i,j), i, j)
    end do
  end do

  !  Initialize Cells.
  do j = 1, JMAX-1
    do i = 1, IMAX-1
      call initialize_cells(Cells(i,j), Points, i, j)
    end do
  end do

  ! Set up secondary areas needed for integration.
  do j = 1, JMAX-1
    do i = 1, IMAX-1
      call set_secondary_areas(Cells(i,j), Points, i, j)
    end do
  end do

  !  Set up Dirichlet condition.
  do j = 1, JMAX
    i = 1
    call set_temperature(Points(i,j), 3. * Points(i,j)%yp + 2.)
    i = IMAX
    call set_temperature(Points(i,j), 3. * Points(i,j)%yp + 2.)
  end do

  do i = 1, IMAX
    j = 1
    call set_temperature(Points(i,j), abs(cos(pi * Points(i,j)%xp)) + 1.)
    j = IMAX
    call set_temperature(Points(i,j), 5. * (sin(pi * Points(i,j)%xp) + 1.))
  end do
  !  End Dirichlet condition.

  !  Set some useful pointers.
  Temperature => Points%T
  tempTemperature => Points%tempT
  !  End set up.

  !  Begin main loop, stop if we hit our mark or after 100,000 iterations.
  do while (residual >= .00001 .and. step <= 100000)
    Temperature = tempTemperature
!    points(2:JMAX-1, 2:IMAX-1)%T = points(2:JMAX-1, 2:IMAX-1)%tempT
    write(*,*), 'step = ', step

    do j = 1, JMAX - 1
      do i = 1, IMAX - 1
!        call update_temperature(Points, i, j)
        call first_derivative(Points, Cells, i, j)
        call second_derivative(Points, Cells, i, j)
      end do
    end do

!    points(2:JMAX-1, 2:IMAX-1)%tempT = points(2:JMAX-1, 2:IMAX-1)%tempT/4.


    points(2:JMAX-1, 2:IMAX-1)%tempT = points(2:JMAX-1, 2:IMAX-1)%d2Td2x + &
                                       points(2:JMAX-1, 2:IMAX-1)%d2Td2y
!    points(2:JMAX-1, 2:IMAX-1)%tempT = points(2:JMAX-1, 2:IMAX-1)%tempT/4.
    points(2:JMAX-1, 2:IMAX-1)%d2Td2x = 0.
    points(2:JMAX-1, 2:IMAX-1)%d2Td2y = 0.

    residual = maxval(abs(Temperature - tempTemperature))
    write(*, *), "residual ", residual

    step = step + 1
  end do
  !  End main loop.

  ! Let's find the last cell to change temperature and write some output.
  open (unit = 1, file = "data.dat")
  do i = 1, IMAX
    do j = 1, JMAX
      Point => Points(i,j)

      if (abs(Point%T - Point%tempT) > maxDiff) then
        maxDiff = abs(Point%T - Point%tempT)
        max_i = i
        max_j = j
      end if

      write (1,'(I5, 5X, F10.8, 5X, F10.8, 5X, F10.8)'), step, Point%x, Point%y, Point%T
    end do
  end do

  write (*,*), "Max residual = ", maxDiff, "At i ", max_i, ", j ", max_j
  write (*,*), "Max temperature = ", maxval(Temperature), " Low temperature = ", minval(Temperature)

  close(1)
  ! End output.

end program heat
