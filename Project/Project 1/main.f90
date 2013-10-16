program heat
  use constants
  use clock
  use GridPointModule
  use GridCellModule
  use UpdateTemperature
  implicit none

  integer :: i, j, step = 1

  type (GridPoint), target, allocatable :: Points(:,:)

!  type (GridCell), pointer :: Cell
  type (GridCell), target, allocatable :: Cells(:,:)

  real*8, pointer :: Temperature(:,:), tempTemperature(:,:)
  real*8 :: residual = 1. ! Arbitrary initial residual.
  real*8, allocatable :: timestep(:,:)

  ! Set up our grid size, grid points, grid cells and our arrays.
  call SetGridSize(101)
  allocate(Points(1:IMAX, 1:JMAX))
  allocate(Cells(1:IMAX-1, 1:JMAX-1))
  allocate(timestep(2:IMAX-1, 2:JMAX-1))

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

  ! Set up Dirichlet condition.
  do j = 1, JMAX
    call set_temperature(Points(1,j), 3. * Points(1,j)%yp + 2.)
    call set_temperature(Points(IMAX,j), 3. * Points(IMAX,j)%yp + 2.)
  end do

  do i = 1, IMAX
    call set_temperature(Points(i,1), abs(cos(pi * Points(i,1)%xp)) + 1.)
    call set_temperature(Points(i,JMAX), 5. * (sin(pi * Points(i,JMAX)%xp) + 1.))
  end do
  ! End Dirichlet condition.

  ! Calculate timesteps and assign secondary volumes.
  do j = 2, JMAX - 1
    do i = 2, IMAX - 1
      Points(i, j)%timestep = ( 0.5 / (2 * alpha) ) * ( ( Cells(i, j)%V ** 2) / &
                 ( (Points(i+1, j)%x - Points(i, j)%x)**2 + ((Points(i, j+1)%y - Points(i, j)%y)**2 )) )
      Points(i, j)%Vol2 = ( ( Cells(i, j)%V + Cells(i - 1, j)%V + Cells(i, j - 1)%V + Cells(i - 1, j - 1)%V ) / 4.)
    end do
  end do

  !  Set some useful pointers.
  Temperature => Points(2:IMAX-1, 2:JMAX-1)%T
  tempTemperature => Points(2:IMAX-1, 2:JMAX-1)%tempT

!  End set up.

!  Begin main loop, stop if we hit our mark or after 1,000,000 iterations.
  call start_clock()
  do while (residual >= .00001 .and. step <= 1000000)

    Points%d2Td2x = 0.
    Points%d2Td2y = 0.

    do j = 1, JMAX - 1
      do i = 1, IMAX - 1
        call first_derivative(Points, Cells, i, j)
        call second_derivative(Points, Cells, i, j)
      end do
    end do

    do j = 2, JMAX - 1
      do i = 2, IMAX - 1
        Points(i, j)%tempT = Points(i,j)%timestep * alpha * ( Points(i, j)%d2Td2x + Points(i, j)%d2Td2y ) / Points(i, j)%Vol2
      end do
    end do

    Temperature = Temperature + tempTemperature
    residual = maxval(abs(tempTemperature))
!    write(*, *), step, "residual ", residual

    step = step + 1
  end do
  call end_clock()
  !  End main loop.

  ! Let's find the last cell to change temperature and write some output.
  open (unit = 1, file = "steady_state.dat")
  do i = 1, IMAX
    do j = 1, JMAX
      write (1,'(F10.7, 5X, F10.7, 5X, F10.7, I5, F10.7)'), Points(i,j)%x, Points(i,j)%y, Points(i,j)%T
    end do
  end do
  close (1)

  write (*,*), "IMAX/JMAX", IMAX, JMAX
  write (*,*), "steps", step
  write (*,*), "residual", maxval(tempTemperature)
  write (*,*), "ij", maxloc(tempTemperature)

  open (unit = 2, file = "other_info.dat")
  write (2,*), "For a ", IMAX, " by ", JMAX, "size grid, we ran for: "
  write (2,*), step, "steps"
  write (2,*), wall_time, "seconds"
  write (2,*)
  write (2,*), "Found max residual of ", maxval(tempTemperature)
  write (2,*), "At ij of ", maxloc(tempTemperature)
  close (2)
  ! End output.

  deallocate(Points)
  deallocate(Cells)

end program heat
