module clock
  integer clock_start,clock_end,clock_max,clock_rate
  real*4 wall_time

contains
  subroutine start_clock()
    ! call system time to determine flow solver wall time
    call system_clock(count_max=clock_max,count_rate=clock_rate)
    call system_clock(clock_start)
  end subroutine start_clock

  subroutine end_clock()
    ! determine total wall time for solver
    call system_clock(clock_end)
    wall_time=float(clock_end-clock_start)/float(clock_rate)
    print*,'solver wall clock time (seconds)',wall_time
  end subroutine end_clock
end module

program heat
  use constants
  use clock
  use GridPointModule
  use GridCellModule
  use UpdateTemperature

  implicit none

  integer :: i, j, max_i = 0, max_j = 0, step = 1

  type (GridPoint), pointer :: Point
  type (GridPoint), target, allocatable :: Points(:,:)
  type (GridPoint), pointer :: innerPoints(:,:)

!  type (GridCell), pointer :: Cell
  type (GridCell), target, allocatable :: Cells(:,:)

  real*8, pointer :: Temperature(:,:), tempTemperature(:,:)
  real*8 :: timestep, maxDiff = 0., residual = 1. ! Arbitrary initial residual.
  real*8 :: tt = 10.

  call start_clock()

  ! Set up our grid size, grid points, grid cells and our arrays.
  call SetGridSize(501)
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
    call set_temperature(Points(1,j), 3. * Points(1,j)%yp + 2.)
    call set_temperature(Points(IMAX,j), 3. * Points(IMAX,j)%yp + 2.)
  end do

  do i = 1, IMAX
    call set_temperature(Points(i,1), abs(cos(pi * Points(i,1)%xp)) + 1.)
    call set_temperature(Points(i,JMAX), 5. * (sin(pi * Points(i,JMAX)%xp) + 1.))
  end do

!  End Dirichlet condition.

! do j = 1, JMAX
!    i = 1
!    call set_temperature(Points(i,j), tt)
!    i = IMAX
!    call set_temperature(Points(i,j), tt)
!  end do
!
!  do i = 1, IMAX
!    j = 1
!    call set_temperature(Points(i,j), tt)
!    j = IMAX
!    call set_temperature(Points(i,j), tt)
!  end do

  !  Set some useful pointers.
  innerPoints => Points(2:IMAX-1, 2:JMAX-1)
  Temperature => innerPoints%T
  tempTemperature => innerPoints%tempT
!  tempTemperature = 0.


  timestep = 17 * ( 0.5 / (2 * alpha) ) * ( ( Cells(IMAX-1, JMAX-1)%V ** 2) / &
  ( (Points(IMAX, JMAX)%x - Points(IMAX-1, JMAX)%x)**2 + ((Points(IMAX, JMAX)%y - Points(IMAX, JMAX-1)%y)**2 )) )
  !  End set up.

  !  Begin main loop, stop if we hit our mark or after 100,000 iterations.
  do while (residual >= .00001 .and. step <= 1000000)
!    write(*,*), 'step = ', step

!    do j = 2, JMAX - 1
!      do i = 2, IMAX - 1
!        call update_temperature(Points, i, j)
!      end do
!    end do
!    innerPoints%tempT = innerPoints%tempT / 4.

    do j = 1, JMAX - 1
      do i = 1, IMAX - 1
        call first_derivative(Points, Cells, i, j)
        call second_derivative(Points, Cells, i, j)
      end do
    end do
!
!    do j = 1, JMAX
!      do i = 1, IMAX
!              write(*,*), i, j, cells(i, j)%dTdx, cells(i, j)%dTdy
!end do
!end do

    do j = 2, JMAX - 1
      do i = 2, IMAX - 1
        Points(i, j)%tempT = timestep * alpha * ( Points(i, j)%d2Td2x + Points(i, j)%d2Td2y ) / &
                            ( ( Cells(i, j)%V + Cells(i - 1, j)%V + Cells(i, j - 1)%V + Cells(i - 1, j - 1)%V ) / 4.)
      end do
    end do
    Temperature = Temperature + tempTemperature
    Points%d2Td2x = 0.
    Points%d2Td2y = 0.

    residual = maxval(abs(tempTemperature))
    write(*, *), step, "residual ", residual

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

      write (1,'(I5, F10.7, 5X, F10.7, 5X, F10.7, I5, F10.7)'), step, Point%x, Point%y, Point%T
    end do
  end do

  write (*,*), "Max residual = ", maxDiff, "At i ", max_i, ", j ", max_j
  write (*,*), "Max temperature = ", maxval(Points%T), " Low temperature = ", minval(Points%T)

  close(1)
  ! End output.
  call end_clock()

end program heat
