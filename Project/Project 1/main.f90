module constants
  real*8, parameter, public :: k = 18.8, rho = 8000., c_p = 500.
  real*8, parameter, public :: pi = 3.141592654, rot = 30.*pi/180.
  integer :: IMAX, JMAX !101 or 501 (2 cases)
contains
  subroutine SetGridSize(length)
    IMAX = length
    JMAX = length
  end subroutine SetGridSize
end module

module GridPointModule
  use constants
  save

  public
  type GridPoint
    integer :: i, j
    real*8 :: x, xp, y, yp, T, tempT
  end type GridPoint

contains
  subroutine initialize(p, i, j)
    save
    type(GridPoint), intent(inout) :: p
    integer :: i, j

    p%i = i
    p%j = j

    p%xp = cos(0.5*pi*real((IMAX-i))/real((IMAX-1)))
    p%yp = cos(0.5*pi*real((JMAX-j))/real((JMAX-1)))

    p%x = p%xp*cos(rot)+(1.-p%yp)*sin(rot)
    p%y = p%yp*cos(rot)+(p%xp)*sin(rot)

    p%T = 3.5
    p%tempT = 3.5

  end subroutine initialize

  subroutine set_temperature(p, iT)
    save
    type(GridPoint), intent(inout) :: p
    real*8, optional :: iT

    p%tempT = iT

  end subroutine set_temperature

  subroutine update_temperature(p)
    ! subroutine should use update gridpoint temperature based on neighbors t
    save
    type(GridPoint), intent(inout) :: p

    p%T = p%T / 2.
    !    if (present(iT)) then
    !      p%T = iT
    !    else
    !      p%T = 3.5
    !    endif
  end subroutine

end module GridPointModule

program heat
  use constants
  use GridPointModule

  integer :: i, j, max_i = 0, max_j = 0, num, step = 1
  type (GridPoint), pointer :: Point, upperPoint, lowerPoint, leftPoint, rightPoint
  type (GridPoint), target, allocatable :: Points(:)
  real*8, pointer :: Temperature(:), tempTemperature(:)
  real*8 :: maxDiff = 0., temp, residual = 999.

  call SetGridSize(101)
  allocate(Points(1: IMAX * IMAX))
  Temperature => Points%T
  tempTemperature => Points%tempT

  !  Initialize grid.
  open (unit = 1, file = "data.dat")
  do i = 1, IMAX
    do j = 1, IMAX
      Point => Points( (i - 1) * IMAX + j)

      call initialize(Point, i, j)
      !      call set_temperature(Point, DFLOAT(INT(RAND(0)*100)))

      !      write (1,'(F10.5, 5X, F10.5, 5X, F10.5)'), Point%x, Point%y, Point%T
    end do
  end do
  !  End initialization.

  !  Set up Dirichlet condition.
  do i = 1, IMAX
    j = 1
    Point => Points( (i - 1) * IMAX + j )
    call set_temperature(Point, abs(cos(pi * Point%xp)) + 1.)

    j = IMAX
    Point => Points( (i- 1) * IMAX + j )
    call set_temperature(Point, 5. * (sin(pi * Point%xp) + 1.) )
  end do

  do j = 1, IMAX
    i = 1
    Point => Points( (i - 1) * IMAX + j )
    call set_temperature(Point, 3. * Point%yp + 2.)

    i = IMAX
    Point => Points( (i - 1) * IMAX + j)
    call set_temperature(Point, 3. * Point%yp + 2.)
  end do

  Temperature = tempTemperature
  !  End Dirichlet condition.

  !  Begin main loop.
  do while (residual > .00001 .and. step <= 15000)
    Temperature = tempTemperature

    write(*,*), 'step = ', step
    !    write (*, *), 'cell ', maxloc(Temperature), 'has max temp = ', maxval(Temperature)

    do i = 2, IMAX - 1
      do j = 2, IMAX - 1
        Point      => Points( (i - 1) * IMAX + (j) )

        temp = 0.
        num = 0

        !       Currently just summing up nearby nodes and averaging.
        !       I probably have the wrong labels here, but the idea should be correct
        !       and as long as I use consistent logic it shouldn't matter?

        if ( ((i - 1) * IMAX + (j - 1)) > 0 .and. ((i - 1) * IMAX + (j - 1)) < IMAX ** 2 ) then
          upperPoint => Points( (i - 1) * IMAX + (j - 1) )
          temp = temp + upperPoint%T
          num = num + 1
        end if

        if (( (i) * IMAX + (j) ) > 0 .and. ( (i) * IMAX + (j) ) < IMAX ** 2 ) then
          rightPoint => Points( (i) * IMAX + (j) )
          temp = temp + rightPoint%T
          num = num + 1
        end if

        if (( (i - 1) * IMAX + (j + 1) ) > 0 .and. ( (i - 1) * IMAX + (j + 1) < IMAX ** 2 )) then
          lowerPoint => Points( (i - 1) * IMAX + (j + 1) )
          temp = temp + lowerPoint%T
          num = num + 1
        end if

        if (( (i - 2) * IMAX + (j) ) > 0 .and. ( (i - 2) * IMAX + (j) < IMAX ** 2 )) then
          leftPoint  => Points( (i - 2) * IMAX + (j) )
          temp = temp + leftPoint%T
          num = num + 1
        end if

        temp = temp + Point%T
        num = num + 1
        temp = temp/dfloat(num)
        call set_temperature(Point, temp)

        !        Currently just divides by two.
        !        call update_temperature(Point)

        !        write (1,'(I5, 5X, F10.8, 5X, F10.8, 5X, F10.5)'), step, Point%x, Point%y, Point%tempT
      end do
    end do

    residual = maxval(abs(Temperature - tempTemperature))
    write(*, *), "residual ", residual

    step = step + 1
  end do

  do i = 2, IMAX - 1
    do j = 2, IMAX - 1
      Point      => Points( (i - 1) * IMAX + (j) )

      if (abs(Point%T - Point%tempT) > maxDiff) then
        maxDiff = abs(Point%T - Point%tempT)
        max_i = i
        max_j = j
      end if

      write (1,'(I5, 5X, F10.8, 5X, F10.8, 5X, F10.8)'), step, Point%x, Point%y, Point%T
    end do
  end do

  write (*,*), maxDiff, max_i, max_j
  !  End main loop.

  close(1)
end program heat
