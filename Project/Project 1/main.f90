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
    real*8 :: x, xp, y, yp, T
  end type GridPoint

contains
  subroutine initialize_points(p, i, j)
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

  end subroutine initialize_points

  subroutine set_temperature(p, iT)
    save
    type(GridPoint), intent(inout) :: p
    real*8, optional :: iT

    p%T = iT

  end subroutine set_temperature

end module GridPointModule

module GridCellModule
  use constants
  use GridPointModule
  save

  public
  type GridCell
    integer :: i, j
    real*8 :: V, T, tempT
  end type GridCell

contains
  subroutine initialize_cells(c, points, i, j)
    save
    type(GridCell), intent(inout) :: c
    type (GridPoint), pointer :: upperLeftPoint, upperRightPoint, lowerLeftPoint, lowerRightPoint
    type (GridPoint), target :: points(1:IMAX, 1:JMAX)
    integer :: i, j
    real*8 :: V

    upperLeftPoint => Points(i,j)
    upperRightPoint => Points(i+1,j)
    lowerLeftPoint => Points(i,j+1)
    lowerRightPoint => Points(i+1,j+1)

    c%i = i
    c%j = j

    c%V = (upperRightPoint%x - upperLeftPoint%x) * (lowerRightPoint%y - lowerLeftPoint%y)

    c%T = (upperRightPoint%T + upperLeftPoint%T + lowerRightPoint%T + lowerLeftPoint%T)/4.
    c%tempT = c%T

  end subroutine initialize_cells

!  subroutine find_neighbor_cells(c, i, j)
!    save
!    type(GridCell), intent(inout) :: c
!
!    if (( i > 0 .and. i <= IMAX ) .and. ( j > 0 .and. j <= JMAX )) then
!      c => Cells(i,j)
!      temp = temp + c%T
!      num = num + 1
!    end if
!  end subroutine

  subroutine update_temperature(c, T)
    save
    type(GridCell), intent(inout) :: c
    real*8 :: T

    c%tempT = T

  end subroutine update_temperature

end module GridCellModule

program heat
  use constants
  use GridPointModule
  use GridCellModule

  integer :: i, j, max_i = 0, max_j = 0, num, step = 1

  type (GridPoint), pointer :: Point
  type (GridPoint), target, allocatable :: Points(:,:)

  type (GridCell), pointer :: Cell, upperCell, lowerCell, leftCell, rightCell
  type (GridCell), target, allocatable :: Cells(:,:)

  real*8, pointer :: Temperature(:,:), tempTemperature(:,:)
  real*8 :: maxDiff = 0., temp, residual = 999.

  ! Set up our grid size, grid points, grid cells and our arrays.
  call SetGridSize(101)
  allocate(Points(1:IMAX, 1:JMAX))
  allocate(Cells(1:IMAX-1, 1:JMAX-1))

  !  Initialize grid.
  do i = 1, IMAX
    do j = 1, JMAX
      Point => Points(i,j)
      call initialize_points(Point, i, j)
    end do
  end do

  !  Set up Dirichlet condition.
  do i = 1, IMAX
    j = 1
    Point => Points(i,j)
    call set_temperature(Point, abs(cos(pi * Point%xp)) + 1.)

    j = IMAX
    Point => Points(i,j)
    call set_temperature(Point, 5. * (sin(pi * Point%xp) + 1.) )
  end do

  do j = 1, JMAX
    i = 1
    Point => Points(i,j)
    call set_temperature(Point, 3. * Point%yp + 2.)

    i = IMAX
    Point => Points(i,j)
    call set_temperature(Point, 3. * Point%yp + 2.)
  end do
  !  End Dirichlet condition.

  !  Initialize Cells.
  do i = 1, IMAX-1
    do j = 1, JMAX-1
      Cell => Cells(i,j)
      call initialize_cells(Cell, Points, i, j)
    end do
  end do

  !  Set some useful pointers.
  Temperature => Cells%T
  tempTemperature => Cells%tempT

  !  End set up.

  !  Begin main loop, stop if we hit our mark or after 100,000 iterations.
  timesteps: do while (residual >= .00001 .and. step <= 100000)
    Temperature = tempTemperature

    write(*,*), 'step = ', step

    i_loop: do i = 2, size(Cells,1) - 1
      j_loop: do j = 2, size(Cells,2) - 1
        Cell => Cells(i,j)

        temp = 0.
        num = 0

        if ( j - 1 > 0 ) then
          upperCell => Cells(i,j-1)
          temp = temp + upperCell%T
          num = num + 1
        end if

        if ( i + 1 <= IMAX ) then
          rightCell => Cells(i+1,j)
          temp = temp + rightCell%T
          num = num + 1
        end if

        if ( j + 1 <= JMAX ) then
          lowerCell => Cells(i,j+1)
          temp = temp + lowerCell%T
          num = num + 1
        end if

        if ( i - 1 > 0 ) then
          leftCell  => Cells(i-1, j)
          temp = temp + leftCell%T
          num = num + 1
        end if

        temp = temp + Cell%T
        num = num + 1
        temp = temp/dfloat(num)
        call update_temperature(Cell, temp)

      end do j_loop
    end do i_loop

    residual = maxval(abs(Temperature - tempTemperature))
    write(*, *), "residual ", residual

    step = step + 1
  end do timesteps
  !  End main loop.

  ! Write some output.
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
