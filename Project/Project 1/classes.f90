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

end module GridPointModule

module GridCellModule
  use constants
  use GridPointModule
  save

  public
  type GridCell
    integer :: i, j, numberOfNeighbors
    real*8 :: V, T, tempT
  end type GridCell

contains
  subroutine initialize_cells(c, points, i, j)
    save
    type(GridCell), intent(inout) :: c
    type(GridPoint), pointer :: p1, p2, p3, p4
    type (GridPoint), target :: points(1:IMAX, 1:JMAX)
    integer :: i, j

    c%i = i
    c%j = j

    p1 => Points(i,j)
    p2 => Points(i+1,j)
    p3 => Points(i,j+1)
    p4 => Points(i+1,j+1)

    c%V = abs( (p1%x *  p2%y - p1%y * p2%x) + &
               (p2%x *  p3%y - p2%y * p3%x) + &
               (p3%x *  p4%y - p3%y * p4%x) + &
               (p4%x *  p1%y - p4%y * p1%x) ) / 2.

    c%T = 3.5
    c%tempT = c%T

    if ( j - 1 > 0 ) then
      c%numberOfNeighbors = c%numberOfNeighbors + 1
    end if

    if ( i + 1 <= IMAX ) then
      c%numberOfNeighbors = c%numberOfNeighbors + 1
    end if

    if ( j + 1 <= JMAX ) then
      c%numberOfNeighbors = c%numberOfNeighbors + 1
    end if

    if ( i - 1 > 0 ) then
      c%numberOfNeighbors = c%numberOfNeighbors + 1
    end if

  end subroutine initialize_cells

  !  subroutine find_neighbor_cells(c, cells, i, j)
  !    save
  !    type (GridCell), intent(inout) :: c
  !    type (GridCell), pointer :: cell
  !    type (GridCell), target :: cells(1:IMAX-1, 1:JMAX-1)
  !
  !    if (( i > 0 .and. i <= IMAX ) .and. ( j > 0 .and. j <= JMAX )) then
  !      c%neighborCells => cells(i:i,j:j)
  !    end if
  !  end subroutine find_neighbor_cells

  subroutine set_temperature(c, T)
    save
    type(GridCell), intent(inout) :: c
    real*8 :: T

    c%tempT = T

  end subroutine set_temperature

  subroutine update_temperature(Cells, i, j)
    save
    type(GridCell), pointer :: c
    type (GridCell), target :: cells(1:IMAX-1, 1:JMAX-1)
    integer :: num

    c => Cells(i,j)
    c%tempT = 0
    c%tempT = c%tempT + c%T

    if (c%numberOfNeighbors == 4) then
      c%tempT = c%tempT + Cells(i,j-1)%T
      c%tempT = c%tempT + Cells(i+1,j)%T
      c%tempT = c%tempT + Cells(i,j+1)%T
      c%tempT = c%tempT + Cells(i-1, j)%T
      num = 5
    else
      num = 1

      if ( j - 1 > 0 ) then
        c%tempT = c%tempT + Cells(i,j-1)%T
        num = num + 1
      end if

      if ( i + 1 <= IMAX ) then
        c%tempT = c%tempT + Cells(i+1,j)%T
        num = num + 1
      end if

      if ( j + 1 <= JMAX ) then
        c%tempT = c%tempT + Cells(i,j+1)%T
        num = num + 1
      end if

      if ( i - 1 > 0 ) then
        c%tempT = c%tempT + Cells(i-1, j)%T
        num = num + 1
      end if
    end if

    c%tempT = c%tempT/dfloat(num)

  end subroutine update_temperature

end module GridCellModule
