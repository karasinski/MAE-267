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
  subroutine initialize_points(p, i, j)
    save
    type(GridPoint), intent(inout) :: p
    integer :: i, j

    p%i = i
    p%j = j

    p%xp = cos(0.5*pi*dfloat((IMAX-i))/dfloat((IMAX-1)))
    p%yp = cos(0.5*pi*dfloat((JMAX-j))/dfloat((JMAX-1)))

    p%x = p%xp*cos(rot)+(1.-p%yp)*sin(rot)
    p%y = p%yp*cos(rot)+(p%xp)*sin(rot)

    p%T = 3.5
    p%tempT = p%T

  end subroutine initialize_points

subroutine set_temperature(p, T)
    save
    type(GridPoint), intent(inout) :: p
    real*8 :: T

    p%tempT = T

  end subroutine set_temperature

  subroutine update_temperature(points, i, j)
    save
    type(GridPoint), pointer :: p
    type (GridPoint), target :: points(1:IMAX, 1:JMAX)
    integer :: num

    p => points(i,j)
    p%tempT = 0
    p%tempT = p%tempT + p%T

    !    if (c%numberOfNeighbors == 4) then
    p%tempT = p%tempT + points(i,j-1)%T
    p%tempT = p%tempT + points(i+1,j)%T
    p%tempT = p%tempT + points(i,j+1)%T
    p%tempT = p%tempT + points(i-1,j)%T
    num = 5
    !    else
    !      num = 1
    !
    !      if ( j - 1 > 0 ) then
    !        c%tempT = c%tempT + Cells(i,j-1)%T
    !        num = num + 1
    !      end if
    !
    !      if ( i + 1 <= IMAX ) then
    !        c%tempT = c%tempT + Cells(i+1,j)%T
    !        num = num + 1
    !      end if
    !
    !      if ( j + 1 <= JMAX ) then
    !        c%tempT = c%tempT + Cells(i,j+1)%T
    !        num = num + 1
    !      end if
    !
    !      if ( i - 1 > 0 ) then
    !        c%tempT = c%tempT + Cells(i-1, j)%T
    !        num = num + 1
    !      end if
    !    end if

    p%tempT = p%tempT/dfloat(num)

  end subroutine update_temperature

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
    p3 => Points(i+1,j+1)
    p4 => Points(i,j+1)

    c%V = abs( (p1%x *  p2%y - p1%y * p2%x) + &
               (p2%x *  p3%y - p2%y * p3%x) + &
               (p3%x *  p4%y - p3%y * p4%x) + &
               (p4%x *  p1%y - p4%y * p1%x) ) / 2.

    c%T = 3.5
    c%tempT = c%T

    c%numberOfNeighbors = 4
    !    if ( j - 1 > 0 ) then
    !      c%numberOfNeighbors = c%numberOfNeighbors + 1
    !    end if
    !
    !    if ( i + 1 <= IMAX ) then
    !      c%numberOfNeighbors = c%numberOfNeighbors + 1
    !    end if
    !
    !    if ( j + 1 <= JMAX ) then
    !      c%numberOfNeighbors = c%numberOfNeighbors + 1
    !    end if
    !
    !    if ( i - 1 > 0 ) then
    !      c%numberOfNeighbors = c%numberOfNeighbors + 1
    !    end if

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

end module GridCellModule
