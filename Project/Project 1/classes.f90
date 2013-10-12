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
    real*8 :: dTdx, dTdy, d2Td2x, d2Td2y
    real*8 :: Ayi, Axi, Ayj, Axj
    real*8 :: Ayi_half, Axi_half, Ayj_half, Axj_half
  end type GridCell

contains
  subroutine initialize_cells(c, points, i, j)
    save
    type (GridCell), intent(inout) :: c
    type (GridPoint), pointer :: p1, p2, p3, p4
    type (GridPoint), target :: points(1:IMAX, 1:JMAX)
    integer :: i, j

    ! These values come in handy later when referencing cells.
    c%i = i
    c%j = j

    ! Find the location of nearby points...
    p1 => Points(i,j)
    p2 => Points(i+1,j)
    p3 => Points(i+1,j+1)
    p4 => Points(i,j+1)

    ! ...to calculate the volume of each cell.
    c%V = abs( (p1%x *  p2%y - p1%y * p2%x) + &
      (p2%x *  p3%y - p2%y * p3%x) + &
      (p3%x *  p4%y - p3%y * p4%x) + &
      (p4%x *  p1%y - p4%y * p1%x) ) / 2.

    ! We set each cell to an initial temperature of 3.5, some
    ! cells will be overwritten when we declare boundary conditions.
    c%T = 3.5
    c%tempT = c%T

    ! These 'area's are used to do trapezoidal counter-clockwise
    ! integration to get the first derivatives in the x/y directions
    ! at the cell-center using Gauss's theorem.
    c%Ayi = p4%y - p1%y
    c%Axi = p4%x - p1%x
    c%Ayj = p2%y - p1%y
    c%Axj = p2%x - p1%x

    ! This sets the number of neighbors for later when we pass out
    ! temperatures.
    !    c%numberOfNeighbors = 4
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

  subroutine set_secondary_areas(c, i, j)
    save
    type (GridCell), intent(inout) :: c
    type (GridCell), target :: cells(1:IMAX, 1:JMAX)
    integer :: i, j

    ! These areas are used in the calculation of fluxes
    ! for the alternate distributive scheme second-
    ! derivative operator.
    c%Ayi_half = ( cells(i+1, j)%Ayi + cells(i, j)%Ayi   ) / 4.
    c%Axi_half = ( cells(i+1, j)%Axi + cells(i, j)%Axi   ) / 4.
    c%Ayj_half = ( cells(i, j)%Ayj   + cells(i, j+1)%Ayj ) / 4.
    c%Axj_half = ( cells(i, j)%Ayj   + cells(i, j+1)%Ayj ) / 4.
  end subroutine
end module GridCellModule

module UpdateTemperature
  use constants
  use GridPointModule
  use GridCellModule

contains
  subroutine first_derivative(p, c, i, j)
    save
    type (GridPoint) :: p(1:IMAX, 1:JMAX)
    type (GridCell)  :: c(1:IMAX-1, 1:JMAX-1)

    ! Trapezoidal counter-clockwise integration to get the first
    ! derivatives in the x/y directions at the cell-center using
    ! Gauss's theorem.

    !These Vs should be V(i+1/2, j+1/2), not clear if this is currently
    !correct or if I need to wiggle it.
    c(i,j)%dTdx = &
      ( (p(i+1, j)%T + p(i+1, j+1)%T ) * c(i+1, j)%Ayi - &
        (p(i, j)%T   + p(i, j+1)%T )   * c(i, j)%Ayi   - &
        (p(i, j+1)%T + p(i+1, j+1)%T ) * c(i, j+1)%Ayj + &
      (  p(i, j)%T   + p(i+1, j)%T )   * c(i, j)%Ayj     &
      ) / ( 2. * c(i,j)%V )

    c(i,j)%dTdy = &
      ( (p(i+1, j)%T + p(i+1, j+1)%T ) * c(i+1, j)%Axi - &
        (p(i, j)%T   + p(i, j+1)%T )   * c(i, j)%Axi   - &
        (p(i, j+1)%T + p(i+1, j+1)%T ) * c(i, j+1)%Axj + &
        (p(i, j)%T   + p(i+1, j)%T )   * c(i, j)%Axj     &
      ) / ( 2. * c(i,j)%V )
  end subroutine

  subroutine second_derivative(c, i, j)
    save
    type (GridCell) :: c(1:IMAX-1, 1:JMAX-1)

    ! Alternate distributive scheme second-derivative operator.
    ! Can definitely clean this up with some thought, but this essentially
    ! just updates the second derivative by adding the first times a constant
    ! during each time step.
    c(i,j+1)%d2Td2x   = c(i,j+1)%d2Td2x +   &
      (  c(i,j)%Axi_half + c(i,j)%Axj_half ) * c(i,j)%dTdx
    c(i+1,j+1)%d2Td2x = c(i+1,j+1)%d2Td2x + &
      ( -c(i,j)%Axi_half + c(i,j)%Axj_half ) * c(i,j)%dTdx
    c(i+1,j)%d2Td2x   = c(i+1,j)%d2Td2x +   &
      ( -c(i,j)%Axi_half - c(i,j)%Axj_half ) * c(i,j)%dTdx
    c(i,j)%d2Td2x     = c(i,j)%d2Td2x +     &
      (  c(i,j)%Axi_half - c(i,j)%Axj_half ) * c(i,j)%dTdx

    c(i,j+1)%d2Td2y   = c(i,j+1)%d2Td2y +   &
      ( -c(i,j)%Ayi_half - c(i,j)%Ayj_half ) * c(i,j)%dTdy
    c(i+1,j+1)%d2Td2y = c(i+1,j+1)%d2Td2y + &
      (  c(i,j)%Ayi_half - c(i,j)%Ayj_half ) * c(i,j)%dTdy
    c(i+1,j)%d2Td2y   = c(i+1,j)%d2Td2y +   &
      (  c(i,j)%Ayi_half + c(i,j)%Ayj_half ) * c(i,j)%dTdy
    c(i,j)%d2Td2y     = c(i,j)%d2Td2y +     &
      ( -c(i,j)%Ayi_half + c(i,j)%Ayj_half ) * c(i,j)%dTdy
  end subroutine

end module UpdateTemperature
