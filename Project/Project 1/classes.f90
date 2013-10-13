module constants
  implicit none
  real*8, parameter, public :: k = 18.8, rho = 8000., c_p = 500.
  real*8, parameter, public :: pi = 3.141592654, rot = 30.*pi/180.
  integer :: IMAX, JMAX !101 or 501 (2 cases)
contains
  subroutine SetGridSize(length)
    integer :: length
    IMAX = length
    JMAX = length
  end subroutine SetGridSize
end module

module GridPointModule
  use constants
  implicit none
  save

  public
  type GridPoint
    integer :: i, j
    real*8 :: x, xp, y, yp, T, tempT, d2Td2x, d2Td2y
  end type GridPoint

contains
  subroutine initialize_points(p, i, j)
    save
    type(GridPoint), intent(inout) :: p
    integer :: i, j

    ! Set the initial locations and initial temperature of
    ! each grid point.
    p%i = i
    p%j = j

    p%xp = cos(0.5*pi*dfloat((IMAX-i))/dfloat((IMAX-1)))
    p%yp = cos(0.5*pi*dfloat((JMAX-j))/dfloat((JMAX-1)))

    p%x = p%xp*cos(rot)+(1.-p%yp)*sin(rot)
    p%y = p%yp*cos(rot)+(p%xp)*sin(rot)

    p%T = 3.5
    p%tempT = p%T

    p%d2Td2x = 0.
    p%d2Td2y = 0.

  end subroutine initialize_points

  subroutine set_temperature(p, T)
    save
    type(GridPoint), intent(inout) :: p
    real*8 :: T

    ! Mostly unnecessary subroutine to set initial temperatures.
    p%tempT = T
    p%T = T

  end subroutine set_temperature

  subroutine update_temperature(points, i, j)
    save
    type(GridPoint), pointer :: p
    type (GridPoint), target :: points(1:IMAX, 1:JMAX)
    integer :: i, j

    ! This is my 'fake' method to transfer heat: sum the
    ! nearby points and divide by five. This is a temporary
    ! function to help test display data until I get the
    ! real derivatives working.
    p => points(i, j)

    p%tempT = points(i, j-1)%T + points(i+1, j)%T + points(i, j+1)%T + points(i-1, j)%T

  end subroutine update_temperature

end module GridPointModule

module GridCellModule
  use constants
  use GridPointModule

  implicit none
  save

  public
  type GridCell
    integer :: i, j, numberOfNeighbors
    real*8 :: V, T, tempT
    real*8 :: dTdx, dTdy
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
    p1 => Points(i, j)
    p2 => Points(i+1, j)
    p3 => Points(i+1, j+1)
    p4 => Points(i, j+1)

    ! ...to calculate the volume of each cell.
    c%V = abs( &
      (p1%x *  p2%y - p1%y * p2%x) + &
      (p2%x *  p3%y - p2%y * p3%x) + &
      (p3%x *  p4%y - p3%y * p4%x) + &
      (p4%x *  p1%y - p4%y * p1%x)   &
      ) / 2.

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

    c%dTdx = 0
    c%dTdy = 0

  end subroutine initialize_cells

  subroutine set_secondary_areas(c, p, i, j)
    save
    type (GridCell), intent(inout) :: c
    type (GridPoint), target :: p(1:IMAX, 1:JMAX)
    integer :: i, j
    real*8 :: Ayi, Axi, Ayj, Axj

    ! These areas are used in the calculation of fluxes
    ! for the alternate distributive scheme second-
    ! derivative operator.
    Ayi(i,j) = p(i, j+1)%y - p(i, j)%y
    Axi(i,j) = p(i, j+1)%x - p(i, j)%x
    Ayj(i,j) = p(i+1, j)%y - p(i, j)%y
    Axj(i,j) = p(i+1, j)%x - p(i, j)%x

    ! No longer using the cells values due to boundary problems... to be thought about.
    c%Ayi_half = ( Ayi(i+1, j) + Ayi(i, j)   ) / 4.
    c%Axi_half = ( Axi(i+1, j) + Axi(i, j)   ) / 4.
    c%Ayj_half = ( Ayj(i, j)   + Ayj(i, j+1) ) / 4.
    c%Axj_half = ( Ayj(i, j)   + Ayj(i, j+1) ) / 4.
  end subroutine
end module GridCellModule

module UpdateTemperature
  use constants
  use GridPointModule
  use GridCellModule

  implicit none

contains
  subroutine first_derivative(p, c, i, j)
    save
    type (GridPoint) :: p(1:IMAX, 1:JMAX)
    type (GridCell), intent(inout)  :: c(1:IMAX-1, 1:JMAX-1)
    real*8 :: Ayi, Axi, Ayj, Axj
    integer :: i, j

    ! Trapezoidal counter-clockwise integration to get the first
    ! derivatives in the x/y directions at the cell-center using
    ! Gauss's theorem.

    ! This is also a mess, so I should cleverly clean it up somehow.

    ! No longer using the cells values due to boundary problems... to be thought about.
    Ayi(i,j) = p(i, j+1)%y - p(i, j)%y
    Axi(i,j) = p(i, j+1)%x - p(i, j)%x
    Ayj(i,j) = p(i+1, j)%y - p(i, j)%y
    Axj(i,j) = p(i+1, j)%x - p(i, j)%x

    c(i, j)%dTdx = + &
      ( ( p(i+1, j)%T + p(i+1, j+1)%T ) * Ayi(i+1, j) - &
        ( p(i,   j)%T + p(i,   j+1)%T ) * Ayi(i,   j) - &
        ( p(i, j+1)%T + p(i+1, j+1)%T ) * Ayj(i, j+1) + &
        ( p(i,   j)%T + p(i+1,   j)%T ) * Ayj(i,   j)   &
      ) / ( 2. * c(i, j)%V )

    c(i, j)%dTdy = - &
      ( ( p(i+1, j)%T + p(i+1, j+1)%T ) * Axi(i+1, j) - &
        ( p(i,   j)%T + p(i,   j+1)%T ) * Axi(i,   j) - &
        ( p(i, j+1)%T + p(i+1, j+1)%T ) * Axj(i, j+1) + &
        ( p(i,   j)%T + p(i+1,   j)%T ) * Axj(i,   j)   &
      ) / ( 2. * c(i ,j)%V )
  end subroutine

  subroutine second_derivative(p, c, i, j)
    save
    type (GridPoint) :: p(1:IMAX, 1:JMAX)
    type (GridCell) :: c(1:IMAX-1, 1:JMAX-1)
    integer :: i, j

    ! Alternate distributive scheme second-derivative operator.
    ! Can definitely clean this up with some thought, but this essentially
    ! just updates the second derivative by adding the first times a constant
    ! during each time step.
!    p(i, j+1)%tempT   = p(i, j+1)%tempT +   &
!      ((  c(i, j)%Axi_half + c(i, j)%Axj_half ) * c(i, j)%dTdx + &
!       ( -c(i, j)%Ayi_half - c(i, j)%Ayj_half ) * c(i, j)%dTdy &
!      )
!
!    p(i+1, j+1)%tempT = p(i+1, j+1)%tempT + &
!      (( -c(i, j)%Axi_half + c(i, j)%Axj_half ) * c(i, j)%dTdx + &
!       (  c(i, j)%Ayi_half - c(i, j)%Ayj_half ) * c(i, j)%dTdy &
!      )
!
!    p(i+1, j)%tempT   = p(i+1, j)%tempT +   &
!      (( -c(i, j)%Axi_half - c(i, j)%Axj_half ) * c(i, j)%dTdx + &
!       (  c(i, j)%Ayi_half + c(i, j)%Ayj_half ) * c(i, j)%dTdy &
!      )
!
!    p(i, j)%tempT     = p(i, j)%tempT +     &
!      ((  c(i, j)%Axi_half - c(i, j)%Axj_half ) * c(i, j)%dTdx + &
!       ( -c(i, j)%Ayi_half + c(i, j)%Ayj_half ) * c(i, j)%dTdy &
!      )

    p(i, j+1)%d2Td2x   = p(i, j+1)%d2Td2x +   &
      (  c(i, j)%Ayi_half + c(i, j)%Ayj_half ) * c(i, j)%dTdx

    p(i+1, j+1)%d2Td2x = p(i+1, j+1)%d2Td2x + &
      ( -c(i, j)%Ayi_half + c(i, j)%Ayj_half ) * c(i, j)%dTdx

    p(i+1, j)%d2Td2x   = p(i+1, j)%d2Td2x +   &
      ( -c(i, j)%Ayi_half - c(i, j)%Ayj_half ) * c(i, j)%dTdx

    p(i, j)%d2Td2x     = p(i, j)%d2Td2x +     &
      (  c(i, j)%Ayi_half - c(i, j)%Ayj_half ) * c(i, j)%dTdx


    p(i,   j+1)%d2Td2y = p(i,   j+1)%d2Td2y +   &
     ( -c(i, j)%Axi_half - c(i, j)%Axj_half ) * c(i, j)%dTdy

    p(i+1, j+1)%d2Td2y = p(i+1, j+1)%d2Td2y + &
     (  c(i, j)%Axi_half - c(i, j)%Axj_half ) * c(i, j)%dTdy

    p(i+1,   j)%d2Td2y = p(i+1,   j)%d2Td2y +   &
     (  c(i, j)%Axi_half + c(i, j)%Axj_half ) * c(i, j)%dTdy

    p(i,     j)%d2Td2y = p(i,     j)%d2Td2y +     &
     ( -c(i, j)%Axi_half + c(i, j)%Axj_half ) * c(i, j)%dTdy

  end subroutine
!
!  ! This finds the difference between two points in the x or y direction.
!  real*8 function dx(p1, p2)
!    type (GridPoint), pointer :: p1, p2
!    dx = p2%x - p1%x
!  end function
!
!  real*8 function dy(p1, p2)
!    type (GridPoint), pointer :: p1, p2
!    dy = p2%y - p1%y
!  end function
!
!  subroutine all_in_one(c, p, i, j)
!    save
!    type (GridPoint) :: p(1:IMAX, 1:JMAX)
!    type (GridPoint) :: a, b, c, d
!    type (GridPoint) :: a_, b_, c_, d_
!!    type (GridCell) :: c(1:IMAX-1, 1:JMAX-1)
!    integer :: i, j, d_i, d_j
!    real*8 :: temp, prim_vol, secon_vol
!    type (GridPoint), pointer :: p1, p2, p3, p4
!
!  Volume(p1, p2, p3, p4) = abs( &
!            (p1%x *  p2%y - p1%y * p2%x) + &
!            (p2%x *  p3%y - p2%y * p3%x) + &
!            (p3%x *  p4%y - p3%y * p4%x) + &
!            (p4%x *  p1%y - p4%y * p1%x)   &
!           ) / 2.
!
!  T_a(i, j) = ( p(i, j)%T + p(i - 1, j)%T + p(i - 1, j - 1)%T + p(i, j - 1)%T ) / 4.
!  T_b(i, j) = ( p(i, j)%T + p(i + 1, j)%T + p(i + 1, j - 1)%T + p(i, j - 1)%T ) / 4.
!
!  ! In text we see this as dT/dy(i, j - 1/2), but we call it as dT(dy, i, j, 0, -1).
!  dT(delta, i, j, d_i, d_j) = ( p(i + d_i, j + d_j)%T * delta(a_, b_) + T_b(i, j) * delta(b_, c_) + &
!                                p(i,             j)%T * delta(c_, d_) + T_a(i, j) * delta(d_, a_) ) / Volume(a_, b_, c_, d_)
!
!  ! In text we see this as the sum of the dT/dx s, but here we call it as Sum_dT(dx, dy, i, j).
!  Sum_dT(dir, delta, i, j) = dT(dir, i, j, 0, -1) * delta(a, b) + dT(dir, i, j,  1, 0) * delta(b, c) + &
!                             dT(dir, i, j, 0,  1) * delta(c, d) + dT(dir, i, j, -1, 0) * delta(d, a)
!
!  ! This is the final function to determine the temperature at point i, j at step n + 1.
!  p(i, j)%T = p(i, j)%T +  step_size * ( alpha/Volume(a, b, c, d) * ( Sum_dT(dx, dy, i, j) - Sum_dT(dy, dx, i, j) ))
!
!  a%x = ( p(i, j)%x + p(i, j - 1)%x + p(i - 1, j)%x + p(i - 1, j - 1)%x ) / 4.
!  a%y = ( p(i, j)%y + p(i, j - 1)%y + p(i - 1, j)%y + p(i - 1, j - 1)%y ) / 4.
!  b%x = ( p(i, j)%x + p(i, j - 1)%x + p(i + 1, j)%x + p(i + 1, j - 1)%x ) / 4.
!  b%y = ( p(i, j)%y + p(i, j - 1)%y + p(i + 1, j)%y + p(i + 1, j - 1)%y ) / 4.
!  c%x = ( p(i, j)%x + p(i, j + 1)%x + p(i + 1, j)%x + p(i + 1, j + 1)%x ) / 4.
!  c%y = ( p(i, j)%y + p(i, j + 1)%y + p(i + 1, j)%y + p(i + 1, j + 1)%y ) / 4.
!  d%x = ( p(i, j)%x + p(i, j + 1)%x + p(i - 1, j)%x + p(i - 1, j + 1)%x ) / 4.
!  d%y = ( p(i, j)%y + p(i, j + 1)%y + p(i - 1, j)%y + p(i - 1, j + 1)%y ) / 4.
!
!  a%x = ( p(i, j)%x + p(i, j - 1)%x + p(i - 1, j)%x + p(i - 1, j - 1)%x ) / 4.
!  a%y = ( p(i, j)%y + p(i, j - 1)%y + p(i - 1, j)%y + p(i - 1, j - 1)%y ) / 4.
!  b%x = ( p(i, j)%x + p(i, j - 1)%x + p(i + 1, j)%x + p(i + 1, j - 1)%x ) / 4.
!  b%y = ( p(i, j)%y + p(i, j - 1)%y + p(i + 1, j)%y + p(i + 1, j - 1)%y ) / 4.
!  c%x = ( p(i, j)%x + p(i, j + 1)%x + p(i + 1, j)%x + p(i + 1, j + 1)%x ) / 4.
!  c%y = ( p(i, j)%y + p(i, j + 1)%y + p(i + 1, j)%y + p(i + 1, j + 1)%y ) / 4.
!  d%x = ( p(i, j)%x + p(i, j + 1)%x + p(i - 1, j)%x + p(i - 1, j + 1)%x ) / 4.
!  d%y = ( p(i, j)%y + p(i, j + 1)%y + p(i - 1, j)%y + p(i - 1, j + 1)%y ) / 4.
!
!  end subroutine

end module UpdateTemperature
