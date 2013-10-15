module constants
  implicit none
  real*8, parameter, public :: k = 18.8, rho = 8000., c_p = 500.
  real*8, parameter, public :: pi = 3.141592654, rot = 30.*pi/180.
  real*8, public :: alpha = k / (c_p * rho)
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
    type (GridPoint), target :: points(1:IMAX, 1:JMAX)
    integer :: i, j

    ! This is my 'fake' method to transfer heat: sum the
    ! nearby points and divide by five. This is a temporary
    ! function to help test display data until I get the
    ! real derivatives working.

    points(i, j)%tempT = points(i, j-1)%T + points(i+1, j)%T + points(i, j+1)%T + points(i-1, j)%T
    !    points(i, j)%tempT = points(i, j)%tempT / 4.

  end subroutine update_temperature

end module GridPointModule

module GridCellModule
  use constants
  use GridPointModule

  implicit none
  save

  public
  type GridCell
    real*8 :: V
    real*8 :: dTdx, dTdy
    real*8 :: yPP, yNP, yNN, yPN, xNN, xPN, xPP, xNP
    real*8 :: Ayi_half, Axi_half, Ayj_half, Axj_half
  end type GridCell

contains
  subroutine initialize_cells(c, points, i, j)
    save
    type (GridCell), intent(inout) :: c
    type (GridPoint), target :: points(1:IMAX, 1:JMAX)
    integer :: i, j

    ! Calculate the volume of each cell.
    c%V = (Points(i+1, j)%x-Points(i, j)%x)*(Points(i, j+1)%y-Points(i, j)%y)

  end subroutine initialize_cells

  subroutine set_secondary_areas(c, p, i, j)
    save
    type (GridCell), intent(inout) :: c
    type (GridPoint), target :: p(1:IMAX, 1:JMAX)
    integer :: i, j
    real*8 :: Ayi, Axi, Ayj, Axj

    ! These 'area's are used to do trapezoidal counter-clockwise
    ! integration to get the first derivatives in the x/y directions
    ! at the cell-center using Gauss's theorem.

    ! These areas are used in the calculation of fluxes
    ! for the alternate distributive scheme second-
    ! derivative operator.
    Ayi(i,j) = ( p(i, j+1)%y - p(i, j)%y )
    Axi(i,j) = ( p(i, j+1)%x - p(i, j)%x )
    Ayj(i,j) = ( p(i+1, j)%y - p(i, j)%y )
    Axj(i,j) = ( p(i+1, j)%x - p(i, j)%x )

    c%Ayi_half = ( Ayi(i+1, j) + Ayi(i, j)   ) / 4.
    c%Axi_half = ( Axi(i+1, j) + Axi(i, j)   ) / 4.
    c%Ayj_half = ( Ayj(i, j)   + Ayj(i, j+1) ) / 4.
    c%Axj_half = ( Axj(i, j)   + Axj(i, j+1) ) / 4.

    c%yPP = (  c%Ayi_half + c%Ayj_half )
    c%yNP = ( -c%Ayi_half + c%Ayj_half )
    c%yNN = ( -c%Ayi_half - c%Ayj_half )
    c%yPN = (  c%Ayi_half - c%Ayj_half )
    c%xNN = ( -c%Axi_half - c%Axj_half )
    c%xPN = (  c%Axi_half - c%Axj_half )
    c%xPP = (  c%Axi_half + c%Axj_half )
    c%xNP = ( -c%Axi_half + c%Axj_half )

  end subroutine
end module GridCellModule

module UpdateTemperature
  use constants
  use GridPointModule
  use GridCellModule

  implicit none
  public

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

    Ayi(i,j) = ( p(i, j+1)%y - p(i, j)%y )
    Axi(i,j) = ( p(i, j+1)%x - p(i, j)%x )
    Ayj(i,j) = ( p(i+1, j)%y - p(i, j)%y )
    Axj(i,j) = ( p(i+1, j)%x - p(i, j)%x )

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
    ! Updates the second derivative by adding the first times a constant
    ! during each time step.

    ! Pass out x second derivatives.
    p(i,  j+1)%d2Td2x = p(i,  j+1)%d2Td2x + c(i, j)%yPP * c(i, j)%dTdx
    p(i+1,j+1)%d2Td2x = p(i+1,j+1)%d2Td2x + c(i, j)%yNP * c(i, j)%dTdx
    p(i+1,  j)%d2Td2x = p(i+1,  j)%d2Td2x + c(i, j)%yNN * c(i, j)%dTdx
    p(i,    j)%d2Td2x = p(i,    j)%d2Td2x + c(i, j)%yPN * c(i, j)%dTdx

    ! Pass out y second derivatives.
    p(i,  j+1)%d2Td2y = p(i,  j+1)%d2Td2y + c(i, j)%xNN * c(i, j)%dTdy
    p(i+1,j+1)%d2Td2y = p(i+1,j+1)%d2Td2y + c(i, j)%xPN * c(i, j)%dTdy
    p(i+1,  j)%d2Td2y = p(i+1,  j)%d2Td2y + c(i, j)%xPP * c(i, j)%dTdy
    p(i,    j)%d2Td2y = p(i,    j)%d2Td2y + c(i, j)%xNP * c(i, j)%dTdy

  end subroutine
end module UpdateTemperature
