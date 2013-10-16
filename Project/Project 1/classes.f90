! Various constants we need in a lot of places and a routine to set
! the size of the grid.
module constants
  implicit none
  real*8, parameter :: CFL = 0.5
  real*8, parameter :: k = 18.8, rho = 8000., c_p = 500.
  real*8, parameter :: pi = 3.141592654, rot = 30.*pi/180.
  real*8  :: alpha = k / (c_p * rho)
  integer :: IMAX, JMAX !101 or 501 (2 cases)

contains
  subroutine SetGridSize(length)
    integer :: length
    IMAX = length
    JMAX = length
  end subroutine SetGridSize
end module

! Prof's clock module.
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

! GridPointModule contains the GridPoint type and routines
! to initialize and set boundary conditions.
module GridPointModule
  use constants
  implicit none

  public
  type GridPoint
    integer :: i, j
    real*8 :: x, xp, y, yp
    real*8 :: T, tempT, d2Td2x, d2Td2y
    real*8 :: timestep, Vol2
  end type GridPoint

contains
  ! Set the initial locations and initial temperature of
  ! each grid point.
  subroutine initialize_points(p, i, j)
    type (GridPoint), intent(inout) :: p
    integer :: i, j

    p%i = i
    p%j = j

    p%xp = cos( 0.5 * pi * dfloat(IMAX-i) / dfloat(IMAX-1) )
    p%yp = cos( 0.5 * pi * dfloat(JMAX-j) / dfloat(JMAX-1) )

    p%x = p%xp * cos( rot ) + ( 1. - p%yp ) * sin( rot )
    p%y = p%yp * cos( rot ) + ( p%xp ) * sin( rot )

    p%T = 3.5
  end subroutine initialize_points

  ! This is used to set the boundary conditions.
  subroutine set_temperature(p, T)
    type (GridPoint), intent(inout) :: p
    real*8 :: T

    p%T = T
  end subroutine set_temperature
end module GridPointModule

module GridCellModule
  use GridPointModule
  implicit none

  public
  type GridCell
    real*8 :: V
    real*8 :: dTdx, dTdy
    real*8 :: yPP, yNP, yNN, yPN, xNN, xPN, xPP, xNP
    real*8 :: Ayi_half, Axi_half, Ayj_half, Axj_half
  end type GridCell

contains
  subroutine initialize_cells(c, points, i, j)
    type (GridCell), intent(inout) :: c
    type (GridPoint), target :: points(1:IMAX, 1:JMAX)
    integer :: i, j

    ! Calculate the volume of each cell.
    c%V = ( Points(i+1, j)%x - Points(i, j)%x) * &
          ( Points(i, j+1)%y - Points(i, j)%y)
  end subroutine initialize_cells

  subroutine set_secondary_areas(c, p, i, j)
    type (GridCell), intent(inout) :: c
    type (GridPoint), target :: p(1:IMAX, 1:JMAX)
    integer :: i, j
    real*8 :: Ayi, Axi, Ayj, Axj

    ! These 'area's are used to do trapezoidal counter-clockwise
    ! integration to get the first derivatives in the x/y directions
    ! at the cell-center using Gauss's theorem.

    Ayi(i,j) = ( p(i, j+1)%y - p(i, j)%y )
    Axi(i,j) = ( p(i, j+1)%x - p(i, j)%x )
    Ayj(i,j) = ( p(i+1, j)%y - p(i, j)%y )
    Axj(i,j) = ( p(i+1, j)%x - p(i, j)%x )

    ! These areas are used in the calculation of fluxes
    ! for the alternate distributive scheme second-
    ! derivative operator.
    c%Ayi_half = ( Ayi(i+1, j) + Ayi(i, j)   ) / 4.
    c%Axi_half = ( Axi(i+1, j) + Axi(i, j)   ) / 4.
    c%Ayj_half = ( Ayj(i, j)   + Ayj(i, j+1) ) / 4.
    c%Axj_half = ( Axj(i, j)   + Axj(i, j+1) ) / 4.

    ! And these are the numbers that actually appear in the equations,
    ! saved here to (hopefully) save a moment or two during iteration.
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
  use GridPointModule
  use GridCellModule

  implicit none
  public

contains
  subroutine first_derivative(p, c, i, j)
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
      ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * Ayi(i+1, j) - &
        ( p(i,   j)%T + p(i,  j+1)%T ) * Ayi(i,   j) - &
        ( p(i, j+1)%T + p(i+1,j+1)%T ) * Ayj(i, j+1) + &
        ( p(i,   j)%T + p(i+1,  j)%T ) * Ayj(i,   j)   &
      ) / ( 2. * c(i, j)%V )

    c(i, j)%dTdy = - &
      ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * Axi(i+1, j) - &
        ( p(i,   j)%T + p(i,  j+1)%T ) * Axi(i,   j) - &
        ( p(i, j+1)%T + p(i+1,j+1)%T ) * Axj(i, j+1) + &
        ( p(i,   j)%T + p(i+1,  j)%T ) * Axj(i,   j)   &
      ) / ( 2. * c(i ,j)%V )
  end subroutine

  subroutine second_derivative(p, c, i, j)
    type (GridPoint) :: p(1:IMAX, 1:JMAX)
    type (GridCell)  :: c(1:IMAX-1, 1:JMAX-1)
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
