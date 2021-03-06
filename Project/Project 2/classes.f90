! Various constants we need in a lot of places and a routine to set
! the size of the grid.
module constants
  implicit none
  real(kind=8), parameter :: CFL = 1.14d0
  real(kind=8), parameter :: k = 18.8d0, rho = 8000.d0, c_p = 500.d0
  real(kind=8), parameter :: pi = 3.141592654d0, rot = 30.d0*pi/180.d0
  real(kind=8), parameter :: alpha = k / (c_p * rho)
  integer :: IMAX, JMAX
  integer :: N, M ! Number of blocks.

contains
  subroutine SetGridSize(length)
    integer :: length
    IMAX = length
    JMAX = length
  end subroutine SetGridSize

  subroutine SetNumberOfBlocks(n_, m_)
    integer :: n_, m_
    N = n_
    M = m_
  end subroutine SetNumberOfBlocks
end module

! Prof's clock module.
module clock
  integer clock_start,clock_end,clock_max,clock_rate
  real(kind=8) wall_time

contains
  subroutine start_clock()
    ! call system time to determine flow solver wall time
    call system_clock(count_max=clock_max,count_rate=clock_rate)
    call system_clock(clock_start)
  end subroutine start_clock

  subroutine end_clock()
    ! determine total wall time for solver
    call system_clock(clock_end)
    wall_time=dfloat(clock_end-clock_start)/dfloat(clock_rate)
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
    real(kind=8) :: x, xp, y, yp
    real(kind=8) :: T, tempT
    real(kind=8) :: timestep, Vol2, const
  end type GridPoint

contains
  ! Set the initial locations and initial temperature of
  ! each grid point.
  subroutine initialize_points(Points)
    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
    type (GridPoint), pointer :: p
    integer :: i, j

    do j = 1, JMAX
      do i = 1, IMAX
        p => Points(i, j)

        p%i = i
        p%j = j

        p%xp = cos( 0.5d0 * pi * dfloat(IMAX-i) / dfloat(IMAX-1) )
        p%yp = cos( 0.5d0 * pi * dfloat(JMAX-j) / dfloat(JMAX-1) )

        p%x = p%xp * cos( rot ) + ( 1.d0 - p%yp ) * sin( rot )
        p%y = p%yp * cos( rot ) + ( p%xp ) * sin( rot )

        p%T = 3.5d0
      end do
    end do
  end subroutine initialize_points

  ! This is used to set the boundary conditions.
  subroutine set_temperature(p, T)
    type (GridPoint), intent(inout) :: p
    real(kind=8) :: T

    p%T = T
  end subroutine set_temperature
end module GridPointModule

module GridCellModule
  use GridPointModule
  implicit none

  public
  type GridCell
    real(kind=8) :: V
    real(kind=8) :: yPP, yNP, yNN, yPN
    real(kind=8) :: xNN, xPN, xPP, xNP
  end type GridCell

contains
  subroutine initialize_cells(Cells, Points)
    type (GridCell), target :: Cells(1:IMAX-1, 1:JMAX-1)
    type (GridPoint) :: Points(1:IMAX, 1:JMAX)
    integer :: i, j

    do j = 1, JMAX-1
      do i = 1, IMAX-1
        ! Calculate the volume of each cell.
        Cells(i, j)%V = ( Points(i+1, j)%xp - Points(i, j)%xp) * &
                        ( Points(i, j+1)%yp - Points(i, j)%yp)
      end do
    end do
  end subroutine initialize_cells

  subroutine set_secondary_areas(Cells, p)
    type (GridCell), pointer :: c
    type (GridCell), target :: Cells(1:IMAX-1, 1:JMAX-1)
    type (GridPoint), target :: p(1:IMAX, 1:JMAX)
    integer :: i, j
    real(kind=8) :: Ayi, Axi, Ayj, Axj
    real(kind=8) :: Ayi_half, Axi_half, Ayj_half, Axj_half

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
    Ayi_half(i,j) = ( Ayi(i+1, j) + Ayi(i, j)   ) * 0.25d0
    Axi_half(i,j) = ( Axi(i+1, j) + Axi(i, j)   ) * 0.25d0
    Ayj_half(i,j) = ( Ayj(i, j)   + Ayj(i, j+1) ) * 0.25d0
    Axj_half(i,j) = ( Axj(i, j)   + Axj(i, j+1) ) * 0.25d0

    do j = 1, JMAX-1
      do i = 1, IMAX-1
        c => Cells(i, j)
        ! And these are the numbers that actually appear in the equations,
        ! saved here to (hopefully) save a moment or two during iteration.
        c%yPP = (  Ayi_half(i,j) + Ayj_half(i,j) )
        c%yNP = ( -Ayi_half(i,j) + Ayj_half(i,j) )
        c%yNN = ( -Ayi_half(i,j) - Ayj_half(i,j) )
        c%yPN = (  Ayi_half(i,j) - Ayj_half(i,j) )

        c%xNN = ( -Axi_half(i,j) - Axj_half(i,j) )
        c%xPN = (  Axi_half(i,j) - Axj_half(i,j) )
        c%xPP = (  Axi_half(i,j) + Axj_half(i,j) )
        c%xNP = ( -Axi_half(i,j) + Axj_half(i,j) )
      end do
    end do
  end subroutine

  subroutine set_constants(Cells, Points)
    type (GridCell), target :: Cells(1:IMAX-1, 1:JMAX-1)
    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
    integer :: i, j

    ! Calculate timesteps and assign secondary volumes.
    do j = 2, JMAX - 1
      do i = 2, IMAX - 1
        ! Calculate the timestep using the CFL method described in class.
        Points(i, j)%timestep = ( ( CFL * 0.5d0 ) / alpha ) * Cells(i, j)%V ** 2 / &
                                ( ( Points(i+1, j)%xp - Points(i, j)%xp )**2 + &
                                  ( Points(i, j+1)%yp - Points(i, j)%yp )**2 )

        ! Calculate the secondary volumes around each point. As we have rectangular points,
        ! these are simply the sum of the surrounding primary cells divied by four.
        Points(i, j)%Vol2 = ( Cells(i, j)%V + Cells(i - 1, j)%V + &
                              Cells(i, j - 1)%V + Cells(i - 1, j - 1)%V ) * 0.25d0

        ! Calculate this constant now so we don't recalculate in the solver loop.
        Points(i, j)%const = ( Points(i, j)%timestep * alpha / Points(i, j)%Vol2 )
      end do
    end do
  end subroutine
end module GridCellModule

module UpdateTemperature
  use GridPointModule
  use GridCellModule

  implicit none
  public

contains
  subroutine derivatives(p, c)
    type (GridPoint), intent(inout) :: p(1:IMAX, 1:JMAX)
    type (GridCell), intent(inout)  :: c(1:IMAX-1, 1:JMAX-1)
    real(kind=8) :: Ayi, Axi, Ayj, Axj
    real(kind=8) :: dTdx, dTdy
    integer :: i, j

    ! Trapezoidal counter-clockwise integration to get the first
    ! derivatives in the x/y directions at the cell-center using
    ! Gauss's theorem.
    Ayi(i,j) = ( p(i, j+1)%y - p(i, j)%y )
    Axi(i,j) = ( p(i, j+1)%x - p(i, j)%x )
    Ayj(i,j) = ( p(i+1, j)%y - p(i, j)%y )
    Axj(i,j) = ( p(i+1, j)%x - p(i, j)%x )

    ! Reset the change in temperature to zero before we begin summing again.
    p%tempT = 0.d0

    do j = 1, JMAX - 1
      do i = 1, IMAX - 1
        dTdx = + 0.5d0 * &
          ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * Ayi(i+1, j) - &
            ( p(i,   j)%T + p(i,  j+1)%T ) * Ayi(i,   j) - &
            ( p(i, j+1)%T + p(i+1,j+1)%T ) * Ayj(i, j+1) + &
            ( p(i,   j)%T + p(i+1,  j)%T ) * Ayj(i,   j)   &
          ) / c(i, j)%V

        dTdy = - 0.5d0 * &
          ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * Axi(i+1, j) - &
            ( p(i,   j)%T + p(i,  j+1)%T ) * Axi(i,   j) - &
            ( p(i, j+1)%T + p(i+1,j+1)%T ) * Axj(i, j+1) + &
            ( p(i,   j)%T + p(i+1,  j)%T ) * Axj(i,   j)   &
          ) / c(i ,j)%V

        ! Alternate distributive scheme second-derivative operator.
        ! Updates the second derivative by adding the first times a constant
        ! during each time step.
        ! Pass out x and y second derivatives contributions.
        p(i+1,  j)%tempT = p(i+1,  j)%tempT + p(i+1,  j)%const * ( c(i, j)%yNN * dTdx + c(i, j)%xPP * dTdy )
        p(i,    j)%tempT = p(i,    j)%tempT + p(i,    j)%const * ( c(i, j)%yPN * dTdx + c(i, j)%xNP * dTdy )
        p(i,  j+1)%tempT = p(i,  j+1)%tempT + p(i,  j+1)%const * ( c(i, j)%yPP * dTdx + c(i, j)%xNN * dTdy )
        p(i+1,j+1)%tempT = p(i+1,j+1)%tempT + p(i+1,j+1)%const * ( c(i, j)%yNP * dTdx + c(i, j)%xPN * dTdy )
      end do
    end do
  end subroutine
end module UpdateTemperature
