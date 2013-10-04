module constants
  real, parameter, public :: k = 18.8, rho = 8000., c_p = 500.
  real, parameter, public :: pi = 3.141592654, rot = 30.*pi/180.
  integer :: IMAX, JMAX !101 or 501 (2 cases)
contains
  subroutine grid_size(length)
    IMAX = length
    JMAX = length
  end subroutine grid_size
end module

module GridPointModule
  use constants

  public
  type GridPoint
    integer :: i, j
    real :: x, xp, y, yp, T
  end type GridPoint

contains
  subroutine initialize(p, i, j)
    type(GridPoint), intent(inout) :: p
    integer :: i, j

    p%i = i
    p%j = j

    p%xp = cos(0.5*pi*(IMAX-real(i))/(IMAX-1))
    p%yp = cos(0.5*pi*(JMAX-real(j))/(JMAX-1))

    p%x = p%xp*cos(rot)+(1.-p%yp)*sin(rot)
    p%y = p%yp*cos(rot)+(p%xp)*sin(rot)

  end subroutine initialize

  subroutine set_temperature(p, iT)
    type(GridPoint), intent(inout) :: p
    real, optional :: iT

    if (present(iT)) then
      p%T = iT
    else
      p%T = 3.5
    endif
  end subroutine set_temperature

  real function kineticenergy(p)
    type(GridPoint), intent(in) :: p
    integer i
    real :: ke = 0.0

    do i=0,2
      ke = ke + 1
    end do
    kineticenergy = p%T * ke
    return
  end function kineticenergy
end module GridPointModule

program heat
  use constants
  use GridPointModule

  integer :: a, b
  type (GridPoint) :: Point
  type (GridPoint), allocatable :: Points(:)

  call grid_size(101)
  allocate(Points(1: IMAX))

  open (unit = 1, file = "data.dat")
  do a = 1, IMAX
    do b = 1, IMAX
      Point = Points((a+b-1))

      call set_temperature(Point, 50.)
      call initialize(Point, a, b)
      write (1,*), Point
      write (1,*)
    end do
  end do
  close(1)



end program heat
