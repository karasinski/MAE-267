module constants
  real, parameter, public :: k = 18.8, rho = 8000., c_p = 500.
  real, parameter, public :: pi = 3.141592654, rot = 30.*pi/180.
  integer :: IMAX, JMAX !101 or 501 (2 cases)
contains
  subroutine calc(length)
    IMAX = length
    JMAX = length
  end subroutine calc
end module

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

subroutine compute()
  use constants

  implicit none
  integer :: i, j
  real :: xp, yp

  call calc(101)
  write (*,*), xp(1), IMAX
end subroutine compute

program heat_conduction
  use clock

  call start_clock()
  call compute()
  call end_clock()
end program

real function xp(i)
  use constants
  implicit none
  integer, intent(in) :: i

  xp = cos(0.5*pi*(IMAX-i)/(IMAX-1))
end function

real function yp(j)
  use constants
  implicit none
  integer, intent(in) :: j

  yp = cos(0.5*pi*(JMAX-j)/(JMAX-1))
end function
