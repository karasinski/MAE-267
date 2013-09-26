program period_of_a_pendulum

! Calculate the period of a pendulum of length L.
!
! Declare variables and constants.
! constants=pi, g (acceleration due to gravity)
! variables=L (length of pendulum)

  implicit none    ! Require all variables to be explicitly declared

  integer :: ierr
  character(1) :: yn
  real :: L, T
  real, parameter :: pi = 3.141592653589793, g = 9.81

  interactive_loop: do

!   Prompt the user for length and read it.

    write (*,*) 'Enter the length of the pendulum in meters.'
    read (*,*,iostat=ierr) L

!   If length could not be read from input,
!   then cycle through the loop.

    if (ierr /= 0) then
      write(*,*) 'Error, invalid input.'
      cycle interactive_loop
    end if

!   Compute the period.

    T = 2 * pi * sqrt(L / g)

!   Write the input variables (L)
!   and output (T) to the screen.

    write (*,'(1x,a6,f6.2,a7,5x,a6,f6.2,a8)') &
      'Length=',L,' meters','Period=',T,' seconds'

    yn = ' '
    yn_loop: do
      write(*,*) 'Perform another calculation? y[n]'
      read(*,'(a1)') yn
      if (yn=='y' .or. yn=='Y') exit yn_loop
      if (yn=='n' .or. yn=='N' .or. yn==' ') exit interactive_loop
    end do yn_loop

  end do interactive_loop

end program period_of_a_pendulum
