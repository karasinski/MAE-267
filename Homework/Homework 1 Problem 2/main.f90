program minimum_tension

  ! Calculate the position to attach a cable to
  ! produce minimum tension.

  implicit none
  real, parameter :: W = 200, lc = 8, lp = 8
  doubleprecision :: min_T, T, opt_d, d = 1.0


  do while (d < 7.0)
    ! Calculate the tension on the cable.
    T = W * lc * lp / (d * sqrt(lp ** 2 - d ** 2))

    ! If this tension is less than the current minimum
    ! tension, note the tension and distance.
    if (T < min_T) then
      min_T = T
      opt_d = d
    end if

    d = d + 0.1
  end do

  write (*,'(1x,a12,f6.2,a5,5x,a12,f8.2,a8)') &
      'Optimal d = ',opt_d,' feet','Minimum T = ',min_T,' Newtons'

end program minimum_tension
