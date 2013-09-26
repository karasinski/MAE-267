program minimum_tension

  ! Calculate the position to attach a cable to
  ! produce minimum tension.

  implicit none
  real, parameter :: g = 9.81
  doubleprecision :: min_T, T, W = 200 * g, lc = 8, lp = 8, opt_d, d = 1.0


  do d = 1, 7, 0.1
    ! Calculate the tension on the cable.
    T = W * lc * lp / (d * sqrt(lp ** 2 - d ** 2))

    ! If this tension is less than the current minimum
    ! tension, note the tension and distance.
    if (T < min_T) then
      min_T = T
      opt_d = d
    end if

  end do

  write (*,'(1x,a12,f6.2,a5,5x,a12,f8.2,a8)') &
      'Optimal d = ',opt_d,' feet','Minimum T = ',min_T,' Newtons'

end program minimum_tension
