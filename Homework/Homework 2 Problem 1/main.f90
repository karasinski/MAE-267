
program linear_least_squares_fit
  implicit none

  integer, parameter :: length = 20
  real, dimension(1:length) :: x, y

  interface
    subroutine calc(x, y)
      real, dimension(:) :: x, y
    end subroutine calc
  end interface

  ! Declare our x values.
  x = (/ -4.91, -3.84, -2.41, -2.62, -3.78, -0.52, -1.83, &
    -2.01, 0.28, 1.08, -0.94, 0.59, 0.69, 3.04, 1.01, 3.60, &
    4.53, 5.13, 4.43, 4.12 /)

  ! Declare our y values.
  y = (/ -8.18, -7.49, -7.11, -6.15, -5.62, -3.30, -2.05, &
    -2.83, -1.16, 0.52, 0.21, 1.73, 3.96, 4.26, 5.75,     &
    6.67, 7.70, 7.31, 9.05, 10.95 /)

  call calc(x, y)

end program

subroutine calc(x, y)
  implicit none

  integer :: length_x, length_y
  real, dimension(:) :: x, y
  real x_sum, y_sum, xy_sum, xx_sum, yy_sum, x_ave, y_ave, m, b, r

  length_x = size(x)
  length_y = size(y)

  if (length_x /= length_y) then
    write (*, *), "The length of your arrays is not the same, this won't work!"
  else
    ! Calculate neccesary values.
    x_sum = sum(x)
    y_sum = sum(y)
    xy_sum = sum(x*y)
    xx_sum = sum(x*x)
    yy_sum = sum(y*y)
    x_ave = sum(x)/length_x
    y_ave = sum(y)/length_y

    ! Calculate fit values.
    m = (xy_sum - x_sum * y_ave) / (xx_sum - x_sum * x_ave)
    b = y_ave - m * x_ave

    ! Calculate the correlation coefficient.
    r = (length_x * xy_sum - x_sum * y_sum) / &
        sqrt((length_x * xx_sum - x_sum ** 2)*(length_x * yy_sum - y_sum ** 2))

    ! Write our calculated values.
    write (*, *), "Slope, m = ", m, " Intercept, b = ", b
    write (*, *), "Correlation coeff, r = ", r
  end if

end subroutine
