module MainRoutines
  use constants
  use GridPointModule
  use GridCellModule
  use UpdateTemperature

contains
  subroutine initialization(Points, Cells)
    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
    type (GridCell),  target :: Cells(1:IMAX-1, 1:JMAX-1)
    integer :: i, j

    !  Initialize grid.
    do j = 1, JMAX
      do i = 1, IMAX
        call initialize_points(Points(i,j), i, j)
      end do
    end do

    !  Initialize Cells.
    do j = 1, JMAX-1
      do i = 1, IMAX-1
        call initialize_cells(Cells(i,j), Points, i, j)
      end do
    end do

    ! Set up secondary areas needed for integration.
    do j = 1, JMAX-1
      do i = 1, IMAX-1
        call set_secondary_areas(Cells(i,j), Points, i, j)
      end do
    end do

    ! Set up Dirichlet condition.
    do j = 1, JMAX
      call set_temperature(Points(1,j), 3. * Points(1,j)%yp + 2.)
      call set_temperature(Points(IMAX,j), 3. * Points(IMAX,j)%yp + 2.)
    end do

    do i = 1, IMAX
      call set_temperature(Points(i,1), abs(cos(pi * Points(i,1)%xp)) + 1.)
      call set_temperature(Points(i,JMAX), 5. * (sin(pi * Points(i,JMAX)%xp) + 1.))
    end do

    ! Calculate timesteps and assign secondary volumes.
    do j = 2, JMAX - 1
      do i = 2, IMAX - 1
        ! Calculate the timestep using the CFL method described in class.
        Points(i, j)%timestep = ( CFL / (2. * alpha) ) * Cells(i, j)%V ** 2 / &
                                ( ( Points(i+1, j)%x - Points(i, j)%x )**2 + &
                                  ( Points(i, j+1)%y - Points(i, j)%y )**2 )

        ! Calculate the secondary volumes around each point. As we have rectangular points,
        ! these are simply the sum of the surrounding primary cells divied by four.
        Points(i, j)%Vol2 = ( Cells(i, j)%V + Cells(i - 1, j)%V + Cells(i, j - 1)%V + Cells(i - 1, j - 1)%V ) / 4.
      end do
    end do
  end subroutine

  subroutine solve(Points, Cells, step)
    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
    type (GridCell),  target :: Cells(1:IMAX-1, 1:JMAX-1)
    real*8, pointer :: Temperature(:,:), tempTemperature(:,:)
    real*8 :: residual = 1. ! Arbitrary initial residual.
    integer :: i, j, step

    !  Set some useful pointers.
    Temperature => Points(2:IMAX-1, 2:JMAX-1)%T
    tempTemperature => Points(2:IMAX-1, 2:JMAX-1)%tempT

    !  Begin main loop, stop if we hit our mark or after 1,000,000 iterations.
    do while (residual >= .00001 .and. step <= 1000000)
      ! Another day...
      step = step + 1

      ! Reset second derivatives to zero before we begin summing again.
      Points%d2Td2x = 0.
      Points%d2Td2y = 0.

      ! Calculate our first and second derivatives for all our points.
      do j = 1, JMAX - 1
        do i = 1, IMAX - 1
          call first_derivative(Points, Cells, i, j)
          call second_derivative(Points, Cells, i, j)
        end do
      end do

      ! Calculate the new temperature for all of our interior points.
      do j = 2, JMAX - 1
        do i = 2, IMAX - 1
          Points(i, j)%tempT = ( Points(i,j)%timestep * alpha / Points(i, j)%Vol2) * &
            ( Points(i, j)%d2Td2x + Points(i, j)%d2Td2y )
        end do
      end do

      ! Update all our temperatures.
      Temperature = Temperature + tempTemperature
      residual = maxval(abs(tempTemperature))
    end do
  end subroutine

  subroutine output(Points, step)
    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
    real*8, pointer :: Temperature(:,:), tempTemperature(:,:)
    integer :: step

    Temperature => Points(2:IMAX-1, 2:JMAX-1)%T
    tempTemperature => Points(2:IMAX-1, 2:JMAX-1)%tempT
    ! Let's find the last cell to change temperature and write some output.
    ! Write down the 'steady state' configuration.
    open (unit = 1, file = "steady_state.dat")
    do i = 1, IMAX
      do j = 1, JMAX
        write (1,'(F10.7, 5X, F10.7, 5X, F10.7, I5, F10.7)'), Points(i,j)%x, Points(i,j)%y, Points(i,j)%T
      end do
    end do
    close (1)

    ! Some output to the screen so we know something happened.
    write (*,*), "IMAX/JMAX", IMAX, JMAX
    write (*,*), "steps", step
    write (*,*), "residual", maxval(tempTemperature)
    write (*,*), "ij", maxloc(tempTemperature)

    ! Write down misc. info asked for by Prof.
    open (unit = 2, file = "other_info.dat")
    write (2,*), "For a ", IMAX, " by ", JMAX, "size grid, we ran for: "
    write (2,*), step, "steps"
    write (2,*), wall_time, "seconds"
    write (2,*)
    write (2,*), "Found max residual of ", maxval(tempTemperature)
    write (2,*), "At ij of ", maxloc(tempTemperature)
    close (2)
    ! End output.
  end subroutine
end module MainRoutines
