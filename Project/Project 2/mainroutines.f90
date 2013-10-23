module MainRoutines
  use constants
  use GridPointModule
  use GridCellModule
  use UpdateTemperature

  implicit none

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

  subroutine identify_grid(x, n_, m_)
    integer :: x, m_, n_
    n_ = x/1000
    m_ = x - n_ * 1000
    write(*, *), "x ", x, " n ", n_, " m ", m_
  end subroutine identify_grid

  subroutine make_blocks(Points)
    integer :: i, j, i_, j_, m_, n_
    integer :: iBound, jBound
    integer :: BlocksFile  = 99   ! Unit for grid files
    integer :: counter = 0
    type (GridPoint) :: Points(1:IMAX, 1:JMAX)
    type (GridPoint), allocatable :: Blocks(:,:,:,:)

    ! Size of each block.
    iBound = 1 + (IMAX - 1) / N
    jBound = 1 + (JMAX - 1) / M

    ! Allocate up a blocks array to store all our blocks.
    ! To be used when we eventually pass out blocks to indivual processors.
    allocate(Blocks(1:M, 1:N, 1:iBound, 1: jBound))
    20     format(10I10)
    open(unit=BlocksFile,file='grids.dat',form='formatted')

    ! Can identify a block based off its unique ID by calling identify grid
    !  call identify_grid(5004, n_, m_)

    do n_ = 1, N
      do m_ = 1, M
        ! Block unique ID.
        write (BlocksFile, *), "Block: ", n_ * 1000 + m_

        ! Each block has neighbors (and BCs).
        if (m_ - 1 > 0) then
          write (BlocksFile, *), "Neighbor 1: ", n_ * 1000 + (m_ - 1)
        else
          write (BlocksFile, *), "Neighbor 1: BOUNDARY"
        end if
        if (m_ + 1 <= M) then
          write (BlocksFile, *), "Neighbor 2: ", n_ * 1000 + (m_ + 1)
        else
          write (BlocksFile, *), "Neighbor 2: BOUNDARY"
        end if
        if (n_ - 1 > 0) then
          write (BlocksFile, *), "Neighbor 3: ", (n_ - 1) * 1000 + m_
        else
          write (BlocksFile, *), "Neighbor 3: BOUNDARY"
        end if
        if (n_ + 1 <= N) then
          write (BlocksFile, *), "Neighbor 4: ", (n_ + 1) * 1000 + m_
        else
          write (BlocksFile, *), "Neighbor 4: BOUNDARY"
        end if

        ! Each block has an orientation. This doesn't do anything.
        write(BlocksFile, *), "Orientation: 1"

        j = 0
        do j_ = 1 + (m_ - 1) * jBound, m_ * jBound
        i = 0
          j = j + 1
          do i_ = 1 + (n_ - 1) * iBound, n_ * iBound
            i = i + 1

            ! If we're pass the number of blocks continue.
            if (i_ > IMAX .or. j_ > JMAX) then
              continue
            else
              ! Hand out points to the blocks.
              Blocks(m_, n_, i, j) = Points(i_, j_)
              write (BlocksFile, 20), i_, j_, i, j

              ! Counter to test that we have the right number of points
              ! when we are done creating all the blocks.
              counter = counter + 1
            end if

          end do
        end do
        write (BlocksFile, *)
      end do
    end do

    ! Did we hit the correct number of points?
    write(*, *), counter, IMAX*JMAX, counter == IMAX*JMAX
    close(BlocksFile)
  end subroutine make_blocks

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
    write(*,*) "Converged."
  end subroutine
end module MainRoutines
