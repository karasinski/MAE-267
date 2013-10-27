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
    call initialize_points(Points)

    !  Initialize Cells.
    call initialize_cells(Cells, Points)

    ! Set up secondary areas needed for integration.
    call set_secondary_areas(Cells, Points)

    ! Calculate constants for integration.
    call set_constants(Cells, Points)

    ! Set up Dirichlet condition.
    ! Consider refactor.
    do j = 1, JMAX
      call set_temperature(Points(1,j), 3.d0 * Points(1,j)%yp + 2.d0)
      call set_temperature(Points(IMAX,j), 3.d0 * Points(IMAX,j)%yp + 2.d0)
    end do

    do i = 1, IMAX
      call set_temperature(Points(i,1), abs(cos(pi * Points(i,1)%xp)) + 1.d0)
      call set_temperature(Points(i,JMAX), 5.d0 * (sin(pi * Points(i,JMAX)%xp) + 1.d0))
    end do
  end subroutine

  subroutine identify_grid(x, n_, m_)
    integer :: x, m_, n_
    n_ = x/1000
    m_ = x - n_ * 1000
    write(*, *), "x ", x, " n ", n_, " m ", m_
  end subroutine identify_grid

  subroutine make_blocks(Points, Blocks)
    integer :: i, j, i_, j_, m_, n_
    integer :: iBound, jBound
    integer :: BlocksFile  = 99   ! Unit for grid files
!    integer :: counter = 0
    type (GridPoint) :: Points(1:IMAX, 1:JMAX)
    type (GridPoint), allocatable :: Blocks(:,:,:,:)

    ! Size of each block.
    iBound = 1 + (IMAX - 1) / N
    jBound = 1 + (JMAX - 1) / M

    ! Allocate up a blocks array to store all our blocks.
    ! To be used when we eventually pass out blocks to indivual processors.
    20     format(10I10)
    open(unit=BlocksFile,file='blocks.dat',form='formatted')

    ! Can identify a block based off its unique ID by calling identify grid
    !  call identify_grid(5004, n_, m_)

    do n_ = 1, N
      do m_ = 1, M
        ! Block unique ID.
        write (BlocksFile, *), n_ * 1000 + m_

        ! Each block has neighbors (and/or BCs).
        ! Here we either write out the neighboring block
        ! or the code for a boundary block, '-1'

        ! Neighbor 1
        if (m_ - 1 > 0) then
          write (BlocksFile, *), n_ * 1000 + (m_ - 1)
        else
          write (BlocksFile, *), -1
        end if
        ! Neighbor 2
        if (m_ + 1 <= M) then
          write (BlocksFile, *), n_ * 1000 + (m_ + 1)
        else
          write (BlocksFile, *), 1
        end if
        ! Neighbor 3
        if (n_ - 1 > 0) then
          write (BlocksFile, *), (n_ - 1) * 1000 + m_
        else
          write (BlocksFile, *), -1
        end if
        ! Neighbor 4
        if (n_ + 1 <= N) then
          write (BlocksFile, *), (n_ + 1) * 1000 + m_
        else
          write (BlocksFile, *), -1
        end if

        ! Each block has an orientation. This doesn't do anything.
        write(BlocksFile, *), 1

        j = 0
        do j_ = 1 + (m_ - 1) * jBound, m_ * jBound
        i = 0
          j = j + 1
          do i_ = 1 + (n_ - 1) * iBound, n_ * iBound
            i = i + 1

            ! If we're passed the number of points continue.
            if (i_ > IMAX .or. j_ > JMAX) then
              continue
            else
              ! Hand out points to the blocks.
              Blocks(m_, n_, i, j) = Points(i_, j_)
              write (BlocksFile, 20), i_, j_, i, j

              ! Counter to test that we have the right number of points
              ! when we are done creating all the blocks.
!              counter = counter + 1
            end if

          end do
        end do
        write (BlocksFile, *)
      end do
    end do

    ! Did we hit the correct number of points?
!    write(*, *), counter, IMAX*JMAX, counter == IMAX*JMAX

    ! Close our file.
    close(BlocksFile)
  end subroutine make_blocks

  subroutine solve(Points, Cells, step)
    type (GridPoint) :: Points(1:IMAX, 1:JMAX)
    type (GridCell)  :: Cells(1:IMAX-1, 1:JMAX-1)
    real(kind=8) :: residual = 1.d0 ! Arbitrary initial residual.
    integer :: step

    !  Begin main loop, stop if we hit our mark or after 1,000,000 iterations.
    do while (residual >= .00001d0 .and. step <= 1000000)
      ! Another day...
      step = step + 1

      ! Calculate our first and second derivatives for all our points.
      call derivatives(Points, Cells)

      ! Calculate the new temperature for all of our interior points.
      Points(2:IMAX-1, 2:JMAX-1)%tempT = Points(2:IMAX-1, 2:JMAX-1)%const * &
                                         ( Points(2:IMAX-1, 2:JMAX-1)%d2Td2x + Points(2:IMAX-1, 2:JMAX-1)%d2Td2y )

      ! Update all our temperatures.
      Points(2:IMAX-1, 2:JMAX-1)%T = Points(2:IMAX-1, 2:JMAX-1)%T + Points(2:IMAX-1, 2:JMAX-1)%tempT
      residual = maxval(abs(Points(2:IMAX-1, 2:JMAX-1)%tempT))
    end do

    if (step <= 1000000) then
      write(*,*) "Converged."
    else
      write(*,*) "Failed to converge."
    end if

  end subroutine
end module MainRoutines
