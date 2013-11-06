module MainRoutines
  use constants
  use GridPointModule
  use GridCellModule
  use BlockModule
  use UpdateTemperature
  use plot3D_module

  implicit none

contains

  subroutine initialization(Blocks)
    type (BlockType) :: Blocks(:)
    !    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
    !    type (GridCell),  target :: Cells(1:IMAX-1, 1:JMAX-1)
    !    integer :: i, j

    ! We first read in the connectivity file.
    call read_configuration_file(Blocks)

    ! We then read in the grid file.
    call read_grid_file(Blocks)

    ! We then read in the initial temperature file.
    call read_temp_file(Blocks)

    !  Initialize the points.
    call initialize_points(Blocks)

    !  Initialize the cells.
    call initialize_cells(Blocks)

    ! Set up secondary areas needed for integration.
    call set_secondary_areas(Blocks)

    ! Calculate constants for integration.
    call set_constants(Blocks)

    ! Calculate proper bounds for each block.
    call set_bounds(Blocks)
  end subroutine

  subroutine solve(Blocks, step)
    type (BlockType), target :: Blocks(:)
    type (GridPoint), pointer :: Points(:,:)
    !    type (GridCell), pointer :: Cells(:,:)
    real(kind=8) :: temp_residual = 1.d0, residual = 1.d0 ! Arbitrary initial residuals.
    integer :: n_, step, max_steps = 10000

    write(*,*), "Start solver"

    !  Begin main loop, stop if we hit our mark or after max_steps iterations.
    do while (residual >= .00001d0 .and. step < max_steps)
      ! Another day...
      step = step + 1

      ! Calculate our first and second derivatives for all our points.
      ! Calculate the new temperature for all of our interior points.
      residual = 0d0
      call derivatives(Blocks)
!      Blocks%Points%T = Blocks%Points%T + Blocks%Points%tempT
      call update_ghosts(Blocks)

      ! Find block with largest residual.
      do n_ = 1, nBlocks
        temp_residual = maxval(abs(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))
        write(*,*), "n ", n_, temp_residual

        if (temp_residual > residual) then
          residual = temp_residual
        end if
      end do
      write(*,*), step, residual
    end do

    ! Check for convergence.
    if (step < max_steps) then
      write(*,*) "Converged."
    else
      write(*,*) "Failed to converge."
    end if

  end subroutine
end module MainRoutines
