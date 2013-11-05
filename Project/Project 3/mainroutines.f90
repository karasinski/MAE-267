module MainRoutines
  use constants
  use GridPointModule
  use GridCellModule
  use BlockModule
  use UpdateTemperature
  use plot3D_module

  implicit none

contains

  subroutine initialization(BlocksCollection)
    type (BlockType) :: BlocksCollection(:,:)
!    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
!    type (GridCell),  target :: Cells(1:IMAX-1, 1:JMAX-1)
!    integer :: i, j

    ! We first read in the connectivity file.
    call read_configuration_file(BlocksCollection)

    ! We then read in the grid file.
    call read_grid_file(BlocksCollection)

    ! We then read in the initial temperature file.
    call read_temp_file(BlocksCollection)

    !  Initialize the cells.
    call initialize_cells(BlocksCollection)

    ! Set up secondary areas needed for integration.
    call set_secondary_areas(BlocksCollection)

    ! Calculate constants for integration.
    call set_constants(BlocksCollection)

    ! Calculate proper bounds for each block.
    call set_bounds(BlocksCollection)
  end subroutine

  subroutine solve(Blocks, step)
    type (BlockType), target :: Blocks(:,:)
    type (GridPoint), pointer :: Points(:,:)
!    type (GridCell), pointer :: Cells(:,:)
    real(kind=8) :: temp_residual = 1.d0, residual = 1.d0 ! Arbitrary initial residuals.
    integer :: step, max_steps = 10
    integer :: m_, n_

    !  Begin main loop, stop if we hit our mark or after max_steps iterations.
    do while (residual >= .00001d0 .and. step < max_steps)
      ! Another day...
      step = step + 1

!      do m_ = 1, size(Blocks, 1)
!        do n_ = 1, size(Blocks, 2)
!          Points => Blocks(m_,n_)%Points(2:Blocks(m_,n_)%iBound-1,2:Blocks(m_,n_)%jBound-1)
!
!          ! Update all our temperatures.
!          ! Need to NOT update dirichlet points.
!          Points%T = Points%T + Points%tempT
!
!          ! update ghost nodes
!          ! Vertical Passing
!          if (m_ - 1 > 0) then
!            Blocks(m_ - 1, n_)%Points(:,2) = Blocks(m_, n_)%Points(:,Blocks(m_, n_)%jBound-1)
!          end if
!
!          if (m_ + 1 <= size(Blocks, 1)) then
!!            write(*,*),m_,n_, "> 0, mbound"
!            Blocks(m_ + 1, n_)%Points(:,2) = Blocks(m_, n_)%Points(:,Blocks(m_, n_)%jBound-1)
!          end if
!
!          ! Horizontal passing.
!          if (n_ - 1 > 0) then
!            Blocks(m_, n_ - 1)%Points(2,:) = Blocks(m_, n_)%Points(Blocks(m_, n_)%iBound-1,:)
!          end if
!
!          if (n_ + 1 <= size(Blocks, 2)) then
!!            write(*,*),m_,n_, "> 0, nbound"
!!            write(*,*), 'c', Blocks(m_, n_ - 1)%Points(1,:)%i
!!            Blocks(m_, n_ + 1)%Points(0,:)%T = 10.
!            Blocks(m_, n_ + 1)%Points(2,:) = Blocks(m_, n_)%Points(Blocks(m_, n_)%iBound-1,:)
!!            write(*,*), 'd', Blocks(m_, n_ - 1)%Points(1,:)%i
!          end if
!
!        end do
!      end do


      ! Calculate our first and second derivatives for all our points.
      ! Calculate the new temperature for all of our interior points.
      residual = 0d0
      do m_ = 1, M
        do n_ = 1, N

!      write(*,*), '--------------------'
!      write(*,*), "m ", m_, " n ", n_
          ! need to make ghost nodes, otherwise this won't work
          ! or use accumulation operator

          call derivatives(Blocks(m_, n_))

          ! Find block with largest residual.

          ! temp_residual = maxval(abs(Blocks(m_,n_)%Points(2:Blocks(m_,n_)%iBound - 1, &
          !                                                 2:Blocks(m_,n_)%jBound - 1)%tempT))

          temp_residual = maxval(abs(Blocks(m_,n_)%Points(2:Blocks(m_,n_)%iBound-1, 2:Blocks(m_,n_)%jBound-1)%tempT))
!          write(*,*), temp_residual

          if (temp_residual > residual) then
            residual = temp_residual
          end if
        end do
      end do

!      write(*,*), '--------------------'
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
