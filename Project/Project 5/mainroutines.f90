module MainRoutines
  use BlockModule
  use plot3D_module
  use GridCreation

  implicit none

  ! Block array for solver.
  type (BlockType), allocatable, target :: Blocks(:)

contains

  ! This subroutine partitions the grid into blocks.
  subroutine initialize_grid(b)
    type (BlockType) :: b(:)

    ! Create blocks and find neighbors.
    call create_blocks(b)

    ! Set up the points in each block.
    call initialize_block_grid(b)

    ! Set initial temperatures.
    call initialize_block_temp(b)

    ! Calculate proper bounds for each block.
    call set_bounds(b)

    write(*,*), 'Processor ', MyID, ' initialized simulation.'
    write(*,*)
  end subroutine

  ! This subroutine reads the block files and initializes the simulation.
  subroutine initialization
    ! We then read in the grid file.
    call read_grid_file(Blocks)

    ! We then read in the initial temperature file.
    call read_temp_file(Blocks)

    !  Initialize the points.
    call initialize_points(Blocks)

    !  Initialize the primary face areas and volumes.
    call initialize_faces_and_volumes(Blocks)

    ! Calculate constants for integration.
    call set_constants(Blocks)
  end subroutine

  ! This is the main solver.
  subroutine solve(Blocks)
    type (BlockType), target :: Blocks(:)
    integer, parameter :: max_steps = 10000
    real(kind=8) :: temp_residual = 1.d0, residual = 1.d0 ! Arbitrary initial residuals.
    real(kind=8) :: residuals(max_steps) = 0.d0
    integer :: n_
    character(2) :: name, str

    !  Begin main loop, stop if we hit our mark or after max_steps iterations.
    do while (residual >= .00001d0 .and. step < max_steps)
      ! Another day...
      step = step + 1

      ! Calculate our first and second derivatives for all our points, and
      ! calculate the new temperature for all of our interior points.
      call derivatives

      call MPI_Barrier(barrier, ierror)

      call update_ghosts

      ! Hold after we send information. This requires that all Procs have
      ! send and received their info for this loop before we continue.
      call MPI_Barrier(barrier, ierror)

      ! Find block with largest residual.
      residual = 0.d0
      do n_ = 1, MyNBlocks
        temp_residual = maxval(abs(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))

        if (temp_residual > residual) then
          residual = temp_residual
        end if
      end do

      call MPI_Barrier(barrier, ierror)
      write(*,*), MyID, step, residual

      residuals(step) = residual
    end do

    ! Write the residual information to output file.
!     write( name, '(i2)' )  MyID
!     read( name, * ) str

!     open(unit = 666, file = 'convergence.dat.p'//str, form='formatted')
!     write(666,*), residuals
!     close(666)

    ! Check for convergence.
    call MPI_Barrier(barrier, ierror)
    if (step < max_steps) then
      write(*,*) "Converged."
    else
      write(*,*) "Failed to converge."
    end if
  end subroutine

  subroutine derivatives
    type (BlockType), pointer :: MyBlock
    type (GridPoint), pointer :: p(:,:)
    real(kind=8) :: dTdx, dTdy
    integer :: i, j, n_

    ! Loop over each block.
    do n_ = 1, MyNBlocks
      MyBlock => Blocks(n_)
      p => MyBlock%Points

      ! Reset the change in temperature to zero before we begin summing again.
      p%tempT = 0.d0

      do j = MyBlock%localJMIN, MyBlock%localJMAX
        do i =  MyBlock%localIMIN, MyBlock%localIMAX

          ! Trapezoidal counter-clockwise integration to get the first
          ! derivatives in the x/y directions at the cell-center using
          ! Gauss's theorem.
          dTdx = + 0.5d0 * &
            ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * p(i+1, j)%Ayi - &
              ( p(i,   j)%T + p(i,  j+1)%T ) * p(i,   j)%Ayi - &
              ( p(i, j+1)%T + p(i+1,j+1)%T ) * p(i, j+1)%Ayj + &
              ( p(i,   j)%T + p(i+1,  j)%T ) * p(i,   j)%Ayj   &
            ) / p(i, j)%V

          dTdy = - 0.5d0 * &
            ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * p(i+1, j)%Axi - &
              ( p(i,   j)%T + p(i,  j+1)%T ) * p(i,   j)%Axi - &
              ( p(i, j+1)%T + p(i+1,j+1)%T ) * p(i, j+1)%Axj + &
              ( p(i,   j)%T + p(i+1,  j)%T ) * p(i,   j)%Axj   &
            ) / p(i ,j)%V

          ! Alternate distributive scheme second-derivative operator. Updates the
          ! second derivative by adding the first times a constant during each time
          ! step. Pass out x and y second derivatives contributions.
          p(i+1,  j)%tempT = p(i+1,  j)%tempT + p(i+1,  j)%const * &
                           ( p(i, j)%yNN * dTdx + p(i, j)%xPP * dTdy )
          p(i,    j)%tempT = p(i,    j)%tempT + p(i,    j)%const * &
                           ( p(i, j)%yPN * dTdx + p(i, j)%xNP * dTdy )
          p(i,  j+1)%tempT = p(i,  j+1)%tempT + p(i,  j+1)%const * &
                           ( p(i, j)%yPP * dTdx + p(i, j)%xNN * dTdy )
          p(i+1,j+1)%tempT = p(i+1,j+1)%tempT + p(i+1,j+1)%const * &
                           ( p(i, j)%yNP * dTdx + p(i, j)%xPN * dTdy )
        end do
      end do

      ! Update temperatures in the block.
      do j = MyBlock%lowJTemp, MyBlock%localJMAX
        do i =  MyBlock%lowITemp, MyBlock%localIMAX
          p(i,j)%T = p(i,j)%T + p(i,j)%tempT
        end do
      end do
    end do
  end subroutine

  ! After each iteration of temperature updates we need to update our ghost nodes.
  subroutine update_ghosts
    type (BlockType), pointer :: b
    integer :: n_, i, j

    ! For each block...
    do n_ = 1, MyNBlocks
      b => Blocks(n_)

      ! Faces
      ! North face ghost nodes
      if (b%northFace%BC == -1) then
        ! Internal boundary

        if (b%proc == b%northFace%neighborProc) then
          ! Our neighbor block is on same proc, we can grab directly.
          do i = 1, iBlockSize
            b%Points(i, jBlockSize+1)%T = Blocks(b%northFace%neighborLocalBlock)%Points(i, 2)%T
          end do
        else
          ! Our neighbor block is on different proc, so it will also need information from
          ! this block. We must send this information now and do receives at the end of this
          ! subroutine.
          ! Nonblocking send.
          ! MPI CALL HERE

          ! Our neighbor block is on different proc, we must receive with MPI.
          ! Blocking receive.
          ! MPI CALL HERE
!           write(*,*), MyID, " North info is on other proc."
        end if
      end if

      ! South face ghost nodes
      if (b%southFace%BC == -1) then
        if (b%proc == b%southFace%neighborProc) then
          ! Internal boundary
          do i = 1, iBlockSize
            b%Points(i, 0)%T = Blocks(b%southFace%neighborLocalBlock)%Points(i, jBlockSize-1)%T
          end do
        else
!           write(*,*), MyID, " South info is on other proc."
        end if
      end if

      ! East face ghost nodes
      if (b%eastFace%BC == -1) then
        if (b%proc == b%eastFace%neighborProc) then
          ! Internal boundary
          do j = 1, jBlockSize
            b%Points(iBlockSize+1, j)%T = Blocks(b%eastFace%neighborLocalBlock)%Points(2, j)%T
          end do
        else
!           write(*,*), MyID, " East info is on other proc."
        end if
      end if

      ! West face ghost nodes
      if (b%westFace%BC == -1) then
        if (b%proc == b%westFace%neighborProc) then
          ! Internal boundary
          do j = 1, jBlockSize
            b%Points(0, j)%T = Blocks(b%westFace%neighborLocalBlock)%Points(iBlockSize-1, j)%T
          end do
        else
!           write(*,*), MyID, " West info is on other proc."
        end if
      end if

      ! Corners
      ! North east corner
      if (b%NECorner%BC == -1) then
        if (b%proc == b%NECorner%neighborProc) then
          b%Points(iBlockSize+1,jBlockSize+1)%T = Blocks(b%NECorner%neighborLocalBlock)%Points(2,2)%T
        else
!           write(*,*), MyID, " NE info is on other proc."
        end if
      end if

      ! South east corner
      if (b%SECorner%BC == -1) then
        if (b%proc == b%SECorner%neighborProc) then
          b%Points(iBlockSize+1,0)%T = Blocks(b%SECorner%neighborLocalBlock)%Points(2,jBlockSize-1)%T
        else
!           write(*,*), MyID, " SE info is on other proc."
        end if
      end if

      ! South west corner
      if (b%SWCorner%BC == -1) then
        if (b%proc == b%SWCorner%neighborProc) then
          b%Points(0,0)%T = Blocks(b%SWCorner%neighborLocalBlock)%Points(iBlockSize-1,jBlockSize-1)%T
        else
!           write(*,*), MyID, " SW info is on other proc."
        end if
      end if

      ! North west corner
      if (b%NWCorner%BC == -1) then
        if (b%proc == b%NWCorner%neighborProc) then
          b%Points(0,jBlockSize+1)%T = Blocks(b%NWCorner%neighborLocalBlock)%Points(iBlockSize-1,2)%T
        else
!           write(*,*), MyID, " NW info is on other proc."
        end if
      end if

    end do

  end subroutine
end module MainRoutines
