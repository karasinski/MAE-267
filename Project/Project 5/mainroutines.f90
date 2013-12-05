module MainRoutines
  use BlockModule
  use plot3D_module
  use GridCreation

  implicit none

  ! Block array for solver.
  type (BlockType), pointer :: Blocks(:)

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
    real(kind=8) :: temp_residual = 1.d0, residual = 1.d0, local_residual = 1.d0
!     real(kind=8) :: residuals(max_steps) = 0.d0
    integer :: n_
    character(2) :: name, str

    ! Processor 0 times our iterations until convergence.
    if (MyID == 0) then
      call start_clock
    end if

    !  Begin main loop, stop if we hit our mark or after max_steps iterations.
    do step = 1, max_steps
      ! Calculate our first and second derivatives for all our points, and
      ! calculate the new temperature for all of our interior points.
      call derivatives
      call update_ghosts
      call cross_proc_talk

      ! Find block with largest residual.
      local_residual = 0.d0
      do n_ = 1, MyNBlocks
        temp_residual = maxval(abs(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))

        if (temp_residual > local_residual) then
          local_residual = temp_residual
        end if
      end do

      call mpi_allreduce(local_residual, residual, 1, MPI_REAL8, MPI_MAX, mpi_comm_world, ierror)
!       residuals(step) = residual
      ! write(*,*), MyID, step, residual, local_residual

      if (residual < .00001d0) exit
    end do

    ! We have converged and all procs have stopped, end clock.
    call MPI_Barrier(barrier, ierror)
    if (MyID == 0) then
      call end_clock
    end if

    if (MyID == 0) then
      write(*,*), step, residual

      ! Write the residual information to output file.
!       write( name, '(i2)' )  mpi_nprocs
!       read( name, * ) str

!       open(unit = 666, file = str//'_convergence.dat', form='formatted')
!       write(666,*), residuals
!       close(666)

      ! Check for convergence.
      if (step < max_steps) then
        write(*,*) "Converged."
      else
        write(*,*) "Failed to converge."
      end if
    end if

    ! Write some results to file/screen.
!     call output(Blocks, residuals)
  end subroutine

  subroutine derivatives
    type (BlockType), pointer :: MyBlock
    type (GridPoint), pointer :: p(:,:)
    type (GridPoint), pointer :: p0, p1, p2, p3
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
          p0 => Blocks(n_)%Points(i, j)
          p1 => Blocks(n_)%Points(i+1, j)
          p2 => Blocks(n_)%Points(i, j+1)
          p3 => Blocks(n_)%Points(i+1, j+1)
          
          ! Trapezoidal counter-clockwise integration to get the first
          ! derivatives in the x/y directions at the cell-center using
          ! Gauss's theorem.
          dTdx = + 0.5d0 * &
            ( ( p1%T + p3%T ) * p1%Ayi - &
              ( p0%T + p2%T ) * p0%Ayi - &
              ( p2%T + p3%T ) * p2%Ayj + &
              ( p0%T + p1%T ) * p0%Ayj   &
            ) / p0%V

          dTdy = - 0.5d0 * &
            ( ( p1%T + p3%T ) * p1%Axi - &
              ( p0%T + p2%T ) * p0%Axi - &
              ( p2%T + p3%T ) * p2%Axj + &
              ( p0%T + p1%T ) * p0%Axj   &
            ) / p0%V

          ! Alternate distributive scheme second-derivative operator. Updates the
          ! second derivative by adding the first times a constant during each time
          ! step. Pass out x and y second derivatives contributions.
          p1%tempT = p1%tempT + p1%const * &
                           ( p0%yNN * dTdx + p0%xPP * dTdy )
          p0%tempT = p0%tempT + p0%const * &
                           ( p0%yPN * dTdx + p0%xNP * dTdy )
          p2%tempT = p2%tempT + p2%const * &
                           ( p0%yPP * dTdx + p0%xNN * dTdy )
          p3%tempT = p3%tempT + p3%const * &
                           ( p0%yNP * dTdx + p0%xPN * dTdy )
        end do
      end do

      ! Update temperatures in the block.
      do j = MyBlock%lowJTemp, MyBlock%localJMAX
        do i =  MyBlock%lowITemp, MyBlock%localIMAX
          p0 => Blocks(n_)%Points(i,j)
          p0%T = p0%T + p0%tempT
        end do
      end do
    end do
  end subroutine

  ! After each iteration of temperature updates we need to update our ghost nodes.
  subroutine update_ghosts
    type (BlockType), pointer :: b
    real(kind=8), pointer :: p1, p2
    integer :: status(MPI_STATUS_SIZE)
    integer :: n_, i, j, tag, destination, source, request
    integer :: i_count = iBlockSize, j_count = jBlockSize
    real(kind = 8) :: i_buffer(iBlockSize), j_buffer(jBlockSize)
    real(kind = 8) :: buffer
    
    ! For each block...
    do n_ = 1, MyNBlocks
      b => Blocks(n_)

      ! North face ghost nodes
      if (b%northFace%BC == INTERNAL_BOUNDARY) then
        ! Our neighbor block is on same proc, we can grab directly.
        do i = 1, iBlockSize
          p1 => Blocks(n_)%Points(i, jBlockSize+1)%T
          p2 => Blocks(b%northFace%neighborLocalBlock)%Points(i, 2)%T
          p1 = p2
        end do
      else if (b%northFace%BC == PROC_BOUNDARY) then
        ! Our neighbor block is on different proc, so it will also need information from
        ! this block. We do a nonblocking send now and a blocking receive later.

        ! Pack our values to send into a buffer.
        do i = 1, iBlockSize
          p1 => Blocks(n_)%Points(i, jBlockSize-1)%T
          i_buffer(i) = p1
        end do

        ! Find the destination.
        destination = b%northFace%neighborProc

        ! Generate a tag unique within the iteration.
        tag = nB

        ! Send everything to the proc and continue without immediate confirmation.
        call MPI_Isend(i_buffer, i_count, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
      end if

      ! South face ghost nodes
      if (b%southFace%BC == INTERNAL_BOUNDARY) then
        do i = 1, iBlockSize
          p1 => Blocks(n_)%Points(i, 0)%T
          p2 => Blocks(b%southFace%neighborLocalBlock)%Points(i, jBlockSize-1)%T
          p1 = p2
        end do
      else if (b%southFace%BC == PROC_BOUNDARY) then
        do i = 1, iBlockSize
          p1 => Blocks(n_)%Points(i, 2)%T
          i_buffer(i) = p1
        end do
        destination = b%southFace%neighborProc
        tag = sB
        call MPI_Isend(i_buffer, i_count, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
      end if

      ! East face ghost nodes
      if (b%eastFace%BC == INTERNAL_BOUNDARY) then
        do j = 1, jBlockSize
          p1 => Blocks(n_)%Points(iBlockSize+1, j)%T
          p2 => Blocks(b%eastFace%neighborLocalBlock)%Points(2, j)%T
          p1 = p2
        end do
      else if (b%eastFace%BC == PROC_BOUNDARY) then
        do j = 1, jBlockSize
          p1 => Blocks(n_)%Points(iBlockSize-1, j)%T
          j_buffer(j) = p1
        end do
        destination = b%eastFace%neighborProc
        tag = eB
        call MPI_Isend(j_buffer, j_count, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
      end if

      ! West face ghost nodes
      if (b%westFace%BC == INTERNAL_BOUNDARY) then
        do j = 1, jBlockSize
          p1 => Blocks(n_)%Points(0, j)%T
          p2 => Blocks(b%westFace%neighborLocalBlock)%Points(iBlockSize-1, j)%T
          p1 = p2
        end do
      else if (b%westFace%BC == PROC_BOUNDARY) then
        do j = 1, jBlockSize
          p1 => Blocks(n_)%Points(2, j)%T
          j_buffer(j) = p1
        end do
        destination = b%westFace%neighborProc
        tag = wB
        call MPI_Isend(j_buffer, j_count, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
      end if

      ! North east corner ghost node
      if (b%NECorner%BC == INTERNAL_BOUNDARY) then
        p1 => Blocks(n_)%Points(iBlockSize+1,jBlockSize+1)%T
        p2 => Blocks(b%NECorner%neighborLocalBlock)%Points(2,2)%T
        p1 = p2
      else if (b%NECorner%BC == PROC_BOUNDARY) then
        destination = b%NECorner%neighborProc
        p1 => Blocks(n_)%Points(iBlockSize-1,jBlockSize-1)%T
        buffer = p1
        tag = nB + eB
        call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
      end if

      ! South east corner ghost node
      if (b%SECorner%BC == INTERNAL_BOUNDARY) then
        p1 => Blocks(n_)%Points(iBlockSize+1,0)%T
        p2 => Blocks(b%SECorner%neighborLocalBlock)%Points(2,jBlockSize-1)%T
        p1 = p2
      else if (b%SECorner%BC == PROC_BOUNDARY) then
        destination = b%SECorner%neighborProc
        p1 => Blocks(n_)%Points(iBlockSize-1,2)%T
        buffer = p1
        tag = sB + eB
        call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
      end if

      ! South west corner ghost node
      if (b%SWCorner%BC == INTERNAL_BOUNDARY) then
        p1 => Blocks(n_)%Points(0,0)%T
        p2 => Blocks(b%SWCorner%neighborLocalBlock)%Points(iBlockSize-1,jBlockSize-1)%T
        p1 = p2
      else if (b%SWCorner%BC == PROC_BOUNDARY) then
        destination = b%SWCorner%neighborProc
        p1 => Blocks(n_)%Points(2,2)%T
        buffer = p1
        tag = sB + wB
        call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
      end if

      ! North west corner ghost node
      if (b%NWCorner%BC == INTERNAL_BOUNDARY) then
        p1 => Blocks(n_)%Points(0,jBlockSize+1)%T
        p2 => Blocks(b%NWCorner%neighborLocalBlock)%Points(iBlockSize-1,2)%T
        p1 = p2
      else if (b%NWCorner%BC == PROC_BOUNDARY) then
        destination = b%NWCorner%neighborProc
        p1 => Blocks(n_)%Points(2,jBlockSize-1)%T
        buffer = p1
        tag = nB + wB
        call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
      end if

    end do
  end subroutine

  subroutine cross_proc_talk
    type (BlockType), pointer :: b
    real(kind=8), pointer :: p1
    integer :: status(MPI_STATUS_SIZE)
    integer :: n_, i, j, tag, destination, source, request
    integer :: i_count = iBlockSize, j_count = jBlockSize
    real(kind = 8) :: i_buffer(iBlockSize), j_buffer(jBlockSize)
    real(kind = 8) :: buffer
    ! For each block...
    do n_ = 1, MyNBlocks
      b => Blocks(n_)
      ! Our neighbor block is on different proc, we must receive with MPI.

      ! South Face Blocking Receive
      if (b%southFace%BC == PROC_BOUNDARY) then
        source = b%southFace%neighborProc
        tag = nB
        call MPI_RECV(i_buffer, i_count, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
        do i = 1, iBlockSize
          p1 => Blocks(n_)%Points(i, 0)%T
          p1 = i_buffer(i)
        end do
      end if

      ! North Face Blocking Receive
      if (b%northFace%BC == PROC_BOUNDARY) then
        source = b%northFace%neighborProc
        tag = sB
        call MPI_RECV(i_buffer, i_count, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
        do i = 1, iBlockSize
          p1 => Blocks(n_)%Points(i, jBlockSize+1)%T
          p1 = i_buffer(i)
        end do
      end if

      ! West Face Blocking Receive
      if (b%westFace%BC == PROC_BOUNDARY) then
        source = b%westFace%neighborProc
        tag = eB
        call MPI_RECV(j_buffer, j_count, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
        do j = 1, jBlockSize
          p1 => Blocks(n_)%Points(0, j)%T
          p1 = j_buffer(j)
        end do
      end if

      ! East Face Blocking Receive
      if (b%eastFace%BC == PROC_BOUNDARY) then
        source = b%eastFace%neighborProc
        tag = wB
        call MPI_RECV(j_buffer, j_count, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
        do j = 1, jBlockSize
          p1 => Blocks(n_)%Points(iBlockSize+1, j)%T
          p1 = j_buffer(j)
        end do
      end if

      ! South West Corner Blocking Receive
      if (b%SWCorner%BC == PROC_BOUNDARY) then
        source = b%SWCorner%neighborProc
        tag = nB + eB
        call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
        p1 => Blocks(n_)%Points(0,0)%T
        p1 = buffer
      end if

      ! North West Corner Blocking Receive
      if (b%NWCorner%BC == PROC_BOUNDARY) then
        source = b%NWCorner%neighborProc
        tag = sB + eB
        call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
        p1 => Blocks(n_)%Points(0,jBlockSize+1)%T
        p1 = buffer
      end if

      ! North East Corner Blocking Receive
      if (b%NECorner%BC == PROC_BOUNDARY) then
        source = b%NECorner%neighborProc
        tag = sB + wB
        call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
        p1 => Blocks(n_)%Points(iBlockSize+1,jBlockSize+1)%T
        p1 = buffer
      end if

      ! South East Corner Blocking Receive
      if (b%SECorner%BC == PROC_BOUNDARY) then
        source = b%SECorner%neighborProc
        tag = nB + wB
        call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
        p1 => Blocks(n_)%Points(iBlockSize+1,0)%T
        p1 = buffer
      end if

    end do
  end subroutine
end module MainRoutines
