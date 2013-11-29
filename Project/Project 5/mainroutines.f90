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
    real(kind=8) :: temp_residual = 1.d0, residual = 1.d0, local_residual = 1.d0
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

      call cross_proc_talk

      ! Hold after we send information. This requires that all Procs have
      ! send and received their info for this loop before we continue.
      call MPI_Barrier(barrier, ierror)

      ! Find block with largest residual.
      local_residual = 0.d0
      do n_ = 1, MyNBlocks
        temp_residual = maxval(abs(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))

        if (temp_residual > local_residual) then
          local_residual = temp_residual
        end if
      end do

      call MPI_Barrier(barrier, ierror)
      call mpi_allreduce(local_residual, residual, 1, MPI_REAL8, MPI_MAX, mpi_comm_world, ierror)

      write(*,*), MyID, step, residual, local_residual

      residuals(step) = residual
    end do

    call MPI_Barrier(barrier, ierror)

    if (MyID == 0) then
      ! Write the residual information to output file.
      write( name, '(i2)' )  mpi_nprocs
      read( name, * ) str

      open(unit = 666, file = str//'_convergence.dat', form='formatted')
      write(666,*), residuals
      close(666)

      ! Check for convergence.
      if (step < max_steps) then
        write(*,*) "Converged."
      else
        write(*,*) "Failed to converge."
      end if
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
    integer :: n_, i, j, tag, destination, source, request
    real(kind = 8) :: buffer
    integer :: status(MPI_STATUS_SIZE)

    ! For each block...
    do n_ = 1, MyNBlocks
      b => Blocks(n_)

      ! Faces
      ! North face ghost nodes
      if (b%northFace%BC == -1) then
        ! Internal boundary

        if (MyID == b%northFace%neighborProc) then
          ! Our neighbor block is on same proc, we can grab directly.
          do i = 1, iBlockSize
            b%Points(i, jBlockSize+1)%T = Blocks(b%northFace%neighborLocalBlock)%Points(i, 2)%T
          end do
        else
          ! Our neighbor block is on different proc, so it will also need information from
          ! this block. We must send this information now and do receives later.
          ! Nonblocking send.

          ! Need to send our North Face to our North Neighbor as its South Face.
          destination = b%northFace%neighborProc
          do i = 1, iBlockSize
            buffer = Blocks(b%southFace%neighborLocalBlock)%Points(i, jBlockSize-1)%T
            tag = nB * 100 + i + b%northFace%neighborBlock * 1000
            call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
          end do
        end if
      end if

      ! South face ghost nodes
      if (b%southFace%BC == -1) then
        if (MyID == b%southFace%neighborProc) then
          ! Internal boundary
          do i = 1, iBlockSize
            b%Points(i, 0)%T = Blocks(b%southFace%neighborLocalBlock)%Points(i, jBlockSize-1)%T
          end do
        else
          destination = b%southFace%neighborProc
          do i = 1, iBlockSize
            buffer = Blocks(b%northFace%neighborLocalBlock)%Points(i, 2)%T
            tag = sB * 100 + i + b%southFace%neighborBlock * 1000
            call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
          end do
        end if
      end if

      ! East face ghost nodes
      if (b%eastFace%BC == -1) then
        if (MyID == b%eastFace%neighborProc) then
          ! Internal boundary
          do j = 1, jBlockSize
            b%Points(iBlockSize+1, j)%T = Blocks(b%eastFace%neighborLocalBlock)%Points(2, j)%T
          end do
        else
          destination = b%eastFace%neighborProc
          do j = 1, jBlockSize
            buffer = Blocks(b%westFace%neighborLocalBlock)%Points(iBlockSize-1, j)%T
            tag = eB * 100 + j + b%eastFace%neighborBlock * 1000
            call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
          end do
        end if
      end if

      ! West face ghost nodes
      if (b%westFace%BC == -1) then
        if (MyID == b%westFace%neighborProc) then
          ! Internal boundary
          do j = 1, jBlockSize
            b%Points(0, j)%T = Blocks(b%westFace%neighborLocalBlock)%Points(iBlockSize-1, j)%T
          end do
        else
          destination = b%westFace%neighborProc
          do j = 1, jBlockSize
            buffer = Blocks(b%eastFace%neighborLocalBlock)%Points(2, j)%T
            tag = wB * 100 + j + b%westFace%neighborBlock * 1000
            call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
          end do
        end if
      end if

      ! Corners
      ! North east corner
      if (b%NECorner%BC == -1) then
        if (MyID == b%NECorner%neighborProc) then
          b%Points(iBlockSize+1,jBlockSize+1)%T = Blocks(b%NECorner%neighborLocalBlock)%Points(2,2)%T
        else! if (b%SWCorner%neighborBlock /= -1) then
          destination = b%NECorner%neighborProc
          buffer = Blocks(b%SWCorner%neighborLocalBlock)%Points(iBlockSize-1,jBlockSize-1)%T
          tag = nB + eB + b%NECorner%neighborBlock * 1000
          call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
        end if
      end if

      ! South east corner
      if (b%SECorner%BC == -1) then
        if (MyID == b%SECorner%neighborProc) then
          b%Points(iBlockSize+1,0)%T = Blocks(b%SECorner%neighborLocalBlock)%Points(2,jBlockSize-1)%T
        else! if (b%NWCorner%neighborBlock /= -1) then
          destination = b%SECorner%neighborProc
          buffer = Blocks(b%NWCorner%neighborLocalBlock)%Points(iBlockSize-1,2)%T
          tag = sB + eB + b%SECorner%neighborBlock * 1000
          call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
        end if
      end if

      ! South west corner
      if (b%SWCorner%BC == -1) then
        if (MyID == b%SWCorner%neighborProc) then
          b%Points(0,0)%T = Blocks(b%SWCorner%neighborLocalBlock)%Points(iBlockSize-1,jBlockSize-1)%T
        else! if (b%NECorner%neighborBlock /= -1) then
          destination = b%SWCorner%neighborProc
          buffer = Blocks(b%NECorner%neighborLocalBlock)%Points(2,2)%T
          tag = sB + wB + b%SWCorner%neighborBlock * 1000
          call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
        end if
      end if

      ! North west corner
      if (b%NWCorner%BC == -1) then
        if (MyID == b%NWCorner%neighborProc) then
          b%Points(0,jBlockSize+1)%T = Blocks(b%NWCorner%neighborLocalBlock)%Points(iBlockSize-1,2)%T
        else! if (b%SECorner%neighborBlock /= -1) then
          destination = b%NWCorner%neighborProc
          buffer = Blocks(b%SECorner%neighborLocalBlock)%Points(2,jBlockSize-1)%T
          tag = nB + wB + b%NWCorner%neighborBlock * 1000
          call MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, request, ierror)
        end if
      end if

    end do

  end subroutine

  subroutine cross_proc_talk
    type (BlockType), pointer :: b
    integer :: n_, i, j, tag, destination, source, request
    real(kind = 8) :: buffer
    integer :: status(MPI_STATUS_SIZE)

    ! For each block...
    do n_ = 1, MyNBlocks
      b => Blocks(n_)

      ! Faces
      ! South Face Blocking Receive
      if (b%southFace%BC == -1) then
        if (MyID /= b%southFace%neighborProc) then
          ! Our neighbor block is on different proc, we must receive with MPI.
          ! Blocking receive.

          source = b%southFace%neighborProc
          do i = 1, iBlockSize
            tag = nB * 100 + i + b%id * 1000
            call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
            b%Points(i, 0)%T = buffer
          end do
        end if
      end if

      ! North Face Blocking Receive
      if (b%northFace%BC == -1) then
        if (MyID /= b%northFace%neighborProc) then
          source = b%northFace%neighborProc
          do i = 1, iBlockSize
            tag = sB * 100 + i + b%id * 1000
            call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
            b%Points(i, jBlockSize+1)%T = buffer
          end do
        end if
      end if

      ! West Face Blocking Receive
      if (b%westFace%BC == -1) then
        if (MyID /= b%westFace%neighborProc) then
          source = b%westFace%neighborProc
          do j = 1, jBlockSize
            tag = eB * 100 + j + b%id * 1000
            call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
            b%Points(0, j)%T = buffer
          end do
        end if
      end if

      ! East Face Blocking Receive
      if (b%eastFace%BC == -1) then
        if (MyID /= b%eastFace%neighborProc) then
          source = b%eastFace%neighborProc
          do j = 1, jBlockSize
            tag = wB * 100 + j + b%id * 1000
            call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
            b%Points(iBlockSize+1, j)%T = buffer
          end do
        end if
      end if

      ! Corners
      ! South West corner
      if (b%SWCorner%BC == -1) then
        if (MyID /= b%SWCorner%neighborProc) then
          source = b%SWCorner%neighborProc
          tag = nB + eB + b%id * 1000
          call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
          if (buffer /= 0.d0) b%Points(0,0)%T = buffer
        end if
      end if

      ! North West corner
      if (b%NWCorner%BC == -1) then
        if (MyID /= b%NWCorner%neighborProc) then
          source = b%NWCorner%neighborProc
          tag = sB + eB + b%id * 1000
          call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
          if (buffer /= 0.d0) b%Points(0,jBlockSize+1)%T = buffer
        end if
      end if

      ! North East corner
      if (b%NECorner%BC == -1) then
        if (MyID /= b%NECorner%neighborProc) then
          source = b%NECorner%neighborProc
          tag = sB + wB + b%id * 1000
          call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
          if (buffer /= 0.d0) b%Points(iBlockSize+1,jBlockSize+1)%T = buffer
        end if
      end if

      ! South East corner
      if (b%SECorner%BC == -1) then
        if (MyID /= b%SECorner%neighborProc) then
          source = b%SECorner%neighborProc
          tag = nB + wB + b%id * 1000
          call MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, status, ierror)
          if (buffer /= 0.d0) b%Points(iBlockSize+1,0)%T = buffer
        end if
      end if

    end do
  end subroutine
end module MainRoutines
