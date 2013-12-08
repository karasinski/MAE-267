module MainRoutines
  use BlockModule
  use plot3D_module
  use GridCreation

  implicit none

  ! Block array for solver.
  type (BlockType), pointer :: Blocks(:)

  ! Iteration lists
  type(LinkedList), pointer :: northLocal, southLocal, eastLocal, westLocal
  type(LinkedList), pointer :: northMPI, southMPI, eastMPI, westMPI
  type(LinkedList), pointer :: neLocal, nwLocal, swLocal, seLocal
  type(LinkedList), pointer :: neMPI, nwMPI, swMPI, seMPI

  ! Linked lists
  type(LinkedList), pointer :: northLocalList, southLocalList, eastLocalList, westLocalList
  type(LinkedList), pointer :: northMPIList, southMPIList, eastMPIList, westMPIList
  type(LinkedList), pointer :: neLocalList, nwLocalList, swLocalList, seLocalList
  type(LinkedList), pointer :: neMPIList, nwMPIList, swMPIList, seMPIList

  ! Temp lists
  type(LinkedList), pointer :: northLocalTemp, southLocalTemp, eastLocalTemp, westLocalTemp
  type(LinkedList), pointer :: northMPITemp, southMPITemp, eastMPITemp, westMPITemp
  type(LinkedList), pointer :: neLocalTemp, nwLocalTemp, swLocalTemp, seLocalTemp
  type(LinkedList), pointer :: neMPITemp, nwMPITemp, swMPITemp, seMPITemp

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

    ! Set up our linked lists.
    call initialize_linked_lists
  end subroutine

  ! This is the main solver.
  subroutine solve(Blocks)
    type (BlockType), target :: Blocks(:)
    real(kind=8) :: temp_residual = 1.d0, residual = 1.d0, local_residual = 1.d0
    real(kind=8) :: residuals(max_steps) = 0.d0, tol = .00001d0
    integer :: n_

    ! Processor 0 times our iterations until convergence.
    if (MyID == 0) then
      call start_clock
    end if

    !  Begin main loop, stop if we hit our mark or after max_steps iterations.
    do step = 1, max_steps
      
      ! Calculate our first and second derivatives for all our points, and
      ! calculate the new temperature for all of our interior points.
      call derivatives

      ! After each iteration of temperature updates we need to update our ghost nodes.
      call update_local_ghosts
      call mpi_sends
      call mpi_receives

      ! Find block with largest residual.
      local_residual = 0.d0

      do n_ = 1, MyNBlocks
        temp_residual = maxval(abs(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))

        if (temp_residual > local_residual) then
          local_residual = temp_residual
        end if
      end do

      call MPI_Barrier(mpi_comm_world, ierror)
      call mpi_allreduce(local_residual, residual, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpi_comm_world, ierror)
      residuals(step) = residual
      ! write(*,*), MyID, step, residual, local_residual

      if (residual < tol) exit
    end do

    ! We have converged and all procs have stopped, end clock.
    call MPI_Barrier(mpi_comm_world)
    if (MyID == 0) then
      call end_clock
    end if

    ! Write some results to file/screen.
    call output(Blocks, residuals)
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

  ! Our neighbor block is on same proc, we can grab directly.
  subroutine update_local_ghosts
    type (BlockType), pointer :: b1
    real(kind=8), pointer :: p1, p2
    integer :: i, j
    
    ! North face ghost nodes
    northLocal => northLocalList
    do
        if (.NOT. associated(northLocal)) exit
        do i = 1, iBlockSize
          b1 => Blocks(northLocal%id)
          p1 => b1%Points(i, jBlockSize+1)%T
          p2 => Blocks(b1%northFace%neighborLocalBlock)%Points(i, 2)%T
          p1 = p2
        end do
        northLocal => northLocal%next
    end do

    ! south face ghost nodes
    southLocal => southLocalList
    do
        if (.NOT. associated(southLocal)) exit
        do i = 1, iBlockSize
          b1 => Blocks(southLocal%id)
          p1 => b1%Points(i, 0)%T
          p2 => Blocks(b1%southFace%neighborLocalBlock)%Points(i, jBlockSize-1)%T
          p1 = p2
        end do
        southLocal => southLocal%next
    end do

    ! east face ghost nodes
    eastLocal => eastLocalList
    do
        if (.NOT. associated(eastLocal)) exit
        do j = 1, jBlockSize
          b1 => Blocks(eastLocal%id)
          p1 => b1%Points(iBlockSize+1, j)%T
          p2 => Blocks(b1%eastFace%neighborLocalBlock)%Points(2, j)%T
          p1 = p2
        end do
        eastLocal => eastLocal%next
    end do

    ! west face ghost nodes
    westLocal => westLocalList
    do
        if (.NOT. associated(westLocal)) exit
        do j = 1, jBlockSize
          b1 => Blocks(westLocal%id)
          p1 => b1%Points(0, j)%T
          p2 => Blocks(b1%westFace%neighborLocalBlock)%Points(iBlockSize-1, j)%T
          p1 = p2
        end do
        westLocal => westLocal%next
    end do

    ! ne corner ghost node
    neLocal => neLocalList
    do
      if (.NOT. associated(neLocal)) exit
      b1 => Blocks(neLocal%id)
      p1 => b1%Points(iBlockSize+1,jBlockSize+1)%T
      p2 => Blocks(b1%NECorner%neighborLocalBlock)%Points(2,2)%T
      p1 = p2
      neLocal => neLocal%next
    end do

    ! se corner ghost node
    seLocal => seLocalList
    do
      if (.NOT. associated(seLocal)) exit
      b1 => Blocks(seLocal%id)
      p1 => b1%Points(iBlockSize+1,0)%T
      p2 => Blocks(b1%SECorner%neighborLocalBlock)%Points(2,jBlockSize-1)%T
      p1 = p2
      seLocal => seLocal%next
    end do

    ! sw corner ghost node
    swLocal => swLocalList
    do
      if (.NOT. associated(swLocal)) exit
      b1 => Blocks(swLocal%id)
      p1 => b1%Points(0,0)%T
      p2 => Blocks(b1%SWCorner%neighborLocalBlock)%Points(iBlockSize-1,jBlockSize-1)%T
      p1 = p2
      swLocal => swLocal%next
    end do

    ! nw corner ghost node
    nwLocal => nwLocalList
    do
      if (.NOT. associated(nwLocal)) exit
      b1 => Blocks(nwLocal%id)
      p1 => b1%Points(0,jBlockSize+1)%T
      p2 => Blocks(b1%NWCorner%neighborLocalBlock)%Points(iBlockSize-1,2)%T
      p1 = p2
      nwLocal => nwLocal%next
    end do

  end subroutine

subroutine mpi_sends
    type (BlockType), pointer :: b1
    real(kind=8), pointer :: p1
    integer :: i, j, tag, destination
    real(kind = 8) :: i_buffer(iBlockSize), j_buffer(jBlockSize)
    real(kind = 8) :: buffer
    
    northMPI => northMPIList
    do
      if (.NOT. associated(northMPI)) exit
      ! Our neighbor block is on different proc, so it will also need information from
      ! this block. We do a nonblocking send now and a blocking receive later.
      b1 => Blocks(northMPI%id)

      ! Pack our values to send into a buffer.
      do i = 1, iBlockSize
        p1 => b1%Points(i, jBlockSize-1)%T
        i_buffer(i) = p1
      end do
    
      ! Find the destination.
      destination = b1%northFace%neighborProc

      ! Generate a tag unique within the iteration.
      tag = nB + b1%northFace%neighborBlock * 1000          

      ! Send everything to the proc and continue without immediate confirmation.
      call MPI_Isend(i_buffer, iBlockSize, MPI_DOUBLE_PRECISION, destination, tag, &
                     mpi_comm_world, request, ierror)

      northMPI => northMPI%next
    end do

    southMPI => southMPIList
    do
      if (.NOT. associated(southMPI)) exit
      b1 => Blocks(southMPI%id)
      do i = 1, iBlockSize
        p1 => b1%Points(i, 2)%T
        i_buffer(i) = p1
      end do
      destination = b1%southFace%neighborProc
      tag = sB + b1%southFace%neighborBlock * 1000          
      call MPI_Isend(i_buffer, iBlockSize, MPI_DOUBLE_PRECISION, destination, tag, &
                     mpi_comm_world, request, ierror)
      southMPI => southMPI%next
    end do

    eastMPI => eastMPIList
    do
      if (.NOT. associated(eastMPI)) exit
      b1 => Blocks(eastMPI%id)
      do j = 1, jBlockSize
        p1 => b1%Points(iBlockSize-1, j)%T
        j_buffer(j) = p1
      end do
      destination = b1%eastFace%neighborProc
      tag = eB + b1%eastFace%neighborBlock * 1000          
      call MPI_Isend(j_buffer, jBlockSize, MPI_DOUBLE_PRECISION, destination, tag, &
                     mpi_comm_world, request, ierror)
      eastMPI => eastMPI%next
    end do

    westMPI => westMPIList
    do
      if (.NOT. associated(westMPI)) exit
      b1 => Blocks(westMPI%id)
      do j = 1, jBlockSize
        p1 => b1%Points(2, j)%T
        j_buffer(j) = p1
      end do
      destination = b1%westFace%neighborProc
      tag = wB + b1%westFace%neighborBlock * 1000          
      call MPI_Isend(j_buffer, jBlockSize, MPI_DOUBLE_PRECISION, destination, tag, &
                     mpi_comm_world, request, ierror)
      westMPI => westMPI%next
    end do

    neMPI => neMPIList
    do
      if (.NOT. associated(neMPI)) exit
      b1 => Blocks(neMPI%id)
      p1 => b1%Points(iBlockSize-1,jBlockSize-1)%T
      buffer = p1
      destination = b1%NECorner%neighborProc
      tag = nB + eB * 10 + b1%NECorner%neighborBlock * 1000
      call MPI_Isend(buffer, 1, MPI_DOUBLE_PRECISION, destination, tag, mpi_comm_world, &
                     request, ierror)
      neMPI => neMPI%next
    end do

    seMPI => seMPIList
    do
      if (.NOT. associated(seMPI)) exit
      b1 => Blocks(seMPI%id)
      p1 => b1%Points(iBlockSize-1,2)%T
      buffer = p1
      destination = b1%SECorner%neighborProc
      tag = sB + eB * 10 + b1%SECorner%neighborBlock * 1000
      call MPI_Isend(buffer, 1, MPI_DOUBLE_PRECISION, destination, tag, mpi_comm_world, &
                     request, ierror)
      seMPI => seMPI%next
    end do

    swMPI => swMPIList
    do
      if (.NOT. associated(swMPI)) exit
      b1 => Blocks(swMPI%id)
      p1 => b1%Points(2,2)%T
      buffer = p1
      destination = b1%SWCorner%neighborProc
      tag = sB + wB * 10 + b1%SWCorner%neighborBlock * 1000
      call MPI_Isend(buffer, 1, MPI_DOUBLE_PRECISION, destination, tag, mpi_comm_world, &
                     request, ierror)
      swMPI => swMPI%next
    end do

    nwMPI => nwMPIList
    do
      if (.NOT. associated(nwMPI)) exit
      b1 => Blocks(nwMPI%id)
      p1 => b1%Points(2,jBlockSize-1)%T
      buffer = p1
      destination = b1%nwCorner%neighborProc
      tag = nB + wB * 10 + b1%NWCorner%neighborBlock * 1000
      call MPI_Isend(buffer, 1, MPI_DOUBLE_PRECISION, destination, tag, mpi_comm_world, &
                     request, ierror)
      nwMPI => nwMPI%next
    end do

  end subroutine


  subroutine mpi_receives
    type (BlockType), pointer :: b1
    real(kind=8), pointer :: p1
    integer :: i, j, tag, source
    real(kind = 8) :: i_buffer(iBlockSize), j_buffer(jBlockSize)
    real(kind = 8) :: buffer

      ! Our neighbor block is on different proc, we must receive with MPI.
      southMPI => southMPIList
      do
        if (.NOT. associated(southMPI)) exit
        b1 => Blocks(southMPI%id)
        source = b1%southFace%neighborProc
        tag = nB + b1%id * 1000
        call MPI_RECV(i_buffer, iBlockSize, MPI_DOUBLE_PRECISION, source, tag, mpi_comm_world, &
                      status, ierror)
        do i = 1, iBlockSize
          p1 => b1%Points(i, 0)%T
          p1 = i_buffer(i)
        end do
        southMPI => southMPI%next
      end do

      northMPI => northMPIList
      do
        if (.NOT. associated(northMPI)) exit
        b1 => Blocks(northMPI%id)
        source = b1%northFace%neighborProc
        tag = sB + b1%id * 1000
        call MPI_RECV(i_buffer, iBlockSize, MPI_DOUBLE_PRECISION, source, tag, mpi_comm_world, &
                      status, ierror)
        do i = 1, iBlockSize
          p1 => b1%Points(i, jBlockSize+1)%T
          p1 = i_buffer(i)
        end do
        northMPI => northMPI%next
      end do

      westMPI => westMPIList
      do
        if (.NOT. associated(westMPI)) exit
        b1 => Blocks(westMPI%id)
        source = b1%westFace%neighborProc
        tag = eB + b1%id * 1000
        call MPI_RECV(j_buffer, jBlockSize, MPI_DOUBLE_PRECISION, source, tag, mpi_comm_world, &
                      status, ierror)
        do j = 1, jBlockSize
          p1 => b1%Points(0, j)%T
          p1 = j_buffer(j)
        end do
        westMPI => westMPI%next
      end do

      eastMPI => eastMPIList
      do
        if (.NOT. associated(eastMPI)) exit
        b1 => Blocks(eastMPI%id)
        source = b1%eastFace%neighborProc
        tag = wB + b1%id * 1000
        call MPI_RECV(j_buffer, jBlockSize, MPI_DOUBLE_PRECISION, source, tag, mpi_comm_world, &
                      status, ierror)
        do j = 1, jBlockSize
          p1 => b1%Points(iBlockSize+1, j)%T
          p1 = j_buffer(j)
        end do
        eastMPI => eastMPI%next
      end do

      swMPI => swMPIList
      do
        if (.NOT. associated(swMPI)) exit
        b1 => Blocks(swMPI%id)
        source = b1%SWCorner%neighborProc
        tag = nB + eB * 10 + b1%id * 1000
        call MPI_RECV(buffer, 1, MPI_DOUBLE_PRECISION, source, tag, mpi_comm_world, status, ierror)
        p1 => b1%Points(0,0)%T
        p1 = buffer
        swMPI => swMPI%next
      end do

      nwMPI => nwMPIList
      do
        if (.NOT. associated(nwMPI)) exit
        b1 => Blocks(nwMPI%id)
        source = b1%NWCorner%neighborProc
        tag = sB + eB * 10 + b1%id * 1000
        call MPI_RECV(buffer, 1, MPI_DOUBLE_PRECISION, source, tag, mpi_comm_world, status, ierror)
        p1 => b1%Points(0,jBlockSize+1)%T
        p1 = buffer
        nwMPI => nwMPI%next
      end do

      neMPI => neMPIList
      do
        if (.NOT. associated(neMPI)) exit
        b1 => Blocks(neMPI%id)
        source = b1%NECorner%neighborProc
        tag = sB + wB * 10 + b1%id * 1000
        call MPI_RECV(buffer, 1, MPI_DOUBLE_PRECISION, source, tag, mpi_comm_world, status, ierror)
        p1 => b1%Points(iBlockSize+1,jBlockSize+1)%T
        p1 = buffer
        neMPI => neMPI%next
      end do

      seMPI => seMPIList
      do
        if (.NOT. associated(seMPI)) exit
        b1 => Blocks(seMPI%id)
        source = b1%SECorner%neighborProc
        tag = nB + wB * 10 + b1%id * 1000
        call MPI_RECV(buffer, 1, MPI_DOUBLE_PRECISION, source, tag, mpi_comm_world, status, ierror)
        p1 => b1%Points(iBlockSize+1,0)%T
        p1 = buffer
        seMPI => seMPI%next
      end do

  end subroutine

  subroutine initialize_linked_lists
    type (BlockType), pointer :: b
    integer :: n_
    
    ! For each block...
    do n_ = 1, MyNBlocks
      b => Blocks(n_)

      ! north Face
      if (b%northFace%BC == INTERNAL_BOUNDARY) then
        if (.not. associated(northLocalList)) then
          allocate(northLocalList)
          northLocalTemp => northLocalList
          nullify(northLocalTemp%next)
          northLocalTemp%id = n_
        else
          allocate(northLocalTemp%next)
          northLocalTemp => northLocalTemp%next
          nullify(northLocalTemp%next)
          northLocalTemp%id = n_
        end if
      else if (b%northFace%BC == PROC_BOUNDARY) then
        if (.not. associated(northMPIList)) then
          allocate(northMPIList)
          northMPITemp => northMPIList
          nullify(northMPITemp%next)
          northMPITemp%id = n_
        else
          allocate(northMPITemp%next)
          northMPITemp => northMPITemp%next
          nullify(northMPITemp%next)
          northMPITemp%id = n_
        end if
      end if

      ! south Face
      if (b%southFace%BC == INTERNAL_BOUNDARY) then
        if (.not. associated(southLocalList)) then
          allocate(southLocalList)
          southLocalTemp => southLocalList
          nullify(southLocalTemp%next)
          southLocalTemp%id = n_
        else
          allocate(southLocalTemp%next)
          southLocalTemp => southLocalTemp%next
          nullify(southLocalTemp%next)
          southLocalTemp%id = n_
        end if
      else if (b%southFace%BC == PROC_BOUNDARY) then
        if (.not. associated(southMPIList)) then
          allocate(southMPIList)
          southMPITemp => southMPIList
          nullify(southMPITemp%next)
          southMPITemp%id = n_
        else
          allocate(southMPITemp%next)
          southMPITemp => southMPITemp%next
          nullify(southMPITemp%next)
          southMPITemp%id = n_
        end if
      end if

      ! east Face
      if (b%eastFace%BC == INTERNAL_BOUNDARY) then
        if (.not. associated(eastLocalList)) then
          allocate(eastLocalList)
          eastLocalTemp => eastLocalList
          nullify(eastLocalTemp%next)
          eastLocalTemp%id = n_
        else
          allocate(eastLocalTemp%next)
          eastLocalTemp => eastLocalTemp%next
          nullify(eastLocalTemp%next)
          eastLocalTemp%id = n_
        end if
      else if (b%eastFace%BC == PROC_BOUNDARY) then
        if (.not. associated(eastMPIList)) then
          allocate(eastMPIList)
          eastMPITemp => eastMPIList
          nullify(eastMPITemp%next)
          eastMPITemp%id = n_
        else
          allocate(eastMPITemp%next)
          eastMPITemp => eastMPITemp%next
          nullify(eastMPITemp%next)
          eastMPITemp%id = n_
        end if
      end if

      ! west Face
      if (b%westFace%BC == INTERNAL_BOUNDARY) then
        if (.not. associated(westLocalList)) then
          allocate(westLocalList)
          westLocalTemp => westLocalList
          nullify(westLocalTemp%next)
          westLocalTemp%id = n_
        else
          allocate(westLocalTemp%next)
          westLocalTemp => westLocalTemp%next
          nullify(westLocalTemp%next)
          westLocalTemp%id = n_
        end if
      else if (b%westFace%BC == PROC_BOUNDARY) then
        if (.not. associated(westMPIList)) then
          allocate(westMPIList)
          westMPITemp => westMPIList
          nullify(westMPITemp%next)
          westMPITemp%id = n_
        else
          allocate(westMPITemp%next)
          westMPITemp => westMPITemp%next
          nullify(westMPITemp%next)
          westMPITemp%id = n_
        end if
      end if

      ! ne Corner
      if (b%NECorner%BC == INTERNAL_BOUNDARY) then
        if (.not. associated(neLocalList)) then
          allocate(neLocalList)
          neLocalTemp => neLocalList
          nullify(neLocalTemp%next)
          neLocalTemp%id = n_
        else
          allocate(neLocalTemp%next)
          neLocalTemp => neLocalTemp%next
          nullify(neLocalTemp%next)
          neLocalTemp%id = n_
        end if
      else if (b%NECorner%BC == PROC_BOUNDARY) then
        if (.not. associated(neMPIList)) then
          allocate(neMPIList)
          neMPITemp => neMPIList
          nullify(neMPITemp%next)
          neMPITemp%id = n_
        else
          allocate(neMPITemp%next)
          neMPITemp => neMPITemp%next
          nullify(neMPITemp%next)
          neMPITemp%id = n_
        end if
      end if

    ! nw Corner
    if (b%NWCorner%BC == INTERNAL_BOUNDARY) then
      if (.not. associated(nwLocalList)) then
        allocate(nwLocalList)
        nwLocalTemp => nwLocalList
        nullify(nwLocalTemp%next)
        nwLocalTemp%id = n_
      else
        allocate(nwLocalTemp%next)
        nwLocalTemp => nwLocalTemp%next
        nullify(nwLocalTemp%next)
        nwLocalTemp%id = n_
      end if
    else if (b%NWCorner%BC == PROC_BOUNDARY) then
      if (.not. associated(nwMPIList)) then
        allocate(nwMPIList)
        nwMPITemp => nwMPIList
        nullify(nwMPITemp%next)
        nwMPITemp%id = n_
      else
        allocate(nwMPITemp%next)
        nwMPITemp => nwMPITemp%next
        nullify(nwMPITemp%next)
        nwMPITemp%id = n_
      end if
    end if

    ! sw Corner
    if (b%swCorner%BC == INTERNAL_BOUNDARY) then
      if (.not. associated(swLocalList)) then
        allocate(swLocalList)
        swLocalTemp => swLocalList
        nullify(swLocalTemp%next)
        swLocalTemp%id = n_
      else
        allocate(swLocalTemp%next)
        swLocalTemp => swLocalTemp%next
        nullify(swLocalTemp%next)
        swLocalTemp%id = n_
      end if
    else if (b%swCorner%BC == PROC_BOUNDARY) then
      if (.not. associated(swMPIList)) then
        allocate(swMPIList)
        swMPITemp => swMPIList
        nullify(swMPITemp%next)
        swMPITemp%id = n_
      else
        allocate(swMPITemp%next)
        swMPITemp => swMPITemp%next
        nullify(swMPITemp%next)
        swMPITemp%id = n_
      end if
    end if

    ! se Corner
    if (b%seCorner%BC == INTERNAL_BOUNDARY) then
      if (.not. associated(seLocalList)) then
        allocate(seLocalList)
        seLocalTemp => seLocalList
        nullify(seLocalTemp%next)
        seLocalTemp%id = n_
      else
        allocate(seLocalTemp%next)
        seLocalTemp => seLocalTemp%next
        nullify(seLocalTemp%next)
        seLocalTemp%id = n_
      end if
    else if (b%seCorner%BC == PROC_BOUNDARY) then
      if (.not. associated(seMPIList)) then
        allocate(seMPIList)
        seMPITemp => seMPIList
        nullify(seMPITemp%next)
        seMPITemp%id = n_
      else
        allocate(seMPITemp%next)
        seMPITemp => seMPITemp%next
        nullify(seMPITemp%next)
        seMPITemp%id = n_
      end if
    end if

    end do
  end subroutine

end module MainRoutines
