MODULE MainRoutines
  USE BlockModule
  USE plot3D_module
  USE GridCreation

  IMPLICIT NONE

  ! Block array for solver.
  TYPE (BlockType), POINTER :: Blocks(:)

  ! Iteration lists
  TYPE(LinkedList), POINTER :: northLocal, southLocal, eastLocal, westLocal
  TYPE(LinkedList), POINTER :: northMPI, southMPI, eastMPI, westMPI
  TYPE(LinkedList), POINTER :: neLocal, nwLocal, swLocal, seLocal
  TYPE(LinkedList), POINTER :: neMPI, nwMPI, swMPI, seMPI

  ! Linked lists
  TYPE(LinkedList), POINTER :: northLocalList, southLocalList, eastLocalList, westLocalList
  TYPE(LinkedList), POINTER :: northMPIList, southMPIList, eastMPIList, westMPIList
  TYPE(LinkedList), POINTER :: neLocalList, nwLocalList, swLocalList, seLocalList
  TYPE(LinkedList), POINTER :: neMPIList, nwMPIList, swMPIList, seMPIList

CONTAINS

  ! This subroutine partitions the grid into blocks.
  SUBROUTINE initialize_grid(b)
    TYPE (BlockType) :: b(:)

    ! Create blocks and find neighbors.
    CALL create_blocks(b)

    ! Set up the points in each block.
    CALL initialize_block_grid(b)

    ! Set initial temperatures.
    CALL initialize_block_temp(b)

    ! Calculate proper bounds for each block.
    CALL set_bounds(b)

    WRITE(*,*), 'Processor ', MyID, ' initialized simulation.'
    WRITE(*,*)
  END SUBROUTINE

  ! This subroutine reads the block files and initializes the simulation.
  SUBROUTINE initialization
    ! We then read in the grid file.
    CALL read_grid_file(Blocks)

    ! We then read in the initial temperature file.
    CALL read_temp_file(Blocks)

    !  Initialize the points.
    CALL initialize_points(Blocks)

    !  Initialize the primary face areas and volumes.
    CALL initialize_faces_and_volumes(Blocks)

    ! Calculate constants for integration.
    CALL set_constants(Blocks)

    ! Set up our linked lists.
    CALL initialize_linked_lists
  END SUBROUTINE

  ! This is the main solver.
  SUBROUTINE solve(Blocks)
    TYPE (BlockType), TARGET :: Blocks(:)
    REAL(KIND=8) :: temp_residual = 1.d0, residual = 1.d0, local_residual = 1.d0
    REAL(KIND=8) :: residuals(max_steps) = 0.d0, tol = .00001d0
    INTEGER :: n_

    ! Processor 0 times our iterations until convergence.
    IF (MyID == 0) THEN
      CALL start_clock
    END IF

    !  Begin main loop, stop if we hit our mark or after max_steps iterations.
    DO step = 1, max_steps

      ! Calculate our first and second derivatives for all our points, and
      ! calculate the new temperature for all of our interior points.
      CALL derivatives

      ! After each iteration of temperature updates we need to update our ghost nodes.
      CALL update_local_ghosts
      CALL mpi_sends
      CALL mpi_receives

      ! Find block with largest residual.
      local_residual = 0.d0

      DO n_ = 1, MyNBlocks
        temp_residual = MAXVAL(ABS(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))

        IF (temp_residual > local_residual) THEN
          local_residual = temp_residual
        END IF
      END DO

      CALL mpi_allreduce(local_residual, residual, 1, MPI_REAL8, MPI_MAX, mpi_comm_world, ierror)
      residuals(step) = residual
      !       write(*,*), MyID, step, residual, local_residual

      IF (residual < tol) EXIT
    END DO

    ! We have converged and all procs have stopped, end clock.
    CALL MPI_Barrier(mpi_comm_world, ierror)
    IF (MyID == 0) THEN
      CALL end_clock
    END IF

    ! Write some results to file/screen.
    CALL output(Blocks, residuals)
  END SUBROUTINE

  SUBROUTINE derivatives
    TYPE (BlockType), POINTER :: MyBlock
    TYPE (GridPoint), POINTER :: p0, p1, p2, p3
    REAL (KIND=8), POINTER :: tempT(:,:)
    REAL (KIND=8) :: dTdx, dTdy
    INTEGER :: i, j, n_

    ! Loop over each block.
    DO n_ = 1, MyNBlocks
      MyBlock => Blocks(n_)
      tempT => MyBlock%Points%tempT

      ! Reset the change in temperature to zero before we begin summing again.
      tempT = 0.d0

      DO j = MyBlock%localJMIN, MyBlock%localJMAX
        DO i =  MyBlock%localIMIN, MyBlock%localIMAX
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
        END DO
      END DO

      ! Update temperatures in the block.
      DO j = MyBlock%lowJTemp, MyBlock%localJMAX
        DO i =  MyBlock%lowITemp, MyBlock%localIMAX
          p0 => Blocks(n_)%Points(i,j)
          p0%T = p0%T + p0%tempT
        END DO
      END DO
    END DO
  END SUBROUTINE

  ! Our neighbor block is on same proc, we can grab directly.
  SUBROUTINE update_local_ghosts
    TYPE (BlockType), POINTER :: b1
    REAL(KIND=8), POINTER :: p1, p2
    INTEGER :: i, j

    ! North face ghost nodes
    northLocal => northLocalList
    DO
      IF (.NOT. ASSOCIATED(northLocal)) EXIT
      DO i = 1, iBlockSize
        b1 => Blocks(northLocal%ID)
        p1 => b1%Points(i, jBlockSize+1)%T
        p2 => Blocks(b1%northFace%neighborLocalBlock)%Points(i, 2)%T
        p1 = p2
      END DO
      northLocal => northLocal%next
    END DO

    ! South face ghost nodes
    southLocal => southLocalList
    DO
      IF (.NOT. ASSOCIATED(southLocal)) EXIT
      DO i = 1, iBlockSize
        b1 => Blocks(southLocal%ID)
        p1 => b1%Points(i, 0)%T
        p2 => Blocks(b1%southFace%neighborLocalBlock)%Points(i, jBlockSize-1)%T
        p1 = p2
      END DO
      southLocal => southLocal%next
    END DO

    ! East face ghost nodes
    eastLocal => eastLocalList
    DO
      IF (.NOT. ASSOCIATED(eastLocal)) EXIT
      DO j = 1, jBlockSize
        b1 => Blocks(eastLocal%ID)
        p1 => b1%Points(iBlockSize+1, j)%T
        p2 => Blocks(b1%eastFace%neighborLocalBlock)%Points(2, j)%T
        p1 = p2
      END DO
      eastLocal => eastLocal%next
    END DO

    ! West face ghost nodes
    westLocal => westLocalList
    DO
      IF (.NOT. ASSOCIATED(westLocal)) EXIT
      DO j = 1, jBlockSize
        b1 => Blocks(westLocal%ID)
        p1 => b1%Points(0, j)%T
        p2 => Blocks(b1%westFace%neighborLocalBlock)%Points(iBlockSize-1, j)%T
        p1 = p2
      END DO
      westLocal => westLocal%next
    END DO

    ! NE corner ghost node
    neLocal => neLocalList
    DO
      IF (.NOT. ASSOCIATED(neLocal)) EXIT
      b1 => Blocks(neLocal%ID)
      p1 => b1%Points(iBlockSize+1,jBlockSize+1)%T
      p2 => Blocks(b1%NECorner%neighborLocalBlock)%Points(2,2)%T
      p1 = p2
      neLocal => neLocal%next
    END DO

    ! SE corner ghost node
    seLocal => seLocalList
    DO
      IF (.NOT. ASSOCIATED(seLocal)) EXIT
      b1 => Blocks(seLocal%ID)
      p1 => b1%Points(iBlockSize+1,0)%T
      p2 => Blocks(b1%SECorner%neighborLocalBlock)%Points(2,jBlockSize-1)%T
      p1 = p2
      seLocal => seLocal%next
    END DO

    ! SW corner ghost node
    swLocal => swLocalList
    DO
      IF (.NOT. ASSOCIATED(swLocal)) EXIT
      b1 => Blocks(swLocal%ID)
      p1 => b1%Points(0,0)%T
      p2 => Blocks(b1%SWCorner%neighborLocalBlock)%Points(iBlockSize-1,jBlockSize-1)%T
      p1 = p2
      swLocal => swLocal%next
    END DO

    ! NW corner ghost node
    nwLocal => nwLocalList
    DO
      IF (.NOT. ASSOCIATED(nwLocal)) EXIT
      b1 => Blocks(nwLocal%ID)
      p1 => b1%Points(0,jBlockSize+1)%T
      p2 => Blocks(b1%NWCorner%neighborLocalBlock)%Points(iBlockSize-1,2)%T
      p1 = p2
      nwLocal => nwLocal%next
    END DO
  END SUBROUTINE

  ! Our neighbor block is on different proc, so it will also need information from
  ! this block. We do a nonblocking send now and a blocking receive later.
  SUBROUTINE mpi_sends
    TYPE (BlockType), POINTER :: b1
    REAL(KIND=8), POINTER :: p1
    INTEGER :: i, j, tag, destination
    REAL(KIND = 8) :: i_buffer(iBlockSize), j_buffer(jBlockSize)
    REAL(KIND = 8) :: buffer

    ! Set pointers to first elements in linked lists.
    northMPI => northMPIList
    southMPI => southMPIList
    eastMPI  => eastMPIList
    westMPI  => westMPIList
    neMPI    => neMPIList
    seMPI    => seMPIList
    swMPI    => swMPIList
    nwMPI    => nwMPIList

    ! North face MPI sends
    DO
      ! If it's associated, continue, otherwise we exit.
      IF (.NOT. ASSOCIATED(northMPI)) EXIT

      ! Identify our block from our linked list.
      b1 => Blocks(northMPI%ID)

      ! Pack our values to send into a buffer.
      DO i = 1, iBlockSize
        p1 => b1%Points(i, jBlockSize-1)%T
        i_buffer(i) = p1
      END DO

      ! Find the destination.
      destination = b1%northFace%neighborProc

      ! Generate a tag unique within the iteration and interproc communication.
      tag = nB

      ! Send everything to the proc and continue without immediate confirmation.
      CALL MPI_Isend(i_buffer, iBlockSize, MPI_REAL8, destination, tag, &
        mpi_comm_world, request, ierror)

      ! Traverse our linked list.
      northMPI => northMPI%next
    END DO

    ! South face MPI sends
    DO
      IF (.NOT. ASSOCIATED(southMPI)) EXIT
      b1 => Blocks(southMPI%ID)
      DO i = 1, iBlockSize
        p1 => b1%Points(i, 2)%T
        i_buffer(i) = p1
      END DO
      destination = b1%southFace%neighborProc
      tag = sB
      CALL MPI_Isend(i_buffer, iBlockSize, MPI_REAL8, destination, tag, &
        mpi_comm_world, request, ierror)
      southMPI => southMPI%next
    END DO

    ! East face MPI sends
    DO
      IF (.NOT. ASSOCIATED(eastMPI)) EXIT
      b1 => Blocks(eastMPI%ID)
      DO j = 1, jBlockSize
        p1 => b1%Points(iBlockSize-1, j)%T
        j_buffer(j) = p1
      END DO
      destination = b1%eastFace%neighborProc
      tag = eB
      CALL MPI_Isend(j_buffer, jBlockSize, MPI_REAL8, destination, tag, &
        mpi_comm_world, request, ierror)
      eastMPI => eastMPI%next
    END DO

    ! West face MPI sends
    DO
      IF (.NOT. ASSOCIATED(westMPI)) EXIT
      b1 => Blocks(westMPI%ID)
      DO j = 1, jBlockSize
        p1 => b1%Points(2, j)%T
        j_buffer(j) = p1
      END DO
      destination = b1%westFace%neighborProc
      tag = wB
      CALL MPI_Isend(j_buffer, jBlockSize, MPI_REAL8, destination, tag, &
        mpi_comm_world, request, ierror)
      westMPI => westMPI%next
    END DO

    ! North East corner MPI sends
    DO
      IF (.NOT. ASSOCIATED(neMPI)) EXIT
      b1 => Blocks(neMPI%ID)
      p1 => b1%Points(iBlockSize-1,jBlockSize-1)%T
      buffer = p1
      destination = b1%NECorner%neighborProc
      tag = nB + eB * 10
      CALL MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, &
        request, ierror)
      neMPI => neMPI%next
    END DO

    ! South East corner MPI sends
    DO
      IF (.NOT. ASSOCIATED(seMPI)) EXIT
      b1 => Blocks(seMPI%ID)
      p1 => b1%Points(iBlockSize-1,2)%T
      buffer = p1
      destination = b1%SECorner%neighborProc
      tag = sB + eB * 10
      CALL MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, &
        request, ierror)
      seMPI => seMPI%next
    END DO

    ! South West corner MPI sends
    DO
      IF (.NOT. ASSOCIATED(swMPI)) EXIT
      b1 => Blocks(swMPI%ID)
      p1 => b1%Points(2,2)%T
      buffer = p1
      destination = b1%SWCorner%neighborProc
      tag = sB + wB * 10
      CALL MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, &
        request, ierror)
      swMPI => swMPI%next
    END DO

    ! North West corner MPI sends
    DO
      IF (.NOT. ASSOCIATED(nwMPI)) EXIT
      b1 => Blocks(nwMPI%ID)
      p1 => b1%Points(2,jBlockSize-1)%T
      buffer = p1
      destination = b1%nwCorner%neighborProc
      tag = nB + wB * 10
      CALL MPI_Isend(buffer, 1, MPI_REAL8, destination, tag, mpi_comm_world, &
        request, ierror)
      nwMPI => nwMPI%next
    END DO
  END SUBROUTINE

  ! Our neighbor block is on different proc, we must receive with MPI.
  SUBROUTINE mpi_receives
    TYPE (BlockType), POINTER :: b1
    REAL(KIND=8), POINTER :: p1
    INTEGER :: i, j, tag, source
    REAL(KIND = 8) :: i_buffer(iBlockSize), j_buffer(jBlockSize)
    REAL(KIND = 8) :: buffer

    ! Set pointers to first elements in linked lists
    southMPI => southMPIList
    northMPI => northMPIList
    westMPI  => westMPIList
    eastMPI  => eastMPIList
    swMPI    => swMPIList
    nwMPI    => nwMPIList
    neMPI    => neMPIList
    seMPI    => seMPIList

    ! South face MPI receives
    DO
      IF (.NOT. ASSOCIATED(southMPI)) EXIT
      b1 => Blocks(southMPI%ID)
      source = b1%southFace%neighborProc
      tag = nB
      CALL MPI_RECV(i_buffer, iBlockSize, MPI_REAL8, source, tag, mpi_comm_world, &
        STATUS, ierror)
      DO i = 1, iBlockSize
        p1 => b1%Points(i, 0)%T
        p1 = i_buffer(i)
      END DO
      southMPI => southMPI%next
    END DO

    ! North face MPI receives
    DO
      IF (.NOT. ASSOCIATED(northMPI)) EXIT
      b1 => Blocks(northMPI%ID)
      source = b1%northFace%neighborProc
      tag = sB
      CALL MPI_RECV(i_buffer, iBlockSize, MPI_REAL8, source, tag, mpi_comm_world, &
        STATUS, ierror)
      DO i = 1, iBlockSize
        p1 => b1%Points(i, jBlockSize+1)%T
        p1 = i_buffer(i)
      END DO
      northMPI => northMPI%next
    END DO

    ! West face MPI receives
    DO
      IF (.NOT. ASSOCIATED(westMPI)) EXIT
      b1 => Blocks(westMPI%ID)
      source = b1%westFace%neighborProc
      tag = eB
      CALL MPI_RECV(j_buffer, jBlockSize, MPI_REAL8, source, tag, mpi_comm_world, &
        STATUS, ierror)
      DO j = 1, jBlockSize
        p1 => b1%Points(0, j)%T
        p1 = j_buffer(j)
      END DO
      westMPI => westMPI%next
    END DO

    ! East face MPI receives
    DO
      IF (.NOT. ASSOCIATED(eastMPI)) EXIT
      b1 => Blocks(eastMPI%ID)
      source = b1%eastFace%neighborProc
      tag = wB
      CALL MPI_RECV(j_buffer, jBlockSize, MPI_REAL8, source, tag, mpi_comm_world, &
        STATUS, ierror)
      DO j = 1, jBlockSize
        p1 => b1%Points(iBlockSize+1, j)%T
        p1 = j_buffer(j)
      END DO
      eastMPI => eastMPI%next
    END DO

    ! South West corner MPI sends
    DO
      IF (.NOT. ASSOCIATED(swMPI)) EXIT
      b1 => Blocks(swMPI%ID)
      source = b1%SWCorner%neighborProc
      tag = nB + eB * 10
      CALL MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, STATUS, ierror)
      p1 => b1%Points(0,0)%T
      p1 = buffer
      swMPI => swMPI%next
    END DO

    ! North West corner MPI sends
    DO
      IF (.NOT. ASSOCIATED(nwMPI)) EXIT
      b1 => Blocks(nwMPI%ID)
      source = b1%NWCorner%neighborProc
      tag = sB + eB * 10
      CALL MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, STATUS, ierror)
      p1 => b1%Points(0,jBlockSize+1)%T
      p1 = buffer
      nwMPI => nwMPI%next
    END DO

    ! North East corner MPI sends
    DO
      IF (.NOT. ASSOCIATED(neMPI)) EXIT
      b1 => Blocks(neMPI%ID)
      source = b1%NECorner%neighborProc
      tag = sB + wB * 10
      CALL MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, STATUS, ierror)
      p1 => b1%Points(iBlockSize+1,jBlockSize+1)%T
      p1 = buffer
      neMPI => neMPI%next
    END DO

    ! South East corner MPI sends
    DO
      IF (.NOT. ASSOCIATED(seMPI)) EXIT
      b1 => Blocks(seMPI%ID)
      source = b1%SECorner%neighborProc
      tag = nB + wB * 10
      CALL MPI_RECV(buffer, 1, MPI_REAL8, source, tag, mpi_comm_world, STATUS, ierror)
      p1 => b1%Points(iBlockSize+1,0)%T
      p1 = buffer
      seMPI => seMPI%next
    END DO
  END SUBROUTINE

  SUBROUTINE initialize_linked_lists
    TYPE (BlockType), POINTER :: b
    INTEGER :: n_

    ! Temp lists
    TYPE(LinkedList), POINTER :: northLocalTemp, southLocalTemp, eastLocalTemp, westLocalTemp
    TYPE(LinkedList), POINTER :: northMPITemp, southMPITemp, eastMPITemp, westMPITemp
    TYPE(LinkedList), POINTER :: neLocalTemp, nwLocalTemp, swLocalTemp, seLocalTemp
    TYPE(LinkedList), POINTER :: neMPITemp, nwMPITemp, swMPITemp, seMPITemp

    ! For each block...
    DO n_ = 1, MyNBlocks
      b => Blocks(n_)

      ! North Face
      ! If this is an internal boundary then we put the
      ! block into the appropriate local list.
      IF (b%northFace%BC == INTERNAL_BOUNDARY) THEN

        ! If this list has not been allocated then we must
        ! first allocate it.
        IF (.NOT. ASSOCIATED(northLocalList)) THEN
          ALLOCATE(northLocalList)
          northLocalTemp => northLocalList
          NULLIFY(northLocalTemp%next)
          northLocalTemp%ID = n_
        ELSE
          ! If the list is already allocated then we can
          ! allocate the next link in the linked list and
          ! assign the block to it.
          ALLOCATE(northLocalTemp%next)
          northLocalTemp => northLocalTemp%next
          NULLIFY(northLocalTemp%next)
          northLocalTemp%ID = n_
        END IF

      ! Otherwise if this is a interproc boundary
      ! then we do a similar thing with the appropriate
      ! MPI list.
      ELSE IF (b%northFace%BC == PROC_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(northMPIList)) THEN
          ALLOCATE(northMPIList)
          northMPITemp => northMPIList
          NULLIFY(northMPITemp%next)
          northMPITemp%ID = n_
        ELSE
          ALLOCATE(northMPITemp%next)
          northMPITemp => northMPITemp%next
          NULLIFY(northMPITemp%next)
          northMPITemp%ID = n_
        END IF
      END IF

      ! South Face
      IF (b%southFace%BC == INTERNAL_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(southLocalList)) THEN
          ALLOCATE(southLocalList)
          southLocalTemp => southLocalList
          NULLIFY(southLocalTemp%next)
          southLocalTemp%ID = n_
        ELSE
          ALLOCATE(southLocalTemp%next)
          southLocalTemp => southLocalTemp%next
          NULLIFY(southLocalTemp%next)
          southLocalTemp%ID = n_
        END IF
      ELSE IF (b%southFace%BC == PROC_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(southMPIList)) THEN
          ALLOCATE(southMPIList)
          southMPITemp => southMPIList
          NULLIFY(southMPITemp%next)
          southMPITemp%ID = n_
        ELSE
          ALLOCATE(southMPITemp%next)
          southMPITemp => southMPITemp%next
          NULLIFY(southMPITemp%next)
          southMPITemp%ID = n_
        END IF
      END IF

      ! East Face
      IF (b%eastFace%BC == INTERNAL_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(eastLocalList)) THEN
          ALLOCATE(eastLocalList)
          eastLocalTemp => eastLocalList
          NULLIFY(eastLocalTemp%next)
          eastLocalTemp%ID = n_
        ELSE
          ALLOCATE(eastLocalTemp%next)
          eastLocalTemp => eastLocalTemp%next
          NULLIFY(eastLocalTemp%next)
          eastLocalTemp%ID = n_
        END IF
      ELSE IF (b%eastFace%BC == PROC_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(eastMPIList)) THEN
          ALLOCATE(eastMPIList)
          eastMPITemp => eastMPIList
          NULLIFY(eastMPITemp%next)
          eastMPITemp%ID = n_
        ELSE
          ALLOCATE(eastMPITemp%next)
          eastMPITemp => eastMPITemp%next
          NULLIFY(eastMPITemp%next)
          eastMPITemp%ID = n_
        END IF
      END IF

      ! West Face
      IF (b%westFace%BC == INTERNAL_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(westLocalList)) THEN
          ALLOCATE(westLocalList)
          westLocalTemp => westLocalList
          NULLIFY(westLocalTemp%next)
          westLocalTemp%ID = n_
        ELSE
          ALLOCATE(westLocalTemp%next)
          westLocalTemp => westLocalTemp%next
          NULLIFY(westLocalTemp%next)
          westLocalTemp%ID = n_
        END IF
      ELSE IF (b%westFace%BC == PROC_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(westMPIList)) THEN
          ALLOCATE(westMPIList)
          westMPITemp => westMPIList
          NULLIFY(westMPITemp%next)
          westMPITemp%ID = n_
        ELSE
          ALLOCATE(westMPITemp%next)
          westMPITemp => westMPITemp%next
          NULLIFY(westMPITemp%next)
          westMPITemp%ID = n_
        END IF
      END IF

      ! NE Corner
      IF (b%NECorner%BC == INTERNAL_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(neLocalList)) THEN
          ALLOCATE(neLocalList)
          neLocalTemp => neLocalList
          NULLIFY(neLocalTemp%next)
          neLocalTemp%ID = n_
        ELSE
          ALLOCATE(neLocalTemp%next)
          neLocalTemp => neLocalTemp%next
          NULLIFY(neLocalTemp%next)
          neLocalTemp%ID = n_
        END IF
      ELSE IF (b%NECorner%BC == PROC_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(neMPIList)) THEN
          ALLOCATE(neMPIList)
          neMPITemp => neMPIList
          NULLIFY(neMPITemp%next)
          neMPITemp%ID = n_
        ELSE
          ALLOCATE(neMPITemp%next)
          neMPITemp => neMPITemp%next
          NULLIFY(neMPITemp%next)
          neMPITemp%ID = n_
        END IF
      END IF

      ! NW Corner
      IF (b%NWCorner%BC == INTERNAL_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(nwLocalList)) THEN
          ALLOCATE(nwLocalList)
          nwLocalTemp => nwLocalList
          NULLIFY(nwLocalTemp%next)
          nwLocalTemp%ID = n_
        ELSE
          ALLOCATE(nwLocalTemp%next)
          nwLocalTemp => nwLocalTemp%next
          NULLIFY(nwLocalTemp%next)
          nwLocalTemp%ID = n_
        END IF
      ELSE IF (b%NWCorner%BC == PROC_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(nwMPIList)) THEN
          ALLOCATE(nwMPIList)
          nwMPITemp => nwMPIList
          NULLIFY(nwMPITemp%next)
          nwMPITemp%ID = n_
        ELSE
          ALLOCATE(nwMPITemp%next)
          nwMPITemp => nwMPITemp%next
          NULLIFY(nwMPITemp%next)
          nwMPITemp%ID = n_
        END IF
      END IF

      ! SW Corner
      IF (b%swCorner%BC == INTERNAL_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(swLocalList)) THEN
          ALLOCATE(swLocalList)
          swLocalTemp => swLocalList
          NULLIFY(swLocalTemp%next)
          swLocalTemp%ID = n_
        ELSE
          ALLOCATE(swLocalTemp%next)
          swLocalTemp => swLocalTemp%next
          NULLIFY(swLocalTemp%next)
          swLocalTemp%ID = n_
        END IF
      ELSE IF (b%swCorner%BC == PROC_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(swMPIList)) THEN
          ALLOCATE(swMPIList)
          swMPITemp => swMPIList
          NULLIFY(swMPITemp%next)
          swMPITemp%ID = n_
        ELSE
          ALLOCATE(swMPITemp%next)
          swMPITemp => swMPITemp%next
          NULLIFY(swMPITemp%next)
          swMPITemp%ID = n_
        END IF
      END IF

      ! SE Corner
      IF (b%seCorner%BC == INTERNAL_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(seLocalList)) THEN
          ALLOCATE(seLocalList)
          seLocalTemp => seLocalList
          NULLIFY(seLocalTemp%next)
          seLocalTemp%ID = n_
        ELSE
          ALLOCATE(seLocalTemp%next)
          seLocalTemp => seLocalTemp%next
          NULLIFY(seLocalTemp%next)
          seLocalTemp%ID = n_
        END IF
      ELSE IF (b%seCorner%BC == PROC_BOUNDARY) THEN
        IF (.NOT. ASSOCIATED(seMPIList)) THEN
          ALLOCATE(seMPIList)
          seMPITemp => seMPIList
          NULLIFY(seMPITemp%next)
          seMPITemp%ID = n_
        ELSE
          ALLOCATE(seMPITemp%next)
          seMPITemp => seMPITemp%next
          NULLIFY(seMPITemp%next)
          seMPITemp%ID = n_
        END IF
      END IF
    END DO
  END SUBROUTINE

END MODULE MainRoutines
