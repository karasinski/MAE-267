! This module contains all the subroutines used to
! create our blocks and initial grid/temp files.
module GridCreation
  use BlockModule
  implicit none

contains
  subroutine create_blocks(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:)
    type (BlockType), pointer :: b
    integer :: n_, proc = 1, iN, iM
    integer :: nBound, sBound, eBound, wBound
    integer :: highI, highJ

    !
    !            |             |
    !            |    North    |
    !          NW|   (n + N)   |NE
    ! (n + N + 1)|             |(n + N + 1)
    ! -------------------------------------
    !            |             |
    !     West   |   Current   |    East
    !   (n - 1)  |     (n)     |  (n + 1)
    !            |             |
    ! -------------------------------------
    !          SW|             |SE
    ! (n - N - 1)|    South    |(n - N + 1)
    !            |   (n - N)   |
    !            |             |
    !

    ! Loop over our M x N blocks and a pack an
    ! array of BlocksCollection(n_).
    n_ = 1
    do iM = 1, M
      do iN = 1, N
        b => BlocksCollection(n_)

        ! Global id.
        b%id = n_

        ! We are only using one processor, so we must
        ! be on that processor.
        b%proc = proc

        ! High/low i/j corresponds to the
        ! global block number start and stop.
        highI = 1 + iN * (iBlockSize - 1)
        b%lowI = highI - (iBlockSize - 1)
        highJ = 1 + iM * (jBlockSize - 1)
        b%lowJ = highJ - (jBlockSize - 1)

        ! Assume all faces are internal and find out who
        ! their neighbor is. Set the neighbor block and set
        ! their processor to 1 (only one processor).
        b%northFace%neighborBlock = n_ + N
        b%northFace%neighborProc = proc
        b%southFace%neighborBlock = n_ - N
        b%southFace%neighborProc = proc
        b%eastFace%neighborBlock = n_ + 1
        b%eastFace%neighborProc = proc
        b%westFace%neighborBlock = n_ - 1
        b%westFace%neighborProc = proc

        ! Assume all corners are internal (BC=-1) and set their neighbor.
        ! Set their processor to 1 (only one processor).
        b%NECorner%BC = -1
        b%NECorner%neighborBlock = n_ + N + 1
        b%NECorner%neighborProc = proc
        b%SECorner%BC = -1
        b%SECorner%neighborBlock = n_ - N + 1
        b%SECorner%neighborProc = proc
        b%SWCorner%BC = -1
        b%SWCorner%neighborBlock = n_ - N - 1
        b%SWCorner%neighborProc = proc
        b%NWCorner%BC = -1
        b%NWCorner%neighborBlock = n_ + N - 1
        b%NWCorner%neighborProc = proc

        ! Determine whether block is on boundary or
        ! if they is internal for each face.
        if (highJ == JMAX) then
          nBound = nB
        else
          nBound = -1
        end if
        if (highI == IMAX) then
          eBound = eB
        else
          eBound = -1
        end if
        if (b%lowJ == 1) then
          sBound = sB
        else
          sBound = -1
        end if
        if (b%lowI == 1) then
          wBound = wB
        else
          wBound = -1
        end if

        ! Bounds and Faces
        ! Set our boundary conditions for each face.
        b%northFace%BC = nBound
        b%southFace%BC = sBound
        b%eastFace%BC = eBound
        b%westFace%BC = wBound

        ! Check if faces are on the boundary and
        ! if so, set them to point to nothing.
        if (b%northFace%BC /= -1) then
          b%northFace%neighborBlock = 0
          b%northFace%neighborProc = 0
        end if

        if (b%southFace%BC /= -1) then
          b%southFace%neighborBlock = 0
          b%southFace%neighborProc = 0
        end if

        if (b%eastFace%BC /= -1) then
          b%eastFace%neighborBlock = 0
          b%eastFace%neighborProc = 0
        end if

        if (b%westFace%BC /= -1) then
          b%westFace%neighborBlock = 0
          b%westFace%neighborProc = 0
        end if

        ! If corners are on boundary then they don't
        ! have associated blocks. Set their processor
        ! to 0 so we know that there's nowhere to look.
        ! North East Corner
        if (b%northFace%BC == nB) then
          b%NECorner%BC = nBound
          b%NECorner%neighborBlock = 0
          b%NECorner%neighborProc = 0
        else if (b%eastFace%BC == eB) then
          b%NECorner%BC = eBound
          b%NECorner%neighborBlock = 0
          b%NECorner%neighborProc = 0
        end if

        ! South East Corner
        if (b%southFace%BC == sB) then
          b%SECorner%BC = sBound
          b%SECorner%neighborBlock = 0
          b%SECorner%neighborProc = 0
        else if (b%eastFace%BC == eB) then
          b%SECorner%BC = eBound
          b%SECorner%neighborBlock = 0
          b%SECorner%neighborProc = 0
        end if

        ! South West Corner
        if (b%southFace%BC == sB) then
          b%SWCorner%BC = sBound
          b%SWCorner%neighborBlock = 0
          b%SWCorner%neighborProc = 0
        else if (b%westFace%BC == wB) then
          b%SWCorner%BC = wBound
          b%SWCorner%neighborBlock = 0
          b%SWCorner%neighborProc = 0
        end if

        ! North West Corner
        if (b%northFace%BC == nB) then
          b%NWCorner%BC = nBound
          b%NWCorner%neighborBlock = 0
          b%NWCorner%neighborProc = 0
        else if (b%westFace%BC == wB) then
          b%NWCorner%BC = wBound
          b%NWCorner%neighborBlock = 0
          b%NWCorner%neighborProc = 0
        end if

        ! Calculate 'true' bounds of blocks.
        ! North face.
        if (b%northFace%BC == -1) then
          b%localJMAX = jBlockSize
        else
          b%localJMAX = jBlockSize - 1
        end if

        ! South face.
        if (b%southFace%BC == -1) then
          b%localJMIN = 0
        else
          b%localJMIN = 1
        end if

        ! East face.
        if (b%eastFace%BC == -1) then
          b%localIMAX = iBlockSize
        else
          b%localIMAX = iBlockSize - 1
        end if

        ! West face.
        if (b%westFace%BC == -1) then
          b%localIMIN = 0
        else
          b%localIMIN = 1
        end if

        ! Weight of each block, used for load balancing.
        ! We add the 'area' of the block to the length of the sides of the block
        ! (communication cost), which we weigh equal to the solver weight.
        b%comm = 0
        if (b%northFace%BC == -1) b%comm = b%comm + b%localJMAX - b%localJMIN
        if (b%southFace%BC == -1) b%comm = b%comm + b%localJMAX - b%localJMIN
        if (b%eastFace%BC == -1) b%comm = b%comm + b%localIMAX - b%localIMIN
        if (b%westFace%BC == -1) b%comm = b%comm + b%localIMAX - b%localIMIN
        b%size = (b%localIMAX - b%localIMIN) * (b%localJMAX - b%localJMIN) + b%comm

        n_ = n_ + 1
      end do
    end do
  end subroutine

  ! Intialize our locations and prime locations.
  subroutine initialize_block_grid(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:)
    type (GridPoint), pointer :: p
    type (BlockType), pointer :: b
    integer :: i, j, n_

    do n_ = 1, nBlocks
      do j = 0, jBlockSize+1
        do i = 0, iBlockSize+1
          b => BlocksCollection(n_)
          p => BlocksCollection(n_)%Points(i, j)

          ! Have to convert from local i/j to global i/j.
          p%xp = cos( 0.5d0 * pi * dfloat(IMAX - (b%lowI + (i - 1 ))) / dfloat(IMAX-1))
          p%yp = cos( 0.5d0 * pi * dfloat(JMAX - (b%lowJ + (j - 1 ))) / dfloat(JMAX-1))

          p%x = p%xp * cos( rot ) + ( 1.d0 - p%yp ) * sin( rot )
          p%y = p%yp * cos( rot ) + p%xp * sin( rot )
        end do
      end do
    end do
  end subroutine

  ! Dirichlet conditions.
  subroutine initialize_block_temp(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:)
    type (BlockType), pointer :: b
    type (GridPoint), pointer :: p(:,:)
    integer :: i, j, n_
    real(kind=8) :: T_0 = 3.5d0

    do n_ = 1, nBlocks
      b => BlocksCollection(n_)
      p => BlocksCollection(n_)%Points

      ! Initialize all points to T_0.
      p%T = T_0

      ! Set the boundary conditions if our blocks are on the boundary.
      if (b%northFace%BC == nB) then
        do i = 1, iBlockSize
          p(i, jBlockSize)%T = 5.d0 * (sin(pi * p(i, jBlockSize)%xp) + 1.d0)
        end do
      end if

      if (b%southFace%BC == sB) then
        do i = 1, iBlockSize
          p(i, 1)%T = abs(cos(pi * p(i,1)%xp)) + 1.d0
        end do
      end if

      if (b%eastFace%BC == eB) then
        do j = 1, jBlockSize
          p(iBlockSize, j)%T = (3.d0 * p(iBlockSize, j)%yp) + 2.d0
        end do
      end if

      if (b%westFace%BC == wB) then
        do j = 1, jBlockSize
          p(1, j)%T = (3.d0 * p(1, j)%yp) + 2.d0
        end do
      end if

    end do
  end subroutine

  recursive subroutine distribute_blocks(BlocksCollection, Procs, method, fudge_factor)
    type (BlockType), target :: BlocksCollection(:)
    type (Proc), allocatable :: Procs(:)
    type (BlockType), pointer :: b
    real(kind=8) :: optimal, rel, fudge_factor
    integer :: method, neighbors, p_, n_, b_, m_, i, sum = 0
    integer :: most_empty_proc = 0, largest_block = 0, largest_block_number = 0
    logical :: unAssignedBlocks = .TRUE., error = .FALSE.

    ! Set starting weights on each proc equal to zero.
    Procs%weight = 0
    Procs%nBlocks = 0

    ! Find the sum of the weights of the blocks.
    ! This does not take into account any interblock communcation,
    ! but is simply a measure of how long the solver should take
    ! per iteration. A better algorithm would attempt to place blocks
    ! near each other in order to reduce (eliminate) communcation cost.

    do n_ = 1, nBlocks
      b => BlocksCollection(n_)
      sum = sum + b%size
    end do
    optimal = dfloat(sum)/dfloat(nProcs)

    ! Log each step for debuggin.
    !      write(*,*), '------------------'
    !      do i = 1, nProcs
    !        write(*, *), 'proc ', i, Procs(i)%weight
    !      end do


    if (method == 1) then
      ! While there are unassigned blocks, distribute
      ! blocks to most empty processor. This method does not
      ! attempt to put nearby blocks on the same processor.
      do while (unAssignedBlocks)
        ! Reset largest block to 0.
        largest_block = 0
        do n_ = 1, nBlocks
          b => BlocksCollection(n_)

          ! Find the block with the most weight.
          if (b%size > largest_block) then
            largest_block = b%size
            largest_block_number = n_
          end if
        end do

        ! Find the proc with the least weight.
        most_empty_proc = minloc(Procs%weight,1)

        ! Add this block to this proc.
        ! Increase number of blocks on this proc by 1.
        Procs(most_empty_proc)%nBlocks = Procs(most_empty_proc)%nBlocks + 1

        ! Increase weight on the proc by the added block's weight.
        Procs(most_empty_proc)%weight = Procs(most_empty_proc)%weight + largest_block

        ! Copy block over to proc.
        Procs(most_empty_proc)%Blocks(Procs(most_empty_proc)%nBlocks) = &
          BlocksCollection(largest_block_number)

        ! Reduce the size of the block in the blockscollection to zero so it is
        ! not assigned to any other procs.
        BlocksCollection(largest_block_number)%size = 0

        ! Check for unassigned blocks.
        unAssignedBlocks = .FALSE.
        do n_ = 1, nBlocks
          b => BlocksCollection(n_)

          ! If a block hasn't been assigned, we need to
          ! continue looping.
          if (b%size > 0) unAssignedBlocks = .TRUE.
        end do
      end do
    else if (method == 2) then
      ! While this processor has less than the optimal weight,
      ! if the largest weight won't push it past optimal,
      ! add the next block's weight.
      !    do while (unAssignedBlocks)

      ! Check to see if proc has less than optimal weight.
      do p_ = 1, nProcs
        do n_ = 1, nBlocks
          b => BlocksCollection(n_)

          ! Find the first nonzero size block.
          if (b%size > 0) then
            largest_block = b%size
            largest_block_number = n_
          end if
        end do

        do while (Procs(p_)%weight + largest_block < 1.05 * optimal)
          ! Reset largest block to zero.
          largest_block = 0

          ! Loop over all blocks.
          do n_ = 1, nBlocks
            b => BlocksCollection(n_)

            ! Find the first nonzero size block.
            if (b%size > 0) then
              largest_block = b%size
              largest_block_number = n_
            end if
          end do

          if (largest_block == 0) then
            write(*,*) "Largest Block is 0, exiting"
            exit
          end if

          ! Add this block to this proc.
          ! Increase number of blocks on this proc by 1.
          Procs(p_)%nBlocks = Procs(p_)%nBlocks + 1

          ! Increase weight on the proc by the added block's weight.
          Procs(p_)%weight = Procs(p_)%weight + largest_block

          ! Copy block over to proc.
          Procs(p_)%Blocks(Procs(p_)%nBlocks) = &
            BlocksCollection(largest_block_number)

          ! Reduce the size of the block in the blockscollection to zero so it is
          ! not assigned to any other procs.
          BlocksCollection(largest_block_number)%size = 0
        end do
      end do
    else if (method == 3) then
      ! Vary this number to find an optimal setting.

      ! Check to see if proc has less than optimal weight.
      do p_ = 1, nProcs

        do n_ = 1, nBlocks
          b => BlocksCollection(n_)

          ! Find the first nonzero size block.
          if (b%size > 0) then
            largest_block = b%size
            largest_block_number = n_
          end if
        end do

        do while (Procs(p_)%weight + largest_block < fudge_factor * optimal)
          ! Reset largest block to zero.
          largest_block = 0

          ! Loop over all blocks.
          do n_ = 1, nBlocks
            b => BlocksCollection(n_)

            ! Find the first nonzero size block.
            if (b%size > 0) then
              largest_block = b%size
              largest_block_number = n_
            end if
          end do

          if (largest_block == 0) then
            write(*,*) "Largest Block is 0, we're done here"
            exit
          end if

          ! Add this block to this proc.
          ! Increase number of blocks on this proc by 1.
          Procs(p_)%nBlocks = Procs(p_)%nBlocks + 1

          ! Increase weight on the proc by the added block's weight.
          Procs(p_)%weight = Procs(p_)%weight + largest_block

          ! Copy block over to proc.
          Procs(p_)%Blocks(Procs(p_)%nBlocks) = &
            BlocksCollection(largest_block_number)

          ! Reduce the size of the block in the blockscollection to zero so it is
          ! not assigned to any other procs.
          BlocksCollection(largest_block_number)%size = 0
        end do
      end do

    else if (method == 4) then
      ! While there are unassigned blocks, distribute
      ! blocks to most empty processor. This method does not
      ! attempt to put nearby blocks on the same processor.

      ! Reset largest block to 0.
      largest_block = 0
      i = 1
      do n_ = 1, N
        do m_ = 1, M
          i = m_ + (n_ - 1) * M
          write(*,*), n_, N, m_, M, i
          b => BlocksCollection(i)

          ! Find the block with the most weight.
          if (b%size > largest_block) then
            largest_block = b%size
          end if


          ! Find the proc with the least weight.
          most_empty_proc = minloc(Procs%weight,1)

          ! Add this block to this proc.
          ! Increase number of blocks on this proc by 1.
          Procs(most_empty_proc)%nBlocks = Procs(most_empty_proc)%nBlocks + 1

          ! Increase weight on the proc by the added block's weight.
          Procs(most_empty_proc)%weight = Procs(most_empty_proc)%weight + largest_block

          ! Copy block over to proc.
          Procs(most_empty_proc)%Blocks(Procs(most_empty_proc)%nBlocks) = &
            BlocksCollection(i)

          ! Reduce the size of the block in the blockscollection to zero so it is
          ! not assigned to any other procs.
          BlocksCollection(i)%size = 0

        end do
      end do

    end if

    ! If we screwed this up we need to terminated execution.
    do n_ = 1, nBlocks
      if (BlocksCollection(n_)%size > 0) then
        write(*,*), "We forgot block ", n_, BlocksCollection(n_)%size
        !        error = .TRUE.
      end if
    end do

    if (error) then
      write(*,*), "Sorry, something went terribly wrong."
      write(*,*), "Program exiting."
      STOP
    end if

    ! Write the total weight, number of procs, ideal weight per proc.
    write(*,*)
    write(*,*), '      Weight  ', 'Ideal Weight/proc'
    write(*,*), sum, optimal
    write(*,*)

    ! Write out the weight on each proc.
    ! Find the relative error.
    rel = 0.d0
    write(*, *), '          proc #      ', "nBlocks  ", "Weight ", "Err"
    do n_ = 1, nProcs
      write(*, *), n_, Procs(n_)%nBlocks, Procs(n_)%weight, dfloat(Procs(n_)%weight) - optimal
      rel = rel + abs( optimal - dfloat(Procs(n_)%weight) )
    end do

    rel = rel / optimal * 100.d0
    write (*,*), "The relative error is ", rel, " percent."

  end subroutine
end module
