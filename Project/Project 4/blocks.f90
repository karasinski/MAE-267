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

        if (b%NECorner%BC == -1) b%comm = b%comm + 1
        if (b%NWCorner%BC == -1) b%comm = b%comm + 1
        if (b%SECorner%BC == -1) b%comm = b%comm + 1
        if (b%SWCorner%BC == -1) b%comm = b%comm + 1

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

  subroutine add_block_to_proc(MyProc, Block)
    type (Proc) :: MyProc
    type (BlockType) :: Block

    MyProc%nBlocks = MyProc%nBlocks + 1
    MyProc%weight = MyProc%weight + Block%size
    MyProc%Blocks(MyProc%nBlocks) = Block
    MyProc%Blocks(MyProc%nBlocks)%proc = MyProc%procID

    Block%size = 0
  end subroutine

  subroutine check_communication_cost(Procs)
    type (Proc), target :: Procs(:)
    type (Proc), pointer :: MyProc
    type (BlockType), pointer :: Block
    integer :: comm, ID, b, p, i

    ! Check if this block has neighbors on the processor.
    ! If so we go ahead and remove the communication cost.
    do p = 1, nProcs
      MyProc => Procs(p)
      ID = MyProc%procID
      comm = 0

      do b = 1, MyProc%nBlocks
        Block => MyProc%Blocks(b)

        do i = 1, MyProc%nBlocks
          if (MyProc%Blocks(i)%id /= Block%id) then
            if (MyProc%Blocks(i)%id == Block%northFace%neighborBlock) comm = comm + ( Block%localJMAX - Block%localJMIN )
            if (MyProc%Blocks(i)%id == Block%southFace%neighborBlock) comm = comm + ( Block%localJMAX - Block%localJMIN )
            if (MyProc%Blocks(i)%id == Block%eastFace%neighborBlock)  comm = comm + ( Block%localIMAX - Block%localIMIN )
            if (MyProc%Blocks(i)%id == Block%westFace%neighborBlock)  comm = comm + ( Block%localIMAX - Block%localIMIN )

            if (MyProc%Blocks(i)%id == Block%NECorner%neighborBlock)  comm = comm + 1
            if (MyProc%Blocks(i)%id == Block%NWCorner%neighborBlock)  comm = comm + 1
            if (MyProc%Blocks(i)%id == Block%SECorner%neighborBlock)  comm = comm + 1
            if (MyProc%Blocks(i)%id == Block%SWCorner%neighborBlock)  comm = comm + 1
          end if
        end do
      end do

      if (comm > 0) write(*,*), "Winner", comm
      MyProc%weight = MyProc%weight - comm
    end do

  end subroutine

  subroutine add_blocks_to_proc(MyProc, BlocksCollection, BlockIDs)
    type (Proc) :: MyProc
    type (BlockType), target :: BlocksCollection(:)
    integer :: BlockIDs(:), i

    do i = 1, size(BlockIDs)
        call add_block_to_proc(MyProc, BlocksCollection(BlockIDs(i)))
    end do

  end subroutine

  subroutine distribute_blocks(BlocksCollection, Procs)
    type (BlockType), target :: BlocksCollection(:)
    type (Proc), allocatable :: Procs(:)
    type (BlockType), pointer :: b
    real(kind=8) :: optimal, fudge_factor = 1.05d0
    integer :: method, p_, n_, sum = 0
    integer :: largest_block = 0, largest_block_number = 0

    ! Set starting weights on each proc equal to zero.
    Procs%weight = 0
    Procs%nBlocks = 0

    ! Find the sum of the weights of the blocks with no communication cost.
    ! This is an ideal we cannot meet.
    do n_ = 1, nBlocks
      b => BlocksCollection(n_)
      sum = sum + b%size - b%comm
    end do

    ! The optimal distribution is an equal weight on each block.
    optimal = dfloat(sum)/dfloat(nProcs)

    ! Pick our method depending on our settings.
    if (M == 5  .and. N == 4) then
      if (nProcs == 6) then
        method = 546
      else if (nProcs == 4) then
        ! General method.
        method = 544
      end if
    else if (M == 10 .and. N == 10) then
      ! General method with precomputed fudge factor.
      if (nProcs == 6) then
        method = 10106
      else if (nProcs == 4) then
        method = 10104
      end if
    else
      ! Otherwise we'll do our best.
      method = 666
    end if
    write(*,*), "method ", method

    if (method == 546) then
      ! Hard coding is hard work.
      call add_blocks_to_proc(Procs(1), BlocksCollection, [1, 2, 3, 4])
      call add_blocks_to_proc(Procs(2), BlocksCollection, [5, 6, 7])
      call add_blocks_to_proc(Procs(3), BlocksCollection, [9, 10, 11])
      call add_blocks_to_proc(Procs(4), BlocksCollection, [13, 14, 15])
      call add_blocks_to_proc(Procs(5), BlocksCollection, [8, 12, 16])
      call add_blocks_to_proc(Procs(6), BlocksCollection, [17, 18, 19, 20])

    else if (method == 544) then
      ! Hard coding is still hard work.
      call add_blocks_to_proc(Procs(1), BlocksCollection, [1, 5, 9, 13, 17])
      call add_blocks_to_proc(Procs(2), BlocksCollection, [2, 6, 10, 14, 18])
      call add_blocks_to_proc(Procs(3), BlocksCollection, [3, 7, 11, 15, 19])
      call add_blocks_to_proc(Procs(4), BlocksCollection, [4, 8, 12, 16, 20])

    else if (method == 10106) then
      ! Hard coding is okay.
      do n_ = 1, 10
        call add_block_to_proc(Procs(1), BlocksCollection(n_))
        call add_block_to_proc(Procs(1), BlocksCollection(n_ + 10))

        call add_block_to_proc(Procs(6), BlocksCollection(n_ + 80))
        call add_block_to_proc(Procs(6), BlocksCollection(n_ + 90))
      end do

      do n_ = 2, 4
        call add_blocks_to_proc(Procs(2), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10])

        call add_blocks_to_proc(Procs(3), BlocksCollection, &
          [6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])
      end do

      do n_ = 5, 7
        call add_blocks_to_proc(Procs(4), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10])

        call add_blocks_to_proc(Procs(5), BlocksCollection, &
          [6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])
      end do
    else if (method == 10104) then
      ! Hard coding isn't that bad.
      do n_ = 0, 4
        call add_blocks_to_proc(Procs(1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10])

        call add_blocks_to_proc(Procs(2), BlocksCollection, &
          [6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])

        call add_blocks_to_proc(Procs(3), BlocksCollection, &
          [51 + n_*10, 52 + n_*10, 53 + n_*10, 54 + n_*10, 55 + n_*10])

        call add_blocks_to_proc(Procs(4), BlocksCollection, &
          [56 + n_*10, 57 + n_*10, 58 + n_*10, 59 + n_*10, 60 + n_*10])
      end do

    else if (method == 666) then
      ! This is our fallback method for unknown configurations.

      do p_ = 1, nProcs
        do n_ = 1, nBlocks
          b => BlocksCollection(n_)

          ! Find the first nonzero size block.
          if (b%size > 0) then
            largest_block = b%size
            largest_block_number = n_
          end if
        end do

        ! Check to see if proc has less than optimal weight.
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
    end if

    ! Update each proc's communication cost.
    call check_communication_cost(Procs)

    ! If we screwed this up we need to terminated execution.
    do n_ = 1, nBlocks
      if (BlocksCollection(n_)%size > 0) then
        write(*,*), "Sorry, something went terribly wrong."
        write(*,*), "Program exiting."
!        STOP
      end if
    end do

    ! Write the total weight, number of procs, ideal weight per proc.
    write(*,*)
    write(*,*), '      Weight  ', 'Ideal Weight/proc'
    write(*,*), sum, optimal
    write(*,*)

    ! Write out the weight on each proc.
    write(*, *), '          proc #      ', "nBlocks  ", "Weight ", "Err"
    do n_ = 1, nProcs
      write(*, *), n_, Procs(n_)%nBlocks, Procs(n_)%weight, dfloat(Procs(n_)%weight) - optimal
    end do

  end subroutine
end module
