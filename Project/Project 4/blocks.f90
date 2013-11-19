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

        ! High/low i/j corresponds to the
        ! global block number start and stop.
        highI = 1 + iN * (iBlockSize - 1)
        b%lowI = highI - (iBlockSize - 1)
        highJ = 1 + iM * (jBlockSize - 1)
        b%lowJ = highJ - (jBlockSize - 1)

        ! Assume all faces are internal and find out who
        ! their neighbor is. Set the neighbor block and set
        ! their processor to 1 for now.
        b%northFace%neighborBlock = n_ + N
        b%southFace%neighborBlock = n_ - N
        b%eastFace%neighborBlock = n_ + 1
        b%westFace%neighborBlock = n_ - 1

        ! Assume all corners are internal (BC=-1) and set their neighbor.
        ! Set their processor to 1 for now.
        b%NECorner%BC = -1
        b%NECorner%neighborBlock = n_ + N + 1
        b%SECorner%BC = -1
        b%SECorner%neighborBlock = n_ - N + 1
        b%SWCorner%BC = -1
        b%SWCorner%neighborBlock = n_ - N - 1
        b%NWCorner%BC = -1
        b%NWCorner%neighborBlock = n_ + N - 1

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
          b%northFace%neighborBlock = -1
          b%northFace%neighborProc = -1
        end if

        if (b%southFace%BC /= -1) then
          b%southFace%neighborBlock = -1
          b%southFace%neighborProc = -1
        end if

        if (b%eastFace%BC /= -1) then
          b%eastFace%neighborBlock = -1
          b%eastFace%neighborProc = -1
        end if

        if (b%westFace%BC /= -1) then
          b%westFace%neighborBlock = -1
          b%westFace%neighborProc = -1
        end if

        ! If corners are on boundary then they don't
        ! have associated blocks. Set their processor
        ! to -1 so we know that there's nowhere to look.
        ! North East Corner
        if (b%northFace%BC == nB) then
          b%NECorner%BC = nBound
          b%NECorner%neighborBlock = -1
          b%NECorner%neighborProc = -1
        else if (b%eastFace%BC == eB) then
          b%NECorner%BC = eBound
          b%NECorner%neighborBlock = -1
          b%NECorner%neighborProc = -1
        end if

        ! South East Corner
        if (b%southFace%BC == sB) then
          b%SECorner%BC = sBound
          b%SECorner%neighborBlock = -1
          b%SECorner%neighborProc = -1
        else if (b%eastFace%BC == eB) then
          b%SECorner%BC = eBound
          b%SECorner%neighborBlock = -1
          b%SECorner%neighborProc = -1
        end if

        ! South West Corner
        if (b%southFace%BC == sB) then
          b%SWCorner%BC = sBound
          b%SWCorner%neighborBlock = -1
          b%SWCorner%neighborProc = -1
        else if (b%westFace%BC == wB) then
          b%SWCorner%BC = wBound
          b%SWCorner%neighborBlock = -1
          b%SWCorner%neighborProc = -1
        end if

        ! North West Corner
        if (b%northFace%BC == nB) then
          b%NWCorner%BC = nBound
          b%NWCorner%neighborBlock = -1
          b%NWCorner%neighborProc = -1
        else if (b%westFace%BC == wB) then
          b%NWCorner%BC = wBound
          b%NWCorner%neighborBlock = -1
          b%NWCorner%neighborProc = -1
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
        ! We calculate the 'area' of the block.
        b%size = (b%localIMAX - b%localIMIN) * (b%localJMAX - b%localJMIN)

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

  subroutine update_neighbor_procs(Procs)
    type (Proc), target :: Procs(:)
    type (Proc), pointer :: MyProc, ThisProc
    type (BlockType), pointer :: MyBlock, ThisBlock
    integer :: b, p, b2, p2

    do p = 1, nProcs
      MyProc => Procs(p)

      ! Set proc ids. (Start at 0.)
      MyProc%procID = p-1

      ! All blocks on this proc are on this proc.
      MyProc%Blocks%proc = MyProc%procID

      ! Default all neighbors to be on this proc as well.
      MyProc%Blocks%northFace%neighborProc = MyProc%procID
      MyProc%Blocks%southFace%neighborProc = MyProc%procID
      MyProc%Blocks%eastFace%neighborProc = MyProc%procID
      MyProc%Blocks%westFace%neighborProc = MyProc%procID

      MyProc%Blocks%NECorner%neighborProc = MyProc%procID
      MyProc%Blocks%NWCorner%neighborProc = MyProc%procID
      MyProc%Blocks%SECorner%neighborProc = MyProc%procID
      MyProc%Blocks%SWCorner%neighborProc = MyProc%procID
    end do

    ! Now we find which procs all our blocks neighbors live on.
    ! We loop through all procs...
    do p = 1, nProcs
      MyProc => Procs(p)

      ! ...we loop through each proc's blocks...
      do b = 1, Procs(p)%nBlocks
        MyBlock => MyProc%Blocks(b)

        ! ...we check each proc...
        do p2 = 1, nProcs
          ThisProc => Procs(p2)

          if (MyProc%procID /= ThisProc%procID) then
            ! ...and compare this blocks id to our faces id.
            ! When we find a match we set our faces proc to the proc this block is on.
            do b2 = 1, Procs(p2)%nBlocks
              ThisBlock => ThisProc%Blocks(b2)

              ! Check each face and corner for a match and assign.
              if (MyBlock%northFace%neighborBlock == ThisBlock%id) MyBlock%northFace%neighborProc = ThisProc%procID
              if (MyBlock%southFace%neighborBlock == ThisBlock%id) MyBlock%southFace%neighborProc = ThisProc%procID
              if (MyBlock%eastFace%neighborBlock  == ThisBlock%id) MyBlock%eastFace%neighborProc  = ThisProc%procID
              if (MyBlock%westFace%neighborBlock  == ThisBlock%id) MyBlock%westFace%neighborProc  = ThisProc%procID

              if (MyBlock%NECorner%neighborBlock == ThisBlock%id)  MyBlock%NECorner%neighborProc  = ThisProc%procID
              if (MyBlock%NWCorner%neighborBlock == ThisBlock%id)  MyBlock%NWCorner%neighborProc  = ThisProc%procID
              if (MyBlock%SECorner%neighborBlock == ThisBlock%id)  MyBlock%SECorner%neighborProc  = ThisProc%procID
              if (MyBlock%SWCorner%neighborBlock == ThisBlock%id)  MyBlock%SWCorner%neighborProc  = ThisProc%procID
            end do

          end if
        end do
      end do
    end do
  end subroutine

  subroutine check_communication_cost(Procs)
    type (Proc), target :: Procs(:)
    type (Proc), pointer :: MyProc
    type (BlockType), pointer :: MyBlock
    integer :: comm, b, p, i
    integer, pointer :: northProc, southProc, eastProc, westProc, NEProc, NWProc, SEProc, SWProc

    ! Check if this block has neighbors on the processor.
    ! If so we go ahead and add the communication cost.
    ! For each proc...
    do p = 1, nProcs
      MyProc => Procs(p)
      MyProc%comm = 0
      comm = 0

      ! ...we loop over each block.
      do b = 1, MyProc%nBlocks
        MyBlock => MyProc%Blocks(b)

        ! Check to see if this block is a neighbor to another block on this processor.
        ! The communication cost is equal to the length of the sides of the block.
        northProc => MyBlock%northFace%neighborProc
        southProc => MyBlock%southFace%neighborProc
        eastProc => MyBlock%eastFace%neighborProc
        westProc => MyBlock%westFace%neighborProc
        NEProc => MyBlock%NECorner%neighborProc
        NWProc => MyBlock%NWCorner%neighborProc
        SEProc => MyBlock%SECorner%neighborProc
        SWProc => MyBlock%SWCorner%neighborProc

        if (northProc /= MyProc%procID .and. northProc /= -1) comm = comm + ( MyBlock%localJMAX - MyBlock%localJMIN )
        if (southProc /= MyProc%procID .and. southProc /= -1) comm = comm + ( MyBlock%localJMAX - MyBlock%localJMIN )
        if (eastProc  /= MyProc%procID .and. eastProc /= -1)  comm = comm + ( MyBlock%localIMAX - MyBlock%localIMIN )
        if (westProc  /= MyProc%procID .and. westProc /= -1)  comm = comm + ( MyBlock%localIMAX - MyBlock%localIMIN )

        if (NEProc /= MyProc%procID .and. NEProc /= -1)       comm = comm + 1
        if (NWProc /= MyProc%procID .and. NWProc /= -1)       comm = comm + 1
        if (SEProc /= MyProc%procID .and. SEProc /= -1)       comm = comm + 1
        if (SWProc /= MyProc%procID .and. SWProc /= -1)       comm = comm + 1
      end do

      ! Set the communication cost.
      MyProc%comm = comm
    end do
  end subroutine

  ! Helper subroutine to allow us to pass a list of ids.
  subroutine add_blocks_to_proc(MyProc, BlocksCollection, BlockIDs)
    type (Proc) :: MyProc
    type (BlockType), target :: BlocksCollection(:)
    integer :: BlockIDs(:), i

    do i = 1, size(BlockIDs)
        call add_block_to_proc(MyProc, BlocksCollection(BlockIDs(i)))
    end do
  end subroutine

  ! Various actions to add a block to a processor.
  subroutine add_block_to_proc(MyProc, Block)
    type (Proc) :: MyProc
    type (BlockType) :: Block

    MyProc%nBlocks = MyProc%nBlocks + 1
    MyProc%weight = MyProc%weight + Block%size

    MyProc%Blocks(MyProc%nBlocks) = Block
    MyProc%Blocks(MyProc%nBlocks)%proc = MyProc%procID

    ! Setting the old block size to zero to note it has been assigned.
    Block%size = 0
  end subroutine

  subroutine distribute_blocks(BlocksCollection, Procs)
    type (BlockType), target :: BlocksCollection(:)
    type (Proc), allocatable :: Procs(:)
    type (BlockType), pointer :: b
    real(kind=8) :: fudge_factor = 1.1d0
    integer :: optimal, method, p_, n_, sum = 0
    integer :: largest_block = 0, largest_block_number = 0

    ! Set starting weights on each proc equal to zero.
    Procs%weight = 0
    Procs%nBlocks = 0

    ! Find the sum of the weights of the blocks with no communication cost.
    ! This is an ideal we cannot meet.
    do n_ = 1, nBlocks
      b => BlocksCollection(n_)
      sum = sum + b%size
    end do

    ! The optimal distribution is an equal weight on each block.
    optimal = int(dfloat(sum)/dfloat(nProcs))

    ! Pick our method depending on our settings.
    ! Hard coded for all the decompositions we're interested in.
    if (M == 5  .and. N == 4) then
      if (nProcs == 6) then
        method = 546
      else if (nProcs == 4) then
        method = 544
      end if
    else if (M == 10 .and. N == 10) then
      if (nProcs == 6) then
        method = 10106
      else if (nProcs == 4) then
        method = 10104
      end if
    else
      ! Otherwise we'll do our best and use an automated method.
      method = 666
    end if

    if (method == 546) then
      ! 5x4 blocks on 6 processors.
      ! Hard coding is hard work.
      call add_blocks_to_proc(Procs(1), BlocksCollection, [1, 2, 3, 4])
      call add_blocks_to_proc(Procs(2), BlocksCollection, [5, 6, 7])
      call add_blocks_to_proc(Procs(3), BlocksCollection, [9, 10, 11])
      call add_blocks_to_proc(Procs(4), BlocksCollection, [13, 14, 15])
      call add_blocks_to_proc(Procs(5), BlocksCollection, [8, 12, 16])
      call add_blocks_to_proc(Procs(6), BlocksCollection, [17, 18, 19, 20])

    else if (method == 544) then
      ! 5x4 blocks on 4 processors. Break the domain into columns.
      ! Hard coding is still hard work.
      call add_blocks_to_proc(Procs(1), BlocksCollection, [1, 5, 9, 13, 17])
      call add_blocks_to_proc(Procs(2), BlocksCollection, [2, 6, 10, 14, 18])
      call add_blocks_to_proc(Procs(3), BlocksCollection, [3, 7, 11, 15, 19])
      call add_blocks_to_proc(Procs(4), BlocksCollection, [4, 8, 12, 16, 20])

    else if (method == 10106) then
      ! 10x10 blocks on 6 processors.
      ! Hard coding is okay.
      do n_ = 0, 2
        call add_blocks_to_proc(Procs(1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10])

        call add_blocks_to_proc(Procs(2), BlocksCollection, &
          [6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])

        call add_blocks_to_proc(Procs(3), BlocksCollection, &
          [31 + n_*10, 32 + n_*10, 33 + n_*10, 34 + n_*10, 35 + n_*10])

        call add_blocks_to_proc(Procs(4), BlocksCollection, &
          [36 + n_*10, 37 + n_*10, 38 + n_*10, 39 + n_*10, 40 + n_*10])
      end do

      do n_ = 0, 3
        call add_blocks_to_proc(Procs(5), BlocksCollection, &
          [61 + n_*10, 62 + n_*10, 63 + n_*10, 64 + n_*10, 65 + n_*10])

        call add_blocks_to_proc(Procs(6), BlocksCollection, &
          [66 + n_*10, 67 + n_*10, 68 + n_*10, 69 + n_*10, 70 + n_*10])
      end do

    else if (method == 10104) then
      ! 10x10 blocks on 4 processors. Break the domain into quadrants.
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
    call update_neighbor_procs(Procs)
    call check_communication_cost(Procs)

    ! If we screwed this up we need to terminated execution.
    do n_ = 1, nBlocks
      if (BlocksCollection(n_)%size > 0) then
        write(*,*), n_
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
    write(*, *), '          proc #     ', "nBlocks   ", "Weight       ", "Comm        ", "Err"
    do n_ = 1, nProcs
      write(*, *), n_, Procs(n_)%nBlocks, Procs(n_)%weight, Procs(n_)%comm, Procs(n_)%weight + Procs(n_)%comm - optimal
    end do

  end subroutine
end module
