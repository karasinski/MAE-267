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
        ! their neighbor is. Set the neighbor block.
        b%northFace%neighborBlock = n_ + N
        b%southFace%neighborBlock = n_ - N
        b%eastFace%neighborBlock = n_ + 1
        b%westFace%neighborBlock = n_ - 1

        ! Assume all corners are internal (BC=-1) and set their neighbor.
        b%NECorner%BC = INTERNAL_BOUNDARY
        b%NECorner%neighborBlock = n_ + N + 1
        b%SECorner%BC = INTERNAL_BOUNDARY
        b%SECorner%neighborBlock = n_ - N + 1
        b%SWCorner%BC = INTERNAL_BOUNDARY
        b%SWCorner%neighborBlock = n_ - N - 1
        b%NWCorner%BC = INTERNAL_BOUNDARY
        b%NWCorner%neighborBlock = n_ + N - 1

        ! Set local block numbers to -1 for now.
        b%northFace%neighborLocalBlock = -1
        b%southFace%neighborLocalBlock = -1
        b%eastFace%neighborLocalBlock = -1
        b%westFace%neighborLocalBlock = - 1

        b%NECorner%neighborLocalBlock = - 1
        b%SECorner%neighborLocalBlock = - 1
        b%SWCorner%neighborLocalBlock = - 1
        b%NWCorner%neighborLocalBlock = - 1

        ! Determine whether block is on boundary or
        ! if they is internal for each face.
        if (highJ == JMAX) then
          nBound = nB
        else
          nBound = INTERNAL_BOUNDARY
        end if
        if (highI == IMAX) then
          eBound = eB
        else
          eBound = INTERNAL_BOUNDARY
        end if
        if (b%lowJ == 1) then
          sBound = sB
        else
          sBound = INTERNAL_BOUNDARY
        end if
        if (b%lowI == 1) then
          wBound = wB
        else
          wBound = INTERNAL_BOUNDARY
        end if

        ! Bounds and Faces
        ! Set our boundary conditions for each face.
        b%northFace%BC = nBound
        b%southFace%BC = sBound
        b%eastFace%BC = eBound
        b%westFace%BC = wBound

        ! Check if faces are on the boundary and
        ! if so, set them to point to nothing.
        if (b%northFace%BC /= INTERNAL_BOUNDARY) then
          b%northFace%neighborBlock = EXTERNAL_BOUNDARY
          b%northFace%neighborProc = EXTERNAL_BOUNDARY
        end if

        if (b%southFace%BC /= INTERNAL_BOUNDARY) then
          b%southFace%neighborBlock = EXTERNAL_BOUNDARY
          b%southFace%neighborProc = EXTERNAL_BOUNDARY
        end if

        if (b%eastFace%BC /= INTERNAL_BOUNDARY) then
          b%eastFace%neighborBlock = EXTERNAL_BOUNDARY
          b%eastFace%neighborProc = EXTERNAL_BOUNDARY
        end if

        if (b%westFace%BC /= INTERNAL_BOUNDARY) then
          b%westFace%neighborBlock = EXTERNAL_BOUNDARY
          b%westFace%neighborProc = EXTERNAL_BOUNDARY
        end if

        ! If corners are on boundary then they don't
        ! have associated blocks. Set their processor
        ! to -1 so we know that there's nowhere to look.
        ! North East Corner
        if (b%northFace%BC == nB) then
          b%NECorner%BC = nBound
          b%NECorner%neighborBlock = EXTERNAL_BOUNDARY
          b%NECorner%neighborProc = EXTERNAL_BOUNDARY
        else if (b%eastFace%BC == eB) then
          b%NECorner%BC = eBound
          b%NECorner%neighborBlock = EXTERNAL_BOUNDARY
          b%NECorner%neighborProc = EXTERNAL_BOUNDARY
        end if

        ! South East Corner
        if (b%southFace%BC == sB) then
          b%SECorner%BC = sBound
          b%SECorner%neighborBlock = EXTERNAL_BOUNDARY
          b%SECorner%neighborProc = EXTERNAL_BOUNDARY
        else if (b%eastFace%BC == eB) then
          b%SECorner%BC = eBound
          b%SECorner%neighborBlock = EXTERNAL_BOUNDARY
          b%SECorner%neighborProc = EXTERNAL_BOUNDARY
        end if

        ! South West Corner
        if (b%southFace%BC == sB) then
          b%SWCorner%BC = sBound
          b%SWCorner%neighborBlock = EXTERNAL_BOUNDARY
          b%SWCorner%neighborProc = EXTERNAL_BOUNDARY
        else if (b%westFace%BC == wB) then
          b%SWCorner%BC = wBound
          b%SWCorner%neighborBlock = EXTERNAL_BOUNDARY
          b%SWCorner%neighborProc = EXTERNAL_BOUNDARY
        end if

        ! North West Corner
        if (b%northFace%BC == nB) then
          b%NWCorner%BC = nBound
          b%NWCorner%neighborBlock = EXTERNAL_BOUNDARY
          b%NWCorner%neighborProc = EXTERNAL_BOUNDARY
        else if (b%westFace%BC == wB) then
          b%NWCorner%BC = wBound
          b%NWCorner%neighborBlock = EXTERNAL_BOUNDARY
          b%NWCorner%neighborProc = EXTERNAL_BOUNDARY
        end if

        ! Calculate 'true' bounds of blocks.
        ! North face.
        if (b%northFace%BC == INTERNAL_BOUNDARY) then
          b%localJMAX = jBlockSize
        else
          b%localJMAX = jBlockSize - 1
        end if

        ! South face.
        if (b%southFace%BC == INTERNAL_BOUNDARY) then
          b%localJMIN = 0
        else
          b%localJMIN = 1
        end if

        ! East face.
        if (b%eastFace%BC == INTERNAL_BOUNDARY) then
          b%localIMAX = iBlockSize
        else
          b%localIMAX = iBlockSize - 1
        end if

        ! West face.
        if (b%westFace%BC == INTERNAL_BOUNDARY) then
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
    type (GridPoint), pointer :: p1
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
          p1 => BlocksCollection(n_)%Points(i, jBlockSize)
          p1%T = 5.d0 * (sin(pi * p(i, jBlockSize)%xp) + 1.d0)
        end do
      end if

      if (b%southFace%BC == sB) then
        do i = 1, iBlockSize
          p1 => BlocksCollection(n_)%Points(i, 1)
          p1%T = abs(cos(pi * p(i,1)%xp)) + 1.d0
        end do
      end if

      if (b%eastFace%BC == eB) then
        do j = 1, jBlockSize
          p1 => BlocksCollection(n_)%Points(iBlockSize, j)
          p1%T = (3.d0 * p(iBlockSize, j)%yp) + 2.d0
        end do
      end if

      if (b%westFace%BC == wB) then
        do j = 1, jBlockSize
          p1 => BlocksCollection(n_)%Points(1, j)
          p1%T = (3.d0 * p(1, j)%yp) + 2.d0
        end do
      end if

    end do
  end subroutine

  subroutine update_neighbor_procs(Procs)
    type (Proc), target :: Procs(:)
    type (Proc), pointer :: MyProc, ThisProc
    type (BlockType), pointer :: MyBlock, ThisBlock
    integer :: b, p, b2, p2

    do p = 1, mpi_nprocs
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
    do p = 1, mpi_nprocs
      MyProc => Procs(p)

      ! ...we loop through each proc's blocks...
      do b = 1, Procs(p)%nBlocks
        MyBlock => MyProc%Blocks(b)

        ! ...we check each proc...
        do p2 = 1, mpi_nprocs
          ThisProc => Procs(p2)

            ! ...and compare this blocks id to our faces id.
            ! When we find a match we set our faces proc to the proc this block is on.
            do b2 = 1, Procs(p2)%nBlocks
              ThisBlock => ThisProc%Blocks(b2)

              ! Check each face and corner for a match and assign.
              if (MyBlock%northFace%neighborBlock == ThisBlock%id) then
                MyBlock%northFace%neighborProc = ThisProc%procID
                MyBlock%northFace%neighborLocalBlock = b2
                
                if (MyProc%procID /= ThisProc%procID) then
                  MyBlock%northFace%BC = PROC_BOUNDARY
                end if
              end if
              if (MyBlock%southFace%neighborBlock == ThisBlock%id) then 
                MyBlock%southFace%neighborProc = ThisProc%procID
                MyBlock%southFace%neighborLocalBlock = b2
                if (MyProc%procID /= ThisProc%procID) then
                  MyBlock%southFace%BC = PROC_BOUNDARY
                end if
              end if
              if (MyBlock%eastFace%neighborBlock  == ThisBlock%id) then 
                MyBlock%eastFace%neighborProc  = ThisProc%procID
                MyBlock%eastFace%neighborLocalBlock = b2
                if (MyProc%procID /= ThisProc%procID) then
                  MyBlock%eastFace%BC = PROC_BOUNDARY
                end if
              end if
              if (MyBlock%westFace%neighborBlock  == ThisBlock%id) then 
                MyBlock%westFace%neighborProc  = ThisProc%procID
                MyBlock%westFace%neighborLocalBlock = b2
                if (MyProc%procID /= ThisProc%procID) then
                  MyBlock%westFace%BC = PROC_BOUNDARY
                end if
              end if

              if (MyBlock%NECorner%neighborBlock == ThisBlock%id)  then 
                MyBlock%NECorner%neighborProc  = ThisProc%procID
                MyBlock%NECorner%neighborLocalBlock = b2
                if (MyProc%procID /= ThisProc%procID) then
                  MyBlock%NECorner%BC = PROC_BOUNDARY
                end if
              end if
              if (MyBlock%NWCorner%neighborBlock == ThisBlock%id)  then
                MyBlock%NWCorner%neighborProc  = ThisProc%procID
                MyBlock%NWCorner%neighborLocalBlock = b2
                if (MyProc%procID /= ThisProc%procID) then
                  MyBlock%NWCorner%BC = PROC_BOUNDARY
                end if
              end if
              if (MyBlock%SECorner%neighborBlock == ThisBlock%id)  then 
                MyBlock%SECorner%neighborProc  = ThisProc%procID
                MyBlock%SECorner%neighborLocalBlock = b2
                if (MyProc%procID /= ThisProc%procID) then
                  MyBlock%SECorner%BC = PROC_BOUNDARY
                end if
              end if
              if (MyBlock%SWCorner%neighborBlock == ThisBlock%id)  then 
                MyBlock%SWCorner%neighborProc  = ThisProc%procID
                MyBlock%SWCorner%neighborLocalBlock = b2
                if (MyProc%procID /= ThisProc%procID) then
                  MyBlock%SWCorner%BC = PROC_BOUNDARY
                end if
              end if
            end do

        end do
      end do
    end do
  end subroutine

  subroutine check_communication_cost(Procs)
    type (Proc), pointer :: Procs(:)
    type (Proc), pointer :: MyProc
    type (BlockType), pointer :: MyBlock
    integer :: comm, b, p
    integer, pointer :: northProc, southProc, eastProc, westProc, NEProc, NWProc, SEProc, SWProc

    ! Check if this block has neighbors on the processor.
    ! If so we go ahead and add the communication cost.
    ! For each proc...
    do p = 1, mpi_nprocs
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
    type (Proc), pointer :: Procs(:)
    type (BlockType), pointer :: b
    real(kind=8) :: fudge_factor = 1.05d0
    integer :: optimal, method = 666, p_, n_, sum = 0
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
    optimal = int(dfloat(sum)/dfloat(mpi_nprocs))

    ! Pick our method depending on our settings.
    ! Hard coded for all the decompositions we're interested in.
    if (M == 10 .and. N == 10) then
      if (mpi_nprocs == 10) then
        method = 101010
      else if (mpi_nprocs == 8) then
        fudge_factor = 10108
      else if (mpi_nprocs == 6) then
        method = 10106
      else if (mpi_nprocs == 4) then
        method = 10104
      else if (mpi_nprocs == 2) then
        method = 10102
      end if
    else
      ! Otherwise we'll do our best and use an automated method.
      method = 666
    end if

    if (method == 101010) then
      ! 10x10 blocks on 10 processors.

      do n_ = 0, 9 
        call add_blocks_to_proc(Procs(n_+1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10, &
           6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])
      end do

    else if (method == 10108) then
      ! 10x10 blocks on 8 processors.

      call add_blocks_to_proc(Procs(1), BlocksCollection, &
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
      call add_blocks_to_proc(Procs(2), BlocksCollection, &
        [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
      call add_blocks_to_proc(Procs(3), BlocksCollection, &
        [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38])
      call add_blocks_to_proc(Procs(4), BlocksCollection, &
        [39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50])
      call add_blocks_to_proc(Procs(5), BlocksCollection, &
        [51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63])
      call add_blocks_to_proc(Procs(6), BlocksCollection, &
        [64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75])
      call add_blocks_to_proc(Procs(7), BlocksCollection, &
        [76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88])
      call add_blocks_to_proc(Procs(8), BlocksCollection, &
        [89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100])

    else if (method == 10106) then
      ! 10x10 blocks on 6 processors.
      do n_ = 0, 2
        call add_blocks_to_proc(Procs(1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10])

        call add_blocks_to_proc(Procs(2), BlocksCollection, &
          [6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])
      end do

      call add_blocks_to_proc(Procs(1), BlocksCollection, &
        [34, 35])

      call add_blocks_to_proc(Procs(2), BlocksCollection, &
        [36, 37])

      call add_blocks_to_proc(Procs(3), BlocksCollection, &
        [31, 32, 33])

      call add_blocks_to_proc(Procs(4), BlocksCollection, &
        [38, 39, 40])

      do n_ = 1, 2
        call add_blocks_to_proc(Procs(3), BlocksCollection, &
          [31 + n_*10, 32 + n_*10, 33 + n_*10, 34 + n_*10, 35 + n_*10])

        call add_blocks_to_proc(Procs(4), BlocksCollection, &
          [36 + n_*10, 37 + n_*10, 38 + n_*10, 39 + n_*10, 40 + n_*10])
      end do

      call add_blocks_to_proc(Procs(3), BlocksCollection, &
        [61, 62, 63])

      call add_blocks_to_proc(Procs(4), BlocksCollection, &
        [68, 69, 70])

      call add_blocks_to_proc(Procs(5), BlocksCollection, &
        [64, 65])

      call add_blocks_to_proc(Procs(6), BlocksCollection, &
        [66, 67])

      do n_ = 1, 3
        call add_blocks_to_proc(Procs(5), BlocksCollection, &
          [61 + n_*10, 62 + n_*10, 63 + n_*10, 64 + n_*10, 65 + n_*10])

        call add_blocks_to_proc(Procs(6), BlocksCollection, &
          [66 + n_*10, 67 + n_*10, 68 + n_*10, 69 + n_*10, 70 + n_*10])
      end do

    else if (method == 10104) then
      ! 10x10 blocks on 4 processors. Break the domain into quadrants.
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

    else if (method == 10102) then
      ! 10x10 blocks on 4 processors. Break the domain into quadrants.
      do n_ = 0, 4
        call add_blocks_to_proc(Procs(1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10, &
           6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])

        call add_blocks_to_proc(Procs(2), BlocksCollection, &
          [51 + n_*10, 52 + n_*10, 53 + n_*10, 54 + n_*10, 55 + n_*10, &
           56 + n_*10, 57 + n_*10, 58 + n_*10, 59 + n_*10, 60 + n_*10])
      end do

    else if (method == 666) then
      ! This is our fallback method for unknown configurations.

      do p_ = 1, mpi_nprocs
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

          call add_block_to_proc(Procs(p_), BlocksCollection(largest_block_number))
        end do
      end do

    end if

    ! Update each proc's communication cost.
    call update_neighbor_procs(Procs)
    call check_communication_cost(Procs)

    ! Write the total weight, number of procs, ideal weight per proc.
    write(*,*)
    write(*,*), '      Weight  ', 'Ideal Weight/proc'
    write(*,*), sum, optimal
    write(*,*)

    ! Write out the weight on each proc.
    write(*, *), '     proc #     ', "nBlocks   ", "   Weight        ", "Comm        ", "  Err"
    do n_ = 1, mpi_nprocs
      write(*, *), n_, Procs(n_)%nBlocks, Procs(n_)%weight, Procs(n_)%comm, &
                  100*real(Procs(n_)%weight + Procs(n_)%comm - optimal)/real(optimal)
    end do

    ! If we screwed this up we need to terminated execution.
    do n_ = 1, nBlocks
      if (BlocksCollection(n_)%size > 0) then
        write(*,*), n_
        write(*,*), "Sorry, something went terribly wrong."
        write(*,*), "Program exiting."
        STOP
      end if
    end do

  end subroutine
end module
