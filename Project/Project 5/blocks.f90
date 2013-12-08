! This module contains all the subroutines used to
! create our blocks and initial grid/temp files.
MODULE GridCreation
  USE BlockModule
  IMPLICIT NONE

CONTAINS
  SUBROUTINE create_blocks(BlocksCollection)
    TYPE (BlockType), TARGET :: BlocksCollection(:)
    TYPE (BlockType), POINTER :: b
    INTEGER :: n_, proc = 1, IN, IM
    INTEGER :: nBound, sBound, eBound, wBound
    INTEGER :: highI, highJ

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
    DO IM = 1, M
      DO IN = 1, N
        b => BlocksCollection(n_)

        ! Global id.
        b%ID = n_

        ! High/low i/j corresponds to the
        ! global block number start and stop.
        highI = 1 + IN * (iBlockSize - 1)
        b%lowI = highI - (iBlockSize - 1)
        highJ = 1 + IM * (jBlockSize - 1)
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
        IF (highJ == JMAX) THEN
          nBound = nB
        ELSE
          nBound = INTERNAL_BOUNDARY
        END IF
        IF (highI == IMAX) THEN
          eBound = eB
        ELSE
          eBound = INTERNAL_BOUNDARY
        END IF
        IF (b%lowJ == 1) THEN
          sBound = sB
        ELSE
          sBound = INTERNAL_BOUNDARY
        END IF
        IF (b%lowI == 1) THEN
          wBound = wB
        ELSE
          wBound = INTERNAL_BOUNDARY
        END IF

        ! Bounds and Faces
        ! Set our boundary conditions for each face.
        b%northFace%BC = nBound
        b%southFace%BC = sBound
        b%eastFace%BC = eBound
        b%westFace%BC = wBound

        ! Check if faces are on the boundary and
        ! if so, set them to point to nothing.
        IF (b%northFace%BC /= INTERNAL_BOUNDARY) THEN
          b%northFace%neighborBlock = EXTERNAL_BOUNDARY
          b%northFace%neighborProc = EXTERNAL_BOUNDARY
        END IF

        IF (b%southFace%BC /= INTERNAL_BOUNDARY) THEN
          b%southFace%neighborBlock = EXTERNAL_BOUNDARY
          b%southFace%neighborProc = EXTERNAL_BOUNDARY
        END IF

        IF (b%eastFace%BC /= INTERNAL_BOUNDARY) THEN
          b%eastFace%neighborBlock = EXTERNAL_BOUNDARY
          b%eastFace%neighborProc = EXTERNAL_BOUNDARY
        END IF

        IF (b%westFace%BC /= INTERNAL_BOUNDARY) THEN
          b%westFace%neighborBlock = EXTERNAL_BOUNDARY
          b%westFace%neighborProc = EXTERNAL_BOUNDARY
        END IF

        ! If corners are on boundary then they don't
        ! have associated blocks. Set their processor
        ! to -1 so we know that there's nowhere to look.
        ! North East Corner
        IF (b%northFace%BC == nB) THEN
          b%NECorner%BC = nBound
          b%NECorner%neighborBlock = EXTERNAL_BOUNDARY
          b%NECorner%neighborProc = EXTERNAL_BOUNDARY
        ELSE IF (b%eastFace%BC == eB) THEN
          b%NECorner%BC = eBound
          b%NECorner%neighborBlock = EXTERNAL_BOUNDARY
          b%NECorner%neighborProc = EXTERNAL_BOUNDARY
        END IF

        ! South East Corner
        IF (b%southFace%BC == sB) THEN
          b%SECorner%BC = sBound
          b%SECorner%neighborBlock = EXTERNAL_BOUNDARY
          b%SECorner%neighborProc = EXTERNAL_BOUNDARY
        ELSE IF (b%eastFace%BC == eB) THEN
          b%SECorner%BC = eBound
          b%SECorner%neighborBlock = EXTERNAL_BOUNDARY
          b%SECorner%neighborProc = EXTERNAL_BOUNDARY
        END IF

        ! South West Corner
        IF (b%southFace%BC == sB) THEN
          b%SWCorner%BC = sBound
          b%SWCorner%neighborBlock = EXTERNAL_BOUNDARY
          b%SWCorner%neighborProc = EXTERNAL_BOUNDARY
        ELSE IF (b%westFace%BC == wB) THEN
          b%SWCorner%BC = wBound
          b%SWCorner%neighborBlock = EXTERNAL_BOUNDARY
          b%SWCorner%neighborProc = EXTERNAL_BOUNDARY
        END IF

        ! North West Corner
        IF (b%northFace%BC == nB) THEN
          b%NWCorner%BC = nBound
          b%NWCorner%neighborBlock = EXTERNAL_BOUNDARY
          b%NWCorner%neighborProc = EXTERNAL_BOUNDARY
        ELSE IF (b%westFace%BC == wB) THEN
          b%NWCorner%BC = wBound
          b%NWCorner%neighborBlock = EXTERNAL_BOUNDARY
          b%NWCorner%neighborProc = EXTERNAL_BOUNDARY
        END IF

        ! Calculate 'true' bounds of blocks.
        ! North face.
        IF (b%northFace%BC == INTERNAL_BOUNDARY) THEN
          b%localJMAX = jBlockSize
        ELSE
          b%localJMAX = jBlockSize - 1
        END IF

        ! South face.
        IF (b%southFace%BC == INTERNAL_BOUNDARY) THEN
          b%localJMIN = 0
        ELSE
          b%localJMIN = 1
        END IF

        ! East face.
        IF (b%eastFace%BC == INTERNAL_BOUNDARY) THEN
          b%localIMAX = iBlockSize
        ELSE
          b%localIMAX = iBlockSize - 1
        END IF

        ! West face.
        IF (b%westFace%BC == INTERNAL_BOUNDARY) THEN
          b%localIMIN = 0
        ELSE
          b%localIMIN = 1
        END IF

        ! Weight of each block, used for load balancing.
        ! We calculate the 'area' of the block.
        b%SIZE = (b%localIMAX - b%localIMIN) * (b%localJMAX - b%localJMIN)

        n_ = n_ + 1
      END DO
    END DO
  END SUBROUTINE

  ! Intialize our locations and prime locations.
  SUBROUTINE initialize_block_grid(BlocksCollection)
    TYPE (BlockType), TARGET :: BlocksCollection(:)
    TYPE (GridPoint), POINTER :: p
    TYPE (BlockType), POINTER :: b
    INTEGER :: i, j, n_

    DO n_ = 1, nBlocks
      DO j = 0, jBlockSize+1
        DO i = 0, iBlockSize+1
          b => BlocksCollection(n_)
          p => BlocksCollection(n_)%Points(i, j)

          ! Have to convert from local i/j to global i/j.
          p%xp = COS( 0.5d0 * pi * DFLOAT(IMAX - (b%lowI + (i - 1 ))) / DFLOAT(IMAX-1))
          p%yp = COS( 0.5d0 * pi * DFLOAT(JMAX - (b%lowJ + (j - 1 ))) / DFLOAT(JMAX-1))

          p%x = p%xp * COS( rot ) + ( 1.d0 - p%yp ) * SIN( rot )
          p%y = p%yp * COS( rot ) + p%xp * SIN( rot )
        END DO
      END DO
    END DO
  END SUBROUTINE

  ! Dirichlet conditions.
  SUBROUTINE initialize_block_temp(BlocksCollection)
    TYPE (BlockType), TARGET :: BlocksCollection(:)
    TYPE (BlockType), POINTER :: b
    TYPE (GridPoint), POINTER :: p(:,:)
    TYPE (GridPoint), POINTER :: p1
    INTEGER :: i, j, n_
    REAL(KIND=8) :: T_0 = 3.5d0

    DO n_ = 1, nBlocks
      b => BlocksCollection(n_)
      p => BlocksCollection(n_)%Points

      ! Initialize all points to T_0.
      p%T = T_0

      ! Set the boundary conditions if our blocks are on the boundary.
      IF (b%northFace%BC == nB) THEN
        DO i = 1, iBlockSize
          p1 => BlocksCollection(n_)%Points(i, jBlockSize)
          p1%T = 5.d0 * (SIN(pi * p(i, jBlockSize)%xp) + 1.d0)
        END DO
      END IF

      IF (b%southFace%BC == sB) THEN
        DO i = 1, iBlockSize
          p1 => BlocksCollection(n_)%Points(i, 1)
          p1%T = ABS(COS(pi * p(i,1)%xp)) + 1.d0
        END DO
      END IF

      IF (b%eastFace%BC == eB) THEN
        DO j = 1, jBlockSize
          p1 => BlocksCollection(n_)%Points(iBlockSize, j)
          p1%T = (3.d0 * p(iBlockSize, j)%yp) + 2.d0
        END DO
      END IF

      IF (b%westFace%BC == wB) THEN
        DO j = 1, jBlockSize
          p1 => BlocksCollection(n_)%Points(1, j)
          p1%T = (3.d0 * p(1, j)%yp) + 2.d0
        END DO
      END IF

    END DO
  END SUBROUTINE

  SUBROUTINE update_neighbor_procs(Procs)
    TYPE (Proc), TARGET :: Procs(:)
    TYPE (Proc), POINTER :: MyProc, ThisProc
    TYPE (BlockType), POINTER :: MyBlock, ThisBlock
    INTEGER :: b, p, b2, p2

    DO p = 1, mpi_nprocs
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
    END DO

    ! Now we find which procs all our blocks neighbors live on.
    ! We loop through all procs...
    DO p = 1, mpi_nprocs
      MyProc => Procs(p)

      ! ...we loop through each proc's blocks...
      DO b = 1, Procs(p)%nBlocks
        MyBlock => MyProc%Blocks(b)

        ! ...we check each proc...
        DO p2 = 1, mpi_nprocs
          ThisProc => Procs(p2)

          ! ...and compare this blocks id to our faces id.
          ! When we find a match we set our faces proc to the proc this block is on.
          DO b2 = 1, Procs(p2)%nBlocks
            ThisBlock => ThisProc%Blocks(b2)

            ! Check each face and corner for a match and assign.
            IF (MyBlock%northFace%neighborBlock == ThisBlock%ID) THEN
              MyBlock%northFace%neighborProc = ThisProc%procID
              MyBlock%northFace%neighborLocalBlock = b2

              IF (MyProc%procID /= ThisProc%procID) THEN
                MyBlock%northFace%BC = PROC_BOUNDARY
              END IF
            END IF
            IF (MyBlock%southFace%neighborBlock == ThisBlock%ID) THEN
              MyBlock%southFace%neighborProc = ThisProc%procID
              MyBlock%southFace%neighborLocalBlock = b2
              IF (MyProc%procID /= ThisProc%procID) THEN
                MyBlock%southFace%BC = PROC_BOUNDARY
              END IF
            END IF
            IF (MyBlock%eastFace%neighborBlock  == ThisBlock%ID) THEN
              MyBlock%eastFace%neighborProc  = ThisProc%procID
              MyBlock%eastFace%neighborLocalBlock = b2
              IF (MyProc%procID /= ThisProc%procID) THEN
                MyBlock%eastFace%BC = PROC_BOUNDARY
              END IF
            END IF
            IF (MyBlock%westFace%neighborBlock  == ThisBlock%ID) THEN
              MyBlock%westFace%neighborProc  = ThisProc%procID
              MyBlock%westFace%neighborLocalBlock = b2
              IF (MyProc%procID /= ThisProc%procID) THEN
                MyBlock%westFace%BC = PROC_BOUNDARY
              END IF
            END IF

            IF (MyBlock%NECorner%neighborBlock == ThisBlock%ID)  THEN
              MyBlock%NECorner%neighborProc  = ThisProc%procID
              MyBlock%NECorner%neighborLocalBlock = b2
              IF (MyProc%procID /= ThisProc%procID) THEN
                MyBlock%NECorner%BC = PROC_BOUNDARY
              END IF
            END IF
            IF (MyBlock%NWCorner%neighborBlock == ThisBlock%ID)  THEN
              MyBlock%NWCorner%neighborProc  = ThisProc%procID
              MyBlock%NWCorner%neighborLocalBlock = b2
              IF (MyProc%procID /= ThisProc%procID) THEN
                MyBlock%NWCorner%BC = PROC_BOUNDARY
              END IF
            END IF
            IF (MyBlock%SECorner%neighborBlock == ThisBlock%ID)  THEN
              MyBlock%SECorner%neighborProc  = ThisProc%procID
              MyBlock%SECorner%neighborLocalBlock = b2
              IF (MyProc%procID /= ThisProc%procID) THEN
                MyBlock%SECorner%BC = PROC_BOUNDARY
              END IF
            END IF
            IF (MyBlock%SWCorner%neighborBlock == ThisBlock%ID)  THEN
              MyBlock%SWCorner%neighborProc  = ThisProc%procID
              MyBlock%SWCorner%neighborLocalBlock = b2
              IF (MyProc%procID /= ThisProc%procID) THEN
                MyBlock%SWCorner%BC = PROC_BOUNDARY
              END IF
            END IF
          END DO

        END DO
      END DO
    END DO
  END SUBROUTINE

  SUBROUTINE check_communication_cost(Procs)
    TYPE (Proc), POINTER :: Procs(:)
    TYPE (Proc), POINTER :: MyProc
    TYPE (BlockType), POINTER :: MyBlock
    INTEGER :: comm, b, p
    INTEGER, POINTER :: northProc, southProc, eastProc, westProc, NEProc, NWProc, SEProc, SWProc

    ! Check if this block has neighbors on the processor.
    ! If so we go ahead and add the communication cost.
    ! For each proc...
    DO p = 1, mpi_nprocs
      MyProc => Procs(p)
      MyProc%comm = 0
      comm = 0

      ! ...we loop over each block.
      DO b = 1, MyProc%nBlocks
        MyBlock => MyProc%Blocks(b)

        ! Check to see if this block is a neighbor to another block on this processor.
        ! The communication cost is equal to the number of points on the sides of the block.
        northProc => MyBlock%northFace%BC
        southProc => MyBlock%southFace%BC
        eastProc  => MyBlock%eastFace%BC
        westProc  => MyBlock%westFace%BC
        NEProc    => MyBlock%NECorner%BC
        NWProc    => MyBlock%NWCorner%BC
        SEProc    => MyBlock%SECorner%BC
        SWProc    => MyBlock%SWCorner%BC

        IF (northProc == PROC_BOUNDARY) comm = comm + iBlockSize
        IF (southProc == PROC_BOUNDARY) comm = comm + iBlockSize
        IF (eastProc  == PROC_BOUNDARY) comm = comm + jBlockSize
        IF (westProc  == PROC_BOUNDARY) comm = comm + jBlockSize
        IF (NEProc == PROC_BOUNDARY)    comm = comm + 1
        IF (NWProc == PROC_BOUNDARY)    comm = comm + 1
        IF (SEProc == PROC_BOUNDARY)    comm = comm + 1
        IF (SWProc == PROC_BOUNDARY)    comm = comm + 1
      END DO

      ! Set the communication cost.
      MyProc%comm = comm
    END DO
  END SUBROUTINE

  ! Helper subroutine to allow us to pass a list of ids.
  SUBROUTINE add_blocks_to_proc(MyProc, BlocksCollection, BlockIDs)
    TYPE (Proc) :: MyProc
    TYPE (BlockType), TARGET :: BlocksCollection(:)
    INTEGER :: BlockIDs(:), i

    DO i = 1, SIZE(BlockIDs)
      CALL add_block_to_proc(MyProc, BlocksCollection(BlockIDs(i)))
    END DO
  END SUBROUTINE

  ! Various actions to add a block to a processor.
  SUBROUTINE add_block_to_proc(MyProc, BLOCK)
    TYPE (Proc) :: MyProc
    TYPE (BlockType) :: BLOCK

    MyProc%nBlocks = MyProc%nBlocks + 1
    MyProc%weight = MyProc%weight + BLOCK%SIZE

    MyProc%Blocks(MyProc%nBlocks) = BLOCK
    MyProc%Blocks(MyProc%nBlocks)%proc = MyProc%procID

    ! Setting the old block size to zero to note it has been assigned.
    BLOCK%SIZE = 0
  END SUBROUTINE

  SUBROUTINE distribute_blocks(BlocksCollection, Procs)
    TYPE (BlockType), TARGET :: BlocksCollection(:)
    TYPE (Proc), POINTER :: Procs(:)
    TYPE (BlockType), POINTER :: b
    REAL(KIND=8) :: fudge_factor = 1.05d0, REAL_LOAD, IDEAL_LOAD
    INTEGER :: optimal, method = 666, p_, n_, LOAD_S = 0
    INTEGER :: LPP = 0, TLPP = 0
    REAL(KIND=8) :: PRED_IDEAL_SPEED = 0, PRED_REAL_SPEED = 0
    INTEGER :: largest_block = 0, largest_block_number = 0

    ! Set starting weights on each proc equal to zero.
    Procs%weight = 0
    Procs%nBlocks = 0

    ! Find the sum of the weights of the blocks with no communication cost.
    ! This is an ideal we cannot meet.
    DO n_ = 1, nBlocks
      b => BlocksCollection(n_)
      LOAD_S = LOAD_S + b%SIZE
    END DO

    ! The optimal distribution is an equal weight on each block.
    optimal = INT(DFLOAT(LOAD_S)/DFLOAT(mpi_nprocs))

    ! Pick our method depending on our settings.
    ! Hard coded for all the decompositions we're interested in.
    IF (M == 10 .and. N == 10) THEN
      IF (mpi_nprocs == 10) THEN
        method = 101010
      ELSE IF (mpi_nprocs == 8) THEN
        method = 10108
      ELSE IF (mpi_nprocs == 6) THEN
        method = 10106
      ELSE IF (mpi_nprocs == 4) THEN
        method = 10104
      ELSE IF (mpi_nprocs == 2) THEN
        method = 10102
      END IF
    ELSE
      ! Otherwise we'll do our best and use an automated method.
      method = 666
    END IF

    IF (method == 101010) THEN
      ! 10x10 blocks on 10 processors.

      DO n_ = 0, 9
        CALL add_blocks_to_proc(Procs(n_+1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5  + n_*10, &
           6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])
      END DO

    ELSE IF (method == 10108) THEN
      ! 10x10 blocks on 8 processors.

      CALL add_blocks_to_proc(Procs(1), BlocksCollection, &
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
      CALL add_blocks_to_proc(Procs(2), BlocksCollection, &
        [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
      CALL add_blocks_to_proc(Procs(3), BlocksCollection, &
        [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37])
      CALL add_blocks_to_proc(Procs(4), BlocksCollection, &
        [38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50])
      CALL add_blocks_to_proc(Procs(5), BlocksCollection, &
        [51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63])
      CALL add_blocks_to_proc(Procs(6), BlocksCollection, &
        [64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75])
      CALL add_blocks_to_proc(Procs(7), BlocksCollection, &
        [76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87])
      CALL add_blocks_to_proc(Procs(8), BlocksCollection, &
        [88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100])

    ELSE IF (method == 10106) THEN
      ! 10x10 blocks on 6 processors.
      DO n_ = 0, 2
        CALL add_blocks_to_proc(Procs(1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10])

        CALL add_blocks_to_proc(Procs(2), BlocksCollection, &
          [6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])
      END DO

      CALL add_blocks_to_proc(Procs(1), BlocksCollection, &
        [34, 35])

      CALL add_blocks_to_proc(Procs(2), BlocksCollection, &
        [36, 37])

      CALL add_blocks_to_proc(Procs(3), BlocksCollection, &
        [31, 32, 33])

      CALL add_blocks_to_proc(Procs(4), BlocksCollection, &
        [38, 39, 40])

      DO n_ = 1, 2
        CALL add_blocks_to_proc(Procs(3), BlocksCollection, &
          [31 + n_*10, 32 + n_*10, 33 + n_*10, 34 + n_*10, 35 + n_*10])

        CALL add_blocks_to_proc(Procs(4), BlocksCollection, &
          [36 + n_*10, 37 + n_*10, 38 + n_*10, 39 + n_*10, 40 + n_*10])
      END DO

      CALL add_blocks_to_proc(Procs(3), BlocksCollection, &
        [61, 62, 63])

      CALL add_blocks_to_proc(Procs(4), BlocksCollection, &
        [68, 69, 70])

      CALL add_blocks_to_proc(Procs(5), BlocksCollection, &
        [64, 65])

      CALL add_blocks_to_proc(Procs(6), BlocksCollection, &
        [66, 67])

      DO n_ = 1, 3
        CALL add_blocks_to_proc(Procs(5), BlocksCollection, &
          [61 + n_*10, 62 + n_*10, 63 + n_*10, 64 + n_*10, 65 + n_*10])

        CALL add_blocks_to_proc(Procs(6), BlocksCollection, &
          [66 + n_*10, 67 + n_*10, 68 + n_*10, 69 + n_*10, 70 + n_*10])
      END DO

    ELSE IF (method == 10104) THEN
      ! 10x10 blocks on 4 processors. Break the domain into quadrants.
      DO n_ = 0, 4
        CALL add_blocks_to_proc(Procs(1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5 + n_*10])

        CALL add_blocks_to_proc(Procs(2), BlocksCollection, &
          [6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])

        CALL add_blocks_to_proc(Procs(3), BlocksCollection, &
          [51 + n_*10, 52 + n_*10, 53 + n_*10, 54 + n_*10, 55 + n_*10])

        CALL add_blocks_to_proc(Procs(4), BlocksCollection, &
          [56 + n_*10, 57 + n_*10, 58 + n_*10, 59 + n_*10, 60 + n_*10])
      END DO

    ELSE IF (method == 10102) THEN
      ! 10x10 blocks on 4 processors. Break the domain into quadrants.
      DO n_ = 0, 4
        CALL add_blocks_to_proc(Procs(1), BlocksCollection, &
          [1 + n_*10, 2 + n_*10, 3 + n_*10, 4 + n_*10, 5  + n_*10, &
           6 + n_*10, 7 + n_*10, 8 + n_*10, 9 + n_*10, 10 + n_*10])

        CALL add_blocks_to_proc(Procs(2), BlocksCollection, &
          [51 + n_*10, 52 + n_*10, 53 + n_*10, 54 + n_*10, 55 + n_*10, &
           56 + n_*10, 57 + n_*10, 58 + n_*10, 59 + n_*10, 60 + n_*10])
      END DO

    ELSE IF (method == 666) THEN
      ! This is our fallback method for unknown configurations.

      DO p_ = 1, mpi_nprocs
        DO n_ = 1, nBlocks
          b => BlocksCollection(n_)

          ! Find the first nonzero size block.
          IF (b%SIZE > 0) THEN
            largest_block = b%SIZE
            largest_block_number = n_
          END IF
        END DO

        ! Check to see if proc has less than optimal weight.
        DO WHILE (Procs(p_)%weight + largest_block < fudge_factor * optimal)
          ! Reset largest block to zero.
          largest_block = 0

          ! Loop over all blocks.
          DO n_ = 1, nBlocks
            b => BlocksCollection(n_)

            ! Find the first nonzero size block.
            IF (b%SIZE > 0) THEN
              largest_block = b%SIZE
              largest_block_number = n_
            END IF
          END DO

          IF (largest_block == 0) THEN
            EXIT
          END IF

          CALL add_block_to_proc(Procs(p_), BlocksCollection(largest_block_number))
        END DO
      END DO

    END IF

    ! Update each proc's communication cost.
    CALL update_neighbor_procs(Procs)
    CALL check_communication_cost(Procs)

    ! Write the total weight, number of procs, ideal weight per proc.
    WRITE(*,*)
    WRITE(*,*), '      Weight  ', 'Ideal Weight/proc'
    WRITE(*,*), LOAD_S, optimal
    WRITE(*,*)

    ! Write out the weight on each proc.
    WRITE(*, *), '     proc #', '     nBlocks', '         LPP', &
    '          % REAL LOAD', '              % IDEAL LOAD'

    DO n_ = 1, mpi_nprocs
      LPP = Procs(n_)%weight + Procs(n_)%comm

      ! Total load for all processors is the sum of the load on each.
      TLPP = TLPP + LPP
    END DO

    DO n_ = 1, mpi_nprocs
      ! Load per processor is the sum of the weight and communication cost.
      LPP = Procs(n_)%weight + Procs(n_)%comm
      
      ! Real load is the load on this proc divided by the total load for all.
      REAL_LOAD = DFLOAT(LPP)/DFLOAT(TLPP)

      ! Ideal load is the load on this proc divided by serial load divided
      ! by the number of processors.
      IDEAL_LOAD = DFLOAT(Procs(n_)%weight)/DFLOAT(LOAD_S)

      WRITE(*, *), n_-1, Procs(n_)%nBlocks, LPP, REAL_LOAD, IDEAL_LOAD 

      IF (REAL_LOAD > PRED_REAL_SPEED) THEN
        PRED_REAL_SPEED = REAL_LOAD
      END IF
      IF (IDEAL_LOAD > PRED_IDEAL_SPEED) THEN
        PRED_IDEAL_SPEED = IDEAL_LOAD
      END IF
    END DO

    WRITE(*,*)
    WRITE(*,*), "Ideal                   ", DFLOAT(mpi_nprocs)
    WRITE(*,*), "Predicted Ideal Speedup ", 1.d0/PRED_IDEAL_SPEED
    WRITE(*,*), "Predicted Real Speedup  ", 1.d0/PRED_REAL_SPEED

    ! If we screwed this up we need to terminated execution.
    DO n_ = 1, nBlocks
      IF (BlocksCollection(n_)%SIZE > 0) THEN
        WRITE(*,*), n_
        WRITE(*,*), "Sorry, something went terribly wrong."
        WRITE(*,*), "Program exiting."
        STOP
      END IF
    END DO

  END SUBROUTINE
END MODULE
