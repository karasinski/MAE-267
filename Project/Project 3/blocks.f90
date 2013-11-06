! This module contains all the subroutines used to
! create our blocks and initial grid/temp files.

module GridCreation
  use BlockModule

  contains
  subroutine create_blocks(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:)
    type (BlockType), pointer :: b

    integer :: n_ = 1, proc = 1, type = 1, iN, iM
    integer :: nBound, sBound, eBound, wBound

    do iM = 1, M
      do iN = 1, N
        b => BlocksCollection(n_)

        b%type = type
        b%proc = proc
        b%highI = 1 + iN * (iBlockSize - 1)
        b%lowI = b%highI - (iBlockSize - 1)
        b%highJ = 1 + iM * (jBlockSize - 1)
        b%lowJ = b%highJ - (jBlockSize - 1)

        ! Default values.
        b%northFace%neighborBlock = 0
        b%northFace%neighborProc = 0
        b%southFace%neighborBlock = 0
        b%southFace%neighborProc = 0
        b%eastFace%neighborBlock = 0
        b%eastFace%neighborProc = 0
        b%westFace%neighborBlock = 0
        b%westFace%neighborProc = 0

        b%NECorner%neighborBlock = 0
        b%NECorner%neighborProc = 0
        b%SECorner%neighborBlock = 0
        b%SECorner%neighborProc = 0
        b%SWCorner%neighborBlock = 0
        b%SWCorner%neighborProc = 0
        b%NWCorner%neighborBlock = 0
        b%NWCorner%neighborProc = 0

        if (b%highJ == JMAX) then
          nBound = nB
        else
          nBound = -1
        end if
        if (b%highI == IMAX) then
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
        b%northFace%BC = nBound
        b%southFace%BC = sBound
        b%eastFace%BC = eBound
        b%westFace%BC = wBound

        if (b%northFace%BC == -1) then
          b%northFace%neighborBlock = n_ + N
          b%northFace%neighborProc = proc
        end if

        if (b%eastFace%BC == -1) then
          b%eastFace%neighborBlock = n_ + 1
          b%eastFace%neighborProc = proc
        end if

        if (b%southFace%BC == -1) then
          b%southFace%neighborBlock = n_ - N
          b%southFace%neighborProc = proc
        end if

        if (b%westFace%BC == -1) then
          b%westFace%neighborBlock = n_ - 1
          b%westFace%neighborProc = proc
        end if

        ! Set corner neighbors and procs
        ! North East Corner
        if (b%northFace%BC == nB) then
          b%NECorner%BC = nBound
        else if (b%eastFace%BC == eB) then
          b%NECorner%BC = eBound
        else
          b%NECorner%BC = -1
          b%NECorner%neighborBlock = n_ + N + 1
          b%NECorner%neighborProc = proc
        end if

        ! South East Corner
        if (b%southFace%BC == sB) then
          b%SECorner%BC = sBound
        else if (b%eastFace%BC == eB) then
          b%SECorner%BC = eBound
        else
          b%SECorner%BC = -1
          b%SECorner%neighborBlock = n_ - N + 1
          b%SECorner%neighborProc = proc
        end if

        ! South West Corner
        if (b%southFace%BC == sB) then
          b%SWCorner%BC = sBound
        else if (b%westFace%BC == wB) then
          b%SWCorner%BC = wBound
        else
          b%SWCorner%BC = -1
          b%SWCorner%neighborBlock = n_ - N - 1
          b%SWCorner%neighborProc = proc
        end if

        ! North West Corner
        if (b%northFace%BC == nB) then
          b%NWCorner%BC = nBound
        else if (b%westFace%BC == wB) then
          b%NWCorner%BC = wBound
        else
          b%NWCorner%BC = -1
          b%NWCorner%neighborBlock = n_ + N - 1
          b%NWCorner%neighborProc = proc
        end if

        n_ = n_ + 1
      end do
    end do
  end subroutine

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

          p%xp = cos( 0.5d0 * pi * dfloat(IMAX - (b%lowI + (i - 1 ))) / dfloat(IMAX-1))
          p%yp = cos( 0.5d0 * pi * dfloat(JMAX - (b%lowJ + (j - 1 ))) / dfloat(JMAX-1))

          p%x = p%xp * cos( rot ) + ( 1.d0 - p%yp ) * sin( rot )
          p%y = p%yp * cos( rot ) + p%xp * sin( rot )
        end do
      end do
    end do
  end subroutine

  subroutine initialize_block_temp(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:)
    type (BlockType), pointer :: b
    type (GridPoint), pointer :: p(:,:)
    integer :: i, j, n_
    real(kind=8) :: T_0 = 3.5d0

    do n_ = 1, nBlocks
      b => BlocksCollection(n_)
      p => BlocksCollection(n_)%Points

      ! Initialize all points to T0.
      p(1:iBlockSize, 1:jBlockSize)%T = T_0

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
end module
