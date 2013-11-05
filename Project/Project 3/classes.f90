! Various constants we need in a lot of places and a routine to set
! the size of the grid.
module constants
  implicit none
  real(kind=8), parameter :: CFL = 1.14d0
  real(kind=8), parameter :: k = 18.8d0, rho = 8000.d0, c_p = 500.d0
  real(kind=8), parameter :: pi = 3.141592654d0, rot = 30.d0*pi/180.d0
  real(kind=8), parameter :: alpha = k / (c_p * rho)
  integer :: IMAX, JMAX
  integer :: N, M ! Number of blocks.
  integer :: iBlockSize, jBlockSize

  integer :: nB = 1
  integer :: eB = 2
  integer :: sB = 3
  integer :: wB = 4
contains
  subroutine SetGridSize(length)
    integer :: length
    IMAX = length
    JMAX = length
  end subroutine SetGridSize

  subroutine SetNumberOfBlocks(m_, n_)
    integer :: n_, m_
    M = m_
    N = n_
    iBlockSize = 1 + (IMAX - 1) / N
    jBlockSize = 1 + (JMAX - 1) / M
  end subroutine SetNumberOfBlocks
end module

! Prof's clock module.
module clock
  integer clock_start,clock_end,clock_max,clock_rate
  real(kind=8) wall_time

contains
  subroutine start_clock()
    ! call system time to determine flow solver wall time
    call system_clock(count_max=clock_max,count_rate=clock_rate)
    call system_clock(clock_start)
  end subroutine start_clock

  subroutine end_clock()
    ! determine total wall time for solver
    call system_clock(clock_end)
    wall_time=dfloat(clock_end-clock_start)/dfloat(clock_rate)
    print*,'solver wall clock time (seconds)',wall_time
  end subroutine end_clock
end module

! GridPointModule contains the GridPoint type and routines
! to initialize and set boundary conditions.
module GridPointModule
  use constants
  implicit none

  public
  type GridPoint
    integer :: i, j
    real(kind=8) :: x, xp, y, yp
    real(kind=8) :: T, tempT
    real(kind=8) :: timestep, Vol2, const
  end type GridPoint

contains
  ! Set the initial locations and initial temperature of
  ! each grid point.
!  subroutine initialize_points(Points)
!    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
!    type (GridPoint), pointer :: p
!    integer :: i, j
!
!    do j = 1, JMAX
!      do i = 1, IMAX
!        p => Points(i, j)
!
!        p%i = i
!        p%j = j
!
!        p%xp = cos( 0.5d0 * pi * dfloat(IMAX-i) / dfloat(IMAX-1) )
!        p%yp = cos( 0.5d0 * pi * dfloat(JMAX-j) / dfloat(JMAX-1) )
!
!        p%x = p%xp * cos( rot ) + ( 1.d0 - p%yp ) * sin( rot )
!        p%y = p%yp * cos( rot ) + ( p%xp ) * sin( rot )
!
!        p%T = 3.5d0
!      end do
!    end do
!  end subroutine initialize_points

  ! This is used to set the boundary conditions.
  subroutine set_temperature(p, T)
    type (GridPoint), intent(inout) :: p
    real(kind=8) :: T

    p%T = T
  end subroutine set_temperature
end module GridPointModule

module GridCellModule
  use GridPointModule
  implicit none

  public
  type GridCell
    real(kind=8) :: V
    real(kind=8) :: yPP, yNP, yNN, yPN
    real(kind=8) :: xNN, xPN, xPP, xNP
  end type GridCell
end module GridCellModule

module BlockModule
  use GridPointModule
  use GridCellModule

  implicit none
  public

  type face_type
    integer :: BC, neighborBlock, neighborProc
  end type

  type BlockType
    ! This sucks.
    type (GridPoint) :: Points(0:103,0:103)
    type (GridCell)  :: Cells(0:102,0:102)
    integer :: iStart, jStart, iBound, jBound

    integer :: type, proc, lowJ, highJ, lowI, highI
    integer :: iLength, jLength
    integer :: localJMIN, localIMIN, localJMAX, localIMAX
    type (face_type) :: northFace, southFace, eastFace, westFace
    type (face_type) :: NECorner, SECorner, SWCorner, NWCorner
  end type BlockType

contains

!  subroutine initialize_blocks(Blocks, Points, Cells)
!    type (BlockType) :: Blocks(:,:)
!    type (GridPoint) :: Points(:,:)
!    type (GridCell)  :: Cells(:,:)
!    integer :: m_, n_, i_, j_, i, j, iBound, jBound
!
!    ! Size of each block.
!    iBound = 1 + (IMAX - 1) / N
!    jBound = 1 + (JMAX - 1) / M
!
!    mloop: do m_ = 1, M
!      nloop: do n_ = 1, N
!        j = 0
!        jloop: do j_ = 1 +(m_ - 1) * ( jBound ), m_ * ( jBound ) + 1
!          i = 0
!          j = j + 1
!          iloop: do i_ = 1 + (n_ - 1) * ( iBound ), n_ * ( iBound ) + 1
!            i = i + 1
!
!            ! If we're passed the number of points continue.
!            if ( (i_ > IMAX .or. j_ > JMAX) .or. (i_ <= 0 .or. j_ <= 0) ) then
!              continue
!!              Blocks(m_, n_)%Points(i, j)%T = 666.d0
!            else
!              ! Hand out points to the blocks.
!!              write(*,*), m_,n_,j_,i_
!              Blocks(m_, n_)%Points(i, j) = Points(i_, j_)
!              Blocks(m_, n_)%iBound = i
!              Blocks(m_, n_)%jBound = j
!            end if
!
!            ! Similarly...
!            ! If we're passed the number of cells continue.
!            if ( (i_ > IMAX-1 .or. j_ > JMAX-1) .or. (i_ <= 0 .or. j_ <= 0) ) then
!              continue
!            else
!              ! Hand out cells to the blocks.
!              Blocks(m_, n_)%Cells(i, j) = Cells(i_, j_)
!            end if
!          end do iloop
!        end do jloop
!      end do nloop
!    end do mloop
!
!    write(*, *), "          m_          ", "n_         ", "iBound      ", "jBound"
!    do m_ = 1, M
!      do n_ = 1, N
!        write(*, *), m_, n_, Blocks(m_,n_)%iBound, Blocks(m_,n_)%jBound
!      end do
!    end do
!!    call set_block_bounds(Blocks)
!  end subroutine initialize_blocks
!
!  subroutine set_block_bounds(Blocks)
!    type (BlockType) :: Blocks(:,:)
!    integer :: m_, n_
!
!    write(*, *), "          m_          ", "n_         ", "iBound      ", "jBound"
!    do m_ = 1, M
!      do n_ = 1, N
!        Blocks(m_,n_)%iBound = (maxval(Blocks(m_,n_)%Points(:,:)%i) - &
!                                minval(Blocks(m_,n_)%Points(:,:)%i,   &
!                                MASK = Blocks(m_,n_)%Points(:,:)%i>0))
!
!        Blocks(m_,n_)%jBound = (maxval(Blocks(m_,n_)%Points(:,:)%j) - &
!                                minval(Blocks(m_,n_)%Points(:,:)%j,   &
!                                MASK = Blocks(m_,n_)%Points(:,:)%j>0))
!        write(*, *), m_, n_, Blocks(m_,n_)%iBound, Blocks(m_,n_)%jBound
!      end do
!    end do
!
!!    write(*,*)
!  end subroutine set_block_bounds

  subroutine create_blocks(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:,:)
    type (BlockType), pointer :: b

    integer :: i = 1, proc = 1, type = 1, iN, iM
    integer :: nBound, sBound, eBound, wBound

    i = 1
    do iM = 1, M
      do iN = 1, N
        ! Defaults
        b => BlocksCollection(iM, iN)

        b%type = type
        b%proc = proc
        b%highI = 1 + iN * (iBlockSize - 1)
        b%lowI = b%highI - (iBlockSize - 1)
        b%highJ = 1 + iM * (jBlockSize - 1)
        b%lowJ = b%highJ - (jBlockSize - 1)
        b%iLength = iBlockSize
        b%jLength = jBlockSize

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
          b%northFace%neighborBlock = i + N
          b%northFace%neighborProc = proc
        else
          b%northFace%neighborBlock = 0
          b%northFace%neighborProc = 0
        end if

        write(*,*), b%northFace%neighborBlock

        if (b%highI == IMAX) then
          b%eastFace%neighborBlock = i + 1
          b%eastFace%neighborProc = proc
        else
          b%eastFace%neighborBlock = 0
          b%eastFace%neighborProc = 0
        end if

        if (b%lowJ == 1) then
          b%southFace%neighborBlock = i - N
          b%southFace%neighborProc = proc
        else
          b%southFace%neighborBlock = 0
          b%southFace%neighborProc = 0
        end if

        if (b%lowI == 1) then
          b%westFace%neighborBlock = i - 1
          b%westFace%neighborProc = proc
        else
          b%westFace%neighborBlock = 0
          b%westFace%neighborProc = 0
        end if

        ! Corners
        if (b%northFace%BC == nB) then
          b%NECorner%BC = nBound
          b%NECorner%neighborBlock = 0
          b%NECorner%neighborProc = 0
        else if (b%northFace%BC == eB) then
          b%NECorner%BC = eBound
          b%NECorner%neighborBlock = 0
          b%NECorner%neighborProc = 0
        else
          b%NECorner%BC = -1
          b%NECorner%neighborBlock = i + N + 1
          b%NECorner%neighborProc = proc
        end if

        if (b%southFace%BC == sB) then
          b%SECorner%BC = sBound
          b%SECorner%neighborBlock = 0
          b%SECorner%neighborProc = 0
        else if (b%southFace%BC == eB) then
          b%SECorner%BC = eBound
          b%SECorner%neighborBlock = 0
          b%SECorner%neighborProc = 0
        else
          b%SECorner%BC = -1
          b%SECorner%neighborBlock = i - N + 1
          b%SECorner%neighborProc = proc
        end if

        if (b%southFace%BC == sB) then
          b%SWCorner%BC = sBound
          b%SWCorner%neighborBlock = 0
          b%SWCorner%neighborProc = 0
        else if (b%southFace%BC == wB) then
          b%SWCorner%BC = wBound
          b%SWCorner%neighborBlock = 0
          b%SWCorner%neighborProc = 0
        else
          b%SWCorner%BC = -1
          b%SWCorner%neighborBlock = i - N - 1
          b%SWCorner%neighborProc = proc
        end if

        if (b%northFace%BC == nB) then
          b%NWCorner%BC = nBound
          b%NWCorner%neighborBlock = 0
          b%NWCorner%neighborProc = 0
        else if (b%northFace%BC == wB) then
          b%NWCorner%BC = wBound
          b%NWCorner%neighborBlock = 0
          b%NWCorner%neighborProc = 0
        else
          b%NWCorner%BC = -1
          b%NWCorner%neighborBlock = i - N - 1
          b%NWCorner%neighborProc = proc
        end if

        i = i + 1
      end do
    end do
  end subroutine

  subroutine initialize_block_grid(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:,:)
    type (GridPoint), pointer :: p
    type (BlockType), pointer :: b
    integer :: i, j, iM, iN

    do iM = 1, M
      do iN = 1, N
        do j = 1, jBlockSize
          do i = 1, iBlockSize

            b => BlocksCollection(iM, iN)
            p => BlocksCollection(iM, iN)%Points(i, j)

!            p%i = b%lowI + i - 1
!            p%j = b%lowJ + j - 1

            p%xp = cos( 0.5d0 * pi * dfloat(IMAX - (b%lowI + (i - 1 ))) / dfloat(IMAX-1))
            p%yp = cos( 0.5d0 * pi * dfloat(JMAX - (b%lowJ + (j - 1 ))) / dfloat(JMAX-1))

            p%x = p%xp * cos( rot ) + ( 1.d0 - p%yp ) * sin( rot )
            p%y = p%yp * cos( rot ) + p%xp * sin( rot )

!            write(*, *), iM, iN, p%i, p%j

          end do
        end do
      end do
    end do
  end subroutine

  subroutine set_secondary_areas(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:,:)
    type (GridCell), pointer :: c
    type (GridPoint), pointer :: p(:,:)
    integer :: i, j,  iM, iN
    real(kind=8) :: Ayi, Axi, Ayj, Axj
    real(kind=8) :: Ayi_half, Axi_half, Ayj_half, Axj_half

    ! These 'area's are used to do trapezoidal counter-clockwise
    ! integration to get the first derivatives in the x/y directions
    ! at the cell-center using Gauss's theorem.
    Ayi(i,j) = ( p(i, j+1)%y - p(i, j)%y )
    Axi(i,j) = ( p(i, j+1)%x - p(i, j)%x )
    Ayj(i,j) = ( p(i+1, j)%y - p(i, j)%y )
    Axj(i,j) = ( p(i+1, j)%x - p(i, j)%x )

    ! These areas are used in the calculation of fluxes
    ! for the alternate distributive scheme second-
    ! derivative operator.
    Ayi_half(i,j) = ( Ayi(i+1, j) + Ayi(i, j)   ) * 0.25d0
    Axi_half(i,j) = ( Axi(i+1, j) + Axi(i, j)   ) * 0.25d0
    Ayj_half(i,j) = ( Ayj(i, j)   + Ayj(i, j+1) ) * 0.25d0
    Axj_half(i,j) = ( Axj(i, j)   + Axj(i, j+1) ) * 0.25d0

    do iM=1, M
      do iN=1, N
        do j = 1, jBlockSize-1
          do i = 1, iBlockSize-1
            p => BlocksCollection(iM, iN)%Points
            c => BlocksCollection(iM, iN)%Cells(i,j)
            ! And these are the numbers that actually appear in the equations,
            ! saved here to (hopefully) save a moment or two during iteration.
            c%yPP = (  Ayi_half(i,j) + Ayj_half(i,j) )
            c%yNP = ( -Ayi_half(i,j) + Ayj_half(i,j) )
            c%yNN = ( -Ayi_half(i,j) - Ayj_half(i,j) )
            c%yPN = (  Ayi_half(i,j) - Ayj_half(i,j) )

            c%xNN = ( -Axi_half(i,j) - Axj_half(i,j) )
            c%xPN = (  Axi_half(i,j) - Axj_half(i,j) )
            c%xPP = (  Axi_half(i,j) + Axj_half(i,j) )
            c%xNP = ( -Axi_half(i,j) + Axj_half(i,j) )
          end do
        end do
      end do
    end do
  end subroutine

  subroutine set_constants(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:,:)
    type (GridCell), pointer :: Cells(:,:)
    type (GridPoint), pointer :: Points(:,:)
    integer :: i, j, iM, iN

    ! Calculate timesteps and assign secondary volumes.
    do iM=1, M
      do iN=1, N
        do j = 1, jBlockSize - 1
          do i = 1, iBlockSize - 1
            Points => BlocksCollection(iM, iN)%Points
            Cells => BlocksCollection(iM, iN)%Cells
            ! Calculate the timestep using the CFL method described in class.
            Points(i, j)%timestep = ( ( CFL * 0.5d0 ) / alpha ) * Cells(i, j)%V ** 2 / &
              ( ( Points(i+1, j)%xp - Points(i, j)%xp )**2 + &
              ( Points(i, j+1)%yp - Points(i, j)%yp )**2 )

            ! Calculate the secondary volumes around each point. As we have rectangular points,
            ! these are simply the sum of the surrounding primary cells divied by four.
            Points(i, j)%Vol2 = ( Cells(i, j)%V + Cells(i - 1, j)%V + &
              Cells(i, j - 1)%V + Cells(i - 1, j - 1)%V ) * 0.25d0

            ! Calculate this constant now so we don't recalculate in the solver loop.
            Points(i, j)%const = ( Points(i, j)%timestep * alpha / Points(i, j)%Vol2 )
          end do
        end do
      end do
    end do
  end subroutine

  subroutine set_bounds(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:,:)
    type (GridCell), pointer :: c(:,:)
    type (GridPoint), pointer :: p1, p2
    integer, pointer :: neighbor
    integer :: i, j, iM, iN, neighborM, neighborN

    ! Calculate timesteps and assign secondary volumes.
    do iM=1, M
      do iN=1, N
!        do i = 1, iBlockSize
          p1 =>BlocksCollection(iM, iN)%Points(i, jBlockSize+1)
          p2 =>BlocksCollection(iM, iN)%Points(i, 2)

          neighbor => BlocksCollection(iM, iN)%northFace%neighborBlock
          neighborM = iM * M + iN - 1
          neighborN = mod(neighbor, N)
          ! Need to convert from neighbor to neighbor_M and neighbor_N

          write(*,*), iM, iN, neighbor, neighborM, neighborN
!          p%x =
!        end do

        do j = 1, jBlockSize

        end do
      end do
    end do
  end subroutine

  subroutine initialize_block_temp(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:,:)
    type (BlockType), pointer :: b
    type (GridPoint), pointer :: p(:,:)
    integer :: i, j, iM, iN
    real(kind=8) :: T_0 = 3.5d0

    do iM = 1, M
      do iN = 1, N
        b => BlocksCollection(iM, iN)
        p => BlocksCollection(iM, iN)%Points

        p(2:iBlockSize-1, 2:jBlockSize-1)%T = T_0

        if (b%northFace%BC == -1) then
          p(1:iBlockSize, jBlockSize)%T = T_0
        else if (b%northFace%BC == nB) then
          do i = 1, iBlockSize
            p(i, jBlockSize)%T = 5.d0 * (sin(pi * p(i, jBlockSize)%xp) + 1.d0)
          end do
        end if

        if (b%eastFace%BC == -1) then
          p(iBlockSize, 1:jBlockSize)%T = T_0
        else if (b%eastFace%BC == eB) then
          do j = 1, jBlockSize
            p(iBlockSize, j)%T = (3.d0 * p(iBlockSize, j)%yp) + 2.d0
          end do
        end if

        if (b%southFace%BC == -1) then
          p(1:iBlockSize, 1)%T = T_0
        else if (b%southFace%BC == sB) then
          do i = 1, iBlockSize
            p(i, 1)%T = abs(cos(pi * p(i,1)%xp)) + 1.d0
          end do
        end if

        if (b%westFace%BC == -1) then
          p(1, 1:jBlockSize)%T = T_0
        else if (b%westFace%BC == wB) then
          do j = 1, jBlockSize
            p(1, j)%T = (3.d0 * p(1, j)%yp) + 2.d0
          end do
        end if
      end do
    end do
  end subroutine

  subroutine initialize_grid(b)
    type (BlockType) :: b(:,:)

    call create_blocks(b)
    call initialize_block_grid(b)
    call initialize_block_temp(b)

    write(*,*),'Initialized Grid'
  end subroutine

  subroutine initialize_cells(BlocksCollection)
    type (BlockType), target :: BlocksCollection(:,:)
    type (GridCell), pointer :: Cells(:,:)
    type (GridPoint), pointer :: Points(:,:)
    integer :: i, j, iM, iN

    do iM=1, M
      do iN=1, N
        Points => BlocksCollection(iM, iN)%Points
        Cells => BlocksCollection(iM, iN)%Cells

        do j = 1, jBlockSize-1
          do i = 1, iBlockSize-1
            ! Calculate the volume of each cell.
            Cells(i, j)%V = ( Points(i+1, j)%xp - Points(i, j)%xp) * &
                            ( Points(i, j+1)%yp - Points(i, j)%yp)
          end do
        end do
      end do
    end do
  end subroutine initialize_cells

end module BlockModule

module UpdateTemperature
  use GridPointModule
  use GridCellModule
  use BlockModule

  implicit none
  public

contains
  subroutine derivatives(MyBlock)
    type (BlockType), target :: MyBlock
    type (GridPoint), pointer :: p(:,:)
    type (GridCell), pointer :: c(:,:)

    real(kind=8) :: Ayi, Axi, Ayj, Axj
    real(kind=8) :: dTdx, dTdy
    integer :: i, j

    ! Trapezoidal counter-clockwise integration to get the first
    ! derivatives in the x/y directions at the cell-center using
    ! Gauss's theorem.
    Ayi(i,j) = ( p(i, j+1)%y - p(i, j)%y )
    Axi(i,j) = ( p(i, j+1)%x - p(i, j)%x )
    Ayj(i,j) = ( p(i+1, j)%y - p(i, j)%y )
    Axj(i,j) = ( p(i+1, j)%x - p(i, j)%x )

    p => MyBlock%Points
    c => MyBlock%Cells

    ! Reset the change in temperature to zero before we begin summing again.
    p%tempT = 0.d0

    do j = 1, jBlockSize-1
      do i = 1, iBlockSize-1

        dTdx = + 0.5d0 * &
          ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * Ayi(i+1, j) - &
            ( p(i,   j)%T + p(i,  j+1)%T ) * Ayi(i,   j) - &
            ( p(i, j+1)%T + p(i+1,j+1)%T ) * Ayj(i, j+1) + &
            ( p(i,   j)%T + p(i+1,  j)%T ) * Ayj(i,   j)   &
          ) / c(i, j)%V

        dTdy = - 0.5d0 * &
          ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * Axi(i+1, j) - &
            ( p(i,   j)%T + p(i,  j+1)%T ) * Axi(i,   j) - &
            ( p(i, j+1)%T + p(i+1,j+1)%T ) * Axj(i, j+1) + &
            ( p(i,   j)%T + p(i+1,  j)%T ) * Axj(i,   j)   &
          ) / c(i ,j)%V

        ! Alternate distributive scheme second-derivative operator.
        ! Updates the second derivative by adding the first times a constant
        ! during each time step.
        ! Pass out x and y second derivatives contributions.
        p(i+1,  j)%tempT = p(i+1,  j)%tempT + p(i+1,  j)%const * ( c(i, j)%yNN * dTdx + c(i, j)%xPP * dTdy )
        p(i,    j)%tempT = p(i,    j)%tempT + p(i,    j)%const * ( c(i, j)%yPN * dTdx + c(i, j)%xNP * dTdy )
        p(i,  j+1)%tempT = p(i,  j+1)%tempT + p(i,  j+1)%const * ( c(i, j)%yPP * dTdx + c(i, j)%xNN * dTdy )
        p(i+1,j+1)%tempT = p(i+1,j+1)%tempT + p(i+1,j+1)%const * ( c(i, j)%yNP * dTdx + c(i, j)%xPN * dTdy )

        write(*,*), i, j, p(i,j)%tempT
      end do
    end do
  end subroutine
end module UpdateTemperature
