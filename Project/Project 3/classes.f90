! Various constants we need in a lot of places and a routine to set
! the size of the grid.
module constants
  implicit none
  real(kind=8), parameter :: CFL = 1.d0
  real(kind=8), parameter :: k = 18.8d0, rho = 8000.d0, c_p = 500.d0
  real(kind=8), parameter :: pi = 3.141592654d0, rot = 30.d0*pi/180.d0
  real(kind=8), parameter :: alpha = k / (c_p * rho)
  integer :: IMAX, JMAX
  integer :: N, M ! Number of blocks.
  integer :: iBlockSize, jBlockSize
  integer :: nBlocks

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
    nBlocks = M * N
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
    real(kind=8) :: Ayi, Axi, Ayj, Axj
    real(kind=8) :: Ayi_half, Axi_half, Ayj_half, Axj_half
  end type GridPoint

contains
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
!        b%iLength = iBlockSize
!        b%jLength = jBlockSize

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
        else
          b%northFace%neighborBlock = 0
          b%northFace%neighborProc = 0
        end if

        if (b%eastFace%BC == -1) then
          b%eastFace%neighborBlock = n_ + 1
          b%eastFace%neighborProc = proc
        else
          b%eastFace%neighborBlock = 0
          b%eastFace%neighborProc = 0
        end if

        if (b%southFace%BC == -1) then
          b%southFace%neighborBlock = n_ - N
          b%southFace%neighborProc = proc
        else
          b%southFace%neighborBlock = 0
          b%southFace%neighborProc = 0
        end if

        if (b%westFace%BC == -1) then
          b%westFace%neighborBlock = n_ - 1
          b%westFace%neighborProc = proc
        else
          b%westFace%neighborBlock = 0
          b%westFace%neighborProc = 0
        end if

        ! Set corner neighbors and procs
        ! North East Corner
        if (b%northFace%BC == nB) then
          b%NECorner%BC = nBound
          b%NECorner%neighborBlock = 0
          b%NECorner%neighborProc = 0
        else if (b%eastFace%BC == eB) then
          b%NECorner%BC = eBound
          b%NECorner%neighborBlock = 0
          b%NECorner%neighborProc = 0
        else
          b%NECorner%BC = -1
          b%NECorner%neighborBlock = n_ + N + 1
          b%NECorner%neighborProc = proc
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
        else
          b%SECorner%BC = -1
          b%SECorner%neighborBlock = n_ - N + 1
          b%SECorner%neighborProc = proc
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
        else
          b%SWCorner%BC = -1
          b%SWCorner%neighborBlock = n_ - N - 1
          b%SWCorner%neighborProc = proc
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

  subroutine set_fluxes(Blocks)
    type (BlockType), target :: Blocks(:)
    type (GridCell), pointer :: c
    type (GridPoint), pointer :: p(:,:)
    integer :: i, j, n_

    do n_=1, nBlocks
      p => Blocks(n_)%Points
      do j = 0, jBlockSize
        do i = 0, iBlockSize
          c => Blocks(n_)%Cells(i,j)
          ! And these are the numbers that actually appear in the equations,
          ! saved here to (hopefully) save a moment or two during iteration.
          c%yPP = (  p(i,j)%Ayi_half + p(i,j)%Ayj_half )
          c%yNP = ( -p(i,j)%Ayi_half + p(i,j)%Ayj_half )
          c%yNN = ( -p(i,j)%Ayi_half - p(i,j)%Ayj_half )
          c%yPN = (  p(i,j)%Ayi_half - p(i,j)%Ayj_half )

          c%xNN = ( -p(i,j)%Axi_half - p(i,j)%Axj_half )
          c%xPN = (  p(i,j)%Axi_half - p(i,j)%Axj_half )
          c%xPP = (  p(i,j)%Axi_half + p(i,j)%Axj_half )
          c%xNP = ( -p(i,j)%Axi_half + p(i,j)%Axj_half )
        end do
      end do
    end do
  end subroutine

  subroutine set_constants(Blocks)
    type (BlockType), target :: Blocks(:)
    type (GridCell), pointer :: Cells(:,:)
    type (GridPoint), pointer :: Points(:,:)
    integer :: i, j, n_

    ! Calculate timesteps and assign secondary volumes.
    do n_ = 1, nBlocks
      do j = 1, jBlockSize
        do i = 1, iBlockSize
          Points => Blocks(n_)%Points
          Cells => Blocks(n_)%Cells
          ! Calculate the timestep using the CFL method described in class.
          Points(i, j)%timestep = ( ( CFL * 0.5d0 ) / alpha ) * 4.d0 * Points(i,j)%Vol2 ** 2 / &
            ( ( Points(i+1, j)%xp - Points(i, j)%xp )**2 + &
              ( Points(i, j+1)%yp - Points(i, j)%yp )**2 )

          ! Calculate this constant now so we don't recalculate in the solver loop.
          Points(i, j)%const = ( Points(i, j)%timestep * alpha / Points(i, j)%Vol2 )
        end do
      end do
    end do
  end subroutine

  ! Set the locations of each grid point.
  subroutine initialize_points(Blocks)
    type (BlockType), target :: Blocks(:)
    type (BlockType), pointer :: b
    type (GridPoint), pointer :: p
    integer :: i, j, n_

    do n_ = 1, nBlocks
      b=> Blocks(n_)
      do j = 0, jBlockSize+1
        do i = 0, iBlockSize+1
          p => Blocks(n_)%Points(i, j)

          p%xp = cos( 0.5d0 * pi * dfloat(IMAX - (i + b%lowI - 1)) / dfloat(IMAX-1))
          p%yp = cos( 0.5d0 * pi * dfloat(JMAX - (j + b%lowJ - 1)) / dfloat(JMAX-1))

!          p%x = p%xp * cos( rot ) + ( 1.d0 - p%yp ) * sin( rot )
!          p%y = p%yp * cos( rot ) + ( p%xp ) * sin( rot )
        end do
      end do
    end do
  end subroutine initialize_points

  subroutine set_bounds(Blocks)
    type (BlockType), target :: Blocks(:)
    type (BlockType), pointer :: b
    type (GridPoint), pointer :: p1, p2
    integer :: i, j, n_, neighbor

    do n_= 1, nBlocks
      b => Blocks(n_)

      if (b%northFace%BC == -1) then
        do i = 1, iBlockSize
          neighbor = b%northFace%neighborBlock
          p1 => b%Points(i, jBlockSize+1)
          p2 => Blocks(neighbor)%Points(i, 2)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        end do

        b%localJMAX = jBlockSize
      else
        b%localJMAX = jBlockSize - 1
      end if

      if (b%eastFace%BC == -1) then
        do j = 1, jBlockSize
          neighbor = b%eastFace%neighborBlock
          p1 => b%Points(iBlockSize+1, j)
          p2 => Blocks(neighbor)%Points(2, j)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        end do

        b%localIMAX = iBlockSize
      else
        b%localIMAX = iBlockSize - 1
      end if

      if (b%southFace%BC == -1) then
        do i = 1, iBlockSize
          neighbor = b%southFace%neighborBlock
          p1 => b%Points(i, 0)
          p2 => Blocks(neighbor)%Points(i, jBlockSize - 1)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        end do

        b%localIMIN = 0
      else
        b%localIMIN = 1
      end if

      if (b%westFace%BC == -1) then
        do j = 1, jBlockSize
          neighbor = b%westFace%neighborBlock
          p1 => b%Points(0, j)
          p2 => Blocks(neighbor)%Points(iBlockSize - 1, j)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        end do

        b%localJMIN = 0
      else
        b%localJMIN = 1
      end if

      ! Set corner points.
      if (b%NECorner%BC == -1) then
        neighbor = b%NECorner%neighborBlock
        p1 => b%Points(iBlockSize+1, jBlockSize+1)
        p2 => Blocks(neighbor)%Points(2, 2)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      if (b%SECorner%BC == -1) then
        neighbor = b%SECorner%neighborBlock
        p1 => b%Points(iBlockSize+1, 0)
        p2 => Blocks(neighbor)%Points(2, jBlockSize-1)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      if (b%SWCorner%BC == -1) then
        neighbor = b%SWCorner%neighborBlock
        p1 => b%Points(0, 0)
        p2 => Blocks(neighbor)%Points(iBlockSize-1, jBlockSize-1)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      if (b%NWCorner%BC == -1) then
        neighbor = b%NWCorner%neighborBlock
        p1 => b%Points(0, jBlockSize+1)
        p2 => Blocks(neighbor)%Points(iBlockSize-1, 2)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if
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
  end subroutine

  subroutine initialize_grid(b)
    type (BlockType) :: b(:)

    call create_blocks(b)
    call initialize_block_grid(b)
    call initialize_block_temp(b)

    write(*,*),'Initialized Grid'
  end subroutine

  subroutine initialize_faces_and_volumes(Blocks)
    type (BlockType), target :: Blocks(:)
    type (GridCell), pointer :: Cells(:,:)
    type (GridPoint), pointer :: p(:,:)
    integer :: i, j, n_

    do n_=1, nBlocks
      p => Blocks(n_)%Points
      Cells => Blocks(n_)%Cells

      do j = 0, jBlockSize
        do i = 0, iBlockSize + 1
          p(i,j)%Ayi = p(i,j+1)%y - p(i,j)%y
          p(i,j)%Axi = p(i,j+1)%x - p(i,j)%x
        end do
      end do

      do j = 0, jBlockSize+1
        do i = 0, iBlockSize
          p(i,j)%Ayj = p(i+1,j)%y - p(i,j)%y
          p(i,j)%Axj = p(i+1,j)%x - p(i,j)%x
        end do
      end do

      do j = 0, jBlockSize
        do i = 0, iBlockSize
          ! Calculate the volume of each cell.
          Cells(i, j)%V = abs( &
            p(i,    j)%x * p(i,  j+1)%y - p(i,    j)%y * p(i,  j+1)%x + &
            p(i,  j+1)%x * p(i+1,j+1)%y - p(i,  j+1)%y * p(i+1,j+1)%x + &
            p(i+1,j+1)%x * p(i+1,  j)%y - p(i+1,j+1)%y * p(i+1,  j)%x + &
            p(i+1,  j)%x * p(i,    j)%y - p(i+1,  j)%y * p(i,    j)%x )
          !            write(*,*) n_, i, j, Cells(i,j)%V

          p(i, j)%Vol2 = ( Cells(i, j)%V + Cells(i + 1, j)%V + &
            Cells(i, j + 1)%V + Cells(i + 1, j + 1)%V ) * 0.25d0

          p(i,j)%Ayi_half = ( p(i+1,j)%Ayi + p(i,j)%Ayi ) * 0.25d0
          p(i,j)%Axi_half = ( p(i+1,j)%Axi + p(i,j)%Axi ) * 0.25d0

          p(i,j)%Ayj_half = ( p(i,j+1)%Ayj + p(i,j)%Ayj ) * 0.25d0
          p(i,j)%Axj_half = ( p(i,j+1)%Axj + p(i,j)%Axj ) * 0.25d0
        end do
      end do
    end do
  end subroutine
end module BlockModule

module UpdateTemperature
  use GridPointModule
  use GridCellModule
  use BlockModule

  implicit none
  public

contains
  subroutine derivatives(Blocks)
    type (BlockType), target :: Blocks(:)
    type (BlockType), pointer :: MyBlock
    type (GridPoint), pointer :: p(:,:)
    type (GridCell), pointer :: c(:,:)
    real(kind=8) :: dTdx, dTdy
    integer :: i, j, n_

    ! Trapezoidal counter-clockwise integration to get the first
    ! derivatives in the x/y directions at the cell-center using
    ! Gauss's theorem.

    do n_ = 1, nBlocks
      MyBlock => Blocks(n_)
      p => MyBlock%Points
      c => MyBlock%Cells

      ! Reset the change in temperature to zero before we begin summing again.
      p%tempT = 0.d0

      do j = MyBlock%localJMIN, MyBlock%localJMAX
        do i =  MyBlock%localIMIN, MyBlock%localIMAX
          dTdx = + 0.5d0 * &
            ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * p(i+1, j)%Ayi - &
            ( p(i,   j)%T + p(i,  j+1)%T ) * p(i,   j)%Ayi - &
            ( p(i, j+1)%T + p(i+1,j+1)%T ) * p(i, j+1)%Ayj + &
            ( p(i,   j)%T + p(i+1,  j)%T ) * p(i,   j)%Ayj   &
            ) / c(i, j)%V

          dTdy = - 0.5d0 * &
            ( ( p(i+1, j)%T + p(i+1,j+1)%T ) * p(i+1, j)%Axi - &
            ( p(i,   j)%T + p(i,  j+1)%T ) * p(i,   j)%Axi - &
            ( p(i, j+1)%T + p(i+1,j+1)%T ) * p(i, j+1)%Axj + &
            ( p(i,   j)%T + p(i+1,  j)%T ) * p(i,   j)%Axj   &
            ) / c(i ,j)%V
!           write(*,*), i, j, dTdx, dTdy
!           write(*,*), i, j, c(i ,j)%V

          ! Alternate distributive scheme second-derivative operator.
          ! Updates the second derivative by adding the first times a constant
          ! during each time step.
          ! Pass out x and y second derivatives contributions.
          p(i+1,  j)%tempT = p(i+1,  j)%tempT + p(i+1,  j)%const * ( c(i, j)%yNN * dTdx + c(i, j)%xPP * dTdy )
          p(i,    j)%tempT = p(i,    j)%tempT + p(i,    j)%const * ( c(i, j)%yPN * dTdx + c(i, j)%xNP * dTdy )
          p(i,  j+1)%tempT = p(i,  j+1)%tempT + p(i,  j+1)%const * ( c(i, j)%yPP * dTdx + c(i, j)%xNN * dTdy )
          p(i+1,j+1)%tempT = p(i+1,j+1)%tempT + p(i+1,j+1)%const * ( c(i, j)%yNP * dTdx + c(i, j)%xPN * dTdy )
!                  write(*,*), i, j, p(i,j)%tempT

          ! Update temperatures.
          p(i,j)%T = p(i,j)%T + p(i,j)%tempT
        end do
      end do
    end do
  end subroutine

  subroutine update_ghosts(Blocks)
    type (BlockType), target :: Blocks(:)
    type (BlockType), pointer :: b
    integer :: n_, i, j

    do n_ = 1, nBlocks
      b => Blocks(n_)

      ! North face ghost nodes
      if (b%northFace%BC == -1) then
        ! Internal boundary
        do i = 1, iBlockSize
          b%Points(i, jBlockSize+1)%T = Blocks(b%northFace%neighborBlock)%Points(i, 2)%T
        end do
      else
        ! Reset to derichlet
        do i = 1, iBlockSize
          b%Points(i, jBlockSize)%T = 5.d0 * (sin(pi * b%Points(i, jBlockSize)%xp) + 1.d0)
        end do
      end if

      ! East face ghost nodes
      if (b%eastFace%BC == -1) then
        ! Internal boundary
        do j = 1, jBlockSize
          b%Points(iBlockSize+1, j)%T = Blocks(b%eastFace%neighborBlock)%Points(2, j)%T
        end do
      else
        ! Reset to derichlet
        do j = 1, jBlockSize
          b%Points(iBlockSize, j)%T = 3.d0 * b%Points(iBlockSize, j)%yp + 2.d0
        end do
      end if

      ! South face ghost nodes
      if (b%southFace%BC == -1) then
        ! Internal boundary
        do i = 1, iBlockSize
          b%Points(i, 0)%T = Blocks(b%southFace%neighborBlock)%Points(i, jBlockSize-1)%T
        end do
      else
        ! Reset to derichlet
        do i = 1, iBlockSize
          b%Points(i, 1)%T = abs(cos(pi * b%Points(i,1)%xp)) + 1.d0
        end do
      end if

      ! West face ghost nodes
      if (b%westFace%BC == -1) then
        ! Internal boundary
        do j = 1, jBlockSize
          b%Points(0, j)%T = Blocks(b%westFace%neighborBlock)%Points(iBlockSize-1, j)%T
        end do
      else
        ! Reset to derichlet
        do j = 1, jBlockSize
          b%Points(1, j)%T = 3.d0 * b%Points(1, j)%yp + 2.d0
        end do
      end if

      ! Corners
      ! North east corner
      if (b%NECorner%BC == -1) then
        b%Points(iBlockSize+1,jBlockSize+1)%T = Blocks(b%NECorner%neighborBlock)%Points(2,2)%T
      end if

      ! South east corner
      if (b%SECorner%BC == -1) then
        b%Points(iBlockSize+1,0)%T = Blocks(b%SECorner%neighborBlock)%Points(2,jBlockSize-1)%T
      end if

      ! South west corner
      if (b%SWCorner%BC == -1) then
        b%Points(0,0)%T = Blocks(b%SWCorner%neighborBlock)%Points(iBlockSize-1,jBlockSize-1)%T
      end if

      ! North west corner
      if (b%NWCorner%BC == -1) then
        b%Points(0,jBlockSize+1)%T = Blocks(b%NWCorner%neighborBlock)%Points(iBlockSize-1,2)%T
      end if
    end do
  end subroutine
end module UpdateTemperature
