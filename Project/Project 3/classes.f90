! Various constants we need in a lot of places and a routine to set
! the size of the grid.
module constants
  implicit none
  real(kind=8), parameter :: CFL = 1.1d0
  real(kind=8), parameter :: k = 18.8d0, rho = 8000.d0, c_p = 500.d0
  real(kind=8), parameter :: pi = 3.141592654d0, rot = 30.d0*pi/180.d0
  real(kind=8), parameter :: alpha = k / (c_p * rho)
  integer :: IMAX, JMAX, N, M ! Number of blocks.
  integer :: iBlockSize, jBlockSize, nBlocks
  integer :: nB = 1, eB = 2, sB = 3, wB = 4
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

module BlockModule
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

  type GridCell
    real(kind=8) :: V
    real(kind=8) :: yPP, yNP, yNN, yPN
    real(kind=8) :: xNN, xPN, xPP, xNP
  end type GridCell

  type Face
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
    type (Face) :: northFace, southFace, eastFace, westFace
    type (Face) :: NECorner, SECorner, SWCorner, NWCorner
  end type BlockType

contains
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
          Points(i, j)%timestep = ( ( CFL * 0.5d0 ) / alpha ) * Points(i,j)%Vol2 ** 2 / &
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

        b%localJMIN = 0
      else
        b%localJMIN = 1
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

        b%localIMIN = 0
      else
        b%localIMIN = 1
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
          Cells(i, j)%V = abs(( p(i+1, j)%xp - p(i, j)%xp) * &
                               ( p(i, j+1)%yp - p(i, j)%yp))
        end do
      end do

      do j = 0, jBlockSize
        do i = 0, iBlockSize
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
