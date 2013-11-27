! Various constants we need in a lot of places, a routine to set
! the size of the grid, and and a routine to set number of blocks.
module constants
  implicit none

  integer, parameter :: IMAX = 101
  integer, parameter :: JMAX = 101
  integer, parameter :: M = 10
  integer, parameter :: N = 10
  integer, parameter :: nProcs = 4
  integer, parameter :: iBlockSize = 1 + (IMAX - 1) / N
  integer, parameter :: jBlockSize = 1 + (JMAX - 1) / M
  integer, parameter :: nBlocks = M * N

  real(kind=8), parameter :: CFL = 1.09d0
  real(kind=8), parameter :: k = 18.8d0, rho = 8000.d0, c_p = 500.d0
  real(kind=8), parameter :: pi = 3.141592654d0, rot = 30.d0*pi/180.d0
  real(kind=8), parameter :: alpha = k / (c_p * rho)
  integer :: nB = 1, eB = 2, sB = 3, wB = 4
  integer :: Internal = -1
  integer :: step = 0

  integer :: MyID, MyNBlocks
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

! Contains derived data types and initialization routines.
module BlockModule
  use constants
  implicit none
  public

  type GridPoint
    real(kind=8) :: x, xp, y, yp
    real(kind=8) :: T, tempT
    real(kind=8) :: const
    real(kind=8) :: Ayi, Axi, Ayj, Axj
    real(kind=8) :: V, Vol2
    real(kind=8) :: yPP, yNP, yNN, yPN
    real(kind=8) :: xNN, xPN, xPP, xNP
  end type GridPoint

  type Neighbor
    integer :: BC, neighborBlock, neighborProc
  end type

  type BlockType
    type (GridPoint) :: Points(0:iBlockSize + 1,0:jBlockSize + 1)
    integer :: id, proc, size
    integer :: lowJ, lowI, lowITemp, lowJTemp
    integer :: localJMIN, localIMIN, localJMAX, localIMAX
    type (Neighbor) :: northFace, southFace, eastFace, westFace
    type (Neighbor) :: NECorner, SECorner, SWCorner, NWCorner
  end type BlockType

  type Proc
    integer :: procID, weight, comm, nBlocks
    ! Number of blocks on each proc should be roughly nBlocks/nProcs.
    type (BlockType) :: Blocks(nBlocks/nProcs + 5)
  end type Proc

contains
  ! Set the true bounds and ghost nodes for each block.
  subroutine set_bounds(Blocks)
    type (BlockType), target :: Blocks(:)
    type (BlockType), pointer :: b
    type (GridPoint), pointer :: p1, p2
    integer :: i, j, n_, neighbor

    do n_= 1, nBlocks
      b => Blocks(n_)

      ! Initialize ghost nodes. If block face is internal
      ! also set different bounds for the solver loop.
      ! North face.
      if (b%northFace%BC == -1) then
        do i = 1, iBlockSize
          neighbor = b%northFace%neighborBlock
          p1 => b%Points(i, jBlockSize+1)
          p2 => Blocks(neighbor)%Points(i, 2)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        end do
      end if

      ! East face.
      if (b%eastFace%BC == -1) then
        do j = 1, jBlockSize
          neighbor = b%eastFace%neighborBlock
          p1 => b%Points(iBlockSize+1, j)
          p2 => Blocks(neighbor)%Points(2, j)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        end do
      end if

      ! South face.
      if (b%southFace%BC == -1) then
        do i = 1, iBlockSize
          neighbor = b%southFace%neighborBlock
          p1 => b%Points(i, 0)
          p2 => Blocks(neighbor)%Points(i, jBlockSize - 1)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        end do
      end if

      ! West face.
      if (b%westFace%BC == -1) then
        do j = 1, jBlockSize
          neighbor = b%westFace%neighborBlock
          p1 => b%Points(0, j)
          p2 => Blocks(neighbor)%Points(iBlockSize - 1, j)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        end do
      end if

      ! Set corner points.
      ! North east corner
      if (b%NECorner%BC == -1) then
        neighbor = b%NECorner%neighborBlock
        p1 => b%Points(iBlockSize+1, jBlockSize+1)
        p2 => Blocks(neighbor)%Points(2, 2)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      ! South east corner
      if (b%SECorner%BC == -1) then
        neighbor = b%SECorner%neighborBlock
        p1 => b%Points(iBlockSize+1, 0)
        p2 => Blocks(neighbor)%Points(2, jBlockSize-1)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      ! South west corner
      if (b%SWCorner%BC == -1) then
        neighbor = b%SWCorner%neighborBlock
        p1 => b%Points(0, 0)
        p2 => Blocks(neighbor)%Points(iBlockSize-1, jBlockSize-1)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      ! North west corner
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

  ! Set the prime locations of each grid point.
  subroutine initialize_points(Blocks)
    type (BlockType), target :: Blocks(:)
    type (BlockType), pointer :: b
    type (GridPoint), pointer :: p
    integer :: i, j, n_

    do n_ = 1, MyNBlocks
      ! Set our lower bound for updating the temperature so
      ! we don't update along the edge.
      Blocks(n_)%lowITemp = Blocks(n_)%localIMIN
      Blocks(n_)%lowJTemp = Blocks(n_)%localJMIN

      if (Blocks(n_)%localIMIN == 1 .and. Blocks(n_)%localJMIN == 1) then
        Blocks(n_)%lowITemp = Blocks(n_)%localIMIN+1
        Blocks(n_)%lowJTemp = Blocks(n_)%localJMIN+1
      else if (Blocks(n_)%localIMIN == 1) then
        Blocks(n_)%lowITemp = Blocks(n_)%localIMIN+1
      else if (Blocks(n_)%localJMIN == 1) then
        Blocks(n_)%lowJTemp = Blocks(n_)%localJMIN+1
      end if
    end do

    do n_ = 1, MyNBlocks
      b=> Blocks(n_)
      do j = 0, jBlockSize+1
        do i = 0, iBlockSize+1
          p => Blocks(n_)%Points(i, j)

          ! Have to convert from i, j to global i, j.
          p%xp = cos( 0.5d0 * pi * dfloat(IMAX - (i + b%lowI - 1)) / dfloat(IMAX-1))
          p%yp = cos( 0.5d0 * pi * dfloat(JMAX - (j + b%lowJ - 1)) / dfloat(JMAX-1))
        end do
      end do
    end do
  end subroutine initialize_points

  subroutine initialize_faces_and_volumes(Blocks)
    type (BlockType), target :: Blocks(:)
    type (GridPoint), pointer :: p(:,:)
    integer :: i, j, n_

    do n_=1, MyNBlocks
      p => Blocks(n_)%Points

      ! Calculate fluxes.
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

      ! Calculate the volumes.
      do j = 0, jBlockSize
        do i = 0, iBlockSize
          p(i, j)%V = abs(( p(i+1, j)%xp - p(i, j)%xp) * &
                          ( p(i, j+1)%yp - p(i, j)%yp))
        end do
      end do

      ! Calculate secondary volumes.
      do j = 0, jBlockSize
        do i = 0, iBlockSize
          p(i,    j)%Vol2 = p(i,    j)%Vol2 + p(i,    j)%V * 0.25d0
          p(i+1,  j)%Vol2 = p(i+1,  j)%Vol2 + p(i+1,  j)%V * 0.25d0
          p(i,  j+1)%Vol2 = p(i,  j+1)%Vol2 + p(i,  j+1)%V * 0.25d0
          p(i+1,j+1)%Vol2 = p(i+1,j+1)%Vol2 + p(i+1,j+1)%V * 0.25d0
        end do
      end do

    end do
  end subroutine

  subroutine set_constants(Blocks)
    type (BlockType), target :: Blocks(:)
    type (GridPoint), pointer :: Points(:,:)
    type (GridPoint), pointer :: p
    integer :: i, j, n_
    real(kind=8) :: timestep
    real(kind=8) :: Ayi_half, Axi_half, Ayj_half, Axj_half

    ! Calculate fluxes.
    Ayi_half(i,j) = ( Points(i+1,j)%Ayi + Points(i,j)%Ayi ) * 0.25d0
    Axi_half(i,j) = ( Points(i+1,j)%Axi + Points(i,j)%Axi ) * 0.25d0
    Ayj_half(i,j) = ( Points(i,j+1)%Ayj + Points(i,j)%Ayj ) * 0.25d0
    Axj_half(i,j) = ( Points(i,j+1)%Axj + Points(i,j)%Axj ) * 0.25d0

    ! Constants used during iteration.
    do n_=1, MyNBlocks
      Points => Blocks(n_)%Points
      do j = 0, jBlockSize
        do i = 0, iBlockSize
          p => Points(i,j)
          ! These are the numbers that actually appear in the equations,
          ! saved here to save a moment or two during iteration.
          p%yPP = (  Ayi_half(i,j) + Ayj_half(i,j) )
          p%yNP = ( -Ayi_half(i,j) + Ayj_half(i,j) )
          p%yNN = ( -Ayi_half(i,j) - Ayj_half(i,j) )
          p%yPN = (  Ayi_half(i,j) - Ayj_half(i,j) )

          p%xNN = ( -Axi_half(i,j) - Axj_half(i,j) )
          p%xPN = (  Axi_half(i,j) - Axj_half(i,j) )
          p%xPP = (  Axi_half(i,j) + Axj_half(i,j) )
          p%xNP = ( -Axi_half(i,j) + Axj_half(i,j) )
        end do
      end do
    end do

    ! Calculate timesteps and assign secondary volumes.
    do n_ = 1, MyNBlocks
      do j = 1, jBlockSize
        do i = 1, iBlockSize
          Points => Blocks(n_)%Points
          ! Calculate the timestep using the CFL method described in class.

          timestep = ( ( CFL * 2.d0 ) / alpha ) * Points(i,j)%Vol2 ** 2 / &
                       ( ( Points(i+1, j)%xp - Points(i-1, j)%xp )**2 + &
                         ( Points(i, j+1)%yp - Points(i, j-1)%yp )**2 )


          ! Calculate this constant now so we don't recalculate in the solver loop.
          Points(i, j)%const = ( timestep * alpha / Points(i, j)%Vol2 )
        end do
      end do
    end do
  end subroutine
end module BlockModule
