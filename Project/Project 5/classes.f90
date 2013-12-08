! Various constants we need in a lot of places, a routine to set
! the size of the grid, and and a routine to set number of blocks.
module constants
  implicit none

  ! We use MPI in many of our routines.
  include "mpif.h"

  integer, parameter :: IMAX = 101
  integer, parameter :: JMAX = 101
  integer, parameter :: N = 10
  integer, parameter :: M = 10
  integer, parameter :: iBlockSize = 1 + (IMAX - 1) / N
  integer, parameter :: jBlockSize = 1 + (JMAX - 1) / M
  integer, parameter :: nBlocks = M * N

  real(kind=8), parameter :: CFL = 1.09d0
  integer, parameter :: max_steps = 100000
  
  real(kind=8), parameter :: k = 18.8d0, rho = 8000.d0, c_p = 500.d0
  real(kind=8), parameter :: pi = 3.141592654d0, rot = 30.d0*pi/180.d0
  real(kind=8), parameter :: alpha = k / (c_p * rho)
  
  integer :: nB = 1, eB = 2, sB = 3, wB = 4
  integer :: INTERNAL_BOUNDARY = -1, EXTERNAL_BOUNDARY = -2, PROC_BOUNDARY = -3
  integer :: step = 0

  ! MPI related variables.
  integer :: MyID, MyNBlocks
  integer :: ierror, mpi_nprocs, request
  integer :: status(MPI_STATUS_SIZE)
end module

! Prof's clock module.
module clock
!   use constants
  real(kind=8) :: start_time, end_time, wall_time

contains
  subroutine timestamp ( )
    character ( len = 8 ) ampm
    integer ( kind = 4 ) d
    character ( len = 8 ) date
    integer ( kind = 4 ) h
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mm
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
      'January  ', 'February ', 'March    ', 'April    ', &
      'May      ', 'June     ', 'July     ', 'August   ', &
      'September', 'October  ', 'November ', 'December ' /)
    integer ( kind = 4 ) n
    integer ( kind = 4 ) s
    character ( len = 10 ) time
    integer ( kind = 4 ) values(8)
    integer ( kind = 4 ) y
    character ( len = 5 ) zone

    call date_and_time ( date, time, zone, values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if ( h < 12 ) then
      ampm = 'AM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Noon'
      else
        ampm = 'PM'
      end if
    else
      h = h - 12
      if ( h < 12 ) then
        ampm = 'PM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Midnight'
        else
          ampm = 'AM'
        end if
      end if
    end if

    write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
      trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )
    return
  end

  subroutine start_clock()
    ! call system time to determine flow solver wall time
    start_time = MPI_Wtime()
    write(*,*) "Start time: "
    call timestamp()
  end subroutine start_clock

  subroutine end_clock()
    ! determine total wall time for solver
    end_time = MPI_Wtime()
    write(*,*) "End time: "
    call timestamp()
    wall_time = end_time - start_time
    write(*,*) "Wall time: ", wall_time
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
    integer :: BC, neighborBlock, neighborLocalBlock, neighborProc
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
    type (BlockType) :: Blocks(nBlocks)
  end type Proc

  type LinkedList
      type(LinkedList),pointer :: next
      integer :: id
  end type LinkedList

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
      if (b%northFace%BC == INTERNAL_BOUNDARY) then
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
      if (b%eastFace%BC == INTERNAL_BOUNDARY) then
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
      if (b%southFace%BC == INTERNAL_BOUNDARY) then
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
      if (b%westFace%BC == INTERNAL_BOUNDARY) then
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
      if (b%NECorner%BC == INTERNAL_BOUNDARY) then
        neighbor = b%NECorner%neighborBlock
        p1 => b%Points(iBlockSize+1, jBlockSize+1)
        p2 => Blocks(neighbor)%Points(2, 2)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      ! South east corner
      if (b%SECorner%BC == INTERNAL_BOUNDARY) then
        neighbor = b%SECorner%neighborBlock
        p1 => b%Points(iBlockSize+1, 0)
        p2 => Blocks(neighbor)%Points(2, jBlockSize-1)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      ! South west corner
      if (b%SWCorner%BC == INTERNAL_BOUNDARY) then
        neighbor = b%SWCorner%neighborBlock
        p1 => b%Points(0, 0)
        p2 => Blocks(neighbor)%Points(iBlockSize-1, jBlockSize-1)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      end if

      ! North west corner
      if (b%NWCorner%BC == INTERNAL_BOUNDARY) then
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
    type (GridPoint), pointer :: p1, p2, p3, p4
    integer :: i, j, n_

    do n_=1, MyNBlocks

      ! Calculate fluxes.
      do j = 0, jBlockSize
        do i = 0, iBlockSize + 1
          p1 => Blocks(n_)%Points(i,j)
          p2 => Blocks(n_)%Points(i,j+1)
          p1%Ayi = p2%y - p1%y
          p1%Axi = p2%x - p1%x
        end do
      end do

      do j = 0, jBlockSize+1
        do i = 0, iBlockSize
          p1 => Blocks(n_)%Points(i,j)
          p2 => Blocks(n_)%Points(i+1,j)
          p1%Ayj = p2%y - p1%y
          p1%Axj = p2%x - p1%x
        end do
      end do

      ! Calculate the volumes.
      do j = 0, jBlockSize
        do i = 0, iBlockSize
          p1 => Blocks(n_)%Points(i,j)
          p2 => Blocks(n_)%Points(i+1,j)
          p3 => Blocks(n_)%Points(i,j+1)
          p1%V = abs(( p2%xp - p1%xp) * &
                     ( p3%yp - p1%yp))
        end do
      end do

      ! Calculate secondary volumes.
      do j = 0, jBlockSize
        do i = 0, iBlockSize
          p1 => Blocks(n_)%Points(i,    j)
          p2 => Blocks(n_)%Points(i+1,  j)
          p3 => Blocks(n_)%Points(i,  j+1)
          p4 => Blocks(n_)%Points(i+1,j+1)
          p1%Vol2 = p1%Vol2 + p1%V * 0.25d0
          p2%Vol2 = p2%Vol2 + p2%V * 0.25d0
          p3%Vol2 = p3%Vol2 + p3%V * 0.25d0
          p4%Vol2 = p4%Vol2 + p4%V * 0.25d0
        end do
      end do
    end do
  end subroutine

  subroutine set_constants(Blocks)
    type (BlockType), target :: Blocks(:)
    type (GridPoint), pointer :: Points(:,:)
    type (GridPoint), pointer :: p0, p1, p2, p3, p4
    integer :: i, j, n_
    real(kind=8) :: timestep
    real(kind=8) :: temp

    ! Constants used during iteration.
    do n_=1, MyNBlocks
      do j = 0, jBlockSize
        do i = 0, iBlockSize
          p1 => Blocks(n_)%Points(i,j)
          p2 => Blocks(n_)%Points(i+1,j)
          p3 => Blocks(n_)%Points(i,j+1)
          ! These are the numbers that actually appear in the equations,
          ! saved here to save a moment or two during iteration.
          p1%yPP = (  ( p2%Ayi + p1%Ayi ) + ( p3%Ayj + p1%Ayj ) ) * 0.25d0
          p1%yNP = ( -( p2%Ayi + p1%Ayi ) + ( p3%Ayj + p1%Ayj ) ) * 0.25d0
          p1%yNN = ( -( p2%Ayi + p1%Ayi ) - ( p3%Ayj + p1%Ayj ) ) * 0.25d0
          p1%yPN = (  ( p2%Ayi + p1%Ayi ) - ( p3%Ayj + p1%Ayj ) ) * 0.25d0

          p1%xNN = ( -( p2%Axi + p1%Axi ) - ( p3%Axj + p1%Axj ) ) * 0.25d0
          p1%xPN = (  ( p2%Axi + p1%Axi ) - ( p3%Axj + p1%Axj ) ) * 0.25d0
          p1%xPP = (  ( p2%Axi + p1%Axi ) + ( p3%Axj + p1%Axj ) ) * 0.25d0
          p1%xNP = ( -( p2%Axi + p1%Axi ) + ( p3%Axj + p1%Axj ) ) * 0.25d0
        end do
      end do
    end do

    ! Calculate timesteps and assign secondary volumes.
    do n_ = 1, MyNBlocks
      do j = 1, jBlockSize
        do i = 1, iBlockSize
          Points => Blocks(n_)%Points
          ! Calculate the timestep using the CFL method described in class.
          p0 => Blocks(n_)%Points(i,   j)
          p1 => Blocks(n_)%Points(i+1, j)
          p2 => Blocks(n_)%Points(i-1, j)
          p3 => Blocks(n_)%Points(i, j+1)
          p4 => Blocks(n_)%Points(i, j-1)

          temp = ( ( p1%xp - p2%xp )**2 + ( p3%yp - p4%yp )**2 )
          if (temp > 0) then 
            timestep = ( ( CFL * 2.d0 ) / alpha ) * p0%Vol2 ** 2 / temp

            ! Calculate this constant now so we don't recalculate in the solver loop.
            p0%const = ( timestep * alpha / p0%Vol2 )
          end if

        end do
      end do
    end do
  end subroutine
end module BlockModule
