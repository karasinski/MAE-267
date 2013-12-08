MODULE constants
  IMPLICIT NONE

  ! We use MPI in many of our routines.
  INCLUDE "mpif.h"

  ! Various constants we need in a lot of places.
  INTEGER, PARAMETER :: IMAX = 501
  INTEGER, PARAMETER :: JMAX = 501
  INTEGER, PARAMETER :: N = 10
  INTEGER, PARAMETER :: M = 10
  INTEGER, PARAMETER :: iBlockSize = 1 + (IMAX - 1) / N
  INTEGER, PARAMETER :: jBlockSize = 1 + (JMAX - 1) / M
  INTEGER, PARAMETER :: nBlocks = M * N

  REAL(KIND=8), PARAMETER :: CFL = 1.0d0
  INTEGER, PARAMETER :: max_steps = 100000

  REAL(KIND=8), PARAMETER :: k = 18.8d0, rho = 8000.d0, c_p = 500.d0
  REAL(KIND=8), PARAMETER :: pi = 3.141592654d0, rot = 30.d0*pi/180.d0
  REAL(KIND=8), PARAMETER :: alpha = k / (c_p * rho)

  INTEGER, PARAMETER :: nB = 1, eB = 2, sB = 3, wB = 4
  INTEGER, PARAMETER :: INTERNAL_BOUNDARY = -1, EXTERNAL_BOUNDARY = -2, PROC_BOUNDARY = -3
  INTEGER :: step = 0

  ! MPI related variables.
  INTEGER :: MyID, MyNBlocks
  INTEGER :: ierror, mpi_nprocs, request
  INTEGER :: STATUS(MPI_STATUS_SIZE)
END MODULE

! Contains derived data types and initialization routines.
MODULE BlockModule
  USE constants
  IMPLICIT NONE
  PUBLIC

  TYPE GridPoint
    REAL(KIND=8) :: x, xp, y, yp
    REAL(KIND=8) :: T, tempT
    REAL(KIND=8) :: const
    REAL(KIND=8) :: Ayi, Axi, Ayj, Axj
    REAL(KIND=8) :: V, Vol2
    REAL(KIND=8) :: yPP, yNP, yNN, yPN
    REAL(KIND=8) :: xNN, xPN, xPP, xNP
  END TYPE GridPoint

  TYPE Neighbor
    INTEGER :: BC, neighborBlock, neighborLocalBlock, neighborProc
  END TYPE

  TYPE BlockType
    TYPE (GridPoint) :: Points(0:iBlockSize + 1,0:jBlockSize + 1)
    INTEGER :: ID, proc, SIZE
    INTEGER :: lowJ, lowI, lowITemp, lowJTemp
    INTEGER :: localJMIN, localIMIN, localJMAX, localIMAX
    TYPE (Neighbor) :: northFace, southFace, eastFace, westFace
    TYPE (Neighbor) :: NECorner, SECorner, SWCorner, NWCorner
  END TYPE BlockType

  TYPE Proc
    INTEGER :: procID, weight, comm, nBlocks
    ! Number of blocks on each proc should be roughly nBlocks/nProcs.
    TYPE (BlockType) :: Blocks(nBlocks)
  END TYPE Proc

  TYPE LinkedList
    TYPE (LinkedList), POINTER :: next
    INTEGER :: ID
  END TYPE LinkedList

CONTAINS
  ! Set the true bounds and ghost nodes for each block.
  SUBROUTINE set_bounds(Blocks)
    TYPE (BlockType), TARGET :: Blocks(:)
    TYPE (BlockType), POINTER :: b
    TYPE (GridPoint), POINTER :: p1, p2
    INTEGER :: i, j, n_, neighbor

    DO n_= 1, nBlocks
      b => Blocks(n_)

      ! Initialize ghost nodes. If block face is internal
      ! also set different bounds for the solver loop.
      ! North face.
      IF (b%northFace%BC == INTERNAL_BOUNDARY) THEN
        DO i = 1, iBlockSize
          neighbor = b%northFace%neighborBlock
          p1 => b%Points(i, jBlockSize+1)
          p2 => Blocks(neighbor)%Points(i, 2)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        END DO
      END IF

      ! South face.
      IF (b%southFace%BC == INTERNAL_BOUNDARY) THEN
        DO i = 1, iBlockSize
          neighbor = b%southFace%neighborBlock
          p1 => b%Points(i, 0)
          p2 => Blocks(neighbor)%Points(i, jBlockSize - 1)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        END DO
      END IF

      ! East face.
      IF (b%eastFace%BC == INTERNAL_BOUNDARY) THEN
        DO j = 1, jBlockSize
          neighbor = b%eastFace%neighborBlock
          p1 => b%Points(iBlockSize+1, j)
          p2 => Blocks(neighbor)%Points(2, j)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        END DO
      END IF

      ! West face.
      IF (b%westFace%BC == INTERNAL_BOUNDARY) THEN
        DO j = 1, jBlockSize
          neighbor = b%westFace%neighborBlock
          p1 => b%Points(0, j)
          p2 => Blocks(neighbor)%Points(iBlockSize - 1, j)

          p1%x = p2%x
          p1%y = p2%y
          p1%T = p2%T
        END DO
      END IF

      ! Set corner points.
      ! North east corner
      IF (b%NECorner%BC == INTERNAL_BOUNDARY) THEN
        neighbor = b%NECorner%neighborBlock
        p1 => b%Points(iBlockSize+1, jBlockSize+1)
        p2 => Blocks(neighbor)%Points(2, 2)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      END IF

      ! South east corner
      IF (b%SECorner%BC == INTERNAL_BOUNDARY) THEN
        neighbor = b%SECorner%neighborBlock
        p1 => b%Points(iBlockSize+1, 0)
        p2 => Blocks(neighbor)%Points(2, jBlockSize-1)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      END IF

      ! South west corner
      IF (b%SWCorner%BC == INTERNAL_BOUNDARY) THEN
        neighbor = b%SWCorner%neighborBlock
        p1 => b%Points(0, 0)
        p2 => Blocks(neighbor)%Points(iBlockSize-1, jBlockSize-1)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      END IF

      ! North west corner
      IF (b%NWCorner%BC == INTERNAL_BOUNDARY) THEN
        neighbor = b%NWCorner%neighborBlock
        p1 => b%Points(0, jBlockSize+1)
        p2 => Blocks(neighbor)%Points(iBlockSize-1, 2)

        p1%x = p2%x
        p1%y = p2%y
        p1%T = p2%T
      END IF
    END DO
  END SUBROUTINE

  ! Set the prime locations of each grid point.
  SUBROUTINE initialize_points(Blocks)
    TYPE (BlockType), TARGET :: Blocks(:)
    TYPE (BlockType), POINTER :: b
    TYPE (GridPoint), POINTER :: p
    INTEGER :: i, j, n_

    DO n_ = 1, MyNBlocks
      ! Set our lower bound for updating the temperature so
      ! we don't update along the edge.
      Blocks(n_)%lowITemp = Blocks(n_)%localIMIN
      Blocks(n_)%lowJTemp = Blocks(n_)%localJMIN

      IF (Blocks(n_)%localIMIN == 1 .and. Blocks(n_)%localJMIN == 1) THEN
        Blocks(n_)%lowITemp = Blocks(n_)%localIMIN+1
        Blocks(n_)%lowJTemp = Blocks(n_)%localJMIN+1
      ELSE IF (Blocks(n_)%localIMIN == 1) THEN
        Blocks(n_)%lowITemp = Blocks(n_)%localIMIN+1
      ELSE IF (Blocks(n_)%localJMIN == 1) THEN
        Blocks(n_)%lowJTemp = Blocks(n_)%localJMIN+1
      END IF
    END DO

    DO n_ = 1, MyNBlocks
      b=> Blocks(n_)
      DO j = 0, jBlockSize+1
        DO i = 0, iBlockSize+1
          p => Blocks(n_)%Points(i, j)

          ! Have to convert from i, j to global i, j.
          p%xp = COS( 0.5d0 * pi * DFLOAT(IMAX - (i + b%lowI - 1)) / DFLOAT(IMAX-1))
          p%yp = COS( 0.5d0 * pi * DFLOAT(JMAX - (j + b%lowJ - 1)) / DFLOAT(JMAX-1))
        END DO
      END DO
    END DO
  END SUBROUTINE initialize_points

  SUBROUTINE initialize_faces_and_volumes(Blocks)
    TYPE (BlockType), TARGET :: Blocks(:)
    TYPE (GridPoint), POINTER :: p1, p2, p3, p4
    INTEGER :: i, j, n_

    DO n_= 1, MyNBlocks

      ! Calculate fluxes.
      ! i Direction
      DO j = 0, jBlockSize
        DO i = 0, iBlockSize + 1
          p1 => Blocks(n_)%Points(i,j)
          p2 => Blocks(n_)%Points(i,j+1)
          p1%Ayi = p2%y - p1%y
          p1%Axi = p2%x - p1%x
        END DO
      END DO

      ! j Direction
      DO j = 0, jBlockSize+1
        DO i = 0, iBlockSize
          p1 => Blocks(n_)%Points(i,j)
          p2 => Blocks(n_)%Points(i+1,j)
          p1%Ayj = p2%y - p1%y
          p1%Axj = p2%x - p1%x
        END DO
      END DO

      ! Calculate the volumes.
      DO j = 0, jBlockSize
        DO i = 0, iBlockSize
          p1 => Blocks(n_)%Points(i,j)
          p2 => Blocks(n_)%Points(i+1,j)
          p3 => Blocks(n_)%Points(i,j+1)
          p1%V = ABS(( p2%xp - p1%xp) * &
            ( p3%yp - p1%yp))
        END DO
      END DO

      ! Calculate secondary volumes.
      DO j = 0, jBlockSize
        DO i = 0, iBlockSize
          p1 => Blocks(n_)%Points(i,    j)
          p2 => Blocks(n_)%Points(i+1,  j)
          p3 => Blocks(n_)%Points(i,  j+1)
          p4 => Blocks(n_)%Points(i+1,j+1)
          p1%Vol2 = p1%Vol2 + p1%V * 0.25d0
          p2%Vol2 = p2%Vol2 + p2%V * 0.25d0
          p3%Vol2 = p3%Vol2 + p3%V * 0.25d0
          p4%Vol2 = p4%Vol2 + p4%V * 0.25d0
        END DO
      END DO
    END DO
  END SUBROUTINE

  SUBROUTINE set_constants(Blocks)
    TYPE (BlockType), TARGET :: Blocks(:)
    TYPE (GridPoint), POINTER :: Points(:,:)
    TYPE (GridPoint), POINTER :: p0, p1, p2, p3, p4
    INTEGER :: i, j, n_
    REAL(KIND=8) :: timestep
    REAL(KIND=8) :: temp

    ! Constants used during iteration.
    DO n_=1, MyNBlocks
      DO j = 0, jBlockSize
        DO i = 0, iBlockSize
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
        END DO
      END DO
    END DO

    ! Calculate timesteps and assign secondary volumes.
    DO n_ = 1, MyNBlocks
      DO j = 1, jBlockSize
        DO i = 1, iBlockSize
          Points => Blocks(n_)%Points

          ! Calculate the timestep using the CFL method described in class.
          p0 => Blocks(n_)%Points(i,   j)
          p1 => Blocks(n_)%Points(i+1, j)
          p2 => Blocks(n_)%Points(i-1, j)
          p3 => Blocks(n_)%Points(i, j+1)
          p4 => Blocks(n_)%Points(i, j-1)

          temp = ( ( p1%xp - p2%xp )**2 + ( p3%yp - p4%yp )**2 )
          IF (temp > 0) THEN
            timestep = ( ( CFL * 2.d0 ) / alpha ) * p0%Vol2 ** 2 / temp

            ! Calculate this constant now so we don't recalculate in the solver loop.
            p0%const = ( timestep * alpha / p0%Vol2 )
          END IF

        END DO
      END DO
    END DO
  END SUBROUTINE
END MODULE BlockModule

MODULE clock
  REAL(KIND=8) :: start_time, end_time, wall_time

CONTAINS
  SUBROUTINE timestamp()
    CHARACTER (LEN = 9), PARAMETER, DIMENSION(12) :: month = (/ &
      'January  ', 'February ', 'March    ', 'April    ', &
      'May      ', 'June     ', 'July     ', 'August   ', &
      'September', 'October  ', 'November ', 'December ' /)

    CHARACTER (LEN = 10) :: time
    CHARACTER (LEN = 8) :: ampm, DATE
    CHARACTER (LEN = 5) :: zone
    INTEGER :: d, h, m, mm, n, s, y
    INTEGER :: values(8)

    CALL DATE_AND_TIME(DATE, time, zone, values)

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    IF (h < 12) THEN
      ampm = 'AM'
    ELSE IF ( h == 12 ) THEN
      IF (n == 0 .and. s == 0) THEN
        ampm = 'Noon'
      ELSE
        ampm = 'PM'
      END IF
    ELSE
      h = h - 12
      IF ( h < 12 ) THEN
        ampm = 'PM'
      ELSE IF ( h == 12 ) THEN
        IF (n == 0 .and. s == 0) THEN
          ampm = 'Midnight'
        ELSE
          ampm = 'AM'
        END IF
      END IF
    END IF

    WRITE ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
      TRIM ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, TRIM ( ampm )
    RETURN
  END

  SUBROUTINE start_clock()
    ! call system time to determine flow solver wall time
    start_time = MPI_Wtime()
    WRITE(*,*) "Start time: ", start_time
    CALL timestamp()
  END SUBROUTINE start_clock

  SUBROUTINE end_clock()
    ! determine total wall time for solver
    end_time = MPI_Wtime()
    WRITE(*,*) "End time: ", end_time
    CALL timestamp()
    wall_time = (end_time - start_time) / MPI_Wtick()
    WRITE(*,*) "Wall time: ", wall_time
  END SUBROUTINE end_clock
END MODULE
