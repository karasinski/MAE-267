MODULE plot3D_module
  USE clock
  USE BlockModule
  IMPLICIT NONE

  INTEGER :: gridUnit  = 30   ! Unit for grid file
  INTEGER :: tempUnit = 21    ! Unit for temp file
  INTEGER :: configUnit = 99  ! Unit for config file
  REAL(KIND=8) :: tRef = 1.d0 ! tRef number
  REAL(KIND=8) :: dum = 0.d0  ! dummy values

CONTAINS
  SUBROUTINE write_configuration_file(Procs)
    TYPE (Proc), TARGET :: Procs(:)
    TYPE (BlockType), POINTER :: BlocksCollection(:)
    TYPE (BlockType), POINTER :: b
    INTEGER :: n_, p_
    CHARACTER(2) :: name1, str1, name2, str2

    WRITE( name1, '(i2)' )  mpi_nprocs
    READ( name1, * ) str1

    ! Make a directory for this run (disabled on vortex)
    ! call execute_command_line ('mkdir -p ' // str1 )

    10 FORMAT(3I5)
    20 FORMAT(33I5)

    DO p_ = 1, mpi_nprocs
      BlocksCollection => Procs(p_)%Blocks

      WRITE( name2, '(i2)' )  Procs(p_)%procID
      READ( name2, * ) str2

      OPEN(UNIT = configUnit, FILE = TRIM(TRIM(str1)//'/configuration_file.dat.p'//str2), FORM='formatted')
      WRITE(configUnit, 10) Procs(p_)%nBlocks, iBlockSize, jBlockSize
      DO n_ = 1, Procs(p_)%nBlocks
        b => BlocksCollection(n_)
        WRITE(configUnit, 20) b%ID, b%proc, b%SIZE, &
          b%lowI, b%lowJ, &
          b%localIMIN, b%localIMAX, b%localJMIN, b%localJMAX, &
          b%northFace%BC, b%northFace%neighborBlock, b%northFace%neighborLocalBlock, b%northFace%neighborProc, &
          b%southFace%BC, b%southFace%neighborBlock, b%southFace%neighborLocalBlock, b%southFace%neighborProc, &
          b%eastFace%BC,  b%eastFace%neighborBlock,  b%eastFace%neighborLocalBlock,  b%eastFace%neighborProc,  &
          b%westFace%BC,  b%westFace%neighborBlock,  b%westFace%neighborLocalBlock,  b%westFace%neighborProc, &
          b%NECorner%BC,  b%NECorner%neighborBlock,  b%NECorner%neighborLocalBlock,  b%NECorner%neighborProc, &
          b%NWCorner%BC,  b%NWCorner%neighborBlock,  b%NWCorner%neighborLocalBlock,  b%NWCorner%neighborProc, &
          b%SWCorner%BC,  b%SWCorner%neighborBlock,  b%SWCorner%neighborLocalBlock,  b%SWCorner%neighborProc, &
          b%SECorner%BC,  b%SECorner%neighborBlock,  b%SECorner%neighborLocalBlock,  b%SECorner%neighborProc
      END DO
      CLOSE(configUnit)
    END DO

  END SUBROUTINE

  SUBROUTINE read_configuration_file(Blocks)
    TYPE (BlockType), POINTER :: Blocks(:)
    TYPE (BlockType), POINTER :: b
    INTEGER :: n_, readnBlocks, readiBlockSize, readjBlockSize
    CHARACTER(2) :: name1, str1, name2, str2

    WRITE( name1, '(i2)' )  mpi_nprocs
    READ( name1, * ) str1
    WRITE( name2, '(i2)' )  MyID
    READ( name2, * ) str2
    OPEN(UNIT = 1, FILE = TRIM(TRIM(str1)//'/configuration_file.dat.p'//str2), STATUS='old')

    10 FORMAT(3I5)
    20 FORMAT(33I5)

    READ(1, 10) readnBlocks, readiBlockSize, readjBlockSize

    ! Set my nblocks.
    MyNBlocks = readnBlocks

    ! Allocate our blocks array.
    ALLOCATE(Blocks(1:MyNBlocks))

    DO n_ = 1, readnBlocks
      b => Blocks(n_)
      READ(1, 20) b%ID, b%proc, b%SIZE, &
        b%lowI, b%lowJ, &
        b%localIMIN, b%localIMAX, b%localJMIN, b%localJMAX, &
        b%northFace%BC, b%northFace%neighborBlock, b%northFace%neighborLocalBlock, b%northFace%neighborProc, &
        b%southFace%BC, b%southFace%neighborBlock, b%southFace%neighborLocalBlock, b%southFace%neighborProc, &
        b%eastFace%BC,  b%eastFace%neighborBlock,  b%eastFace%neighborLocalBlock,  b%eastFace%neighborProc,  &
        b%westFace%BC,  b%westFace%neighborBlock,  b%westFace%neighborLocalBlock,  b%westFace%neighborProc, &
        b%NECorner%BC,  b%NECorner%neighborBlock,  b%NECorner%neighborLocalBlock,  b%NECorner%neighborProc, &
        b%NWCorner%BC,  b%NWCorner%neighborBlock,  b%NWCorner%neighborLocalBlock,  b%NWCorner%neighborProc, &
        b%SWCorner%BC,  b%SWCorner%neighborBlock,  b%SWCorner%neighborLocalBlock,  b%SWCorner%neighborProc, &
        b%SECorner%BC,  b%SECorner%neighborBlock,  b%SECorner%neighborLocalBlock,  b%SECorner%neighborProc
    END DO

    CLOSE(1)
    !     write(*,*), 'Processor ', MyID, ' read configuration file.'
  END SUBROUTINE

  SUBROUTINE read_grid_file(Blocks)
    TYPE (BlockType) :: Blocks(:)
    INTEGER :: n_, i, j
    INTEGER :: readnBlocks, readiBlockSize, readjBlockSize
    CHARACTER(2) :: name1, str1, name2, str2

    ! Format statements
    10     FORMAT(I10)
    20     FORMAT(10I10)
    30     FORMAT(10E20.8)

    ! Open files
    WRITE( name1, '(i2)' )  mpi_nprocs
    READ( name1, * ) str1
    WRITE( name2, '(i2)' )  MyID
    READ( name2, * ) str2
    OPEN(UNIT=gridUnit,FILE=TRIM(TRIM(str1)//'/grid.dat.p'//str2), STATUS='old')

    ! Read grid file
    READ(gridUnit,10) readnBlocks
    READ(gridUnit,20) (readiBlockSize, readjBlockSize, i=1, readnBlocks)
    DO n_ = 1, readnBlocks
      READ(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                        ((Blocks(n_)%Points(i,j)%y,i=0,iBlockSize+1),j=0,jBlockSize+1)
    END DO

    ! Close file
    CLOSE(gridUnit)
    !     write(*,*), 'Processor ', MyID, ' read grid file.'
  END SUBROUTINE

  SUBROUTINE read_temp_file(Blocks)
    TYPE (BlockType) :: Blocks(:)
    INTEGER :: n_, i, j
    INTEGER :: readnBlocks, readiBlockSize, readjBlockSize
    CHARACTER(2) :: name1, str1, name2, str2

    ! Format statements
    10     FORMAT(I10)
    20     FORMAT(10I10)
    30     FORMAT(10E20.8)

    ! Open files
    WRITE( name1, '(i2)' )  mpi_nprocs
    READ( name1, * ) str1
    WRITE( name2, '(i2)' )  MyID
    READ( name2, * ) str2
    OPEN(UNIT=tempUnit,FILE=TRIM(TRIM(str1)//'/temp.dat.p'//str2), STATUS='old')

    ! Read temperature file
    READ(tempUnit,10) readnBlocks
    READ(tempUnit,20) (readiBlockSize, readjBlockSize, n_=1, readnBlocks)

    DO n_ = 1, readnBlocks
      READ(tempUnit,30) tRef,dum,dum,dum
      READ(tempUnit,30) ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                        ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                        ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                        ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1)
    END DO

    ! Close file
    CLOSE(tempUnit)

    !     write(*,*), 'Processor ', MyID, ' read temperature file.'
  END SUBROUTINE

  ! Plot3D routine to output grid and temperature files in a
  ! machine readable format.
  SUBROUTINE plotProcs(Procs)
    TYPE (Proc), TARGET :: Procs(:)
    TYPE (BlockType), POINTER :: Blocks(:)
    INTEGER :: globn, n_, p_, i, j
    CHARACTER(2) :: name1, str1, name2, str2

    ! Format statements
    10     FORMAT(I10)
    20     FORMAT(10I10)
    30     FORMAT(10E20.8)

    globn = 1
    DO p_ = 1, mpi_nprocs
      ! Convert integer to string for filename concat.
      WRITE( name1, '(i2)' )  mpi_nprocs
      READ( name1, * ) str1
      WRITE( name2, '(i2)' )  Procs(p_)%procID
      READ( name2, * ) str2

      ! Open files
      OPEN(UNIT=gridUnit,FILE=TRIM(TRIM(str1)//'/grid.dat.p'//str2),FORM='formatted')
      OPEN(UNIT=tempUnit,FILE=TRIM(TRIM(str1)//'/temp.dat.p'//str2),FORM='formatted')

      ! Write to grid file
      WRITE(gridUnit,10) Procs(p_)%nBlocks
      WRITE(gridUnit,20) (iBlockSize, jBlockSize, n_=1, Procs(p_)%nBlocks)

      Blocks => Procs(p_)%Blocks
      DO n_ = globn, globn + Procs(p_)%nBlocks - 1
        WRITE(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                           ((Blocks(n_)%Points(i,j)%y,i=0,iBlockSize+1),j=0,jBlockSize+1)
      END DO

      ! Write to temperature file
      WRITE(tempUnit,10) Procs(p_)%nBlocks
      WRITE(tempUnit,20) (iBlockSize, jBlockSize, n_=1, Procs(p_)%nBlocks)

      DO n_ = globn, globn + Procs(p_)%nBlocks - 1
        WRITE(tempUnit,30) tRef,dum,dum,dum
        WRITE(tempUnit,30) ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                           ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                           ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                           ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1)
      END DO

      ! Close files
      CLOSE(gridUnit)
      CLOSE(tempUnit)
    END DO
  END SUBROUTINE plotProcs

  ! Plot3D routine to output grid and temperature files.
  SUBROUTINE plot3D(Blocks)
    TYPE (BlockType) :: Blocks(:)
    INTEGER :: n_, i, j
    CHARACTER(2) :: name1, str1, name2, str2

    ! Format statements
    10     FORMAT(I10)
    20     FORMAT(10I10)
    30     FORMAT(10E20.8)

    WRITE( name1, '(i2)' )  mpi_nprocs
    READ( name1, * ) str1
    WRITE( name2, '(i2)' )  MyID
    READ( name2, * ) str2

    ! Open files
    OPEN(UNIT=gridUnit,FILE=TRIM(TRIM(str1)//'/final_grid.dat.p'//str2), FORM='formatted')
    OPEN(UNIT=tempUnit,FILE=TRIM(TRIM(str1)//'/final_temp.dat.p'//str2), FORM='formatted')

    ! Write to grid file
    WRITE(gridUnit,10) MyNBlocks
    WRITE(gridUnit,20) (iBlockSize, jBlockSize, n_=1, MyNBlocks)

    DO n_ = 1, MyNBlocks
      WRITE(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%y,i=1,iBlockSize),j=1,jBlockSize)
    END DO

    ! Write to temperature file
    WRITE(tempUnit,10) MyNBlocks
    WRITE(tempUnit,20) (iBlockSize, jBlockSize, n_=1, MyNBlocks)

    DO n_ = 1, MyNBlocks
      WRITE(tempUnit,30) tRef,dum,dum,dum
      WRITE(tempUnit,30) ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize)
    END DO

    ! Close files
    CLOSE(gridUnit)
    CLOSE(tempUnit)
  END SUBROUTINE plot3D

  ! Some output so we know something happened.
  SUBROUTINE output(Blocks, residuals)
    TYPE (BlockType) :: Blocks(:)
    REAL(KIND=8) :: residuals(:)
    INTEGER :: n_, max_n = 1
    INTEGER :: i, j, max_i, max_j
    REAL(KIND=8) :: temp_residual = 0.d0, residual = 0.d0
    CHARACTER(2) :: name1, str1, name2, str2

    IF (MyID == 0) THEN
      ! Write the residual information to output file.
      WRITE( name1, '(i2)' )  mpi_nprocs
      READ( name1, * ) str1
      WRITE( name2, '(i2)' )  MyID
      READ( name2, * ) str2

      OPEN(UNIT = 666, FILE = TRIM(TRIM(str1)//'/convergence.dat.p'//str2), FORM='formatted')
      WRITE(666,*), residuals
      CLOSE(666)

      ! Check for convergence.
      IF (step < max_steps) THEN
        WRITE(*,*) "Converged."
      ELSE
        WRITE(*,*) "Failed to converge."
      END IF

      WRITE(*,*)
    END IF

    IF ( step > 0) THEN
      ! Loop over blocks, find largest residual.
      DO n_ = 1, MyNBlocks
        temp_residual = MAXVAL(ABS(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))

        IF (temp_residual > residual) THEN
          max_n = n_
          residual = temp_residual
        END IF
      END DO

      CALL MPI_Barrier(mpi_comm_world, ierror)
      CALL MPI_Bcast(wall_time, 1, MPI_REAL8, 0, mpi_comm_world, ierror)
      IF (residual == MINVAL(residuals, DIM=1, mask=(residuals>0))) THEN
        WRITE (*,*), "steps ", step
        WRITE (*,*), "residual ", residual
        residual = 0.d0

        ! Loop over points in block, find largest residual.
        DO j = Blocks(max_n)%localJMIN, Blocks(max_n)%localJMAX
          DO i =  Blocks(max_n)%localIMIN, Blocks(max_n)%localIMAX
            temp_residual = ABS(Blocks(max_n)%Points(i, j)%tempT)

            IF (temp_residual > residual) THEN
              residual = temp_residual
              max_i = i
              max_j = j
            END IF
          END DO
        END DO

        ! Readjust to global i, j
        max_i = Blocks(max_n)%lowI + max_i - 2
        max_j = Blocks(max_n)%lowJ + max_j - 2

        WRITE (*,*), "proc ", MyID, " has ij ", max_i, max_j
      END IF
    END IF
  END SUBROUTINE
END MODULE plot3D_module
