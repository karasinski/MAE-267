module plot3D_module
  use clock
  use BlockModule
  implicit none

  integer :: gridUnit  = 30   ! Unit for grid file
  integer :: tempUnit = 21    ! Unit for temp file
  integer :: configUnit = 99  ! Unit for config file
  real(kind=8) :: tRef = 1.d0 ! tRef number
  real(kind=8) :: dum = 0.d0  ! dummy values

contains
  subroutine write_configuration_file(Procs)
    type (Proc), target :: Procs(:)
    type (BlockType), pointer :: BlocksCollection(:)
    type (BlockType), pointer :: b
    integer :: n_, p_
    character(2) :: name1, str1, name2, str2

    write( name1, '(i2)' )  mpi_nprocs
    read( name1, * ) str1

    ! Make a directory for this run
    call execute_command_line ('mkdir -p ' // str1 )

    10 format(3I5)
    20 format(33I5)

    do p_ = 1, mpi_nprocs
      BlocksCollection => Procs(p_)%Blocks

      write( name2, '(i2)' )  Procs(p_)%procID
      read( name2, * ) str2

      open(unit = configUnit, file = trim(trim(str1)//'/configuration_file.dat.p'//str2), form='formatted')
      write(configUnit, 10) Procs(p_)%nBlocks, iBlockSize, jBlockSize
      do n_ = 1, Procs(p_)%nBlocks
        b => BlocksCollection(n_)
        write(configUnit, 20) b%id, b%proc, b%size, &
          b%lowI, b%lowJ, &
          b%localIMIN, b%localIMAX, b%localJMIN, b%localJMAX, &
          b%northFace%BC, b%northFace%neighborBlock, b%northFace%neighborLocalBlock, b%northFace%neighborProc, &
          b%southFace%BC, b%southFace%neighborBlock, b%southFace%neighborLocalBlock, b%southFace%neighborProc, &
          b%eastFace%BC,  b%eastFace%neighborBlock,  b%eastFace%neighborLocalBlock,  b%eastFace%neighborProc,  &
          b%westFace%BC,  b%westFace%neighborBlock,  b%westFace%neighborLocalBlock,  b%westFace%neighborProc, &
          b%NECorner%BC, b%NECorner%neighborBlock, b%NECorner%neighborLocalBlock,  b%NECorner%neighborProc, &
          b%NWCorner%BC, b%NWCorner%neighborBlock, b%NWCorner%neighborLocalBlock,  b%NWCorner%neighborProc, &
          b%SWCorner%BC, b%SWCorner%neighborBlock, b%SWCorner%neighborLocalBlock,  b%SWCorner%neighborProc, &
          b%SECorner%BC, b%SECorner%neighborBlock, b%SECorner%neighborLocalBlock,  b%SECorner%neighborProc
      end do
      close(configUnit)
    end do

  end subroutine

  subroutine read_configuration_file(Blocks)
    type (BlockType), pointer :: Blocks(:)
    type (BlockType), pointer :: b
    integer :: n_, readnBlocks, readiBlockSize, readjBlockSize
    character(2) :: name1, str1, name2, str2

    write( name1, '(i2)' )  mpi_nprocs
    read( name1, * ) str1
    write( name2, '(i2)' )  MyID
    read( name2, * ) str2
    open(unit = 1, file = trim(trim(str1)//'/configuration_file.dat.p'//str2), status='old')

    10 format(3I5)
    20 format(33I5)

    read(1, 10) readnBlocks, readiBlockSize, readjBlockSize

    ! Set my nblocks.
    MyNBlocks = readnBlocks

    ! Allocate our blocks array.
    allocate(Blocks(1:readnBlocks))

    do n_ = 1, readnBlocks
      b => Blocks(n_)
      read(1, 20) b%id, b%proc, b%size, &
        b%lowI, b%lowJ, &
        b%localIMIN, b%localIMAX, b%localJMIN, b%localJMAX, &
        b%northFace%BC, b%northFace%neighborBlock, b%northFace%neighborLocalBlock, b%northFace%neighborProc, &
        b%southFace%BC, b%southFace%neighborBlock, b%southFace%neighborLocalBlock, b%southFace%neighborProc, &
        b%eastFace%BC,  b%eastFace%neighborBlock,  b%eastFace%neighborLocalBlock,  b%eastFace%neighborProc,  &
        b%westFace%BC,  b%westFace%neighborBlock,  b%westFace%neighborLocalBlock,  b%westFace%neighborProc, &
        b%NECorner%BC, b%NECorner%neighborBlock, b%NECorner%neighborLocalBlock,  b%NECorner%neighborProc, &
        b%NWCorner%BC, b%NWCorner%neighborBlock, b%NWCorner%neighborLocalBlock,  b%NWCorner%neighborProc, &
        b%SWCorner%BC, b%SWCorner%neighborBlock, b%SWCorner%neighborLocalBlock,  b%SWCorner%neighborProc, &
        b%SECorner%BC, b%SECorner%neighborBlock, b%SECorner%neighborLocalBlock,  b%SECorner%neighborProc
    end do

    close(1)
    write(*,*), 'Processor ', MyID, ' read configuration file.'
  end subroutine

  subroutine read_grid_file(Blocks)
    type (BlockType) :: Blocks(:)
    integer :: n_, i, j
    integer :: readnBlocks, readiBlockSize, readjBlockSize
    character(2) :: name1, str1, name2, str2

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    ! Open files
    write( name1, '(i2)' )  mpi_nprocs
    read( name1, * ) str1
    write( name2, '(i2)' )  MyID
    read( name2, * ) str2
    open(unit=gridUnit,file=trim(trim(str1)//'/grid.dat.p'//str2), status='old')

    ! Read grid file
    read(gridUnit,10) readnBlocks
    read(gridUnit,20) (readiBlockSize, readjBlockSize, i=1, readnBlocks)
    do n_ = 1, readnBlocks
      read(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                        ((Blocks(n_)%Points(i,j)%y,i=0,iBlockSize+1),j=0,jBlockSize+1)
    end do

    ! Close file
    close(gridUnit)
    write(*,*), 'Processor ', MyID, ' read grid file.'
  end subroutine

  subroutine read_temp_file(Blocks)
    type (BlockType) :: Blocks(:)
    integer :: n_, i, j
    integer :: readnBlocks, readiBlockSize, readjBlockSize
    character(2) :: name1, str1, name2, str2

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    ! Open files
    write( name1, '(i2)' )  mpi_nprocs
    read( name1, * ) str1
    write( name2, '(i2)' )  MyID
    read( name2, * ) str2
    open(unit=tempUnit,file=trim(trim(str1)//'/temp.dat.p'//str2), status='old')

    ! Read temperature file
    read(tempUnit,10) readnBlocks
    read(tempUnit,20) (readiBlockSize, readjBlockSize, n_=1, readnBlocks)

    do n_ = 1, readnBlocks
      read(tempUnit,30) tRef,dum,dum,dum
      read(tempUnit,30) ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                        ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                        ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                        ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1)
    end do

    ! Close file
    close(tempUnit)

    write(*,*), 'Processor ', MyID, ' read temperature file.'
  end subroutine

  ! Plot3D routine to output grid and temperature files in a
  ! machine readable format.
  subroutine plotProcs(Procs)
    type (Proc), target :: Procs(:)
    type (BlockType), pointer :: Blocks(:)
    integer :: globn, n_, p_, i, j
    character(2) :: name1, str1, name2, str2

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    globn = 1
    do p_ = 1, mpi_nprocs
      ! Convert integer to string for filename concat.
      write( name1, '(i2)' )  mpi_nprocs
      read( name1, * ) str1
      write( name2, '(i2)' )  Procs(p_)%procID
      read( name2, * ) str2

      ! Open files
      open(unit=gridUnit,file=trim(trim(str1)//'/grid.dat.p'//str2),form='formatted')
      open(unit=tempUnit,file=trim(trim(str1)//'/temp.dat.p'//str2),form='formatted')

      ! Write to grid file
      write(gridUnit,10) Procs(p_)%nBlocks
      write(gridUnit,20) (iBlockSize, jBlockSize, n_=1, Procs(p_)%nBlocks)

      Blocks => Procs(p_)%Blocks
      do n_ = globn, globn + Procs(p_)%nBlocks - 1
!        write(*,*), Blocks(n_)%id
        write(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                           ((Blocks(n_)%Points(i,j)%y,i=0,iBlockSize+1),j=0,jBlockSize+1)
      end do

      ! Write to temperature file
      write(tempUnit,10) Procs(p_)%nBlocks
      write(tempUnit,20) (iBlockSize, jBlockSize, n_=1, Procs(p_)%nBlocks)

      do n_ = globn, globn + Procs(p_)%nBlocks - 1
        write(tempUnit,30) tRef,dum,dum,dum
        write(tempUnit,30) ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                           ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                           ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1), &
                           ((Blocks(n_)%Points(i,j)%T,i=0,iBlockSize+1),j=0,jBlockSize+1)
      end do

      ! Close files
      close(gridUnit)
      close(tempUnit)
    end do
  end subroutine plotProcs

  ! Plot3D routine to output grid and temperature files.
  subroutine plot3D(Blocks)
    type (BlockType) :: Blocks(:)
    integer :: n_, i, j
    character(2) :: name1, str1, name2, str2

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    write( name1, '(i2)' )  mpi_nprocs
    read( name1, * ) str1
    write( name2, '(i2)' )  MyID
    read( name2, * ) str2

    ! Open files
    open(unit=gridUnit,file=trim(trim(str1)//'/final_grid.dat.p'//str2), form='formatted')
    open(unit=tempUnit,file=trim(trim(str1)//'/final_temp.dat.p'//str2), form='formatted')

    ! Write to grid file
    write(gridUnit,10) MyNBlocks
    write(gridUnit,20) (iBlockSize, jBlockSize, n_=1, MyNBlocks)

    do n_ = 1, MyNBlocks
      write(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%y,i=1,iBlockSize),j=1,jBlockSize)
    end do

    ! Write to temperature file
    write(tempUnit,10) MyNBlocks
    write(tempUnit,20) (iBlockSize, jBlockSize, n_=1, MyNBlocks)

    do n_ = 1, MyNBlocks
      write(tempUnit,30) tRef,dum,dum,dum
      write(tempUnit,30) ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize)
    end do

    ! Close files
    close(gridUnit)
    close(tempUnit)
  end subroutine plot3D

  ! Some output so we know something happened.
  subroutine output(Blocks, residuals)
    type (BlockType) :: Blocks(:)
    real(kind=8) :: residuals(:)
    integer :: n_, max_n = 1
    integer :: i, j, max_i, max_j
    real(kind=8) :: temp_residual = 0.d0, residual = 0.d0
    character(2) :: name1, str1, name2, str2

    if (MyID == 0) then
      ! Write the residual information to output file.
      write( name1, '(i2)' )  mpi_nprocs
      read( name1, * ) str1
      write( name2, '(i2)' )  MyID
      read( name2, * ) str2

      open(unit = 666, file = trim(trim(str1)//'/convergence.dat.p'//str2), form='formatted')
      write(666,*), residuals
      close(666)

      ! Check for convergence.
      if (step < max_steps) then
        write(*,*) "Converged."
      else
        write(*,*) "Failed to converge."
      end if
    
      write(*,*)
    end if

    if ( step > 0) then
      ! Loop over blocks, find largest residual.
      do n_ = 1, MyNBlocks
        temp_residual = maxval(abs(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))

        if (temp_residual > residual) then
          max_n = n_
          residual = temp_residual
        end if
      end do

      call MPI_Barrier(mpi_comm_world, ierror)
      call MPI_Bcast(wall_time, 1, MPI_DOUBLE_PRECISION, 0, mpi_comm_world, ierror)
      if (residual == minval(residuals, dim=1, mask=(residuals>0))) then
        write (*,*), "steps ", step
        write (*,*), "residual ", residual
        residual = 0.d0

        ! Loop over points in block, find largest residual.
        do j = Blocks(max_n)%localJMIN, Blocks(max_n)%localJMAX
          do i =  Blocks(max_n)%localIMIN, Blocks(max_n)%localIMAX
            temp_residual = abs(Blocks(max_n)%Points(i, j)%tempT)

            if (temp_residual > residual) then
              residual = temp_residual
              max_i = i
              max_j = j
            end if
          end do
        end do

        ! Readjust to global i, j
        max_i = Blocks(max_n)%lowI + max_i - 2
        max_j = Blocks(max_n)%lowJ + max_j - 2

        write (*,*), "proc ", MyID, " has ij ", max_i, max_j
      end if
    end if
  end subroutine
end module plot3D_module
