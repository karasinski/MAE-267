module plot3D_module
  use clock
  use BlockModule
  implicit none

  integer :: gridUnit  = 30   ! Unit for grid file
  integer :: tempUnit = 21    ! Unit for temp file
  integer :: configUnit = 99    ! Unit for temp file
  real(kind=8) :: tRef = 1.d0 ! tRef number
  real(kind=8) :: dum = 0.d0  ! dummy values

contains
  subroutine write_configuration_file(Procs)
    type (Proc), target :: Procs(:)
    type (BlockType), pointer :: BlocksCollection(:)
    type (BlockType), pointer :: b
    integer :: n_, p_
    character(2) :: name, str

    10 format(3I5)
    20 format(33I5)

    do p_ = 1, nProcs
      BlocksCollection => Procs(p_)%Blocks

      write( name, '(i2)' )  p_
      read( name, * ) str

      open(unit = configUnit, file = 'configuration_file.dat.p'//str, form='formatted')
      write(configUnit, 10) M*N, iBlockSize, jBlockSize
      do n_ = 1, Procs(p_)%nBlocks
        b => BlocksCollection(n_)
        write(configUnit, 20) n_, b%proc, b%size, &
                     b%lowI, b%lowJ, &
                     b%localIMIN, b%localIMAX, b%localJMIN, b%localJMAX, &
                     b%northFace%BC, b%northFace%neighborBlock, b%northFace%neighborProc, &
                     b%eastFace%BC, b%eastFace%neighborBlock, b%eastFace%neighborProc, &
                     b%southFace%BC, b%southFace%neighborBlock, b%southFace%neighborProc, &
                     b%westFace%BC, b%westFace%neighborBlock, b%westFace%neighborProc, &
                     b%NECorner%BC, b%NECorner%neighborBlock, b%NECorner%neighborProc, &
                     b%NWCorner%BC, b%NWCorner%neighborBlock, b%NWCorner%neighborProc, &
                     b%SWCorner%BC, b%SWCorner%neighborBlock, b%SWCorner%neighborProc, &
                     b%SECorner%BC, b%SECorner%neighborBlock, b%SECorner%neighborProc
      end do
      close(configUnit)
    end do

  end subroutine

  subroutine read_configuration_file(Blocks)
    type (BlockType), target :: Blocks(:)
    type (BlockType), pointer :: b
    integer :: n_, nFile
    integer :: readnBlocks, readiBlockSize, readjBlockSize

    open(unit = 1, file = 'configuration_file.dat', status='old')

    10 format(3I5)
    20 format(33I5)

    read(1, 10) readnBlocks, readiBlockSize, readjBlockSize
    do n_ = 1, nBlocks
      b => Blocks(n_)
      read(1, 20) nFile, b%proc, b%size, &
                  b%lowI, b%lowJ, &
                  b%localIMIN, b%localIMAX, b%localJMIN, b%localJMAX, &
                  b%northFace%BC, b%northFace%neighborBlock, b%northFace%neighborProc, &
                  b%eastFace%BC, b%eastFace%neighborBlock, b%eastFace%neighborProc, &
                  b%southFace%BC, b%southFace%neighborBlock, b%southFace%neighborProc, &
                  b%westFace%BC, b%westFace%neighborBlock, b%westFace%neighborProc, &
                  b%NECorner%BC, b%NECorner%neighborBlock, b%NECorner%neighborProc, &
                  b%NWCorner%BC, b%NWCorner%neighborBlock, b%NWCorner%neighborProc, &
                  b%SWCorner%BC, b%SWCorner%neighborBlock, b%SWCorner%neighborProc, &
                  b%SECorner%BC, b%SECorner%neighborBlock, b%SECorner%neighborProc
    end do

    close(1)
    write(*,*), 'Read configuration file'
  end subroutine

  subroutine read_grid_file(Blocks)
    type (BlockType) :: Blocks(:)
    integer :: n_, i, j
    integer :: readnBlocks, readiBlockSize, readjBlockSize

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    ! Open files
    open(unit=gridUnit,file='i_grid.dat', status='old')

    ! Read grid file
    read(gridUnit,10) readnBlocks
    read(gridUnit,20) (readiBlockSize, readjBlockSize, i=1, nBlocks)
    do n_ = 1, nBlocks
      read(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=1,iBlockSize),j=1,jBlockSize), &
                        ((Blocks(n_)%Points(i,j)%y,i=1,iBlockSize),j=1,jBlockSize)
    end do

    ! Close file
    close(gridUnit)
    write(*,*), 'Read grid file'
  end subroutine

  subroutine read_temp_file(Blocks)
    type (BlockType) :: Blocks(:)
    integer :: n_, i, j
    integer :: readnBlocks, readiBlockSize, readjBlockSize

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    ! Open files
    open(unit=tempUnit,file='i_temp.dat', status='old')

    ! Read temperature file
    read(tempUnit,10) readnBlocks
    read(tempUnit,20) (readiBlockSize, readjBlockSize, n_=1, nBlocks)

    do n_ = 1, nBlocks
      read(tempUnit,30) tRef,dum,dum,dum
      read(tempUnit,30) ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                        ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                        ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize), &
                        ((Blocks(n_)%Points(i,j)%T,i=1,iBlockSize),j=1,jBlockSize)
    end do

    ! Close file
    close(tempUnit)

    write(*,*), 'Read initial temperature file'
  end subroutine

  ! Plot3D routine to output grid and temperature files in a
  ! machine readable format.
  subroutine plotProcs(Procs)
    type (Proc), target :: Procs(:)
    type (BlockType), pointer :: Blocks(:)
    character(2) :: name, str
    integer :: globn, n_, p_, i, j

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    globn = 1
    do p_ = 1, nProcs
      ! Convert integer to string for filename concat.
      write(name,'(i2)')  p_
      read (name,*) str

      ! Open files
      open(unit=gridUnit,file='grid.dat.p'//str,form='formatted')
      open(unit=tempUnit,file='temp.dat.p'//str,form='formatted')

      ! Write to grid file
      write(gridUnit,10) Procs(p_)%nBlocks
      write(gridUnit,20) (iBlockSize, jBlockSize, n_=1, Procs(p_)%nBlocks)

      Blocks => Procs(p_)%Blocks
      do n_ = globn, globn + Procs(p_)%nBlocks - 1
!        write(*,*), Blocks(n_)%id
        write(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=1,iBlockSize),j=1,jBlockSize), &
                           ((Blocks(n_)%Points(i,j)%y,i=1,iBlockSize),j=1,jBlockSize)
      end do

      ! Write to temperature file
      write(tempUnit,10) Procs(p_)%nBlocks
      write(tempUnit,20) (iBlockSize, jBlockSize, n_=1, Procs(p_)%nBlocks)

      do n_ = globn, globn + Procs(p_)%nBlocks - 1
        write(tempUnit,30) tRef,dum,dum,dum
        write(tempUnit,30) ((real(p_),i=1,iBlockSize),j=1,jBlockSize), &
                           ((real(p_),i=1,iBlockSize),j=1,jBlockSize), &
                           ((real(p_),i=1,iBlockSize),j=1,jBlockSize), &
                           ((real(p_),i=1,iBlockSize),j=1,jBlockSize)

                           ! Should be Blocks(n_)%Points(i,j)%T rather than real(p_), but this
                           ! makes it easy to see which processors have which blocks.
      end do

      ! Close files
      close(gridUnit)
      close(tempUnit)
    end do

  end subroutine plotProcs


  ! Plot3D routine to output grid and temperature files in a
  ! machine readable format.
  subroutine plot3D(Blocks, name)
    type (BlockType) :: Blocks(:)
    character :: name
    integer :: n_, i, j

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    ! Open files
    open(unit=gridUnit,file=name//'_grid.dat',form='formatted')
    open(unit=tempUnit,file=name//'_temp.dat',form='formatted')

    ! Write to grid file
    write(gridUnit,10) nBlocks
    write(gridUnit,20) (iBlockSize, jBlockSize, n_=1, nBlocks)

    do n_ = 1, nBlocks
      write(gridUnit,30) ((Blocks(n_)%Points(i,j)%x,i=1,iBlockSize),j=1,jBlockSize), &
                         ((Blocks(n_)%Points(i,j)%y,i=1,iBlockSize),j=1,jBlockSize)
    end do

    ! Write to temperature file
    write(tempUnit,10) nBlocks
    write(tempUnit,20) (iBlockSize, jBlockSize, n_=1, nBlocks)

    do n_ = 1, nBlocks
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
  subroutine output(Blocks)
    type (BlockType) :: Blocks(:)
    integer :: n_, max_n = 1
    integer :: i, j, max_i, max_j
    real(kind=8) :: temp_residual, residual = 0.d0

    ! Write down misc. info asked for by Prof.
    if ( step > 0) then
      open (unit = 2, file = "info.dat")
      write (2,*), wall_time, " seconds"
      write (*,*), "steps ", step
      write (2,*), "steps ", step

      ! Loop over blocks, find largest residual.
      do n_ = 1, nBlocks
        temp_residual = maxval(abs(Blocks(n_)%Points(2:iBlockSize-1, 2:jBlockSize-1)%tempT))

        if (temp_residual > residual) then
          max_n = n_
          residual = temp_residual
        end if
      end do

      write (*,*), "residual ", residual
      write (2,*), "residual ", residual
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

      write (*,*), "ij ", max_i, max_j
      write (2,*), "ij ", max_i, max_j

      close(2)
    end if

  end subroutine
end module plot3D_module
