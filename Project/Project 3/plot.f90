module plot3D_module
  use constants
  use clock
  use GridPointModule
  use BlockModule
  implicit none

  integer :: gridUnit  = 30   ! Unit for grid file
  integer :: tempUnit = 21    ! Unit for temp file
  real(kind=8) :: tRef = 1    ! tRef number
  real(kind=8) :: dum = 0.d0  ! dummy values
  integer :: nBlocks = 1      ! default number of blocks

contains
  ! Plot3D routine to output grid and temperature files in a
  ! machine readable format.
  subroutine plot3D(Blocks)
    implicit none

    type (BlockType) :: Blocks(:,:)
    integer :: m_, n_, M, N
    integer :: i, j
    integer :: iBound, jBound

    ! Read M and N from Blocks.
    M = size(Blocks, 1)
    N = size(Blocks, 2)

    ! Set number of blocks.
    nBlocks = N * M

    ! Size of each block.
    iBound = 1 + (IMAX - 1) / N
    jBound = 1 + (JMAX - 1) / M

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    ! Open files
    open(unit=gridUnit,file='grid.dat',form='formatted')
    open(unit=tempUnit,file='temperature.dat',form='formatted')

    ! Write to grid file
    write(gridUnit,10) nBlocks
    m_ = 1
    write(gridUnit,20) ((iBound,jBound, m_=1, M), n_=1, N)
    do m_ = 1, M
      do n_ = 1, N
        write(gridUnit,30) ((Blocks(m_,n_)%Points(i,j)%x,i=1,iBound),j=1,jBound), &
                           ((Blocks(m_,n_)%Points(i,j)%y,i=1,iBound),j=1,jBound)
      end do
    end do

    ! Write to temperature file
    write(tempUnit,10) nBlocks
    m_ = 1
    write(tempUnit,20) ((iBound,jBound, m_=1, M), n_=1, N)
    do m_ = 1, M
      do n_ = 1, N
        write(tempUnit,30) tRef,dum,dum,dum
        write(tempUnit,30) ((Blocks(m_,n_)%Points(i,j)%T,i=1,iBound),j=1,jBound), &
                           ((Blocks(m_,n_)%Points(i,j)%T,i=1,iBound),j=1,jBound), &
                           ((Blocks(m_,n_)%Points(i,j)%T,i=1,iBound),j=1,jBound), &
                           ((Blocks(m_,n_)%Points(i,j)%T,i=1,iBound),j=1,jBound)
      end do
    end do


    ! Write to grid file
!    write(gridUnit,10) nBlocks
!    m_ = 1
!    write(gridUnit,20) ((Blocks(m_,n_)%iBound,Blocks(m_,n_)%jBound, m_=1, M), n_=1, N)
!    do m_ = 1, M
!      do n_ = 1, N
!        write(gridUnit,30) ((Blocks(m_,n_)%Points(i,j)%x,i=1,Blocks(m_,n_)%iBound),j=1,Blocks(m_,n_)%jBound), &
!                           ((Blocks(m_,n_)%Points(i,j)%y,i=1,Blocks(m_,n_)%iBound),j=1,Blocks(m_,n_)%jBound)
!      end do
!    end do
!
!    ! Write to temperature file
!    write(tempUnit,10) nBlocks
!    m_ = 1
!    write(tempUnit,20) ((Blocks(m_,n_)%iBound,Blocks(m_,n_)%jBound, m_=1, M), n_=1, N)
!    do m_ = 1, M
!      do n_ = 1, N
!        write(tempUnit,30) tRef,dum,dum,dum
!        write(tempUnit,30) ((Blocks(m_,n_)%Points(i,j)%T,i=1,Blocks(m_,n_)%iBound),j=1,Blocks(m_,n_)%jBound), &
!                           ((Blocks(m_,n_)%Points(i,j)%T,i=1,Blocks(m_,n_)%iBound),j=1,Blocks(m_,n_)%jBound), &
!                           ((Blocks(m_,n_)%Points(i,j)%T,i=1,Blocks(m_,n_)%iBound),j=1,Blocks(m_,n_)%jBound), &
!                           ((Blocks(m_,n_)%Points(i,j)%T,i=1,Blocks(m_,n_)%iBound),j=1,Blocks(m_,n_)%jBound)
!      end do
!    end do

    ! Close files
    close(gridUnit)
    close(tempUnit)
  end subroutine plot3D
! Basic info output.
subroutine output(Blocks, step)
  type (BlockType) :: Blocks(:,:)
  integer :: step
  integer :: m_,n_,max_m,max_n
  real(kind=8) :: temp_residual, residual

  ! Some output to the screen so we know something happened.
  do m_ = 1, size(Blocks, 1)
    do n_ = 1, size(Blocks, 2)
      temp_residual = maxval(abs(Blocks(m_,n_)%Points(2:Blocks(m_,n_)%iBound - 1, 2:Blocks(m_,n_)%iBound - 1)%tempT))

      if (temp_residual > residual) then
        max_m = m_
        max_n = n_
        residual = temp_residual
      end if

    end do
  end do

  if ( step > 0) then
    write (*,*), "steps", step
    write (*,*), "residual", residual
    write (*,*), "mn", max_m, max_n, "ij", maxloc(abs(Blocks(max_m,max_n)%Points(2:IMAX - 1, 2:JMAX - 1)%tempT))

    ! Write down misc. info asked for by Prof.
    open (unit = 2, file = "info.dat")
    write (2,*), "For a ", IMAX, " by ", JMAX, "size grid, we ran for: "
    write (2,*), step, "steps"
    write (2,*), wall_time, "seconds"
    write (2,*)
    write (2,*), "Found max residual of ", residual
    write (2,*), "mn",max_m, max_n,"ij", maxloc(abs(Blocks(max_m,max_n)%Points(2:IMAX - 1, 2:JMAX - 1)%tempT))
    close (2)
  end if
end subroutine
end module plot3D_module
