module plot3D_module
  use constants
  use clock
  use GridPointModule
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

    type (GridPoint), allocatable :: Blocks(:,:,:,:)
    integer :: m_, n_, M, N
    integer :: i, j, iBound, jBound

    ! Read M and N from Blocks.
    M = size(Blocks, 1)
    N = size(Blocks, 2)

    ! Set number of blocks.
    nBlocks = N * M

    ! Size of each block.
    iBound = size(Blocks, 3)
    jBound = size(Blocks, 4)

    ! Format statements
    10     format(I10)
    20     format(10I10)
    30     format(10E20.8)

    ! Open files
    open(unit=gridUnit,file='grid.dat',form='formatted')
    open(unit=tempUnit,file='temperature.dat',form='formatted')

    ! Clever scheme to find the true number of nodes in each block.
    ! May want this later, though it has no use now.
    !
    !    do m_ = 1, M
    !      do n_ = 1, N
    !        write(*, *), (maxval(Blocks(m_,n_,:,:)%i)-minval(Blocks(m_,n_,:,:)%i, MASK = Blocks(m_,n_,:,:)%i>0)),  &
    !                     (maxval(Blocks(m_,n_,:,:)%j)-minval(Blocks(m_,n_,:,:)%j, MASK = Blocks(m_,n_,:,:)%j>0))
    !      end do
    !    end do

    ! Write to grid file
    write(gridUnit,10) nBlocks
    m_ = 1
    write(gridUnit,20) ((iBound,jBound, m_=1, M), n_=1, N)
    do m_ = 1, M
      do n_ = 1, N
        write(gridUnit,30) ((Blocks(m_,n_,i,j)%xp,i=1,iBound),j=1,jBound), &
                           ((Blocks(m_,n_,i,j)%yp,i=1,iBound),j=1,jBound)
      end do
    end do

    ! Write to temperature file
    write(tempUnit,10) nBlocks
    m_ = 1
    write(tempUnit,20) ((iBound,jBound, m_=1, M), n_=1, N)
    do m_ = 1, M
      do n_ = 1, N
        write(tempUnit,30) tRef,dum,dum,dum
        write(tempUnit,30) ((Blocks(m_,n_,i,j)%T,i=1,iBound),j=1,jBound), &
                           ((Blocks(m_,n_,i,j)%T,i=1,iBound),j=1,jBound), &
                           ((Blocks(m_,n_,i,j)%T,i=1,iBound),j=1,jBound), &
                           ((Blocks(m_,n_,i,j)%T,i=1,iBound),j=1,jBound)
      end do
    end do

    ! Close files
    close(gridUnit)
    close(tempUnit)
  end subroutine plot3D

  ! Handcrafted output routine to give us some info.
  subroutine output(Points, step)
    type (GridPoint), target :: Points(1:IMAX, 1:JMAX)
    real(kind=8), pointer :: Temperature(:,:), tempTemperature(:,:)
    integer :: step
    !    integer :: i, j

    Temperature => Points(2:IMAX-1, 2:JMAX-1)%T
    tempTemperature => Points(2:IMAX-1, 2:JMAX-1)%tempT
    ! Let's find the last cell to change temperature and write some output.
    ! Write down the 'steady state' configuration.
    !    open (unit = 1, file = "steady_state.dat")
    !    do i = 1, IMAX
    !      do j = 1, JMAX
    !        write (1,'(F10.7, 5X, F10.7, 5X, F10.7, I5, F10.7)'), Points(i,j)%x, Points(i,j)%y, Points(i,j)%T
    !      end do
    !    end do
    !    close (1)

    ! Some output to the screen so we know something happened.
    write (*,*), "IMAX/JMAX", IMAX, JMAX
    write (*,*), "steps", step
    write (*,*), "residual", maxval(tempTemperature)
    write (*,*), "ij", maxloc(tempTemperature)

    ! Write down misc. info asked for by Prof.
    open (unit = 2, file = "info.dat")
    write (2,*), "For a ", IMAX, " by ", JMAX, "size grid, we ran for: "
    write (2,*), step, "steps"
    write (2,*), wall_time, "seconds"
    write (2,*)
    write (2,*), "Found max residual of ", maxval(tempTemperature)
    write (2,*), "At ij of ", maxloc(tempTemperature)
    close (2)
    ! End output.
  end subroutine
end module plot3D_module
