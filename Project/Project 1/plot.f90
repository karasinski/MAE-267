!-------------------
!
! MODULE:		plot3D_module
! AUTHOR:		Todd Beeby
! CLASS: 		MAE 267 Fall 2013
! ASSIGNMENT:	Project 1
! DESCRIPTION:	This module creates a grid and temperature file in
!               the plot3D format for outr steady state solution
!
!-------------------

module plot3D_module
    use constants
    use GridPointModule
    implicit none

    ! Variables
    integer :: gridUnit  = 30   ! Unit for grid file
    integer :: tempUnit = 21    ! Unit for temp file
    real*8 :: tRef = 1          ! tRef number
    real*8 :: dum = 0.          ! dummy values
    integer :: nBlocks = 1      ! number of blocks

    contains
    subroutine plot3D(Points)
        implicit none

        type (GridPoint) :: Points(1:IMAX, 1:JMAX)
        integer :: i, j

        ! Format statements
        10     format(I10)
        20     format(10I10)
        30     format(10E20.8)

        ! Open files
        open(unit=gridUnit,file='grid.dat',form='formatted')
        open(unit=tempUnit,file='temp.dat',form='formatted')

        ! Write to grid file
        write(gridUnit,10) nBlocks
        write(gridUnit,20) IMAX,JMAX
        write(gridUnit,30) ((Points(i,j)%x,i=1,IMAX),j=1,JMAX), ((Points(i,j)%y,i=1,IMAX),j=1,JMAX)

        ! Write to temp file
        write(tempUnit,10) nBlocks
        write(tempUnit,20) IMAX,JMAX
        write(tempUnit,30) tRef,dum,dum,dum
        write(tempUnit,30) ((Points(i,j)%T,i=1,IMAX),j=1,JMAX), ((Points(i,j)%T,i=1,IMAX),j=1,JMAX), &
                           ((Points(i,j)%T,i=1,IMAX),j=1,JMAX), ((Points(i,j)%T,i=1,IMAX),j=1,JMAX)

        ! Close files
        close(gridUnit)
        close(tempUnit)
    end subroutine plot3D
end module plot3D_module
