      PROGRAM CALCPIP
      IMPLICIT NONE

      include "mpif.h"
      REAL*8  :: A,AK,B,BK,H,PI,SUBPI
      INTEGER :: K,MYID,N,NK,NPROCS
      INTEGER :: IERROR,TAG,STATUS
      REAL*8  :: time

      ! INITIALIZE MPI
      CALL MPI_Init(IERROR)

      ! DETERMINE MY PROCESSOR ID
      CALL MPI_Comm_rank(MPI_COMM_WORLD,MYID,IERROR)

      ! FIND OUT HOW MANY PROCESSORS ARE USED
      CALL MPI_Comm_size(MPI_COMM_WORLD,NPROCS,IERROR)

      IF(MYID == 0) THEN
        !Start the clock.
        time = MPI_Wtime()

        N = 60000
      END IF

      ! BROADCAST THE NUMBER OF SUB-INTERVALS
      CALL MPI_Bcast(N,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)

      A = 0.0   !DEFINE INTERVAL START
      B = 1.0   !DEFINE INTERVAL STOP
      H  = (B-A)/REAL(N)

      ! N INTERVALS MUST BE EVENLY DIVISIBLE BY NPROCS
      NK = N/NPROCS
      AK = A + REAL(MYID)*REAL(NK)*H
      BK = AK + REAL(NK)*H

      ! COMPUTE LOCAL INTEGRAL
      CALL simpson(AK,BK,NK,SUBPI)
 
      ! SET UP A MASTER-SLAVE RELATIONSHIP WHERE THE MASTER
      ! IS RESPONSIBLE FOR ACCUMULATING THE SUB-INTEGRALS
      ! AND WRITING OUT THE ANSWER
      IF(MYID == 0) THEN
        ! SUM UP THE INTEGRALS FROM THE OTHER PROCESSORS
        PI = SUBPI
        ! ADD THE SUBPI'S FROM THE OTHER PROCESSORS
        DO K = 1,NPROCS-1
           CALL MPI_Recv(SUBPI,1,MPI_DOUBLE_PRECISION,K,K, 
     &                   MPI_COMM_WORLD,STATUS,IERROR)
           PI = PI + SUBPI
        END DO
        PRINT *,'PI = ',PI
      ELSE
        ! SEND THE INTEGRAL TO THE MASTER
        CALL MPI_Send(SUBPI,1,MPI_DOUBLE_PRECISION,0,MYID,
     &                MPI_COMM_WORLD,IERROR)
      END IF
      ! Syncronize processors
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

      IF(MYID == 0) THEN
        time = MPI_Wtime() - time
        write(*,*), "Elapsed time is ", time
        write(*,*), NPROCS, " processors"
        write(*,*), "Or ", (time/REAL(NPROCS)), " seconds per processor"
      END IF
      ! TERMINATE MPI
 1000 CALL MPI_Finalize(IERROR)

      STOP
      END
