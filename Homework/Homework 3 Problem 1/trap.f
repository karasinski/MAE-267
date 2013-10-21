      SUBROUTINE TRAP(A,B,N,INTEGRAL)
      !CALCULATE INTEGRAL OF FUNCTION WITH EVEN INTERVAL
      IMPLICIT NONE
      REAL*8  :: A,B,F,H,INTEGRAL,X
      INTEGER :: N,I
      F(X) = 4.0/(1.+X**2)        !FUNCTION DEFINITION
      H(A,B,N) = (B-A)/REAL(N)    !INTERVAL DEFINITION
      INTEGRAL = (F(A)+F(B))/2.   !INITIALIZE INTEGRAL
      DO I = 1,N-1
         X = A + REAL(I)*H(A,B,N)
         INTEGRAL = INTEGRAL + F(X)
      END DO
      INTEGRAL = H(A,B,N)*INTEGRAL
      RETURN
      END

      subroutine simpson(a,b,n,integral)
      implicit none
      real*8 f, a, b, integral, s
      real*8 h, x
      integer nint
      integer n, i
      f(x) = 4.0/(1.+x**2)
      ! if n is odd we add +1 to make it even
      ! PURPOSEFUL mixed mode arithmetic
      if((n/2)*2 /= n) n=n+1
      ! loop over n (number of intervals)
      s = 0.0
      h = (b-a)/real(n)
      do i=2, n-2, 2
         x = a+real(i)*h
         s = s + 2.0*f(x) + 4.0*f(x+h)
      end do
      integral = (s + f(a) + f(b) + 4.0*f(a+h))*h/3.0
      return
      end subroutine simpson