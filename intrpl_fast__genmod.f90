        !COMPILER-GENERATED INTERFACE MODULE: Fri Dec  4 10:33:41 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTRPL_FAST__genmod
          INTERFACE 
            SUBROUTINE INTRPL_FAST(L,X,Y,N,I1,I2,U,V)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: L
              REAL(KIND=8), INTENT(IN) :: X(L)
              REAL(KIND=8), INTENT(IN) :: Y(L)
              INTEGER(KIND=4) :: I1
              INTEGER(KIND=4) :: I2
              REAL(KIND=8), INTENT(IN) :: U(N)
              REAL(KIND=8), INTENT(OUT) :: V(N)
            END SUBROUTINE INTRPL_FAST
          END INTERFACE 
        END MODULE INTRPL_FAST__genmod
