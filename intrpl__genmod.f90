        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar  3 16:31:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTRPL__genmod
          INTERFACE 
            SUBROUTINE INTRPL(L,X,Y,N,U,V)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: L
              REAL(KIND=8) :: X(L)
              REAL(KIND=8) :: Y(L)
              REAL(KIND=8) :: U(N)
              REAL(KIND=8) :: V(N)
            END SUBROUTINE INTRPL
          END INTERFACE 
        END MODULE INTRPL__genmod
