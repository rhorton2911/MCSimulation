        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb 14 12:23:02 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TESTTICS__genmod
          INTERFACE 
            SUBROUTINE TESTTICS(SDCSOB,STATEBASIS)
              USE STATE_CLASS
              USE SDCS_MODULE
              TYPE (BASIS_SDCS) :: SDCSOB
              TYPE (BASIS_STATE) :: STATEBASIS
            END SUBROUTINE TESTTICS
          END INTERFACE 
        END MODULE TESTTICS__genmod
