        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb 14 11:39:58 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POPULATESTATES__genmod
          INTERFACE 
            SUBROUTINE POPULATESTATES(STATEBASIS,TCSBASIS)
              USE TOTALCS_MODULE
              USE STATE_CLASS
              TYPE (BASIS_TOTALCS), INTENT(IN) :: TCSBASIS
              TYPE (BASIS_STATE), INTENT(INOUT) :: STATEBASIS
            END SUBROUTINE POPULATESTATES
          END INTERFACE 
        END MODULE POPULATESTATES__genmod
