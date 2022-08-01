        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar  3 16:31:42 2022
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
