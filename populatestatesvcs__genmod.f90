        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar  3 16:31:42 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POPULATESTATESVCS__genmod
          INTERFACE 
            SUBROUTINE POPULATESTATESVCS(STATEBASIS,TCSBASIS)
              USE TOTALCS_MODULE
              USE STATE_CLASS
              TYPE (BASIS_STATE), INTENT(INOUT) :: STATEBASIS
              TYPE (BASIS_TOTALCS), INTENT(IN) :: TCSBASIS
            END SUBROUTINE POPULATESTATESVCS
          END INTERFACE 
        END MODULE POPULATESTATESVCS__genmod
