        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb 14 11:39:58 2022
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
