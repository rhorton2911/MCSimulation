        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb 14 12:23:02 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FINEENERGYGRID__genmod
          INTERFACE 
            SUBROUTINE FINEENERGYGRID(STATEBASIS,TCSBASIS,ENFINE,NENFINE&
     &)
              USE TOTALCS_MODULE
              USE STATE_CLASS
              TYPE (BASIS_STATE), INTENT(IN) :: STATEBASIS
              TYPE (BASIS_TOTALCS), INTENT(IN) :: TCSBASIS
              REAL(KIND=8) :: ENFINE(10000)
              INTEGER(KIND=4) :: NENFINE
            END SUBROUTINE FINEENERGYGRID
          END INTERFACE 
        END MODULE FINEENERGYGRID__genmod
