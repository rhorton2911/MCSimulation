        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb 14 12:23:02 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TESTINPUTDATA__genmod
          INTERFACE 
            SUBROUTINE TESTINPUTDATA(TCSBASIS,STATEBASIS,ELDCS,SDCSBASIS&
     &,ENFINE,NENFINE)
              USE SDCS_MODULE
              USE DCS_MODULE
              USE STATE_CLASS
              USE TOTALCS_MODULE
              INTEGER(KIND=4) :: NENFINE
              TYPE (BASIS_TOTALCS) :: TCSBASIS
              TYPE (BASIS_STATE) :: STATEBASIS
              TYPE (BASIS_DCS) :: ELDCS
              TYPE (BASIS_SDCS) :: SDCSBASIS
              REAL(KIND=8), INTENT(IN) :: ENFINE(NENFINE)
            END SUBROUTINE TESTINPUTDATA
          END INTERFACE 
        END MODULE TESTINPUTDATA__genmod
