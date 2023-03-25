        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar  3 16:31:42 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_CSATEIN__genmod
          INTERFACE 
            SUBROUTINE GET_CSATEIN(STATEBASIS,EINCIDENT,TCS)
              USE TOTALCS_MODULE
              USE STATE_CLASS
              TYPE (BASIS_STATE), INTENT(IN) :: STATEBASIS
              REAL(KIND=8), INTENT(IN) :: EINCIDENT
              TYPE (TOTALCS), INTENT(OUT) :: TCS
            END SUBROUTINE GET_CSATEIN
          END INTERFACE 
        END MODULE GET_CSATEIN__genmod
