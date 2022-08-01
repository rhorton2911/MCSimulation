        !COMPILER-GENERATED INTERFACE MODULE: Mon Feb 14 11:39:58 2022
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
