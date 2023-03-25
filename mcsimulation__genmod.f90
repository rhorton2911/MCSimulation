        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar  3 16:31:45 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MCSIMULATION__genmod
          INTERFACE 
            SUBROUTINE MCSIMULATION(DATA_INPUT,STATEBASIS,SDCSBASIS,    &
     &ELDCS,DICSBASIS,VARPS,DATAMC,MIN_EIN)
              USE MC
              USE PS_MODULE
              USE TOTALCS_MODULE
              USE DCS_MODULE
              USE SDCS_MODULE
              USE STATE_CLASS
              USE INPUT_DATA
              TYPE (INPUT), INTENT(IN) :: DATA_INPUT
              TYPE (BASIS_STATE), INTENT(IN) :: STATEBASIS
              TYPE (BASIS_SDCS) :: SDCSBASIS
              TYPE (BASIS_DCS), INTENT(IN) :: ELDCS
              TYPE (BASIS_DICS), INTENT(IN) :: DICSBASIS
              TYPE (PSVAR), INTENT(IN) :: VARPS
              TYPE (SIMDATA) :: DATAMC
              REAL(KIND=8) :: MIN_EIN
            END SUBROUTINE MCSIMULATION
          END INTERFACE 
        END MODULE MCSIMULATION__genmod
