        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar  3 16:31:45 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COLLISIONSIMULATION__genmod
          INTERFACE 
            SUBROUTINE COLLISIONSIMULATION(DATA_INPUT,STATEBASIS,       &
     &SDCSBASIS,DICSBASIS,PARTICLEBASIS,MINEXCEN,PARTNUM,COLL,SIMINDEX, &
     &DATASIM,IONOP,BMODE,VARPS,ELDCS)
              USE DCS_MODULE
              USE PS_MODULE
              USE MC
              USE TOTALCS_MODULE
              USE SDCS_MODULE
              USE STATE_CLASS
              USE INPUT_DATA
              TYPE (INPUT), INTENT(IN) :: DATA_INPUT
              TYPE (BASIS_STATE), INTENT(IN) :: STATEBASIS
              TYPE (BASIS_SDCS) :: SDCSBASIS
              TYPE (BASIS_DICS), INTENT(IN) :: DICSBASIS
              TYPE (PARTICLE) :: PARTICLEBASIS(0:1000)
              REAL(KIND=8) :: MINEXCEN
              INTEGER(KIND=4), INTENT(IN) :: PARTNUM
              INTEGER(KIND=4), INTENT(IN) :: COLL
              INTEGER(KIND=8) :: SIMINDEX
              TYPE (SIMDATA), INTENT(INOUT) :: DATASIM
              CHARACTER(LEN=60), INTENT(IN) :: IONOP
              LOGICAL(KIND=4), INTENT(IN) :: BMODE
              TYPE (PSVAR), INTENT(IN) :: VARPS
              TYPE (BASIS_DCS), INTENT(IN) :: ELDCS
            END SUBROUTINE COLLISIONSIMULATION
          END INTERFACE 
        END MODULE COLLISIONSIMULATION__genmod
