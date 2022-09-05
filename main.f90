program main

  use numbers
  use input_data         ! definitions for input data type, construct routine, object data_in
  use state_class        ! defines states (and basis of them) with operations on them
  use totalcs_module     ! reading totalcs files
  use sdcs_module !Handles operations on sdcs
  use Ps_module			 ! handles all Ps Formation Benchmark subroutines
  use dcs_module         ! deals with elastic DCS
  use mc



 
  implicit none
 
  logical:: ex, useState
  integer:: nfile, n
  type(basis_totalcs):: tcsbasis
  character(len=40):: filename
  integer:: Nmaxall
  type(basis_state):: statebasis
  type(basis_dcs):: eldcs    ! keep elastic DCS here
  type(basis_dics):: dicsBasis !Stores dissociative ionisation CS
  type(basis_sdcs):: sdcsBasis !Stores single differential CS
  type(simdata), dimension(:), allocatable:: datamcArray    !stores results of each run of mcsimulation
  real(dp), dimension(10000):: enfinetmp
  real(dp), dimension(:), allocatable:: enfine, tempArray
  real(dp), allocatable, dimension(:):: TICS
  integer:: Nenfine, jj, numInel
  real(dp):: eIonBenchmark, intVal, deltaE, cutoffEn
  real(dp):: min_ein
  real(dp):: percentage
  real(dp):: pEl, pIon   !Scaling percentages in error analysis
  real(dp), dimension(:), allocatable:: pInel
  type(totalcs):: tcs
  ! Random Number Generation
  integer :: values(1:8), k, ii
  integer, dimension(:), allocatable :: seed

  type(PsVar):: VarPs ! For Ps Benchmark Simulation
  
  ! Initialise Random Number Seed
  call date_and_time(values=values)
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = values(8)
  call random_seed(put=seed)

  !------------------------------------------------------------------------
  ! Input data: construct routine for input type.
  inquire(file='data.in',exist=ex)
  if (ex) then
     nfile = 10 
     print*, 'data.in found'
  else
     print*,'data.in not found, will use standard input'
     nfile = 5  
  endif
  call readin( data_in, 10,1 )
  !-------------------------------------------------------------------------

  call initialisePsVar(VarPs)                   ! Initialise simulation input variables for Ps Benchmark Simulation
  if(data_in%benchmark) then                    ! True if Ps Benchmark Simulation is to be done
     print*, 'Ps Benchmark Mode Selected'
     Nmaxall = 4                                ! There are 4 possible options (states): Ps,el,exc,ion
     call new_basis(statebasis,Nmaxall,3)       ! Creates an array of 4 states, tcstype hardcoded to 3	 	
     call populatePsstates(statebasis,VarPs)    ! define states
     min_ein = statebasis%b(2)%enex
  else if (data_in%benchmarkOp .eq. 1) then
     !Use cross sections presented in Garvey and Green (1976) to benchmark
     !the simulation.
     print*, "Electron-H2 Benchmark Mode Selected"

     !Read all totalcs files and set up states, as their data is needed to create enfine 
     call new_totalcsbasis(tcsbasis,data_in%Ntcs,data_in%tcstype)
     do n=1,data_in%Ntcs
        filename = TRIM(data_in%filename_tcs(n))
        print*,'filename=', filename
        call  read_totalcs(tcsbasis%b(n),data_in%tcstype,filename)
        print*, 'tcs%Nmax = ', tcsbasis%b(n)%Nmax
	!    call print_tcs(tcsbasis%b(n))
     enddo
     Nmaxall = tcsbasis%b(data_in%Ntcs)%Nmax
     call new_basis(statebasis,Nmaxall,data_in%tcstype)
     ! define states
     call populatestates(statebasis,tcsbasis)

     call fineenergygrid(statebasis,tcsbasis,enfinetmp,Nenfine)
     print*,  'Nenfine =', Nenfine 
     allocate(enfine(Nenfine))
     enfine(1:Nenfine) = enfinetmp(1:Nenfine)
     do n=1,-Nenfine
        print*,  '-->> enfine() =', n, enfine(n)
     enddo
  
     !Clear tcsbasis and statebasis
     call destruct_totalcsbasis(tcsbasis)
     !call destruct_basis(statebasis)

     !Construct SDCS using analytical formula, considering only one
     !ionisation continuum.
     call createAnalyticSdcs(sdcsBasis) 
     !Calculate TICS as a function of incident energy
     allocate(TICS(sdcsBasis%Nein)) 
     eIonBenchmark = 16.0 
     do ii = 1, sdcsBasis%Nein
        intVal = 0.0
        cutoffEn = (sdcsBasis%ein(ii)-eIonBenchmark)/2.0
        do jj = 1, sdcsBasis%Ncgsdcs-1
           if ((sdcsBasis%ecg(jj) .lt. cutoffEn) .and. (sdcsBasis%ein(ii) .gt. eIonBenchmark)) then
              deltaE = sdcsBasis%ecg(jj+1) - sdcsBasis%ecg(jj)
              intVal = intVal + 0.5*(sdcsBasis%cs(jj,ii)+sdcsBasis%cs(jj+1,ii))*deltaE
           end if
        end do
        TICS(ii) = intVal
     end do 

     !Read in table of benchmark fit parameters for H2 excitations and 
     !single state treatment of dissociative ionisation.
     !Also include approximate CS for v=1, v=2 ground state
     !excitations.
     min_ein = statebasis%b(2)%enex
     print*, "BEFORE: ", statebasis%b(2)%enex
     allocate(tempArray(sdcsBasis%Nein))
     tempArray(:) = sdcsBasis%ein(:)
     call readTcsBenchmark(tcsbasis, statebasis, sdcsBasis%Nein, tempArray, TICS)
     deallocate(tempArray)     
     print*, "AFTER: ", statebasis%b(2)%enex
     !Issue: simulation.f90 cuts energy off at statebasis%b(2)%enex. This is set to lower
     !       than lowest point on fine energy grid just above. This causes the selectAngle
     !       subroutine to crash. Need to change it so that enfine starts at altered value
     !       of statebasis%b(2)%enex. 


     !Read in dcs, these do not affect results to benchmark.
     !call make_elastic_dcs(eldcs,size(vcsBasis%ein),vcsBasis%ein)
     call make_elastic_dcs(eldcs,Nenfine,enfine)

     !Fills an array of dcs type in the basis_dcs defined over a fine
     !energy grid. Also finds average scattering angle.
     call fillDcsBasis(eldcs)
  else 
     print*, 'Default Mode Selected'

     ! Read all totalcs files	 
     call new_totalcsbasis(tcsbasis,data_in%Ntcs,data_in%tcstype)
     do n=1,data_in%Ntcs
        filename = TRIM(data_in%filename_tcs(n))
        print*,'filename=', filename
        call  read_totalcs(tcsbasis%b(n),data_in%tcstype,filename)
        print*, 'tcs%Nmax = ', tcsbasis%b(n)%Nmax
	!    call print_tcs(tcsbasis%b(n))
     enddo

     !Set-up the states data structure - from the totalcs file with the largest energy
     Nmaxall = tcsbasis%b(data_in%Ntcs)%Nmax
     call new_basis(statebasis,Nmaxall,data_in%tcstype)
     ! define states
     call populatestates(statebasis,tcsbasis)
     if (data_in%vcsOp .eq. 1) then
        !Read in vcs for elastic and excitation collisions
        call populateStatesVcs(statebasis,tcsbasis)
        call readVcsPseudo(statebasis)
     end if

		 useState = .false.
     if (data_in%numStatesIn .ne. 0) then
				do ii = 1, statebasis%n
			     useState = .false.
					 jj=1
			     do while ((useState .eqv. .false.) .and. (jj .le. data_in%numStatesIn))
				      if (statebasis%b(ii)%stlabel .eq. trim(data_in%statesToUse(jj))) then
								 print*, statebasis%b(ii)%stlabel
                 useState = .true.
					    end if
              jj = jj+1
			     end do
					   
					 if (useState .eqv. .false.) then
							statebasis%b(ii)%cs(:) = 0.0_dp
				   end if
	      end do
		 end if

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        <<<<<<<<-----------------
     !!!!!!!!!!!!!!!!!!!!!!!!!Return to this for testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        <<<<<<<<-----------------
     !Cross section interpolation for ionising pseudostates gives negative values above 500 ev!!!
     !These values come from the totalcs_... (no .txt) files in DATA_1

     !do ii=1, statebasis%n
     !   print*, statebasis%b(ii)%stlabel, statebasis%b(ii)%v
     !   write(6,*) (statebasis%b(ii)%cs(jj), jj=1, size(statebasis%b(ii)%ein)) 

     !   if (statebasis%b(ii)%stlabel(1:3) .eq. 'ION') then
     !      !do jj = 1, size(statebasis%b(ii)%ein)
     !      !   if ((statebasis%b(ii)%cs(jj) .lt. 0.0d0) .and. (statebasis%b(ii)%ein(jj) .lt. 500.0d0)) then
     !      !      print*, "En, CS: ", statebasis%b(ii)%ein(jj), statebasis%b(ii)%cs(jj)
     !      !   end if
     !      !end do
     !   end if
     !end do
     !stop
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		 print*, "MAKE FINE ENERGY GRID"

     call fineenergygrid(statebasis,tcsbasis,enfinetmp,Nenfine)
     print*,  'Nenfine =', Nenfine 
     allocate(enfine(Nenfine))
     enfine(1:Nenfine) = enfinetmp(1:Nenfine)
     do n=1,-Nenfine
        print*,  '-->> enfine() =', n, enfine(n)
     enddo

     call  intp_cs(statebasis,Nenfine,enfine)
     call make_elastic_dcs(eldcs,size(statebasis%b(1)%ein),statebasis%b(1)%ein)
     
     !Fills an array of dcs type in the basis_dcs defined over a fine
     !energy grid. Also finds average scattering angle.
     call fillDcsBasis(eldcs)

     !Reads in inelastic differential cross sections.
     call make_inelastic_dcs(eldcs,statebasis)

     !Reads in dissociative ionisation cross sections
     call readDics("../data-1/DissIonTotal.csv", dicsBasis)

     !Creates a set of dissociation fractions for dissociative
     !ionisation defined over a fine energy grid. 
     call fillDicsBasis(statebasis ,dicsBasis, Nenfine, enfine)

     !Reads in single differential cross section (SDCS) data.
     call readSdcs(data_in,sdcsBasis)

     !Interpolates SDCS to two dimensional energy grid    
     call intpSdcsModFine(data_in,sdcsBasis,Nenfine,enfine)

     if (data_in%sdcsScaleOp .eq. 1) then
        call scaleSdcs(sdcsBasis,stateBasis)
     end if
     if (data_in%distortSdcsOp .eq. 1) then
        call distortSdcs(sdcsBasis,stateBasis)
     end if

     !Construct analytical SDCS distribution
     call constructSDCSDist(sdcsBasis)

     print*
     print*, 'Finished data proccessing'
     print*

     if (data_in%debugOp .eq. 1) then
        print*, 'Testing Input Data'
        call testInputData(statebasis, eldcs, sdcsBasis, enfine, Nenfine)
     end if
     min_ein = statebasis%b(2)%enex
  end if

  !Get number of inelastic cross sections
  numInel = 0
  do ii = 2, statebasis%n
     if (statebasis%b(ii)%en .lt. 0.0d0) then
        numInel = numInel + 1
     end if
  end do


  pEl = 0.0_dp
  pIon = 0.0_dp
  deallocate(seed)
  allocate(datamcArray(data_in%totalVariations))
  do ii = 1, data_in%totalVariations
     !____________Scaling inputs for uncertainty propagation___________!
     !-finds a percentage (positive or negative) to (increase/decrease) the
     ! magnitude of the cross sections by. Percentage taken from normal 
     ! distribution with mean 0% and standard deviation taken to be the 
     ! estimated uncertainty in each type of cross section.
     !-runs full simulation many times with different scaled inputs, spread of
     ! outputs can be used to determine uncertainty in converged output values of the simulation
     ! without knowing the explicit functional relationship between inputs and outputs.
     ! This is the 'total monte carlo' method of uncertainty propagation.

     !Use a different seed for each variation of the simulation to avoid 
     !introducing spurious correlations between the simulations. This is 
     !important for error analysis.

     call date_and_time(values=values)
     call random_seed(size=k)
     allocate(seed(1:k))
     seed(:) = values(8)
     call random_seed(put=seed)
 
     if (data_in%totalVariations .gt. 1) then
        percentage = 0.0d0
        call sampleNormal(percentage, 0.0d0, 5.0d0) 
        pEl = percentage
        !pEl = 5.0d0
        print*, percentage
        if (data_in%elScale .eq. 1) then
           call scaleElCs(statebasis,percentage)
        end if 

        if (data_in%inelScale .eq. 1) then
           allocate(pInel(numInel))
           call sampleNormal(percentage, 0.0d0, 5.0d0) 
           do jj =2, statebasis%n
              if (statebasis%b(jj)%en .lt. 0.0d0) then
                 !call sampleNormal(percentage, 0.0d0, 11.0d0) 
                 pInel(jj-1) = percentage
                 !print*, "pInel, jj-1: ", percentage, jj-1
                 call scaleInelCs(statebasis,percentage,jj)
              end if
           end do
        end if
    
        call sampleNormal(percentage, 0.0d0, 5.0d0) 
        pIon = percentage
        print*, percentage
        if (data_in%ionScale .eq. 1) then
           call scaleIonCs(statebasis,percentage)
        end if
    end if

    !----------------------------  Run MonteCarlo Simulation ----------------
    !Each variation is run data_in%totalSims/data_in%totalVariations times.
     call mcsimulation(data_in,statebasis,sdcsBasis,eldcs,dicsBasis,VarPs,datamcArray(ii),min_ein)
    !----------------------------  END MonteCarlo Simulation ----------------

    !-Restore input data to what it was before modification-----------------
    !Scaling back down is given by scaling by percentage p' = -p/(1+p/100), where
    !p is the original percentage (e.g. p=5.0) scaled up.
    if (data_in%totalVariations .gt. 1) then
       if (data_in%elScale .eq. 1) then
          call scaleElCs(statebasis, -pEl/(1.0d0 + (pEl/100.0d0)))
       end if
       if (data_in%inelScale .eq. 1) then
          do jj = 1, numInel 
             call scaleInelCs(statebasis, -pInel(jj)/(1.0d0 + (pInel(jj)/100.0d0)),jj+1)
          end do
          deallocate(pInel)
       end if
       if (data_in%ionScale .eq. 1) then
          call scaleIonCs(statebasis, -pIon/(1.0d0 + (pIon/100.0d0)))
       end if
    end if 

    deallocate(seed)  
  end do

  !------------Process the results of the data_input%totalVariations number of runs--------!
  if (data_in%totalVariations .gt. 1) then
     print*, "Propagating Errors"
     !Errors are being propagated using fast total monte carlo method
     call calculateErrors(datamcArray,data_in%totalVariations,data_in)
  end if
  if (data_in%totalVariations .eq. 1) then
     open(70,file="variance.txt")
     write(70,*) datamcArray(1)%mssecE - datamcArray(1)%msecE**2
     write(70,*) datamcArray(1)%msexcite - datamcArray(1)%mexcite**2
     write(70,*) datamcArray(1)%msdissociations - datamcArray(1)%mdissociations**2
     write(70,*) datamcArray(1)%msgen - datamcArray(1)%mgen**2
     write(70,*) datamcArray(1)%msinERad - datamcArray(1)%minERad**2
     write(70,*) datamcArray(1)%msB1SuExc - datamcArray(1)%mB1SuExc**2
     write(70,*) datamcArray(1)%msC1PuExc - datamcArray(1)%mC1PuExc**2
     write(70,*) datamcArray(1)%mssingletExc - datamcArray(1)%msingletExc**2
     write(70,*) datamcArray(1)%mstripletExc - datamcArray(1)%mtripletExc**2
     write(70,*) datamcArray(1)%msSingIonPair - datamcArray(1)%mSingIonPair**2
     write(70,*) datamcArray(1)%msTripIonPair - datamcArray(1)%mTripIonPair**2 
     close(70)
  end if
  
  !Clean up memory
  deallocate(enfine) 
  deallocate(datamcArray)
  call delete_totalcs(tcs)   
  call destruct_basis(statebasis)
  call destruct_totalcsbasis(tcsbasis)
  call destruct_basis_dcs(eldcs)
  call destruct_input(data_in)
  call destructDicsBasis(dicsBasis)
  call destructSdcsBasis(sdcsBasis)

  stop
end program main
!
!
subroutine fineenergygrid(statebasis,tcsbasis,enfine,Nenfine)

 use state_class        ! defines states (and basis of them) with operations on them
 use totalcs_module     !  reading totalcs files

 implicit none
 type(basis_state), intent(in):: statebasis
 type(basis_totalcs), intent(in):: tcsbasis
 real(dp), dimension(10000):: enfine
! real(dp), dimension(:), pointer:: enfine
 integer:: Nenfine

 real(dp):: ein_min, ein_max, en, en1, step1, step2, estart
 integer:: i
 
 !  fine energy mesh
 !  lowest energy, largest energy
 !  step: 0.1 eV up to 25 eV
 !  step: 1 eV up to largest energy
 
 print*
 print*, 'Make a fine energy grid'

 ein_min = statebasis%b(2)%enex
 ein_max = tcsbasis%b(tcsbasis%n)%en_incident
 print*, 'ein_min, ein_max', ein_min, ein_max
 en = 0d0
 step1 = 0.1
 en1 = 25.0
 step2 = 1.0

 do 
    if(en .gt. ein_min) exit
    en = en + step1
 enddo
 estart = en - step1

 en = estart - step1 
 i = 0
 do 
    if(en .gt. en1) exit
    en = en + step1
    i = i + 1
    enfine(i) = en
!    print*,'enfine=', i, en
 enddo
 do 
    if(en .gt. ein_max) exit
    en = en + step2
    i = i + 1
    enfine(i) = en
!    print*,'--->>> enfine=', i, en
 enddo
 Nenfine = i

end subroutine fineenergygrid




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: testInputData
!Purpose: runs checks on the input data needed for debugging
!Date last modified: 11/02/2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testInputData(statebasis, eldcs, sdcsBasis, enfine, Nenfine)
    use input_data         ! definitions for input data type, construct routine, object data_in
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     ! reading totalcs files
    use sdcs_module !Handles operations on sdcs
    use Ps_module			 ! handles all Ps Formation Benchmark subroutines
    use dcs_module         ! deals with elastic DCS
    use mc

    implicit none

    integer:: i,n, nst                       
    character(len=40):: filename                   
    character(len=10):: filenamest, enString       
    character(len=4):: stnumber                    
    integer:: Nenfine
    type(basis_state):: statebasis                                   
    type(basis_dcs):: eldcs    ! keep elastic DCS here               
    type(basis_sdcs):: sdcsBasis !Stores single differential CS      
    real(dp), dimension(Nenfine), intent(in):: enfine                       
    real(dp), dimension(12)::enArray                                   
    integer:: ii, jj, enSize 
    real(dp):: en_incident, dissEn 
    type(totalcs):: tcs        
    type(dcs):: dcsaten       
    type(sdcs):: sdcsaten      
    logical:: writeData       
    real(dp), dimension(:), allocatable:: boundX1SgCs, boundB1SuCs, boundC1PuCs, bounda3SgCs, boundc3PuCs, boundEF1SgCs, totalC1PuDiss
     real(dp), dimension(:), allocatable:: totalB1SuDiss,totalb3SuDiss,totala3SgDiss,totalc3PuDiss,totalEF1SgDiss
    real(dp), dimension(:,:), allocatable:: forWriting 

     nst = 3
     do n=1,-statebasis%b(nst)%inum
        print*, n, statebasis%b(nst)%ein(n),statebasis%b(nst)%cs(n)
     enddo

     !------- Tests -----
     i = 1
     if(i .eq. 1) then
        print*, 'test 1'
        en_incident = enfine(1) + 2.0
        call get_csatein(statebasis,en_incident,tcs)
        print*, 'test 2'
        !call print_tcs(tcs)
     endif
     !--------------------------
     i = 1
     if(i .eq. 1) then
        print*, 'test 3'
        do n=1,2 !statebasis%n
           !      if(statebasis%b(n)%idf .ne. 0) then
           write(stnumber,'(i4.4)') n
           filenamest ="st_"//stnumber
           !call print_state(statebasis%b(n),filenamest)
           !      endif
        enddo
     endif
     !---------------------------
     i = 1
     if(i .eq. 1) then
        print*, 'test 4'
        filename = 'dcs_20.0eV'
        en_incident = 21.0
        call get_dcsatein(eldcs,en_incident,dcsaten,1)  !Testing elastic dcs
        !call print_dcs(dcsaten,filename)
     endif
     !----------------------------
     !Write totalcs to file for plotting 
     call get_csatein(stateBasis, 50.0_dp, tcs)
     call write_tcs(tcs)
     !--------------------------- 
     !Write totalcs and symmetrised totalcs to file
     !for comparison at a particular incident energy.
     en_incident = 50.0 !eV
     call get_csatein(statebasis,en_incident,tcs)
     call write_tcs(tcs)
     call symTcs(tcs,en_incident)
     !Write incident energy to a string for use in filename
     write(enString,'(I0)') int(tcs%en_incident)

     filename = 'tcs_sym_'//TRIM(enString)//'eV.txt'
     open(70, file=filename)
     write(70, *) "energy(eV)        CS(a.u)"
     do ii = 1, tcs%Nmax
        write(70,*) tcs%en(ii),  tcs%CS(ii)
     end do     
     close(70)

     !Write total excitation CS to file for plotting.
     call writeExcitationCs(stateBasis)
     call writeFullCs(stateBasis)
     call writeElCs(stateBasis) 
     call writeElCsInState(stateBasis)


     !Test sdcs correct
     call get_sdcsatein(sdcsBasis,500.0_dp, sdcsaten)
     print*, "Test SDCS"
     !call print_sdcs(sdcsaten)
     print*, "End Test SDCS"
     call get_sdcsatein(sdcsBasis, 300.0_dp, sdcsaten)

     !Compare integrated SDCS and totalcs as a consistency check
     call testSdcs(sdcsBasis,enfine(1))
     call testTcs(stateBasis,enfine,Nenfine)
     !call printSdcsAtEn(sdcsaten) 
     call destructSdcsAtEin(sdcsaten)

     !Print SDCS to file at several different energies
     enArray = (/16.0_dp,18.0_dp,20.0_dp,30.0_dp,40.0_dp,50.0_dp,100.0_dp,150.0_dp,200.0_dp,300.0_dp,400.0_dp,500.0_dp/)
     do ii=1, size(enArray)
        call get_sdcsatein(sdcsBasis,enArray(ii),sdcsaten)
        call printSdcsAtEn(sdcsaten)
     end do
     call get_sdcsatein(sdcsBasis, 300.0_dp, sdcsaten)

     !Test analytical SDCS distribution
     en_incident = 100.0_dp
     call get_csatein(stateBasis,en_incident,tcs)
     call getSdcsDistAtEin(sdcsBasis,en_incident,sdcsaten,tcs)
     call printSdcsDistAtEn(sdcsaten) 

     !Print SDCS and Distribution treatment total probabilities
     !to a file for plotting
     call testTics(sdcsBasis,stateBasis)

     !Print SDCS to a single file for all incident energies for creation
     !of a 3D plot, used for testing and debugging.
     call printSDCSAll(sdcsBasis,enfine(1))
     !call printIntpSdcs(sdcsBasis) 
 
    !----------Comparison of SDCS with totalCS treatments------------
     !Compare ejection probabilities, cross sections and mean 
     !ejection energy of different ionisation treatments
     en_incident = 300_dp
     call calcProbs(sdcsBasis,stateBasis,15.0_dp)
     if (data_in%stateIonop .eq. 1) then
        call symTcs(tcs,en_incident)              
     end if    

     call get_sdcsatein(sdcsBasis, 300.0_dp, sdcsaten)

     !______________________Test vibrationally resolved totalcs_______
     if (data_in%vcsOp .eq. 1) then
        !Sum all CS in totalcs and compare to regular total cross section
        call testTcsVib(stateBasis)
     end if

     if (data_in%vcsop .eq. 1) then
        !Test dissociation fractions by printing them to file for given states
        !Compare these files with input data

        do ii = 1, statebasis%n
           if ((statebasis%b(ii)%stlabel .eq. 'X1Sg') .or. (statebasis%b(ii)%stlabel .eq. 'a3Sg') &
              .or. (statebasis%b(ii)%stlabel .eq. 'B1Su') .or. (statebasis%b(ii)%stlabel .eq. 'C1Pu') & 
              .or. (statebasis%b(ii)%stlabel .eq. 'b3Su') .or. (statebasis%b(ii)%stlabel .eq. 'c3Pu') &
              .or. (statebasis%b(ii)%stlabel .eq. 'EF1Sg')) then

              if (statebasis%b(ii)%v .ge. 0) then
                 filename = trim(statebasis%b(ii)%stlabel)//'DFIntp'
                 if (statebasis%b(ii)%v .eq. 0) then
                    open(70, file=trim(filename), position='append')
                    write(70,*) 'v       DF(a.u)'
                 else 
                    open(70, file=trim(filename), position='append')
                 end if
                 !print*, statebasis%b(ii)%v, statebasis%b(ii)%df(5), statebasis%b(ii)%stlabel   
                 write(70,*) statebasis%b(ii)%v, statebasis%b(ii)%df(5)    !DF constant across energy range

                 close(70)
              end if
           end if
        end do 


        enSize = size(statebasis%b(1)%ein)
        allocate(boundX1SgCs(enSize)) 
        allocate(boundB1SuCs(enSize))
        allocate(boundC1PuCs(enSize))
        allocate(bounda3SgCs(enSize))
        allocate(boundc3PuCs(enSize))
        allocate(boundEF1SgCs(enSize))
        boundX1SgCs(:) = 0.0d0
        boundB1SuCs(:) = 0.0d0
        boundC1PuCs(:) = 0.0d0
        bounda3SgCs(:) = 0.0d0
        boundc3PuCs(:) = 0.0d0
        boundEF1SgCs(:) = 0.0d0

        !Print excitation cross sections to file for comparison with input data
        writeData = .false.
        do ii = 1, statebasis%n
           if ((statebasis%b(ii)%stlabel .eq. 'X1Sg') .or. (statebasis%b(ii)%stlabel .eq. 'a3Sg') &
              .or. (statebasis%b(ii)%stlabel .eq. 'B1Su') .or. (statebasis%b(ii)%stlabel .eq. 'C1Pu') & 
              .or. (statebasis%b(ii)%stlabel .eq. 'b3Su') .or. (statebasis%b(ii)%stlabel .eq. 'c3Pu') &
              .or. (statebasis%b(ii)%stlabel .eq. 'EF1Sg')) then
              if ((statebasis%b(ii)%v .ge. 0) .and. (statebasis%b(ii)%v .le. 12)) then 
                 write(enString,'(I0)') statebasis%b(ii)%v
                 filename = 'cs'//trim(statebasis%b(ii)%stlabel)//'vf'//trim(enString)//'Intp'
                 writeData = .true.
              else if (statebasis%b(ii)%v .eq. -1) then
                 filename = 'cs'//trim(statebasis%b(ii)%stlabel)//'DissIntp'
                 writeData = .true.
              end if     
     
              if (writeData .eqv. .true.) then
                 open(70, file=trim(filename))
                 write(70,*) 'Energy(eV),       CS(a.u)'
                 do jj = 1, statebasis%b(ii)%inum
                    write(70,*) statebasis%b(ii)%ein(jj), statebasis%b(ii)%cs(jj)  
                 end do

                 close(70)
              end if
              writeData = .false.
 
              dissEn= statebasis%b(ii)%enex - statebasis%b(ii)%dissThresh
              if ((dissEn .lt. 0.0d0) .or. (statebasis%b(ii)%stlabel .eq. 'X1Sg')) then
                 if (statebasis%b(ii)%stlabel .eq. 'X1Sg') then 
                    boundX1SgCs(:) = boundX1SgCs(:) + statebasis%b(ii)%cs(:)
                 else if (statebasis%b(ii)%stlabel .eq. 'B1Su') then
                    boundB1SuCs(:) = boundB1SuCs(:) + statebasis%b(ii)%cs(:)
                 else if (statebasis%b(ii)%stlabel .eq. 'C1Pu') then
                    boundC1PuCs(:) = boundC1PuCs(:) + statebasis%b(ii)%cs(:)
                 else if (statebasis%b(ii)%stlabel .eq. 'a3Sg') then
                    bounda3SgCs(:) = bounda3SgCs(:) + statebasis%b(ii)%cs(:)
                 else if (statebasis%b(ii)%stlabel .eq. 'c3Pu') then
                    boundc3PuCs(:) = boundc3PuCs(:) + statebasis%b(ii)%cs(:)
                 else if (statebasis%b(ii)%stlabel .eq. 'EF1Sg') then
                    boundEF1SgCs(:) = boundEF1SgCs(:) + statebasis%b(ii)%cs(:)
                 end if 
              end if 
           end if
        end do
        allocate(forWriting(enSize,6)) 
        forWriting(:,1) = boundX1SgCs(:)
        forWriting(:,2) = boundB1SuCs(:)
        forWriting(:,3) = boundC1PuCs(:)
        forWriting(:,4) = bounda3SgCs(:)
        forWriting(:,5) = boundc3PuCs(:)
        forWriting(:,6) = boundEF1SgCs(:)

        open(70,file='BoundTotal.txt')
        write(70,*) "Ein(ev)               X1Sg                 B1Su                C1Pu                a3Sg                 c3Pu                 EF1Sg"
           do ii =1, enSize 
              write(70,*) statebasis%b(1)%ein(ii), (forWriting(ii,jj), jj=1, 6) 
           end do
        close(70)


	!Test total dissociative excitation cross sections for specific states of interest, calculate 
	!by summing vibrational pseudostates and writing to file for plotting
	allocate(totalC1PuDiss(size(statebasis%b(1)%ein)))
	allocate(totalB1SuDiss(size(statebasis%b(1)%ein)))
	allocate(totalb3SuDiss(size(statebasis%b(1)%ein)))
	allocate(totala3SgDiss(size(statebasis%b(1)%ein)))
	allocate(totalc3PuDiss(size(statebasis%b(1)%ein)))
	allocate(totalEF1SgDiss(size(statebasis%b(1)%ein)))
	totalC1PuDiss(:)=0.0_dp
	totalB1SuDiss(:)=0.0_dp
	totalb3SuDiss(:)=0.0_dp
	totala3SgDiss(:)=0.0_dp
	totalc3PuDiss(:)=0.0_dp
	totalEF1SgDiss(:)=0.0_dp

	if (data_in%vcsOp .eq. 1) then
	   do ii =1, statebasis%n
	      if ((statebasis%b(ii)%enex - statebasis%b(ii)%dissThresh) .gt. 0.0_dp) then 
	         if (statebasis%b(ii)%stlabel .eq. 'C1Pu') then
	            totalC1PuDiss(:) = totalC1PuDiss(:) + statebasis%b(ii)%cs(:) 
		 end if
	         if (statebasis%b(ii)%stlabel .eq. 'B1Su') then
	            totalB1SuDiss(:) = totalB1SuDiss(:) + statebasis%b(ii)%cs(:) 
		 end if
	         if (statebasis%b(ii)%stlabel .eq. 'b3Su') then
	            totalb3SuDiss(:) = totalb3SuDiss(:) + statebasis%b(ii)%cs(:) 
		 end if
	         if (statebasis%b(ii)%stlabel .eq. 'a3Sg') then
	            totala3SgDiss(:) = totala3SgDiss(:) + statebasis%b(ii)%cs(:) 
		 end if
	         if (statebasis%b(ii)%stlabel .eq. 'c3Pu') then
	            totalc3PuDiss(:) = totalc3PuDiss(:) + statebasis%b(ii)%cs(:) 
		 end if
	         if (statebasis%b(ii)%stlabel .eq. 'EF1Sg') then
	            totalEF1SgDiss(:) = totalEF1SgDiss(:) + statebasis%b(ii)%cs(:) 
		 end if
	      end if 
	   end do
	end if

	open(70, file='DissTotal.txt')
	write(70,*) "Total dissociation cross section obtained by suming pseudostates as a function of projectile energy"
	write(70,*) "EIN,C1Pu,B1Su,b3Su,a3Sg,EF1Sg"
	do ii = 1, SIZE(statebasis%b(1)%ein)
	   write(70,*) statebasis%b(1)%ein(ii), totalC1PuDiss(ii), totalB1SuDiss(ii), totalb3SuDiss(ii), totala3SgDiss(ii), totalc3PuDiss(ii), totalEF1SgDiss(ii)
	end do
	close(70)

	!Write dissociation fractions of specific states to file for comparison with input
	open(71,file='C1PuDissFrac.txt')
	open(72,file='B1SuDissFrac.txt')
	if (data_in%vcsop .eq. 1) then
     do ii = 1, statebasis%n
        if (statebasis%b(ii)%stlabel .eq. 'C1Pu' ) then
					 if (statebasis%b(ii)%enex - statebasis%b(ii)%dissThresh .lt. 0.0_dp) then
    	        write(71,*) statebasis%b(ii)%v, statebasis%b(ii)%df(1)
				   end if
        end if 
        if (statebasis%b(ii)%stlabel .eq. 'B1Su' ) then
					 if (statebasis%b(ii)%enex - statebasis%b(ii)%dissThresh .lt. 0.0_dp) then
    	        write(72,*) statebasis%b(ii)%v, statebasis%b(ii)%df(1)
				   end if
        end if 
   	end do
  end if
  close(71)
  close(72)


        deallocate(totalC1PuDiss,totalB1SuDiss,totalb3SuDiss,totala3SgDiss,totalc3PuDiss,totalEF1SgDiss)
        deallocate(forWriting)
        deallocate(boundX1SgCs) 
        deallocate(boundB1SuCs)
        deallocate(boundC1PuCs)
        deallocate(bounda3SgCs)
        deallocate(boundc3PuCs)
        deallocate(boundEF1SgCs) 
     end if

end subroutine testInputData





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: testTics
!Purpose: calculates total ionisation cross section using three 
!         different ionisation treatments and writes them to a file
!         for comparison. Used for debugging pruposes.
!Date last modified: 26/09/2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testTics(sdcsob,stateBasis)
    use state_class
    use totalcs_module
    use sdcs_module

    implicit none
    type(basis_state)::stateBasis
    type(basis_sdcs)::sdcsob
    type(sdcs)::sdcsaten
    type(sdcs)::sdcsatenDist
    type(totalcs)::tcs
    integer::ii,jj,kk
    real(dp)::eIn,ticsSdcs,ticsState,ticsSdcsDist
    real(dp):: change, eIon, cutoffEn


    eIon = 15.96632  !Ionisation energy of ground state H2
    open(22,file="ticsTest.txt")
    write(22,*) "EIn (eV)     TotalCS        SDCS       SDCSDist" 
    do ii = 1, sdcsob%Nein-1
       eIn = sdcsob%ein(ii)

       ticsSdcs = 0.0_dp
       ticsState = 0.0_dp
       ticsSdcsDist = 0.0_dp 

       if (eIn .gt. eIon) then
          call get_csatein(stateBasis,eIn,tcs)
          call getSdcsDistAtEIn(sdcsob,eIn,sdcsatenDist,tcs)
          call get_sdcsatein(sdcsob,eIn,sdcsaten)
    
          cutoffEn = (eIn-eIon)/2.0 
          !Integrate over SDCS from E_min to E/2
          jj=1
          do while (sdcsaten%ecg(jj) .le. cutoffEn)
             change = 0.5*(sdcsaten%CS(jj+1)+sdcsaten%CS(jj))*(sdcsaten%ecg(jj+1)-sdcsaten%ecg(jj))
             ticsSdcs = ticsSdcs + change        
             jj=jj+1 
          end do 

          !Integrate over SDCS from E_min to E/2
          jj=1
          do while (sdcsatenDist%ecg(jj) .le. cutoffEn)
             change = 0.5*(sdcsatenDist%CS(jj+1)+sdcsatenDist%CS(jj))*(sdcsatenDist%ecg(jj+1)-sdcsatenDist%ecg(jj))
             ticsSdcsDist = ticsSdcsDist + change        
             jj=jj+1 
          end do 

          do kk =1, tcs%Nmax
             if (tcs%en(kk) .gt. 0.0) then
                ticsState = ticsState + tcs%cs(kk)
             end if
          end do
       end if


       write(22,*) ein, ticsState, ticsSdcs, ticsSdcsDist 
       call destructSdcsAtEIn(sdcsaten)
       call destructSdcsAtEIn(sdcsatenDist)
       call delete_totalcs(tcs)
    end do
    close(22)

end subroutine testTics
