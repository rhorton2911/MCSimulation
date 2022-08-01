module mc 
    use numbers

    public::  update_energy,init_simdata,init_particle
	
    type, public:: particle
        real(dp),dimension(0:30000):: energy ! The energy at each collision
	real(dp),dimension(0:30000):: enArr ! The energy target energy at each collision
	integer:: gen ! Which generation the particle was created in
	integer:: colls ! Total number of collisions the particle undergoes
	integer,dimension(0:30000):: state ! The state selected at each collision
	logical,dimension(0:30000):: diss	! True if dissociation occured at that collision
	!real(dp),dimension(0:70000):: costheta ! Cosine of the scattering angle at each collision 
        !real(dp),dimension(0:70000):: phi     !Azimuthal scattering angle randomly generated after each collision	
	logical,dimension(0:30000):: Ps	! True if Ps Formation occured at that collision (only used for Ps Formation Benchmark Simulation)
        !real(dp),dimension(0:70000)::path !Stores path length after each collisions	
                 

        !Stores position and time of particle at each collision. First element is point of creation.
        real(dp),dimension(0:30000)::x
        real(dp),dimension(0:30000)::y
        real(dp),dimension(0:30000)::z
        real(dp),dimension(0:30000)::time
 
        !Stores unit vectors in coordinate system attached to each collision,
        !where the new zHat is along the final particle trajectory.
        real(dp), dimension(3)::xHat
        real(dp), dimension(3)::yHat
        real(dp), dimension(3)::zHat
    end type particle
	
    type, public:: simdata ! Data for the simulation of ONE incident particle
        integer:: secE ! total number of secondary electrons created in this simulation (num of ionisations)		
	integer:: excite ! total number of excitation collisions created in this simulation
	integer:: elastic ! total number of elastic collisions created in this simulation
	integer:: dissociations ! total number of dissociations in this simulation
	integer:: gen ! total number of sec E generations created in this simulation
	integer:: genIndex ! Current generation 
        real(dp):: inERad  !Final radial distance of the incident electron
        integer:: B1SuExc  !total number of excitations of the B1Su state of H2 in this simulation
        integer:: C1PuExc     !total number of excitations of the C1Pu state of H2 in this simulation
        integer:: singletExc  !total number of singlet H2 excitations in this simulation
        integer:: tripletExc !total number of triplet H2 excitations in this simulation
	integer,dimension(1000):: ePerGen ! Number of electrons in each generation
	real(dp),dimension(1000):: enPerGen ! Total energy of electrons in each generation
	integer,dimension(0:1000):: dissPerGen ! Dissociation count per generation
	!integer,dimension(0:1000):: ionPerGen ! number of ionisations in each gen 
	integer,dimension(0:1000):: excPerGen ! number of excitation collisions in each gen
	integer,dimension(0:1000):: elPerGen ! number of elastic collisions in each gen 
	integer,dimension(0:1000):: collPerGen ! number of collisions in each gen 
	       
        !Data recorded for program consistency checks
        real(dp), dimension(0:1000):: ejEnGrid !Grid of ejection energies in ionisation
        integer,dimension(0:1000):: numEj !Number of ionisations at a given energy in above grid 	
				
	! For Ps Benchmark Simulation
	logical:: PsFormed ! True if Positronium is Formed (exit condition for simulation loop)
	integer:: numPsFormed
	real(dp) :: duration,tResTemp,tRes

        !Newer energy deposition parameters
        real(dp):: W   !Energy per ion pair 
        real(dp):: singIonPair  !Singlet excitations per ion pair 
        real(dp):: tripIonPair  !Triplet excitations per ion pair
        real(dp):: B1SuIonPair  !B1Su excitations per ion pair
        real(dp):: C1PuIonPair  !C1Pu excitations per ion pair
        integer:: groundv1 !Number of X1Sg(v=1) excitations
        integer:: groundv2 !Number of X1Sg(v=2) excitations
        real(dp):: groundRatio  !Ratio of X1Sg(v=2) to (v=1) excitations
        

        !Dissociation parameters
        integer:: singDiss !dissociations due to singlet excitations
        integer:: tripDiss !dissociations due to triplet excitations
        integer:: C1PuDiss !dissociations due to C1Pu excitations
        integer:: B1SuDiss !dissociations due to B1Su excitations
        integer:: b3SuDiss !dissociations due to b3Su excitations
        real(dp)::dissHeat  !Amount of kinetic energy released by dissociation
        real(dp)::PDissEn  !Energy of primary electron deposited as heat through dissociation
        real(dp)::PVibEn   !Energy of primary electron lost through vibrational excitation

        !Data used for error analysis, collected over a large number of runs
        !of the simulation
        integer:: numSims   !number of simulations run
        real(dp):: msecE     !mean values of simulation outputs
        real(dp):: mW
        real(dp):: mexcite 
        real(dp):: mdissociations
        real(dp):: mgen 
        real(dp):: minERad 
        real(dp):: mB1SuExc 
        real(dp):: mC1PuExc
        real(dp):: msingletExc  
        real(dp):: mtripletExc  
        real(dp):: mSingIonPair
        real(dp):: mTripIonPair
        real(dp):: mB1SuIonPair
        real(dp):: mC1PuIonPair
  
        !Mean square values of simulation outputs
        real(dp):: mssecE 
        real(dp):: msW
        real(dp):: msexcite 
        real(dp):: msdissociations
        real(dp):: msgen 
        real(dp):: msinERad 
        real(dp):: msB1SuExc
        real(dp):: msC1PuExc 
        real(dp):: mssingletExc 
        real(dp):: mstripletExc 		
        real(dp):: msSingIonPair
        real(dp):: msTripIonPair
        real(dp):: msB1SuIonPair
        real(dp):: msC1PuIonPair
    end type simdata


	
contains
	
    subroutine init_particle(self,initE,gen,xIn,yIn,zIn,tIn) ! Initialise particle values
        implicit none
        type(particle),intent(inout):: self
	real(dp),intent(in):: initE
	integer,intent(in):: gen
        real(dp)::xIn, yIn, zIn, tIn !Initial particle coordinates
		
	self%energy(:) = 0.0
	self%state(:) = 0
	self%diss(:) = .FALSE.
	!self%costheta(:) = 0.0
	self%Ps(:) = .FALSE.
        !self%path(:) = 0.0
			
	self%energy(0) = initE
	self%gen = gen		
	self%colls = 0

        self%x(0) = xIn
        self%y(0) = yIn
        self%z(0) = zIn
        self%time(0) = tIn

        self%xHat(1) = 1.0
        self%xHat(2) = 0.0
        self%xHat(3) = 0.0
        self%yHat(1) = 0.0
        self%yHat(2) = 1.0
        self%xHat(3) = 0.0
        self%zHat(1) = 0.0
        self%zHat(2) = 0.0
        self%zHat(3) = 1.0
    end subroutine init_particle
	
	
    subroutine init_simdata(self) ! Initialise data for ONE simulation
        implicit none
	type(simdata), intent(inout):: self
	integer:: i
        real(dp):: en	
	
	self%secE = 0
	self%gen = 0
	self%excite = 0
	self%elastic = 0
	self%dissociations = 0
		
	self%PsFormed = .false.
	self%numPsFormed = 0
	self%duration = 0.0
	self%tResTemp = 0.0
	self%tRes = 0.0
                
        self%B1SuExc = 0
        self%C1PuExc = 0
        self%singletExc = 0
        self%tripletExc = 0		

        self%singDiss = 0 
        self%tripDiss = 0 
        self%C1PuDiss = 0 
        self%B1SuDiss = 0 
        self%b3SuDiss = 0 
        self%dissHeat = 0.0d0
        self%PDissEn = 0.0d0
        self%PVibEn = 0.0d0

        do i = 1,1000
	   self%ePerGen(i) = 0
	   self%enPerGen(i) = 0.0
	end do
		
	do i = 0,1000
	   self%dissPerGen(i) = 0
	   !self%ionPerGen(i) = 0
	   self%excPerGen(i) = 0
	   self%elPerGen(i) = 0
	   self%collPerGen(i) = 0
	end do

        en = 0.0
        do i =1, 1000
           self%ejEnGrid(i) = en
           en = en + 0.5
        end do
        self%numEj(:) = 0
        self%W = 0.0d0 
        self%singIonPair = 0.0d0
        self%tripIonPair = 0.0d0
        self%B1SuIonPair = 0.0d0
        self%C1PuIonPair = 0.0d0
        self%groundv1 = 0
        self%groundv2 = 0
        self%groundRatio = 0.0d0

        self%msecE = 0.0d0
        self%mW = 0.0d0
        self%mexcite = 0.0d0
        self%mdissociations = 0.0d0
        self%mgen = 0.0d0
        self%minERad = 0.0d0
        self%mB1SuExc = 0.0d0
        self%mC1PuExc = 0.0d0
        self%msingletExc = 0.0d0
        self%mtripletExc = 0.0d0
        self%mSingIonPair = 0.0d0
        self%mTripIonPair = 0.0d0
        self%mB1SuIonPair = 0.0d0
        self%mC1PuIonPair = 0.0d0

        self%mssecE = 0.0d0
        self%msW = 0.0d0
        self%msexcite = 0.0d0
        self%msdissociations = 0.0d0
        self%msgen = 0.0d0
        self%msinERad = 0.0d0
        self%msB1SuExc = 0.0d0
        self%msC1PuExc = 0.0d0
        self%mssingletExc = 0.0d0
        self%mstripletExc = 0.0d0
        self%msSingIonPair = 0.0d0
        self%msTripIonPair = 0.0d0
        self%msB1SuIonPair = 0.0d0
        self%msC1PuIonPair = 0.0d0
    end subroutine init_simdata


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: recordEjEn
    !Purpose: records the ejection energy chosen using one of several
    !         distributions for ionisation.
    !Date last modified: 15/07/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine recordEjEn(energy, datasim)
        implicit none
        type(simdata)::datasim
        real(dp):: energy
        logical:: found
        integer:: ii

        found = .false.
        !Find energy range in grid chosen energy falls into
        
        ii = 1
        do while (.not. found)
           if (energy .gt. datasim%ejEnGrid(ii)) then
              ii = ii + 1
           else 
              found = .true.
           end if 
        end do
            
        !Increment number in a particular range
        datasim%numEj(ii) = datasim%numEj(ii) +1
    end subroutine recordEjEn

	
    subroutine update_energy(statebasis,sdcsBasis,dcsBasis,stateNum,tcs,particlebasis,partNum,cosangle,coll,ionop,enlossop,datasim,bmode,VarPs)
        use state_class
        use totalcs_module
        use Ps_module
        use sdcs_module
        use input_data
        use dcs_module
        implicit none
		
	type(basis_state),intent(in):: statebasis
	type(state):: selstate
        type(basis_sdcs)::sdcsBasis
        type(sdcs)::sdcsAtEIn
	integer, intent(in):: stateNum
	type(particle),dimension(0:1000),intent(inout):: particlebasis
	integer,intent(in):: partNum,coll, enlossop
	character (len=60),intent(in):: ionop ! Ionistation option
	type(totalcs):: tcs ! used to calculate the mean excitation energy
	type(simdata), intent(inout):: datasim		
	real(dp):: enex,en ! excitation and ionisation energy (eV)
	real(dp):: meanexc ! mean excitation energy (eV)
	real(dp):: elEnergyLoss	  ! energy lost by incident particle (eV)
        real(dp):: inElEnergyLoss   !energy loss in inelastic scattering, includes internal energy
        real(dp):: ejEn, secE, primaryE, eIon, ejEnSelected, eIn,cosangle
        logical:: indist
	! For Ps Benchmark Simulation
	type(PsVar),intent(in):: VarPs 
	logical,intent(in):: bmode
        ! For debugging
        integer:: stateNumDebug
        logical:: debug
        type(totalcs)::tcsDebug
        type(basis_dcs)::dcsBasis
        debug = .false.
        
        eIon = 15.96632_dp !Ionisation energy of ground state H2 (eV)	
        eIn = particlebasis(partNum)%energy(coll)

	! The state that has been selected
        !Note that tcs stores energy in a.u, not eV
        en = 0.0_dp
        if (data_in%stateIonop .eq. 1) then
           !Use symmetrised pseudostates
           if (tcs%en(stateNum)*data_in%eV .lt. 0.0_dp) then
              selstate = statebasis%b(stateNum) 
   	      enex =  selstate%enex
    	      en = selstate%en
           else if (tcs%en(stateNum)*data_in%eV .gt. 0.0_dp) then
              en = tcs%en(stateNum)*data_in%eV
           end if
        else if (data_in%stateIonop .eq. 0) then 
           selstate = statebasis%b(stateNum) 
           enex =  selstate%enex
           en = selstate%en
        end if
        call getIndist(bmode,indist) 	
	!print*, 'particle',partNum,'gen',particlebasis(partNum)%gen,'coll',coll
	!print*, 'enex',enex
	!print*, 'preCollEnergy',particlebasis(partNum)%energy(coll)
	
	!print*, 'el cs', tcs%CS(1)
	!print*, 'exc cs', tcs%CS(2)
	!print*, 'ion cs', tcs%CS(3)
	!print*, 'Ps cs', tcs%CS(4)
        particlebasis(partNum)%enArr(coll) = en

	if (stateNum == 1) then ! elastic collision					
	   datasim%elastic = datasim%elastic + 1	! Update total elastic scattering count for this simulation
	   datasim%elPerGen(particlebasis(partNum)%gen) = datasim%elPerGen(particlebasis(partNum)%gen) + 1 ! Update total elastic scattering count for this generation
	   call elasticScattering(particlebasis(partNum)%energy(coll),cosangle,elEnergyLoss,bmode,dcsBasis) 
	   particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - elEnergyLoss 
	
	else if((en < 0.0_dp) .or. (statebasis%b(stateNum)%ion .eqv. .false.)) then ! excitation collision 
	   datasim%excite = datasim%excite + 1 
	   datasim%excPerGen(particlebasis(partNum)%gen) = datasim%excPerGen(particlebasis(partNum)%gen) + 1 ! Update total excitation count for this generation
           call update_H2Data(datasim, stateBasis, stateNum)
			
	   if(bmode) then
	      call elasticScattering(particlebasis(partNum)%energy(coll),cosangle,elEnergyLoss,bmode,dcsBasis) 
	      particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - enex - elEnergyLoss
	   else
              if (enlossop .eq. 1) then
                 !Use only excitation energy of selected state                 
	         particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - enex
              else if (enlossop .eq. 2) then
                 !Use inelastic formula 
                 call inelasticScattering(particlebasis(partNum)%energy(coll),enex,cosangle,inElEnergyLoss,bmode) 	
                 !Formula used for inelastic scattering includes loss
                 !to internal energy.
                 particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - inElEnergyLoss
              end if
	   end if	
	
        else if(((stateNum .gt. 1) .and. (en > 0.0_dp)) .and. (statebasis%b(stateNum)%ion .eqv. .true.)) then ! ionisation collision 
	   datasim%secE = datasim%secE + 1 ! Increment number of secondary electrons created
	   datasim%gen = particlebasis(partNum)%gen + 1 ! Keep track of the number of generations that have been created
           datasim%ePerGen(particlebasis(partNum)%gen + 1) = datasim%ePerGen(particlebasis(partNum)%gen+1) + 1 ! Keep track of number of electrons in the next generation
			
	   if(bmode) then	! Only for Ps Benchmark Simulation
	      particlebasis(partNum)%energy(coll + 1) = VarPs%Q * (particlebasis(partNum)%energy(coll) - VarPs%e_ion)
	      !print*, 'Q',VarPs%Q
	   else ! For Default Simulation
              if(data_in%stateIonop .eq. 1) then
                 !Symmetrise totalcs

                 !tcs object energy grid in Hartees, need to convert
                 enex = tcs%en(stateNum)*data_in%eV
              end if

	      if(ionop.EQ.'mean') then ! mean excitation energy will be used
	         call meanexcenergy(tcs,statebasis,meanexc)

                 !Note; don't need to deal with indistinguishibility yet, as meanexc is
                 !      always below (eIn-eIon)/2.0
                 !ejEnSelected = meanexc
                 !Deal with particle indistinguishibility		
                 !if (indist) then
                 !   if (ejEnSelected .gt. (eIn-eIon)/2.0) then
                    !Above cutoff energy, energy selected is that of incident e        
                 !      secE = eIn-ejEnSelected
                 !      primaryE = ejEnSelected-eIon
                 !   else 
                 !      secE = ejEnSelected-eIon
                 !      primaryE = eIn-ejEnselected
                 !   end if
                 !end if

                 !particlebasis(partNum)%energy(coll+1) = primaryE 
		 !call init_particle(particlebasis(datasim%secE),secE,particlebasis(partNum)%gen + 1)			
		 !datasim%enPerGen(particlebasis(partNum)%gen + 1) = datasim%enPerGen(particlebasis(partNum)%gen+1) + secE

	         particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - meanexc
	         ! sec E created with energy equal to the mean exc energy minus the ionisation potential, its generation is the one after the current generation
		 call init_particle(particlebasis(datasim%secE),meanexc-eIon,particlebasis(partNum)%gen + 1,particlebasis(partNum)%x(coll),particlebasis(partNum)%y(coll),&
                                    particlebasis(partNum)%z(coll),particlebasis(partNum)%time(coll))		
		 datasim%enPerGen(particlebasis(partNum)%gen + 1) = datasim%enPerGen(particlebasis(partNum)%gen+1) + (meanexc-eIon)
			
		 ! if((partNum.eq.0) .and.(coll.eq.0)) then
		   ! call meanexcenergy(tcs,statebasis,meanexc)
		 ! end if
                 !ejEnSelected = meanexc - eIon
	      else if(ionop .EQ. 'state_distribution') then ! regular excitation energy will be used 
                 !Particle loses energy of pseudostate, subject to indistinguishability
                 !print*, 'enex', enex
                 !Use alternate method to select ejection energy from CCC CS

                 !Get pseudostates at the given incident energy


                 !Symmetrise them

                 !Select positive energy state.




                     
                 if (enlossop .eq. 1) then
                    !Deal with indistinguishibility
                    ejEnSelected = enex - eIon
                    !print*, "IN:", eIn, "SEL:",ejEnSelected
                    if (indist) then
                       if (ejEnSelected .gt. (eIn-eIon)/2.0) then
                          !Above cutoff energy, energy selected is that of incident e       
                          secE = eIn-ejEnSelected-eIon 
                          primaryE = ejEnSelected
                       else 
                          secE = ejEnSelected
                          primaryE = eIn-ejEnselected-eIon
                       end if
                    else
                       secE = enex-eIon
                       primaryE = eIn-enex
                    end if 

	            particlebasis(partNum)%energy(coll + 1) = primaryE
	            !particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - enex
                 else if (enlossop .eq. 2) then
                    !Include recoil for inelastic collisions

                    !Deal with indistinguishibility
                    ejEnSelected = enex - eIon
                    if (indist) then
                       if (ejEnSelected .gt. (eIn-eIon)/2.0) then
                          !Above cutoff energy, energy selected is that of incident e        
                          call inelasticScattering(particlebasis(partNum)%energy(coll),eIn-enex,cosangle,inElEnergyLoss,bmode)
                          secE = eIn-enex
                          primaryE = eIn-inElEnergyLoss 
                       else 
                          call inelasticScattering(particlebasis(partNum)%energy(coll),enex,cosangle,inElEnergyLoss,bmode)
                          secE = enex
                          primaryE = eIn-inElEnergyLoss
                       end if
                    else
                       call inelasticScattering(particlebasis(partNum)%energy(coll),enex,cosangle,inElEnergyLoss,bmode)
                       secE = enex-eIon
                       primaryE = eIn-inElEnergyLoss
                    end if
                    particlebasis(partNum)%energy(coll+1) = primaryE
                    !call inelasticScattering(particlebasis(partNum)%energy(coll),enex,particlebasis(partNum)%costheta(coll),inElEnergyLoss,bmode)
                    !particlebasis(partNum)%energy(coll+1) = particlebasis(partNum)%energy(coll) - inElEnergyLoss
                 end if
	         call init_particle(particlebasis(datasim%secE),secE,particlebasis(partNum)%gen + 1,particlebasis(partNum)%x(coll),particlebasis(partNum)%y(coll),&
                                    particlebasis(partNum)%z(coll), particlebasis(partNum)%time(coll))
                 call recordEjEn(secE,datasim)
	         datasim%enPerGen(particlebasis(partNum)%gen + 1) = datasim%enPerGen(particlebasis(partNum)%gen+1) + secE

                 !ejEnSelected = en
	         ! sec E created with energy equal to positive cont. ionisation energy, its generation is the one after the current generation
	         !call init_particle(particlebasis(datasim%secE),en,particlebasis(partNum)%gen + 1)
                 !call recordEjEn(en,datasim)
	         !datasim%enPerGen(particlebasis(partNum)%gen + 1) = datasim%enPerGen(particlebasis(partNum)%gen+1) + en
	      else if(ionop .EQ. 'SDCS') then !Distribution from SDCS
                 !Instead fit all SDCS in preprocessing ig sdcsOp .eq. 1 
                 !if(data_in%SdcsOp .eq. 0) then
                 !   call get_sdcsatein(sdcsBasis,particlebasis(partNum)%energy(coll),sdcsAtEIn)
                 !else if (data_in%SdcsOp .eq. 1) then 
                 !   call get_sdcsateinfit(sdcsAtEIn,sdcsBasis,particlebasis(partNum)%energy(coll),stateBasis,data_in)
                 !May need to update for more general non-H2 case
                 !end if
 
                 call get_sdcsatein(sdcsBasis,particlebasis(partNum)%energy(coll),sdcsAtEIn)


                 if (debug) then
                    allocate(tcsDebug%CS(sdcsAtEIn%Ncgsdcs))
                    allocate(tcsDebug%en(sdcsAtEIn%Ncgsdcs))
                    tcsDebug%cs(:) = sdcsAtEIn%CS(:)
                    tcsDebug%en(:) = sdcsAtEIn%ecg(:) 
                    call selectstate(tcsDebug,particlebasis(partNum),coll,stateNumDebug)
                    ejEn = tcsDebug%en(stateNumDebug)
                    deallocate(tcsDebug%CS,tcsDebug%en)
                 else 
                    !Choose ejection energy, run monte carlo
                    call selectEjEn(sdcsAtEIn,ejEn) 
                 end if

                 !Deal with indistinguishability
                 ejEnSelected = ejEn	        
                 if (indist)  then
                    if (ejEnSelected .gt. (eIn-eIon)/2.0) then
                       !Above cutoff energy, energy selected is that of incident e        
                       secE = eIn-ejEnSelected-eIon
                       primaryE = ejEnSelected
                    else 
                       secE = ejEnSelected
                       primaryE = eIn-ejEnSelected-eIon
                    end if
                 else
                    secE = ejEn
                    primaryE = eIn-ejEnSelected-eIon
                 end if
                  
                 particlebasis(partNum)%energy(coll+1) = primaryE
	         call init_particle(particlebasis(datasim%secE),secE,particlebasis(partNum)%gen+1,particlebasis(partNum)%x(coll),particlebasis(partNum)%y(coll),&
                                    particlebasis(partNum)%z(coll), particlebasis(partNum)%time(coll))
                 datasim%enPerGen(particlebasis(partNum)%gen + 1) = datasim%enPerGen(particlebasis(partNum)%gen+1) + secE

                 !particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - (ejEn + eIon)				
	         !call init_particle(particlebasis(datasim%secE),ejEn,particlebasis(partNum)%gen + 1)
                 !datasim%enPerGen(particlebasis(partNum)%gen + 1) = datasim%enPerGen(particlebasis(partNum)%gen+1) + ejEn
                 call destructSdcsAtEIn(sdcsAtEIn)
                 call recordEjEn(secE,datasim)
                 !call recordEjEn(ejEn,datasim) !Record ejection energy chosen
              else if (ionop .eq. 'SDCSDist') then
                 !Use analytical distribution from Dalgarno paper
                 call getSdcsDistAtEIn(sdcsBasis,eIn,sdcsAtEIn,tcs)
                 call selectEjEn(sdcsAtEIn, ejEn) 
  
                 !print*, "EJEN: ", ejEn
                 !Deal with indistinguishability
                 !Use same method as symmetrised SDCS for consistency
                 ejEnSelected = ejEn	        
                 if (indist)  then
                    if (ejEnSelected .gt. (eIn-eIon)/2.0) then
                       !Above cutoff energy, energy selected is that of incident e        
                       secE = eIn-ejEnSelected-eIon
                       primaryE = ejEnSelected
                    else 
                       secE = ejEnSelected
                       primaryE = eIn-ejEnSelected-eIon
                    end if
                 else
                    secE = ejEn
                    primaryE = eIn-ejEnSelected-eIon
                 end if
 
                 particlebasis(partNum)%energy(coll+1) = primaryE
	         call init_particle(particlebasis(datasim%secE),secE,particlebasis(partNum)%gen + 1,particlebasis(partNum)%x(coll),particlebasis(partNum)%y(coll), &
                                    particlebasis(partNum)%z(coll), particlebasis(partNum)%time(coll))
                 datasim%enPerGen(particlebasis(partNum)%gen + 1) = datasim%enPerGen(particlebasis(partNum)%gen+1) + secE
 
                 call destructSdcsAtEIn(sdcsAtEIn)
                 call recordEjEn(secE,datasim)
              end if
           end if	
           !particlebasis(partNum)%energy(coll+1) = primaryE
           !call init_particle(particlebasis(datasim%secE),secE,particlebasis(partNum)%gen+1)               
           !datasim%enPerGen(particlebasis(partNum)%gen + 1) = datasim%enPerGen(particlebasis(partNum)%gen+1) + secE
           !call recordEjEn(secE,datasim) !Record ejection energy chosen
        else
	   !write to an error file
        end if	
        !print*, 'postCollEnergy',particlebasis(partNum)%energy(coll + 1) 	
    end subroutine update_energy


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: update_position
    !Purpose: updates the position of the scattered particle and
    !         orientation of the coordinate system.
    !Date last modified: 02/02/2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine update_position(particleIn,rad,costheta,phi,coll,datasim)
        implicit none
        type(simdata)::datasim
        type(particle)::particleIn
        real(dp)::rad,costheta,phi,theta
        real(dp)::xVal,yVal,zVal
        real(dp), dimension(3,3)::rotMat
        integer::coll

        theta = ACOS(costheta)       

        !Construct rotation matrix from theta, phi using analytic formulae
        rotMat(1,1) = costheta*cos(phi)
        rotMat(1,2) = -sin(theta)
        rotMat(1,3) = costheta*sin(phi)
        rotMat(2,1) = sin(theta)*cos(phi)
        rotMat(2,2) = costheta
        rotMat(2,3) = sin(theta)*sin(phi)
        rotMat(3,1) = -sin(phi)
        rotMat(3,2) = 0.0
        rotMat(3,3) = cos(phi)

        !Rotate coordinate system. 
        particleIn%zHat = MATMUL(rotMat,particleIn%zHat)
        particleIn%yHat = MATMUL(rotMat,particleIn%yHat) 
        particleIn%xHat = MATMUL(rotMat,particleIn%xhat)

        !Calculate new position
        particleIn%x(coll+1) = particleIn%x(coll) + rad*particleIn%zHat(1)
        particleIn%y(coll+1) = particleIn%y(coll) + rad*particleIn%zHat(2)
        particleIn%z(coll+1) = particleIn%z(coll) + rad*particleIn%zHat(3)

        xVal = particleIn%x(coll+1)
        yVal = particleIn%y(coll+1)
        zVal = particleIn%z(coll+1)

        !Update incident particle distance from origin
        datasim%inERad = SQRT(xVal**2+yVal**2+zVal**2)

    end subroutine update_position




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: update_H2Data
    !Purpose: updates how many excitations of states of interest
    !         of H2 have occured based on the selected state label
    !         and data.
    !Date last modified: 09/08/2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine update_H2Data(datasim, stateBasis, stateNum)
        use state_class
        implicit none
        type(simdata)::datasim
        type(basis_state)::stateBasis
        integer:: stateNum
 
        !if (stateNum .eq. 6) then    !Old version based on 1-1 stateNum e state correspondence
        if (stateBasis%b(stateNum)%stlabel .eq. 'B1Su') then
           if ((int(stateBasis%b(stateNum)%spin) .eq. 0) .and. (stateBasis%b(stateNum)%ipar .eq. -1) &
              .and. (stateBasis%b(stateNum)%m .eq. 0)) then 
              datasim%B1SuExc = datasim%B1SuExc +1
           else
              print*, "ERROR: invalid state, stopping: update_H2Data"
              stop
           end if
        !else if ((stateNum .eq. 8) .or. (stateNum .eq. 9)) then
        else if (stateBasis%b(stateNum)%stlabel .eq. 'C1Pu') then
           if ((int(stateBasis%b(stateNum)%spin) .eq. 0) .and. (stateBasis%b(stateNum)%ipar .eq. -1) &
              .and. (abs(stateBasis%b(stateNum)%m) .eq. 1)) then 
              datasim%C1PuExc = datasim%C1PuExc +1
           else
              print*, "ERROR: invalid state, stopping: update_H2Data"
              stop
           end if
        else if (statebasis%b(stateNum)%stlabel .eq. 'X1Sg') then
           if (statebasis%b(stateNum)%resolved) then
              if (statebasis%b(stateNum)%v .eq. 1) then
                 datasim%groundv1 = datasim%groundv1 + 1
              else if (statebasis%b(stateNum)%v .eq. 2 ) then
                 datasim%groundv2 = datasim%groundv2 + 1
              end if
           end if
        end if
 
        if (int(statebasis%b(stateNum)%spin) .eq. 0) then
           datasim%singletExc = datasim%singletExc + 1 
        else if (int(statebasis%b(stateNum)%spin) .eq. 1) then
           datasim%tripletExc = datasim%tripletExc + 1
        end if 
 
    end subroutine update_H2Data


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: getIndist
    !Purpose: indentifies whether or not scattered particles are
    !         indistinguishible
    !Date last modified: 20/09/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getIndist(bmode,indist)
        implicit none
        logical, intent(inout)::indist
        logical, intent(in):: bmode

        indist = .false.
        if (.not. bmode) then
           indist = .true.
        end if
    end subroutine

	
    subroutine dissociation(part,statebasis,stateNum,dicsBasis,tcs,datasim,data_input,coll,partNum)
        use input_data
        use state_class
	use totalcs_module	
	implicit none
		
        type(input)::data_input
	type(particle),intent(inout):: part
	type(totalcs),intent(in):: tcs
        type(basis_state),intent(in)::statebasis
        integer,intent(in):: stateNum, coll, partNum
        type(basis_dics),intent(in)::dicsBasis
        type(state)::selstate
        type(simdata), intent(inout):: datasim		
        logical:: diss	! True if dissociation occurs
	real(dp):: df, randNum, en, enIncident, dissEn
                
        dissEn = 0.0d0	
        enIncident = part%energy(coll)
	en=0.0_dp

        !Find energy of selected state
        if (data_in%stateIonop .eq. 1) then
           en = tcs%en(stateNum)  
        else if (data_in%stateIonop .eq. 0) then
           selstate = statebasis%b(stateNum)
           en = selstate%en
        end if
 
        ! Assume dissociation does not occur until proven otherwise
	diss = .FALSE.		
        df = 0.0   

        if (data_input%benchmarkOp .eq. 1) then
           !Dissociation not considered
        else 
  	   ! Get the dissociation fraction for selected state at current energy
           if (data_input%dicsop .eq. 1) then
              df = tcs%df(stateNum)
           else if (data_input%dicsop .eq. 2) then
              if ((en <0.0_dp) .and. (statebasis%b(stateNum)%ion .eqv. .false.) ) then	
                 !excitation collision
	         df = tcs%df(stateNum)
              else if ((en >0.0_dp) .and. (statebasis%b(stateNum)%ion .eqv. .true.)) then
                 !Ionisation collision uses dissociative 
                 !ionisation cross sections.            
                 call getDissIonDf(dicsBasis, enIncident, df)
              end if
           end if		
        end if


        ! Random Number is selected between 0 and 1
	call RANDOM_NUMBER(randNum)
		
        ! Sample to see if dissociation occurs
	if(randNum .lt. df) then
           diss = .TRUE.
           part%diss(coll) = .true.
        else
           part%diss(coll) = .false. 
        end if		
	
        if (diss .eqv. .true.) then
           if (int(statebasis%b(stateNum)%spin) .eq. 1) then
              datasim%tripDiss = datasim%tripDiss + 1 
           else if (int(statebasis%b(stateNum)%spin) .eq. 0) then
              datasim%singDiss = datasim%singDiss + 1
           end if

           if (statebasis%b(stateNum)%stlabel .eq. 'C1Pu') then
              datasim%C1PuDiss = datasim%C1PuDiss + 1
	      !print*, "BOOGALOO"
           else if (statebasis%b(stateNum)%stlabel .eq. 'B1Su') then
              datasim%B1SuDiss = datasim%B1SuDiss + 1
           else if (statebasis%b(stateNum)%stlabel .eq. 'b3Su') then
              datasim%b3SuDiss = datasim%b3SuDiss + 1
           end if  
	end if
  
        !If enex - dissThresh > 0, this is a dissociative pseudostate, record energy release
        !print*, "ENEX, THRESH: ", statebasis%b(stateNum)%enex, statebasis%b(stateNum)%dissThresh
        if ((statebasis%b(stateNum)%resolved) .and. (statebasis%b(stateNum)%enex - statebasis%b(stateNum)%dissThresh .gt. 0.0d0)) then
           dissEn= statebasis%b(stateNum)%enex - statebasis%b(stateNum)%dissThresh
           datasim%dissHeat = datasim%dissHeat + dissEn
           if (partNum .eq. 0) then
              !Update primary energy deposited as heat through dissociation
              datasim%PDissEn = datasim%PDissEn + dissEn
           end if
        else if ((statebasis%b(stateNum)%resolved) .and. (statebasis%b(stateNum)%enex - statebasis%b(stateNum)%dissThresh .lt. 0.0d0))  then
           !Record energy of bound excitation caused by primary
           if (partNum .eq. 0) then
              datasim%PVibEn = datasim%PVibEn + statebasis%b(stateNum)%enex 
           end if
        end if
        !end if
	
        !print*, 'df = ', df
	!print*, 'randNum = ', randNum
	!print*, 'diss = ', diss
        !print*, 'colls = ', part%colls, coll
	part%diss(part%colls) = diss		
	if(diss .eqv. .true.) then ! Update how many dissociations have occured for this generation.
	   datasim%dissPerGen(part%gen) = datasim%dissPerGen(part%gen) + 1
	   datasim%dissociations = datasim%dissociations + 1
	end if
	!print*, 'diss tally ',datasim%dissPerGen(part%gen)
    end subroutine dissociation

	
    subroutine collectdata(datasim,maxgen,datamc,data_input)
        use input_data
        implicit none
	integer,intent(in):: maxgen
        type(input):: data_input
	type(simdata),intent(in):: datasim
	type(simdata),intent(inout):: datamc
	integer:: i

	datamc%secE = datamc%secE + datasim%secE
	datamc%excite = datamc%excite + datasim%excite
	datamc%elastic = datamc%elastic + datasim%elastic
	datamc%gen = datamc%gen + datasim%gen
	datamc%dissociations = datamc%dissociations + datasim%dissociations
        datamc%inERad = datamc%inERad + datasim%inERad
      
        datamc%B1SuExc = datamc%B1SuExc + datasim%B1SuExc	
        datamc%C1PuExc = datamc%C1PuExc + datasim%C1PuExc
        datamc%singletExc = datamc%singletExc + datasim%singletExc
        datamc%tripletExc = datamc%tripletExc + datasim%tripletExc 
        !New parameters
        if (datasim%secE .gt. 0) then
           datamc%W = datamc%W + data_input%energyeV/dble(datasim%secE)
           datamc%singIonPair = datamc%singIonPair + dble(datasim%singletExc)/dble(datasim%secE)
           datamc%tripIonPair = datamc%tripIonPair + dble(datasim%tripletExc)/dble(datasim%secE)
           datamc%B1SuIonPair = datamc%B1SuIonPair + dble(datasim%B1SuExc)/dble(datasim%secE)
           datamc%C1PuIonPair = datamc%C1PuIonPair + dble(datasim%C1PuExc)/dble(datasim%secE)
        end if    
        !Else, do nothing, no ionisation occured, add zero.
        
        datamc%groundv1 = datamc%groundv1 + datasim%groundv1
        datamc%groundv2 = datamc%groundv2 + datasim%groundv2     
        if (datasim%groundv1 .gt. 0) then 
           datamc%groundRatio = datamc%groundRatio + float(datasim%groundv2)/float(datasim%groundv1)
        end if

        datamc%singDiss = datamc%singDiss + datasim%singDiss 
        datamc%tripDiss = datamc%tripDiss + datasim%tripDiss
        datamc%C1PuDiss = datamc%C1PuDiss + datasim%C1PuDiss
        datamc%B1SuDiss = datamc%B1SuDiss + datasim%B1SuDiss
        datamc%b3SuDiss = datamc%b3SuDiss + datasim%b3SuDiss
        datamc%dissHeat = datamc%dissHeat + datasim%dissHeat
        datamc%PDissEn = datamc%PDissEn + datasim%PDissEn
        datamc%PVibEn = datamc%PVibEn + datasim%PVibEn

	do i=1,maxgen
	   datamc%ePerGen(i) = datamc%ePerGen(i) + datasim%ePerGen(i) 
	   datamc%enPerGen(i) = datamc%enPerGen(i) + datasim%enPerGen(i)
        end do	 
	 
	do i=0,maxgen
	   datamc%dissPerGen(i) = datamc%dissPerGen(i) + datasim%dissPerGen(i) 
	   datamc%excPerGen(i) = datamc%excPerGen(i) + datasim%excPerGen(i) 
	   datamc%elPerGen(i) = datamc%elPerGen(i) + datasim%elPerGen(i) 
	   datamc%collPerGen(i) = datamc%collPerGen(i) + datasim%collPerGen(i) 
	end do	 
	 
	if(datasim%PsFormed) then
	   datamc%numPsFormed = datamc%numPsFormed + 1
	end if
	
        !Consistency check data
        do i =1, 1000
           datamc%numEj(i) = datamc%numEj(i) + datasim%numEj(i) 
        end do

        !Update mean and mean square of recorded values, these are used 
        !to calculate statistical uncertainty in the 
        !output of the simulation
        datamc%msecE = datamc%msecE + dble(datasim%secE)/dble(datamc%numSims)
        datamc%mexcite = datamc%mexcite + dble(datasim%excite)/dble(datamc%numSims)
        datamc%mdissociations = datamc%mdissociations + dble(datasim%dissociations)/dble(datamc%numSims)
        datamc%mgen = datamc%mgen + dble(datasim%gen)/dble(datamc%numSims)
        datamc%minERad = datamc%minERad + dble(datasim%inERad)/dble(datamc%numSims)
        datamc%mB1SuExc = datamc%mB1SuExc + dble(datasim%B1SuExc)/dble(datamc%numSims)
        datamc%mC1PuExc = datamc%mC1PuExc + dble(datasim%C1PuExc)/dble(datamc%numSims)
        datamc%msingletExc = datamc%msingletExc + dble(datasim%singletExc)/dble(datamc%numSims)
        datamc%mtripletExc = datamc%mtripletExc + dble(datasim%tripletExc)/dble(datamc%numSims)

        if (datasim%secE .gt. 0) then
           datamc%mW = datamc%mW + (dble(data_input%energyeV)/dble(datasim%secE))/dble(datamc%numSims)
           datamc%mSingIonPair = datamc%mSingIonPair + dble(datasim%singletExc)/dble(datasim%secE*datamc%numSims) 
           datamc%mTripIonPair = datamc%mTripIonPair + dble(datasim%tripletExc)/dble(datasim%secE*datamc%numSims)
           datamc%mB1SuIonPair = datamc%mB1SuIonPair + dble(datasim%B1SuExc)/dble(datasim%secE*datamc%numSims)
           datamc%mC1PuIonPair = datamc%mC1PuIonPair + dble(datasim%C1PuExc)/dble(datasim%secE*datamc%numSims)
        end if

        datamc%mssecE = datamc%mssecE + (dble(datasim%secE)**2)/dble(datamc%numSims)
        datamc%msexcite = datamc%msexcite + (dble(datasim%excite)**2)/dble(datamc%numSims)
        datamc%msdissociations = datamc%msdissociations + (dble(datasim%dissociations)**2)/dble(datamc%numSims)
        datamc%msgen = datamc%msgen + (dble(datasim%gen)**2)/dble(datamc%numSims)
        datamc%msinERad = datamc%msinERad + (dble(datasim%inERad)**2)/dble(datamc%numSims)
        datamc%msB1SuExc = datamc%msB1SuExc + (dble(datasim%B1SuExc)**2)/dble(datamc%numSims)
        datamc%msC1PuExc = datamc%msC1PuExc + (dble(datasim%C1PuExc)**2)/dble(datamc%numSims)
        datamc%mssingletExc = datamc%mssingletExc + (dble(datasim%singletExc)**2)/dble(datamc%numSims)
        datamc%mstripletExc = datamc%mstripletExc + (dble(datasim%tripletExc)**2)/dble(datamc%numSims)

        if (datasim%secE .gt. 0) then
           datamc%msW = datamc%msW + (dble(data_input%energyeV)/dble(datasim%secE))**2/dble(datamc%numSims)
           datamc%msSingIonPair = datamc%msSingIonPair + ((dble(datasim%singletExc)/dble(datasim%secE))**2)/dble(datamc%numSims) 
           datamc%msTripIonPair = datamc%msTripIonPair + ((dble(datasim%tripletExc)/dble(datasim%secE))**2)/dble(datamc%numSims)
           datamc%msB1SuIonPair = datamc%msB1SuIonPair + ((dble(datasim%B1SuExc)/dble(datasim%secE))**2)/dble(datamc%numSims)
           datamc%msC1PuIonPair = datamc%msC1PuIonPair + ((dble(datasim%C1PuExc)/dble(datasim%secE))**2)/dble(datamc%numSims)
        end if
    end subroutine collectdata 

	
    subroutine printsim(simIndex,particlebasis,datasim,bmode)	
        implicit none
	integer,intent(in):: simIndex
	type(particle),dimension(0:1000),intent(in):: particlebasis
	type(simdata),intent(in):: datasim
	logical,intent(in):: bmode ! True if Ps Benchmark Simulation is being done
	integer:: partNum,coll
	
	write(60,*) '-------------------------------- NEW SIMULATION ---------------------------'
	write(60,*) 'simIndex',simIndex
	do partNum=0, datasim%secE   
	   if((.not. bmode) .or. (bmode.and.(partNum.eq.0))) then
	      write(60,*) '--------NEW PARTICLE--------'
	      write(60,*) 'particle',partNum,'gen',particlebasis(partNum)%gen
	      do coll=0,particlebasis(partNum)%colls
	         write(60,*) '-------'
		 write(60,*) 'coll',coll,'state',particlebasis(partNum)%state(coll)
		 write(60,*) 'precollE',particlebasis(partNum)%energy(coll),'postE', particlebasis(partNum)%energy(coll+1)	
		 if(bmode) then
		    write(60,*) 'Ps Formation = ', particlebasis(partNum)%Ps(coll)
		 else
		    write(60,*) 'diss = ', particlebasis(partNum)%diss(coll)
	         end if
	      end do
	   end if
        end do
	write(60,*) 'ions',datasim%secE
	write(60,*) 'elastic', datasim%elastic
	write(60,*) 'excite',datasim%excite		
    end subroutine printsim


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: writePathToFile
    !Purpose: write the path of the incident particle and all secondaries
    !         to a csv file for plotting
    !Date last modified: 21/01/2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine writePathToFile(simIndex,particlebasis,datasim,bmode)
        use input_data
        implicit none
         
        integer,intent(in):: simIndex
        type(particle),dimension(0:1000),intent(in):: particlebasis
        type(simdata),intent(in):: datasim
        logical,intent(in):: bmode ! True if Ps Benchmark Simulation is being done
        character(len=20)::simString, enString,xString,yString,zString,stateNumString
        character(len=20)::partNumString, timeString
        character(len=2)::inString, dissString 
        integer::ii,jj

        !Write particle path to file in order of creation. Collision event at which incident particle first
        !creates a secondary corresponds to event where first secondary in file is created.
        if (bmode) then


        else 
           write(simString,'(I0)') simIndex
           open(70, file="sim"//trim(simString)//".csv")

           write(70,*) datasim%secE+1 !Total number of electrons
           write(70,*) "type,","Incident/Second,","EnSelected,","x(m),","y(m),","z(m),","Diss,","StateNum,","partNum,","Time(s)"               
           do ii = 0, datasim%secE
              if (ii .eq. 0) then
                 inString = "i"
              else
                 inString = "s"
              end if
              inString = trim(inString)//","
  
                
              do jj = 0, particlebasis(ii)%colls
                 if (particlebasis(ii)%diss(jj)) then
                    dissString = "y"
                 else
                    dissString = "n"
                 end if 
                 dissString = trim(dissString)//","

                 write(enString,'(e13.5)') particlebasis(ii)%enArr(jj)
                 enString = trim(enString)//","
                 write(xString,'(e13.5)') particlebasis(ii)%x(jj)
                 xString = trim(xString)//","
                 write(yString,'(e13.5)') particlebasis(ii)%y(jj)
                 yString = trim(yString)//","
                 write(zString,'(e13.5)') particlebasis(ii)%z(jj)
                 zString = trim(zString)//","
                 write(stateNumString, '(i4)') particlebasis(ii)%state(jj)
                 stateNumString = trim(stateNumString)//","
                 write(partNumString, '(i4)') ii
                 partNumString = trim(partNumString)//","
                 write(timeString, '(e13.5)') particlebasis(ii)%time(jj)
                 timeString = trim(timeString)

                 write(70,'(a,a,a,a,a)',advance="no") adjustl("e,"),adjustl(inString),adjustl(enString),adjustl(xString),adjustl(yString)
                 write(70,'(a,a,a,a,a)') adjustl(zString),adjustl(dissString),adjustl(stateNumString),adjustl(partNumString),adjustl(timeString)
              end do
           end do
           close(70)
        end if
    end subroutine writePathToFile 

	
    subroutine simulationresults(datamc,totalSims,en_incident,maxgen,ionop,bmode,runTime,data_input)	
        use input_data
	implicit none

        type(input),intent(in)::data_input
	type(simdata),intent(in):: datamc
	integer,intent(in):: totalSims
	real(dp),intent(in):: en_incident
	integer,intent(in):: maxgen
	character (len=60),intent(in):: ionop ! Ionistation option
	logical,intent(in):: bmode ! True if Ps Benchmark Simulation is being done
	integer:: i
        real(dp)::runTime
	real(dp):: aveIons,aveExcites,aveElastics, avediss ! Average number of each collision type
	real(dp):: aveGens ! Average number of generations created
        real(dp):: aveB1SuExc, aveC1PuExc, aveSingletExc, aveTripletExc  !Average number of different H2 excitations
        real(dp):: aveGroundv1, aveGroundv2, aveGroundRatio
	real(dp),dimension(1000):: aveEperGen ! Average number of electrons per generation
	real(dp),dimension(1000):: aveEnperGen ! Average number of electrons per generation
	real(dp),dimension(0:1000):: avedissPerGen,aveexcPerGen,aveelPerGen,avecollPerGen ! Average number of dissociations per generation
	real(dp):: fSurv ! Survival Fraction (only for Ps Benchmark Simulation)
	real(dp):: aveTRes, eIon
        real(dp):: wInv !1/(mean energy per ion pair)
        real(dp), dimension(0:1000):: aveNumEj
        character (len=3):: enString
        !Total number of dissociations, etc over all secondary generations
        real(dp):: aveIonTotal, aveDissTotal, aveExcTotal, aveRad
        real(dp):: aveElTotal, aveCollTotal, aveNumTotal, aveEnTotal

        !Number of dissociation per ion pair: total, singlet and triplet excitation contributions, triplet b contribution
        real(dp):: aveDITot, aveDISing, aveDITri, aveDIb                 
        real(dp):: aveDIPri, aveDISec  !Average number of ionisations per ion pair, incident and seconday contributions
        real(dp):: aveDissHeat !Average amount of energy given to fragments due to direct dissociative excitation (by all electrons). (excitation energy - dissociation threshold)
        real(dp):: aveB1SuDiss, aveC1PuDiss !Average number of singlet B and C dissociations (both predissociation and radiative decay) 
        real(dp):: avePDissEn    !Average amount of primary energy lost through direct dissociative excitation.
        real(dp):: avePVibEn  !Average energy of primary electron lost through vibrational excitation 

 
        eIon = 15.96632 !Ionisation energy of ground state H2 (eV)	
	
        !Initialises arrays to zero
        aveEperGen(:) = 0
        aveEnperGen(:) = 0
	

        open(unit=70,file='output.txt')
	
	aveIons = float(datamc%secE)/float(totalSims)
	aveExcites = float(datamc%excite)/float(totalSims)
	aveElastics = float(datamc%elastic)/float(totalSims) 
	!aveEperGen = float(datamc%ePerGen(1))/float(totalSims)
	aveGens = float(datamc%gen)/float(totalSims)
	avediss = float(datamc%dissociations)/float(totalSims)

        aveB1SuExc = float(datamc%B1SuExc)/float(totalSims)
        aveC1PuExc = float(datamc%C1PuExc)/float(totalSims)
        aveSingletExc = float(datamc%singletExc)/float(totalSims)
        aveTripletExc = float(datamc%tripletExc)/float(totalSims)
        aveGroundv1= float(datamc%groundv1)/float(totalSims)
        aveGroundv2= float(datamc%groundv2)/float(totalSims)      
        aveGroundRatio= datamc%groundRatio/float(totalSims)
  
        !Calculate mean penetration distance of incident electron
        aveRad = datamc%inERad/float(totalSims)
	
        !Calculate (1/mean energy per ion pair)
        if (en_incident .gt. eIon) then
           wInv = aveIons/en_incident
        else
           wInv = 0.0
        end if

	write(70,*) 'program runtime(s)',runTime	
	write(70,*) 'inital particle energy (eV)', en_incident
	write(70,*) 'total number of simulations', totalSims
	if(bmode .eqv. .false.) then ! ionisation option only applies to default simulation
		write(70,*) 'ionisation treatment: ', ionop
	end if
	! write(70,*) 'totalIons',datamc%secE,'aveIons',aveIons
	! write(70,*) 'totalExcites',datamc%excite,'aveExcites',aveExcites
	! write(70,*) 'totalElastics',datamc%elastic,'aveElastics',aveElastics
	! write(70,*) 'totaldiss',datamc%dissociations,'avediss',avediss
	
	write(70,*) 'ave collisional ionisations',aveIons
	write(70,*) 'ave collisional excitations',aveExcites
	write(70,*) 'ave elastic collisions',aveElastics
	if(.NOT. bmode) then
		write(70,*) 'ave dissociations',avediss
	end if
	
	if(bmode) then ! Ps Formation Benchmark Run
	   fSurv = 100 * (1 - float(datamc%numPsFormed)/float(totalSims))
	   write(70,*) 'fSurv: ', fSurv
	   write(70,*) 'datamc%numPsFormed: ', datamc%numPsFormed
	   aveTRes = datamc%tRes / (float(totalSims) - float(datamc%numPsFormed))
	   write(70,*) 'TRes: ', aveTRes
	else ! Default Simulation Run
           write(70,*) 'vibrationally resolved: ', data_input%vcsOp
	   write(70,*) 'highest number of generations', maxgen
	   write(70,*) 'ave number of generations', aveGens 
           write(70,*) 'reciprocal of mean energy per ion pair:', wInv		
           write(70,*) 'average incident particle final radius:', aveRad
           write(70,*) 'average B1Su excitations per ion pair:', aveB1SuExc/aveIons
           write(70,*) 'average C1Pu excitations per ion pair:', aveC1PuExc/aveIons
           write(70,*) 'average singlet excitations per ion pair:', aveSingletExc/aveIons
           write(70,*) 'average triplet excitations per ion pair:', aveTripletExc/aveIons
           write(70,*) 'average X1Sg(v=1) excitations:', aveGroundv1
           write(70,*) 'average X1Sg(v=2) excitations:', aveGroundv2
           write(70,*) 'average ratio of X1Sg(v=2) to (v=1) excitations', aveGroundRatio	

	   do i=1,maxgen
     	      aveEperGen(i) = float(datamc%ePerGen(i)) / float(totalSims)
	      aveEnperGen(i) = datamc%enPerGen(i) / float(datamc%ePerGen(i))	
	   end do
		
           do i=0,maxgen
	      write(70,*) '------------------------------------------------------------'
	      write(70,*) 'generation:',i
	      if(i .gt. 0) then					
	         write(70,*) 'ave number of electrons:',aveEperGen(i)!,datamc%ePerGen(i)!,&
	         write(70,*) 'ave energy of electrons (eV):',aveEnperGen(i)
	      end if 
	      avedissPerGen(i) = float(datamc%dissPerGen(i)) / float(totalSims)
	      aveexcPerGen(i) = float(datamc%excPerGen(i)) / float(totalSims)
	      aveelPerGen(i) = float(datamc%elPerGen(i)) / float(totalSims)
	      avecollPerGen(i) = float(datamc%collPerGen(i)) / float(totalSims)
	      write(70,*) 'ave number of ionisation collisions:',aveEperGen(i+1)				
	      write(70,*) 'ave number of excitation collisions:', aveexcPerGen(i)				
	      write(70,*) 'ave number of elastic collisions:', aveelPerGen(i)				
	      write(70,*) 'ave number of collisions:',avecollPerGen(i)
	      write(70,*) 'ave number of dissociations:', avedissPerGen(i)
	   end do	 
  
           aveNumTotal = 0.0 
           aveEnTotal = 0.0
           aveIonTotal = 0.0
           aveDissTotal = 0.0
           aveExcTotal = 0.0
           aveElTotal = 0.0
           aveCollTotal = 0.0

           !Print results for secondary electrons overall 
           do i =1,maxgen
              aveNumTotal = aveNumTotal + aveEPerGen(i)
              aveEnTotal = aveEnTotal + aveEnPerGen(i) 
              aveIonTotal = aveIonTotal + aveEPerGen(i+1)
              aveDissTotal = aveDissTotal + float(datamc%dissPerGen(i))
              aveExcTotal = aveExcTotal + float(datamc%excPerGen(i))
              aveElTotal = aveElTotal + float(datamc%elPerGen(i))
              aveCollTotal = aveCollTotal + float(datamc%collPerGen(i))
           end do	
                 

           aveDissTotal = aveDissTotal/float(totalSims)
           aveExcTotal = aveExcTotal/float(totalSims)
           aveElTotal = aveElTotal/float(totalSims) 
           aveCollTotal = aveCollTotal/float(totalSims)
           write(70,*) '------------------------------------------------'
           write(70,*) 'secondary total:'
           write(70,*) 'ave number of electrons:', aveNumTotal 
           write(70,*) 'ave energy of electrons (eV):', aveEnTotal
	   write(70,*) 'ave number of ionisation collisions:',aveIonTotal		
           write(70,*) 'ave number of excitation collisions:', aveExcTotal				
	   write(70,*) 'ave number of elastic collisions:', aveElTotal		
	   write(70,*) 'ave number of collisions:',aveCollTotal
	   write(70,*) 'ave number of dissociations:', aveDissTotal
       
           !Average number of other dissociative effects
           aveDITot = avediss/aveIons 
           aveDISec = aveDissTotal/aveIons  !Definition: average number of secondary dissociations/average of total ionisation (including both secondary and incident)
           aveDIPri = abs(aveDITot - aveDISec) !Definition: same, but average number of incident dissociations in the numerator
           aveDISing = (float(datamc%singDiss)/float(totalSims))/aveIons
           aveDITri = (float(datamc%tripDiss)/float(totalSims))/aveIons
           aveDIb = (float(datamc%b3SuDiss)/float(totalSims))/aveIons
           aveDissHeat = datamc%dissHeat/float(totalSims)
           aveB1SuDiss = (float(datamc%B1SuDiss)/float(totalSims))/aveIons
           aveC1PuDiss = (float(datamc%C1PuDiss)/float(totalSims))/aveIons
           avePDissEn = (datamc%PDissEn/float(totalSims))/en_incident 
           avePVibEn = (datamc%PVibEn/float(totalSims))/en_incident

           write(70,*) '_________________________________________________'
           write(70,*) 'dissociative effects:'
           write(70,*) 'ave dissociations per ion pair (total): ', aveDITot
           write(70,*) 'ave dissociations per ion pair (primary): ', aveDIPri
           write(70,*) 'ave dissociations per ion pair (secondary): ', aveDISec
           write(70,*) 'ave dissociations per ion pair (singlet): ', aveDISing
           write(70,*) 'ave dissociations per ion pair (triplet): ', aveDITri
           write(70,*) 'ave dissociations per ion pair (b3Su): ', aveDIb
           write(70,*) 'ave dissociations per ion pair (B1Su): ', aveB1SuDiss
           write(70,*) 'ave dissociations per ion pair (C1Pu): ', aveC1PuDiss
           write(70,*) 'ave kinetic energy released through dissociation (eV): ', aveDissHeat 
           write(70,*) 'ave fraction of primary energy deposited as heat through dissociation: ', avePDissEn
           write(70,*) 'ave fraction of primary energy lost through vibrational excitation: ', avePVibEn
	end if
	
	
	!!!!!!!!!!!!!!!!!!!!!!!! NICE FORMAT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(.NOT. bmode) then ! Outputs only for Default Sim
 	   write(70,*) '------------------------------------------------'
	   write(70,*) '------------------------------------------------'
	   write(70,*) 'ions: total, gen0,gen1,etc'
	   write(70,*) aveIons
	   do i=0,maxgen-1
	      write(70,*) aveEperGen(i+1)				
           end do
        	
	
	   write(70,*) 'exc: total, gen0,gen1,etc'
	   write(70,*) aveExcites
	   do i=0,maxgen
	      write(70,*) aveexcPerGen(i)				
	   end do	
		
           write(70,*) 'el: total, gen0,gen1,etc'
	   write(70,*) aveElastics
	   do i=0,maxgen
	      write(70,*) aveelPerGen(i)			
	   end do
		
	   write(70,*) 'colls: gen0,gen1,etc'
	   do i=0,maxgen
	      write(70,*) avecollPerGen(i)			
	   end do
		
	   write(70,*) 'diss: total, gen0,gen1,etc'
	   write(70,*) avediss
	   do i=0,maxgen
              write(70,*) avedissPerGen(i)	
	   end do 
		
	   write(70,*) 'enpergen: gen0,gen1,etc'
	   do i=0,maxgen -1
	      write(70,*) aveEnperGen(i+1)
	   end do
	end if
        close(70)		

        !Writes consistency check data to separate file for plotting
        if ((.not. bmode) .and. (data_in%debugOp .eq. 1)) then
           if ((ionop .eq. "SDCS") .or. (ionop .eq. "state_distribution")) then
              aveNumEj(:) = float(datamc%numEj(:))/float(totalSims) 
              write(enString,'(I0)') int(data_input%energyeV)
              if (ionop .eq. 'SDCS') then
                open(70, file = "numEjAtEnSDCS"//trim(enString)//".txt")
              else if (ionop .eq. "state_distribution") then
                 open(70, file = "numEjAtEnDist"//trim(enString)//".txt")
              end if
            
              write(70,*) "Ejection Energy (ev)           Ave Number Ejected"
              do i = 1,1000
                 write(70,*) datamc%ejEnGrid(i), aveNumEj(i)
              end do
              close(70)
           end if
        end if  

    end subroutine simulationresults	
	
    subroutine meanexcendist(statebasis)
        use totalcs_module
	use state_class
	implicit none
	type(basis_state),intent(in):: statebasis
	type(totalcs):: tcs 
	real(dp):: meanexc,en_incident
		
	print*, 'creating secondary electron energy distribution (file=secEenergydist)'
		
	en_incident = 15.96632 !(eV)
		
	open(unit=90,file='secEenergydist')
	!write(90,*) 'en_incident	','meanexc	','meansecE'
	do while(en_incident .LE. 700)
	   call get_csatein(statebasis,en_incident,tcs) ! Create totalcs for the current energy !call print_tcs(tcs) 
	   call meanexcenergy(tcs,statebasis,meanexc)
	   write(90,*) en_incident, meanexc-15.96632
	   en_incident = en_incident + 0.001
	end do
	close(90)
    end subroutine meanexcendist
	
    subroutine elasticScattering(en_incident,costheta, elEnergyLoss,bmode,dcsBasis)
        use input_data 
        use dcs_module 
	implicit none
	real(dp),intent(in):: en_incident	  	  ! initial energy of incident particle (eV)
	real(dp),intent(out):: elEnergyLoss	  ! energy lost by incident particle (eV)
	real(dp):: mParticle,mTarget 			  ! mass of particle and target (kg)
	real(dp),intent(in):: costheta			  ! cosine of scattering angle (RADIANS)
	real(dp):: th,randNum
	real(dp):: PI
	logical,intent(in):: bmode
        type(basis_dcs)::dcsBasis
        type(dcs)::dcsEn
		
	PI = 4.D0*DATAN(1.D0)
		
	if(.NOT. bmode) then ! Run Default Simulation
	   mParticle = 9.10938356E-31		! mass of electron (kg)
	   !mParticle = 9.1093835D-31
	   mTarget = 2 * 1.6726219E-27		! mass of H2 molecule (kg), assuming it is 2 * mass of proton			
	
           if (data_in%momOp .eq. 0) then
              !use inelastic scattering energy loss formula	
	      elEnergyLoss = (2*mParticle*mTarget)/((mParticle+mTarget)*(mParticle+mTarget))*en_incident*(1-costheta)
              !print*, 'theta',theta*PI/180,'factor',(4*mParticle*mTarget)/((mParticle+mTarget)*(mParticle+mTarget))*(sin((theta*PI/180)/2))*(sin((theta*PI/180)/2))
	      !elEnergyLoss = (3.6E-4) * en_incident
           else if (data_in%momOp .eq. 1) then
              !Use momentum transfer cross sections to get average energy loss
              call get_dcsatein(dcsBasis,en_incident,dcsEn,1)    !Elastic scattering is a 1->1 transition
              elEnergyLoss = 2*((mParticle*mTarget)/(mParticle+mTarget)**2)*en_incident*dcsEn%momtcs/dcsEn%intCs
              call delete_dcs(dcsEn) 
           end if
        else
	   ! Random Number is selected between 0 and 1
	   call RANDOM_NUMBER(randNum)
	   th = 180.0 * randNum
	   mParticle = 9.10938356E-31		! mass of positron (kg)
	   !mParticle = 9.11E-31		! mass of positron (kg)
	   !mTarget = 1.6735575E-27		! mass of atomic H (kg)
	   mTarget = 1.6726219E-27		! mass of H atom (kg), assuming it is mass of proton
	   elEnergyLoss = (4*mParticle*mTarget)/((mParticle+mTarget)*(mParticle+mTarget))*en_incident*(sin(th/2))*(sin(th/2))			
	end if
    end subroutine elasticScattering
!
!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: inelasticScattering
    !Purpose: calculates the recoil energy loss in an inelastic collision
    !Date last modified: 29/02/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine inelasticScattering(en_incident, enex, costheta, energyLoss, bmode )
        real(dp),intent(in):: en_incident	  	  ! initial energy of incident particle (eV)
        real(dp),intent(out):: energyLoss	  ! energy lost by incident particle (eV)
	real(dp):: mParticle,mTarget 			  ! mass of particle and target (kg)
	real(dp),intent(in):: costheta			  ! cosine of scattering angle (RADIANS)
        logical::bmode
        real(dp):: enex  !Excitation/internal energy
        real(dp)::en, sqrtVal

        en = en_incident

        !When generalising, masses will be input variables.
        if (bmode) then
           !Positronium benchmark mode
           mTarget = 2 * 1.6726219E-27   ! mass of H2 molecule (kg), assuming it is 2 * mass of proton			
           mParticle = 9.10938356E-31   ! mass of positron (kg)
           energyLoss = en*((2*mTarget*mParticle)/((mTarget+mParticle)**2) + &
                            enex/(en*(1+(mParticle/mTarget))) - &
                            (2*mParticle*costheta/(mParticle+mTarget))*sqrt((mTarget/(mTarget+mParticle))**2 &
                            - enex/(en*(1+mParticle/mTarget))))
        else
           !Default simulation 
           mTarget = 2 * 1.6726219E-27   ! mass of H2 molecule (kg), assuming it is 2 * mass of proton			
           mParticle = 9.10938356E-31   ! mass of electron (kg) 
           sqrtVal =  (mTarget/(mTarget+mParticle))**2 - enex/(en*(1.0+mParticle/mTarget))
           if (sqrtVal .gt. 0.0) then
              energyLoss = en*((2.0*mTarget*mParticle)/((mTarget+mParticle)**2) + &
                            enex/(en*(1.0+(mParticle/mTarget))) - &
                            (2.0*mParticle*costheta/(mParticle+mTarget))*sqrt((mTarget/(mTarget+mParticle))**2 &
                            - enex/(en*(1.0+mParticle/mTarget)))) 
           else
              print*, "SQRTVAL: ", sqrtVal
              print*, "ENEX, EIN: ", enex, en
              print*, "ERROR: ratio excitationEn/eIncident .gt. 1/(1+mParticle/mTarget), stopping"
              stop
           end if
        end if
    end subroutine inelasticScattering

    subroutine selectstate(tcs,part,coll,stateNum)
        use input_data
        use totalcs_module 
           
        implicit none 
        integer, intent(in):: coll
        type(totalcs), intent(in):: tcs
        type(particle), intent(inout):: part ! The current particle
        integer,intent(out):: stateNum
        real(dp), dimension(:),allocatable:: stateprob
        real(dp):: randNum
        integer:: i

        allocate(stateprob(tcs%Nmax))

        ! stateprob is cumulative probabilty for each state
        call normalisetcs(tcs,stateprob,tcs%Nmax)
           
        ! Random Number is selected between 0 and 1
        call RANDOM_NUMBER(randNum)
           
        ! Determines which state (the 'i'th process in the array) is selected
        i = 1
        do while (stateprob(i) .le. randNum)
           i = i + 1		
        end do    

        ! Ensures that the state is valid
        ! Due to numerical accuracy,stateprob(tcs%Nmax) may be less than 1.0
        if(i .lt. tcs%Nmax) then
           stateNum = i
        else
           stateNum = tcs%Nmax
        end if

        part%state(coll) = stateNum ! records the state chosen at this collision for this particle
        if (tcs%en(i) .gt. 0.0) then
           !print*, tcs%en(i) 
        end if	


        !print*, 'randNum',randNum
        ! ! Be careful printing this next line if stateNum=1
        if(stateNum .eq. 1) then
        !print*,'stateprob(stateNum)',stateprob(stateNum)
        else
        !print*,'stateprob(stateNum)',stateprob(stateNum),'stateprob(stateNum-1)',stateprob(stateNum-1)
        end if
        !print*, 'stateNum =',stateNum

        !Write function to sample path length at given energy
        deallocate(stateprob) !deallocate variable allocated in 'normalisetcs' 	
    end subroutine selectstate


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: selectPath
    !Purpose: selects path length from TCS at given energy.
    !Date last modified: 29/01/2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine selectPath(path,stateBasis,eIncident)
        use state_class
        use totalcs_module
        use input_data
        implicit none
        real(dp)::eIncident, randVal, sigma
        real(dp),intent(out):: path
        type(basis_state)::stateBasis
        type(totalcs)::tcs
 
        randVal = 0.0
 
        !Get total cross section at given energy.
        call get_csatein(stateBasis,eIncident,tcs)
        sigma = SUM(tcs%cs)*(data_in%bohrRadius**2)  !Convert to SI
 
        !Draw path length from exponential distribution
        if(data_in%radop .eq. 1) then
           call RANDOM_NUMBER(randVal)
           if (.not.(randVal .gt. 0.0)) then
              do while (.not.(randVal .gt. 0.0))
                 call RANDOM_NUMBER(randVal)
              end do
           end if
           path = (-1.0)*(1/(data_in%density*sigma))*log(randVal)
        else if (data_in%radop .eq. 0) then
           path = 1/(data_in%density*sigma) 
        end if
        
        call delete_totalcs(tcs)
  
    end subroutine selectPath

  
    SUBROUTINE timeElapsed(energy,crossSection,datasim)
        use Ps_module	
	IMPLICIT NONE
	REAL(dp), INTENT(IN) :: energy
	REAL(dp) :: m_p, eVToJoules, constant, con2, m
	real(dp), intent(in) :: crossSection
	!REAL*8, INTENT(OUT) :: duration
	type(simdata),intent(inout):: datasim
	
	m_p = 	9.109E-31 
	m = 9.1093897
	m_p = 9.10938356E-31 ! mass of positron
	eVToJoules = 1.60218E-19
		
	! IF(stateNum .EQ. 4) THEN
		! crossSection = VarPs%cs_Ps
	! ELSE IF(stateNum .EQ. 1) THEN
		! crossSection = VarPs%cs_el
	! ELSE IF(stateNum .EQ. 3) THEN
		! crossSection = VarPs%cs_ion
	! ELSE IF(stateNum .EQ. 2) THEN
		! crossSection = VarPs%cs_exc
	! ELSE 
		! !crossSection = 0.0
	! END IF
	
	! Convert Cross section from A^2 to m^2
	!crossSection = crossSection * 10**(-20)
	
	! constant is equal to eVToJoules/m_p
	constant = 1.7565E11
	
	con2 = 1!*10**(14)
	
	!velocity = (2 * energy * eVToJoules / m_p) ** (0.5)
	
	!duration = 1 / (crossSection * velocity)	
	
	!duration = ((m_p)**(1/2)) / (crossSection * ((2 * energy * eVToJoules)**(1/2)))
	!duration = ((m_p)**(1/2)) / (((2 * energy * eVToJoules)**(1/2)))
	
	!duration = (m**0.5)/((2*energy * 1.60218)**(0.5) * crossSection) / 10 ! in terms of s*m^-3 * 10^15
	!duration = ((m/(2*energy*1.60218))**(0.5))/(crossSection*10*3.2) ! in terms of s*m^-3 * 10^15
	
	! THis is what I had before the update
	!velocity = (2 * energy * constant) ** (0.5)
	!duration = ((m/(2*energy*1.602176565))**(0.5))/(crossSection*10) ! in terms of s*m^-3 * 10^15
	
	
	
	! Update 03/07/2019 !
	
	datasim%duration = 10E15 * sqrt(m_p / (2*energy*eVToJoules)) * 10E20 / crossSection ! puts it into 10^15 form?
	datasim%duration = sqrt(m_p / (2*energy*eVToJoules)) * 10E20 / crossSection ! unaltered
	datasim%duration = 10E-15 * sqrt(m_p / (2*energy*eVToJoules)) * 10E20 / crossSection ! altered?
	!print*, 'duration', datasim%duration
	
	!write(2,*) 'vel', velocity, 'cs', crossSection, 'energy', energy, 'duration', duration	
    END SUBROUTINE timeElapsed


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: selectAngle
    !Purpose: selects a scattering angle (cos(theta)) from cross sections
    !         differential in solid angle (dcs)
    !Date last modified: 06/07/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE selectAngle(eldcs,eIncident,costheta)
        use dcs_module
        implicit none
        type(basis_dcs), intent(in)::eldcs
        real(dp), intent(in)::eIncident 
        real(dp), intent(out)::costheta 
        type(dcs)::dcsBelow, dcsAbove
        integer::ii, numPoints
        real(dp)::total,  randNum
        real(dp), allocatable, dimension(:)::sums !array of partial sums 
        real(dp), allocatable, dimension(:)::dcsArray !array of dcs at incident energy
        real(dp)::intp1, intp2, intp3 !Used to break up long interpolation formula
   
        !Finds dcs for above and below given incident energy
        ii =1
        do while( (ii .le. eldcs%Ndcs) .and. (eldcs%ein(ii) .lt. eIncident) )
           ii = ii +1
        end do
        call copyDcs(dcsBelow,eldcs%dcsbasis(ii-1))
        call copyDcs(dcsAbove,eldcs%dcsbasis(ii))
        !!dcsBelow = eldcs%dcsbasis(ii-1)    
        !!dcsAbove = eldcs%dcsbasis(ii)
   
        allocate(sums(dcsBelow%Nth), dcsArray(dcsBelow%Nth))
      
        do ii = 1, dcsBelow%Nth
           !Interpolates to find cross sections at given energy.
           intp1 = (eIncident-dcsBelow%en_incident)
           intp2 = (dcsAbove%dcs(ii)-dcsBelow%dcs(ii))
           intp3 = (dcsAbove%en_incident-dcsBelow%en_incident)
           dcsArray(ii) = dcsBelow%dcs(ii) + intp1*(intp2/intp3)
        end do
   
        !Computes total cross section
        call getSum(dcsArray, dcsBelow%Nth, total) 
   
        !Constructs distribution.      
        sums = 0.d0 !Initialises array
        numPoints = dcsBelow%Nth !Number of cross sections, also num angle points
        sums(1) = (dcsArray(1)/total)
        do ii = 2, numPoints
           sums(ii) = sums(ii-1) + (dcsArray(ii)/total)
           !print*, "SUM, ANGLE: ", sums(ii), dcsBelow%theta(ii)
        end do
        !sums(:) = sums(:)*abs(sin(dcsBelow%theta(:)))
   
        !Generates random number
        call RANDOM_NUMBER(randNum)
        ii = 1          
        do while (sums(ii) .lt. randNum)
           ii = ii + 1
        end do
   
        !Selected angle
        costheta = cos(dcsBelow%theta(ii)) 
   
        call delete_dcs(dcsBelow)
        call delete_dcs(dcsAbove)
        deallocate(sums, dcsArray)
   
    END SUBROUTINE selectAngle



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: getSum
    !Purpose: sums over differential cross sections in imported array
    !Date last modified: 27/03/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getSum(dcsArrIn, Ndcs, total)
        implicit none
        real(dp), dimension(:), intent(in)::dcsArrIn
        real(dp)::total
        integer::ii, Ndcs
   
        !Sums differential cross sections
        total = 0
        do ii = 1, Ndcs
           total = total + dcsArrIn(ii)
        end do   
    end subroutine getSum



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: selectEjEn
    !Purpose: selects an ejection energy for secondary
    !         electrons from sdcs, defined over 2D
    !         fine energy grid.
    !Date last modified: 06/05/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine selectEjEn(sdcsAtEIn, ejEn)
        use sdcs_module
        implicit none
        type(sdcs)::sdcsAtEIn
        real(dp)::ejEn !Selected ejection energy
        real(dp), allocatable, dimension(:)::enProb
        real(dp)::randNum !Uniform random number
        integer::enNum, ii
   
        allocate(enProb(sdcsAtEIn%Ncgsdcs))    
   
        !Construct cumulative distribution
        call normalisesdcs(sdcsAtEIn, enProb, sdcsAtEIn%Ncgsdcs)
   
        !Select energy
        call RANDOM_NUMBER(randNum)
   
        !Find selected state
        ii = 1
        do while (enProb(ii) .lt. randNum)
           ii = ii +1
        end do
   
        !Account for rare case when eIncident -eIon < deltaE, where deltaE is the
        !grid step size. 
        if (int(enProb(1)) .eq. 2) then
           !Normalisesdcs sets enProb(1) to 2 when this occurs       
           ejEn = sdcsAtEIn%ecg(2) - sdcsAtEIn%ecg(1) !Ensures no further events can occur. 
        else
           !Check validity
           if (ii .lt. sdcsAtEIn%Ncgsdcs) then
              enNum = ii
           else
              enNum = sdcsAtEIn%Ncgsdcs
           end if
           ejEn = sdcsAtEIn%ecg(enNum)
        end if 
        deallocate(enProb) 
    end subroutine selectEjEn


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: sampleNormal
    !Purpose: samples normal probability distribution with given mean 
    !         mu and standard deviation sigma
    !Date last modified: 14/06/2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sampleNormal(xToGet, mu, sigma)
        implicit none
        real(dp), intent(inout):: xToGet
        real(dp):: mu, sigma   !Mean and standard deviation
        real(dp), dimension(:), allocatable:: xVec, PVec, func
        real(dp):: dx, randVal
        integer:: ii, numX
        integer:: istart, istop 
        logical:: found
  
        dx = 0.1
        !x domain will be from -5*sigma to 5*sigma
        numX = int(10.0d0*sigma/dx) +1
        allocate(xVec(numX))
        xVec(1) = -5.0d0*sigma
        do ii =2, numX
           xVec(ii) = xVec(ii-1) + dx
        end do
 
        allocate(PVec(numX), func(numX))

        PVec(:) = 0.0d0
        func(:) = 0.0d0
        !Use composite integration rule to form cumulative probability distribution
        func(:) = (1.0d0/(sqrt(2.0d0*4.0d0*atan(1.0d0))*sigma))*exp(-(xVec(:)-mu)**2/(2.0d0*sigma**2))
        istart = 1 
        istop = numX+1
        PVec(1) = 0.0d0    !P(-infinity) = 0, cumulative distribution over the reals
        do ii=istart+2, istop-2,2
           !even points: Simpson 3/8 rule:
           PVec(ii+1)=PVec(ii-2)+(func(ii-2)+3.0d0*func(ii-1)+3.0d0*func(ii)+func(ii+1))*0.375d0*dx
           !odd points: Simpson 1/3 rule:
           PVec(ii)=PVec(ii-2)+(func(ii-2)+4.0d0*func(ii-1)+func(ii))*dx/3.0d0
        enddo
 
        call random_number(randVal)
 
        found = .false.
        ii = 1
        do while( .not. found) 
           if (PVec(ii) .ge. randVal) then
              found = .true.
              xToGet = xVec(ii)
           end if
           ii = ii+1
        end do

        if (.not. found) then
           print*, "ERROR: value not selected from normal distribution. Check function 'sampleNormal' in monteCarlo.f90"
           stop
        end if
        deallocate(xVec)
    end subroutine sampleNormal



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: calculateErrors
    !Purpose: calculates the errors and other useful values and prints
    !         them to a file. 
    !Date last modified: 15/06/2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calculateErrors(datamcArray, numVars, data_input)
        use input_data
        implicit none
        integer:: numVars   !Number of variations of the input parameters
        type(input)::data_input
        type(simdata), dimension(numVars):: datamcArray
        integer:: numSimsVal, ii
        !Mean values of the outputs of each run of the simulation
        real(dp):: mObIons, mObExcites, mObElastics, mObGens, mObDiss, mObB1SuExc, mObC1PuExc
        real(dp):: mObSingletExc, mObTripletExc, mObRad
        real(dp):: mObSingletPerIonPair, mObTripletPerIonPair 
        real(dp):: mObB1SuIonPair, mObC1PuIonPair
        real(dp):: mObW
        !Mean squares of the outputs of each run of the simulation
        real(dp):: msObIons, msObExcites, msObElastics, msObGens, msObDiss, msObB1SuExc, msObC1PuExc
        real(dp):: msObSingletExc, msObTripletExc, msObRad
        real(dp):: msObSingletPerIonPair, msObTripletPerIonPair
        real(dp):: msObB1SuIonPair, msObC1PuIonPair
        real(dp):: msObW
        !Observed ncertainties(standard deviation) in the outputs of each run of the simulation
        real(dp):: uIons, uExcites, uElastics, uGens, uDiss, uB1SuExc, uC1PuExc
        real(dp):: uSingletExc, uTripletExc, uRad
        real(dp):: uB1SuPerIonPair, uC1PuPerIonPair
        real(dp):: uSingletPerIonPair, uTripletPerIonPair
        real(dp):: uB1SuIonPair, uC1PuIonPair
        real(dp):: uW
        !Final unceratinties after subtractin statistical part
        real(dp):: uIonsFin, uExcitesFin, uGensFin, uDissFin, uB1SuExcFin, uC1PuExcFin
        real(dp):: uSingletExcFin, uTripletExcFin, uRadFin
        real(dp):: uB1SuPerIonPairFin, uC1PuPerIonPairFin
        real(dp):: uSingletPerIonPairFin, uTripletPerIonPairFin
        real(dp):: uB1SuIonPairFin, uC1PuIonPairFin
        real(dp):: uWFin
        !Statistical uncertainties 
        real(dp):: uStatIons, uStatExcites, uStatGens, uStatDiss
        real(dp):: uStatB1SuExc, uStatC1PuExc, uStatSingletExc, uStatTripletExc
        real(dp):: uStatRad, uStatSingletPerIonPair, uStatTripletPerIonPair
        real(dp):: uStatB1SuIonPair, uStatC1PuIonPair
        real(dp):: uStatW
 
        !Find the mean and standard deviation of each of the calculated
        !quantities in each of the runs
        mObIons = 0.0d0
        mObW = 0.0d0
        mObExcites = 0.0d0
        mObElastics = 0.0d0
        mObGens = 0.0d0
        mObDiss = 0.0d0
        mObB1SuExc = 0.0d0
        mObC1PuExc = 0.0d0
        mObSingletExc = 0.0d0
        mObTripletExc = 0.0d0
        mObRad = 0.0d0
        mObB1SuIonPair = 0.0d0
        mObC1PuIonPair = 0.0d0

        msObIons = 0.0d0
        msObW = 0.0d0
        msObExcites = 0.0d0
        msObElastics = 0.0d0
        msObGens = 0.0d0
        msObDiss = 0.0d0
        msObB1SuExc = 0.0d0
        msObC1PuExc = 0.0d0
        msObSingletExc = 0.0d0
        msObTripletExc = 0.0d0
        msObRad = 0.0d0
        msObSingletPerIonPair =0.0_dp
        msObTripletPerIonPair =0.0_dp
        msObB1SuIonPair = 0.0d0
        msObC1PuIonPair = 0.0d0

        mObIons = 0.0_dp
        mObW = 0.0_dp
        mObExcites = 0.0_dp
        mObElastics = 0.0_dp
        mObGens = 0.0_dp
        mObDiss = 0.0_dp
        mObB1SuExc = 0.0_dp
        mObC1PuExc = 0.0_dp
        mObSingletExc = 0.0_dp
        mObTripletExc = 0.0_dp
        mObRad = 0.0_dp
        mObSingletPerIonPair =0.0_dp
        mObTripletPerIonPair =0.0_dp
        mObB1SuIonPair = 0.0_dp
        mObC1PuIonPair = 0.0_dp

        uStatIons = 0.0d0
        uStatW = 0.0d0
        uStatExcites = 0.0d0
        !uStatElastics = uStat
        uStatGens = 0.0d0
        uStatDiss = 0.0d0
        uStatB1SuExc = 0.0d0
        uStatC1PuExc = 0.0d0
        uStatSingletExc = 0.0d0
        uStatTripletExc = 0.0d0
        uStatRad = 0.0d0
        uStatSingletPerIonPair = 0.0d0
        uStatTripletPerIonPair = 0.0d0
        uStatB1SuIonPair = 0.0d0
        uStatC1PuIonPair = 0.0d0

        numSimsVal = 0
        do ii= 1, data_in%totalVariations
           numSimsVal = datamcArray(ii)%numSims
           !print*, mIons, datamcArray(ii)%secE, dble(datamcArray(ii)%secE)/dble(numSimsVal)           
           !print*, "SIMS: ", dble(numSimsVal) 
           !print*, "W: ", datamcArray(ii)%W
           !print*, ii

           !Sum statistical uncertainties 
           uStatIons = uStatIons + (datamcArray(ii)%mssecE - datamcArray(ii)%msecE**2)/dble(numSimsVal) 
           uStatW = uStatW + (datamcArray(ii)%msW - datamcArraY(ii)%mW**2)/dble(numSimsVal)
           uStatExcites = uStatExcites + (datamcArray(ii)%msexcite - datamcArray(ii)%mexcite**2)/dble(numSimsVal)
           !uStatElastics = uStatElastics + (datamcArray(ii)%ms - datamcArray(ii)%m  ) 
           uStatGens = uStatGens + (datamcArray(ii)%msgen - datamcArray(ii)%mgen**2)/dble(numSimsVal) 
           uStatDiss = uStatDiss + (datamcArray(ii)%msdissociations - datamcArray(ii)%mdissociations**2)/dble(numSimsVal) 
           uStatB1SuExc = uStatB1SuExc + (datamcArray(ii)%msB1SuExc - datamcArray(ii)%mB1SuExc**2)/dble(numSimsVal) 
           uStatC1PuExc = uStatC1PuExc + (datamcArray(ii)%msC1PuExc - datamcArray(ii)%mC1PuExc**2)/dble(numSimsVal) 
           uStatSingletExc = uStatSingletExc + (datamcArray(ii)%mssingletExc - datamcArray(ii)%msingletExc**2)/dble(numSimsVal) 
           uStatTripletExc = uStatTripletExc + (datamcArray(ii)%mstripletExc - datamcArray(ii)%mtripletExc**2)/dble(numSimsVal) 
           uStatRad = uStatRad + (datamcArray(ii)%msinERad - datamcArray(ii)%minERad**2)/dble(numSimsVal) 
           uStatSingletPerIonPair = uStatSingletPerIonPair + (datamcArray(ii)%msSingIonPair - datamcArray(ii)%mSingIonPair**2)/dble(numSimsVal)
           uStatTripletPerIonPair = uStatTripletPerIonPair + (datamcArray(ii)%msTripIonPair - datamcArray(ii)%mTripIonPair**2)/dble(numSimsVal) 
           uStatB1SuIonPair = uStatB1SuIonPair + (datamcArray(ii)%msB1SuIonPair - datamcArray(ii)%mB1SuIonPair**2)/dble(numSimsVal)
           uStatC1PuIonPair = uStatC1PuIonPair + (datamcArray(ii)%msC1PuIonPair - datamcArray(ii)%mC1PuIonPair**2)/dble(numSimsVal)

           !Calculate observed mean square and mean
           mObIons = mObIons + dble(datamcArray(ii)%secE)/dble(numSimsVal)
           mObW = mObW + dble(datamcArray(ii)%W)/dble(numSimsVal) 
           mObExcites = mObExcites + dble(datamcArray(ii)%excite)/dble(numSimsVal)
           mObElastics = mObElastics + dble(datamcArray(ii)%elastic)/dble(numSimsVal)
           mObGens = mObGens + dble(datamcArray(ii)%gen)/dble(numSimsVal)
           mObDiss = mObDiss + dble(datamcArray(ii)%dissociations)/dble(numSimsVal)
           mObB1SuExc = mObB1SuExc + dble(datamcArray(ii)%B1SuExc)/dble(numSimsVal)
           mObC1PuExc = mObC1PuExc + dble(datamcArray(ii)%C1PuExc)/dble(numSimsVal)
           mObSingletExc = mObSingletExc + dble(datamcArray(ii)%singletExc)/dble(numSimsVal)
           mObTripletExc = mObTripletExc + dble(datamcArray(ii)%tripletExc)/dble(numSimsVal)
           mObRad = mObRad + dble(datamcArray(ii)%inERad)/dble(numSimsVal)
           mObSingletPerIonPair = mObSingletPerIonPair + dble(datamcArray(ii)%singIonPair)/dble(numSimsVal) 
           mObTripletPerIonPair = mObTripletPerIonPair + dble(datamcArray(ii)%tripIonPair)/dble(numSimsVal)
           mObB1SuIonPair = mObB1SuIonPair + dble(datamcArray(ii)%B1SuIonPair)/dble(numSimsVal)
           mObC1PuIonPair = mObC1PuIonPair + dble(datamcArray(ii)%C1PuIonPair)/dble(numSimsVal)
 
           msObIons = msObIons + (dble(datamcArray(ii)%secE)/dble(numSimsVal))**2
           msObW = msObW + (dble(datamcArray(ii)%W)/dble(numSimsVal))**2 
           msObExcites = msObExcites + (dble((datamcArray(ii)%excite))/dble(numSimsVal))**2
           msObElastics =  msObElastics + (dble((datamcArray(ii)%elastic))/dble(numSimsVal))**2
           msObGens = msObGens + (dble((datamcArray(ii)%gen))/dble(numSimsVal))**2
           msObDiss = msObDiss + (dble((datamcArray(ii)%dissociations))/dble(numSimsVal))**2
           msObB1SuExc = msObB1SuExc + (dble((datamcArray(ii)%B1SuExc))/dble(numSimsVal))**2
           msObC1PuExc = msObC1PuExc + (dble((datamcArray(ii)%C1PuExc))/dble(numSimsVal))**2
           msObSingletExc = msObSingletExc + (dble((datamcArray(ii)%singletExc))/dble(numSimsVal))**2
           msObTripletExc = msObTripletExc + (dble((datamcArray(ii)%tripletExc))/dble(numSimsVal))**2
           msObRad = msObRad + (dble((datamcArray(ii)%inERad))/dble(numSimsVal))**2
           msObSingletPerIonPair = msObSingletPerIonPair + (dble(datamcArray(ii)%singIonPair)/dble(numSimsVal))**2
           msObTripletPerIonPair = msObTripletPerIonPair + (dble(datamcArray(ii)%tripIonPair)/dble(numSimsVal))**2
           msObB1SuIonPair = msObB1SuIonPair + (dble(datamcArray(ii)%B1SuIonPair)/dble(numSimsVal))**2
           msObC1PuIonPair = msObC1PuIonPair + (dble(datamcArray(ii)%C1PuIonPair)/dble(numSimsVal))**2 
        end do

        print*, "W recorded" 
        print*, msObW
        print*, mObW


        !Better to use formulas from sokolnikoff
        !Uncertainties in outputs due to varying input values:
        !- Var(X) approx = (1/(N-1))*(sum(X_j^2) - (1/N)*(sum(X_j)^2))
        !- E(X) approx = (1/N)*sum(X_j)
      
        !Same, but from paper:
        ! - sigma_{obs}^2 = (1/(N-1))*sum((X_j-XBar)^2)
        !These are actually the same formula.
        !For large N this is essentially: mean(X^2) - (mean(X))^2

        uStatIons = uStatIons/dble(data_in%totalVariations)
        uStatW = uStatW/dble(data_in%totalVariations)
        uStatExcites = uStatExcites/dble(data_in%totalVariations)
        uStatGens = uStatGens/dble(data_in%totalVariations)
        uStatDiss = uStatDiss/dble(data_in%totalVariations)
        uStatB1SuExc = uStatB1SuExc/dble(data_in%totalVariations)
        uStatC1PuExc = uStatC1PuExc/dble(data_in%totalVariations)
        uStatSingletExc = uStatSingletExc/dble(data_in%totalVariations)
        uStatTripletExc = uStatTripletExc/dble(data_in%totalVariations)
        uStatRad = uStatRad/dble(data_in%totalVariations)
        uStatSingletPerIonPair = uStatSingletPerIonPair/dble(data_in%totalVariations)
        uStatTripletPerIonPair = uStatTripletPerIonPair/dble(data_in%totalVariations)
        uStatB1SuIonPair = uStatB1SuIonPair/dble(data_in%totalVariations)
        uStatC1PuIonPair = uStatC1PuIonPair/dble(data_in%totalVariations)   

        print*, "UStat"
        print*, sqrt(uStatIons)
        print*, sqrt(uStatW)
        print*, sqrt(uStatExcites)
        print*, sqrt(uStatGens)
        print*, sqrt(uStatDiss)
        print*, sqrt(uStatB1SuExc)
        print*, sqrt(uStatC1PuExc)
        print*, sqrt(uStatSingletExc)
        print*, sqrt(uStatTripletExc)
        print*, sqrt(uStatRad)
        print*, sqrt(uStatSingletPerIonPair)
        print*, sqrt(uStatTripletPerIonPair)
        print*, sqrt(uStatB1SuIonPair)
        print*, sqrt(uStatC1PuIonPair)

        mObIons = mObIons/dble(data_input%totalVariations)
        mObW = mObW/dble(data_input%totalVariations)
        mObExcites = mObExcites/dble(data_input%totalVariations)
        mObElastics = mObElastics/dble(data_input%totalVariations) 
        mObGens = mObGens/dble(data_input%totalVariations) 
        mObDiss = mObDiss/dble(data_input%totalVariations)
        mObB1SuExc = mObB1SuExc/dble(data_input%totalVariations) 
        mObC1PuExc = mObC1PuExc/dble(data_input%totalVariations)
        mObSingletExc = mObSingletExc/dble(data_input%totalVariations)
        mObTripletExc = mObTripletExc/dble(data_input%totalVariations) 
        mObRad = mObRad/dble(data_input%totalVariations) 
        mObSingletPerIonPair = mObSingletPerIonPair/dble(data_input%totalVariations)                     
        mObTripletPerIonPair = mObTripletPerIonPair/dble(data_input%totalVariations)                     
        mObB1SuIonPair = mObB1SuIonPair/dble(data_input%totalVariations)
        mObC1PuIonPair = mObC1PuIonPair/dble(data_input%totalVariations)
 
        msObIons = msObIons/dble(data_input%totalVariations) 
        msObW = msObW/dble(data_input%totalVariations)
        msObExcites = msObExcites/dble(data_input%totalVariations) 
        msObElastics = msObElastics/dble(data_input%totalVariations) 
        msObGens = msObGens/dble(data_input%totalVariations)
        msObDiss = msObDiss/dble(data_input%totalVariations)
        msObB1SuExc = msObB1SuExc/dble(data_input%totalVariations) 
        msObC1PuExc = msObC1PuExc/dble(data_input%totalVariations) 
        msObSingletExc = msObSingletExc/dble(data_input%totalVariations) 
        msObTripletExc = msObTripletExc/dble(data_input%totalVariations) 
        msObRad = msObRad/dble(data_input%totalVariations) 
        msObSingletPerIonPair = msObSingletPerIonPair/dble(data_input%totalVariations)
        msObTripletPerIonPair = msObTripletPerIonPair/dble(data_input%totalVariations)
        msObB1SuIonPair = msObB1SuIonPair/dble(data_input%totalVariations)
        msObC1PuIonPair = msObC1PuIonPair/dble(data_input%totalVariations)


        uIons = sqrt((msObIons-mObIons**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))   
        uW = sqrt((msObW - mObW**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uExcites = sqrt((msObExcites-mObExcites**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uElastics = sqrt((msObElastics-mObElastics**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uGens = sqrt((msObGens-mObGens**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uDiss = sqrt((msObDiss-mObDiss**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uB1SuExc = sqrt((msObB1SuExc-mObB1SuExc**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uC1PuExc = sqrt((msObC1PuExc-mObC1PuExc**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uSingletExc = sqrt((msObSingletExc-mObSingletExc**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uTripletExc = sqrt((msObTripletExc-mObTripletExc**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uRad = sqrt((msObRad-mObRad**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uSingletPerIonPair = sqrt((msObSingletPerIonPair - mObSingletPerIonPair**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uTripletPerIonPair = sqrt((msObTripletPerIonPair - mObTripletPerIonPair**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))
        uB1SuIonPair = sqrt((msObB1SuIonPair-mObB1SuIonPair**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1))) 
        uC1PuIonPair = sqrt((msObC1PuIonPair-mObC1PuIonPair**2)*(dble(data_input%totalVariations)/dble(data_input%totalVariations-1)))

        print*, "U without stat"
        print*, uIons 
        print*, uW
        print*, uExcites 
        print*, uGens 
        print*, uDiss 
        print*, uB1SuExc 
        print*, uC1PuExc 
        print*, uSingletExc 
        print*, uTripletExc 
        print*, uRad 
        print*, uSingletPerIonPair
        print*, uTripletPerIonPair
        print*, uB1SuIonPair
        print*, uC1PuIonPair

        !Subtract statistical uncertainty to get final result. This is
        !sigma(input) = sqrt( sigma(observed)^2 - sigma(statistical)^2)
        uIonsFin = sqrt(uIons**2 - uStatIons)   !uStatIons actually stores mean variance
        uWFin = sqrt(uW**2 - uStatW)
        uExcitesFin = sqrt(uExcites**2 - uStatExcites)
        uGensFin = sqrt(uGens**2 - uStatGens)
        uDissFin = sqrt(uDiss**2 - uStatDiss)
        uB1SuExcFin = sqrt(uB1SuExc**2 - uStatB1SuExc)
        uC1PuExcFin = sqrt(uC1PuExc**2 - uStatC1PuExc)
        uSingletExcFin = sqrt(uSingletExc**2 - uStatSingletExc)
        uTripletExcFin = sqrt(uTripletExc**2 - uStatTripletExc)
        uRadFin = sqrt(uRad**2 - uStatRad)
        uSingletPerIonPairFin = sqrt(uSingletPerIonPair**2 - uStatSingletPerIonPair)  
        uTripletPerIonPairFin = sqrt(uTripletPerIonPair**2 - uStatTripletPerIonPair)
        uB1SuIonPairFin = sqrt(uB1SuIonPair**2 - uStatB1SuIonPair)
        uC1PuIonPairFin = sqrt(uC1PuIonPair**2 - uStatC1PuIonPair) 

        !Calculate any additional uncertainties depending on the above values.
        !uWInv = uW/((data_input%energy*data_input%eV/mObIons)**2)
        !uB1SuPerIonPair = (mObB1SuExc/mObIons)*sqrt((uB1SuExc/mObB1SuExc)**2 + (uIons/mObIons)**2) 
        !uC1PuPerIonPair = (mObC1PuExc/mObIons)*sqrt((uC1PuExc/mObC1PuExc)**2 + (uIons/mObIons)**2)

        !Print these uncertainties to file.
        open(70, file='diagnostics.txt')
        !write(70,*) "Mean Variances"
        !write(70,*) "_______________________________________________________"
        !write(70,*) "Ions: ",  uStatIons
        !write(70,*) "Excites: ", uStatExcites 
        !write(70,*) "Gens: ", uStatGens 
        !write(70,*) "Dissociations: ", uStatDiss 
        !write(70,*) "B1SuExc: ", uStatB1SuExc 
        !write(70,*) "C1PuExc: ", uStatC1PuExc 
        !write(70,*) "singletExc: ", uStatSingletExc
        !write(70,*) "tripletExc: ", uStatTripletExc
        !write(70,*) "rad: ", uStatRad 
        !write(70,*) "singletPerIonPair: ", uStatSingletPerIonPair
        !write(70,*) "tripletPerIonPair: ", uStatTripletPerIonPair

        write(70,*) "Incident energy (eV): ", data_input%energy*data_input%eV
        write(70,*) "Mean Values"
        write(70,*) "-----------------------------------------------------"
        write(70,*) "mean number of ionisations: ", mObIons
        write(70,*) "mean energy per ion pair: ", mObW
        write(70,*) "mean number of excitations: ", mObExcites
!        write(70,*) "mean number of elastic scattering events: ", mObElastics
        write(70,*) "mean number of generations: ", mObGens
        write(70,*) "mean number of dissociations: ", mObDiss
        write(70,*) "mean number of B1Su excitations: ", mObB1SuExc
        write(70,*) "mean number of C1Pu excitations: ", mObC1PuExc
        write(70,*) "mean number of singlet excitations: ", mObSingletExc
        write(70,*) "mean number of triplet excitations: ", mObTripletExc
        write(70,*) "mean final radius: ", mObRad

        write(70,*) "Uncertainties"
        write(70,*) "-----------------------------------------------------"
        write(70,*) "u(ionisations): ", uIonsFin
        write(70,*) "u(mean energy per ion pair): ", uWFin
        write(70,*) "u(num excitations): ", uExcitesFin
!        write(70,*) "u(num elastic colls): ", uElasticsFin
        write(70,*) "u(num generations): ", uGensFin
        write(70,*) "u(num dissociations): ", uDissFin
       ! write(70,*) "u(1/mean energy per ion pair):", uWInvFin
        write(70,*) "u(B1SuPerIonPair):", uB1SuPerIonPairFin
        write(70,*) "u(C1PuPerIonPair):", uC1PuPerIonPairFin
        write(70,*) "u(num B1Su excitations): ", uB1SuExcFin
        write(70,*) "u(num C1Pu excitations): ", uC1PuExcFin
        write(70,*) "u(num singlet excitations): ", uSingletExcFin
        write(70,*) "u(num triplet excitations): ", uTripletExcFin
        write(70,*) "u(num singlet excitations per ion pair): ", uSingletPerIonPairFin
        write(70,*) "u(num triplet excitations per ion pair): ", uTripletPerIonPairFin
        write(70,*) "u(num B1Su excitations per ion pair): ", uB1SuIonPairFin
        write(70,*) "u(num C1Pu excitations per ion pair): ", uC1PuIonPairFin
        write(70,*) "u(final radius): ", uRadFin

        write(70,*) "Observed Standard Deviations"
        write(70,*) "_______________________________________________________" 
        write(70,*) "sigmaOb(ionisations): ", uIons
        write(70,*) "sigmaOb(mean energy per ion pair): ", uW
        write(70,*) "sigmaOb(num excitations): ", uExcites
!        write(70,*) "sigmaOb(num elastic colls): ", uElastics
        write(70,*) "sigmaOb(num generations): ", uGens
        write(70,*) "sigmaOb(num dissociations): ", uDiss
        !write(70,*) "sigmaOb(1/mean energy per ion pair):", uWInv
        write(70,*) "sigmaOb(B1SuPerIonPair):", uB1SuPerIonPair
        write(70,*) "sigmaOb(C1PuPerIonPair):", uC1PuPerIonPair
        write(70,*) "sigmaOb(num B1Su excitations): ", uB1SuExc
        write(70,*) "sigmaOb(num C1Pu excitations): ", uC1PuExc
        write(70,*) "sigmaOb(num singlet excitations): ", uSingletExc
        write(70,*) "sigmaOb(num triplet excitations): ", uTripletExc
        write(70,*) "sigmaOb(num singlet excitations per ion pair): ", uSingletPerIonPair
        write(70,*) "sigmaOb(num triplet excitations per ion pair): ", uTripletPerIonPair
        write(70,*) "sigmaOb(num B1Su excitations per ion pair): ", uB1SuIonPair
        write(70,*) "sigmaOb(num C1Pu excitations per ion pair): ", uC1PuIonPair
        write(70,*) "sigmaOb(final radius): ", uRad

        close(70)

    end subroutine calculateErrors
end module mc
