subroutine mcsimulation(data_input,statebasis,sdcsBasis,eldcs,dicsBasis,VarPs,datamc,min_ein)

    use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data	! contains input such as incident energy, benchmark mode
	  use ParticleDynamics	
    use Ps_module
    use dcs_module         ! deals with elastic DCS
    use OMP_LIB
 
    implicit none
    
    ! Input Parameters
    real(dp):: en_incident
	  real(dp), dimension(3)::Elec
    type(input), intent(in):: data_input
    type(basis_state),intent(in):: statebasis
    type(basis_dcs),intent(in):: eldcs    ! keep elastic DCS here
    type(basis_dics),intent(in):: dicsBasis !Dissociative ionisation CS and DF  
    type(basis_sdcs)::sdcsBasis
    !!type(basis_vcs)::vcsBasis
                                 
    ! Monte Carlo Parameters
    integer:: simIndex,totalSims
    type(simdata):: datamc
    integer:: maxgen 
    real(dp)::runTime
    integer:: t1, t2, clockRate !For measuring runtime
    
    ! Individual Simulation Parameters  
    type(particle),dimension(0:1000):: particlebasis
    integer:: partNum,coll
    type(simdata):: datasim
    character (len=60):: ionop ! Ionistation option
    logical:: diss
    real(dp):: cutoffEn
    real(dp):: min_ein    

    ! For Ps Benchmark Simulation
    type(PsVar),intent(in):: VarPs 
    logical:: bmode ! Defines whether benchmark (true) or default sim (false) is being run
   
    !real(dp), external :: OMP_GET_WTIME 
    
    !---------------- Monte Carlo Simulation ------------------------------------------------------------------
    print*, 'commencing mcsimulation', data_input%energyeV
    totalSims = int(real(data_input%totalSims)/real(data_input%totalVariations))
    en_incident = data_input%energyeV
    bmode = data_input%benchmark
    call init_simdata(datamc) 
    maxgen = 0 
    !ionop = 'mean' ! uses the mean excitation energy to calculate energy loss of incident particle and energy of ejected electron
    ionop = 'positive ionisation energy' 
    ionop = data_input%ionop
    diss = .false.
	  Elec = (/1.0E-8, 1.0E-8, 1.0E-8/) !Hard code electric field for testing purposes
   
    
    datamc%numSims = totalSims
    open(unit=60,file='particletrace.txt') ! Output of each collision for every particle

    call system_clock(t1, clockRate)

    !$OMP PARALLEL DO private(particleBasis,diss,datasim,partNum,coll) 
    do simIndex = 1, totalSims
       ! initialise all energy values in the simulation to 0 so that the previous simulation is forgotten
       !call clearparticlebasis(particlebasis)		
	
       ! Initialise the number of electrons created in each generation to 0
       !call clearsimdata(datasim)
	
       !----------- Simulation starting at initial particle -------------------------------------
       partNum = 0 
       call init_simdata(datasim)
       call init_particle(particlebasis(0),en_incident,0,0.0_dp,0.0_dp,0.0_dp,0.0_dp)	
	
       !print*, '---------------------------------------'
       !print*, '---------------------------------------'
       !print*, '---------------------------------------'
       !print*, '---------------------------------------'
       !print*, "Thread Num: ", OMP_GET_THREAD_NUM()
       !print*, "SIMINDEX: ", simIndex
	
       do while(((partNum.LE.datasim%secE).and.(.not.bmode))  .or.  ((bmode).and.(partNum.eq.0)))	! Tracking of secondary electrons is turned off for the Ps Benchmark (bmode=true)	
          coll = 0	
      	  if(bmode .eqv. .true.) then ! Ps Formation Benchmark Simulation will be run instead of default simulation
   	     do while((particlebasis(partNum)%energy(coll).GE. 6.80) .AND. (datasim%PsFormed .EQV. .FALSE.)) 			
   	        particlebasis(partNum)%colls = coll ! Updates the amount of collisions the particle has gone through
   		datasim%collPerGen(particlebasis(partNum)%gen) = datasim%collPerGen(particlebasis(partNum)%gen) + 1 ! Updates collisions per generation				

		call collisionsimulation(data_input,statebasis,sdcsBasis,dicsBasis,particlebasis,partNum,coll,simIndex,datasim,ionop,bmode,VarPs,eldcs,Elec)
		!call timeElapsed(e_energy(x), scattEvent, cs_Ps, cs_el, cs_ion, cs_exc, duration)		
		!tResTemp = tResTemp + duration
				
		coll = coll + 1 ! Increment to consider the next collision					
	     end do			
	  else ! Run default simulation
             cutoffEn = min_ein  !Cutoff energy for the simulation
             !print*, "Cutoff energy: ", cutoffEn 
	     do while((particlebasis(partNum)%energy(coll).GE. cutoffEn) .AND. (particlebasis(partNum)%energy(coll) .LE. 700)) 				                               
	        particlebasis(partNum)%colls = coll ! Updates the amount of collisions the particle has gone through	
		datasim%collPerGen(particlebasis(partNum)%gen) = datasim%collPerGen(particlebasis(partNum)%gen) + 1	! Updates collisions per generation			
		call collisionsimulation(data_input,statebasis,sdcsBasis,dicsBasis,particlebasis,partNum,coll,simIndex,datasim,ionop,bmode,VarPs,eldcs,Elec)
		coll = coll + 1 ! Increment to consider the next collision		
	     end do
	  end if
		
	  ! Calculates the highest generation number reached in the Monte Carlo simulation
	  if(maxgen .LT. particlebasis(partNum)%gen) then
	     maxgen = particlebasis(partNum)%gen
	  end if
		
	  partNum = partNum + 1	! Increment to consider the next particle
       end do	 
 
       !!!!!!!!$OMP CRITICAL
       !Wait for threads to catch up to each other before updating data
       !print*, "COLLECTING DATA" 
       call collectdata(datasim,maxgen,datamc,data_input) ! Update Monte Carlo Parameters after each simulation	 
       !!!!!!!!!$OMP END CRITICAL	 
 
       !Collect spatial data for a few simulations (number given in input file)
       if (simIndex .le. data_in%numToWrite) then
          call writePathToFile(simIndex,particlebasis,datasim,bmode)
       end if

       if(simIndex .eq. 1) then
          call printsim(simIndex,particlebasis,datasim,bmode) ! Outputs trace for first simulation
       end if
       !Make arrays in particle type allocatable, then deallocate and end of simulation, reallocate as needed during the simulation.
       !deallocate(particlebasis)
    end do
    !$OMP END PARALLEL DO

    close(60) ! Close particletrace.txt

    call system_clock(t2, clockRate)
    runTime = real(t2-t1)/real(clockRate)
 
    call simulationresults(datamc,totalSims,en_incident,maxgen,ionop,bmode,runTime,data_input) 
    !call meanexcendist(statebasis) ! Do this if you want to produce a file that shows the distribution of sec e energies
 
    print*, 'Finished mcsimulation'
end subroutine mcsimulation



subroutine collisionsimulation(data_input,statebasis,sdcsBasis,dicsBasis,particlebasis,partNum,coll,simIndex,datasim,ionop,bmode,VarPs,eldcs,Elec)
    use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use sdcs_module
    use totalcs_module     !  reading totalcs files
    use input_data			! contains input such as incident energy, benchmark mode
    use mc				! contains subroutines for monte carlo simulation
	  use ParticleDynamics
    use Ps_module
    use dcs_module         ! deals with elastic DCS
    use OMP_LIB	
    use dcs_module

    implicit none
    type(input), intent(in)::data_input	 
    type(particle),dimension(0:1000):: particlebasis
    integer,intent(in):: partNum,coll
    integer:: simIndex
    type(simdata),intent(inout):: datasim
    type(basis_state),intent(in):: statebasis
    type(basis_dics),intent(in):: dicsBasis 
    character (len=60),intent(in):: ionop
    type(basis_dcs),intent(in):: eldcs    ! keep elastic DCS here
    type(basis_sdcs)::sdcsBasis
	  real(dp), dimension(3),intent(in)::Elec
    !type(basis_vcs)::vcsBasis
    type(totalcs):: tcs
    integer:: stateNum
    type(dcs):: dcsaten
    real(dp):: cosangle, randNum, phi, path	
    real(dp):: mass

    ! For Ps Benchmark Simulation
    type(PsVar),intent(in):: VarPs 
    logical,intent(in):: bmode ! Defines whether benchmark (true) or default sim (false) is being run
  
    mass = 9.11e-31  !Electron mass in SI units
 
    if(bmode .eqv. .false.) then ! Run Default Simulation	
       call get_csatein(statebasis,particlebasis(partNum)%energy(coll),tcs) ! Create totalcs for the current energy !call print_tcs(tcs) 
       !SYMMETRISE ONLY IN UPDATE_ENERGY, RESELECT POSITIVE ENERGY STATE
       if (data_input%stateIonop .eq. 1) then
          call symTcs(tcs,particlebasis(partNum)%energy(coll))                
       end if    
       call selectstate(tcs,particlebasis(partNum),coll,stateNum) ! Select a state based on totalcs	

       !call get_dcsatein(eldcs,particlebasis(partNum)%energy(coll),dcsaten)
       !particlebasis(partNum)%theta(coll) = dcsaten%avang
       if ( data_input%angleop .eq. 1 ) then
          !Use mean scattering angle
          call getAvAngle(eldcs,particlebasis(partNum)%energy(coll),cosangle)                                           
          !particlebasis(partNum)%costheta(coll) = cosangle ! Set the angle for the current collision
          !print*, particlebasis(partNum)%energy(coll)
       else if ( data_input%angleop .eq. 2 ) then
          !Select cos(angle) from differential cross sections(DCS)
          call selectAngle(eldcs,particlebasis(partNum)%energy(coll),cosangle)
          !particlebasis(partNum)%costheta(coll) = cosangle
       end if
       !Randomly generate an azimuthal scattering angle.
       call RANDOM_NUMBER(randNum)
       !particlebasis(partNum)%phi(coll) = randNum*4d0*atan(1.d0)                
       phi = randNum*4d0*atan(1.d0)   

       !Generate path length at given energy
       call selectPath(path,stateBasis,particleBasis(partNum)%energy(coll))
       !call update_position(particlebasis(partNum,rad,costheta,phi,coll)   !Update particle position

       call E_Field(particlebasis(partNum),path,cosangle,phi,coll,datasim,Elec,particleBasis(partNum)%energy(coll),statebasis)
       call update_energy(statebasis,sdcsBasis,eldcs,stateNum,tcs,particlebasis,partNum,cosangle,coll,ionop,data_input%enlossop,datasim,bmode,VarPs) ! Update the energy of the particle	


       !call testInputData(statebasis, eldcs, sdcsBasis, statebasis%b(1)%ein, size(statebasis%b(1)%ein))
	     !Where should E_Field go exactly?
	   
	     !The below line is commented out as it is the time taken in the case of no E field 
	     !deltaT = path*SQRT(mass/(2*particleBasis(partNum)%energy(coll)*data_in%evToSi)) 
       !particleBasis(partNum)%time(coll+1) = particlebasis(partNum)%time(coll) + deltaT
       call dissociation(particlebasis(partNum),statebasis,stateNum,dicsBasis,tcs,datasim,data_input,coll,partNum)

       ! Clean up memory
       call delete_totalcs(tcs)
       call delete_dcs(dcsaten)
    else ! Run Ps Benchmark Simulation	
       call create_totalcs(tcs,particlebasis(partNum)%energy(coll),VarPs) ! Create totalcs for the current energy !call print_tcs(tcs) 		
       call selectstate(tcs,particlebasis(partNum),coll,stateNum) ! Select a state based on totalcs	
	
       if(stateNum .eq. 4) then ! State 4 is Ps formation (simulation stops when ps is formed)
          datasim%PsFormed = .true.
	  particlebasis(partNum)%Ps(coll) = .true.
	  datasim%tRes = datasim%tRes + datasim%tResTemp
	  print*, datasim%tRes
       else 
          call update_energy(statebasis,sdcsBasis,eldcs,stateNum,tcs,particlebasis,partNum,cosangle,coll,ionop,data_input%enlossop,datasim,bmode,VarPs) ! Update the energy of the particle
	        call timeElapsed(particlebasis(partNum)%energy(coll),tcs%CS(stateNum),datasim)
	        datasim%tResTemp = datasim%tResTemp + datasim%duration
	        !print*, 'TResTemp', datasim%tResTemp
       end if		
	
       !print*, stateNum, tcs%CS(1), tcs%CS(2), tcs%CS(3), tcs%CS(4)

       ! Clean up memory
       call delete_totalcs(tcs) 	
    end if
	
    if(particlebasis(partNum)%energy(coll + 1)  .lt. 0) then ! Check energy loss is valid
       print*, 'invalid energy loss! (simulation.f90)'
       print*, 'At collision ',coll,'for particle ',partNum,' in sim ',simIndex,'the postCollEnergy was (eV) ',particlebasis(partNum)%energy(coll + 1)
    end if						
end subroutine collisionsimulation
