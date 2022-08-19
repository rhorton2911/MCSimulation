module AnalyticScattering
	use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS
	
	public :: ReidRampTCS, update_energy_analytic
	
	contains
	
	subroutine ReidRampTCS(statebasis,eincident,tcs)
	
		type(totalcs), intent(out):: tcs
		real(dp), intent(in):: eincident
		type(basis_state), intent(in):: statebasis
                integer::i,Nmax,tcstype

		
		if(eincident .lt. 0.2) then
			
			Nmax = 1 !Only elastic scattering possible
		
		else 
		
			Nmax = 2 !Toy inelastic channel open
		
		end if

                tcstype = statebasis%tcstype                
		call new_totalcs(tcs,Nmax,tcstype) !Create tcs structure
		
		!Populate parameters of tcs
		
		tcs%en_incident = eincident 
		tcs%tcstype = statebasis%tcstype
		
		do i=1,Nmax
			tcs%nf(i) = i
			tcs%ni(i) = 1  
			if(i .eq. 1) then !Populate total cross sections for interaction tydpe (GIVEN HERE IN AMSTRONG, MUST BE CONVERTED TO ATOMIC UNITS AT SOME STAGE)
				tcs%CS(i) = 6.0d0
			else if(i .eq. 2) then
				tcs%CS(i) = 10.0d0 * (eincident - 0.2d0)
			end if
			
			if(i .eq. 1) then !Populate energy of states
			tcs%en(i) = 0.0d0  !Energy of state with reference to some reference energy, set to zero for elastic, 0.2 for inelastic
			else if (i .eq. 2) then
			tcs%en(i) = 0.2d0
			end if
			tcs%df(i) = 0.0d0 !Dissociative fraction, zero for toy model
		enddo
	
	end subroutine

        subroutine update_energy_analytic(particlebasis, tcs, stateNum, partNum, datasim, coll, statebasis)
		use state_class
        use totalcs_module
        use Ps_module
        use sdcs_module
        use input_data
        use dcs_module
        implicit none

		type(basis_state),intent(in):: statebasis
		integer, intent(in):: stateNum
		type(particle),dimension(0:1000),intent(inout):: particlebasis
		integer,intent(in):: partNum, coll
		type(simdata), intent(inout):: datasim
		type(totalcs):: tcs
		real(dp):: elEnergyLoss	  ! energy lost by incident particle (eV)

	elEnergyLoss = 0.0d0

	if(stateNum==1) then !Elastic collision
		datasim%elastic = datasim%elastic + 1	! Update total elastic scattering count for this simulation
		datasim%elPerGen(particlebasis(partNum)%gen) = datasim%elPerGen(particlebasis(partNum)%gen) + 1 ! Update total elastic scattering count for this generation
		particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - elEnergyLoss !How should elEnergyLoss be calculated in analytic model?
	else if(stateNum==2) then !Toy inelastic model, energy decreases by 0.2
		 datasim%excite = datasim%excite + 1
	     datasim%excPerGen(particlebasis(partNum)%gen) = datasim%excPerGen(particlebasis(partNum)%gen) + 1 ! Update total excitation count for this generation
		 particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll) - 0.2d0

	end if

	end subroutine
	
end module AnalyticScattering
