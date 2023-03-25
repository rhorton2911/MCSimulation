module Ps_module

        use numbers

	type, public:: PsVar
		! Parameters for Cross Sections
		real(dp) :: e_Ps, p_Ps, lamda_Ps, a_Ps, cs_Ps
		real(dp) :: e_el, p_el, lamda_el, a_el, cs_el
		real(dp) :: e_ion, p_ion, lamda_ion, a_ion, cs_ion
		real(dp) :: e_exc, p_exc, lamda_exc, a_exc, cs_exc	
		real(dp) :: Q
		
	end type PsVar
	
	
	contains
	
	
	subroutine initialisePsVar(VarPs)
		implicit none 
		type(PsVar),intent(inout):: VarPs
		
		VarPs%Q = 1.0 ! Fraction of energy that the positron retains after ionisation
		
		! UNITS
		! Threshold energy (eV) ! p (unitless) ! lamda (eV) ! 'A' (Angstrom^2)
		
		! Ps
		VarPs%e_Ps = 6.8  !13.6 - 6.8 = 6.8 
		VarPs%p_Ps = 0.5 
		VarPs%lamda_Ps = 7.0
		VarPs%a_Ps = 1.0
		
		! Elastic		
		VarPs%e_el = -20.0 
		VarPs%p_el = 2.0
		VarPs%lamda_el = 20.0
		VarPs%a_el = 0.4 * VarPs%a_Ps
		
		! Ionisation
		VarPs%e_ion = 13.6
		VarPs%p_ion = 1.0
		VarPs%lamda_ion = 30.0 ! Varies on benchmarks. A = 30
		VarPs%a_ion = 0.3 ! Varies on benchmarks. A = 0.3
		
		! Excitation
		VarPs%e_exc = 10.0 ! Varies on benchmarks. A = N/A, B = 10, C = 2		
		VarPs%p_exc = 1.5		 
		VarPs%lamda_exc = 12.0 ! Varies on benchmarks. A = N/A B = 12 C = 10
		VarPs%a_exc = 0.5 ! Varies on benchmarks. A = N/A B = 0.5 C = 1.0
		
		! Initialise to 0 for debugging
		VarPs%cs_Ps = 0d0
		VarPs%cs_el = 0d0
		VarPs%cs_ion = 0d0
		VarPs%cs_exc = 0d0
	end subroutine initialisePsVar
		
		
	subroutine create_totalcs(self,en_incident,VarPs)
		use totalcs_module 
		implicit none
		type(totalcs),intent(inout):: self
		real(dp),intent(in):: en_incident 
		type(PsVar):: VarPs
		integer:: Nmax
		
		Nmax = 4
		
		self%en_incident = en_incident	! Set the incident energy for the tcs object
		self%Nmax = Nmax ! There are 4 possible states
		self%tcstype = 3 ! Used for this analytical positron
		
		if(Nmax .ne. 0)  then
			allocate(self%nf(Nmax),self%ni(Nmax),self%BICS(Nmax),self%CS(Nmax),self%en(Nmax),self%df(Nmax))
		endif
		
		! initialise cs for all states to be 0
		self%CS(1) = 0d0
		self%CS(2) = 0d0
		self%CS(3) = 0d0
		self%CS(4) = 0d0
		
		!print*,'en_incident,VarPs%e_el,VarPs%p_el,VarPs%lamda_el,VarPs%a_el',en_incident,VarPs%e_el,VarPs%p_el,VarPs%lamda_el,VarPs%a_el
		
		self%CS(1) = CS_SPWR(en_incident,VarPs%e_el,VarPs%p_el,VarPs%lamda_el,VarPs%a_el)		! Elastic Scattering 
		self%CS(2) = CS_SPWR(en_incident,VarPs%e_exc,VarPs%p_exc,VarPs%lamda_exc,VarPs%a_exc)	! Excitation
		self%CS(3) = CS_SPWR(en_incident,VarPs%e_ion,VarPs%p_ion,VarPs%lamda_ion,VarPs%a_ion)	! Ionisation
		self%CS(4) = CS_SEXP(en_incident,VarPs%e_Ps,VarPs%p_Ps,VarPs%lamda_Ps,VarPs%a_Ps)		! Ps Formation
		
		!print*, 'creating tcs for Ps Benchmark Simulation, cs are below:'
		!print*, self%CS(1),self%CS(2),self%CS(3),self%CS(4)
		
	end subroutine create_totalcs
	
	! Used to calculate A*S_exp function found in the positronium paper
	function CS_SEXP(particleEnergy, eThresh, p , lamda, a) result(crossSection)
                implicit none
		real(dp), INTENT(IN) :: particleEnergy, eThresh, p , lamda, a
		real(dp):: crossSection
		real(dp) :: x, e	

		! Need to fix this and find actual function
		e = 2.718281828459
		x = particleEnergy - eThresh

		! Follows the S_exp function given in the paper
		if(x >= 0) then		  
			crossSection = a * (((e/lamda)**p) * (x**p) * e**(- (p * x) / lamda))
		else		  
			crossSection = 0d0			
		end if		

	end function CS_SEXP


	! Used to calculate A*S_pwr function found in the positronium paper
	function CS_SPWR(particleEnergy, eThresh, p , lamda, a) result(crossSection)
                implicit none
		real(dp),intent(in) :: particleEnergy, eThresh, p , lamda, a
		real(dp):: crossSection
		real(dp):: x, e	

		! Need to fix this and find actual function
		e = 2.718281828459
		x = particleEnergy - eThresh
		
		! Follows the S_pwr function given in the paper
		if(x >= 0) then		  
			crossSection = a * (lamda**p) * ((p+(1.0))**(p+(1.0))) * (x / ((x + (p*lamda))**(p+(1.0))))
		else		  
			crossSection = 0.0			
		end if
		
		!print*, 'calculated cross section: ',crossSection		

	end function CS_SPWR
	
	subroutine populatePsstates(statebasis,VarPs)
		use state_class
		implicit none
		type(basis_state), intent(inout):: statebasis
		type(PsVar),intent(in):: VarPs
		
		statebasis%b(1)%enex = 0d0 ! Elastic Scattering,energy loss must be calculated as a function of scattering angle
		statebasis%b(2)%enex = VarPs%e_exc ! Excitation energy loss equal to threshold energy
		statebasis%b(3)%enex = VarPs%e_ion! Ionisation energy loss equal to threshold energy
		statebasis%b(4)%enex = 0d0 ! Ps Formation stops the simulation. Excitation energy is not relevant
		
		! Only purpose of setting en (ionisation energy) is so that the general update_energy subroutine knows which one is ionisation
		statebasis%b(1)%en = -1.0! Elastic Scattering,energy loss must be calculated as a function of scattering angle
		statebasis%b(2)%en = -1.0 ! Excitation energy loss equal to threshold energy
		statebasis%b(3)%en = 1.0! Ionisation energy loss equal to threshold energy
		statebasis%b(4)%en = -1.0! Ps Formation stops the simulation. Excitation energy is not relevant
		
	end subroutine populatePsstates
	
	
	! SUBROUTINE timeElapsed(energy,stateNum,datasim)

	! IMPLICIT NONE
	! REAL, INTENT(IN) :: energy
	! INTEGER, INTENT(IN) :: stateNum
	! type(PsVar) :: VarPs
	! REAL*8 :: m_p, eVToJoules, crossSection, velocity, constant, con2, m
	! !REAL*8, INTENT(OUT) :: duration
	! type(simdata),intent(inout):: datasim
	
	! m_p = 	9.109E-31 
	! m = 9.1093897
	
	! m_p = 9.10938356E-31 ! mass of positron
	
	! eVToJoules = 1.60218E-19
		
	! IF(stateNum .EQ. 4) THEN
		! crossSection = varPs%cs_Ps
	! ELSE IF(stateNum .EQ. 1) THEN
		! crossSection = varPs%cs_el
	! ELSE IF(stateNum .EQ. 3) THEN
		! crossSection = varPs%cs_ion
	! ELSE IF(stateNum .EQ. 2) THEN
		! crossSection = varPs%cs_exc
	! ELSE 
		! crossSection = 0.0
	! END IF
	
	! ! Convert Cross section from A^2 to m^2
	! !crossSection = crossSection * 10**(-20)
	
	! ! constant is equal to eVToJoules/m_p
	! constant = 1.7565E11
	
	! con2 = 1!*10**(14)
	
	! !velocity = (2 * energy * eVToJoules / m_p) ** (0.5)
	
	! !duration = 1 / (crossSection * velocity)	
	
	! !duration = ((m_p)**(1/2)) / (crossSection * ((2 * energy * eVToJoules)**(1/2)))
	! !duration = ((m_p)**(1/2)) / (((2 * energy * eVToJoules)**(1/2)))
	
	! !duration = (m**0.5)/((2*energy * 1.60218)**(0.5) * crossSection) / 10 ! in terms of s*m^-3 * 10^15
	! !duration = ((m/(2*energy*1.60218))**(0.5))/(crossSection*10*3.2) ! in terms of s*m^-3 * 10^15
	
	! ! THis is what I had before the update
	! !velocity = (2 * energy * constant) ** (0.5)
	! !duration = ((m/(2*energy*1.602176565))**(0.5))/(crossSection*10) ! in terms of s*m^-3 * 10^15
	
	
	
	! ! Update 03/07/2019 !
	
	! datasim%duration = sqrt(m_p / (2*energy*eVToJoules)) * 10E20 / crossSection
	
	! !write(2,*) 'vel', velocity, 'cs', crossSection, 'energy', energy, 'duration', duration
	

	
	
! END SUBROUTINE timeElapsed


end module Ps_module

