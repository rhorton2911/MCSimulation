module AnalyticScattering
	use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS
	
	public :: update_energy_analytic
	
	contains
	subroutine update_energy_analytic(particlebasis, tcs, stateNum, partNum, datasim, coll, statebasis, costheta)
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
                real(dp),intent(in)::costheta
		real(dp):: elEnergyLoss	  ! energy lost by incident particle (eV)
                real(dp):: theta, mParticle, mTarget


                mParticle = 9.10938356E-31
                mTarget = 4.0d0 * 1.6726219E-27
                theta = acos(costheta)

                elEnergyLoss = (4.0d0*mParticle*mTarget)/((mParticle+mTarget)*(mParticle+mTarget))*(sin(theta/2))*(sin(theta/2)) !Taken from original Reid paper

	if(stateNum==1) then !Elastic collision
		datasim%elastic = datasim%elastic + 1	! Update total elastic scattering count for this simulation
		datasim%elPerGen(particlebasis(partNum)%gen) = datasim%elPerGen(particlebasis(partNum)%gen) + 1 ! Update total elastic scattering count for this generation
		particlebasis(partNum)%energy(coll + 1) = particlebasis(partNum)%energy(coll+1)* (1.0d0 - elEnergyLoss) !How should elEnergyLoss be calculated in analytic model?
	else if(stateNum==2) then !Toy inelastic model, energy decreases by 0.2
		 datasim%excite = datasim%excite + 1
	     datasim%excPerGen(particlebasis(partNum)%gen) = datasim%excPerGen(particlebasis(partNum)%gen) + 1 ! Update total excitation count for this generation
		 particlebasis(partNum)%energy(coll + 1) = (particlebasis(partNum)%energy(coll+1) - 0.2d0)*(1.0d0 - elEnergyLoss)

	end if

	end subroutine

end module AnalyticScattering
