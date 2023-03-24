module ParticleDynamics
	use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS
	
public :: E_Field, timeTaken, DriftVelocityConvergence

contains

subroutine E_Field(particleIn,rad,costheta,phi,coll,datasim,Elec,energy,statebasis,tcs,ScatteringModel)
       
    use numbers
    use AnalyticScattering
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS

	   implicit none
	    type(basis_state),intent(in):: statebasis
        type(simdata)::datasim
        type(particle),intent(inout)::particleIn
		real(dp)::energy, energyJ, energyInit, energyJInit, energyOut
        real(dp)::rad, path, costheta, phi, theta, randNum, energySampled
        real(dp)::xVal,yVal,zVal
		real(dp)::velx, vely, velz
		real(dp)::time, mass, charge, magVel
		real(dp)::deltaT
		real(dp), dimension(3)::Elec !E field, defined in input
		real(dp), dimension(3)::Acceleration, Velocity
        real(dp), dimension(3,3)::rotMat
        integer::coll, ScatteringModel
        type(totalcs),intent(out):: tcs

		!Figure out values that will be needed for calculation, acceleration in each direction due to E field, initial velocities
		theta = ACOS(costheta) 

		!Accelerations in x, y, z directions (lab frame)
		mass = 9.10938356E-31 !mass of an electron in kg
		charge = 1.60217663E-19 !charge of an electron in C
		energyJ = energy * 1.60218E-19
                energyInit = energy
                energyJInit = energyJ
           
		Acceleration(1) = -(Elec(1)*charge)/mass
		Acceleration(2) = -(Elec(2)*charge)/mass
		Acceleration(3) = -(Elec(3)*charge)/mass

                !Test acceleration is what is expected
                if (mod(coll,1000) .eq. 0) then
                        !print*, 'Acceleration at collision ', coll, 'is: ', Acceleration
                end if


		!print*, 'Acceleration(1) is: ', Acceleration(1)
		
		!Velocities in x, y, z directions after collision
		
		magVel = sqrt((2.0*energyJ)/mass)
		!print*,'Magnitude of velocity at collision ', coll, 'is: ', magVel
		
		Velocity(1) = magVel * cos(phi)*sin(theta)
		Velocity(2) = magVel * sin(phi)*sin(theta)
		Velocity(3) = magVel * costheta
		particleIn%abs_vel(coll) = magVel

                !PathConvergence uses the mean free path to find convergent tcs, energy
                if(mod(coll, 10) .eq. 0) then
                !print*, 'Made it to coll: ', coll  
                end if
                call PathConvergence(rad,time,datasim,statebasis,Velocity,Acceleration,tcs,ScatteringModel,energyOut,energyJ, coll)
                if(mod(coll, 10) .eq. 0) then
                !print*, 'rad calculated by PathConvergence is: ', rad    
                end if    
                !Accurate path now generated in PathConvergence, selectPath merely used as first guess
                !call selectPath(path,statebasis,energyOut,ScatteringModel)
		!Calculate new time taken, which will be the time returned to the simulation
		call timeTaken(rad, Velocity, Acceleration, time)

		!print *, 'Time taken for collision ', coll, 'is: ', time
		!Update position, position from origin, energy, for new collision
		velx = Velocity(1) + Acceleration(1)*time
		vely = Velocity(2) + Acceleration(2)*time
		velz = Velocity(3) + Acceleration(3)*time
		
		magVel = sqrt(velx*velx + vely*vely + velz*velz)
		
		energyJ = 0.5*mass*magVel*magVel 
                !Need two energies, one which measures average energy and one which updates to simulation     
                call EnergySampling(particleIn, coll, Velocity, Acceleration, time, energySampled, mass) !Sample a random point along the path to take energy from, for calculation of average energy
                particleIn%energySampled(coll) = energySampled
  
        !Construct rotation matrix from theta, phi using analytic formulae
		!Unnecessary in E_Field, rotation elements included in velocity calculations
		
        !rotMat(1,1) = costheta*cos(phi)
        !rotMat(1,2) = -sin(theta)
        !rotMat(1,3) = costheta*sin(phi)
        !rotMat(2,1) = sin(theta)*cos(phi)
        !rotMat(2,2) = costheta
        !rotMat(2,3) = sin(theta)*sin(phi)
        !rotMat(3,1) = -sin(phi)
        !rotMat(3,2) = 0.0
        !rotMat(3,3) = cos(phi)

        !Rotate coordinate system. 
        particleIn%zHat = MATMUL(rotMat,particleIn%zHat)
        particleIn%yHat = MATMUL(rotMat,particleIn%yHat) 
        particleIn%xHat = MATMUL(rotMat,particleIn%xhat)

        !Calculate new position !!!NEED TO ADJUST FOR CURVED PATH!!!
		
		particleIn%x(coll+1) = particleIn%x(coll) + Velocity(1)*time + 0.5*Acceleration(1)*time*time
		particleIn%y(coll+1) = particleIn%y(coll) + Velocity(2)*time + 0.5*Acceleration(2)*time*time
		particleIn%z(coll+1) = particleIn%z(coll) + Velocity(3)*time + 0.5*Acceleration(3)*time*time

        !Below lines commented out, were originally from straight-trajectory "update_position"
		
		!particleIn%x(coll+1) = particleIn%x(coll) + rad*particleIn%zHat(1)
        !particleIn%y(coll+1) = particleIn%y(coll) + rad*particleIn%zHat(2)
        !particleIn%z(coll+1) = particleIn%z(coll) + rad*particleIn%zHat(3)

        xVal = particleIn%x(coll+1)
        yVal = particleIn%y(coll+1)
        zVal = particleIn%z(coll+1)

        !Update incident particle distance from origin
        datasim%inERad = SQRT(xVal**2+yVal**2+zVal**2)
		
		!Update particle energy - added to current coll value so update_energy functions properly
		energy = energyJ * 6.242E+18
                particleIn%energy(coll+1) = energy
                
        !Perform tcs 1 more time so that selectstate is given tcs at proper collision energy		 
        if(ScatteringModel .eq. 1) then !Use MCCC data
                call get_csatein(statebasis,energy,tcs) ! Create totalcs for the current energy !call print_tcs(tcs) 
        else if(ScatteringModel .eq. 2) then !Use Reid Ramp model
                call ReidRampTCS(statebasis,energy,tcs)
        end if

		!Export time to main simulation

                !Update velocities in particle structure
                !Multiply by 0.5 again so that the velocity is halfway along the path (average)

                particleIn%velx(coll) = (particleIn%x(coll+1) - particleIn%x(coll))/time
                particleIn%vely(coll) = (particleIn%y(coll+1) - particleIn%y(coll))/time 
                particleIn%velz(coll) = (particleIn%z(coll+1) - particleIn%z(coll))/time

		particleIn%time(coll+1) = particleIn%time(coll) + time

end subroutine E_Field
	
	
subroutine timeTaken(path, Velocity, Acceleration, time)

    use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS

	implicit none
	
	real(dp)::posx, posy, posz, deltax, deltay, deltaz
	real(dp),intent(inout)::time
	real(dp)::path, magVel, magAcc
	real(dp), dimension(3),intent(inout)::Acceleration, Velocity
	real(dp), dimension(3)::VelocityInt
	real(dp)::lower_bound, increment, integral
	real(dp)::arg1, arg2, arsinh1, arsinh2
	
	!Calculate a, b and c constants for integral from velocity, acceleration vectors
	
	!Calculate magAcc and magVel required for lower_bound
	
	!Calculate lower bound to integrate above by assuming particle travels in straight line perfectly aligned with E field
	
    !Do loop continues until calculated path exceeds given path
	integral = 0.0

	increment = 1.0E-11
	time = increment
	VelocityInt = Velocity
	posx = 0.0
	posy = 0.0
	posz = 0.0
	!print*, 'Before do loop, path is: ', path
	do while(integral .lt. path)
		
		deltax = abs(VelocityInt(1) * increment + 0.5 * Acceleration(1) * increment * increment)
		deltay = abs(VelocityInt(2) * increment + 0.5 * Acceleration(2) * increment * increment)
		deltaz = abs(VelocityInt(3) * increment + 0.5 * Acceleration(3) * increment * increment)
		
		integral = integral + sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz)
		
		VelocityInt(1) = VelocityInt(1) + Acceleration(1) * increment
		VelocityInt(2) = VelocityInt(2) + Acceleration(2) * increment
		VelocityInt(3) = VelocityInt(3) + Acceleration(3) * increment
		
		
		
		time = time + increment
		!print *, "Current integral in loop is: ", integral
	end do
	
	!Export time taken
	time = time - increment

end subroutine timeTaken

subroutine DriftVelocityConvergence(W, ConvergenceFlag, particleIn, coll)

	real(dp),dimension(3),intent(inout) :: W
	real(dp) :: W_new
	real(dp) :: W_old
	real(dp) :: Difference
	integer :: i
	integer,intent(inout) :: ConvergenceFlag
	type(particle),intent(inout)::particleIn
	integer,intent(in)::coll

	W_new = 0.0
	W_old = W(3)

	do i = coll-999, coll

		W_new = W_new + particleIn%velz(i)

	end do

	W_new = W_new/1000.0d0

	Difference = abs((W_old - W_new)/W_new)
        
        !print*, 'At collision ', coll, 'difference in drifts is: ', Difference


	if(Difference .lt. 0.50) then
		ConvergenceFlag = 1
	else
		ConvergenceFlag = 0
	end if

	W(3) = W_new

end subroutine DriftVelocityConvergence

subroutine PathConvergence(rad,time,datasim,statebasis,Velocity,Acceleration,tcs,ScatteringModel,energyOut,energyJ, coll)

	use numbers
    use AnalyticScattering
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS


	real(dp), dimension(3)::Acceleration, VelocityInt
	real(dp), dimension(3),intent(inout)::Velocity
	type(basis_state),intent(in):: statebasis
	type(totalcs),intent(inout):: tcs
    type(simdata),intent(in)::datasim
	integer::ScatteringModel, testCounter, coll
	real(dp)::energy,sigma,lambda_new,lambda_old,diff, magVel, mass, velx, vely, velz
        real(dp),intent(out)::energyOut
	real(dp),intent(out)::time
	real(dp),intent(inout)::rad
        real(dp)::energyJ, energyInit, energyFinal, E_grid, sigma_tot, counter, timeIncrement, randVal

	mass = 9.10938356E-31 !mass of an electron in kg
        
        !Generate total cross section as usual
        energy = energyJ * 6.242E+18
        randVal = 0.0
        !print*, 'energyJ is: ', energyJ
        energyInit = energy
        !print*, 'Velocity given to PathConvergence is: ', Velocity
	if(ScatteringModel .eq. 1) then !Use MCCC data
		call get_csatein(statebasis,energy,tcs) ! Create totalcs for the current energy !call print_tcs(tcs)
                sigma = SUM(tcs%cs)*(data_in%bohrRadius**2)
	else if(ScatteringModel .eq. 2) then !Use Reid Ramp model
		call ReidRampTCS(statebasis,energy,tcs)
                sigma = SUM(tcs%cs)*1.0E-20
	end if
        !print*, 'sigma old is: ', sigma
        !print*, 'sum of cross sections is: ', SUM(tcs%cs)
	!Calculate mean free path for first guess
        lambda_old = 1/(data_in%density*sigma)
        if(mod(coll, 10) .eq. 0) then
           !print*, 'Initial guess for path is: ', lambda_old
        end if

	diff = 1.0
        testCounter = 0
        timeIncrement = 1.0E-11
	!Run do loop to calculate convergent value of sigma and lambda
	do while(diff.gt.0.001)
                testCounter = testCounter + 1
                if(testCounter .gt. 10) then
                   print*, 'At collision: ', coll, 'testCounter is at: ', testCounter
                end if
		!At given guess for lambda, call timeTaken
		call timeTaken(lambda_old, Velocity, Acceleration, time)
                !print*, 'Time returned from timeTaken is: ', time
                !print*, 'Velocity returned from timeTaken is: ', Velocity

                velx = Velocity(1) + Acceleration(1)*time
		vely = Velocity(2) + Acceleration(2)*time
		velz = Velocity(3) + Acceleration(3)*time

                !print*, 'Velocity given to energyFinal is: ', velx, vely, velz

		magVel = sqrt(velx*velx + vely*vely + velz*velz)

                !print*, 'magVel is: ', magVel

		energyJ = 0.5*mass*magVel*magVel
		
		energyFinal = energyJ * 6.242E+18
		!print*, 'energyFinal is: ', energyFinal
		!Take new cross section at energy halfway along the path
                !print*, 'energyInit is: ', energyInit
		sigma_tot = 0.0d0
		counter = 0.0d0
                
                if(energyInit.le.energyFinal) then
                E_grid = energyInit
		do while(E_grid.le.energyFinal)
		
			if(ScatteringModel .eq. 1) then !Use MCCC data
				call get_csatein(statebasis,E_grid,tcs) ! Create totalcs for the current energy !call print_tcs(tcs) 
				sigma = SUM(tcs%cs)*(data_in%bohrRadius**2)  !Convert to SI
			else if(ScatteringModel .eq. 2) then !Use Reid Ramp model
				call ReidRampTCS(statebasis,E_grid,tcs)
				sigma = SUM(tcs%cs)*1.0E-20
			end if
			sigma_tot = sigma_tot + sigma
			counter = counter + 1.0d0
			!Calculate next E_grid in terms of time step
                        !VelocityInt(1) = VelocityInt(1) + Acceleration(1) * timeIncrement 
                        !VelocityInt(2) = VelocityInt(2) + Acceleration(2) * timeIncrement
                        !VelocityInt(3) = VelocityInt(3) + Acceleration(3) * timeIncrement

                        !magVel = sqrt(VelocityInt(1)**2 + VelocityInt(2)**2 + VelocityInt(3)**2)
                        !E_grid = 0.5*mass*magVel*magVel
                        
                        !Calculate next E_grid in terms of E increment
                        E_grid = E_grid + 0.0001
		end do
                end if 
                
                if(energyInit.gt.energyFinal) then
                E_grid = energyFinal
                 do while(E_grid.lt.energyInit)

                        if(ScatteringModel .eq. 1) then !Use MCCC data
                                call get_csatein(statebasis,E_grid,tcs) ! Create totalcs for the current energy !call print_tcs(tcs)
                                sigma = SUM(tcs%cs)*(data_in%bohrRadius**2)  !Convert to SI
                        else if(ScatteringModel .eq. 2) then !Use Reid Ramp model
                                call ReidRampTCS(statebasis,E_grid,tcs)
                                sigma = SUM(tcs%cs)*1.0E-20
                        end if
                        sigma_tot = sigma_tot + sigma
                        counter = counter + 1.0d0
                        
                         !Calculate next E_grid in terms of time step
                        !VelocityInt(1) = VelocityInt(1) + Acceleration(1) * timeIncrement
                        !VelocityInt(2) = VelocityInt(2) + Acceleration(2) * timeIncrement
                        !VelocityInt(3) = VelocityInt(3) + Acceleration(3) * timeIncrement

                        !magVel = sqrt(VelocityInt(1)**2 + VelocityInt(2)**2 + VelocityInt(3)**2)
                        !E_grid = 0.5*mass*magVel*magVel

                        !Calculate E_grid in terms of E step
                        E_grid = E_grid + 0.0001

                end do
                end if

		
		sigma = sigma_tot/counter
                !sigma = sigma_tot
                !print*, 'sigma new is: ', sigma
		!Calculate new lambda based on new TCS
		lambda_new = 1/(data_in%density*sigma)

		!Calculate relative difference in lambda's
                diff = abs((lambda_new - lambda_old)/(lambda_new))
                if(testCounter .gt. 10) then
                print*, 'lambda_old is: ', lambda_old
                print*, 'lambda_new is: ', lambda_new
                print*, 'Diff is: ', diff
                end if
                !Set old lambda for next iteration to be current new lambda
		lambda_old = lambda_new
	end do
        
        energyOut = energy
        !print *, 'Number of lambda loops is: ', testCounter

        !Accurate sigma now known, sample path length from this sigma (bypassing selectPath)
        call RANDOM_NUMBER(randVal)
           rad = (-1.0)*(1/(data_in%density*sigma))*log(randVal)
           
           if(randVal .gt. 0.99) then
             ! print*, 'Test counter number for path convergence loops: ', testCounter
           end if

end subroutine PathConvergence

subroutine EnergySampling(particleIn, coll, Velocity, Acceleration, time, energySampled, mass)
use numbers
    use AnalyticScattering
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS

  
    type(particle),intent(inout)::particleIn
    integer::coll
  real(dp), dimension(3)::Acceleration, Velocity, VelocitySampled
  real(dp)::time, energySampled, mass
  real(dp)::randNum, magVel, checkPoint, timeSample

!We sample the energy of the particle at every time interval of 10^7 s


  checkPoint = REAL(CEILING(particleIn%time(coll)*1.0E+7)) * 1.0E-7
  if(mod(coll, 100) .eq. 0) then
     !print*, 'checkPoint is: ', checkPoint, 'Time at end is: ', particleIn%time(coll) + time
  end if
  energySampled = 0.0d0

  if(particleIn%time(coll) + time .gt. checkPoint) then
     timeSample = checkPoint - particleIn%time(coll)

     VelocitySampled(1) = Velocity(1) + Acceleration(1)*timeSample
  VelocitySampled(2) = Velocity(2) + Acceleration(2)*timeSample
  VelocitySampled(3) = Velocity(3) + Acceleration(3)*timeSample
  magVel = sqrt(VelocitySampled(1)*VelocitySampled(1) + VelocitySampled(2)*VelocitySampled(2) + VelocitySampled(3)*VelocitySampled(3))

     

  !call RANDOM_NUMBER(randNum)
  !VelocitySampled(1) = Velocity(1) + Acceleration(1)*time*randNum
  !VelocitySampled(2) = Velocity(2) + Acceleration(2)*time*randNum
  !VelocitySampled(3) = Velocity(3) + Acceleration(3)*time*randNum
  !magVel = sqrt(VelocitySampled(1)*VelocitySampled(1) + VelocitySampled(2)*VelocitySampled(2) + VelocitySampled(3)*VelocitySampled(3))
  energySampled = 0.5*mass*magVel*magVel
  energySampled = energySampled * 6.242E+18 !Convert sampled energy to eV
  end if

  !IF THE ENERGY IS ZERO IT IS SKIPPED OVER IN FINAL CALCULATION OF MEAN ENERGY

  

end subroutine EnergySampling
  
end module ParticleDynamics
