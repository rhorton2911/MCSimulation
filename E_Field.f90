
subroutine E_Field(particleIn,rad,costheta,phi,coll,datasim,Elec,energy)
       
	   use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS

	   implicit none
        type(simdata)::datasim
        type(particle),intent(inout)::particleIn
		real(dp)::energy
        real(dp)::rad, path, costheta, phi, theta
        real(dp)::xVal,yVal,zVal
		real(dp)::velx, vely, velz
		real(dp)::time, mass, charge, magVel
		real(dp),intent(out)::deltaT
		real(dp), dimension(3)::Elec !E field, defined in input
		real(dp), dimension(3)::Acceleration, Velocity
        real(dp), dimension(3,3)::rotMat
        integer::coll

		!Figure out values that will be needed for calculation, acceleration in each direction due to E field, initial velocities
		theta = ACOS(costheta) 
		!Accelerations in x, y, z directions (lab frame)
		mass = 9.10938356E-31 !mass of an electron in kg
		charge = 1.60217663E-19 !charge of an electron in C
		Acceleration(1) = (Elec(1)*charge)/mass
		Acceleration(2) = (Elec(2)*charge)/mass
		Acceleration(3) = (Elec(3)*charge)/mass
		
		!Velocities in x, y, z directions after collision
		
		Velocity(1) = sqrt((2*energy)/mass) * costheta*sin(phi)
		Velocity(2) = sqrt((2*energy)/mass) * sin(theta)*sin(phi)
		Velocity(3) = sqrt((2*energy)/mass) * cos(phi)

        !Begin by taking 'rad' calculated in selectPath, and figuring out time taken
		time = timeTaken(rad, Velocity, Acceleration)
		
		
		
		!Next, establish what kinetic energy would be half way along this path using classical kinematics
		velx = Velocity(1) + Acceleration(1)*time*0.5
		vely = Velocity(2) + Acceleration(2)*time*0.5
		velz = Velocity(3) + Acceleration(3)*time*0.5
		
		magVel = sqrt(velx*velx + vely*vely + velz*velz)
		
		energy = 0.5*mass*magVel*magVel
		
		!Then, calculate path again using new approximated energy
		call selectPath(path,stateBasis,energy,coll)
		
		
		!Calculate new time taken, which will be the time returned to the simulation
		time = timeTaken(path, Velocity, Acceleration)
		
		
		!Update position, position from origin, energy, for new collision
		velx = Velocity(1) + Acceleration(1)*time
		vely = Velocity(2) + Acceleration(2)*time
		velz = Velocity(3) + Acceleration(3)*time
		
		magVel = sqrt(velx*velx + vely*vely + velz*velz)
		
		energy = 0.5*mass*magVel*magVel      

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
		
		particleIn%energy(coll) = energy
		
		!Export time to main simulation
		deltaT = time

end subroutine E_Field
	
	
subroutine timeTaken(path, Velocity, Acceleration)

    use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS

	implicit none
	
	real(dp)::a, b, c
	real(dp),intent(out)::time
	real(dp)::path, magVel, magAcc
	real(dp), dimension(3),intent(inout)::Acceleration, Velocity
	real(dp)::lower_bound, increment, integral
	real(dp)::arg1, arg2, arsinh1, arsinh2
	
	!Calculate a, b and c constants for integral from velocity, acceleration vectors
	a = Acceleration(1)*Acceleration(1) + Acceleration(2)*Acceleration(2) + Acceleration(3)*Acceleration(3)
	b = 2*Velocity(1)*Acceleration(1) + 2*Velocity(2)*Acceleration(2) + 2*Velocity(3)*Acceleration(3)
	c = Velocity(1)*Velocity(1) + Velocity(2)*Velocity(2) + Velocity(3)*Velocity(3)
	
	!Calculate magAcc and magVel required for lower_bound
	magAcc = sqrt(Acceleration(1)*Acceleration(1) + Acceleration(2)*Acceleration(2) + Acceleration(3)*Acceleration(3))
	magVel = sqrt(Velocity(1)*Velocity(1) + Velocity(2)*Velocity(2) + Velocity(3)*Velocity(3))
	
	!Calculate lower bound to integrate above by assuming particle travels in straight line perfectly aligned with E field
	lower_bound = (-magVel + sqrt(magVel*magVel + 2*magAcc*path))/(2*magAcc)
	
    !Do loop continues until calculated path exceeds given path
	integral = 0.0
	time = lower_bound
	increment = 1.0*10E-11
	do while(integral .lt. path)
		arg1 = (2*a*time + b)/(sqrt(4*a*c - b*b))
		arg2 = b/(sqrt(4*a*c - b*b))
		arsinh1 = log(arg1 + sqrt(arg1**2 + 1))
		arsinh2 = log(arg2 + sqrt(arg2**2 + 1))
		integral = ((4*a*c - b*b)* (arsinh1 - arsinh2))/(8*(a**(3/2))) + ((2*a*time + b) * sqrt(time*(a*time + b) + c) - b*sqrt(c))/(4*a)
		time = time + increment
	end do
	
	!Export time taken
	time = time - increment

end subroutine timeTaken

