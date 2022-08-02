module ParticleDynamics
	use numbers
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    use sdcs_module
    use mc					! contains subroutines for monte carlo simulation
    use input_data			! contains input such as incident energy, benchmark mode
    use Ps_module
    use dcs_module         ! deals with elastic DCS
	
public :: E_Field, timeTaken

contains

subroutine E_Field(particleIn,rad,costheta,phi,coll,datasim,Elec,energy,statebasis)
       
	use numbers
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
		real(dp)::energy, energyJ
        real(dp)::rad, path, costheta, phi, theta
        real(dp)::xVal,yVal,zVal
		real(dp)::velx, vely, velz
		real(dp)::time, mass, charge, magVel
		real(dp)::deltaT
		real(dp), dimension(3)::Elec !E field, defined in input
		real(dp), dimension(3)::Acceleration, Velocity
        real(dp), dimension(3,3)::rotMat
        integer::coll

		!Figure out values that will be needed for calculation, acceleration in each direction due to E field, initial velocities
		theta = ACOS(costheta) 
		!Accelerations in x, y, z directions (lab frame)
		mass = 9.10938356E-31 !mass of an electron in kg
		charge = 1.60217663E-19 !charge of an electron in C
		energyJ = energy * 1.60218E-19
		Acceleration(1) = (Elec(1)*charge)/mass
		Acceleration(2) = (Elec(2)*charge)/mass
		Acceleration(3) = (Elec(3)*charge)/mass
		
		!print*, 'Acceleration(1) is: ', Acceleration(1)
		
		!Velocities in x, y, z directions after collision
		
		magVel = sqrt((2.0*energyJ)/mass)
		!print*,'Magnitude of velocity at collision ', coll, 'is: ', magVel
		
		Velocity(1) = magVel * costheta*sin(phi)
		Velocity(2) = magVel * sin(theta)*sin(phi)
		Velocity(3) = magVel * cos(phi)
		
		!Update particle velocity for given collision generation
		
		particleIn%velx(coll) = Velocity(1)
		particleIn%vely(coll) = Velocity(2)
		particleIn%velz(coll) = Velocity(3)
		particleIn%abs_vel(coll) = magVel

        !Begin by taking 'rad' calculated in selectPath, and figuring out time taken
		print*, 'Before timeTaken call', coll
		call timeTaken(rad, Velocity, Acceleration, time)
		print*, 'After timeTaken call', coll
		
		
		
		!Next, establish average kinetic energy of initial and final positions
		velx = Velocity(1) + Acceleration(1)*time
		vely = Velocity(2) + Acceleration(2)*time
		velz = Velocity(3) + Acceleration(3)*time
		
		magVel = sqrt(velx*velx + vely*vely + velz*velz)
		
		energyJ = (energyJ + 0.5*mass*magVel*magVel)/2.0
		
		energy = energyJ * 6.242E+18
		
	
		
		!Then, calculate path again using new approximated energy
		call selectPath(path,statebasis,energy)
		
		
		!Calculate new time taken, which will be the time returned to the simulation
		call timeTaken(path, Velocity, Acceleration, time)
		print *, 'Time taken for collision ', coll, 'is: ', time
		
		!Update position, position from origin, energy, for new collision
		velx = Velocity(1) + Acceleration(1)*time
		vely = Velocity(2) + Acceleration(2)*time
		velz = Velocity(3) + Acceleration(3)*time
		
		magVel = sqrt(velx*velx + vely*vely + velz*velz)
		
		energyJ = 0.5*mass*magVel*magVel      

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
		particleIn%energy(coll) = energy
		
		!Export time to main simulation
		deltaT = time
		particleIn%time(coll+1) = particleIn%time(coll) + deltaT

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
	print*, 'Before do loop, path is: ', path
	do while(integral .lt. path)
		
		deltax = abs(VelocityInt(1) * increment + Acceleration(1) * increment * increment)
		deltay = abs(VelocityInt(2) * increment + Acceleration(2) * increment * increment)
		deltaz = abs(VelocityInt(3) * increment + Acceleration(3) * increment * increment)
		
		integral = integral + sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz)
		
		VelocityInt(1) = VelocityInt(1) + Acceleration(1) * increment
		VelocityInt(2) = VelocityInt(2) + Acceleration(2) * increment
		VelocityInt(3) = VelocityInt(3) + Acceleration(3) * increment
		
		
		
		time = time + increment
		print *, "Current integral in loop is: ", integral
	end do
	print*, 'After do loop'
	
	!Export time taken
	time = time - increment

end subroutine timeTaken
end module ParticleDynamics
