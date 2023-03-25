!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module: lossFunc.f90
!Purpose: contains subroutines related to the energy loss function L(E),
!         a parameter used in modelling energy deposition, constructed 
!         from scattering cross sections. Old code no longer in use.
!Author: Reese Horton, Curtin University, Student ID: 19456155
!Date last modified: 09/02/2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module lossFunc




contains







     !______________________Get energy loss function and print to file_____!
     call lossFunction(statebasis,eldcs,sdcsbasis,data_in,enfine, Nenfine, lossFunc, lossFuncLen)  
     call lossFunctionPseudo(statebasis,eldcs,sdcsbasis,data_in,enfine, Nenfine, lossFuncPseudo, lossFuncLen) 
     open(70, file="lossFunc.txt") 
     write(70,*) "Energy(eV)       L(E) (SDCS)(cm^2eV)     L(E) (Psuedo)(cm^2eV)"
     do ii = 1, lossFuncLen
        write(70,*) enfine(ii), lossFunc(ii), lossFuncPseudo(ii)
     end do 
     close(70)
     deallocate(lossFuncPseudo,lossFunc)

     !____________________Get separate parts of L(E) and print to file_______!
     call lossFunctionSplit(statebasis,eldcs,sdcsbasis,data_in,enfine, Nenfine, lossFuncBound, lossFuncIon, lossFuncLen) 
     open(70, file="lossFuncSplit.txt") 
     write(70,*) "Energy(eV)       L(E) (Bound)(cm^2eV)     L(E) (Ion)(cm^2eV)"
     do ii = 1, lossFuncLen
        write(70,*) enfine(ii), lossFuncBound(ii), lossFuncIon(ii)
     end do 
     close(70)
     deallocate(lossFuncBound, lossFuncIon)


     call lossFunctionSplitPseudo(statebasis,eldcs,sdcsbasis,data_in,enfine, Nenfine, lossFuncBound, lossFuncIon, lossFuncLen) 
     open(70, file="lossFuncSplitPseudo.txt") 
     write(70,*) "Energy(eV)       L(E) (Bound)(cm^2eV)     L(E) (Ion)(cm^2eV)"
     do ii = 1, lossFuncLen
        write(70,*) enfine(ii), lossFuncBound(ii), lossFuncIon(ii)
     end do 
     close(70)
     deallocate(lossFuncBound, lossFuncIon)





   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !subroutine: lossFunction
   !Purpose: calculates and returns the energy loss distribution
   !         function L(E) used in the continuous slowing down
   !         approximation. Used to check validity of interpolated
   !         cross sections.
   !Date last modified: 19/03/2021
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine lossFunction(statebasis, eldcs, sdcsBasis, data_input, enfine, Nenfine, lossFunc, lossFuncLen)
       use state_class        ! defines states (and basis of them) with operations on them
       use totalcs_module     !  reading totalcs files
       use sdcs_module
       use input_data			! contains input such as incident energy, benchmark mode
       use dcs_module         ! deals with elastic DCS
       implicit none
       
       type(input), intent(in):: data_input
       type(basis_state),intent(in):: statebasis
       type(basis_dcs),intent(in):: eldcs    ! keep elastic DCS here
       type(basis_sdcs)::sdcsBasis
       real(dp), dimension(:), allocatable:: enfine
       integer:: Nenfine, ii, jj
       integer, intent(out):: lossFuncLen
       real(dp), dimension(:), allocatable:: lossFunc 
       !real(dp), dimension(:), allocatable:: sigmaSp
       real(dp):: massProj
       real(dp):: massTarget 
       type(totalcs)::tcs
       type(dcs)::dcsEn
       type(sdcs)::sdcsEn
       real(dp):: eIncident, eIon, intValTop, intValBot
       real(dp):: enCurr, sigmaMt, elCollLoss
       real(dp):: excEnSum, sdcsInt, cutoffEn 

       !sigmaSp is the stopping power cross section defined in 
       !Fursa, Zammit, et al. 2017- Electron mass stopping power in H2.
       !allocate(sigmaSp(Nenfine))

       !Definition for lossFunc = L(E) taken from: Dalgarno, Yan and Liu,
       !1999 - Electron energy deposition in a gas mxture of atomica and 
       !molecuar hydrogen and helium. Eqn 3.
       if (allocated(lossFunc)) then
          deallocate(lossFunc)
       end if

       allocate(lossFunc(Nenfine-1))
       lossFuncLen = Nenfine-1

       massProj = 9.11e-31  !Electron mass

       massTarget = (2.0)*1.6726219E-27      !Mass of H2
       eIon = 15.96632 !Ionisation energy of ground state H2 (eV)	
       do ii = 1, Nenfine-1
          eIncident = enfine(ii)
          if (eIncident .gt. statebasis%b(2)%enex) then
             !print*, eIncident
             call get_csatein(statebasis,eIncident, tcs)
             call get_dcsatein(eldcs, eIncident, dcsEn,1) 
             call get_sdcsatein(sdcsBasis, eIncident,sdcsEn) 
 
             sigmaMt = dcsEn%momtcs              !Momentum transfer CS
             !Sum over cross section weighted bound excitation energies
             excEnSum = 0.0
             jj = 1
             do while((tcs%en(jj+1) .lt. 0.0) .and. (jj .lt. tcs%Nmax-1)) 
                excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV
                jj = jj+1
             end do
             !excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV

             !Term used by dalgarno to model energy loss in elastic collisions
             !of fast electrons with ambient electrons.
             !ne =     !Fast electron number density
             !vEl = sqrt(2*eIncident/massProj)  !Fast electron velocity
             !Ee = 8.62e-3  !Value used by dalgarno
             !elCollLoss = ((2.0*1e-4)*(ne**0.97)/(eIncident)**0.44)*((eIncident-Ee)/(eIncident-0.53*Ee))**2.36       !Set to zero for initial testing
             !elCollLoss = elCollLoss/(data_input%density*vEl)
             elCollLoss = 0.0

             !Integrate over SDCS.
             sdcsInt = 0.0
             jj = 1
             cutoffEn = (eIncident -eIon)/2.0
             enCurr = sdcsEn%ecg(jj)
             do while((enCurr .le. cutoffEn) .and. (jj .le. sdcsEn%Ncgsdcs-1))
                !Use trapezoidal integration rule
                intValTop = (eIon + sdcsEn%ecg(jj+1))*sdcsEn%CS(jj+1)
                intValBot = (eIon + sdcsEn%ecg(jj))*sdcsEn%CS(jj)
                sdcsInt = sdcsInt + ((intValTop + intValBot)/2.0)*(sdcsEn%ecg(jj+1)-sdcsEn%ecg(jj))

                jj = jj+1
                enCurr = sdcsEn%ecg(jj)
             end do
 
             lossFunc(ii) = 2*(massProj*eIncident/massTarget)*sigmaMt + elCollLoss + excEnSum + sdcsInt 
             lossFunc(ii) = lossFunc(ii)*(data_input%bohrRadius*100)**2  !Print result in cm^2*eV
          end if
       end do 
   end subroutine lossFunction



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !subroutine: lossFunctionPseudo
   !Purpose: calculates and returns the energy loss distribution
   !         function L(E) used in the continuous slowing down
   !         approximation. Used to check validity of interpolated
   !         cross sections. Uses pseudostates in stead of SDCS.
   !Date last modified: 21/03/2021
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine lossFunctionPseudo(statebasis, eldcs, sdcsBasis, data_input, enfine, Nenfine, lossFunc, lossFuncLen)
       use state_class        ! defines states (and basis of them) with operations on them
       use totalcs_module     !  reading totalcs files
       use sdcs_module
       use input_data			! contains input such as incident energy, benchmark mode
       use dcs_module         ! deals with elastic DCS
       implicit none
       
       type(input), intent(in):: data_input
       type(basis_state),intent(in):: statebasis
       type(basis_dcs),intent(in):: eldcs    ! keep elastic DCS here
       type(basis_sdcs)::sdcsBasis
       real(dp), dimension(:), allocatable:: enfine
       integer:: Nenfine, ii, jj
       integer, intent(out):: lossFuncLen
       real(dp), dimension(:), allocatable:: lossFunc 
       !real(dp), dimension(:), allocatable:: sigmaSp
       real(dp):: massProj
       real(dp):: massTarget 
       type(totalcs)::tcs
       type(dcs)::dcsEn
       !type(sdcs)::sdcsEn
       real(dp):: eIncident, eIon, intValTop, intValBot
       real(dp):: enCurr, sigmaMt, elCollLoss
       real(dp):: excEnSum, sdcsInt, cutoffEn 

       !sigmaSp is the stopping power cross section defined in 
       !Fursa, Zammit, et al. 2017- Electron mass stopping power in H2.
       !allocate(sigmaSp(Nenfine))

       !Definition for lossFunc = L(E) taken from: Dalgarno, Yan and Liu,
       !1999 - Electron energy deposition in a gas mxture of atomica and 
       !molecuar hydrogen and helium. Eqn 3.
       if (allocated(lossFunc)) then
          deallocate(lossFunc)
       end if
       allocate(lossFunc(Nenfine-1))
       lossFuncLen = Nenfine-1

       massProj = 9.11e-31  !Electron mass
       massTarget = (2.0)*1.6726219E-27      !Mass of H2
       eIon = 15.96632 !Ionisation energy of ground state H2 (eV)	

       do ii = 1, Nenfine-1
          eIncident = enfine(ii)
          if (eIncident .gt. statebasis%b(2)%enex) then
             !print*, eIncident
             call get_csatein(statebasis,eIncident, tcs)
             call get_dcsatein(eldcs, eIncident, dcsEn,1) 
             !call get_sdcsatein(sdcsBasis, eIncident,sdcsEn) 
 
             sigmaMt = dcsEn%momtcs              !Momentum transfer CS
             !Sum over cross section weighted bound excitation energies
             excEnSum = 0.0
             jj = 1
             do while(jj .lt. tcs%Nmax-1)
                excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV
                jj = jj+1
             end do
             excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV

             !Term used by dalgarno to model energy loss in elastic collisions
             !of fast electrons with ambient electrons.
             !ne =     !Fast electron number density
             !vEl = sqrt(2*eIncident/massProj)  !Fast electron velocity
             !Ee = 8.62e-3  !Value used by dalgarno
             !elCollLoss = ((2.0*1e-4)*(ne**0.97)/(eIncident)**0.44)*((eIncident-Ee)/(eIncident-0.53*Ee))**2.36       !Set to zero for initial testing
             !elCollLoss = elCollLoss/(data_input%density*vEl)
             elCollLoss = 0.0

             !Integrate over SDCS.
             !sdcsInt = 0.0
             !jj = 1
             !cutoffEn = (eIncident -eIon)/2.0
             !enCurr = sdcsEn%ecg(jj)
             !do while((enCurr .le. cutoffEn) .and. (jj .le. sdcsEn%Ncgsdcs-1))
                !Use trapezoidal integration rule
             !   intValTop = (eIon + sdcsEn%ecg(jj+1))*sdcsEn%CS(jj+1)
             !   intValBot = (eIon + sdcsEn%ecg(jj))*sdcsEn%CS(jj)
             !   sdcsInt = sdcsInt + ((intValTop + intValBot)/2.0)*(sdcsEn%ecg(jj+1)-sdcsEn%ecg(jj))

             !   jj = jj+1
             !   enCurr = sdcsEn%ecg(jj)
             !end do
 
             lossFunc(ii) = 2*(massProj*eIncident/massTarget)*sigmaMt + elCollLoss + excEnSum 
             lossFunc(ii) = lossFunc(ii)*(data_input%bohrRadius*100)**2  !Print result in cm^2*eV
          end if
       end do 

   end subroutine lossFunctionPseudo



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !subroutine: lossFunctionSplit
   !Purpose: calculates and returns the energy loss distribution
   !         function L(E) used in the continuous slowing down
   !         approximation. Returns bound excitation part and 
   !         continuum part. Used to check validity of interpolated
   !         cross sections. Uses SDCS.
   !Date last modified: 24/03/2021
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine lossFunctionSplit(statebasis, eldcs, sdcsBasis, data_input, enfine, Nenfine, lossFuncBound, lossFuncIon, lossFuncLen)
       use state_class        ! defines states (and basis of them) with operations on them
       use totalcs_module     !  reading totalcs files
       use sdcs_module
       use input_data			! contains input such as incident energy, benchmark mode
       use dcs_module         ! deals with elastic DCS
       implicit none
       
       type(input), intent(in):: data_input
       type(basis_state),intent(in):: statebasis
       type(basis_dcs),intent(in):: eldcs    ! keep elastic DCS here
       type(basis_sdcs)::sdcsBasis
       real(dp), dimension(:), allocatable:: enfine
       integer:: Nenfine, ii, jj
       integer, intent(out):: lossFuncLen
       real(dp), dimension(:), allocatable:: lossFuncBound, lossFuncIon 
       !real(dp), dimension(:), allocatable:: sigmaSp
       real(dp):: massProj
       real(dp):: massTarget 
       type(totalcs)::tcs
       type(dcs)::dcsEn
       type(sdcs)::sdcsEn
       real(dp):: eIncident, eIon, intValTop, intValBot
       real(dp):: enCurr, sigmaMt, elCollLoss
       real(dp):: excEnSum, sdcsInt, cutoffEn 

       !sigmaSp is the stopping power cross section defined in 
       !Fursa, Zammit, et al. 2017- Electron mass stopping power in H2.
       !allocate(sigmaSp(Nenfine))

       !Definition for lossFunc = L(E) taken from: Dalgarno, Yan and Liu,
       !1999 - Electron energy deposition in a gas mxture of atomica and 
       !molecuar hydrogen and helium. Eqn 3.
       if (allocated(lossFuncBound)) then
          deallocate(lossFuncBound)
       end if
       if (allocated(lossFuncIon)) then
          deallocate(lossFuncIon)
       end if


       allocate(lossFuncBound(Nenfine-1), lossFuncIon(Nenfine-1))
       lossFuncLen = Nenfine-1

       massProj = 9.11e-31  !Electron mass
       massTarget = (2.0)*1.6726219E-27      !Mass of H2
       eIon = 15.96632 !Ionisation energy of ground state H2 (eV)	

       do ii = 1, Nenfine-1
          eIncident = enfine(ii)
          if (eIncident .gt. statebasis%b(2)%enex) then
             !print*, eIncident
             call get_csatein(statebasis,eIncident, tcs)
             call get_dcsatein(eldcs, eIncident, dcsEn,1) 
             call get_sdcsatein(sdcsBasis, eIncident,sdcsEn) 
 
             sigmaMt = dcsEn%momtcs              !Momentum transfer CS
             !Sum over cross section weighted bound excitation energies
             excEnSum = 0.0
             jj = 1
             do while((tcs%en(jj+1) .lt. 0.0) .and. (jj .lt. tcs%Nmax-1))
                excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV
                jj = jj+1
             end do
             !excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV

             !Term used by dalgarno to model energy loss in elastic collisions
             !of fast electrons with ambient electrons.
             !ne =     !Fast electron number density
             !vEl = sqrt(2*eIncident/massProj)  !Fast electron velocity
             !Ee = 8.62e-3  !Value used by dalgarno
             !elCollLoss = ((2.0*1e-4)*(ne**0.97)/(eIncident)**0.44)*((eIncident-Ee)/(eIncident-0.53*Ee))**2.36       !Set to zero for initial testing
             !elCollLoss = elCollLoss/(data_input%density*vEl)
             elCollLoss = 0.0

             !Integrate over SDCS.
             sdcsInt = 0.0
             jj = 1
             cutoffEn = (eIncident -eIon)/2.0
             enCurr = sdcsEn%ecg(jj)
             do while((enCurr .le. cutoffEn) .and. (jj .le. sdcsEn%Ncgsdcs-1))
                !Use trapezoidal integration rule
                intValTop = (eIon + sdcsEn%ecg(jj+1))*sdcsEn%CS(jj+1)
                intValBot = (eIon + sdcsEn%ecg(jj))*sdcsEn%CS(jj)
                sdcsInt = sdcsInt + ((intValTop + intValBot)/2.0)*(sdcsEn%ecg(jj+1)-sdcsEn%ecg(jj))

                jj = jj+1
                enCurr = sdcsEn%ecg(jj)
             end do
 
             lossFuncBound(ii) = excEnSum 
             lossFuncBound(ii) = lossFuncBound(ii)*(data_input%bohrRadius*100)**2  !Print result in cm^2*eV
             lossFuncIon(ii) = sdcsInt
             lossFuncIon(ii) = lossFuncIon(ii)*(data_input%bohrRadius*100)**2  !Print result in cm^2*eV
          end if
       end do 

   end subroutine lossFunctionSplit





   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !subroutine: lossFunctionSplitPseudo
   !Purpose: calculates and returns the energy loss distribution
   !         function L(E) used in the continuous slowing down
   !         approximation. Returns bound excitation part and 
   !         continuum part. Used to check validity of interpolated
   !         cross sections. Uses pseudostates in stead of SDCS.
   !Date last modified: 24/03/2021
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine lossFunctionSplitPseudo(statebasis, eldcs, sdcsBasis, data_input, enfine, Nenfine, lossFuncBound, lossFuncIon, lossFuncLen)
       use state_class        ! defines states (and basis of them) with operations on them
       use totalcs_module     !  reading totalcs files
       use sdcs_module
       use input_data			! contains input such as incident energy, benchmark mode
       use dcs_module         ! deals with elastic DCS
       implicit none
       
       type(input), intent(in):: data_input
       type(basis_state),intent(in):: statebasis
       type(basis_dcs),intent(in):: eldcs    ! keep elastic DCS here
       type(basis_sdcs)::sdcsBasis
       real(dp), dimension(:), allocatable:: enfine
       integer:: Nenfine, ii, jj
       integer, intent(out):: lossFuncLen
       real(dp), dimension(:), allocatable:: lossFuncBound, lossFuncIon 
       !real(dp), dimension(:), allocatable:: sigmaSp
       real(dp):: massProj
       real(dp):: massTarget 
       type(totalcs)::tcs
       type(dcs)::dcsEn
       type(sdcs)::sdcsEn
       real(dp):: eIncident, eIon, intValTop, intValBot
       real(dp):: enCurr, sigmaMt, elCollLoss
       real(dp):: excEnSum, sdcsInt, cutoffEn 

       !sigmaSp is the stopping power cross section defined in 
       !Fursa, Zammit, et al. 2017- Electron mass stopping power in H2.
       !allocate(sigmaSp(Nenfine))

       !Definition for lossFunc = L(E) taken from: Dalgarno, Yan and Liu,
       !1999 - Electron energy deposition in a gas mxture of atomica and 
       !molecuar hydrogen and helium. Eqn 3.
       if (allocated(lossFuncBound)) then
          deallocate(lossFuncBound)
       end if
       if (allocated(lossFuncIon)) then
          deallocate(lossFuncIon)
       end if

       allocate(lossFuncBound(Nenfine-1), lossFuncIon(Nenfine-1))
       lossFuncLen = Nenfine-1

       massProj = 9.11e-31  !Electron mass
       massTarget = (2.0)*1.6726219E-27      !Mass of H2
       eIon = 15.96632 !Ionisation energy of ground state H2 (eV)	

       do ii = 1, Nenfine-1
          eIncident = enfine(ii)
          if (eIncident .gt. statebasis%b(2)%enex) then
             !print*, eIncident
             call get_csatein(statebasis,eIncident, tcs)
             call get_dcsatein(eldcs, eIncident, dcsEn,1) 
             call get_sdcsatein(sdcsBasis, eIncident,sdcsEn) 
 
             sigmaMt = dcsEn%momtcs              !Momentum transfer CS
             !Sum over cross section weighted bound excitation energies
             excEnSum = 0.0
             jj = 1
             do while((tcs%en(jj+1) .lt. 0.0) .and. (jj .lt. tcs%Nmax-1))
                excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV
                jj = jj+1
             end do
             !excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV

             !Term used by dalgarno to model energy loss in elastic collisions
             !of fast electrons with ambient electrons.
             !ne =     !Fast electron number density
             !vEl = sqrt(2*eIncident/massProj)  !Fast electron velocity
             !Ee = 8.62e-3  !Value used by dalgarno
             !elCollLoss = ((2.0*1e-4)*(ne**0.97)/(eIncident)**0.44)*((eIncident-Ee)/(eIncident-0.53*Ee))**2.36       !Set to zero for initial testing
             !elCollLoss = elCollLoss/(data_input%density*vEl)

             elCollLoss = 0.0

             lossFuncBound(ii) = excEnSum 
             lossFuncBound(ii) = lossFuncBound(ii)*(data_input%bohrRadius*100)**2  !Print result in cm^2*eV

             excEnSum = 0.0
             if (jj .lt. tcs%Nmax-1) then
                do while((tcs%en(jj+1) .gt. 0.0) .and. (jj .lt. tcs%Nmax-1))
                   excEnSum = excEnSum + tcs%cs(jj+1)*abs((tcs%en(jj+1)-tcs%en(1)))*data_input%eV
                   jj = jj+1
                end do
                excEnSum = excEnSum + tcs%cs(jj+1)*abs(tcs%en(jj+1)-tcs%en(1))*data_input%eV
             end if

             lossFuncIon(ii) = excEnSum
             lossFuncIon(ii) = lossFuncIon(ii)*(data_input%bohrRadius*100)**2  !Print result in cm^2*eV
          end if
       end do 

   end subroutine lossFunctionSplitPseudo

end module lossFunc













