!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!File: sdcsOld.f90
!Purpose: stores old/redundant sdcs subroutines. 
!Author: Reese Horton, Curtin University, ID: 19456155
!Date last modified: 10/02/2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




module sdcsOld






contains




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: printSdcsAtEnFit
  !Purpose: prints cross sections in imported sdcs type to a file
  !         for plotting. Assumes sdcs are fitted from totalcs
  !Date last modified: 08/12/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine printSdcsAtEnFit(self)
      implicit none
      type(sdcs)::self
      integer:: ii
      character (len=3)::enString
 
      write(enString,'(I0)') int(self%en_incident)
      open(70, file="sdcsFit"//trim(enString)//".txt")
      write(70,*) "Energy(eV)          CS(a.u)"
      do ii = 1, self%Ncgsdcs
         write(70,*) self%ecg(ii), self%CS(ii)
      end do
      close(70)    

  end subroutine printSdcsAtEnFit

























    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: printDiffTcs
    !Purpose: numerically differentiates TCS at certain incident
    !          energies and prints to a file.
    !Date last modified: 01/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printDiffTcs(stateBasis)
        use state_class
        implicit none
        type(basis_state)::stateBasis
        integer::ii
        real(dp), dimension(6)::enArray
        real(dp)::eIn
        type(totalcs)::tcs 

        enArray = (/17.0,20.0,50.0,100.0,300.0,500.0/) !eV
        do ii = 1,6
           eIn = enArray(ii)
           call get_csatein(stateBasis,eIn,tcs)
           call symTcs(tcs,eIn)       
           call diffTcs(tcs)
           call printDiff(tcs,eIn)
           call delete_totalcs(tcs) 
        end do
     end subroutine printDiffTcs








    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: diffTcs
    !Purpose: numerically differentiates tcs to product approximate SDCS
    !         stores these SDCS in tcs object's cs array.
    !Date last modified: 01/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine diffTcs(tcs)
        use input_data
        implicit none
        type(totalcs)::tcs
        integer::ii,kk
        real(dp), allocatable, dimension(:)::tempArray

        allocate(tempArray(tcs%Nmax))
        tempArray(:) = 0.0
        tcs%en(:) = tcs%en(:)*data_in%eV

 
        do ii=2, tcs%Nmax
           !Need to avoid degenerate states
           if ((tcs%en(ii) .gt. 0.0) .and. (tcs%en(ii)-tcs%en(ii-1) .gt. 0.00000001)) then
              !Need to correct sign due to decreasing pseudostate CS
              tempArray(ii) = (-1)*(tcs%cs(ii)-tcs%cs(ii-1))/(tcs%en(ii)-tcs%en(ii-1))
           else
              !Find previous state with a different energy, in case of degenerate
              !states.
              kk = 1
              do while (tcs%en(ii) - tcs%en(ii-kk) .lt. 0.00000001)
                 kk = kk+1
              end do
              !Need to correct sign due to decreasing pseudostate CS
              tempArray(ii) = (-1)*(tcs%cs(ii)-tcs%cs(ii-kk))/(tcs%en(ii)-tcs%en(ii-kk))
           end if
        end do 
        tempArray(1) = tempArray(2)
        deallocate(tcs%cs)
        allocate(tcs%cs(tcs%Nmax))
        tcs%cs(:) = tempArray(:)
        deallocate(tempArray) 

    end subroutine diffTcs






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: printDiff
    !Purpose: prints the numerically differentiated pseudostate cross 
    !         sections to a file for plotting. 
    !         Should in principle be similar to SDCS
    !Date last modified: 01/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printDiff(tcs,eIn)
        implicit none
        type(totalcs)::tcs
        real(dp)::eIn
        character (len=3)::enString
        integer::ii

        write(enString,'(I0)') int(eIn)
        open(70, file="tcsDiff"//trim(enString)//"eV.txt")
        write(70,*) "EOut             CS"
        do ii=1, tcs%Nmax
           write(70,*) tcs%en(ii), tcs%cs(ii)
        end do
        close(70)

    end subroutine printDiff





!Functions for construction of sdcs from fitted pseudostate excitation CS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: get_sdcsateinfit
    !Purpose: returns a set of sdcs at the given incident energy by
    !         using pseudostate cross sections. These are fitted 
    !         to an exponential curve, then differentiated to give
    !         SDCS.
    !Date last modified: 07/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine get_sdcsateinfit(sdcsaten,sdcsBasis,eIncident,stateBasis,data_input)
        use state_class
        use input_data
        use totalcs_module 
        implicit none
        type(basis_sdcs)::sdcsBasis
        type(sdcs)::sdcsaten
        type(input)::data_input
        type(basis_state)::stateBasis
        real(dp)::eIncident, hartreeToEV
        type(totalcs)::tcs
        real(dp), allocatable, dimension(:)::enInArray       
        real(dp)::cutoffEn, deltaE, deltaCS, diffEn, eIon
        integer::cutoffInd, jj,kk
        logical::found 

        if (associated(sdcsaten%cs)) then
           deallocate(sdcsaten%cs)
        end if
        if (associated(sdcsaten%ecg)) then
           deallocate(sdcsaten%ecg)
        end if

        call get_csatein(stateBasis,eIncident,tcs)
        allocate(enInArray(sdcsBasis%Ncgsdcs))
        hartreeToEV = data_input%eV
        enInArray(:) = sdcsBasis%ecg(:)/hartreeToEV
        call fitTotalCs(tcs,enInArray,sdcsBasis%Ncgsdcs)
        !diffTcs numerically differentiates totalcs and stores the result in 
        !in the object's cross section array.
        deallocate(enInArray)
        call diffTcs(tcs)

        allocate(sdcsaten%cs(sdcsBasis%Ncgsdcs),sdcsaten%ecg(sdcsBasis%Ncgsdcs))
        sdcsaten%cs(:) = tcs%cs(:)
        sdcsaten%ecg(:) = tcs%en(:)
        sdcsaten%Ncgsdcs = tcs%Nmax

        call delete_totalcs(tcs)

        eIon = 15.96632  !Ionisation energy of ground state H2
        !Note: Need to check energy of each grid point above E/2 
        !      and perform linear interpolation over appropriate
        !      values below E/2 to get correct symmetry.
        cutoffEn = (eIncident-eIon)/2.0
        cutoffInd=1
        found = .false.
        do while (.not. found)
           if(sdcsaten%ecg(cutoffInd) .gt. cutoffEn) then
              found = .true.
           end if
           cutoffInd=cutoffInd+1
        end do
        cutoffInd = cutoffInd-1



        if (cutoffEn .gt. sdcsaten%ecg(1)) then
           kk = cutoffInd
           do while (sdcsaten%ecg(kk) .le. (eIncident-eIon))
              jj = cutoffInd-1
              found = .false.
              do while ((jj .ge. 2) .and. (found .eq. .false.))
                 if (sdcsaten%ecg(kk)-cutoffEn .lt. cutoffEn-sdcsaten%ecg(jj)) then
                    found = .true.
                 end if
                 jj = jj-1
              end do
              if (found .eq. .true.) then
                 !Linear interpolation over appropriate values below E/2
                 deltaE = sdcsaten%ecg(jj+1) - sdcsaten%ecg(jj)
                 deltaCS = sdcsaten%cs(jj+1) - sdcsaten%cs(jj)
                 diffEn = sdcsaten%ecg(kk)-cutoffEn - (cutoffEn-sdcsaten%ecg(jj+1))
                 ! '-' corrects sign of gradient, CS is decreasing
                 sdcsaten%cs(kk) = sdcsaten%cs(jj+1) - (deltaCS/deltaE)*diffEn
              end if
              kk = kk+1
           end do
        end if
    end subroutine get_sdcsateinfit




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: fitSdcsAll
    !Purpose: replaces all SDCS in sdcsBasis with SDCS fitted from totalcs
    !         pseudostate cross section values.
    !Date last modified: 15/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fitSdcsAll(sdcsBasis,stateBasis,data_input)
        use state_class
        use totalcs_module
        use input_data
        implicit none
        type(basis_sdcs)::sdcsBasis
        type(basis_state)::stateBasis
        type(totalcs)::tcs
        type(input)::data_input
        type(sdcs)::sdcsAtEIn
        integer::ii
        real(dp), allocatable, dimension(:,:)::tempBasis
        real(dp), allocatable, dimension(:)::enOutArray

        allocate(tempBasis(sdcsBasis%Ncgsdcs,sdcsBasis%Nein+1))

        do ii = 1, sdcsBasis%Nein
           call get_sdcsateinfit(sdcsAtEIn,sdcsBasis,sdcsBasis%ein(ii),stateBasis,data_input) 
           tempBasis(:,ii) = sdcsAtEIn%cs(:)  
           call destructSdcsAtEIn(sdcsAtEIn)
        end do
        sdcsBasis%cs(:,:) = tempBasis(:,:)
        deallocate(tempBasis)
    end subroutine fitSdcsAll











end module sdcsOld













