module sdcs_module

  use numbers

  type, public:: sdcs
     ! private
     integer:: tcstype
     real(dp):: en_incident  ! incident energy 
     real(dp):: en_ionization  ! en_ionization
     real(dp):: en_total  ! total energy 
     integer:: Ncgsdcs  ! number of energy grid points for SDCS
     integer:: nFinal !Index of final non-zero CS defined on common grid 
     real(dp), pointer, dimension(:):: CS => NULL()
     real(dp), pointer, dimension(:):: ecg => NULL()    ! common grid for ejected electron energies
  end type sdcs

  type, public:: basis_sdcs
     !     private
     integer:: tcstype
     real(dp), dimension(:,:), pointer:: CS => NULL()
     real(dp), pointer, dimension(:):: ein => NULL()    ! incident energies
     real(dp), pointer, dimension(:):: ecg => NULL()    ! common grid for ejected electron energies
     !Below needed for file IO due to varying energy grid sizes in files
     real(dp), pointer, dimension(:,:):: readEOut => NULL() 
     integer, pointer, dimension(:):: nBelowInFile => NULL() !Number of CS below enfine(1) in a file

     integer:: Ninsdcs   ! the last point on the fine energy grid where SDCS is available
     integer:: Ninsdcs_first ! the first point on the fine energy grid where SDCS is available 
     integer:: Ncgsdcs  ! number of ejection energy grid points for SDCS 
     integer:: Nsdcsatein !Number of incident energy grid points
     integer:: Ncgbelow !Number of ejection energy grid points below 4.4 eV

     real(dp), pointer, dimension(:,:):: CSDist => NULL()
     real(dp), pointer, dimension(:)::ecgDist => NULL()
     integer:: NcgsdcsDist     

     integer:: Nein 
     integer:: einLast !the last point on the fine energy grid of incident energies
     integer:: einFirst !first point on fine energy grid of incident energies
     real(dp):: en_ionization

  end type basis_sdcs
  !  

contains



  subroutine make_sdcs(self,Nenfine,enfine)

    use input_data         ! definitions for input data type, construct routine, object data_in

    implicit none

    type(basis_sdcs):: self
    integer, intent(in):: Nenfine
    real(dp), dimension(Nenfine), intent(in):: enfine

    integer:: Nlarge
    parameter(Nlarge = 1000)    
    real(dp):: eioniz
    real(dp), dimension(Nlarge,data_in%Nsdcs):: cs,en
    real(dp), dimension(data_in%Nsdcs):: ein
    integer, dimension(data_in%Nsdcs):: inumsdcs
    character(len=40):: tmpstr
    integer:: incident_energy, i, j, ns, Ncg, jm
    real(dp):: tmp1, tmp2, etotd2, step
    real(dp), dimension(:,:), allocatable:: cscg
    real(dp), dimension(Nenfine):: venfine
    ! to be placed in the type sdcs
    integer:: Ninsdcs, Ninsdcs_first, Ncgsdcs
    real(dp), dimension(:,:), allocatable:: sdcs
    real(dp), dimension(:), allocatable:: ensdcs
    real(dp), dimension(:), allocatable:: ecg, vcg
    !
    
    eioniz = self%en_ionization
    print*,'eioniz = self%en_ionization =' , eioniz

    print*, 'Read SDCS files: Nsdcs =', data_in%Nsdcs
    do i=1,data_in%Nsdcs

       tmpstr = TRIM(data_in%filename_sdcs(i)(9:))
       !     print*, '1. tmpstr = ', tmpstr
       ns = SCAN(tmpstr,'_')
       tmpstr = tmpstr(1:ns-1)
       read(tmpstr,*) incident_energy
!            print*, '2. tmpstr = ', tmpstr
       ein(i) = incident_energy

       open(27,file=data_in%filename_sdcs(i),action='read')
       read(27,*)
       j = 0
       do 
          j = j + 1
          read(27,*,END=127,ERR=127) en(j,i), tmp1, tmp2, cs(j,i)
          !        print*, en(j,i)
       enddo
127    continue
       close(27)
       inumsdcs(i) = j - 1
       print*, data_in%filename_sdcs(i), ein(i), inumsdcs(i) 

    enddo

    ! interpolate these SDCS to a standard grid: [0,1]

    ! set a common grid
    Ncg = 200  ! number of points in the common grid
    allocate(ecg(Ncg),vcg(Ncg))
    allocate(cscg(Ncg,data_in%Nsdcs))
    step = 1d0/dble(Ncg)
    print*,'SDCS: set common grid with number of points Ncg and step =', Ncg,step
    do i=1,Ncg
       ecg(i) = i*step
    enddo
    self%Ncgsdcs = Ncg
    allocate(self%ecg(Ncg))
    self%ecg(1:Ncg) = ecg(1:Ncg)

    ! interpolate to the common grid for SDCS: ecg  
    print*,' interpolate to the common grid for SDCS: ecg  '
    do i=1,data_in%Nsdcs
       etotd2 = (ein(i) - eioniz)/2.0 /data_in%eV
       jm = inumsdcs(i)
       en(1:jm,i) = en(1:jm,i) /etotd2

       call INTRPL(jm,en(1:jm,i),cs(1:jm,i),Ncg,ecg,vcg)

       cscg(1:Ncg,i) = vcg(1:Ncg)
    enddo

    ! for each point on the common grid interpolate over incident energies to the common incident energy grid
    print*,'for each point on the common grid interpolate over incident energies to the common incident energy grid'
    Ncgsdcs = Ncg
    do i=1,Nenfine
       if(enfine(i) .gt. ein(data_in%Nsdcs) ) exit
    enddo
    Ninsdcs =  i   ! it is OK to go a bit up... it is a fine grid...
    print*, 'Ninsdcs=', Ninsdcs, enfine(i)
! find the first point on the fine energy grid the first (lowest) SDCS is available
    do i=Ninsdcs,1,-1
       if(enfine(i) .lt. ein(1) ) exit
    enddo
    Ninsdcs_first = i ! it is OK to go a bit down... it is a fine grid...
    print*, 'Ninsdcs_first=', Ninsdcs_first, enfine(i), ein(1)
    
    self%Ninsdcs = Ninsdcs   ! this is the last point on the fine energy grid where SDCS are available 
    self%Ninsdcs_first = Ninsdcs_first ! this is the first point on the fine energy grid where SDCS is available 
    
    allocate(ensdcs(Ninsdcs))
    allocate(self%ein(Ninsdcs), self%CS(Ncgsdcs,Ninsdcs))
    ensdcs(1:Ninsdcs) = enfine(1:Ninsdcs)
    self%ein(1:Ninsdcs) = ensdcs(1:Ninsdcs)


    do i=1,Ncg
       jm = data_in%Nsdcs
       call INTRPL(jm,ein(1:jm),cscg(i,1:jm),Ninsdcs,ensdcs,venfine)
       self%CS(i,1:Ninsdcs) = venfine(1:Ninsdcs)
    enddo

    deallocate(ecg,cscg,vcg,ensdcs)

  end subroutine make_sdcs
  !
  !

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: get_sdcsatein
  !Purpose: returns a set of sdcs at the requested incident
  !         energy.
  !Important: eincident is assumed to be greater than the ionisation
  !           energy of the target state. It is the user's job
  !           to ensure this occurs.
  !Date last modified: 24/02/2022
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine get_sdcsatein(sdcsob,eincident,sdcsen)

    implicit none

    type(basis_sdcs), intent(in):: sdcsob
    real(dp), intent(in):: eincident
    type(sdcs):: sdcsen
    real(dp):: e1, e2, w1, w2
    integer:: i, j, j1, j2, Ncgsdcs, cutoffInd
    integer:: numPoints, ii, kk, jj
    real(dp):: deltaE, CS1, CS2, en1, en2, weight1, weight2
    real(dp):: eIon, cutoffEn, diffEn, deltaCS
    logical:: found, belowMax

    sdcsen%tcstype = sdcsob%tcstype
    sdcsen%Ncgsdcs = sdcsob%Ncgsdcs
    sdcsen%en_incident = eincident
    sdcsen%en_ionization = sdcsob%en_ionization
    sdcsen%en_total = eincident + sdcsob%en_ionization

    Ncgsdcs = sdcsen%Ncgsdcs

    !Determine size of uniform energy grid
    numPoints = 5000
    eIon = 15.96632  !Ionisation energy of ground state H2
    if (eincident .lt. eIon) then
       print*, "ERROR: request for sdcs below ionisation energy, stopping"
       stop
    end if
    cutoffEn = (eIncident-eIon)/2.0
    deltaE = (sdcsob%ecg(sdcsob%Ncgsdcs) - sdcsob%ecg(1))/real(numPoints)
    !print*, sdcsob%ecg(sdcsob%ncgsdcs), eIncident, deltaE
    
    !CutoffEn/deltaE > 2.0 required for algorithm to work
    if (int(cutoffEn/deltaE) .lt. 4) then
       !Handle case where eIncident approx equal to eIon
       deltaE = (eIncident-eIon)/4.0
       numPoints = int(sdcsob%ecg(sdcsob%Ncgsdcs)/deltaE)
    end if

    if ( associated(sdcsen%ecg) ) deallocate(sdcsen%ecg)
    if ( associated(sdcsen%cs) ) deallocate(sdcsen%cs)
        
    !allocate(sdcsen%ecg(Ncgsdcs),sdcsen%CS(Ncgsdcs))

    !Copy ejection energy grid
    !sdcsen%ecg(1:Ncgsdcs) = sdcsob%ecg(1:Ncgsdcs)

    !Create common energy grid of around 1000 points across 0-700eV.
    !Constant spacing required for sampling to work correctly
    allocate(sdcsen%ecg(numPoints),sdcsen%CS(numPoints))
    sdcsen%ecg(:)= 0.0
    sdcsen%CS(:) = 0.0
    sdcsen%ecg(1) = sdcsob%ecg(1)
    do ii=2,numPoints
       sdcsen%ecg(ii) = sdcsen%ecg(ii-1) + deltaE
    end do
    sdcsen%Ncgsdcs = numPoints


    !Fix below, use symmetrisation code from intpCS
    !Check CS are available adjacent to requested incident energy
    if(eincident-1e-6 .gt. sdcsob%ein(sdcsob%Ninsdcs)) then
       print*,'sdcs.f90: get_sdcsatein(): eincident is larger than the largest available incident energy for SDCS'
       print*, 'eincident, sdcsob%ein(sdcsob%Ninsdcs))', eincident, sdcsob%ein(sdcsob%Ninsdcs)
       print*, 'will return zero for SDCS'
       sdcsen%CS(1:Ncgsdcs) = 0d0
       return
    endif
    if(eincident+1e-6 .lt. sdcsob%ein(sdcsob%Ninsdcs_first)) then
       print*,'sdcs.f90: get_sdcsatein(): eincident is smaller than the smallest available incident energy for SDCS'
       print*, 'eincident, sdcsob%ein(sdcsob%Ninsdcs_first))', eincident, sdcsob%ein(sdcsob%Ninsdcs_first)
       print*, 'will return zero for SDCS'
       sdcsen%CS(1:Ncgsdcs) = 0d0
       return
    endif
    
    !Find sdcs arrays at incident energies adjacent to given energy
    j = 1
    do i=1,sdcsob%Ninsdcs
       if(eincident .lt. sdcsob%ein(i)) exit
       j = j + 1
    enddo
    if(j .eq. 1) then
       print*,'sdcs.f90: get_sdcsatein(): j=1'
       print*, ' could not find the required two energies'
       stop
    endif

    j1 = j -1
    j2 = j

    !print*, 'eincident = ', eincident
    !print*, 'j1, j2:', j1, j2
    !print*, 'sdcsob%ein(j2) = ', sdcsob%ein(j2)

   !Calculate interpolation weights
    e1 = sdcsob%ein(j1)
    e2 = sdcsob%ein(j2)
    w1 = (eincident - e1)/(e2 - e1) 
    w2 = (e2 - eincident)/(e2 - e1)  

    if (eIncident .ge. eIon) then
       !Faster intrpolation to energy grid with uniform spacing
       kk = 1
       belowMax = .true.
       ii = 2
       cutoffInd = int(cutoffEn/deltaE)+1
       do while(sdcsob%ecg(ii) .le. sdcsen%ecg(cutoffInd)) !ii .le. cutoffInd
          !Including condition for kk in this is something of a hack.
          !Needed to get program working.
          if (belowMax .and. (sdcsob%ecg(ii) .gt. sdcsen%ecg(kk))) then
             do while(belowMax .and. (sdcsob%ecg(ii) .gt. sdcsen%ecg(kk)))
                !Interpolate       
                en1 = sdcsob%ecg(ii-1)
                en2 = sdcsob%ecg(ii)
                weight1 = (sdcsen%ecg(kk)-en1)/(en2-en1)
                weight2 = (en2-sdcsen%ecg(kk))/(en2-en1)
                CS1 = weight1*sdcsob%CS(ii-1,j1) + weight2*sdcsob%CS(ii,j1) 
                CS2 = weight1*sdcsob%CS(ii-1,j2) + weight2*sdcsob%CS(ii,j2)  
                sdcsen%CS(kk) = w1*CS1 + w2*CS2 
 
                if (kk .lt. sdcsen%Ncgsdcs) then
                   kk = kk+1
                else
                   belowMax = .false.
                end if
             end do
          end if
          ii = ii+1
       end do

       !Symmetrise about point E/2 
       ii = cutoffInd !Point to symmetrise about
       do while (sdcsen%ecg(ii) .le. (eIncident-eIon))
          jj = cutoffInd-1
          found = .false.
          do while ((jj .ge. 2) .and. (found .eqv. .false.))
             if (sdcsen%ecg(ii)-cutoffEn .lt. cutoffEn-sdcsen%ecg(jj)) then
                found = .true.
             end if
             jj = jj-1
          end do
          !Linear interpolation over appropriate values below E/2
          deltaE = sdcsen%ecg(jj+1) - sdcsen%ecg(jj)
          deltaCS = sdcsen%CS(jj+1) - sdcsen%CS(jj)
          diffEn = sdcsen%ecg(ii)-cutoffEn - (cutoffEn-sdcsen%ecg(jj+1))
          !print*, "DIFFEN: ", diffEn, "DELTAE: ", deltaE
          ! '-' corrects sign of gradient, CS is decreasing
          sdcsen%CS(ii) = sdcsen%CS(jj+1) - (deltaCS/deltaE)*diffEn
          !print*, "CS: ", sdcsen%CS(ii)
          ii = ii+1
       end do
    else
       !Simpler interpolation, without worrying about perfect symmetry
       !Not actually sampled at these energies, so can use obect grid instead
       !Still gives correct result for TICS at low energies
       deallocate(sdcsen%CS,sdcsen%ecg)
       allocate(sdcsen%CS(sdcsob%Ncgsdcs),sdcsen%ecg(sdcsob%Ncgsdcs))
       sdcsen%ecg(:)=0.0
       sdcsen%CS(:) = 0.0
       sdcsen%Ncgsdcs= sdcsob%Ncgsdcs
       sdcsen%ecg(:) = sdcsob%ecg(:)   
       sdcsen%CS(1:Ncgsdcs) = w1*sdcsob%CS(1:Ncgsdcs,j1) + w2*sdcsob%CS(1:Ncgsdcs,j2)
    end if

    !Old interpolation when common grid was used in both objects
    !sdcsen%CS(1:Ncgsdcs) = w1*sdcsob%CS(1:Ncgsdcs,j1) + w2*sdcsob%CS(1:Ncgsdcs,j2)

    !Find index of last non-zero sdcs value on common grid
    i = 1
    found = .false.
    do while ((i .le. sdcsen%Ncgsdcs) .and. (.not. found))
       !Find first zero value within tolerance
       !Smallest CS  order of magnitude 10^-4
       if (sdcsen%CS(i) .lt. 0.000001) then
          found = .true.
       end if
       i = i +1
    end do
    i = i-1
    sdcsen%nFinal = i
  end subroutine get_sdcsatein
!
!




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: testSdcs
  !Purpose: integrates sdcs from E_0 to E_tot/2 and prints them to file 
  !         for comparison with total cross sections as a consistency 
  !         check. Values should be the same.
  !Date last modified: 13/07/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine testSdcs(sdcsBasis,minE)
      implicit none
      integer:: ii
      type(basis_sdcs)::sdcsBasis
      type(sdcs)::sdcsAtEn
      real(dp):: deltaE, change, minE
      real(dp):: eIncident, csInt, eIon, cutoffEn
     
      eIon = 15.96632  !Ionisation energy of ground state H2
      eIncident = minE
      deltaE = 2.0
      change = 0.0

      open(27, file="sdcsIntTest.txt")
      write(27,*)  "Energy(ev),            IntegratedCS(a.u)"
      do while (eIncident .lt. 500.05)
         if (eIncident .gt. eIon) then
            call get_sdcsatein(sdcsBasis, eIncident, sdcsAtEn) 
            cutoffEn = (eIncident-eIon)/2.0

           !Integrate over SDCS from E_min to E/2
            csInt = 0.0
            ii=1
            do while (sdcsAtEn%ecg(ii) .le. cutoffEn)
               change = 0.5*(sdcsAtEn%CS(ii+1)+sdcsAtEn%CS(ii))*(sdcsAtEn%ecg(ii+1)-sdcsAtEn%ecg(ii))
               csInt = csInt + change        
               ii=ii+1 
            end do 
         else if( eIncident .lt. eIon) then
            !SDCS are all zero below the ionisation energy.
            csInt = 0.0_dp
         end if

         !Print to file at current energy 
         write(27,*)  eIncident, csInt
         eIncident = eIncident + deltaE
         call destructSdcsAtEIn(sdcsAtEn)
      end do
      close(27)
  end subroutine testSdcs


  subroutine print_sdcs(self)

    implicit none
    
    type(sdcs), intent(in):: self
    integer:: i

    print*,'# SDCS'
    print*,'# tcstype =', self%tcstype
    print*,'# en_incident =', self%en_incident
    print*,'# en_total =', self%en_total
    print*,'# en_ionization =', self%en_ionization
    print*,'#  i      ecg(i)     ecg(i)*en_total (eV)   self%CS(i) (a.u.)'
    do i=1,self%Ncgsdcs
       print'(i5,1P,2E15.5,5X,E15.5)', i, self%ecg(i), self%ecg(i)*self%en_total, self%CS(i)
    enddo
 
  end subroutine print_sdcs



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: printSdcsSet
   !Purpose: prints SDCS to a file for plotting for a particular set
   !         of incident energies.
   !Date last modified: 01/12/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine printSdcsSet(sdcsBasis)
       implicit none
       type(basis_sdcs)::sdcsBasis
       type(sdcs)::sdcsAtEn
       integer::ii
       real(dp),dimension(6)::enArray
       real(dp)::eIn 

       enArray = (/17.0,20.0,50.0,100.0,300.0,500.0/)
       do ii =1,6
          eIn = enArray(ii)
          call get_sdcsatein(sdcsBasis,eIn,sdcsAtEn)
          call printSdcsAtEn(sdcsAtEn) 
          call destructSdcsAtEIn(sdcsAtEn) 
       end do

   end subroutine printSdcsSet



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: printSdcsAtEn
  !Purpose: prints cross sections in imported sdcs type to a file
  !         for plotting. Used to check validity of interpolation.
  !Date last modified: 16/07/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine printSdcsAtEn(self)
      implicit none
      type(sdcs)::self
      integer:: ii
      character (len=3)::enString
 
      write(enString,'(I0)') int(self%en_incident)
      open(70, file="sdcs"//trim(enString)//".txt")
      write(70,*) "Energy(eV)          CS(a.u)"
      do ii = 1, self%Ncgsdcs
         write(70,*) self%ecg(ii), self%CS(ii)
      end do
      close(70)    

  end subroutine printSdcsAtEn




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: printSDCSAll
   !Purpose: prints SDCS over entire energy range (eMin-500 eV) at all
   !         incident energies to a file. Output used to create a
   !         3D plot of SDCS cross sections.
   !Date last modified: 30/11/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine printSDCSAll(sdcsbasis,eMin)
       implicit none
       type(basis_sdcs)::sdcsbasis
       integer::kk
       type(sdcs)::sdcsatein
       real(dp)::eIncident, eMin, eIon

       eIon = 15.96632_dp  !Ionisation energy of ground state H2
 
       open(12, file="sdcsSurface.txt") 
       write(12,*) "EIN          EOUT           SDCS"
       eIncident = eMin
       do while (eIncident .lt. 500.0_dp) 
          !SDCS are only non-zero above ionisation energy of the target
          if (eIncident .gt. eIon) then
             call get_sdcsatein(sdcsbasis,eIncident,sdcsatein)       
             kk = 1   
             do while (sdcsatein%ecg(kk) .lt. 500.0_dp)
                write(12,*) eIncident, sdcsatein%ecg(kk), sdcsatein%CS(kk) 
                kk = kk + int(real(sdcsatein%Ncgsdcs)/100.0_dp)
             end do
             call destructSdcsAtEIn(sdcsatein)
          end if
          eIncident = eIncident + 5.0_dp
       end do
       close(12)

   end subroutine printSDCSAll



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: printIntpSdcs
    !Purpose: calls a function in a loop to print interpolated  SDCS to a file
    !         as a function of incident energy, over a range of ejection
    !         energies.
    !Date last modified: 01/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printIntpSdcs(sdcsBasis)
        implicit none
        type(basis_sdcs)::sdcsBasis
        real(dp), dimension(6)::outArray
        real(dp)::eOut
        integer::kk

        outArray = (/0.0,1.0,20.0,50.0,100.0,300.0/)

        do kk =1,6
           eOut = outArray(kk)
           call printSdcsAtEOut(sdcsBasis,eOut) 
        end do 
     end subroutine printIntpSdcs







    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: printSdcsAtEOut
    !Purpose:prints the SDCS for given ejection energy at all incident energies
    !Date last modified: 01/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printSdcsAtEOut(sdcsBasis,eOut)
        implicit none
        type(basis_sdcs)::sdcsBasis
        integer::ii, kk
        real(dp)::eIncident, eOut
        type(sdcs)::sdcsaten 
        character (len=3)::enString

        write(enString,'(I0)') int(eOut)
        open(12, file="sdcsInit"//trim(enString)//"eV.txt")
        ii=1
        write(12,*) "EIncident            SDCS"
        eIncident = 0.0
        do while (eIncident .lt. 500.0)
           eIncident = sdcsBasis%eIn(ii)
           call get_sdcsatein(sdcsBasis,eIncident,sdcsaten)
           kk = 1
           !Find index corresponding to given ejection energy
           do while (sdcsaten%ecg(kk) .lt. eOut)
              kk = kk+1
           end do

           write(12,*) eIncident, sdcsaten%CS(kk)
           call destructSdcsAtEIn(sdcsaten)
           ii = ii+1
        end do 
        close(12)

    end subroutine printSdcsAtEOut




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: printInputSdcs
    !Purpose: calls a function in a loop to print CCC SDCS to a file
    !         as a function of incident energy, over a range of ejection
    !         energies.
    !Date last modified: 01/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printInputSdcs(sdcsBasis)
        implicit none
        type(basis_sdcs)::sdcsBasis
        real(dp), dimension(6)::outArray
        real(dp)::eOut
        integer::kk

        outArray = (/0.0,1.0,20.0,50.0,100.0,300.0/)

        do kk =1,6
           eOut = outArray(kk)
           call printInputSdcsAtEOut(sdcsBasis,eOut) 
        end do 
     end subroutine printInputSdcs


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: printInputSdcsAtEOut
    !Purpose: prints the input SDCS generated by the CCC code for a 
    !         given ejection energy at all incident energies.
    !Note: correct use requires calling while sdcsBasis object still 
    !      stores cross sections read in from a file.
    !Date last modified: 01/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printInputSdcsAtEOut(sdcsBasis,eOut)
        implicit none
        type(basis_sdcs)::sdcsBasis
        integer::ii, kk
        real(dp)::eIncident, eOut, eIon
        character (len=3)::enString


        eIon = 15.96632 !Ionisation energy of ground state H2		

        write(enString,'(I0)') int(eOut)
        open(70, file="cccSdcsAtEOut"//trim(enString)//"eV.txt")
        ii=1
        write(70,*) "EIncident            SDCS"
        eIncident = 0.0
        do while (eIncident .lt. 500.0)
           eIncident = sdcsBasis%eIn(ii)
           if (eIncident-eIon .gt. eOut) then
              kk = 1
              !Find index corresponding to given ejection energy
              do while (sdcsBasis%readEOut(ii,kk) .lt. eOut)
                 kk = kk+1
              end do

              write(70,*) eIncident, sdcsBasis%CS(ii,kk)
           else
              write(70,*) eIncident, 0.0
           end if 
           ii = ii+1
        end do 
        close(70)

    end subroutine printInputSdcsAtEOut



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: readSdcs
  !Purpose: reads single differential cross sections from 
  !         csv files listed in input file.
  !Date last modified:29/04/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readSdcs(data_input, sdcsBasis)
      use input_data

      implicit none
      type(basis_sdcs)::sdcsBasis
      type(input)::data_input
      integer:: largeN, ii, jj, numSDCS
      real(dp), allocatable, dimension(:)::tempArray
      integer:: nBelowCount

      !Reallocate arrays to fit input data
      largeN = 3000 
      if (associated(sdcsBasis%ein)) then
         deallocate(sdcsBasis%ein)
      end if
      allocate(sdcsBasis%ein(data_input%Nsdcs))

      if(associated(sdcsBasis%ecg)) then
         deallocate(sdcsBasis%ecg)
      end if
      allocate(sdcsbasis%ecg(largeN))

      if (associated(sdcsBasis%CS)) then
         deallocate(sdcsBasis%CS)
      end if
      allocate(sdcsBasis%CS(data_input%Nsdcs,largeN))

      if (associated(sdcsBasis%readEOut)) then
         deallocate(sdcsBasis%readEOut)
      end if

      !Find number of SDCS, should be same for all files.
      !Unfortunately, fortran IO forces use of a goto,
      !place marker as close as possible to 'read' call
      open(30, file=data_input%filename_sdcs(1),action='read')
      read(30,*)
      read(30,*)
      do ii=1, largeN
         read(30,*,ERR=70,END=70)
      end do
70 continue
      close(30)
      numSDCS=ii
      allocate(sdcsBasis%readEOut(data_input%Nsdcs, numSDCS))
      allocate(tempArray(numSDCS)) !Used for file IO only

      !Initialise arrays to zero
      sdcsBasis%readEOut(:,:) = 0.0
      sdcsBasis%CS(:,:) = 0.0

      jj=0
      !Read all new sdcs files 
      do ii = 1, data_input%Nsdcs
         tempArray(:) = 0.0
         nBelowCount = 0
         open(30, file=data_input%filename_sdcs(ii),action='read')
         read(30,*)
         read(30,*) sdcsBasis%ein(ii)
         jj=0
         do while (jj .le. largeN)
            jj = jj +1
            read(30,*, ERR=50, END=50) sdcsBasis%readEOut(ii,jj),  tempArray(jj)
            !File units are in rydbergs
            sdcsBasis%readEOut(ii,jj) = sdcsBasis%readEOut(ii,jj)*13.6 
            sdcsBasis%CS(ii,jj) = tempArray(jj)
            sdcsBasis%CS(ii,jj) = sdcsBasis%CS(ii,jj)/13.6  
         end do
50 continue
         close(30)
         !sdcsBasis%nBelowInFile(ii) = nBelowCount
      end do
      deallocate(tempArray)
      !Should be same number of CS in all files
      sdcsBasis%Nsdcsatein = jj-1 
      sdcsBasis%Nein = data_input%Nsdcs
  end subroutine readSdcs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: intpSdcsModFine
  !Purpose: interpolates sdcs cross sections to a common 
  !         energy grid for each incident energy. Overall
  !         a 2D cubic interpolation. Uses a fine grid.
  !Note: assumes number of grid points is a multiple of 1000
  !Date last modified: 09/02/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Note on subroutine: Interpolates over *patches* of the fine energy grid, 
      !  then combines them. This allows for the necessary fineness at 
      !  low energies, while still using INTRPL
  subroutine intpSdcsModFine(data_input,sdcsBasis, Nenfine, enfine)
      use input_data
      implicit none
      type(basis_sdcs)::sdcsBasis
      type(input)::data_input
      integer::Nenfine, ii, jj, numEnIn, kk, counter,ll
      integer::enJustBelow, minEInIndex
      real(dp), dimension(:)::enfine
      real(dp), allocatable, dimension(:,:)::newBasis
      real(dp), allocatable, dimension(:,:)::tempBasis, transposeBasis
      real(dp), allocatable, dimension(:)::tempCS
      real(dp), allocatable, dimension(:)::tempCSAtEIn, tempEOut, tempIntpCs
      real(dp), allocatable, dimension(:)::enGridBelow
      real(dp):: currEn, eIon, enSpacing, cutoffEn, enLast
      logical::found
      real(dp)::deltaE, deltaCS, diffEn, minEIn 
      real(dp), dimension(3)::tempCSBelow,tempEIn
      !Used to find size of array for ejection energies below enfine(1)
      integer:: numSDCS, cutoffInd, numThousands

      numEnIn = sdcsBasis%Nein

      !Create common energy grid and initialise new basis
      if (associated(sdcsBasis%ecg)) then
         deallocate(sdcsBasis%ecg)
      end if 

      !Finer energy grid for SDCS
      enSpacing = 700.0/10000.0 !spacing of uniform SDCS grid (eV)
      sdcsBasis%Ncgsdcs = int(700.0/enSpacing)
      allocate(sdcsBasis%ecg(sdcsBasis%Ncgsdcs))
      !Use new, uniform energy grid  
      sdcsBasis%ecg(1) = sdcsBasis%readEOut(1,1)
      do ii = 2, sdcsBasis%Ncgsdcs 
         sdcsBasis%ecg(ii) = sdcsBasis%ecg(ii-1) + enSpacing
      end do
 
      !Create new basis and intialise to zero for safety
      allocate(newBasis(Nenfine,sdcsBasis%Ncgsdcs))
      newBasis(:,:) = 0.0

      !Define temporary array to store interpolated SDCS values
      allocate(tempBasis(data_input%Nsdcs,sdcsBasis%Ncgsdcs))
      tempBasis(:,:) = 0.0       
 
      !Each row of the 2D array is SDCS at a particular
      !incident electron energy.
      eIon = 15.96632  !Ionisation energy of ground state H2
      minEIn = sdcsBasis%ein(1)
      sdcsBasis%einFirst = 1
      sdcsBasis%einLast = 1
      enLast = sdcsBasis%ein(data_input%Nsdcs)
      !At each available incident energy, interpolate over ejection energy grid.
      do ii = 1, data_input%Nsdcs
         currEn = sdcsBasis%ein(ii)
         numSDCS = size(sdcsBasis%readEOut(ii,:))

         !Find the index of the point on the new energy grid where the 
         !energy is equal to (eIncident - eIon), interpolate only to here.
         counter=1
         do while (sdcsBasis%ecg(counter)  .lt. currEn-eIon)
            counter= counter+1
         end do
 
         !Interpolate over patches of the grid until complete. Done because the INTRPL
         !function can only handle 1000 data points at a time.  
         allocate(tempEOut(numSDCS-1),tempCSAtEIn(numSDCS-1))
         tempEOut(:) = sdcsBasis%readEOut(ii,1:numSDCS-1)
         tempCSAtEIn(:) = sdcsBasis%cs(ii,1:numSDCS-1)
         numThousands = int(real(sdcsBasis%Ncgsdcs)/1000.0)
         do ll = 1, numThousands
            allocate(tempIntpCs(1000))
            tempIntpCs(:) = 0.0
            call INTRPL(numSDCS-1,tempEOut,tempCSAtEin,1000,sdcsBasis%ecg((ll-1)*1000+1:ll*1000),tempIntpCs)
            tempBasis(ii,(ll-1)*1000+1:ll*1000) = tempIntpCs(:)
            deallocate(tempIntpCs)
         end do
         !INTRPL subroutine extrapolates cross sections to energies above eIncident-eIon
         !set CS at these energies to zero
         tempBasis(ii,counter:sdcsBasis%Ncgsdcs) = 0.0
         deallocate(tempEOut,tempCSAtEIn)

         !Symmetrise about E/2.
         !Note: Need to check energy of each grid point above E/2 
         !      and perform linear interpolation over appropriate
         !      values below E/2 to get correct symmetry.
         cutoffEn = (currEn-eIon)/2.0
         cutoffInd=1
         found = .false.
         do while (.not. found)
            if(sdcsBasis%ecg(cutoffInd) .gt. cutoffEn) then
               found = .true.
            end if
            cutoffInd=cutoffInd+1
         end do
         cutoffInd = cutoffInd-1

         if (cutoffEn .gt. sdcsBasis%ecg(1)) then
            kk = cutoffInd
            do while (sdcsBasis%ecg(kk) .le. (currEn-eIon))
               jj = cutoffInd-1
               found = .false.
               do while ((jj .ge. 2) .and. (found .eqv. .false.))
                  if (sdcsBasis%ecg(kk)-cutoffEn .lt. cutoffEn-sdcsBasis%ecg(jj)) then
                     found = .true.
                  end if
                  jj = jj-1
               end do
               if (found .eqv. .true.) then
                  !Linear interpolation over appropriate values below E/2
                  deltaE = sdcsBasis%ecg(jj+1) - sdcsBasis%ecg(jj)
                  deltaCS = tempBasis(ii,jj+1) - tempBasis(ii,jj)
                  diffEn = sdcsBasis%ecg(kk)-cutoffEn - (cutoffEn-sdcsBasis%ecg(jj+1))
                  ! '-' corrects sign of gradient, CS is decreasing
                  tempBasis(ii,kk) = tempBasis(ii,jj+1) - (deltaCS/deltaE)*diffEn
               end if
               !Else, don't symmetrise, as selected energy < eIon when this occurs, so can 
               !safely ignore particle, as it will be removed from array in next iteration.
               kk = kk+1
            end do
         end if
      end do

      !Repeat process for new grid of incident energies.
      do ii = 1, sdcsBasis%Ncgsdcs
         !Create array of cross sections at given ejection energy
         allocate(tempCS(data_input%Nsdcs))
         tempCS(:) = tempBasis(:,ii)
         !Interpolate over cross sections
         call INTRPL(data_input%Nsdcs,sdcsBasis%ein,tempCS,Nenfine,enfine,newBasis(:,ii))
         deallocate(tempCS)
      end do
      
      !Need to do energies below lowest file incident energy separately,
      !as the INTPL subroutine attempts to extrapolate. 
      enJustBelow = 1
      minEInIndex = 1
      do ii =1,Nenfine
         if (enfine(ii) .lt. eIon) then 
            !SDCS are zero below ionisation energy by definition
            newBasis(ii,:) = 0.0
            enJustBelow = enJustBelow + 1
         end if
         if (enfine(ii) .lt. sdcsBasis%ein(1)) then
            minEInIndex= minEInIndex+1
         end if
      end do
      enJustBelow = enJustBelow-1

      !Procedure: need to interpolate from zero for energies below minimum incident
      !           in input files.
      !Create array of cs at energy minEIn
      !Looping over ejection energies, create array(len=3) of cs at a given ejection energy
      !As above, interpolate this array over incident energies in the range eJustBelow-minEIn
      allocate(enGridBelow(minEInIndex-enJustBelow+2))
      enGridBelow(:) = enfine(enJustBelow-1:minEInIndex) !start at two points below eIon, CS=0 at both
      !For each ejection energy, interpolate over incident energy grid
      do ii=1, sdcsBasis%Ncgsdcs
         !Interpolate from zero, starting at two points below eIon on fine energy grid 
         tempCSBelow = (/0.0_dp,0.0_dp,newBasis(minEInIndex,ii)/)    
         tempEIn = (/enfine(enJustBelow-1),enfine(enJustBelow),enfine(minEInIndex)/)  ! sdcsBasis%ein(1)/) !Is this third value right?
         call INTRPL(3,tempEIn,tempCSBelow,minEInIndex-enJustBelow+2,enGridBelow,newBasis(enJustBelow-1:minEInIndex,ii)) 
         !print*, "eout, cs: ", sdcsBasis%ecg(ii), newBasis(minEInIndex,ii), newBasis(minEInIndex+1,ii), newBasis(minEInIndex+2,ii), newBasis(minEInIndex+5,ii), newBasis(minEInIndex+20,ii), newBasis(minEInIndex+50,ii)
         !print*, "eOut: ", sdcsBasis%ecg(ii)
         !print*, newBasis(enJustBelow-1:minEInIndex, ii)
      end do 
      deallocate(tempBasis, enGridBelow)
 
      !Final basis is the transpose of that calculated, in order to match
      !indexing scheme used in previously written functions 
      allocate(transposeBasis(sdcsBasis%Ncgsdcs,Nenfine),sdcsBasis%CS(sdcsBasis%Ncgsdcs,Nenfine)) 
      transposeBasis = transpose(newBasis)
      sdcsBasis%CS(:,:) = transposeBasis(:,:)  
      deallocate(transposeBasis,newBasis)

      !Replace ein with grid of incident energies
      deallocate(sdcsBasis%ein)
      allocate(sdcsBasis%ein(Nenfine))
      sdcsBasis%ein(:) = enfine(:)       
      sdcsBasis%Nein = Nenfine

      !Find index of point on incident energry grid above which
      !cross sections are no longer available
      do ii = 1, Nenfine
         if (sdcsBasis%ein(ii) .lt. enlast) then
            sdcsBasis%einLast = sdcsBasis%einLast +1
         end if
      end do
      sdcsBasis%einLast = sdcsBasis%einLast -1

      !Finish interpolation
      sdcsBasis%Nein = sdcsBasis%einLast - sdcsBasis%einFirst +1
      sdcsBasis%Ninsdcs = sdcsBasis%einLast
      sdcsBasis%Ninsdcs_first = sdcsBasis%einFirst 
      sdcsBasis%Nsdcsatein = sdcsBasis%Ncgsdcs
  end subroutine intpSdcsModFine





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !subroutine: scaleSdcs
    !Purpose: scales SDCS at a particular incident energy so 
    !         that the calculated TICS agrees with that
    !         in the totalCS files. 
    !Needed for: needed because the input SDCS come from an 
    !            interpolation procedure which involves scaling
    !            values. Sometimes this is not perfect and
    !            may need correcting.
    !Warning: use sparingly, as this will obscure major issues 
    !         with the SDCS input files if there are any.
    !Date last modified: 10/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine scaleSdcs(sdcsBasis,stateBasis)
        use state_class
        use totalcs_module
        implicit none
        integer:: ii
        real(dp):: ratio, tics, eIon
        real(dp), allocatable,dimension(:)::intSdcs
        type(totalcs)::tcs
        type(basis_sdcs)::sdcsBasis
        type(basis_state)::stateBasis

        allocate(intSdcs(sdcsBasis%Nein))
        call getIntSdcs(sdcsBasis,intSdcs)

        !Write to file for debugging
        open(70, file="intTest1.txt")
        eIon = 15.96632  !Ionisation energy of ground state H2
        do ii=1, sdcsBasis%Nein
           if (sdcsBasis%ein(ii) .gt. eIon) then
              call get_csatein(stateBasis,sdcsBasis%ein(ii),tcs)
              call getTics(tcs,tics) !Get accurate total ionisation CS
              !Find the integrated SDCS and totalCS at a given energy
              if (intSDCS(ii) .gt. 0.0) then
                 ratio = tics/intSdcs(ii)
                 write(70,*) sdcsBasis%ein(ii), tics, intSdcs(ii)
              else
                 !Total ionization CS is zero
                 ratio = 1.0
              end if

              !Scale SDCS at EIN by ratio totalCS(EIN)/INTSDCS(EIN)
              sdcsBasis%cs(:,ii) = sdcsBasis%cs(:,ii)*ratio
              call delete_totalcs(tcs)
           else
              !Do nothing cross sections set to zero in interpolation
           end if
        end do 
        close(70)


        deallocate(intSdcs)
    end subroutine scaleSdcs




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: distortSdcs
    !Purpose: distorts SDCS over energy range in such a way that 
    !         total ionisation cross section is conserved. Used to
    !         check sensitivity of results to SDCS values.
    !Note: done by multiplying by a 'distortion function' defined 
    !      over same energy grid as SDCS.
    !Date last modified: 11/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine distortSdcs(sdcsBasis,stateBasis)
        use totalcs_module
        use state_class
        implicit none
        type(basis_sdcs)::sdcsBasis
        type(basis_state)::stateBasis
        type(totalcs)::tcs
        real(dp), allocatable, dimension(:)::distortFunc
        real(dp)::eIon, tics
        real(dp), allocatable, dimension(:)::intSdcs
        integer::ii,kk

        eIon = 15.96632 !Ionisation energy of ground state H2 (eV)	
        !Perform distortion for all incident energies capable of causing
        !ionisation
        do kk = 1, sdcsBasis%Nein
              if (sdcsBasis%ein(kk) .gt. eIon) then
                 allocate(distortFunc(sdcsBasis%Ncgsdcs))
                 do ii = 1, sdcsBasis%Ncgsdcs
                    !Can choose any function, doesn't matter which, try several
                    distortFunc(ii) = sin(sdcsBasis%ecg(ii)/sdcsBasis%ecg(sdcsBasis%Ncgsdcs))
                    !distortFunc(ii) = exp((-5)*sdcsBasis%ecg(ii)/sdcsBasis%ecg(sdcsBasis%Ncgsdcs))
                    !distortFunc(ii) = log(1.0+sdcsBasis%ecg(ii)/sdcsBasis%ecg(sdcsBasis%Ncgsdcs))
                 end do
                 do ii=1, sdcsBasis%Ncgsdcs
                    sdcsBasis%cs(ii,kk) = sdcsBasis%cs(ii,kk)*distortFunc(ii)
                 end do
                 deallocate(distortFunc)
              end if
        end do
 
        allocate(intSdcs(sdcsBasis%Nein))
        call getIntSdcs(sdcsBasis,intSdcs)
        do kk=1,sdcsBasis%Nein
           !Get TICS from distorted SDCS and totalCS, scale by their ratio to 
           !preserve TICS.
           if (sdcsBasis%ein(kk) .gt. eIon) then
              call get_csatein(stateBasis,sdcsBasis%ein(kk),tcs)
              tics = 0.0
              call getTics(tcs,tics) 
              sdcsBasis%cs(:,kk) = sdcsBasis%cs(:,kk)*(tics/intSdcs(kk))         
              call delete_totalcs(tcs) 
           end if 
        end do
        deallocate(intSdcs)
    end subroutine distortSdcs



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: getIntSdcs
    !Purpose: returns an array of integrated sdcs, equal to the 
    !         the total ionisation cross section over the incident 
    !         energy grid.
    !Date last modified: 10/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine getIntSdcs(sdcsBasis, intSdcs)
        implicit none
        type(basis_sdcs)::sdcsBasis
        real(dp), allocatable, dimension(:)::intSdcs
        integer:: ii, jj
        real(dp)::intVal, deltaE, cutoffEn, eIon


        eIon = 15.96632 !Ionisation energy of ground state H2 (eV)	

        do ii = 1, sdcsBasis%Nein
           intVal = 0.0
           cutoffEn = (sdcsBasis%ein(ii)-eIon)/2.0
           do jj = 1, sdcsBasis%Ncgsdcs-1
              if ((sdcsBasis%ecg(jj) .lt. cutoffEn) .or. (sdcsBasis%ein(ii) .lt. 20.0)) then
                 deltaE = sdcsBasis%ecg(jj+1) - sdcsBasis%ecg(jj)
                 intVal = intVal + 0.5*(sdcsBasis%cs(jj,ii)+sdcsBasis%cs(jj+1,ii))*deltaE
              end if
           end do
           intSdcs(ii) = intVal
        end do 
    end subroutine getIntSdcs
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Submodule: normalisesdcs
  !Purpose: normalises single differential cross section set
  !         and constructs probability distribution
  !Date last modified: 06/05/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine normalisesdcs(sdcsAtEin, enProb,numCS)
     implicit none
     type(sdcs)::sdcsAtEin
     integer:: numCS, ii
     real(dp), dimension(numCS)::enProb
     real(dp)::total
   
     !numCS = sdcsAtEIn%Ncgsdcs 
     !allocate(enProb(numCS))
     total = SUM(sdcsAtEIn%CS)

     if (total .gt. 0.0) then
        enProb(1) = sdcsAtEin%CS(1)/total
        do ii = 2,numCS
           enProb(ii) = enProb(ii-1) + (sdcsAtEIn%CS(ii)/total)
        end do
     else
        enProb(:) = 0.0
        enProb(1) = 2.01
     end if
  end subroutine normalisesdcs






  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: destructSdcsBasis
  !Purpose: deallocates all memory used in storing SDCS
  !Date last modified: 30/04/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine destructSdcsBasis(sdcsBasis)
     implicit none
     type(basis_sdcs)::sdcsBasis
   
     if (associated(sdcsBasis%CS)) then
        deallocate(sdcsBasis%CS)
     end if
     if (associated(sdcsBasis%ein)) then
        deallocate(sdcsBasis%ein)
     end if 
     if (associated(sdcsBasis%ecg)) then
        deallocate(sdcsBasis%ecg)
     end if
     if (associated(sdcsBasis%readEOut)) then
        deallocate(sdcsBasis%readEOut)
     end if
     if (associated(sdcsBasis%nBelowInFile)) then
        deallocate(sdcsBasis%nBelowInFile)
     end if
     if (associated(sdcsBasis%CSDist)) then
        deallocate(sdcsBasis%CSDist)
     end if
     if (associated(sdcsBasis%ecgDist)) then
        deallocate(sdcsBasis%ecgDist)
     end if
  end subroutine destructSdcsBasis


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: destructSdcsAtEIn
   !Purpose: destructs an sdcs object at a specific energy
   !Date last modified: 28/05/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine destructSdcsAtEin(sdcsAtEIn)
       implicit none
       type(sdcs)::sdcsAtEIn
     
       if (associated(sdcsAtEIn%CS)) then
          deallocate(sdcsAtEIn%CS)
       end if
       if (associated(sdcsAtEIn%ecg)) then
          deallocate(sdcsAtEIn%ecg)
       end if

   end subroutine destructSdcsAtEIn



!Code above used for numerical SDCS calculated using MCCC method.
!Code below is for approximate SDCS function of opal, peterson and beaty (1971)
!This function is also used by Dalgarno et al (1998). 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: getSdcsDistAtEIn
   !Purpose: calculates SDCS at given incident energy using analytical
   !         distribution.
   !Note: as with version of this function above, creates a new energy 
   !      grid each time to avoid issues at incident energies close to 
   !      to the ionisation potential.
   !Important: eincident is assumed to be greater than the ionisation 
   !           potential of the target state. It is the user's job to ensure 
   !           that this occurs.
   !Date last modified: 24/09/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine getSdcsDistAtEIn(sdcsob,eincident,sdcsen,tcsob)
       use totalcs_module
       use state_class
       implicit none

       type(basis_sdcs), intent(in):: sdcsob
       real(dp), intent(in):: eincident
       type(sdcs):: sdcsen
       type(totalcs)::tcsob
       integer:: Ncgsdcs, cutoffInd
       integer:: numPoints, ii, jj
       real(dp):: deltaE
       real(dp):: eIon, cutoffEn, diffEn, deltaCS, eOut, eBar
       logical:: found

		
       sdcsen%tcstype = sdcsob%tcstype
       sdcsen%Ncgsdcs = sdcsob%Ncgsdcs
       sdcsen%en_incident = eincident
       sdcsen%en_ionization = sdcsob%en_ionization
       sdcsen%en_total = eincident + sdcsob%en_ionization

       Ncgsdcs = sdcsen%Ncgsdcs

       !Determine size of uniform energy grid
       numPoints = 1000
       eIon = 15.96632  !Ionisation energy of ground state H2
       if (eincident .lt. eIon) then
          print*, "ERROR: request for sdcs below ionisation energy, stopping"
          stop
       end if

       cutoffEn = (eIncident-eIon)/2.0
       deltaE = (sdcsob%ecgDist(sdcsob%NcgsdcsDist) - sdcsob%ecgDist(1))/real(numPoints)
       !CutoffEn/deltaE > 2.0 required for algorithm to work
       if (int(cutoffEn/deltaE) .lt. 4) then
          !Handle case where eIncident approx equal to eIon
          deltaE = (eIncident-eIon)/4.0
          numPoints = int((sdcsob%ecgDist(sdcsob%NcgsdcsDist)-sdcsob%ecgDist(1))/deltaE)
       end if

       if ( associated(sdcsen%ecg) ) deallocate(sdcsen%ecg)
       if ( associated(sdcsen%cs) ) deallocate(sdcsen%cs)
         
       !allocate(sdcsen%ecg(Ncgsdcs),sdcsen%CS(Ncgsdcs))

       !Copy ejection energy grid
       !sdcsen%ecg(1:Ncgsdcs) = sdcsob%ecg(1:Ncgsdcs)

       !Create common energy grid of around 1000 points across 0-700eV.
       !Constant spacing required for sampling to work correctly
       allocate(sdcsen%ecg(numPoints),sdcsen%CS(numPoints))
       sdcsen%ecg(:)= 0.0
       sdcsen%CS(:) = 0.0
       sdcsen%ecg(1) = sdcsob%ecgDist(1)
       do ii=2,numPoints
          sdcsen%ecg(ii) = sdcsen%ecg(ii-1) + deltaE
       end do
       sdcsen%Ncgsdcs = numPoints

       !Use analytical formula instead of interpolation
       !Energy parameter for H2, will need updating in general case.
       eBar = 8.3 !eV 
       do ii = 1, sdcsen%Ncgsdcs
          eOut = sdcsen%ecg(ii)
          sdcsen%CS(ii) = 1/(1+(eOut/eBar)**2.1) !Formula
       end do

       !Normalise calculated distribution
       call normaliseSDCSDist(tcsob,sdcsen,eIncident)

       !Symmetrise analytical SDCS about midpoint
       cutoffInd = int(cutoffEn/deltaE)+1
       ii = cutoffInd !Point to symmetrise about
       do while (sdcsen%ecg(ii) .le. (eIncident-eIon))
          jj = cutoffInd-1
          found = .false.
          do while ((jj .ge. 2) .and. (found .eqv. .false.))
             if (sdcsen%ecg(ii)-cutoffEn .lt. cutoffEn-sdcsen%ecg(jj)) then
                found = .true.
             end if
             jj = jj-1
          end do
          !Linear interpolation over appropriate values below E/2
          deltaE = sdcsen%ecg(jj+1) - sdcsen%ecg(jj)
          deltaCS = sdcsen%CS(jj+1) - sdcsen%CS(jj)
          diffEn = sdcsen%ecg(ii)-cutoffEn - (cutoffEn-sdcsen%ecg(jj+1))
          !print*, "DIFFEN: ", diffEn, "DELTAE: ", deltaE
          ! '-' corrects sign of gradient, CS is decreasing
          sdcsen%CS(ii) = sdcsen%CS(jj+1) - (deltaCS/deltaE)*diffEn
          !print*, "CS: ", sdcsen%CS(ii)
          ii = ii+1
       end do

       !!NOTE: sdcsen is returned by the function, memory deallocation 
       !       must occur in wherever this function is called.

   end subroutine getSdcsDistAtEIn






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: constructSdcsDist
    !Purpose: creates array of SDCS defined over 2D energy 
    !         grid using an analytical distribution.
    !Date last modified: 22/10/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine constructSdcsDist(sdcsob)
        use totalcs_module
        use state_class
        implicit none
        type(basis_sdcs)::sdcsob
        integer::ii, numPoints
        real(dp):: eBar !Energy parameter from fit to experiment
        real(dp):: eOut, deltaE

        !Determine size of uniform energy grid
        numPoints = 1000
        deltaE = (sdcsob%ecg(sdcsob%Ncgsdcs) - sdcsob%ecg(1))/real(numPoints)

        if (associated(sdcsob%CSDist)) then
           deallocate(sdcsob%CSDist)
        end if
        if (associated(sdcsob%ecgDist)) then
           deallocate(sdcsob%ecgDist)
        end if
        allocate(sdcsob%CSDist(numPoints,sdcsob%Nein),sdcsob%ecgDist(numPoints))
        sdcsob%NcgsdcsDist = numPoints  
        
        sdcsob%ecgDist(1) = sdcsob%ecg(1)
        do ii = 2, numPoints
           sdcsob%ecgDist(ii) = sdcsob%ecgDist(ii-1) + deltaE 
        end do
 
        !Energy parameter for H2, will need updating in general case.
        eBar = 8.3 !eV 
        do ii = 1, sdcsob%NcgsdcsDist
           eOut = sdcsob%ecgDist(ii)
           sdcsob%CSDist(ii,:) = 1/(1+(eOut/ebar)**2.1) !Formula
        end do

    end subroutine constructSdcsDist






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: normaliseSDCSDist
    !Purpose: calculates normalisation constants for the analytical SDCS
    !         distribution at the given incident energy
    !Date last modified; 22/09/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine normaliseSDCSDist(tcsob,sdcsen,eIncident)
        use totalcs_module
        use state_class
        implicit none
        type(totalcs)::tcsob !Contains total ionisation CS
        type(sdcs)::sdcsen
        integer:: jj 
        real(dp)::intDist !Integral over analytical distribution
        real(dp)::eIon, tics, deltaE, eIncident
        
        eIon = 15.96632  !H2 ionisation energy (eV)
        call getTics(tcsob,tics)       
        if (eIncident .gt. eIon) then
           intDist=0.0
           jj=2
           do while (sdcsen%ecg(jj) .le. (eIncident-eIon)/2.0)
              deltaE = sdcsen%ecg(jj)-sdcsen%ecg(jj-1)
              intDist = intDist + sdcsen%CS(jj)*deltaE 
              jj = jj+1
           end do
           !Normalises distribution
           if ( intDist .lt. 0.0001 ) then
              sdcsen%CS(:) = 0.0
           else
              sdcsen%CS(:) = sdcsen%CS(:)*tics/intDist 
           end if
        else
           sdcsen%CS(:) = 0.0
        end if
       
    end subroutine normaliseSDCSDist


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: printSdcsDistAtEn
  !Purpose: prints cross sections in imported sdcs type to a file
  !         for plotting. Used to check validity of interpolation.
  !Date last modified: 23/09/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine printSdcsDistAtEn(self)
      implicit none
      type(sdcs)::self
      integer:: ii
      character (len=3)::enString
 
      write(enString,'(I0)') int(self%en_incident)
      open(70, file="sdcsDist"//trim(enString)//".txt")
      write(70,*) "Energy(eV)          CS(a.u)"
      do ii = 1, self%Ncgsdcs
         write(70,*) self%ecg(ii), self%CS(ii)
      end do
      close(70)    

  end subroutine printSdcsDistAtEn



!Additional functions for comparison of two SDCS treatments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: calcProbs
    !Purpose: calculate and compare probabilities of selecting given
    !         ejection energies between ionisation treatments.
    !         Comparison done by writing probabilities to file for 
    !         plotting.
    !Note: mostly used for debugging purposes
    !Date last modified: 12/10/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calcProbs(sdcsBasis,stateBasis,eOut)
        use totalcs_module
        use state_class
        use input_data
        implicit none
        type(basis_sdcs)::sdcsBasis
        type(basis_state)::stateBasis
        real(dp)::eOut !Given ejection energy
        integer:: ii,jj
        real(dp):: sdcsProb, sdcsDistProb, cccDistProb, eIncident
        real(dp):: eIon
        type(sdcs)::sdcsEn, sdcsDistEn
        type(totalcs)::tcsEn
        real(dp), allocatable, dimension(:)::enProb       
        character (len=3)::enString

        write(enString,'(I0)') int(eOut)
        open(20,file="ejEnProb"//trim(enString)//".txt") 
        write(20,*) "eIncident (eV)        sdcs         cccDist"
        
        eIon = 15.96632  !H2 ionisation energy (eV)
        !Use sdcs common energy grid for convenience
        do ii=1,sdcsBasis%Nein-1
           eIncident = sdcsBasis%ein(ii)

           !One probability for each of the three ionisation treatments
           sdcsProb = 0.0
           sdcsDistProb = 0.0
           cccDistProb = 0.0
           if (eOut .lt. eIncident-eIon) then
              !Get distributions for each treatment
              call get_sdcsatein(sdcsBasis,eIncident,sdcsEn)
              call get_csatein(stateBasis,eIncident,tcsEn)
              call getSdcsDistAtEIn(sdcsBasis,eIncident,sdcsDistEn,tcsEn)

              !Energy converted to a.u. in get_csatein for some reason,
              !convert back to atomic units
              !Unfortunate use of global variable required
              tcsEn%en(:) = tcsEn%en(:)*data_in%eV

              !Loop through cross sections to calculate probabilities
              allocate(enProb(sdcsEn%Ncgsdcs)) 
              call normalisesdcs(sdcsEn,enProb,sdcsEn%Ncgsdcs)
              jj=1
              do while (sdcsEn%ecg(jj) .lt. eOut)
                 jj = jj+1
              end do
              jj=jj-1
              sdcsProb = enProb(jj)
              deallocate(enProb)   

              allocate(enProb(sdcsDistEn%Ncgsdcs))
              call normalisesdcs(sdcsDistEn,enProb,sdcsDistEn%Ncgsdcs)
              jj=1
              do while(sdcsDistEn%ecg(jj) .lt. eOut)
                 jj= jj+1
              end do
              jj=jj-1
              sdcsDistProb = enProb(jj)
              deallocate(enProb)
                 
              allocate(enProb(tcsEn%Nmax))
              call normalisetcs(tcsEn,enProb,tcsEn%Nmax)     
              jj= 1
              do while((jj .le. tcsEn%Nmax-1) .and. (tcsEn%en(jj) .lt. eOut))
                 jj=jj+1
              end do
              jj = jj-1
              cccDistProb = enProb(jj) 
              deallocate(enProb)

              !Write probabilities to file
              !write(20,*) eIncident, sdcsProb,  sdcsDistProb, cccDistProb
              write(20,*) eIncident, sdcsProb,  cccDistProb
              call delete_totalcs(tcsEn)
              call destructSdcsAtEIn(sdcsEn)
              call destructSdcsAtEIn(sdcsDistEn)
           else
              !write(20,*) eIncident, sdcsProb,  sdcsDistProb, cccDistProb
              write(20,*) eIncident, sdcsProb,  cccDistProb
           end if 
        end do
        close(20)
    end subroutine calcProbs



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: compareCS
    !Purpose: writes cross sections for different ionisation treatments to 
    !         file for comparison
    !Date last modified: 15/10/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compareCS(sdcsBasis,stateBasis,eIncident) 
        use totalcs_module
        use state_class
        use input_data
        implicit none
        type(basis_sdcs)::sdcsBasis
        type(basis_state)::stateBasis
        type(sdcs)::sdcsEn
        type(totalcs)::tcsEn
        real(dp)::eIncident
        character (len=3)::enString
        character (len=22)::sdcsEnString,sdcsCSString,tcsEnString,tcsEnCSString
        integer::ii

        write(enString,'(I0)') int(eIncident)
        open(20,file="compCS"//trim(enString)//".txt")

        call get_sdcsatein(sdcsBasis,eIncident,sdcsEn)
        call get_csatein(stateBasis,eIncident,tcsEn)
        do ii=1,sdcsEn%Ncgsdcs
           write(sdcsEnString,*) sdcsEn%ecg(ii)
           write(sdcsCSString,*) sdcsEn%CS(ii)

           if (ii .le. tcsEn%Nmax) then
              write(tcsEnString,*) tcsEn%en(ii)
              write(tcsEnCSString,*) tcsEn%cs(ii)
              write(20,*) sdcsEnString(1:5)//sdcsEnString(18:22), sdcsCSString(1:5)//sdcsCSString(18:22), tcsEnString(1:5)//tcsEnString(18:22), tcsEnCSString(1:5)//tcsEnCSString(18:22)
           else if (ii .gt. tcsEn%Nmax) then
              write(20,*) sdcsEnString(1:5)//sdcsEnString(18:22), sdcsCSString(1:5)//sdcsCSString(18:22),0.0000, 0.0000
           end if
        end do  
        call delete_totalcs(tcsEn)
        call destructSdcsAtEIn(sdcsEn)
        close(20)

    end subroutine compareCS




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: compareCSMod
    !Purpose: writes cross sections for different ionisation treatments to 
    !         file for comparison.
    !Note: modified version of above, rewritten to allow for quicker debugging.
    !      Done to meet a deadline.
    !Date last modified: 15/10/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compareCSMod(sdcsBasis,stateBasis,eIncident) 
        use totalcs_module
        use state_class
        use input_data
        implicit none
        type(basis_sdcs)::sdcsBasis
        type(basis_state)::stateBasis
        type(sdcs)::sdcsEn
        type(totalcs)::tcsEn
        real(dp)::eIncident
        character (len=3)::enString
        integer::ii

        write(enString,'(I0)') int(eIncident)
        open(20,file="compCSSDCS"//trim(enString)//".txt") 
        open(50,file="compCSCCC"//trim(enString)//".txt")

        call get_sdcsatein(sdcsBasis,eIncident,sdcsEn)
        call get_csatein(stateBasis,eIncident,tcsEn)
    
        if (data_in%stateIonop .eq. 1) then
           call symtcs(tcsEn, eIncident)
        end if 

        do ii=1,sdcsEn%Ncgsdcs
           if (ii .le. tcsEn%Nmax) then
              write(20,*) sdcsEn%ecg(ii), sdcsEn%CS(ii)
              !Note: tcs values in other units, need to convert to
              ! eV and a_0^2
              write(50,*) tcsEn%en(ii)*data_in%eV, tcsEn%cs(ii)
              print*, tcsEn%en(ii), tcsEn%cs(ii)
           else if (ii .gt. tcsEn%Nmax) then
              write(20,*) sdcsEn%ecg(ii), sdcsEn%CS(ii)
              write(50,*) 0.0000, 0.0000
           end if
        end do  
        call delete_totalcs(tcsEn)
        call destructSdcsAtEIn(sdcsEn)
        close(20)
        close(50)

    end subroutine compareCSMod







    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: compareMeanEjEn
    !Purpose: calculates and compares the mean ejection energies
    !         of the two ionisatoin treatments as a function of 
    !         incident energy
    !Date last modified: 15/10/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compareMeanEjEn(sdcsBasis,stateBasis) 
        use totalcs_module
        use state_class
        use input_data
        implicit none
        type(basis_sdcs)::sdcsBasis
        type(basis_state)::stateBasis
        type(sdcs)::sdcsEn
        type(totalcs)::tcsEn
        integer::ii,jj
        real(dp)::deltaE,eCutoff,eIon, eIncident
        real(dp)::meanSDCS, meanCCC, intSDCS, intCCC


        open(20,file="meanEjEn.txt")
        eIon = 15.96632  !H2 ionisation energy (eV)
        do ii=1,sdcsBasis%Nein-2
           eIncident = sdcsBasis%ein(ii)
           if (eIncident .gt. eIon) then
              eCutoff = (eIncident-eIon)/2.0
              call get_sdcsatein(sdcsBasis,eIncident,sdcsEn)
              call get_csatein(stateBasis,eIncident,tcsEn)
 
              !Calculate mean up to E/2
              !Calculate SDCS mean via integration
              jj=2
              meanSDCS = 0.0
              intSDCS = 0.0
              do while (sdcsEn%ecg(jj) .lt. eCutoff)
                 deltaE = sdcsEn%ecg(jj)-sdcsEn%ecg(jj-1)
                 meanSDCS = meanSDCS + sdcsEn%ecg(jj)*sdcsEn%CS(jj)*deltaE
                 intSDCS = intSDCS + sdcsEn%CS(jj)*deltaE
                 jj=jj+1
              end do
              if (.not. (intSDCS .gt. 0.0)) then
                 !intSDCS = sdcsEn%CS(1)*(sdcsEn%ecg(2)-sdcsEn%ecg(1))
                 meanSDCS = 0.0
              else
                 meanSDCS = meanSDCS/intSDCS
              end if 

              jj=2
              meanCCC = 0.0
              intCCC = 0.0
              do while((tcsEn%en(jj) .lt. eCutoff) .and. (jj .lt. tcsEn%Nmax-1))
                 if (tcsEn%en(jj) .ge. 0.0 ) then
                    meanCCC = meanCCC + tcsEn%cs(jj)*tcsEn%en(jj)
                    intCCC = intCCC + tcsEn%cs(jj)
                 end if
                 jj=jj+1
              end do
              if (.not. (intCCC .gt. 0.0)) then
                 intCCC = tcsEn%cs(1)
              end if
              meanCCC = (meanCCC/intCCC)*data_in%eV

              write(20,*) eIncident, meanSDCS, meanCCC 

              call delete_totalcs(tcsEn)
              call destructSdcsAtEIn(sdcsEn)
           else
              write(20,*) eIncident, 0.0, 0.0
           end if
        end do  

        close(20)
    end subroutine compareMeanEjEn



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: intSdcsIn
    !Purpose: integrates SDCS from input files to get total ionisation
    !         cross section(TICS) and writes results to file
    !         for plotting.
    !Date last modified: 08/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine intSdcsIn(data_input,sdcsBasis)
        use input_data
        implicit none
        type(input)::data_input
        type(basis_sdcs)::sdcsBasis
        integer::ii,kk
        real(dp)::totalCSVal, deltaE, currEn
        real(dp)::cutoffEn, eIon


        eIon = 15.96632  !Ionisation energy of ground state H2
        open(70,file='intSdcsIn.txt')
        do ii = 1,data_input%Nsdcs
           currEn = sdcsBasis%ein(ii)
           cutoffEn = (currEn-eIon)/2.0       
     
           !Integrate SDCS, energies and SDCS start at zero
           totalCSVal = sdcsBasis%readEOut(ii,1)*sdcsBasis%cs(ii,1)
           do kk = 2,size(sdcsBasis%readEOut(ii,:)) 
              if ((sdcsBasis%readEOut(ii,kk) .lt. cutoffEn) .and. (sdcsBasis%readEOut(ii,kk) .gt. 0.0)) then
                 deltaE = sdcsBasis%readEOut(ii,kk) - sdcsBasis%readEOut(ii,kk-1)
                 totalCSVal = totalCSVal + 0.5*(sdcsBasis%cs(ii,kk)+sdcsBasis%cs(ii,kk-1))*deltaE
              end if
           end do
           write(70,*) currEn, totalCSVal
        end do
        close(70) 
    end subroutine intSdcsIn



   !_____________________Functions for constructing analytical SDCS___________________________!
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: createAnalyticSdcs
   !Purpose: fills the sdcsBasis with analytical Sdcs given in Garvey and
   !         Green (1976).
   !Date last modified: 26/03/2021
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine createAnalyticSdcs(sdcsBasis)
       use input_data 
       implicit none
       type(basis_sdcs)::sdcsBasis
       integer:: ii, jj, numEIn, numEOut
       real(dp):: currEIn, currEOut, enSpacing
       real(dp):: Ae, Te, gammaE, I, sigma 
 
       if (associated(sdcsBasis%CS)) then
          deallocate(sdcsBasis%CS)
       end if
       if (associated(sdcsBasis%ein)) then
          deallocate(sdcsBasis%ein)
       end if
       if (associated(sdcsBasis%ecg)) then
          deallocate(sdcsBasis%ecg)
       end if

       enSpacing = 700.0/10000.0 !spacing of uniform SDCS grid (eV)
       sdcsBasis%Ncgsdcs = int(700.0/enSpacing)
       allocate(sdcsBasis%ecg(sdcsBasis%Ncgsdcs))
       allocate(sdcsBasis%ein(sdcsBasis%Ncgsdcs))
       !Use new, uniform energy grid  
       sdcsBasis%ecg(1) = 0.0
       sdcsBasis%ein(1) = 0.5114
       do ii = 2, sdcsBasis%Ncgsdcs 
          sdcsBasis%ecg(ii) = sdcsBasis%ecg(ii-1) + enSpacing
          sdcsBasis%ein(ii) = sdcsBasis%ein(ii-1) + enSpacing
       end do
       sdcsBasis%Nein = SIZE(sdcsBasis%ein)

       numEIn = sdcsBasis%Ncgsdcs
       numEOut = sdcsBasis%Ncgsdcs
       allocate(sdcsBasis%CS(numEIn, numEOut))


       I = 16.0   !Ionisation potential(eV)
       sigma = (1e-16)*((0.01)**2)/(data_in%bohrRadius**2)

       do ii = 1, numEIn
          currEIn = sdcsBasis%ein(ii)
          do jj = 1, numEOut
             currEOut = sdcsBasis%ecg(jj)

             !Factors used in analytical formula
             Ae = (2.871/currEIn)*sigma*log(currEIn/0.5114)
             Te = 1.870-1000.0/(currEIn+3*I)
             gammaE = 7.070*currEIn/(currEIn-7.700)
                    
             sdcsBasis%CS(jj,ii) = Ae*(gammaE**2/((currEOut-Te)**2 + gammaE**2))
          end do
       end do

       sdcsBasis%Ninsdcs = sdcsBasis%Ncgsdcs 
       sdcsBasis%Ninsdcs_first = 1


   end subroutine createAnalyticSdcs 





end module sdcs_module
