subroutine populatestates(statebasis,tcsbasis)

  use numbers
  use input_data         ! definitions for input data type, construct routine, object data_in
  use state_class        ! defines states (and basis of them) with operations on them
  use totalcs_module     !  reading totalcs files

  implicit none


  integer:: Nlarge
  parameter(Nlarge = 1000)
  type(basis_totalcs), intent(in):: tcsbasis
  type(basis_state), intent(inout):: statebasis
  type(state):: psFormation
  integer:: Nmaxall
  logical:: ex
  integer:: i
  integer:: N, M,ipar
  logical:: ion
  real(dp):: S, twoelen, enex, en, eneV, threshold, ground
  integer:: inum
  real(dp), dimension(tcsbasis%n):: ar_ein,ar_cs
  integer:: Nmax, l, j, idf, iac, inext, io
  character(len=10):: stlabel
  real(dp), dimension(Nlarge):: vdf,eindf
  character(len=10):: filenamest
  character(len=4):: stnumber
  character(len=50):: filename

  print*
  print*, 'Populate states'

  inum = tcsbasis%n
!  print*, 'inum = ', inum

  ! read States_list
  
  if (data_in%posmode .eq. 1) then
    filename = '../pos-H2_DCS/states.core_parts'
  else
    filename = TRIM(data_in%DATApath)//'/'//'States_list'
  end if
  
  inquire(file=filename,exist=ex)
  if (ex) then
     print*, 'States_list found: intialize states'
  else
     print*,'States_list not found, stopping' 
     stop
  endif
  open(25,file=filename,action='read')

  Nmaxall = tcsbasis%b(inum)%Nmax   ! get it from the totalcs file with the largest energy
  if (data_in%posmode .eq. 1) then
    Nmaxall = 0
    do 
      read(25, *, iostat=io)
      if (io/=0) exit
      Nmaxall = Nmaxall + 1
    end do
    rewind(25)
    Nmaxall = Nmaxall - 1
  end if
  

  print*, 'Nmaxall = ', Nmaxall
  read(25,*)
  ground = 0.0
  do i=1, Nmaxall 
     if (data_in%posmode .eq. 1) then
       read(25, *) N, ipar, M, S, twoelen, stlabel
       print *, stlabel
      if (i .eq. 1) then 
        ground = twoelen
      endif
       enex = (twoelen - ground)*13.6*2
       eneV = enex - 15.96632
       en = eneV/(2*13.6)
     else
      read(25,*) N, M,ipar,S, stlabel,twoelen, enex, en, eneV
     end if
     l = -1
     ion = .false.
     if (eneV .gt. 0.0d0) then
        ion = .true.
     end if
     !print *, size(statebasis%b)
     print *, 'calling init state'
     call init_state(statebasis%b(i), trim(stlabel),eneV,enex,l,M,ipar,S,ion)
     print *, 'finished init state'
  enddo

  if (data_in%posmode .eq. 1) then
    statebasis%b(:)%inum = inum
  end if


  close(25)
  !---------------------
  ! populate cross sections from totalcs files
  print*, 'write cross sections from totalcs files to states '
  do n=1,Nmaxall
    !print*, statebasis%b(n)%inum 
     do i=1,inum
        Nmax = tcsbasis%b(i)%Nmax
        if(n .le. Nmax) then
           ar_ein(i) = tcsbasis%b(i)%en_incident
           ar_cs(i) = tcsbasis%b(i)%cs(n)
           !          print*, '----->>> ', tcsbasis%b(i)%en_incident, tcsbasis%b(i)%cs(n)
        else
           ar_ein(i) = tcsbasis%b(i)%en_incident
           ar_cs(i) = 0d0
        endif
     enddo
     l = 1
     if( n .eq. 1 ) l = 0  ! do not add the threshild energy for the elastic channel
     call set_cs(statebasis%b(n),inum,ar_ein,ar_cs,l) ! set CS and add a point at the threshold energy
     !print*, statebasis%b(n)%inum 

  enddo
!
!-----------------------------------------------------
! read diss. fractions and write them into states
  print*, 'Read dissociation fractions and write them to states,   number of states with diss. fractions: Ndf =', data_in%Ndf
  do i=1,data_in%Ndf
     open(27,file=data_in%filename_df(i),action='read')
     j = 0
     do 
        j = j + 1
        read(27,*,END=127,ERR=127) eindf(j),vdf(j)
     enddo
127  continue
     close(27)
     idf = j - 1
     
     ! find state with this label
     do n=1,Nmaxall
!        print*, TRIM(data_in%filename_df(i))
        j = LEN(TRIM(data_in%DATApath))+ 1
        print *, statebasis%b(n)%stlabel, ' COMPARE WITH ', data_in%filename_df(i)(j+4:j+10)
        !print *, data_in%filename_df(i)
        !print*, LEN(trim(adjustl(statebasis%b(n)%stlabel))), LEN(trim(adjustl(data_in%filename_df(i)(j+4:j+10))))
        if(trim(adjustl(statebasis%b(n)%stlabel(1:4))) .eq. trim(adjustl(data_in%filename_df(i)(j+4:j+10)))) then 
           print*,'found for diss.fractions:', n,statebasis%b(n)%stlabel  
           exit
        else if(trim(adjustl(statebasis%b(n)%stlabel(1:5))) .eq. trim(adjustl(data_in%filename_df(i)(j+4:j+10)))) then   
           print*,'found for diss.fractions:', n,statebasis%b(n)%stlabel  
           exit
        endif
     enddo
     if(n .eq. Nmaxall+1) then
        print*,'tcstostate.f90: populatestate() 1: could not find a state with this label:',  data_in%filename_df(i)(j+4:j+10)
        stop
     endif
     
     ! for states that have nonzero M  ( P,D, ...) the next state, n+1, 
     ! should also have the same label and quantum numbers (but -M), 
     ! need to write into that state too
     ! Check that the next state has the same label
     inext = 0
     if(statebasis%b(n)%stlabel .eq. statebasis%b(n+1)%stlabel) then
        inext = 1
     endif
     
     do j=0,inext
        if (associated(statebasis%b(n+j)%eindf)) deallocate(statebasis%b(n+j)%eindf)
        if (associated(statebasis%b(n+j)%df)) deallocate(statebasis%b(n+j)%df)
        allocate(statebasis%b(n+j)%df(idf),statebasis%b(n+j)%eindf(idf))

        call set_df(statebasis%b(n+j),idf,eindf(1:idf),vdf(1:idf))
     enddo

  enddo
!-----------------------------------------------------

  if (data_in%posmode .ne. 1) then
   ! read AN CS files write them into states
     print*, 'read AN cross sections from additional CS files'
     do i=1,data_in%Naddics
        open(27,file=data_in%filename_addics(i),action='read')
        do j=1,6
           read(27,*)
        enddo
        j = 2   ! start from 2 and set CS at these points to zero to avoid problems with cubic splines 
        do 
           j = j + 1
           read(27,*,END=137,ERR=137) eindf(j),vdf(j)
        enddo
   137  continue
        close(27)
    
        iac = j - 1  ! number of energies
    
        eindf(1) = eindf(3) - 0.02
        vdf(1) = 0d0
        eindf(2) = eindf(3) - 0.01
        vdf(2) = 0d0
   
        ! find state with this label
        do n=1,Nmaxall
           j = LEN(TRIM(data_in%DATApath)) + 1
           if(statebasis%b(n)%stlabel .eq. TRIM(data_in%filename_addics(i)(j+5:j+9)) .or. statebasis%b(n)%stlabel .eq. TRIM(data_in%filename_addics(i)(j+5:j+8)) ) then 
              print*,'found for add. CS:', n,statebasis%b(n)%stlabel  
              exit
           endif
        enddo
        if(n .eq. Nmaxall+1) then
           print*,'tcstostate.f90: populatestate() 2: could not find a state with this label:',  data_in%filename_addics(i)(5:9)
           stop
        endif
        if( n .eq. 1) then ! elastic channel vdf(1) = vdf(3)
           vdf(2) = vdf(3)
        endif
   
        ! for states that have nonzero M  ( P,D, ...) the next state, n+1, 
        ! should also have the same label and quantum numbers (but -M), 
        ! need to write into that state too
        ! Check that the next state has the same label
        inext = 0
        if(statebasis%b(n)%stlabel .eq. statebasis%b(n+1)%stlabel) then
           inext = 1
        endif
        
        do j=0,inext
           if (associated(statebasis%b(n+j)%ein)) deallocate(statebasis%b(n+j)%ein)
           if (associated(statebasis%b(n+j)%cs)) deallocate(statebasis%b(n+j)%cs)
           allocate(statebasis%b(n+j)%cs(iac),statebasis%b(n+j)%ein(iac))
   !        call set_csad(statebasis%b(n+j),iac,eindf(1:iac),vdf(1:iac)/dble(1+inext))  ! divide by 2 cross section for non M=0 states 
           statebasis%b(n+j)%inum = iac 
           statebasis%b(n+j)%ein(1:iac) = eindf(1:iac)
           statebasis%b(n+j)%cs(1:iac) = vdf(1:iac)/dble(1+inext)
   
           statebasis%b(n+j)%enex = statebasis%b(n)%enex
           statebasis%b(n+j)%en = statebasis%b(n)%en      
    
           !statebasis%b(n+j)%enex = eindf(3)   ! reset excitation energy to the first entry in add. CS files
           !statebasis%b(n+j)%en = eindf(3) + statebasis%b(1)%en  !Modify ionisation energy
        enddo
   
        !ground state energy of H2+
        threshold=-0.57644d0
        do j=1, size(statebasis%b)
           if (statebasis%b(j)%en .gt. threshold) then
              statebasis%b(j)%stlabel="ION-H2"
           end if
        end do
   
   
   ! ------- Test ------
        l = 1
        !Note: data_in is (unfortunately) a global variable defined in input_data.f90
        if ((l .eq. 1) .and. (data_in%debugOp .eq. 1)) then
           write(stnumber,'(i4.4)') n
           filenamest ="oldst_"//stnumber     
           call print_state(statebasis%b(n),filenamest)
        endif
     enddo
  else
     !Read in scaling factor
  ! allocate(statebasis%b(statebasis%n)%ein(statebasis%b(statebasis%n)%inum)


     call readPsFormCs(psFormation, statebasis%b(statebasis%n)%ein, statebasis%b(statebasis%n)%inum, "../data-2/PsTCS.txt")
     
!     print*, 'psF%inum', psFormation%inum
     psFormation%inum = inum
!     print*, 'stlabel: ', psFormation%stlabel
!     print*, 'en: ', psFormation%en
!     print*, 'enex: ', psFormation%enex
!     print*, 'el, m: ', psFormation%l, psFormation%m
!     print*, 'ipar, spin, inum: ', psFormation%ipar, psFormation%spin, psFormation%inum
!     print*, 'cs: ', psFormation%cs

!     print*, 'ALLOCATE psFormation%df'
     allocate(psFormation%df(inum), psFormation%eindf(inum), psFormation%enEin(inum))
     psFormation%df(:) = 0
     psFormation%enEin(:) = 0
     psFormation%df(:) = psFormation%eindf(:)
!     print*, 'enEin: ', psFormation%enEin
!     print*, 'df, eindf: ', psFormation%df, psFormation%eindf

!     print*, 'resolved, v, dissThresh: ', psFormation%resolved, psFormation%v, psFormation%dissThresh
!     print*, 'ion: ', psFormation%ion
!     print*, 'psform: ', psFormation%psFormation

     call addState(psFormation, statebasis)
!     stop 
  end if
!  print*, 'statebasis%n', statebasis%n 

!  do i=1, statebasis%n
    ! print*, statebasis%b(i)%inum
    ! if (associated(statebasis%b(i)%cs)) print*, "CS"
    ! if (associated(statebasis%b(i)%ein)) print*, "EIN"
    ! if (associated(statebasis%b(i)%df)) print*, "DF"
    ! if (associated(statebasis%b(i)%eindf)) print*, "EINDF"
!  end do

end subroutine populatestates





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Function: populateStatesVcs
!Purpose: uses vibrational cross section data and data for electronic states
!         summed over vibrational states to create a set of vibrationally 
!         resolved 'states' to use in the simulation. 
!Date last modified: 02/08/2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine populateStatesVcs(statebasis,tcsbasis)

    use numbers
    use input_data         ! definitions for input data type, construct routine, object data_in
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files

    implicit none

    integer:: Nlarge
    parameter(Nlarge = 1000)
    type(basis_totalcs), intent(in):: tcsbasis
    type(basis_state), intent(inout):: statebasis
    type(state), allocatable, dimension(:):: tempbasis, leftovers
    logical, allocatable, dimension(:):: resolved, ignoreState 
    real(dp), allocatable, dimension(:)::newGrid
    integer:: ii,jj,kk, numEnergies, n, j, totalStates
    integer:: counter,  Nmax, numElastic, numVibStatesTot
    integer, allocatable, dimension(:):: numVibStates  !Stores number of vibrational levels per electronic state
    real(dp), allocatable, dimension(:)::einArray, csArray, intpCs, stateEn, dfArray
    integer, allocatable, dimension(:):: vArray
    integer:: vVal, inext
    integer:: inum, Nmaxall, vfVal, gridLength, halfInt, nfile, numLines
    character(len=10)::str1, str2, str3, str4, str5, str6, str7, str8, str9, str10
    character(len=10)::str11, str12  !Used for file IO
    character(len=10)::enString
    character(len=11)::threshEn
    real(dp)::deltaE, eMin, eMax, csVal, maxEnInVcs, minEnInVcs, einVal, dissDiff, dissThreshGround
    real(dp)::groundEn, en1, en2
    character(len=80)::filename
    logical:: ex

    do j =1, size(statebasis%b)
       if ((statebasis%b(j)%stlabel .eq. '?') .or. (statebasis%b(j)%stlabel(1:3) .eq. 'ION')) then
          statebasis%b(j)%v = -2
       end if
    end do

    inum = tcsbasis%n
    Nmaxall = tcsbasis%b(inum)%Nmax   ! get length from the totalcs file with the largest energy

    !Create new data structures
    !numVibStatesTot = data_in%numVib*statebasis%n
    numVibStatesTot = data_in%Nvcs + data_in%NvcsNew
    allocate(tempbasis(numVibStatesTot))  
    allocate(resolved(statebasis%n))
    allocate(ignoreState(statebasis%n))
    allocate(numVibStates(statebasis%n))
    resolved(:) = .false.
    numVibStates(:) = 0  !Number of vibrational levels per electronic state 
    ignoreState(:) = .false.

!    !Find maximum energy in input files
!    filename = data_in%filename_vcsNew(1)
!    open(27,file=filename)
!    eMax = 0d0 !maximum energy in input file
!    do ii = 1,10
!       read(27,*)
!    end do

!     ii = 1
!     do 
!        ii = ii +1
!        read(27,*,END=120,ERR=122) eMax, csVal
!     end do
! 122 print*, "Error in VCS File IO, function=populateStatesVcs, stopping."
!     close(27)
!     stop
! 
! 120 continue
!     close(27)
    
    !Create new common energy grid for vibrationally resolved state type
    eMin = 0d0                 !Minimum energy in input files   
    eMax = tcsbasis%b(tcsbasis%n)%en_incident !Largest totalcs energy
    deltaE = 0.5      !Make grid quite fine
    gridLength = int((eMax-eMin)/deltaE)+1
    allocate(newGrid(gridLength))
    newGrid(1) = 0d0
    do ii = 2, gridLength
       newGrid(ii) = newGrid(ii-1) + deltaE
    end do

    maxEnInVcs = 0.0
    minEnInVcs = 0.0
    
    !Read in ground electonic state energies for H2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    filename = data_in%filename_vcsEn(1)
    open(23,file=filename,action='read')
    do ii = 1, 5
       read(23,*)
    end do
    ii = 0
    do  
      read(23,*,END=70,ERR=70) vfVal, en1, en2
      ii = ii+1
    end do 
70  continue
    close(23)
    numLines = ii
    allocate(stateEn(numLines))

    open(23,file=filename,action='read') 
    do ii = 1, 5
       read(23,*)
    end do
    do ii = 1,numLines
       read(23,*) vfVal, en1, stateEn(ii)  !stateEn in eV 
    end do
    close(23)
    groundEn = stateEn(1)  !Ground state energy in eV


    !Read in dissociation threshold for ground electronic state from file
    filename=data_in%filename_vcsEn(2)
    dissDiff=0.0_dp
    dissThreshGround=0.0_dp
    open(28,file=filename)
		print*, filename
    read(28,*)
    read(28,*)
    read(28,*) vfVal, en1, en2, dissDiff 
    close(28)
    dissThreshGround = abs(dissDiff) 

    !Read in vcs for elastic (electronic) scattering (X1Sg -> X1Sg) (vi -> vf)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'read vibrational excitation cross sections (elastic)'
    do ii = 1, data_in%Nvcs

       filename = data_in%filename_vcs(ii)
       ! open vcs file: filename
       inquire(file=filename,exist=ex)
       if (ex) then
          nfile = 100 
       else
          print*, filename, 'is not found, stopping'
          stop
       endif
       open(nfile,file=filename,action='read')
    
       !Ignore first two lines 
       read(nfile,*)
       read(nfile,*)
       read(nfile,*)   vfVal !Final vibrational state
       read(nfile,*)
       !Find number of data points in file
       jj = 0
       do
          read(nfile,*, ERR=110,END=110) einVal, csVal
          if (jj .eq. 0) then
             minEnInVcs = einVal
          end if
          maxEnInVcs = einVal
          jj = jj + 1
       enddo
110    continue

       Nmax = jj
       allocate(csArray(Nmax),einArray(Nmax))
       rewind(nfile)
       read(nfile,*)
       read(nfile,*)
       read(nfile,*)
       read(nfile,*)
       do jj=1,Nmax
          read(nfile,*) einArray(jj), csArray(jj)
       enddo
       close(nfile)

       !Find state with label given in filename
       do n=1,Nmaxall
          j = LEN(TRIM(data_in%DATApath)) + 1
          if(statebasis%b(n)%stlabel .eq. TRIM(data_in%filename_vcs(ii)(j+8:j+11)) ) then 
             print*,'found for VCS (index, state):', n,statebasis%b(n)%stlabel  
             exit
          endif
       enddo
       if(n .eq. Nmaxall+1) then
          print*, data_in%filename_VcsNew(ii)
          print*,'tcstostate.f90: populateStateVcs() 2: could not find a state with this label:',  data_in%filename_Vcs(ii)(j+8:j+11)
          stop
       endif
       resolved(n) = .true.   !State with given label is now vibrationally resolved
       numVibStates(n) = numVibStates(n) + 1

       inext = 0
       if(statebasis%b(n)%stlabel .eq. statebasis%b(n+1)%stlabel) then
          inext = 1
          !Vibrational CS files are summed over values of M, need to avoid double counting
          ignoreState(n+inext) = .true.   
       endif 


       !Create a new state type with this electronic label and vf
       call copy_state(tempbasis(ii),statebasis%b(n))        
       tempbasis(ii)%v = vfVal     !vf in file 
       tempbasis(ii)%en = stateEn(vfVal+1) !State energy (possibly negative) in eV
       tempbasis(ii)%dissThresh= dissThreshGround !Dissociation threshold (positive), measure from assumed initial state

       !Interpolate cross sections to new energy grid 
       allocate(intpCs(size(newGrid)))
       !Interpolate over each half of the new grid due to INTRPL's restrictions on array sizes
       halfInt = int(real(size(newGrid))/2.0)
       numEnergies = size(einArray)
       call INTRPL(numEnergies,einArray,csArray,halfInt,newGrid(1:halfInt),intpCs(1:halfInt))
       call INTRPL(numEnergies,einArray,csArray,size(newGrid)-halfInt,newGrid(halfInt+1:size(newGrid)),intpCs(halfInt+1:size(newGrid)))

       !Set cross sections at energies outside bounds of input file to zero 
       do kk=1, size(newGrid)
          if (newGrid(kk) .gt. maxEnInVcs) then
             intpCs(kk) = 0d0
          else if (newGrid(kk) .lt. minEnInVcs) then
             intpCs(kk) = 0d0
          end if
       end do
      ! print*, "Interpolated: ", vfVal, stateEn(vfVal+1) 
      ! do kk = 1, size(newGrid)
      !    print*, newGrid(kk), intpCs(kk)
      ! end do

       !Copy cross sections to state type
       deallocate(tempbasis(ii)%ein,tempbasis(ii)%cs)
       allocate(tempbasis(ii)%ein(size(newGrid)),tempbasis(ii)%cs(size(newGrid)))
       allocate(dfArray(size(newGrid)))
       tempbasis(ii)%ein(:) = newGrid(:)
       tempbasis(ii)%cs(:) = intpCs(:)
       tempbasis(ii)%inum = size(newGrid)
       tempbasis(ii)%resolved = .true.
       dfArray(:) = 0.0
       call set_df(tempbasis(ii), size(newGrid), newGrid, dfArray) 
       !Include new excitation energies for given v state
       tempbasis(ii)%enex = stateEn(tempbasis(ii)%v+1) - groundEn
       !print*, "ENEX: ", tempbasis(ii)%enex
       deallocate(intpCs,csArray,eInArray,dfArray)
    end do
    numElastic = ii-1

    !Read in vcs files for bound electronic excitations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !Read vcs files and write them into statebasis
    print*, 'read vibrational excitation cross sections (inelastic)'
    do ii = 1, data_in%NvcsNew
       allocate(einArray(Nlarge),csArray(Nlarge))     

       !Read in file data
       filename = data_in%filename_vcsNew(ii)
       open(27,file=filename,action='read')
       do jj=1,4   !Ignore first six lines
          read(27,*)
       end do
       !Note that 'read' uses both spaces and commas as delimiters
       read(27,*) str1, str2, str3, str4, str5, str6, str7, str8, str9, str10, str11, str12
       if (str12(2:4) .eq. 'dis') then
          read(27,*)
          read(27,*)
          read(27,*) str1, str2, threshEn !Read threshold energy 
       else
          read(27,*)
          read(27,*) str1, str2, str3, enString !Read state energy
          read(27,*) str1,  threshEn !Read threshold energy 
       end if
       do jj = 8, 9
          read(27,*)
       end do     
       
       if (len(TRIM(str12)).eq. 5) then
          read(str12(4:4),*) vfVal
       else if (len(TRIM(str12)).eq. 6) then
          read(str12(4:5),*) vfVal
       else 
          vfVal = -1 !Dissociative state
       end if

       !jj = 2   ! start from 2 and set CS at these points to zero to avoid problems with cubic splines 
       jj = 0
       do 
          jj = jj + 1
          read(27,*,END=137,ERR=137) einArray(jj), csArray(jj)
          if (jj .eq. 1) then
             minEnInVcs = einArray(jj)
          end if
       enddo
137    continue
       close(27)
   
       numEnergies = jj-1 ! number of energies
       maxEnInVcs = einArray(numEnergies)  !highest energy in input file

       !Used to make interpolation smoother  
       !numEnergies = jj
       !einArray(1) = einArray(3) - 0.02
       !csArray(1) = 0d0
       !einArray(2) = einArray(3) - 0.01
       !csArray(2) = 0d0

       !Find state with given electronic label 
       do n=1,Nmaxall
          j = LEN(TRIM(data_in%DATApathVcs)) + 1
          if((statebasis%b(n)%stlabel .eq. TRIM(data_in%filename_vcsNew(ii)(j+13:j+16))) .or. (statebasis%b(n)%stlabel .eq. TRIM(data_in%filename_vcsNew(ii)(j+13:j+17))) ) then 
             print*,'found for VCS (index, state):', n,statebasis%b(n)%stlabel  
             exit
          endif
       enddo
       if(n .eq. Nmaxall+1) then
          print*, data_in%filename_VcsNew(ii)
          print*,'tcstostate.f90: populateStateVcs() 2: could not find a state with this label:',  data_in%filename_VcsNew(ii)(j+13:j+16)
          stop
       endif
       resolved(n) = .true.   !State with given label is now vibrationally resolved
       numVibStates(n) = numVibStates(n) + 1

       inext = 0
       if(statebasis%b(n)%stlabel .eq. statebasis%b(n+1)%stlabel) then
          inext = 1
          !Vibrational CS files are summed over values of M, need to avoid double counting
          ignorestate(n+inext) = .true.   
       endif 

       !Create a new state type with this electronic label and vf
       call copy_state(tempbasis(ii+numElastic),statebasis%b(n))        
       tempbasis(ii+numElastic)%v = vfVal     !vf in file 

       !Calculate new state energy.
       if (vfVal .ge. 0) then
          read(enString,*) tempbasis(ii+numElastic)%en
          tempbasis(ii+numElastic)%en = tempbasis(ii+numElastic)%en*data_in%eV  !Convert to eV, file value in Hartees
       else if (vfVal .eq. -1) then
          !Use threshold energy as excitation energy temporarily, given in eV
          read(threshEn,*) tempbasis(ii+numElastic)%enex  
          tempbasis(ii+numElastic)%dissThresh = tempbasis(ii+numElastic)%enex 
          tempbasis(ii+numElastic)%en = tempbasis(1)%enex + groundEn 
       end if

       allocate(intpCs(size(newGrid)))
       !Calculate new set of excitation CS over new energy grid.
       !Interpolate over each half of the new grid due to INTRPL's restrictions on array sizes
       !halfInt = int(real(size(newGrid))/2.0)
       call INTRPL(numEnergies,einArray,csArray,size(newGrid),newGrid,intpCs)
       !call INTRPL(numEnergies,einArray,csArray,size(newGrid)-halfInt,newGrid(halfInt+1:size(newGrid)),intpCs(halfInt+1:size(newGrid)))

       !Some files only have data up to certain energies, above these the CS are essentially zero.
       do kk=1, size(newGrid)
          if (newGrid(kk) .gt. maxEnInVcs) then
             intpCs(kk) = 0d0
          else if (newGrid(kk) .lt. minEnInVcs) then
             intpCs(kk) = 0d0
          end if
       end do

       !Copy cross sections to state type
       deallocate(tempbasis(ii+numElastic)%ein,tempbasis(ii+numElastic)%cs)
       allocate(tempbasis(ii+numElastic)%ein(size(newGrid)),tempbasis(ii+numElastic)%cs(size(newGrid)))
       allocate(dfArray(size(newGrid)))
       tempbasis(ii+numElastic)%ein(:) = newGrid(:)
       tempbasis(ii+numElastic)%cs(:) = intpCs(:)
       tempbasis(ii+numElastic)%inum = size(newGrid)
       tempbasis(ii+numElastic)%resolved = .true.
       dfArray(:) = 0.0   !Vibrationally resolved fractions read in below
       call set_df(tempbasis(ii+numElastic), size(newGrid), newGrid, dfArray) 
       !tempbasis(ii+numElastic)%idf = tempbasis(ii+numElastic)%idf -1 !set_df adds 1 to account for idf = 0 case

       !Include new excitation energies for given v state. Don't do so for dissociative states (v=-1)
       if (tempbasis(ii+numElastic)%v .ge. 0) then
          tempbasis(ii+numElastic)%enex = tempbasis(ii+numElastic)%en - groundEn  !Assuming excitation from ground n=0,vi=0
       end if
       deallocate(einArray,csArray,intpCs,dfArray)
    end do

    !Recombine states to produce vibrationally resolved statebasis
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    jj = 0
    do ii = 1, size(resolved)
       if (.not. resolved(ii)) then
          jj = jj+1
       end if
    end do


    !Interpolate cs for remaining (pseudo)states to new grid and copy 
    allocate(intpCs(size(newGrid)), leftovers(jj))
    jj = 1
    do ii = 1, size(resolved)
       if ((.not. resolved(ii)) .and. (.not. ignoreState(ii))) then
          halfInt = int(real(size(newGrid))/2.0)
          call INTRPL(statebasis%b(ii)%inum,statebasis%b(ii)%ein,statebasis%b(ii)%cs,halfInt,newGrid(1:halfInt),intpCs(1:halfInt))
          call INTRPL(statebasis%b(ii)%inum,statebasis%b(ii)%ein,statebasis%b(ii)%cs,size(newGrid)-halfInt,newGrid(halfInt+1:size(newGrid)),intpCs(halfInt+1:size(newGrid)))
  
          do kk = 1, size(newGrid)
             !If new grid energy less than lowest in original grid, set CS to zero
             if (newGrid(kk) .lt. statebasis%b(ii)%ein(1)) then
                intpCs(kk) = 0d0
             else if (newGrid(kk) .gt. statebasis%b(ii)%ein(statebasis%b(ii)%inum)) then
                !If new grid energy greater than max in orginal grid, set to zero
                intpCs(kk) = 0d0
             end if
          end do
          call copy_state(leftovers(jj),statebasis%b(ii))
          deallocate(leftovers(jj)%cs,leftovers(jj)%ein)
          allocate(leftovers(jj)%cs(size(newGrid)),leftovers(jj)%ein(size(newGrid)))
          allocate(dfArray(size(newGrid)))
          leftovers(jj)%ein(:) = newGrid(:)
          leftovers(jj)%cs(:) = intpCs(:)
          leftovers(jj)%inum = size(newGrid)
          !Interpolate old dissociation fractions to new energy grid and copy
          !print*, "Fix dissociation fraction issues, tcstostate.f90"
          !print*, "Issue 1: eindf not associated, need eindf and ein to be the same"
          !print*, "eindf is not allocated for some states, these states have idf=0 by default"
          !print*, "These are the first 5 electronic states, which are either ground or triplet states"
          !print*, "If a state is triplet, its dissociation fraction is set to 1 in get_csatein"
          !print*, "Do this setting df for triplet states to 1 now. Have an array of ones over the entire newGrid"
          !stop
         
          if (statebasis%b(ii)%idf .gt. 0) then
             call INTRPL(statebasis%b(ii)%idf,statebasis%b(ii)%eindf,statebasis%b(ii)%df,halfInt,newGrid(1:halfInt),dfArray(1:halfInt))
             call INTRPL(statebasis%b(ii)%idf,statebasis%b(ii)%eindf,statebasis%b(ii)%df,size(newGrid)-halfInt,newGrid(halfInt+1:size(newGrid)),dfArray(halfInt+1:size(newGrid)))
          else if (statebasis%b(ii)%idf .eq. 0) then
             !Found triplet state with no dissociation fractions defined.
             !Allocate arrays and set them all equal to 1 over newGrid
             dfArray(:) = 1.0 
          end if

          call set_df(leftovers(jj), size(newGrid), newGrid, dfArray) 
          deallocate(dfArray)
          jj = jj +1
       end if 
    end do 
    jj = jj-1 
    deallocate(intpCs)

    !Copy over threshold energies for dissociative excitation
    do ii = 1, size(tempbasis) 
       if (tempbasis(ii)%v .eq. -1) then
          do jj = 1, size(tempbasis)
             if (tempbasis(jj)%stlabel .eq. tempbasis(ii)%stlabel) then
                tempbasis(jj)%dissThresh = tempbasis(ii)%dissThresh
             end if
          end do
       end if
    end do
 
    !Deconstruct old basis and copy new one to statebasis 
    call destruct_basis(statebasis)
    totalStates = numVibStatesTot + jj
    call new_basis(statebasis,numVibStatesTot+jj,data_in%tcstype)
    ii = 1
    !jj = 1
    !kk = 1
    counter = 1

    !Order states are copied in doesn't matter, as they are sorted by energy later.
    do ii = 1, size(tempbasis)
       !Ignore dissociative states for now. Treated later by reading in dissociative pseudostates.
       call copy_state(statebasis%b(counter), tempbasis(ii))
       counter = counter + 1
    end do
    do ii = 1, jj   !Copy over non-vibrationally resolved states
       call copy_state(statebasis%b(counter), leftovers(ii))
       counter = counter+1
    end do

    !Copy over dissociation thresholds 
    do ii = 1, statebasis%n
      do jj = 1, statebasis%n
         if ((statebasis%b(ii)%stlabel .eq. statebasis%b(jj)%stlabel) .and. (statebasis%b(jj)%v .eq. -1)) then 
            statebasis%b(ii)%dissThresh = statebasis%b(jj)%dissThresh  
         end if
      end do
    end do

    !Read vibrationally resolved dissociation fractions into new basis
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ii = 1, data_in%NDissFrac  
       allocate(vArray(Nlarge),dfArray(Nlarge))     

       !Read in file data
       filename = data_in%filename_dfVib(ii)
       open(27,file=filename,action='read')
       do jj=1,8   !Ignore first eight lines
          read(27,*)
       end do

       jj = 0
       do 
          jj = jj + 1
          read(27,*,END=120,ERR=120) vArray(jj), dfArray(jj)
       enddo
120    continue
       close(27)

       do vVal=0, jj-2
          !Search for state with given state label and vibrational level
          do n=1, numVibStatesTot
             j = LEN(TRIM(data_in%DATApathDf)) + 1
             if((statebasis%b(n)%stlabel .eq. TRIM(data_in%filename_dfVib(ii)(j+22:j+25)) .or. statebasis%b(n)%stlabel .eq. TRIM(data_in%filename_dfVib(ii)(j+22:j+26))) &
                .and. (statebasis%b(n)%v .eq. vVal)) then 
                print*,'found for dissociation fraction (index, state, v):', n,statebasis%b(n)%stlabel, vVal 
                exit
             endif
          enddo
          if(n .eq. numVibStatesTot+1) then
             print*, data_in%filename_dfVib(ii)
             print*,'tcstostate.f90: populateStateVcs() 2: could not find a state with this label and vib state:',  data_in%filename_dfVib(ii)(j+22:j+25), vVal
             stop
          endif

          !Write dissociation fraction to state
          statebasis%b(n)%df(:) = dfArray(vVal+1)
       end do
   
       deallocate(vArray, dfArray)
    end do
 
    !First numVibStatesTot states are vibrationally resolved. 
    !All triplet states eventually dissociate, set df to 1.
    do ii = 1, numVibStatesTot
       if ((statebasis%b(ii)%v .eq. -1) .or. (int(statebasis%b(ii)%spin) .eq. 1)) then
          statebasis%b(ii)%df(:) = 1.0
       end if
    end do
 

 
    !do ii =1, totalStates
    !      print*, "State: ", ii
    !      print*, statebasis%b(ii)%stlabel, statebasis%b(ii)%v
    !      print*, statebasis%b(ii)%en, statebasis%b(ii)%enex
    !end do
    !stop


    !Sort newly filled basis by excitation energy
    !call sort_by_energy(statebasis)    
    call sortEnergies(statebasis)

    !call sort_by_energy(statebasis)
    !Clean up memory
    if (allocated(leftovers)) then
       do ii = 1, size(leftovers)
          call destruct_state(leftovers(ii))
       end do
       deallocate(leftovers)
    end if

    if (allocated(tempbasis)) then
       do ii = 1, size(tempbasis)
          call destruct_state(tempbasis(ii))
       end do
       deallocate(tempbasis)
    end if


    deallocate(resolved,newGrid,stateEn,numVibStates,ignoreState)
end subroutine populateStatesVcs



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: readVcsPseudo
!Purpose: reads in the vibrational pseudostates modelling excitations 
!         to the vibrational continuum of a particular state.
!Date last modified: 03/09/2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readVcsPseudo(statebasis)
    use numbers
    use input_data         ! definitions for input data type, construct routine, object data_in
    use state_class        ! defines states (and basis of them) with operations on them
    use totalcs_module     !  reading totalcs files
    implicit none
    type(basis_state)::statebasis !Basis of states produced by populateStatesVcs
    type(basis_state)::tempbasis
    type(state), dimension(:), allocatable:: newstates
    integer:: ii, jj, vVal, newStateInd, largeN, maxNewStates, vibInd, numDissOld
    integer:: tempIndex
    character(len=141):: line
    logical:: header, cutoff, dataline, foundState
    real(dp):: currEn  !Energy of current input file
    real(dp), dimension(:), allocatable:: csArray, einArray

    largeN = 20000   !Choose a very large number for maximum possible number of new states

    allocate(newstates(largeN))

    !Use intAN excitation cross sections
    !column labelled 'threshold' contains the threshold energy of each state.
    !This threshold energy is the excitation energy. To obtain the absolute energy,
    !add the (negative) ground state energy to the threshold energy.

    !Read totalcs files for vibrational pseudostate cross sections. Each file is
    !for a particular incident particle energy.
    maxNewStates = 0
    ii = 0 
    jj = 0
    print*, "Totalcs Vibrational Pseudostates, NVcsPseudo = ", data_in%NVcsPseudo
    if (data_in%dissOn .eq. 1) then
       do ii = 1, data_in%NVcsPseudo
          open(70,file=data_in%filename_VcsPseudo(ii),action='read')
          read(data_in%filename_VcsPseudo(ii)(57:67), *) currEn       
          print*, ii, data_in%filename_VcsPseudo(ii)
 
          header = .false.
          cutoff = .false.
          dataline = .false. 
          vVal = 0
          newStateInd = 1
          do while (.not. cutoff)
             read(70,'(A)') line 
             !print*, line
             if (line(3:4) .eq. '--') then
                !Found line denoting excitation CS to new electronic state, next line contains useful data 
                if (header .eqv. .false.) then
                   !Found 'header' for a block of cross section, contains data about initial and final state
                   header = .true. 
                else if (header .eqv. .true.) then
                   !Found end of header
                   header = .false.
                   dataline = .true. 
                end if
             else if ((line(1:1) .eq. '') .and. (line(2:2) .eq. '')) then  !Reached blank line between cross section blocks
                if (dataline .eqv. .true.) then
                   dataline = .false.
                end if
             end if

             if ((header .eqv. .false.) .and. (dataline .eqv. .true.)) then
                !Read in cross sections until they become zero
                if (line(6:6) .eq. '<') then    !'<' denotes dissociative pseudostate
                   if (line(8:8) .eq. '>') then
                      !v in single digits
                      read(line(7:7),*) vVal
                   else if (line(9:9) .eq. '>') then
                      read(line(7:8),*) vVal
                   else if (line(9:9) .ne. '>') then
                      read(line(7:9),*) vVal
                   end if

                   if (ii .eq. 1) then
                      !Reading first file, create new state with given label
                      foundState = .false.
                      jj = 1
                      do while (.not. foundState)
                         if (trim(adjustl(line(1:5))) .eq. trim(adjustl(statebasis%b(jj)%stlabel))) then
                            call copy_state(newstates(newStateInd), statebasis%b(jj))
                            newstates(newStateInd)%v = vVal
                            newstates(newStateInd)%df(:) = 1.0  !All continuum vibrational pseudostates are dissociative by definition
                            allocate(csArray(data_in%NVcsPseudo),einArray(data_in%NVcsPseudo))
                            csArray(:) = 0.0d0
                            einArray(:) = 0.0d0
                            call set_cs(newstates(newStateInd),data_in%NVcsPseudo,einArray,csArray,0)
                            deallocate(csArray, einArray)

                            read(line(131:141), *) newstates(newStateInd)%enex
                            newstates(newStateInd)%en = statebasis%b(1)%en + newstates(newStateInd)%enex 
                            newstates(newStateInd)%resolved = .true.

                            foundState = .true.
                            maxNewStates = maxNewStates + 1
                         else if (jj .eq. statebasis%n) then
                            print*, "ERROR: could not find state with given electronic state label, stopping" 
                            print*, "Label: ", line(1:5)
                            stop
                         end if
                         jj = jj+1
                      end do           
                   else
                      !Find state with this label
                      jj = 1
                      foundState = .false.

                      do while(.not. foundState)
                         if ((trim(adjustl(line(1:5))) .eq. trim(adjustl(newstates(jj)%stlabel))) .and. (vVal .eq. newstates(jj)%v)) then
                            foundState = .true.
                         else if (jj .eq. maxNewStates) then
                            print*, "ERROR: could not find state with given electronic and vibrational state label, stopping" 
                            print*, "Label, v: ", line(1:5), vVal
                            stop
                         end if
                         jj = jj+1
                      end do
                      newStateInd = jj-1
                   end if

                   read(line(92:102), *) newstates(newStateInd)%cs(ii)
                   newstates(newStateInd)%ein(ii) = currEn
                   !print*, "State, v, ein, cs: ", newstates(newStateInd)%stlabel, vVal, currEn, newstates(newStateInd)%cs(ii)

                   newStateInd = newStateInd + 1    !If ii .ne. 1, then newStateInd will be reset anyway, so this covers both cases
                end if

             else if ((header .eqv. .true.) .and. (line(3:3) .ne. '-')) then
                !Number of digits in incident energy (e.g. E=1, E=20, E=100)  affects position of initial vibrational number
                if (int(currEn) .lt. 10) then
                   vibInd = 49
                else if ((int(currEn) .ge. 10) .and. (int(currEn) .lt. 100)) then
                   vibInd = 50
                else if ((int(currEn) .ge. 100) .and. (int(currEn) .lt. 1000)) then
                   vibInd = 51
                else
                   print*, "ERROR: incident energy greater than 1000eV, stopping."
                   stop
                end if
                if (line(vibInd:vibInd) .ne. '0') then
                   !Found cross sections for scattering on excited states, not using at this stage
                   cutoff = .true.
                end if
             end if 
          end do

          close(70)
       end do
    end if    

    numDissOld = 0
    do ii = 1, statebasis%n
       if (statebasis%b(ii)%v .eq. -1) then
          numDissOld = numDissOld + 1
       end if
    end do

    !Merge old and new state bases
    allocate(tempbasis%b(statebasis%n+maxNewStates-numDissOld))
    tempbasis%n = statebasis%n+maxNewStates - numDissOld
    tempIndex = 1
    do ii = 1, statebasis%n
       !Ignore previously existing dissociative states, now modelled with pseudostates.
       if (statebasis%b(ii)%v .ne. -1) then
          !If non-dissociative, copy over old state
          call copy_state(tempbasis%b(tempIndex),statebasis%b(ii))
          tempIndex = tempIndex + 1
       end if
    end do
    do ii = 1, maxNewStates 
       call copy_state(tempbasis%b(ii+statebasis%n-numDissOld),newstates(ii))
       call destruct_state(newstates(ii))
    end do
    deallocate(newstates) 
    call destruct_basis(statebasis) 
    
    call new_basis(statebasis,tempbasis%n,data_in%tcstype)
    do ii=1, tempbasis%n
       call copy_state(statebasis%b(ii), tempbasis%b(ii))
    end do 
    call destruct_basis(tempbasis) 

    !Sort all states by energy
    call sortEnergies(statebasis)

		!Used for testing purposes to switch off excitations of specific states
	!	do ii = 1, statebasis%n
	!		 !Switch off dissociative excitations
	!		 if (statebasis%b(ii)%stlabel .eq. 'C1Pu') then
	!				if( statebasis%b(ii)%enex - statebasis%b(ii)%dissThresh .lt. 0.0_dp) then
	!				   statebasis%b(ii)%cs(:) = 0.0_dp
	!			  end if
	!		 end if
	!		 if (statebasis%b(ii)%stlabel .eq. 'B1Su') then
	!				if( statebasis%b(ii)%enex - statebasis%b(ii)%dissThresh .lt. 0.0_dp) then
	!				   statebasis%b(ii)%cs(:) = 0.0_dp
	!			  end if
	!		 end if
	!  end do

end subroutine readVcsPseudo






!
!
! get CS for all states at the requested incident energy
subroutine get_csatein(statebasis,eincident,tcs)
  use numbers
  use input_data         ! definitions for input data type, construct routine, object data_in
  use state_class        ! defines states (and basis of them) with operations on them
  use totalcs_module     !  reading totalcs files
  implicit none

  type(basis_state), intent(in):: statebasis
  real(dp), intent(in):: eincident
  type(totalcs), intent(out):: tcs

  integer:: inum, j, i, j1, j2, Nmax, tcstype
  real(dp):: e1, e2, w1, w2, cs1, cs2, df1, df2
  real(dp), dimension(statebasis%n):: ar_cs, ar_df

  !print*, 'ein (getcsatein) = ', eincident 

  inum = statebasis%b(1)%inum  ! as they are the same for all states
  e1 = statebasis%b(1)%ein(1)
  e2 = statebasis%b(1)%ein(inum)
  if(eincident .gt. e2) then
     print*, e2
     print*, "ERROR: incident energy greater than maximum energy, stopping"
     print*, "get_csatein"
     stop
  else if (eincident .lt. e1) then
     ! set all cross sections to zero
     Nmax = 1  !Only elastic scattering possible
  else
     j = 1
     do i= 1,inum      !   tests to be put on not exiting ....
        !print*, i, statebasis%b(1)%ein(i) 
        !print*, eincident, statebasis%b(1)%ein(i), i
        if(eincident .lt. statebasis%b(1)%ein(i)) exit
        j = j + 1
     enddo
     j1 = j - 1
     j2 = j
     !print*, 'eincident = ', eincident
     !print*, 'j1, j2:', j1, j2
     !print*, 'statebasis%b(1)%ein(j2) = ', statebasis%b(1)%ein(j2)    
    
     e1 = statebasis%b(1)%ein(j1)
     e2 = statebasis%b(1)%ein(j2)
     w1 = (eincident - e1)/(e2 - e1) 
     w2 = (e2 - eincident)/(e2 - e1)

     Nmax = 0
     do i=1,statebasis%n
        !print*,'i, eincident, statebasis%b(i)%enex', i, eincident, statebasis%b(i)%enex
        if(eincident .lt. statebasis%b(i)%enex) exit
        Nmax = i 
        cs1 = statebasis%b(i)%cs(j1)    
        cs2 = statebasis%b(i)%cs(j2)    
        ar_cs(i) = w1*cs1 + w2*cs2


!    diss. fractions are 0 for all positive energy states
!                        1 for all negative energy triplet states
!                        1 for high-lying singlet states
!                        to be determined form diss. fraction files for states in the list
        if(statebasis%b(i)%en .ge. 0) then
           ar_df(i) = 0d0
        else
           if(statebasis%b(i)%spin .gt. 0) then ! spin=1, triplet states ony for H2 !!! might need a fix for a general case
              ar_df(i) = 1d0
           else   ! spin=0, singlet states for H2
              ar_df(i) = 1d0
              if( statebasis%b(i)%idf .ne. 0) then
                 df1 = statebasis%b(i)%df(j1)    
                 df2 = statebasis%b(i)%df(j2)
                 ar_df(i) = w1*df1 + w2*df2    
              endif
              if(i .eq. 1) then
                 ar_df(i) = 0d0
              end if
           endif
        endif
     enddo
     if(Nmax .eq. 0) then
        print*,'tcstostate.f90: get_csatein(): Nmax=', 0
        print*, 'stopping'
        stop
     endif
     !print*, 'Nmax= ', Nmax
  end if

  tcstype = statebasis%tcstype
!  make an instance of tcs object: 
  call new_totalcs(tcs,Nmax,tcstype)

  tcs%en_incident = eincident
  tcs%tcstype = statebasis%tcstype
  do i=1,Nmax
     tcs%nf(i) = i
     tcs%ni(i) = 1  ! change later on to a more general case?
     tcs%CS(i) =  ar_cs(i) ! update later on... for intial and final states
     tcs%en(i) = statebasis%b(i)%en/data_in%eV  ! ionization energy in a.u.
     tcs%df(i) = ar_df(i)
  enddo

end subroutine get_csatein


