!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!File: vcsOld.f90
!Purpose: contains subroutines for the now redundant 'vcs' type. 
!Author: Reese Horton, Curtin University, ID: 19456155
!Date last modified: 10/02/2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module vcsOld


  !Stores cross sections for a particular vibrational transition as a function 
  !of incident particle energy. Mainly used for file IO. 
  type, public::vcs

     character*40:: filename   !Name of input file, used for file IO only.
     real(dp):: en_incident
     integer::Nein !Number of energy grid points
     integer:: nvf  !number of final vibrational states
     real(dp), pointer, dimension(:) :: en => NULL()   !excitation en array, also used for file IO.  
     integer, pointer, dimension(:):: vf !Final vibrational state in vi =0 -> vf transition, file IO
     real(dp), pointer, dimension(:) :: cs => NULL()   !excitation cross section
  end type

  !Stores a vcs type for vibrational transitions at different incident energies.
  type, public::basis_vcs
     type(vcs), pointer, dimension(:):: b => NULL()
     real(dp), pointer, dimension(:):: ein => NULL()
     integer:: numVcs
  end type



contains



  !Subroutines relating to the vibrational cross sections and the various available
  !states.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine:readVcs 
  !Purpose: reads in vibrational excitation cross sections from input 
  !         files
  !Date last modified: 18/12/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readVcs(vcsIn,filename)
      implicit none
      type(vcs)::vcsIn
      character*40, intent(in):: filename
      logical:: ex
      integer:: nfile, i, Nmax, vfVal
      real(dp)::einVal, csVal

      ! open vcs file: filename
      inquire(file=filename,exist=ex)
      if (ex) then
         nfile = 100 
      else
         print*, filename, 'is not found, stopping'
         stop
      endif
      open(nfile,file=filename)
    
      vcsIn%filename = filename

      !Ignore first two lines 
      read(nfile,*)
      read(nfile,*)
      read(nfile,*)   vfVal !Final vibrational state
      read(nfile,*)
      allocate(vcsIn%vf(1))
      vcsIn%vf(1) = vfVal
      i = 0
      do
         read(nfile,*, ERR=110,END=110) einVal, csVal
         i = i + 1
      enddo
110   continue
      Nmax = i
 
      vcsIn%Nein = i
      allocate(vcsIn%en(i),vcsIn%cs(i)) 
      rewind(nfile)
      read(nfile,*)
      read(nfile,*)
      read(nfile,*)
      read(nfile,*)
      do i=1,Nmax
         read(nfile,*) einVal, csVal 
         vcsIn%en(i) = einVal
         vcsIn%cs(i) = csVal
      enddo

      close(nfile)
  end subroutine readVcs



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: populateVcs
   !Purpose: fills vcsBasis with vcs at various energies 
   !         on the fine energy grid of encident energies.
   !         End results is a set of vcs types, usable for 
   !         sampling vibrational excitations at a given 
   !         incident energy.
   !Date last modified: 18/12/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine populateVcs(vcsBasis,enfine,Nenfine,data_input)
       use input_data
       implicit none 
       type(input)::data_input
       type(basis_vcs)::vcsBasis
       type(vcs), allocatable, dimension(:)::tempBasis
       integer:: ii, jj, kk, Nenfine, numBelow, lenVal
       integer::halfInt
       real(dp), allocatable, dimension(:),intent(in)::enfine
       real(dp), allocatable, dimension(:)::intpCsVals, newGrid, excEnArray
       real(dp), allocatable, dimension(:)::tempEnArray, tempCsArray
       real(dp)::currEn, deltaEMin

 
       !Construct array of excitation energies of vibrational states
       allocate(excEnArray(data_input%Nvcs))
       do ii = 2, data_input%Nvcs
          !ii =1 corresponds to vf = 0 in vcsBasis
          excEnArray(ii) = vcsBasis%b(ii)%en(1)
       end do
       excEnArray(1) = 0.0 !vi=0 -> vf=0 transition has zero excitation energy

       !Reallocate arrays from input files to account for the first value in 
       !each file for vf > 0 being the excitation energy.
       do ii = 2, data_input%Nvcs
          lenVal = size(vcsBasis%b(ii)%en)
          allocate(tempEnArray(lenVal),tempCsArray(lenVal))
          tempEnArray(:) = vcsBasis%b(ii)%en(:)
          tempCsArray(:) = vcsBasis%b(ii)%cs(:)
          deallocate(vcsBasis%b(ii)%en,vcsBasis%b(ii)%cs)
          allocate(vcsBasis%b(ii)%en(lenVal-1),vcsBasis%b(ii)%cs(lenVal-1))
          vcsBasis%b(ii)%en(:) = tempEnArray(2:lenVal)
          vcsBasis%b(ii)%cs(:) = tempCSArray(2:lenVal)
          deallocate(tempEnArray,tempCsArray)
       end do

       !Construct new energy grid in the domain [0,enfine(max)] eV
       deltaEMin = vcsBasis%b(1)%en(2) - vcsBasis%b(1)%en(1)
       numBelow = int(enfine(1)/deltaEMin)
       allocate(newGrid(Nenfine+numBelow)) 
       newGrid(1) = vcsBasis%b(1)%en(1)

       do ii = 2, numBelow
          newGrid(ii) = newGrid(ii-1) + deltaEMin
       end do
       newGrid(numBelow+1:Nenfine+numBelow) = enfine(:) 

       allocate(tempBasis(Nenfine+numBelow))
  
       do ii = 1, Nenfine+numBelow
          currEn = newGrid(ii)
          allocate(tempBasis(ii)%vf(data_input%Nvcs)) !Number of vibrational states = # of input files
          allocate(tempBasis(ii)%cs(data_input%Nvcs))
          tempBasis(ii)%en_incident = currEn
       end do

       allocate(intpCsVals(Nenfine+numBelow))
       intpCsVals(:) = 0.0
       do jj = 1, data_input%Nvcs
          !Interpolate to find excitation cross sections for each vf
          halfInt = int(real(Nenfine+numBelow)/2.0)
          call INTRPL(size(vcsBasis%b(jj)%en),vcsBasis%b(jj)%en,vcsBasis%b(jj)%cs,halfInt,newGrid(1:halfInt),intpCsVals(1:halfInt))
          call INTRPL(size(vcsBasis%b(jj)%en),vcsBasis%b(jj)%en,vcsBasis%b(jj)%cs,Nenfine+numBelow-halfInt,newGrid(halfInt+1:Nenfine+numBelow),intpCsVals(halfInt+1:Nenfine+numBelow))
          tempBasis(:)%vf(jj) = jj-1              
          tempBasis(:)%cs(jj) = intpCsVals(:)
          do kk = 1, Nenfine+numBelow
             if (newGrid(kk) .gt. vcsBasis%b(jj)%en(size(vcsBasis%b(jj)%en))) then
                tempBasis(kk)%cs(jj) = 0.0
             else if (newGrid(kk) .lt. vcsBasis%b(jj)%en(1)) then
                tempBasis(kk)%cs(jj) = 0.0
             end if
          end do
          intpCsVals(:) = 0.0
       end do     

       deallocate(intpCsVals)
       call delete_vcsbasis(vcsBasis)

       !Overwrite vcsBasis 
       vcsBasis%numVcs = Nenfine+numBelow
       allocate(vcsBasis%b(vcsBasis%numVcs),vcsBasis%ein(numBelow+Nenfine))
       vcsBasis%ein(:) = newGrid(:)
       do ii = 1, vcsBasis%numVcs
          allocate(vcsBasis%b(ii)%cs(data_input%Nvcs))
          allocate(vcsBasis%b(ii)%vf(data_input%Nvcs))
          allocate(vcsBasis%b(ii)%en(data_input%Nvcs))
          vcsBasis%b(ii)%vf(:) = tempBasis(ii)%vf(:)
          vcsBasis%b(ii)%cs(:) = tempBasis(ii)%cs(:) 
          vcsBasis%b(ii)%en(:) = excEnArray(:)
          vcsBasis%b(ii)%en_incident = newGrid(ii)
          !Same number of final vibrational states for all energies
          vcsBasis%b(ii)%nvf = data_input%Nvcs
       end do
       do ii = 1, Nenfine+numBelow
          call delete_vcs(tempBasis(ii))
       end do
       deallocate(tempBasis,newGrid,excEnArray)
   end subroutine populateVcs




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: newVcsBasis
   !Purpose: constructs a new basis_vcs type, allocating 
   !         memory for data structure
   !Date last modified: 18/12/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine newVcsBasis(vcsBasis,nVcsIn)
       implicit none
       type(basis_vcs)::vcsBasis
       integer:: nVcsIn

       if (associated(vcsBasis%b)) then
          deallocate(vcsBasis%b)
       end if
       vcsBasis%numVcs = nVcsIn
       allocate(vcsBasis%b(nVcsIn))
   end subroutine newVcsBasis



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: delete_vcs
   !Purpose: deallocates memorey for the imported vcs data type
   !Date last modified: 18/12/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine delete_vcs(vcsIn)
       implicit none
       type(vcs)::vcsIn

       if (associated(vcsIn%en)) then
          deallocate(vcsIn%en)
       end if
       if (associated(vcsIn%cs)) then
          deallocate(vcsIn%cs)
       end if
       if (associated(vcsIn%vf)) then
          deallocate(vcsIn%vf)
       end if
   end subroutine delete_vcs


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: delete_vcsbasis
   !Purpose: deallocates the imported basis_vcs type
   !Date last modified: 18/12/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine delete_vcsbasis(vcsBasis)
       implicit none
       type(basis_vcs)::vcsBasis
       integer::ii       

       if (associated(vcsBasis%b)) then
          do ii = 1, vcsBasis%numVcs
             call delete_vcs(vcsBasis%b(ii))
          end do
       end if
       if (associated(vcsBasis%ein)) then
          deallocate(vcsBasis%ein)
       end if
   end subroutine delete_vcsbasis


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: get_vcsatein
   !Purpose: constructs a vcs type at the requested incident
   !         energy using the imported basis_vcs type.
   !Date last modified: 18/12/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_vcsatein(vcsBasis, eIncident, vcsIn)
       implicit none
       type(basis_vcs)::vcsBasis
       type(vcs)::vcsIn, vcsBelow, vcsAbove
       real(dp)::eIncident, weight1, weight2
       integer:: ii, below, after
       logical:: found

       if (associated(vcsIn%cs)) then
          deallocate(vcsIn%cs)
       end if
       if (associated(vcsIn%vf)) then
          deallocate(vcsIn%vf)
       end if
       if (associated(vcsIn%en)) then
          deallocate(vcsIn%en)
       end if 

       !Find vcs types at energies adjacent to eIncident
       below = 1
       do while ((vcsBasis%b(below)%en_incident .lt. eIncident) .and. (below .lt. vcsBasis%numVcs))
          below = below +1
       end do
   
   
       found = .false.
       if ((below .ge. 2) .and. (below .le. vcsBasis%numVcs)) then
          after = below
          below = below -1
          found = .true.
       else if (below .lt. 2) then
          !Incident energy is below lowest available.
          found = .false.
       end if
      
       if (found) then
          vcsIn%nvf = vcsBasis%b(after)%nvf
          allocate(vcsIn%vf(vcsBasis%b(after)%nvf),vcsIn%cs(vcsBasis%b(after)%nvf),vcsIn%en(vcsBasis%b(after)%nvf)) 
          vcsAbove = vcsBasis%b(after)
          vcsBelow = vcsBasis%b(below)
          !Interpolate cross sections
          weight1 = (eIncident - vcsBelow%en_incident)/(vcsAbove%en_incident-vcsBelow%en_incident)
          weight2 = (vcsAbove%en_incident-eIncident)/(vcsAbove%en_incident-vcsBelow%en_incident)
          vcsIn%cs(:) = weight1*vcsBelow%cs(:) + weight2*vcsAbove%cs(:) 
          vcsIn%en(:) = vcsBelow%en(:) !Excitation energies are always the same.
          vcsIn%en_incident = eIncident
       else 
          !Energy requested is invalid, either below 0.01 eV cutoff or above maximum
          print*, "ERROR: vcs incident energy below minimum, stopping"
          stop
       end if

   end subroutine get_vcsatein


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: printVcsAtEn
   !Purpose: prints vibrational excitation cross sections
   !         at given incident energy to a file for plotting
   !Date last modified: 18/12/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine printVcsAtEn(vcsaten)
       implicit none
       type(vcs)::vcsaten
       integer:: ii       
       character (len=3)::enString
 
       write(enString,'(I0)') int(vcsaten%en_incident)
       open(70, file="vcs"//trim(enString)//".txt")
       write(70,*) "Energy(eV)          CS(a.u)"
       do ii = 1, vcsaten%nvf
         write(70,*) vcsaten%vf(ii), vcsaten%CS(ii)
       end do
       close(70)    

   end subroutine printVcsAtEn




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: printVcsState
   !Purpose: prints vibrational excitation cross sections
   !         for the given transition over the range of
   !         available energies.
   !Date last modified: 18/12/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine printVcsState(vcsBasis, vfIn)
       implicit none
       type(basis_vcs)::vcsBasis
       integer::vfIn !Final vibrational state
       type(vcs)::vcsCurr
       integer::ii
       character (len=1)::stateString
 
       write(stateString,'(I0)') int(vfIn)
       open(70, file='vcsvf'//trim(stateString)//'.txt')
       write(70, *) "Energy (eV)            CS (a.u)"
       do ii = 1, vcsBasis%numVcs
          vcsCurr = vcsBasis%b(ii)
          write(70,*) vcsCurr%en_incident,  vcsCurr%cs(vfIn+1)
       end do
       close(70)
   end subroutine printVcsState











end module vcsOld
