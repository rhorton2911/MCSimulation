module totalcs_module

  use numbers

  type, public:: totalcs
     ! private
     integer:: tcstype
     character(len=40):: filename
     real(dp):: en_incident  ! energy 
     integer:: Nmax  ! number of states
     integer, pointer, dimension(:) :: nf => NULL(), ni => NULL()
     real(dp), pointer, dimension(:) :: BICS => NULL(), CS => NULL(), en => NULL(), df => NULL()
  end type totalcs


  type, public:: basis_totalcs
     !     private
     type(totalcs), pointer, dimension(:) :: b  => NULL()
     integer:: n
     integer:: tcstype
  end type basis_totalcs
  ! 


  !Stores cross sections and dissociation fractions 
  !for dissociative ionisation. 
  !Cross sections converted to a.u in preprocessing.
  type, public:: basis_dics
     integer:: nDics !Number of array elements for below
     real(dp), pointer, dimension(:):: enIncident => NULL() !Incident energy (eV)
     real(dp), pointer, dimension(:):: dics => NULL() !cross sections (a_0^2)
     real(dp), pointer, dimension(:):: dfIon => NULL()!Dissociation fractions
  end type


contains


!Dissociative ionisation cross section functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: readDics
   !Purpose: reads in dissociative ionisation cross sections
   !         from a file.
   !Date last modified: 28/03/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine readDics(fileName, dicsBasis)
       implicit none
       type(basis_dics), intent(inout):: dicsBasis
       character(len=26), intent(in)::fileName  
       logical::exists
       integer:: fileNum, ii, nMax
       real(dp):: en, dicsVal, bohrSquare
 
       bohrSquare = 27.9841E-22

       !Checks to see if file exists and opens it   
       inquire(file=fileName,exist=exists)
       if (exists) then
          fileNum = 20 
       else
          print*, fileName, 'is not found, stopping'
          stop
       endif
       open(fileNum, file = fileName,action='read')

       read(fileNum, *) !Ignores first line 
      
       !Counts number of lines in file
       ii = 0
       do
          read(fileNum,*, ERR=105, END=105) en, dicsVal
          ii = ii+1
       end do
105 continue        

       nMax = ii

       allocate(dicsBasis%enIncident(nMax), dicsBasis%dics(nMax))
       rewind(fileNum)  
       read(fileNum,*)

       do ii=1, nMax
          read(fileNum, *) en, dicsVal
          dicsBasis%enIncident(ii) = en !Energy in (eV)
          dicsBasis%dics(ii) = dicsVal/bohrSquare !Convert to a.u
       end do
              
       close(fileNum)
   end subroutine readDics



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: fillDicsBasis
   !Purpose: calculates dissociation fractions for dissociative
   !         ionisation and stores them in basis_dics. Also
   !         reallocates energy array and cross section array
   !         to match fine energy grid
   !Date last modified: 28/03/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fillDicsBasis(stateBasis,basisToFill,Nenfine,enfine)
       use state_class
       implicit none
       type(basis_state), intent(in)::stateBasis
       type(basis_dics), intent(inout)::basisToFill
       type(totalcs)::tcsAtEn
       integer:: Nenfine, ii
       real(dp), allocatable, dimension(:)::enfine !Fine energy grid
       real(dp), allocatable, dimension(:)::dicsArray
       real(dp):: dfIonVal !Dissociation fraction
       real(dp):: tics !Total ionisatino cross section

       if (associated(basisToFill%dfIon)) then
          deallocate(basisToFill%dfIon)
       end if
       
       allocate(dicsArray(Nenfine), basisToFill%dfIon(Nenfine) )

       !Interpolates dics to fine energy grid.
       call intpDics(dicsArray, basisToFill, Nenfine, enfine)
    
       !Reallocate to match fine energy grid 
       deallocate(basisToFill%dics, basisToFill%enIncident)
       allocate(basisToFill%dics(Nenfine), basisToFill%enIncident(Nenfine))

       !Calculates dissociation fractions and stores values in basis
       !Values for which df are defined lie within bounds of energy grid.
       !As per Dmitry: 
          !-below first point on old grid make all values in new grid
          ! equal to that of the first point on the old grid.
       do ii=1, Nenfine-1
          call get_csatein(stateBasis,enfine(ii),tcsAtEn)
          tics = 0
          call getTics(tcsAtEn,tics)
          if( tics .gt. 0.d0 ) then
             dfIonVal = dicsArray(ii)/tics
          else
             !energy less than required to cause ionisation
             dfIonVal = 0.d0
          end if
          basisToFill%dfIon(ii) = dfIonVal
          basisToFill%dics(ii) = dicsArray(ii)
          basisToFill%enIncident(ii) = enfine(ii)
          call delete_totalcs(tcsAtEn)
       end do
       basisToFill%nDics = Nenfine 

       !Deallocate memory
       deallocate(dicsArray) 
   end subroutine fillDicsBasis



   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: intpDics
   !Purpose: linearly interpolates dissociative ionisation cross sections 
   !         to fine energy grid. Sets dics at energies less than ionisation 
   !         potential (smallest energy from input file) to zero.
   !Date last modified: 28/03/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine intpDics(dicsArray, basisDics, Nenfine, enfine) 
       implicit none
       real(dp), allocatable, dimension(:)::enfine, dicsArray
       type(basis_dics), intent(in)::basisDics
       integer::Nenfine, ii, jj
       real(dp):: energy, dicsAbove, dicsBelow
       real(dp):: intp1, intp2, intp3 !Breaks up long interpolation formula

    
       !Note: dissociative ionisation CS are available for energies 
       !      higher than upper bound on fine energy grid, thus go to this.
       do ii=1, Nenfine
          energy = enfine(ii)

          if( energy .gt. basisDics%enIncident(1) ) then
             !Finds dics above and below given energy in grid
             jj=1
             do while(basisDics%enIncident(jj) .lt. energy)
                jj = jj +1
             end do
             dicsBelow = basisDics%dics(jj-1)
             dicsAbove = basisDics%dics(jj)
             !Performs linear interpolation to grid energy
             intp1 = enfine(ii) - basisDics%enIncident(jj-1) 
             intp2 = dicsAbove - dicsBelow
             intp3 = basisDics%enIncident(jj) - basisDics%enIncident(jj-1)
             dicsArray(ii) = dicsBelow + intp1*(intp2/intp3)
          else 
             !Point is below ionisation energy, cross section is zero
             dicsArray(ii) = 0.0
          end if
       end do

   end subroutine intpDics



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: getDissIonDf
   !Purpose: interpolates to find the dissociative ionisation
   !         dissociation fraction at the given incident
   !         energy.
   !Date last modified:28/03/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine getDissIonDf(dicsBasis,enIn,df)
       implicit none
       real(dp)::df, enIn
       type(basis_dics), intent(in)::dicsBasis
       real(dp)::intp1, intp2, intp3 !Breaks up long formula
       integer::ii
       real(dp)::dicsBelow, dicsAbove              

       ii=1
       do while (dicsBasis%enIncident(ii) .lt. enIn)
          ii = ii +1
       end do
       dicsBelow = dicsBasis%dics(ii-1)
       dicsAbove = dicsBasis%dics(ii)
 
       intp1 = enIn - dicsBasis%enIncident(ii-1)
       intp2 = dicsAbove-dicsBelow
       intp3 = dicsBasis%enIncident(ii) - dicsBasis%enIncident(ii-1)
       df = dicsBelow + intp1*(intp2/intp3)
   end subroutine getDissIonDf




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: printDics
   !Purpose: prints interpolated dissociative ionisation 
   !         cross sections to file to check for consistency
   !Date last modified: 18/08/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine printDics(dicsBasis)
       implicit none
       type(basis_dics)::dicsBasis
       integer:: ii

       open(10,file="dicsTest.txt")
       do ii = 1, dicsBasis%nDics-1
          write(10,*,err=12) dicsBasis%enIncident(ii), dicsBasis%dics(ii)
       end do
12     close(10)
   end subroutine printDics




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: destructDicsBasis
   !Purpose: deallocates the contents of a basis_dics type
   !Date last modified: 28/03/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine destructDicsBasis(basisDics)
       implicit none
       type(basis_dics)::basisDics

       if (associated(basisDics%enIncident)) then
          deallocate(basisDics%enIncident)
       end if
       if (associated(basisDics%dics)) then
          deallocate(basisDics%dics)
       end if
       if (associated(basisDics%dfIon)) then
          deallocate(basisDics%dfIon)
       end if
   end subroutine destructDicsBasis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
  subroutine read_totalcs(self,tcstype,filename)
    implicit none
    type(totalcs), intent(inout):: self
    integer, intent(in):: tcstype
    character(len=40), intent(in):: filename
    logical:: ex
    integer:: nfile, i, Nmax, jf, ji
    character(len=3):: str
    real(dp)::  BICS, pwBICS, pwICS, extICS, en
    character(len=20)::  str1, str2, str3, str4

    self%tcstype = tcstype

    ! open totalcs file: filename
    inquire(file=filename,exist=ex)
    if (ex) then
       nfile = 100 
!       print*, 'found: ', filename
    else
       print*, filename, 'is not found, stoping'
       stop
    endif
    open(nfile,file=filename,action='read')
    
    self%filename = filename
!    nch = LEN(filename)
!   nchi = SCAN(filename,'_') + 1
!   stren = filename(nchi:nch)
!   read(stren,*), self%en_incident  !  fix it from the filename
!    print*, '--> energy = ', stren,self%en_incident 

    
    read(nfile,*)
    read(nfile,*)
    read(nfile,*)   str1, str2, str3, str4, en
!    print*, 'en=', en
    self%en_incident = en
    read(nfile,*)
    i = 0
    do
       read(nfile,*, ERR=110,END=110) jf,str,ji,BICS,pwBICS,pwICS,extICS,en
       i = i + 1
    enddo
110 continue
    Nmax = i

        
    self%Nmax = i
    allocate(self%nf(Nmax),self%ni(Nmax),self%BICS(Nmax),self%CS(Nmax),self%en(Nmax),self%df(Nmax))
    
    rewind(nfile)
    read(nfile,*)
    read(nfile,*)
    read(nfile,*)
    read(nfile,*)
    do i=1,Nmax
       read(nfile,*) jf,str,ji,BICS,pwBICS,pwICS,extICS,en
       self%nf(i) = jf
       self%ni(i) = ji
       self%BICS(i) = BICS
       self%CS(i) = extICS
       self%en(i) = en
       self%df(i) = 0  ! just to set it to a value
    enddo

    close(nfile)
    
  end subroutine read_totalcs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: readPsTcs
  !Purpose: reads in positron scattering electron excitation cross 
  !         section file, columns are state excitation cross sections
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readPsCs(tcsbasis,tcstype,filename)
    implicit none
    type(basis_totalcs), intent(inout):: tcsbasis
    integer, intent(in):: tcstype
    character(len=40), intent(in):: filename
    integer:: ii, lines=0, IERR=0, numcols
    character(len=50) :: a
    real(dp), dimension(:), allocatable:: scalefact
    real(dp):: en

    open(60,filename)
    read(60, *) numcols
    read(60, *)
    do while (IERR == 0)
      read(60, *, iostat=IERR) a
      lines = lines + 1
    end do
    lines = lines - 1
    IERR = 0
    rewind(60)

    !Read in scale factor for multiplying cross sections (2nd column)
    allocate(scalefact(lines))
    read(60,*) 
    read(60,*) 
    do ii = 1, lines
       read(60,*) en, scalefact(ii)
    end do

    rewind(60)

    call new_totalcsbasis(tcsbasis,lines,tcstype)
    do ii = 1, lines
       call new_totalcs(tcsbasis%b(ii), numcols, tcstype)
    end do

    read(60,*)
    do ii = 1, lines
       read(60,*) tcsbasis%b(ii)%energy, (tcsbasis%b(ii)%cs(jj), jj=1, numcols)
       tcsbasis%b(ii)%cs(:) = tcsbasis%b(ii)%cs(:)*scalefact(:) 
    end do

    close(60)

  end subroutine readPsCs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: readPsFormCs
  !Purpose: reads in positronium formation cross sections into
  !         the input state type
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readPsFormCs(self,filename)
    implicit none
    type(state):: self
    character(len=40), intent(in):: filename
    integer:: ii, lines=0, IERR=0, numcols
    character(len=50) :: a

    open(60,filename)
    read(60, *) 
    read(60, *) 
    read(60, *)
    do while (IERR == 0)
      read(60, *, iostat=IERR) a
      lines = lines + 1
    end do
    lines = lines - 1
    IERR = 0
    rewind(60)

    eneV = -1.0_dp
    enex = -1.0_dp

    call init_state(self,"PsFormation",eneV,enex,-1,-1,-1,-1,1)
    self%psFormation = .true.
    
    do ii = 1, lines
       read(60,*) self%en(ii), self%cs(ii)
    end do
    
    close(60)

   end subroutine readPsFormCs



  
 !
  subroutine new_totalcs(self,Nmax,tcstype)
    implicit none
    type(totalcs), intent(inout):: self
    integer, intent(in):: Nmax, tcstype

    self%Nmax = Nmax
    self%tcstype = tcstype
    ! create array of n states
    if(Nmax .ne. 0)  then
       allocate( self%nf(Nmax),self%ni(Nmax),self%BICS(Nmax),self%CS(Nmax),self%en(Nmax),self%df(Nmax))
    endif

  end subroutine new_totalcs
!
  subroutine delete_totalcs(self)
    implicit none
    type(totalcs), intent(inout):: self

    if(associated(self%CS)) then
       deallocate( self%nf,self%ni,self%BICS,self%CS,self%en,self%df)
    endif

  end subroutine delete_totalcs
 !
  subroutine new_totalcsbasis(self,n,tcstype)
    implicit none
    type(basis_totalcs), intent(inout):: self
    integer, intent(in):: n, tcstype

    self%n = n
    self%tcstype = tcstype
    ! create array of n states
    if(n .ne. 0)  then
       allocate( self%b(n))
    endif

  end subroutine new_totalcsbasis
  !
  subroutine destruct_totalcsbasis(self)
    implicit none
    type(basis_totalcs), intent(inout):: self
    integer:: n
    integer:: i
    integer:: stat

    if (associated(self%b)) then
       n= self%n
       do i=1,n
          !       print*,'dealocate: i=', i, ', size=', SIZE(self%b(i)%ein)
          if (associated(self%b(i)%nf)) deallocate(self%b(i)%nf)
          if (associated(self%b(i)%ni)) deallocate(self%b(i)%ni)
          if (associated(self%b(i)%BICS)) deallocate(self%b(i)%BICS)
          if (associated(self%b(i)%CS)) deallocate(self%b(i)%CS)
          if (associated(self%b(i)%en)) deallocate(self%b(i)%en)
          if (associated(self%b(i)%en)) deallocate(self%b(i)%df)
       enddo
       deallocate(self%b, STAT=stat)
    end if

  end subroutine destruct_totalcsbasis
!
 subroutine copy_tcs(tcs_l,tcs_r)
    implicit none
    type(totalcs), intent(out):: tcs_l
    type(totalcs), intent(in):: tcs_r
    integer:: i2

    tcs_l%tcstype = tcs_r%tcstype
    tcs_l%en_incident = tcs_r%en_incident
    tcs_l%Nmax = tcs_r%Nmax

    i2 = tcs_r%Nmax
    if ( associated(tcs_l%nf) ) then
       deallocate(tcs_l%nf)
       deallocate(tcs_l%ni)
       deallocate(tcs_l%BICS)
       deallocate(tcs_l%CS)
       deallocate(tcs_l%en)
       deallocate(tcs_l%df)
    endif
    allocate( tcs_l%nf(1:i2), tcs_l%ni(1:i2), tcs_l%BICS(1:i2), tcs_l%CS(1:i2) , tcs_l%en(1:i2), tcs_l%df(1:i2))
    tcs_l%nf(1:i2) = tcs_r%nf(1:i2)
    tcs_l%ni(1:i2) = tcs_r%ni(1:i2)
    tcs_l%BICS(1:i2) = tcs_r%BICS(1:i2)
    tcs_l%CS(1:i2) = tcs_r%CS(1:i2)
    tcs_l%en(1:i2) = tcs_r%en(1:i2)
    tcs_l%df(1:i2) = tcs_r%df(1:i2)
  end subroutine copy_tcs
!
  subroutine print_tcs(tcs)
    implicit none
    type(totalcs), intent(in):: tcs
    integer:: i
    
    print*, tcs%tcstype
    print*, tcs%en_incident
    print*, tcs%Nmax
    print*,'  Nf   Ni      en               CS             df'
      do i=1,tcs%Nmax
       print'(2i5,1P,5E15.5)', tcs%nf(i), tcs%ni(i), tcs%en(i), tcs%CS(i), tcs%df(i)       
    enddo
  end subroutine print_tcs
!
!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: write_tcs
  !Purpose: writes tcs to file for plotting
  !Date last modified: 31/03/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine write_tcs(tcs)
      implicit none
      type(totalcs)::tcs
      integer:: ii
      character(len=20)::filename, enString
      
      !Write incident energy to a string for use in filename
      write(enString,'(I0)') int(tcs%en_incident)

      filename = 'tcs_'//TRIM(enString)//'eV.txt'
      open(70, file=filename)
      write(70, *) "energy(eV)        CS(a.u)"
      do ii = 1, tcs%Nmax
         write(70,*) tcs%en(ii),  tcs%CS(ii)
      end do     
      close(70)

  end subroutine write_tcs



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: writeExcitationCs
  !Purpose: writes total excitation CS to file as a function
  !         of energy for all energies from 0 to 500 eV.
  !Date last modified: 25/06/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeExcitationCs(stateBasis)
      use state_class
      implicit none
      type(basis_state):: stateBasis
      integer:: ii, jj, numPoints
      type(totalcs):: tcs
      real(dp), dimension(:), allocatable:: energies
      real(dp):: deltaE, totExcCs
      logical:: foundMax       

      deltaE = 2.50_dp
      numPoints = int(500.0_dp/deltaE)+1
      allocate(energies(numPoints))
      energies(1) = deltaE
      do ii= 2, numPoints
         energies(ii) = energies(ii-1) + deltaE
      end do

      open(70,file='totExcitationCs.txt')
      write(70,*) 'energy(eV)       cs(a.u)'
      do ii = 1, numPoints
         call get_csatein(stateBasis,energies(ii), tcs)
         jj = 2  !Ignore elastic CS
         totExcCs = 0.0_dp
         foundMax = .false.
         if ((jj .eq. tcs%Nmax) .or. (energies(ii) .lt. stateBasis%b(1)%ein(1))) then
            !second condition covers case where energies(ii) is less than minimum excitation energy
            foundMax = .true.
         end if

         do while (.not. foundMax)
            if (tcs%en(jj) .lt. 0.0_dp) then 
               totExcCs = totExcCs + tcs%cs(jj)
               jj = jj + 1
            else if (tcs%en(jj) .ge. 0.0_dp) then
               foundMax = .true.
            end if

            if (jj .eq. tcs%Nmax) then
               foundMax = .true.
            end if
         end do
         write(70,*) energies(ii), totExcCs
         call delete_totalcs(tcs) 
      end do
      close(70)
  end subroutine writeExcitationCs



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: writeFullCs
  !Purpose: writes the total interaction cross section to file
  !         as a function of incident energy.
  !Date last modified: 25/06/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeFullCs(stateBasis)
      use state_class
      implicit none
      type(basis_state):: stateBasis
      integer:: ii, jj, numPoints
      type(totalcs):: tcs
      real(dp), dimension(:), allocatable:: energies
      real(dp):: deltaE, totCs
       
      deltaE = 2.50d0
      numPoints = int(500.0d0/deltaE)+1
      allocate(energies(numPoints))
      energies(1) = deltaE
      do ii= 2, numPoints
         energies(ii) = energies(ii-1) + deltaE
      end do

      open(70,file='totInteractionCs.txt')
      write(70,*) 'energy(eV)       cs(a.u)'
      do ii = 1, numPoints
         call get_csatein(stateBasis,energies(ii), tcs)
         jj = 1  !Ignore elastic CS
         totCs = 0.0d0
         do while (jj .le. tcs%Nmax) 
            totCs = totCs + tcs%cs(jj)
            jj = jj + 1
         end do
         write(70,*) energies(ii), totCs
         call delete_totalcs(tcs) 
      end do
      close(70)
  end subroutine writeFullCs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: writeElCs
  !Purpose: writes the elastic cross section to file
  !         as a function of incident energy.
  !Date last modified: 25/06/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeElCs(stateBasis)
      use state_class
      implicit none
      type(basis_state):: stateBasis
      integer:: ii, numPoints
      type(totalcs):: tcs
      real(dp), dimension(:), allocatable:: energies
      real(dp):: deltaE
       
      deltaE = 2.50d0
      numPoints = int(500.0d0/deltaE)+1
      allocate(energies(numPoints))
      energies(1) = deltaE
      do ii= 2, numPoints
         energies(ii) = energies(ii-1) + deltaE
      end do

      open(70,file='totElCs.txt')
      write(70,*) 'energy(eV)       cs(a.u)'
      do ii = 1, numPoints
         call get_csatein(stateBasis,energies(ii), tcs)
         write(70,*) energies(ii), tcs%cs(1)
         call delete_totalcs(tcs) 
      end do
      close(70)
  end subroutine writeElCs




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: writeElCsInState
  !Purpose: writes the elastic cross section to file
  !         as a function of incident energy.
  !Date last modified: 25/06/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeElCsInState(stateBasis)
      use state_class
      implicit none
      type(basis_state):: stateBasis
      integer:: ii
       
      open(70,file='elCsStates.txt')
      write(70,*)  "Energy (eV)         CS(a.u)"
      do ii= 1, stateBasis%b(1)%inum
         write(70,*) stateBasis%b(1)%ein(ii), stateBasis%b(1)%cs(ii)
      end do
      close(70)
  end subroutine writeElCsInState



  subroutine set_cs_tcs(tcs,Nmax,csin)
    implicit none
    type(totalcs), intent(inout):: tcs
    integer, intent(in):: Nmax
    real(dp), dimension(Nmax)::  csin

    integer:: i

    if(Nmax .gt. tcs%Nmax) then
       print*, 'totalcs.f90: set_cs_tcs(): Nmax .gt. tcs%Nmax', Nmax, tcs%Nmax
       stop
       ! should not happen this in our case, but might happen in general, do something?
    endif
    ! assume also that energy grid is the same...
    tcs%Nmax = Nmax
    do i=1,Nmax
       tcs%CS(i) = csin(i)
    enddo
  end subroutine set_cs_tcs
  
  subroutine normalisetcs(tcs,stateprob,numStates)
	implicit none 
	type(totalcs), intent(in):: tcs
	integer:: numStates,i
	real(dp), dimension(numStates),intent(out):: stateprob
	real(dp):: tcsSum

	! allocate the size of the array to the number of state in tcs 
	!numStates = tcs%Nmax
	!allocate(stateprob(numStates)) ! do i need to do some sort of deallocate statement as well?
		
	! Find normalised cross sections
	tcsSum = sum(tcs%CS)
	stateprob(1) = tcs%CS(1) / tcsSum
	!print*, tcs%nf(1),stateprob(1)  ! state number and cumulative probabilty
	do i=2,numStates
		stateprob(i) = stateprob(i-1) + (tcs%CS(i) / tcsSum)		
                !print*, stateprob(i), tcs%en(i)
		!print*, tcs%nf(i),stateprob(i)  ! state number and cumulative probabilty
        enddo
  end subroutine normalisetcs
!
!  
  
  subroutine meanexcenergy(tcs,statebasis,meanexc)
	use state_class
	type(totalcs),intent(in):: tcs
	type(basis_state),intent(in):: statebasis
	real(dp),intent(out):: meanexc
	integer:: i, poscontstate
	logical:: statefound
	real(dp):: csIon, totalCsIon 
	
	! Loop through and determine which state is the first continuum state
	! check if en is positive
	statefound = .false.
	i = 1
	do while((i .le. tcs%Nmax) .and. (statefound .eqv. .false.))
		!print*, statebasis%b(i)%en
		if(statebasis%b(i)%en .ge. 0) then
			statefound = .true.
			poscontstate = i
		end if
		i = i + 1
	end do
	
	meanexc = 0.0
	if(statefound) then
		!print*, 'first pos state is', poscontstate
		csIon = 0.0
		totalCsIon = 0.0
		do i = poscontstate, tcs%Nmax
			!csIon = csIon + exciteE(n)*cs(n)
			csIon = csIon + statebasis%b(i)%enex*tcs%CS(i)
			totalCsIon = totalCsIon + tcs%CS(i)
		end do
		meanexc = csIon / totalCsIon
		!print*, 'meanexc', meanexc, 'secEen',meanexc-15.96632
	else
		print*, 'no positive state found (meanexcenergy subrtn)'
	end if
	
  
  end subroutine meanexcenergy
 


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: testTcs
   !Purpose: prints totalcs to file as a function of incident energy
   !         for comparison with integrated SDCS as a consistency 
   !         check.
   !Date last modified: 13/07/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine testTcs(stateBasis,enfine,Nenfine)
       use state_class
       implicit none
       type(basis_state)::stateBasis
       real(dp):: eIncident, deltaE, tics
       type(totalcs)::tcsAtEn
       integer::Nenfine
       real(dp), dimension(:)::enfine

       !use smallest energy on fine energy grid
       eIncident = enfine(1)
       deltaE = 2.0
       open(27, file="totalcsIntTest.txt")
       write(27,*)  "Energy(ev),           totalIonisationCS(a.u)"
       do while (eIncident .lt. enfine(Nenfine))
          call get_csatein(stateBasis,eIncident,tcsAtEn)
          call getTics(tcsAtEn,tics)
 
          write(27,*) eIncident, tics
          eIncident = eIncident + deltaE
          call delete_totalcs(tcsAtEn)
       end do
       close(27)
   end subroutine testTCS




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: testTcsVib
   !Purpose: prints totalcs to file as a function of incident energy
   !         for comparison with integrated SDCS as a consistency 
   !         check. Different filename to above function
   !Date last modified: 13/07/2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine testTcsVib(stateBasis)
       use state_class
       implicit none
       type(basis_state)::stateBasis
       real(dp):: eIncident, deltaE, tics
       type(totalcs)::tcsAtEn

       eIncident = stateBasis%b(2)%enex
       deltaE = 0.05
       open(27, file="totalcsIntTestVib.txt")
       write(27,*)  "Energy(ev),           totalIonisationCS(a.u)"
       do while (eIncident .lt. 500.05)
          call get_csatein(stateBasis,eIncident,tcsAtEn)
          call getTics(tcsAtEn,tics)
 
          write(27,*) eIncident, tics
          eIncident = eIncident + deltaE
          call delete_totalcs(tcsAtEn)
       end do
       close(27)
   end subroutine testTcsVib





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: getTics
  !Purpose: gets the total ionisation cross section from the imported
  !         totalcs type
  !Date last modified: 28/03/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getTics(tcsIn, tics)
      implicit none
      type(totalcs), intent(in)::tcsIn
      real(dp)::tics
      integer:: ii

      tics = 0.d0
      do ii=1,tcsIn%Nmax 
         if (tcsIn%en(ii) .gt. 0.d0) then
            tics = tics + tcsIn%CS(ii)
         end if
      end do
  end subroutine getTics





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: symTcs
    !Purpose: symmetrises positive energy pseudostates in 
    !         distribution treatment 
    !NOTE: tcs energy measured in HARTEES not rydbergs. To 
    !      convert to eV need to multiply by 2*(13.6).
    !Date last modified: 02/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine symTcs(tcs,eIncident)
        use input_data
        implicit none
        type(totalcs)::tcs
        real(dp)::eIncident, eIon, cutoffEn, eV
        real(dp), allocatable, dimension(:)::enNew, csNew
        integer::zeroInd, cutoffInd, total
        integer::ii, kk

        !Find the indices corresponding to 0 and E/2        
        eIon = 15.96632  !Ionisation energy of ground state H2 (eV)
        cutoffEn = (eIncident - eIon)/2.0 !eV

        !Only need to symmetrise if projectile actually has energy 
        !required to cause ionisation. Some tcs don't have positive energy
        !states close to eIncident=eIon also.
        if ((eIncident .gt. eIon) .and. (tcs%en(tcs%Nmax)*data_in%eV .gt. 0.0)) then
           cutoffInd = 0
           zeroInd = 1
           !Compare energies in eV
           do while(tcs%en(zeroInd)*data_in%eV .lt. 0.0)
              zeroInd=zeroInd+1
           end do
           zeroInd = zeroInd-1

           cutoffInd = zeroInd
           do while((tcs%en(cutoffInd)*data_in%eV .lt. cutoffEn) &
                    .and. (cutoffInd .lt. tcs%Nmax))
              cutoffInd = cutoffInd + 1
           end do
           cutoffInd = cutoffInd -1 
           if (cutoffInd .eq. tcs%Nmax -1) then
              cutoffInd = tcs%Nmax
           end if
       
           total = zeroInd+2*(cutoffInd-zeroInd)
           allocate(enNew(total),csNew(total))
       
           enNew(1:zeroInd) = tcs%en(1:zeroInd)
           enNew(zeroInd+1:cutoffInd) = tcs%en(zeroInd+1:cutoffInd)
           csNew(1:zeroInd) = tcs%cs(1:zeroInd)
           csNew(zeroInd+1:cutoffInd) = tcs%cs(zeroInd+1:cutoffInd)
           kk=0
           !Store energies in tcs in a.u, as per previous convention
           eV = data_in%eV
           do ii = cutoffInd+1,total
              enNew(ii) = cutoffEn/eV + (cutoffEn/eV-tcs%en(cutoffInd-kk))
              csNew(ii) = tcs%cs(cutoffInd-kk)
              kk=kk+1
           end do

           !Copy over new symmetrised arrays
           deallocate(tcs%en,tcs%cs)
           allocate(tcs%en(total),tcs%cs(total))
           tcs%en(:) = enNew(:)
           tcs%cs(:) = csNew(:)
           tcs%Nmax = total
           deallocate(enNew,csNew)
        end if
    end subroutine symTcs



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: fitTotalCs
    !Purpose: interpolates pseudostate excitation cross sections
    !         to an exponentially decaying function, using given 
    !         energy grid. Grid energy assumed to be in Hartrees.
    !Date last modified:15/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fitTotalCs(tcs,enInArray,enLen)
        implicit none
        integer:: ii, counter, numPoints, enLen
        type(totalcs), intent(inout)::tcs
        real(dp)::a,b, currEn
        !a and b are interpolation paramters to find
        real(dp), allocatable, dimension(:)::lnArray, enArray               
        real(dp), allocatable, dimension(:)::enInArray, yVec
        real(dp), allocatable, dimension(:,:)::xMat
        real(dp), allocatable, dimension(:)::work
        integer::info

        a=0.0_dp
	      b=0.0_dp
        !Find the non-zero values and energies
        numPoints = 0
        do ii =1,tcs%Nmax
           if ((tcs%cs(ii) .gt. 0.0) .and. (tcs%en(ii) .gt. 0.0)) then
              numPoints = numPoints + 1
           end if
        end do   
        allocate(lnArray(numPoints),enArray(numPoints))
        lnArray(:) = 0.0
        enArray(:) = 0.0
        counter =1
        do ii =1,tcs%Nmax
           if ((tcs%cs(ii) .gt. 0.0) .and. (tcs%en(ii) .gt. 0.0)) then
              lnArray(counter) = log(tcs%cs(ii))
              enArray(counter) = tcs%en(ii)
              counter = counter +1
           end if
        end do   
        counter = counter-1

        if (counter .gt. 0) then
           !Alternate least squares fit, using DGELS least squares 
           !subroutine from SCALAPACK numerical linear algbra package.
           allocate(yVec(counter),xMat(counter,2),work(4))  
           xMat(:,1) = enArray(1:counter)
           xMat(:,2) = 1.0
           yVec(:) = lnArray(1:counter)
           info = 0
           call DGELS('N',counter,2,1,xMat,counter,yVec,counter,work,4,info) 
           a = yVec(1)
           b = yVec(2)
           deallocate(yVec,xMat,work)
        else
           !Do nothing, as there are no positive energy cross sections

        end if           

        deallocate(lnArray, enArray)
        deallocate(tcs%en, tcs%cs)

        allocate(tcs%en(enLen), tcs%cs(enLen))
        tcs%en(:) = enInArray(:)
        tcs%Nmax = enLen

        if (counter .gt. 0) then
           do ii=1, enLen
              currEn = enInArray(ii)
              tcs%cs(ii) = exp(a)*exp(b*currEn)
           end do
        else
           !There are no positive energy states
           tcs%cs(:) = 0.0
        end if

    end subroutine fitTotalCs



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: printFitTcs
    !Purpose: prints curve fitted to pseudostate excitation cross
    !         sections to file for plotting
    !Date last modified:08/12/2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine printFitTcs(tcs)
        implicit none
        type(totalcs)::tcs
        integer::ii
        character (len=3)::enString

        write(enString,'(I0)') int(tcs%en_incident)
        open(70, file="tcsFitted"//trim(enString)//"eV.txt")
        do ii =1, tcs%Nmax
           write(70,*) tcs%en(ii), tcs%cs(ii)
        end do
        close(70)
    end subroutine printFitTcs


  !____________________Read in and use benchmark excitation CS_________________!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: readTcsBenchmark
  !Purpose: reads in table of fit parameters for excitation cross 
  !         sections and dissociative ionisation, taken from paper
  !         by Garvey and Green. Also reads in defining quantum 
  !         numbers for the considered states. 
  !Date last modified: 09/02/2022
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readTcsBenchmark(tcsbasis,statebasis,Nenfine, enfine, TICS)
      use state_class
      use input_data
      implicit none
      type(basis_totalcs)::tcsbasis
      type(basis_state)::statebasis
      type(state)::tempState
      integer::Nenfine, ii, jj, numStates
      real(dp), dimension(:), allocatable:: TICS !Total ionization cross section to model ionisation.
      real(dp), dimension(:), allocatable::enfine, tempCS
      real(dp), dimension(:), allocatable::stateEnArray, AMat, omegaMat, gammaMat, vMat
      real(dp):: currEn, q, A, omega, gammaval, v, I
      real(dp)::  spin, stateEn
      integer:: m, ipar, remainder, numPatches
      character(len=10):: stlabel
      character(len=10):: ionString

      !Save a copy of the CCC elastic scattering CS, as
      !Garvey and Green do not provide one.
      call copy_state(tempState,statebasis%b(1))

      !Interpolate elastic scattering CS to new energy grid.
      allocate(tempCS(Nenfine))
      numPatches = int(real(Nenfine)/1000)
      do ii = 1, numPatches
         call INTRPL(tempState%inum,tempState%ein,tempState%cs,1000,enfine((ii-1)*1000+1:ii*1000),tempCS((ii-1)*1000+1:ii*1000))  
      end do
      remainder = Nenfine - numPatches*1000
      call INTRPL(tempState%inum,tempState%ein,tempState%cs,remainder,enfine(Nenfine-remainder:Nenfine), tempCS(Nenfine-remainder:Nenfine))
      deallocate(tempState%CS, tempState%ein)
      allocate(tempState%CS(Nenfine),tempState%ein(Nenfine))
      tempState%inum = Nenfine
      tempState%CS(:) = tempCS(:)
      tempState%ein(:) = enfine(:)
      deallocate(tempCS)
      tempState%CS(:) = 0.0


      if(associated(tcsbasis%b)) then
         call destruct_totalcsbasis(tcsbasis)
      end if 
      if(associated(statebasis%b)) then
         call destruct_basis(statebasis)
      end if
  
      allocate(tcsbasis%b(Nenfine))
 
      !Open the file, count the number of states  
      open(70, file="../DATA_1/benchmarkTcs.txt",action='read')
      read(70,*) numStates
      read(70,*)     

      allocate(statebasis%b(numStates+2))
      allocate(AMat(numStates), stateEnArray(numStates), omegaMat(numStates))
      allocate(gammaMat(numStates), vMat(numStates))     


      I = 16.0   !Ionisation potential(eV)
      do ii = 1, numStates
        !Read in parameters for state
        read(70,*) stlabel, stateEnArray(ii), m, ipar, spin, AMat(ii), omegaMat(ii), gammaMat(ii), vMat(ii)    
        call init_state(statebasis%b(ii+1), stlabel, stateEnArray(ii) - I, stateEnArray(ii),-1, m, ipar, spin, .false.)    
      end do         
      
      !Introduce a positive energy state modelling ionisation.
      ionString = "ION-H2"
      call init_state(statebasis%b(numStates+2),ionString, 30.0_dp, I, 0, 0, 0, 0.0_dp, .true.)  !State quantum numbers undefined, use 0.
      close(70)


      !Initialise each state's ein array.
      do ii = 1, numStates+2
         if (associated(statebasis%b(ii)%ein)) then
            deallocate(statebasis%b(ii)%ein)
         end if
 
         allocate(statebasis%b(ii)%ein(Nenfine))
         statebasis%b(ii)%inum = Nenfine
         statebasis%b(ii)%ein(:) = enfine(:)      
      end do

      q=(6.514e-14)*(0.01**2)/(data_in%bohrRadius**2)
      do ii=1, Nenfine
         currEn = enfine(ii)
         !Create a new totalcs at each incident energy.
         call new_totalcs(tcsbasis%b(ii), numStates, 1)  
         tcsbasis%b(ii)%en_incident = currEn
         tcsbasis%b(ii)%Nmax = numStates+1
         
         tcsbasis%b(ii)%en(:) = stateEnArray(:)
         tcsbasis%b(ii)%df(:) = 0.0 !Benchmark paper did not explicitly consider dissociation 

         !Fit CS using formula
         do jj = 1, numStates
            stateEn = stateEnArray(jj)
            A = AMat(jj)
            omega = omegaMat(jj)
            gammaval = gammaMat(jj)
            v = vMat(jj)
           
            if( stateEn .lt. currEn ) then
               tcsbasis%b(ii)%CS(jj) = (q*A/(stateEn**2))*((stateEn/currEn)**(omega))*(1-(stateEn/currEn)**gammaval)**v 
            else if (stateEn .gt. currEn) then
               !Incident particle has insufficient energy to excite state.
               tcsbasis%b(ii)%CS(jj) = 0.0
            end if
         end do
      end do     

      !Copy cross sections to state basis and deallocate
      do ii= 1, numStates
         if (associated(statebasis%b(ii+1)%CS)) then
            deallocate(statebasis%b(ii+1)%CS)
         end if
         allocate(statebasis%b(ii+1)%CS(Nenfine))

         do jj = 1, size(statebasis%b(ii+1)%CS)
            !values sroted in transpose order
            statebasis%b(ii+1)%CS(jj) =  tcsbasis%b(jj)%CS(ii)
         end do         
      end do
      statebasis%n = numStates+2

      !Set initial state CS to interpolated elastic scattering CS.
      call copy_state(statebasis%b(1), tempState) 
      
      if (associated(tempState%ein)) deallocate(tempState%ein)
      if (associated(tempState%CS)) deallocate(tempState%CS)
      if (associated(tempState%df)) deallocate(tempState%df)
      if (associated(tempState%eindf)) deallocate(tempState%eindf)

      !Set final state CS to total ionisation cross section/
      if (associated(statebasis%b(numStates+2)%CS)) then
         deallocate(statebasis%b(numStates+2)%CS)
      end if
      allocate(statebasis%b(numStates+2)%CS(Nenfine))
      statebasis%b(numStates+2)%CS(:) = TICS(:) 
      call destruct_totalcsbasis(tcsbasis)
      !Sort states by magnitude of excitation energy
      call sortEnergies(stateBasis)

      deallocate(AMat, stateEnArray, omegaMat, gammaMat, vMat)

  end subroutine readTcsBenchmark






end module totalcs_module
