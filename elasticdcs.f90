module dcs_module

  use numbers

  !Differential cross sections for a particular transition
  type, public:: dcs
     ! private
     integer:: tcstype
     real(dp):: en_incident  ! incident energy 
     real(dp):: avang  ! average scattering angle theta 
     integer:: Nth  ! number of angle grid points for DCS 
     real(dp), pointer, dimension(:):: DCS => NULL()   ! DCS     
     real(dp), pointer, dimension(:):: p => NULL()  ! distribution
     real(dp), pointer, dimension(:):: theta => NULL()  ! common grid for DCS
     real(dp):: integratedCs !Cross sections integrated over solid angle 
     real(dp):: momtcs   !Momentum transfer cross section at this incident energy
     real(dp):: intcs !Cross sections integrated over solid angle (actually used) 
     integer:: nf, ni, stateNum
  end type dcs

  type, public:: basis_dcs
     !     private
     integer:: tcstype
     real(dp), dimension(:,:), pointer:: DCS => NULL()
     real(dp), dimension(:,:), pointer:: p => NULL()
     real(dp), pointer, dimension(:):: avang => NULL()    ! average  angle
     real(dp), pointer, dimension(:):: ein => NULL()    ! incident energies
     real(dp), pointer, dimension(:):: theta => NULL()    ! common grid for angles
     !Basis of dcs type on common fine energy grid
     type(dcs), pointer, dimension(:):: dcsBasis => NULL() 
     integer:: Nth  ! number of angle grid points for DCS 
     integer:: Nindcs   ! the last point on the fine energy grid where DCS is available ky-->> to be done
     integer:: Nindcs_first ! the first point on the fine energy grid where DCS is availabele -->> to be done
     integer:: Ndcs  ! number of energy grid points for DCS, must be the same as Nenfine at the final step
     type(basis_dcsstate), pointer, dimension(:):: inelDcs => NULL()
  end type basis_dcs
  !

  !Stores differential cross sections for a particular stateNum
  !over a range of incident energies.
  type, public:: basis_dcsstate
     integer:: stateNum
     real(dp), dimension(:), pointer:: ein => NULL() 
     type(dcs), dimension(:), pointer:: dcsArray => NULL()
     integer:: numDcs
     integer:: minIndex !Index of first dcs above minimum excitation energy
     real(dp):: minEn    !Minimum energy at which state can be excited
  end type  


contains

  subroutine make_elastic_dcs(self,Nenfine,enfine)

    use input_data         ! definitions for input data type, construct routine, object data_in

    implicit none

    type(basis_dcs):: self
    integer, intent(in):: Nenfine
    real(dp), dimension(Nenfine), intent(in):: enfine

    integer:: Nang   ! maximum number of theta points
    parameter(Nang = 201)    

    integer:: INT
    real(dp):: FINT
    COMMON/INTEG/ INT,FINT

    character*20::  str1,str2,str3,str4,str5
    real(dp):: ein, vint, ONEpi, TWOpi, th, cs, avtheta, ang, thang
    real(dp), dimension(Nang):: tmptheta, dcs, y
    integer:: i, n, Ndcs, tcstype, im, itm, ien, halfInt
    real(dp), dimension(data_in%Ndcs):: einar
    real(dp), dimension(data_in%Ndcs+10):: yen, xen
    real(dp), dimension(Nenfine):: ven
    real(dp), dimension(:), allocatable:: theta, yth, v
    real(dp), dimension(:,:), allocatable:: dcsar, par
 
    halfInt = 0

    ONEpi = dacos(-1d0)
    TWOpi = 2d0 * ONEpi
    INT = 1   ! to get integration done in itrpl()

    self%tcstype = data_in%tcstype

    allocate(self%ein(Nenfine))
    self%ein(1:Nenfine) = enfine(1:Nenfine)

!  make common grid for all DCS
    print*, 'Elastic DCS'
    print*, 'Make common theta grid'
    i = 0
    ang = -0.5d0
    do 
       ang = ang + 0.5d0
       i = i + 1
       if(ang-1.0d-6 .gt. 180.0d0) exit
    enddo
    im = i - 1

    self%Nth = im
    print*, ' Nth = ', self%Nth

    allocate(self%theta(im), self%dcs(im,Nenfine) , self%p(im,Nenfine), self%avang(Nenfine))
    allocate(theta(im), yth(im), v(im))
    allocate(dcsar(im,data_in%Ndcs), par(im,data_in%Ndcs))

    ang = -0.5d0
    do i=1,im
       ang = ang + 0.5d0
       thang = (ang / 180.0d0) * ONEpi
       self%theta(i) = thang
       theta(i) = thang
!       print*, i, ang, thang
    enddo
!
!!$--------------------------------
!
    print*, '   incident energy   theta-av       theta-av(deg)     ICS  '

    do n=1,data_in%Ndcs

!   Read  DCS files
       open(10,file=data_in%filename_dcs(n),action='read')
       print*, 'Open DCS file ', data_in%filename_dcs(n)

       if (n .gt. data_in%Ndcs2) then
          !Read the first type of elastic dcs file
          read(10,*) str1,str2,str3,str4,ein
          einar(n) = ein
          do i=1,5
             read(10,*)
          enddo

          do i=1,Nang
             read(10,*,ERR=35,END=35) th, dcs(i)
             tmptheta(i) = (th/180.0d0) * ONEpi
!            print*, i, tmptheta(i), dcs(i)
          enddo
35        itm = i -1
          !     print*, 'itm =', itm
       else if (n .le. data_in%Ndcs2) then
          !Read second type of dcs file, reading only elastic dcs for now
          read(10,*)
          read(10,*) str1, str2, str3, str4, str5, ein
          einar(n) = ein
          read(10,*)

          do i=1,Nang
             read(10,*,ERR=36,END=36) th, str1, str2, str3, dcs(i)
             tmptheta(i) = (th/180.0d0) * ONEpi
!            print*, i, tmptheta(i), dcs(i)
          enddo
36        itm = i -1
       end if

       close(10)

! interpolate to the common theta grid standard for all DCS
       do i=1,itm
          y(i) = TWOpi * ( dcs(i) * sin(tmptheta(i)) )
!          print*, i, tmptheta(i), y(i)
       enddo
!       print*, itm, im
       INT = 1   ! to get integration done in itrpl()
       call INTRPL(itm,tmptheta,y,im,theta,v)
       cs = FINT
!          print*, 'cs = ', cs
       par(1:im,n) = v(1:im) / cs

       call INTRPL(itm,tmptheta,dcs,im,theta,v)
       dcsar(1:im,n) = v(1:im)

!!$ obtain the average angle
       do i=1,itm
          y(i) = TWOpi * ( dcs(i) * sin(tmptheta(i)) )*tmptheta(i)
       enddo
       call INTRPL(itm,tmptheta,y,im,theta,v)
       avtheta = FINT / cs
       print'("e, avtheta, avtheta(deg),ICS: ",1P,5E13.5)', ein, avtheta, avtheta*180d0/ONEpi, cs
       
    enddo  ! end of redaing and dealing with this elasctic DCS file
!
!
!!$ interpolate to the extended energy grid, from 0 to 4.4 eV, followed by the 
!   fine energy grid above 4.4 eV. Using old variable names.
    self%Ndcs = Nenfine
    
    do n=1,data_in%Ndcs          
       xen(n) = einar(n)
    enddo

    do i=1,self%Nth   ! for each angle on the common theta grid

       do n=1,data_in%Ndcs          
          yen(n) = par(i,n)
       enddo
       ien = data_in%Ndcs          
       !Interpolate over patches of the array and join
       halfInt = nint(real(Nenfine)/2.0)-1
       call INTRPL(ien,xen,yen,halfInt,enfine(1:halfInt),ven(1:halfInt))
       call INTRPL(ien,xen,yen,Nenfine-halfInt,enfine(halfInt+1:Nenfine),ven(halfInt+1:Nenfine))
       self%p(i,1:Nenfine) = ven(1:Nenfine)

       do n=1,data_in%Ndcs          
          yen(n) = dcsar(i,n)
       enddo
       ien = data_in%Ndcs          
       halfInt = nint(real(Nenfine)/2.0)-1
       call INTRPL(ien,xen,yen,halfInt,enfine(1:halfInt),ven(1:halfInt))
       call INTRPL(ien,xen,yen,Nenfine-halfInt,enfine(halfInt+1:Nenfine),ven(halfInt+1:Nenfine))
       self%dcs(i,1:Nenfine) = ven(1:Nenfine)
    enddo

!!$ obtain the average angle once again, this time on the fine energy grid
    do n=1,Nenfine
       do i=1,self%Nth
          yth(i) = self%p(i,n)*self%theta(i)
       enddo
       call INTRPL(self%Nth,self%theta,yth,self%Nth,self%theta,v)
       avtheta = FINT 
       self%avang(n) = avtheta
!       print'("2. ",1P,5E15.5)', enfine(n), avtheta, avtheta*180d0/ONEpi
    enddo


!!$  open(20,file='tmpfile')
!!$  
!!$  do i=1,itm
!!$     write(20,'(1P,10E15.5)') theta(i)*180d0/ONEpi, p(i), uav(i), dcs(i)
!!$  enddo
!!$  close(20)
  
    deallocate(theta,yth,v,dcsar,par)

  end subroutine make_elastic_dcs



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: make_inelastic_dcs
  !Purpose: read in dcs for inelastic scattering events and store in 
  !         dcs data structures.
  !Date last modified: 05/03/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_inelastic_dcs(dcsBasis,stateBasis)
      use input_data    !Contains array of dcs filenames
      use state_class
      implicit none
      type(basis_dcs)::dcsBasis
      type(dcs), dimension(:), allocatable::newDcsArray
      type(basis_state)::stateBasis
      integer:: ii, numFiles, numLines, numStates, maxStates
      integer:: firstSpace, colSpace, counter, ind, jj, kk, space
      integer:: ioStatVar, l, enInd, stateInd, numFilesCor, ll
      integer:: ind1, ind2
      real(dp):: w1, w2
      logical:: found
      real(dp)::thetaVal, dcsVal, currEn, elasticVal, currEnPrev
      character*10000::line
      real(dp), dimension(:), allocatable::thetaArray, dcsCorrArray, thetaArray1, intpDCS
      real(dp), dimension(:,:), allocatable::dcsArray
      character*30::str1, str2, str3, str4, str5, str6, str7, str8, str9
 
      if (associated(dcsBasis%inelDcs)) then
         do ii = 1, size(dcsBasis%inelDcs)   
            if (associated(dcsBasis%ineldcs(ii)%ein)) then
               deallocate(dcsBasis%ineldcs(ii)%ein)
            end if

            if (associated(dcsBasis%ineldcs(ii)%dcsArray)) then
               call freeDcsBasis(dcsBasis%ineldcs(ii)%dcsArray,dcsBasis%ineldcs(ii)%numDcs)
               deallocate(dcsBasis%ineldcs(ii)%dcsArray)
            end if
         end do
      end if

      !!!!!!!!!!!!!!!Reads in regular inelastic dcs files!!!!!!!!!!!

      numFiles = data_in%Ndcsinel
      !Find largest stateNum available in input files
      maxStates = 0   
      numStates = 0
      do ii=1, numFiles
         call getNumStatesDcs(data_in%filename_dcsinel(ii),numStates,space)
         if (numStates .gt. maxStates) then
            maxStates = numStates
         end if
      end do

      !Determine number of points on energy grid by using number of files
      allocate(dcsBasis%inelDcs(numStates))
      do ii = 1, numStates
         allocate(dcsBasis%inelDcs(ii)%dcsArray(numFiles),dcsBasis%inelDcs(ii)%ein(numFiles))
      end do      
     
      !Set minimum excitation energies
      do ii = 1, maxStates
         dcsBasis%inelDcs(ii)%minEn = stateBasis%b(ii)%enex 
      end do    

      currEnPrev = 0.0
      !Loop through all files
      do ii = 1, numFiles
         !For a particular incident energy read in dcs for each available state
         !i.e. Each column in the file fills a single dcs type.

         !Find number of states included in this particular file, i.e. count columns.
         numStates = 0
         call getNumStatesDcs(data_in%filename_dcsinel(ii),numStates,space)
         call getNumLines(data_in%filename_dcsinel(ii),numLines)

         open(70,file=data_in%filename_dcsinel(ii),action='read')
         !Read incident energy from first line
         read(70,*) str1, str2, str3, str4, currEn, str6, str7, str8, str9 

         do kk = 2, 2*numStates+space+1
            !Ignore first few lines in input file. See file format.
            read(70,*)
         end do

         !Read in all columns, making use of iostat to prevent crashing.
         numLines = numLines - 2*numStates-space-1
         allocate(thetaArray(numLines), dcsArray(numLines,numStates))
         do jj = 1, numLines
            read(70,*, iostat=ioStatVar) thetaArray(jj),  (dcsArray(jj,l),l=1,numStates)
         end do

         if (currEnPrev .lt. dcsBasis%inelDcs(numStates)%minEn) then
            dcsBasis%inelDcs(numStates)%minIndex = ii
         end if

         !Loop through dcs array and store values
         if( ioStatVar .eq. 0) then
            do jj = 1, numStates 
               dcsBasis%inelDcs(jj)%ein(ii) = currEn

               !Initialise dcs type for given eIncident and stateNum
               if (.not. (associated(dcsBasis%inelDcs(jj)%dcsArray(ii)%DCS))) then
                  allocate(dcsBasis%inelDcs(jj)%dcsArray(ii)%DCS(numLines))
                  dcsBasis%inelDcs(jj)%dcsArray(ii)%DCS(:) = 0.0
               end if
               if(.not. (associated(dcsBasis%inelDcs(jj)%dcsArray(ii)%theta))) then
                  allocate(dcsBasis%inelDcs(jj)%dcsArray(ii)%theta(numLines))
                  dcsBasis%inelDcs(jj)%dcsArray(ii)%theta(:) = 0.0
               end if
               if(.not. (associated(dcsBasis%inelDcs(jj)%dcsArray(ii)%p))) then
                  allocate(dcsBasis%inelDcs(jj)%dcsArray(ii)%p(numLines))
                  dcsBasis%inelDcs(jj)%dcsArray(ii)%p(:) = 0.0  !Set in get_dcsatein
               end if

               dcsBasis%inelDcs(jj)%dcsArray(ii)%Nth = numLines
               dcsBasis%inelDcs(jj)%dcsArray(ii)%en_incident = currEn 

               !Copy values to arrays in dcs type. 
               dcsBasis%inelDcs(jj)%dcsArray(ii)%theta(:) = thetaArray(:)
               dcsBasis%inelDcs(jj)%dcsArray(ii)%DCS(:) = dcsArray(:,jj)
               dcsBasis%inelDcs(jj)%dcsArray(ii)%stateNum = jj          
            end do
            do jj = numStates+1, maxStates 
               if (.not. (associated(dcsBasis%inelDcs(jj)%dcsArray(ii)%DCS))) then
                  allocate(dcsBasis%inelDcs(jj)%dcsArray(ii)%DCS(numLines))
                  dcsBasis%inelDcs(jj)%dcsArray(ii)%DCS(:) = 0.0
               end if
               if(.not. (associated(dcsBasis%inelDcs(jj)%dcsArray(ii)%theta))) then
                  allocate(dcsBasis%inelDcs(jj)%dcsArray(ii)%theta(numLines))
                  dcsBasis%inelDcs(jj)%dcsArray(ii)%theta(:) = 0.0
               end if
               if(.not. (associated(dcsBasis%inelDcs(jj)%dcsArray(ii)%p))) then
                  allocate(dcsBasis%inelDcs(jj)%dcsArray(ii)%p(numLines))
                  dcsBasis%inelDcs(jj)%dcsArray(ii)%p(:) = 0.0  !Set in get_dcsatein
               end if
 
               dcsBasis%inelDcs(jj)%dcsArray(ii)%Nth = numLines
               dcsBasis%inelDcs(jj)%dcsArray(ii)%en_incident = currEn 

               !Copy values to arrays in dcs type. 
               !Leave values initialised to zero, as channel is closed.
               dcsBasis%inelDcs(jj)%dcsArray(ii)%stateNum = jj          
            end do
         else if (ioStatVar .ne. 0) then
            print*, "ERROR: error occured while reading file: ", data_in%filename_dcsinel(ii), "ABORT"
            stop
         end if         

         currEnPrev = currEn

         close(70)
         deallocate(thetaArray, dcsArray)
      end do         
      allocate(thetaArray1(size(dcsBasis%inelDcs(2)%dcsArray(1)%theta)))
      thetaArray1(:) = dcsBasis%inelDcs(2)%dcsArray(1)%theta(:)

      !!!!!!!!!!Reads in corrected cross sections for excitationt to b1Su state!!!!!!!!!!!!!!
       
      numFilesCor = data_in%Ndcsinelcor
      allocate(newDcsArray(numFilesCor))
      do ii=1, numFilesCor         
         call getNumLines(data_in%filename_dcsinelcor(ii), numLines)

         open(70,file=data_in%filename_dcsinelcor(ii),action='read')
         read(70,*)  str1, str2, str3, str4, currEn, str6, str7, str8
         do jj=2, 3
            !Ignore first few lines
            read(70,*)  
         end do
         
         allocate(thetaArray(numLines-3),dcsCorrArray(numLines-3))
         do jj =1,numLines-3
            read(70,*) thetaArray(jj), elasticVal, dcsCorrArray(jj)
         end do

         !Interpolate data to angular grid used in standard input files.
         allocate(intpDCS(size(thetaArray1)))
         intpDCS(:) = 0.0
         ll=1
         do kk=1,size(thetaArray1)
            found = .false.
            do while((thetaArray(ll) .lt. thetaArray1(kk)) .and. (ll .le. size(thetaArray))) 
               ll= ll+1
            end do
            if (ll .eq. 1) then
               !Only occurs if thetaArray1(1) < thetaArray(1). Since they should be the same,
               !this can only occur as a result of numerical error, in which case just use
               !cs from file.
               w1 = 1.0 
               w2 = 0.0
               ind1 = ll 
               ind2 = ll+1
            else if (ll .eq. size(thetaArray)) then
               !Same as for the above case, but at the end of the array.
               w1 = 1.0
               w2 = 0.0 
               ind1 = ll
               ind2 = ll-1
            else
               w1 = (thetaArray1(kk) - thetaArray(ll-1))/(thetaArray(ll-1)-thetaArray(ll))
               w2 = (thetaArray(ll)- thetaArray1(kk))/ (thetaArray(ll-1)-thetaArray(ll))
               ind1 = ll-1
               ind2 = ll
            end if
            intpDCS(kk) = w1*dcsCorrArray(ind1) + w2*dcsCorrArray(ind2)
         end do

         !Create new dcs and insert into dcsArray in dcsBasis.
         allocate(newDcsArray(ii)%theta(size(thetaArray1)), newDcsArray(ii)%DCS(size(thetaArray1)), newDcsArray(ii)%p(size(thetaArray1)))
         newDcsArray(ii)%theta(:) = thetaArray1(:)
         newDcsArray(ii)%DCS(:) = intpDCS(:)
         newDcsArray(ii)%stateNum = 2    !Index of b3Su state.
         newDcsArray(ii)%en_incident = currEn  
         newDcsArray(ii)%p(:) = 0.0 !Set later in get_dcsatein
         newDcsArray(ii)%Nth = size(thetaArray1)
 
         deallocate(thetaArray, dcsCorrArray,intpDcs)
         close(70)
      end do
      deallocate(thetaArray1)
      call placeDcs(newDcsArray,numFilesCor,dcsBasis) 
      call freeDcsArray(newDcsArray,numFilesCor)
 

  end subroutine make_inelastic_dcs





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: getNumStatesDcs
  !Purpose: reads dcs input file and finds the number of states in the file
  !         for which dcs are available. ie. The number of columns bar the 
  !         first.
  !Date last modified: 07/03/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getNumStatesDcs(filename,numStates,space)
      implicit none
      character(len=80), intent(in)::filename
      integer, intent(out)::numStates, space
      integer::ii, charNum, counter
      character*1000::line
      character*6::lineString, thetaString, unitString, angleString
      integer:: unitCounter

      open(50,file=filename,action='read')
   
      thetaString = "#theta"
      unitString = "#units"
      angleString = "#angle"

      counter = 1
      read(50, *) line
      lineString = line(1:6)
      do while(.not. ((lineString .eq. thetaString) .or. (lineString .eq. angleString)))
         read(50,*) line
         lineString = line(1:6)
         counter = counter +1
         if (lineString .eq. unitString) then
            unitCounter = counter
         end if
      end do

      !Number of states is given by the number of lines of data
      !printed to the top of the input file. In general, two lines
      !per states, followed by four lines of additional info before
      !data actually starts.
      numStates = (unitCounter -1)/2

      !Space is the number of lines between end of info about each
      !transition and beginning of data columns. Usually 4 or 5 lines
      space = counter - unitCounter
      close(50) !Close file to avoid memory leaks
  end subroutine getNumStatesDcs



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: getNumLines
   !Purpose: gets the number of lines in the input file
   !Date last modified: 10/03/2021
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine getNumLines(filename,numLines)
       implicit none
       integer, intent(out):: numLines
       character(len=*)::filename
       integer::counter
       integer:: ii     
 
       counter = 0
       open(35,file=filename,action='read')
       do 
          read(35,*, END=20) 
          counter = counter+1
       end do
20     continue
       numLines = counter
       close(35)

   end subroutine getNumLines






    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Subroutine: placeDcs
    !Purpose: places the imported array of dcs types in the dcsbasis in  
    !         positions ordered by incident energy.
    !Date last modified: 12/03/2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine placeDcs(dcsIn, inLen, dcsBasis) 
        implicit none
        type(dcs), dimension(:), allocatable::dcsIn, tempArray
        type(basis_dcs)::dcsBasis
        real(dp):: currEn, tolerance
        logical::found
        integer:: enIndex, inLen, basisSize, numCommon
        integer::ii, jj, kk, ll
        
        tolerance = 0.00001

        basisSize = size(dcsBasis%inelDcs(2)%dcsArray)
        !Copy over the current array
        allocate(tempArray(basisSize))
        do ii = 1, basisSize   
           call copyDcs(tempArray(ii), dcsBasis%inelDcs(2)%dcsArray(ii)) 
        end do
        
        !Deallocate and reallocate with larger size.           
        deallocate(dcsBasis%inelDcs(2)%dcsArray, dcsBasis%inelDcs(2)%ein)
        
        !Count number of common energies
        numCommon = 0
        jj = 1
        kk = 1
        do while( (jj .le. basisSize) .and. (kk .le. inLen))
           !print*, "basisSize, inLen, jj, kk: ", basisSize, inLen, jj, kk
           if ((abs(tempArray(jj)%en_incident - dcsIn(kk)%en_incident)) .lt. tolerance) then
              !If same energy, use values from file
              numCommon = numCommon +1
              jj = jj+1 
              kk = kk+1
           else if (tempArray(jj)%en_incident .gt. dcsIn(kk)%en_incident) then
              kk = kk+1 
           else if (tempArray(jj)%en_incident .lt. dcsIn(kk)%en_incident) then
              jj = jj+1
           end if
        end do
        allocate(dcsBasis%inelDcs(2)%dcsArray(basisSize+inLen-numCommon),dcsBasis%inelDcs(2)%ein(basisSize+inLen-numCommon))
        dcsBasis%inelDcs(2)%numDcs = basisSize+inLen-numCommon

        !Copy back and order.
        jj = 1
        kk = 1
        ii = 1
        do while((jj .le. basisSize) .and. (kk .le. inLen))
           !print*, "basisSize, inLen, jj, kk: ", basisSize, inLen, jj, kk
           !print*, tempArray(kk)%en_incident, dcsIn(kk)%en_incident
           if ((abs(tempArray(jj)%en_incident - dcsIn(kk)%en_incident)) .lt. tolerance) then
              !If same energy, use values from file
              call copyDcs(dcsBasis%inelDcs(2)%dcsArray(ii), dcsIn(kk))
              dcsBasis%inelDcs(2)%ein(ii) = dcsIn(kk)%en_incident
              jj = jj+1 
              kk = kk+1
           else if (tempArray(jj)%en_incident .gt. dcsIn(kk)%en_incident) then
              call copyDcs(dcsBasis%inelDcs(2)%dcsArray(ii), dcsIn(kk))
              dcsBasis%inelDcs(2)%ein(ii) = dcsIn(kk)%en_incident
              kk = kk+1 
           else if (tempArray(jj)%en_incident .lt. dcsIn(kk)%en_incident) then
              call copyDcs(dcsBasis%inelDcs(2)%dcsArray(ii), tempArray(jj))
              dcsBasis%inelDcs(2)%ein(ii) = tempArray(jj)%en_incident
              jj = jj+1
           end if
           ii = ii+1
        end do

        !Copy over remaining dcs
        if (jj .le. basisSize) then
           do while(jj .le. basisSize) 
              call copyDcs(dcsBasis%inelDcs(2)%dcsArray(ii), tempArray(jj))
              dcsBasis%inelDcs(2)%ein(ii) = tempArray(jj)%en_incident
              ii = ii+1
              jj = jj+1
           end do
        else if (kk .le. inLen) then
           do while(kk .le. inLen)
              call copyDcs(dcsBasis%inelDcs(2)%dcsArray(ii), dcsIn(kk))
              dcsBasis%inelDcs(2)%ein(ii) = dcsIn(kk)%en_incident
              ii = ii+1
              kk = kk+1
           end do
        end if
        call freeDcsArray(tempArray,basisSize)
    end subroutine placeDcs













  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: copyDcs
  !Purpose: copies the contents of the imported dcs type to a new instance
  !         of that type.
  !Date last modififed: 02/02/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine copyDcs(newDcs, dcsToCopy)
      implicit none
      type(dcs):: newDcs, dcsToCopy

      allocate(newDcs%DCS(size(dcsToCopy%DCS)),newDcs%p(size(dcsToCopy%p)),newDcs%theta(size(dcsToCopy%theta)))
      newDcs%tcstype = dcsToCopy%tcstype
      newDcs%en_incident= dcsToCopy%en_incident
      newDcs%avang = dcsToCopy%avang
      newDcs%Nth = dcsToCopy%Nth
      newDcs%DCS(:) = dcsToCopy%DCS(:)
      newDcs%p(:) = dcsToCopy%p(:)
      newDcs%theta(:) = dcsToCopy%theta(:)
      newDcs%integratedCs = dcsToCopy%integratedCs
      newDcs%momtcs = dcsToCopy%momtcs
      newDcs%intcs = dcsToCopy%intcs

  end subroutine copyDcs


 


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Subroutine: get_dcsatein
 !Purpose: returns a dcs type at the given incident 
 !         energy for a given transition.
 !Date last modified: 05/03/2021
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine get_dcsatein(dcsob,eincident,dcsen,stateNum)
    implicit none

    type(basis_dcs), intent(in):: dcsob
    real(dp), intent(in):: eincident
    type(dcs):: dcsen
    real(dp):: e1, e2, w1, w2, cs1, cs2
    integer:: i, kk, j, j1, j2, Nth, Ndcs
    real(dp)::leftVal, rightVal, leftMom, rightMom
    integer::stateNum   !Index of the final state selected.


    if ( associated(dcsen%DCS) ) deallocate(dcsen%DCS)
    if ( associated(dcsen%p) ) deallocate(dcsen%p)
    if ( associated(dcsen%theta) ) deallocate(dcsen%theta)

    !Treat H2 1->1 transition as before.
    if (stateNum .eq. 1) then
       Nth = dcsob%Nth
       dcsen%tcstype = dcsob%tcstype
       dcsen%Nth = dcsob%Nth
       dcsen%en_incident = eincident
       dcsen%stateNum = stateNum

       allocate(dcsen%DCS(Nth), dcsen%p(Nth), dcsen%theta(Nth))

       !Search dcs object for array of dcs at various energies for given stateNum

       dcsen%theta(1:Nth) = dcsob%theta(1:Nth)

       j = 1
       do i=1,dcsob%Ndcs
          if(eincident .lt. dcsob%ein(i)) then
             exit
          end if
          j = j + 1
       enddo
       if(j .eq. 1) then
          print*,'elasticdcs.f90: get_dcsatein(): j=1'
          print*, ' could not find the required two energies'
          stop
       endif

       j1 = j -1
       j2 = j

   !    print*, 'eincident = ', eincident
   !    print*, 'j1, j2:', j1, j2
   !    print*, 'dcsob%ein(j2) = ', dcsob%ein(j2)
    

       e1 = dcsob%ein(j1)
       e2 = dcsob%ein(j2)
       w1 = (eincident - e1)/(e2 - e1) 
       w2 = (e2 - eincident)/(e2 - e1)  

       dcsen%DCS(1:Nth) = w1*dcsob%DCS(1:Nth,j1) + w2*dcsob%DCS(1:Nth,j2)
       dcsen%p(1:Nth) = w1*dcsob%p(1:Nth,j1) + w2*dcsob%p(1:Nth,j2)
  

       dcsen%avang = w1*dcsob%avang(j1) + w2*dcsob%avang(j2)

       dcsen%momtcs = 0.0
       dcsen%intcs = 0.0
       !Calculate momentum transfer cross section and integrated dcs
       do i=1,dcsen%Nth-1
          !Use trapezoidal integration rule
          leftMom = dcsen%DCS(i)*(1-cos(dcsen%theta(i)))*sin(dcsen%theta(i))
          rightMom = dcsen%DCS(i+1)*(1-cos(dcsen%theta(i+1)))*sin(dcsen%theta(i+1))
          dcsen%momtcs = dcsen%momtcs + ((leftMom+rightMom)/2.0)*(dcsen%theta(i+1)-dcsen%theta(i))
 
          leftVal = dcsen%DCS(i)*sin(dcsen%theta(i))
          rightVal = dcsen%DCS(i+1)*sin(dcsen%theta(i+1))  
          dcsen%intcs = dcsen%intcs + ((leftVal+rightVal)/2.0)*(dcsen%theta(i+1)-dcsen%theta(i))
        end do
     else if (stateNum .gt. 1) then
        if (eincident .lt. dcsob%inelDcs(stateNum)%minEn) then
           print*, "ERROR: attempt to get dcs for excitation to stateNum = ", stateNum, " at eincident below minimum excitation energy"
           print*, "STOPPING"
           stop
        end if

        !Use inelastic dcs for given state, taken from array of dcs in dcsob.
        Nth = dcsob%ineldcs(stateNum)%dcsArray(1)%Nth   !Common angular grid used
        dcsen%tcstype = dcsob%tcstype
        dcsen%Nth = dcsob%ineldcs(stateNum)%dcsArray(1)%Nth
        dcsen%en_incident = eincident
        dcsen%stateNum = stateNum
        Ndcs = dcsob%ineldcs(stateNum)%numDcs

        allocate(dcsen%DCS(Nth), dcsen%p(Nth), dcsen%theta(Nth))

        !Interpolate dcs for given state to given energy
        dcsen%theta(1:Nth) = dcsob%ineldcs(stateNum)%dcsArray(dcsob%ineldcs(stateNum)%minIndex)%theta(1:Nth)

        j = 1
        do i=1, Ndcs
           if(eincident .lt. dcsob%ineldcs(stateNum)%ein(i)) then
              exit
           end if
           j = j + 1
        enddo
        if(j .eq. 1) then
           print*,'elasticdcs.f90: get_dcsatein(): j=1'
           print*, ' could not find the required two energies'
           stop
        endif

        j1 = j -1
        j2 = j

        e1 = dcsob%ein(j1)
        e2 = dcsob%ein(j2)
        w1 = (eincident - e1)/(e2 - e1) 
        w2 = (e2 - eincident)/(e2 - e1)  

        dcsen%DCS(1:Nth) = w1*dcsob%ineldcs(stateNum)%dcsArray(j1)%DCS(1:Nth) + w2*dcsob%ineldcs(stateNum)%dcsArray(j2)%DCS(1:Nth)
        dcsen%p(1:Nth) = w1*dcsob%ineldcs(stateNum)%dcsArray(j1)%p(1:Nth) + w2*dcsob%ineldcs(stateNum)%dcsArray(j2)%p(1:Nth)
 
        dcsen%avang = w1*dcsob%ineldcs(stateNum)%dcsArray(j1)%avang + w2*dcsob%ineldcs(stateNum)%dcsArray(j2)%avang

        dcsen%momtcs = 0.0
        dcsen%intcs = 0.0
        !Calculate momentum transfer cross section and integrated dcs
        do i=1,Nth-1
           !Use trapezoidal integration rule
           leftMom = dcsen%DCS(i)*(1-cos(dcsen%theta(i)))*sin(dcsen%theta(i))
           rightMom = dcsen%DCS(i+1)*(1-cos(dcsen%theta(i+1)))*sin(dcsen%theta(i+1))
           dcsen%momtcs = dcsen%momtcs + ((leftMom+rightMom)/2.0)*(dcsen%theta(i+1)-dcsen%theta(i))
 
           leftVal = dcsen%DCS(i)*sin(dcsen%theta(i))
           rightVal = dcsen%DCS(i+1)*sin(dcsen%theta(i+1))  
           dcsen%intcs = dcsen%intcs + ((leftVal+rightVal)/2.0)*(dcsen%theta(i+1)-dcsen%theta(i))
         end do
         dcsen%p(:) = dcsen%DCS(:)/dcsen%intCs
     end if

  end subroutine get_dcsatein


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: fillDcsBasis
  !Purpose: fills the array of dcs type to match the number of points on 
  !         the fine energy grid.
  !Date last modified: 05/03/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fillDcsBasis(dcsob)
      implicit none
      type(basis_dcs)::dcsob !contains basis to fill
      integer:: ii 

      if ( associated( dcsob%dcsBasis ) ) then
         call freeDcsBasis( dcsob%dcsBasis, dcsob%Ndcs )
      end if
      allocate(dcsob%dcsBasis(dcsob%Ndcs))  
 
      do ii = 1, dcsob%Ndcs-1
         call get_dcsatein(dcsob,dcsob%ein(ii),dcsob%dcsBasis(ii),1) 
      end do 
  end subroutine fillDcsBasis





  subroutine print_dcs(self,filename)

    implicit none
    
    type(dcs), intent(in):: self
    character*40, intent(in):: filename
    integer:: i
    real(dp):: ONEpi, th

    ONEpi = dacos(-1d0)

    open(10,file=filename)
    print*, 'open DCS file: ', filename
    
    write(10,'("# DCS")')
    write(10,'("# tcstype =",i5)') self%tcstype
    write(10,'("# stateNum =", i5)') self%stateNum
    write(10,'("# en_incident =",1P,E15.5)') self%en_incident
    write(10,'("# average angle: avang =",1P,E15.5)') self%avang
    write(10,'("#  i      theta      theta(deg)     DCS (a.u.)")')

    do i=1,self%Nth
       th = self%theta(i)
       write(10,'(i5,1P,2E15.5,5X,E15.5)'), i, th, th*(180.0d0/ONEpi), self%DCS(i)
    enddo
    
    close(10)

  end subroutine print_dcs




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: getAvAngle
  !Purpose: performs a linear interpolation to get cosine of average 
  !         scattering angle at the given incident energy
  !         from preprocessed dcs types.
  !Date last modified: 05/03/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getAvAngle(dcsob, eIncident, cosang)
     implicit none
     type(basis_dcs):: dcsob
     real(dp):: avang !Average scattering angle to find
     real(dp):: cosang
     !For interpolation
     real(dp):: ang1, ang2
     real(dp):: eIncident
     !Energies just below and above eIncident in find energy grid.
     real(dp):: e1, e2 
     integer:: ii

     ii=1
     do while( dcsob%ein(ii) .lt. eIncident )
        ii = ii + 1
     end do
   
     e1 = dcsob%ein(ii-1)
     e2 = dcsob%ein(ii)

     ang1 = dcsob%dcsBasis(ii-1)%avang
     ang2 = dcsob%dcsBasis(ii)%avang 
  
     !Performs linear interpolation to get ave scattering angle
     avang = ang1 + ((ang2-ang1)/(e2-e1))*(eIncident-e1) 	
     cosang = cos(avang)

  end subroutine getAvAngle


 subroutine delete_dcs(self)

    implicit none

    type(dcs):: self
    
    if ( associated(self%DCS) ) deallocate(self%DCS)
    if ( associated(self%p) ) deallocate(self%p)
    if ( associated(self%theta) ) deallocate(self%theta)

  end subroutine delete_dcs



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: destruct_basis_dcs
  !Purpose: deallocates a 'basis_dcs' type to clean up memory
  !Date last modified: 06/02/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine destruct_basis_dcs(self)
      implicit none
      type(basis_dcs) ::self
      if ( associated(self%DCS)) then
         deallocate(self%DCS)
      endif
      if ( associated(self%p))then
         deallocate(self%p)
      endif
      if ( associated(self%avang)) then
         deallocate(self%avang)
      endif
      if ( associated(self%ein)) then
         deallocate(self%ein)
      endif
      if (associated(self%theta)) then
         deallocate(self%theta)
      endif
      if (associated(self%dcsBasis)) then
         call freeDcsBasis(self%dcsBasis, self%Ndcs)
      end if
  end subroutine



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: freeDcsBasis
  !Purpose: frees the array of dcs type
  !Date last modified: 05/03/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine freeDcsBasis(dcsBasis, Ndcs)
      implicit none
      integer:: Ndcs, ii
      type(dcs), pointer, dimension(:):: dcsBasis
     
      
      if (associated(dcsBasis)) then
         do ii = 1, Ndcs
            call delete_dcs(dcsBasis(ii))
         end do
      end if

      deallocate(dcsBasis)
  end subroutine freeDcsBasis 



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: freeDcsArray
  !Purpose: frees the array of dcs type
  !Date last modified: 12/03/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine freeDcsArray(dcsArray, Ndcs)
      implicit none
      integer:: Ndcs, ii
      type(dcs), dimension(:), allocatable:: dcsArray
     
      
      if (allocated(dcsArray)) then
         do ii = 1, Ndcs
            call delete_dcs(dcsArray(ii))
         end do
      end if

      deallocate(dcsArray)
  end subroutine freeDcsArray 

end module dcs_module
