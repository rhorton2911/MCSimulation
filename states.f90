module state_class

  use numbers

  public::  new_basis, init_state, copy_state, basis_size, get_energy, get_energy_ex, get_ang_mom, get_ang_mom_proj, get_parity, get_spin, get_cs_array, get_ein_array, get_cs_i, get_ein_i, print_energy, sort_by_energy, set_cs, set_df, print_state


  type, public:: state
!     private
     character(len=10):: stlabel    ! state label or '?'
     real(dp):: en    ! ionization energy   in eV
     real(dp):: enex  ! excitation energy   in eV
     integer:: l    ! angular momentum 
     integer:: m    ! ang. mom. projection m  
     integer:: ipar ! parity
     real(dp):: spin  ! spin 
     integer:: inum ! number of energy points for cross sections
     real(dp), pointer, dimension(:) :: cs => NULL(), ein => NULL()   ! in a.u. : cs and ein (ionization energy) 
     !Energy dependent excitation energy (eV) defined on grid ein (above). 
     real(dp), pointer, dimension(:) :: enEin => NULL()   
     integer:: idf ! number of energy points for diss. fractions
     real(dp), pointer, dimension(:) :: df => NULL(), eindf => NULL()
     !Vibrational part of state.
     logical:: resolved    !True if state is a vibrational level of a vibrationally resolved electronic state, false if only electronic
     integer:: v !vibrational quantum number
     real(dp):: dissThresh   !threshold energy (eV) for dissociative excitation(from ground) of state with given electronic label
     logical:: ion   !Tracks whether state corresponds to ionisation (i.e. Is the state bound?)
     logical:: psFormation  !Tracks if state corresponds to postronium formation
     
  end type state
  !  
  !  
  type, public:: basis_state
     !     private
     type(state), pointer, dimension(:) :: b  => NULL()
     integer:: n
     integer:: tcstype   ! type of totalcs files: 1 for electrons, 2 for other primary particles
  end type basis_state
  !  
contains
  !
  subroutine init_state(self,stlabel,en,enex,l,m,ipar,spin,ion)
    implicit none
    type(state), intent(inout):: self
    character(len=10), intent(in):: stlabel
    real(dp), intent(in):: en,enex
    integer, intent(in):: l, m, ipar
    real(dp), intent(in):: spin
    logical:: ion

    self%stlabel = stlabel
    self%en = en
    self%enex = enex
    self%l = l
    self%m = m
    self%ipar = ipar
    self%spin = spin
    self%inum = 0
    self%idf = 0
    self%dissThresh = 0.0d0
    self%ion = ion
    self%resolved = .false. 
    self%psFormation = .false.

  end subroutine init_state
  !
   subroutine set_cs(self,inum,ar_ein,ar_cs,iadd)
    implicit none
    type(state), intent(inout):: self
    integer, intent(in):: inum
    real(dp), dimension(inum), intent(in):: ar_ein, ar_cs
    integer:: iadd

    if(iadd .ne. 0) then
       iadd = 1   ! add an energy point equal to the excitation energy
    endif

    self%inum = inum + iadd
    ! create arrays ein and cs
    if(associated(self%ein)) then
       deallocate(self%ein, self%cs)
    endif
    allocate(self%ein(inum+iadd),self%cs(inum+iadd))
    self%ein(1+iadd:inum+iadd) = ar_ein(1:inum)
    self%cs(1+iadd:inum+iadd) = ar_cs(1:inum)

    if(iadd .eq. 0) return
    if( .not. (self%enex .gt. 0.0d0)) return  ! not for elastic scattering, self%enex always positive
! here assuming that iadd=1
    self%ein(1) = self%enex 
    self%cs(1) = 0.0

    call placeen(1,self%ein(1),self%cs(1),self%inum, self%ein, self%cs)

  end subroutine set_cs
!
  subroutine set_df(self,idf,ar_ein,ar_df)
    implicit none
    type(state), intent(inout):: self
    integer, intent(in):: idf
    real(dp), dimension(idf), intent(in):: ar_ein, ar_df
    integer:: j, iadd

    iadd = 0
    if (idf .eq. 0) then
       iadd = 1
    end if
    self%idf = idf + iadd

    ! create arrays diss. fractions
    if(associated(self%df)) then
       deallocate(self%df)
    endif
    allocate(self%df(idf+iadd),self%eindf(idf+iadd))
    self%df(1+iadd:idf+iadd) = ar_df(1:idf)
    self%eindf(1+iadd:idf+iadd) = ar_ein(1:idf)
    do j=1,iadd
       self%eindf(j) = ar_ein(1) - 0.1*real(iadd-j+1)
       self%df(j) = ar_df(1)
    enddo

  end subroutine set_df
  !
  subroutine copy_state(state_l,state_r)
    type(state), intent(out):: state_l
    type(state), intent(in):: state_r
    integer:: i2

    state_l%stlabel = state_r%stlabel
    state_l%en = state_r%en
    state_l%enex = state_r%enex
    state_l%l = state_r%l
    state_l%m = state_r%m
    state_l%ipar = state_r%ipar
    state_l%spin = state_r%spin
    state_l%inum = state_r%inum
    state_l%v = state_r%v
    state_l%dissThresh = state_r%dissThresh
    state_l%idf = state_r%idf
    state_l%ion = state_r%ion
    state_l%resolved = state_r%resolved

    i2 = state_r%inum

    if ( associated(state_l%ein) ) deallocate(state_l%ein)
    if ( associated(state_l%cs) ) deallocate(state_l%cs)
    allocate( state_l%ein(1:i2) )
    allocate( state_l%cs(1:i2) )
    state_l%ein(1:i2) = state_r%ein(1:i2)
    state_l%cs(1:i2) = state_r%cs(1:i2)

    i2 = state_r%idf
    if ( associated(state_l%eindf) ) deallocate(state_l%eindf)
    if ( associated(state_l%df) ) deallocate(state_l%df)
    allocate( state_l%eindf(1:i2) )
    allocate( state_l%df(1:i2) )
    if (associated(state_r%eindf)) then
       state_l%eindf(1:i2) = state_r%eindf(1:i2)
    else
       deallocate(state_l%eindf)
    end if
    if (associated(state_r%df)) then
       state_l%df(1:i2) = state_r%df(1:i2)
    else
       deallocate(state_l%df)
    end if

  end subroutine copy_state
  !
  subroutine new_basis(self,n,tcstype)
    implicit none
    type(basis_state), intent(inout):: self
    integer, intent(in):: n, tcstype

    self%tcstype = tcstype
    self%n = n
    ! create array of n states
    if(n .ne. 0)  then
       allocate( self%b(n))
    endif

  end subroutine new_basis
  !
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: destruct_state
  !Purpose: deallocates all memory in the input state type
  !Date last modified: 20/08/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine destruct_state(self)
      implicit none
      type(state)::self

      if (associated(self%cs)) deallocate(self%cs)
      if (associated(self%ein)) deallocate(self%ein)
      if (associated(self%eindf)) deallocate(self%eindf)
      if (associated(self%df)) deallocate(self%df)
      if (associated(self%enEin)) deallocate(self%enEin)
     
  end subroutine destruct_state



  subroutine destruct_basis(self)
    implicit none
    type(basis_state), intent(inout):: self
    integer:: n
    integer:: i
    integer:: stat

    n= self%n
    do i=1,n
       !       print*,'dealocate: i=', i, ', size=', SIZE(self%b(i)%ein)
       if (associated(self%b(i)%cs)) deallocate(self%b(i)%cs)
       if (associated(self%b(i)%ein)) deallocate(self%b(i)%ein)
       if (associated(self%b(i)%eindf)) deallocate(self%b(i)%eindf)
       if (associated(self%b(i)%df)) deallocate(self%b(i)%df)
       if (associated(self%b(i)%enEin)) deallocate(self%b(i)%enEin)
    enddo
    deallocate(self%b, STAT=stat)

  end subroutine destruct_basis
  !
  function basis_size(self)
    implicit none
    integer:: basis_size
    type(basis_state), intent(in):: self

    basis_size = self%n
  end function basis_size
  !
  function get_ang_mom(self)
    implicit none
    integer:: get_ang_mom
    type(state), intent(in):: self

    get_ang_mom = self%l
  end function get_ang_mom
  !
  function get_ang_mom_proj(self)
    implicit none
    integer:: get_ang_mom_proj
    type(state), intent(in):: self

    get_ang_mom_proj = self%m
  end function get_ang_mom_proj
  !
  function get_energy(self)
    implicit none
    real(dp):: get_energy
    type(state), intent(in):: self

    get_energy = self%en
  end function get_energy
!  
  function get_energy_ex(self)
    implicit none
    real(dp):: get_energy_ex
    type(state), intent(in):: self
    
    get_energy_ex = self%enex
  end function get_energy_ex
  !
  function get_parity(self)
    implicit none
    integer:: get_parity
    type(state), intent(in):: self

    get_parity = self%ipar
  end function get_parity
  !
  function get_spin(self)
    implicit none
    real(dp):: get_spin
    type(state), intent(in):: self

    get_spin = self%spin
  end function get_spin
  !
  function get_cs_array(self)
    implicit none
    real(dp), pointer, dimension(:):: get_cs_array
    type(state), intent(in):: self

    get_cs_array => self%cs
  end function get_cs_array
  !
  function get_ein_array(self)
    implicit none
    real(dp), pointer, dimension(:):: get_ein_array
    type(state), intent(in):: self

    get_ein_array => self%ein
  end function get_ein_array
  !
  function get_cs_i(self,i)
    implicit none
    real(dp):: get_cs_i
    type(state), intent(in):: self
    integer, intent(in):: i 

    if( i > self%inum .or. i < 1) then
       get_cs_i = 0d0
       print*, 'WARNING: states.f90: get_cs_i(): trying to access the value out of range'
    else
       get_cs_i = self%cs(i)
    endif

  end function get_cs_i
  !
  function get_ein_i(self,i)
    implicit none
    real(dp):: get_ein_i
    type(state), intent(in):: self
    integer, intent(in):: i 

    if( i > self%inum .or. i < 1) then
       get_ein_i = 0d0
       print*, 'WARNING: states.f90: get_cs_i(): trying to access the value out of range'
    else
       get_ein_i = self%ein(i)
    endif

  end function get_ein_i
  !    
  subroutine  print_energy(self)
    use input_data
    implicit none
    type(basis_state), intent(in):: self
    integer:: n
    real(dp):: ioniz_en, exit_en 

    write(*,'("    N     M  parity    spin   exitation en.   ionization en.")')
    do n=1,self%n
       exit_en = (self%b(n)%en-self%b(1)%en)*data_in%eV  
       ioniz_en = (self%b(n)%en)*data_in%eV  
       write(*,'(i6,2i5,F5.1,2F17.5)') n, self%b(n)%m, self%b(n)%ipar, self%b(n)%spin, exit_en, ioniz_en
    enddo
    print*

  end subroutine print_energy

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: addState
  !Purpose: takes in the input state 'self' and adds it to the list of
  !         states in the statebasis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine addState(self,statebasis)
     implicit none
     type(state):: self
     type(basis_state):: statebasis
     integer:: ii, numStates, tcstype
     type(state), dimension(:), allocatable:: temp

     allocate(temp(statebasis%n))
     tcstype = statebasis%tcstype

     do ii = 1, statebasis%n
        call copy_state(temp(ii), statebasis%b(ii))
     end do
     numStates = statebasis%n

     call destruct_basis(statebasis)
     call new_basis(statebasis,numStates+1,tcstype)
     do ii = 1, numStates
        call copy_state(statebasis%b(ii), temp(ii))
        call destruct_state(temp(ii))
     end do
     call copy_state(statebasis%b(numStates+1), self)
     deallocate(temp)

  end subroutine addState




  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Shellsort algorithm
  subroutine sort_by_energy(self)
    implicit none
    type(basis_state), intent(inout):: self
    integer:: gap, i, j, N
    type(state):: Tmp

    N = self%n
    gap = N/2
    do 
       if(gap .le. 0) exit
       do i=gap+1,N
          call copy_state(Tmp,self%b(i))
          do j=i,gap+1,-gap
             if(Tmp%en .lt. self%b(j-gap)%en) then
                call copy_state(self%b(j),self%b(j-gap))
             else
                exit
             endif
             call copy_state(self%b(j-gap),Tmp)
          enddo
       enddo
       if ( gap .eq. 2) then
          gap = 1
       else
          gap = nint( real(gap)/2.2d0 )
       endif
    enddo

  end subroutine sort_by_energy
!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: sortEnergies
  !Purpose: replacement subroutine for sort_by_energy, which 
  !         doesn't sort correctly. 
  !Note: calls a recursive mergesort algorithm to organise data 
  !Date last modified: 14/01/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sortEnergies(self)
      implicit none
      type(basis_state), intent(inout)::self

      call mergeSort(self,1,self%n)
  end subroutine sortEnergies


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: mergeSort
   !Purpose: sorts statebasis data structure.
   !Date last modified: 14/01/2021
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   recursive subroutine mergeSort(self,leftInd,rightInd)
       implicit none
       type(basis_state)::self
       integer::leftInd, rightInd, midInd

       if (leftInd .lt. rightInd) then
          midInd = (leftInd + rightInd)/2 
       
          call mergeSort(self,leftInd,midInd)
          call mergeSort(self,midInd+1,rightInd)

          call mergeArray(self,leftInd,midInd,rightInd)
       end if
   end subroutine mergeSort


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Subroutine: mergeArray
   !Purpose: merges the two halves of the array
   !Date last modified: 14/01/2021
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mergeArray(self,leftInd,midInd,rightInd)
       implicit none
       type(basis_state)::self
       integer:: leftInd, rightInd, midInd
       integer:: ii, jj, kk, ll, mm
       type(state), allocatable, dimension(:)::tempArray

       allocate(tempArray(rightInd-leftInd+1))
       ii = leftInd
       jj = midInd +1
       kk = 1


       do while( (ii .le. midInd) .and. (jj .le. rightInd)) 
          if (self%b(ii)%enex .le. self%b(jj)%enex) then
             call copy_state(tempArray(kk), self%b(ii))
             ii = ii+1
          else 
             call copy_state(tempArray(kk), self%b(jj))
             jj = jj+1
          end if
          kk = kk+1
       end do
       
       do ll = ii, midInd
          call copy_state(tempArray(kk),self%b(ll))
          kk = kk+1
       end do
       do mm = jj, rightInd
          call copy_state(tempArray(kk),self%b(mm))
          kk=kk+1
       end do
       do kk = leftInd, rightInd
          call copy_state(self%b(kk),tempArray(kk-leftInd+1))
       end do

       deallocate(tempArray)

   end subroutine mergeArray




  subroutine print_state(self,filename)
    implicit none
    type(state), intent(inout):: self
    character(len=10), intent(in):: filename
    integer:: i, imax, imin

    open(15,file=filename)
    print*,'open ', filename

    write(15,*) self%stlabel
    write(15,*) "energy ionization= ", self%en
    write(15,*) "energy excitation= ",self%enex
    write(15,'("m, ipar, spin= ",2i5,F5.1)') self%m, self%ipar, self%spin 
    write(15,*) 'inum= ', self%inum
    write(15,*) 'idf= ', self%idf
    write(15,*) 'ion= ', self%ion
    write(15,*) '     i       ein             cs            eindf      diss.fraction '
    if(self%inum .eq. self%idf) then
       do i=1,self%inum
          write(15,'(i8,1P,3E15.5)') i, self%ein(i), self%cs(i), self%df(i)
       enddo
    else
       if(self%idf .eq. 0) then
          do i=1,self%inum
             write(15,'(i8,1P,2E15.5)') i, self%ein(i), self%cs(i)
          enddo
       else          
          imax = max(self%inum,self%idf)
          imin = min(self%inum,self%idf)
          do i=1,imin
             write(15,'(i8,1P,4E15.5)') i, self%ein(i), self%cs(i), self%eindf(i), self%df(i)
          enddo
          if(self%inum .eq. imax) then
             do i=imin+1,imax
                write(15,'(i8,1P,2E15.5,2A15)') i, self%ein(i), self%cs(i), "!", "!"
             enddo
          else
             do i=imin+1,imax
                write(15,'(i8,2A15,1P,2E15.5)') i,  "!", "!", self%eindf(i), self%df(i)
             enddo
          endif
       endif
    endif
    close(15)
  end subroutine print_state
!
!
! interpolate cross sections to the fine energy grid
  subroutine  intp_cs(self,Nenfine,enfine)
    implicit none
    type(basis_state), intent(in):: self
    integer, intent(in):: Nenfine
    real(dp), dimension(Nenfine), intent(in):: enfine
        
    integer:: n, Nmax,im, Nenmaxold, j
    real(dp), dimension(Nenfine):: csnew, dfnew
    real(dp), dimension(:), allocatable:: eold, csold, edf, dfold

    print*, 'Interpolate CS to a fine energy grid'
    Nmax = self%n

!   find the maxium number of energies for the old energy grid 
    Nenmaxold = 0
    do n=1,Nmax
       Nenmaxold = max(Nenmaxold,self%b(n)%inum)
    enddo
    print*,'Nenmaxold =', Nenmaxold 
    allocate(eold(Nenmaxold), csold(Nenmaxold))

    do n=1,Nmax
       im = self%b(n)%inum
       eold(1:im) = self%b(n)%ein(1:im)
       csold(1:im) = self%b(n)%cs(1:im)
!       print*, 'CS', n, im
       if(n .eq. -2) then
          do j=1,im
             print*, 'intp_cs():', n, eold(j), csold(j)
          enddo
       endif
       call INTRPL(im,eold(1:im),csold(1:im),Nenfine,enfine,csnew)
       do j=1,Nenfine ! below the first point on the old energy grid make all CS on the new grid equal to zero
          if(enfine(j) .le. eold(1)) then
             csnew(j) = 0d0
          else
             exit
          endif
       enddo
       do j=1,Nenfine ! above the last point on the old energy grid make all CS on the new grid equal to zero
          if(enfine(j)- 1e-6 .gt. eold(im)) then
             csnew(j) = 0d0
          endif
       enddo
       
!       statetoprint = 'state1'
!       call print_state(self%b(n),statetoprint)
       j = 0
       call set_cs(self%b(n),Nenfine,enfine,csnew,j)

       ! diss fractions
       im = self%b(n)%idf
       if(im .eq. 0) cycle
!       print*, 'DF', n, im
       allocate(edf(im),dfold(im))
       eold(1:im) = self%b(n)%eindf(1:im)
       dfold(1:im) = self%b(n)%df(1:im)
       do j=1,-im
          print*, 'intp_cs():', n, self%b(n)%stlabel, eold(j), dfold(j)
       enddo
       call INTRPL(im,eold(1:im),dfold(1:im),Nenfine,enfine,dfnew)
       do j=1,Nenfine ! below the first point on the old energy grid make all diss.fractions on the new grid equal to the value at the first point on the old grid
          if(enfine(j) .le. eold(1)) then
             dfnew(j) = dfold(1)
          else
             exit
          endif
       enddo
       deallocate(edf,dfold)

       call set_df(self%b(n),Nenfine,enfine,dfnew)
    enddo
    
    deallocate(eold, csold)

    !
    print*,'finish intp_cs'

  end subroutine intp_cs
!



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: scaleElCs
  !Purpose: scales elastic scattering cross sections by given
  !         percentage. Takes in processed data stored in stateBasis 
  !Date last modified: 10/06/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine scaleElCs(statebasis,percent)
      implicit none
      type(basis_state):: statebasis
      real(dp):: percent

      !Only first state in basis corresponds to elastic scattering,
      !this is the ground state.
      statebasis%b(1)%cs(:) = ((100.0d0+percent)/100.d0)*statebasis%b(1)%cs(:) 
  end subroutine scaleElCs



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: scaleInelCs
  !Purpose: scales inelastic scattering cross sections by given
  !         percentage. Takes in processed data stored in stateBasis 
  !Date last modified: 10/06/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine scaleInelCs(statebasis,percent,ii)
      implicit none
      type(basis_state):: statebasis
      real(dp):: percent
      integer:: ii
   
      !Start at 2, ignore first (ground) state corresponding to elastic scattering
      !do ii=2, statebasis%n
      !   if (statebasis%b(ii)%en .lt. 0.0d0) then
            !Scale bound state
      statebasis%b(ii)%cs(:) = ((100.0d0+percent)/100.d0)*statebasis%b(ii)%cs(:) 
      !   end if
      !end do
  end subroutine scaleInelCs



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: scaleIonCs
  !Purpose: scales ionisation cross sections by given
  !         percentage. Takes in processed data stored in stateBasis 
  !Date last modified: 10/06/2021
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine scaleIonCs(statebasis,percent)
      implicit none
      type(basis_state):: statebasis
      real(dp):: percent
      integer:: ii
   
      !Start at 2, ignore first (ground) state corresponding to elastic scattering
      do ii=2, statebasis%n
         if (statebasis%b(ii)%en .gt. 0.0d0) then
            !Scale bound state
            statebasis%b(ii)%cs(:) = ((100.0d0+percent)/100.d0)*statebasis%b(ii)%cs(:) 
         end if
      end do
  end subroutine scaleIonCs


end module state_class





subroutine placeen(ix,xv, yv, n, x, ar)
  implicit none
  integer, parameter:: dp =selected_real_kind(15,307)
  real(dp), intent(in):: xv,yv
  integer, intent(in):: n
  real(dp), dimension(n), intent(inout):: x,ar
  integer:: i, ip, ix
  real(dp), dimension(n):: tmpar

  do i=1,n
     if(xv .lt. x(i)) exit
  enddo
  ip = i-1  ! where to place xv
!  print*, '--->>>> ip =', ip
  if(ip .eq. 0 ) return ! all done

  tmpar(1:ix-1) = x(1:ix-1)
  tmpar(ix:ip-1) = x(ix+1:ip)
  tmpar(ip) = xv
  tmpar(ip+1:n) = x(ip+1:n)

  x(:) = tmpar(:)

  tmpar(1:ix-1) = ar(1:ix-1)
  tmpar(ix:ip-1) = ar(ix+1:ip)
  tmpar(ip) = yv
  tmpar(ip+1:n) = ar(ip+1:n)

  ar(:) = tmpar(:)



end subroutine placeen
