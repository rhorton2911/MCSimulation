module sturmian_class
 
  integer, parameter:: Nonrel=0, Rel=1

  public::  new_basis, copy_basis, init_function, copy, multiply_byscalar, construct, construct_spheroidal, destruct, combine_basis, basis_size, get_k, get_ang_mom, get_ang_mom_proj, get_max_L, value, get_minf, get_maxf, fpointer, me_1el, ovlp3_nr, get_energy, print_energy, sort_by_energy, print_wf, ovlp, set_k, set_ang_mom_prog, get_alpha, integral_oidJ
  
  private:: construct_nr, destruct_nr, value_nr
!  
  type, public:: sturmian_nr
     private
     real(dp):: en  ! enrgy of one-electron state, for Sturmians en is set to zero: en=0.0
     integer:: l  ! angular momentum of sturmian function
     integer:: m  ! projection m  
     integer:: k  ! order of sturmian function, or principal q. number, or redundant for one-e.states 
     integer:: minf, maxf
     real(dp), pointer, dimension(:) :: f => NULL(), g => NULL()
     real(dp):: alpha

  end type sturmian_nr 
!  
!  
  type, public:: basis_sturmian_nr
!     private
     type(sturmian_nr), pointer, dimension(:) :: b  => NULL()
     integer:: n
     real(dp), pointer, dimension(:,:) :: ortint => NULL(), ham1el => NULL()
  end type basis_sturmian_nr
!  
!
  interface init_function
     module procedure init_function_nr,  init_function_nr_m, init_function_nr_st, init_function_nr_st_m, init_function_oid
  end interface
!
  interface copy
     module procedure copy_nr
  end interface
!
  interface multiply_byscalar
     module procedure multiply_byscalar_nr
  end interface
!
interface new_basis
     module procedure new_basis_nr
  end interface
!     
interface new_basis_cont
     module procedure  new_basis_proj
  end interface
!
interface copy_basis
     module procedure copy_basis_nr
  end interface
!     
  interface construct
     module procedure construct_nr,  construct_wflocalpot_nr, construct_all_nr, construct_all_alpha_nr
  end interface
!
  interface construct_spheroidal
     module procedure construct_spheroidal_basis
  end interface
!  
  interface destruct
     module procedure destruct_nr,  destruct_nr_bf
  end interface
!
  interface combine_basis
     module procedure combine_basis_nr
  end interface
!
  interface basis_size
     module procedure basis_size_nr
  end interface
!
  interface value
     module procedure  value_nr
  end interface
!
  interface get_energy
     module procedure  get_energy_nr
  end interface
!
  interface get_ang_mom
     module procedure get_ang_mom_nr
  end interface
!
  interface set_ang_mom
     module procedure set_ang_mom_nr
  end interface
!
  interface get_ang_mom_proj
     module procedure get_ang_mom_proj_nr
  end interface
!
  interface set_ang_mom_proj
     module procedure set_ang_mom_proj_nr
  end interface
!
  interface get_max_L
     module procedure get_max_L
  end interface
!
  interface get_k
     module procedure get_k_nr
  end interface
!
  interface set_k
     module procedure set_k_nr
  end interface
!
  interface fpointer
     module procedure  fpointer_nr
  end interface
!
  interface gpointer
     module procedure  gpointer_nr
  end interface
!
  interface get_minf
     module procedure  get_minf_nr
  end interface
!
  interface get_maxf
     module procedure  get_maxf_nr
  end interface
!
  interface print_energy
     module procedure  print_energy_nr
  end interface
!
  interface print_wf
     module procedure  print_wf_nr
  end interface
!
  interface sort_by_energy
     module procedure  sort_by_energy_nr
  end interface
!
!
  interface gsort
     module procedure gsort_nr
  end interface
!
  interface me_1el
     module procedure  me_1el_nr
  end interface
!
!
  interface ovlp
     module procedure  ovlp_nr
  end interface
!
  interface write_sturmian
     module procedure write_sturmians
  end interface
!  
contains
  !
  subroutine init_function_nr(self,l,k,i1,i2,temp,nr,alpha)
    type(sturmian_nr), intent(inout):: self
    integer, intent(in):: l, k, i1, i2, nr
    real(dp), dimension(nr), intent(in):: temp
    real(dp):: alpha

    self%en = 0d0
    self%l = l
    self%m = 0
    self%k = k
    self%minf = i1
    self%maxf = i2
    ! create array for one-electron function
    allocate(self%f(i2))
    self%f(1:i1) = 0.0d0
    self%f(i1:i2) = temp(i1:i2)
    self%alpha = 0d0
    self%alpha = alpha
  end subroutine init_function_nr
!
  subroutine init_function_nr_m(self,l,m,k,alpha,i1,i2,temp,nr)
    implicit none
    type(sturmian_nr), intent(inout):: self
    integer, intent(in):: l, m, k, i1, i2, nr
    real(dp), intent(in) :: alpha
    real(dp), dimension(nr), intent(in):: temp

    self%en = 0d0
    self%l = l
    self%m = m
    self%k = k
    self%alpha = alpha
    self%minf = i1
    self%maxf = i2
    ! create array for one-electron function
    allocate(self%f(i2))
    self%f(1:i1) = 0.0d0
    self%f(i1:i2) = temp(i1:i2)
  end subroutine init_function_nr_m
!
  subroutine init_function_nr_st(self,en,l,k,i1,i2,temp,nr)
    type(sturmian_nr), intent(inout):: self
    real(dp), intent(in):: en
    integer, intent(in):: l, k, i1, i2, nr
    real(dp), dimension(nr), intent(in):: temp

    self%en = en
    self%l = l
    self%m = 0
    self%k = k
    self%minf = i1
    self%maxf = i2
    ! create array for one-electron function
    allocate(self%f(i2))
    self%f(1:i1) = 0.0d0
    self%f(i1:i2) = temp(i1:i2)
  end subroutine init_function_nr_st
!
  subroutine init_function_nr_st_m(self,en,l,m,k,i1,i2,temp,nr)
    type(sturmian_nr), intent(inout):: self
    real(dp), intent(in):: en
    integer, intent(in):: l, m, k, i1, i2, nr
    real(dp), dimension(nr), intent(in):: temp

    self%en = en
    self%l = l
    self%m = m
    self%k = k
    self%minf = i1
    self%maxf = i2
    ! create array for one-electron function
    allocate(self%f(i2))
    self%f(1:i1) = 0.0d0
    self%f(i1:i2) = temp(i1:i2)
  end subroutine init_function_nr_st_m
!
  subroutine init_function_oid(self, k,l,m, alpha, i1,i2,fVec,gVec)
    type(sturmian_nr), intent(inout):: self
    integer, intent(in) :: k,l,m, i1,i2
    real(dp), intent(in) :: alpha
    real(dp), dimension(i1:i2), intent(in) :: fVec, gVec

    self%k = k ; self%l = l ; self%m = m
    self%alpha = alpha
    self%minf = i1 ; self%maxf = i2
    allocate( self%f(1:i2), self%g(1:i2) )
    self%f(1:i1) = 0d0 ; self%g(1:i1) = 0d0
    self%f(i1:i2) = fVec(:) ; self%g(i1:i2) = gVec(:)
  end subroutine init_function_oid
!
  subroutine copy_nr(state_l,state_r)
    type(sturmian_nr), intent(out):: state_l
    type(sturmian_nr), intent(in):: state_r
    integer:: i1,i2

    state_l%en = state_r%en
    state_l%l = state_r%l
    state_l%m = state_r%m
    state_l%k = state_r%k
    state_l%alpha = state_r%alpha 

    i1 = state_r%minf
    i2 = state_r%maxf
    state_l%minf = i1
    state_l%maxf = i2

    if ( associated(state_l%f) ) deallocate(state_l%f)
    allocate( state_l%f(1:i2) )
    state_l%f(1:i1) = 0d0
    state_l%f(i1:i2) = state_r%f(i1:i2)

    if ( associated(state_l%g) ) deallocate(state_l%g)
    if ( associated(state_r%g) ) then
       allocate( state_l%g(1:i2) )
       state_l%g(1:i1) = 0d0
       state_l%g(i1:i2) = state_r%g(i1:i2)
    end if

  end subroutine copy_nr
!
!
  subroutine multiply_byscalar_nr(scalar,state)
    real(dp), intent(in):: scalar
    type(sturmian_nr), intent(inout):: state
    integer:: i1,i2

    i1 = state%minf
    i2 = state%maxf
    if(associated(state%f)) then
       state%f(i1:i2) = scalar*state%f(i1:i2) 
    endif
  end subroutine multiply_byscalar_nr
!
!
  subroutine new_basis_nr(self,n)
    implicit none

    type(basis_sturmian_nr), intent(inout):: self
    integer, intent(in):: n

    self%n = n
    ! create array of n one-electron functions
    if(n .ne. 0)  then
       allocate( self%b(n), self%ortint(n,n), self%ham1el(n,n) )
       self%ortint(:,:) = 0.0d0
       self%ham1el(:,:) = 0d0
    endif

  end subroutine new_basis_nr
!
!
  subroutine new_basis_proj(self,n)
    implicit none

    type(basis_sturmian_nr), intent(inout):: self
    integer, intent(in):: n
    

    self%n = n
    ! create array of n one-electron functions
    if(n .ne. 0)  then
       allocate( self%b(n) )
    end if

  end subroutine new_basis_proj
!
!
  subroutine copy_basis_nr(basis_l,basis_r)
    use MPI_module
    type(basis_sturmian_nr), intent(inout):: basis_l,basis_r
    integer :: n

    if(basis_l%n .ne. basis_r%n) then
       if (myid == 0) print*, 'sturmian.f90: Can not copy basis of unequal size'
       stop
    endif

    if(basis_l%n .eq. 0) then
       if (myid == 0) print*, 'sturmian.f90: basis size is zero'
       stop
    endif
 
    n = basis_l%n
    do i = 1, n
       call copy(basis_l%b(i), basis_r%b(i))
    end do
    basis_l%ortint(:,:) =  basis_r%ortint(:,:)
    basis_l%ham1el(:,:) = basis_r%ham1el(:,:)

  end subroutine copy_basis_nr
!
  subroutine combine_basis_nr(self,basis_F5,basis_data,basis_type)
!!$ Mark: This subroutine is used to combine two bases together. 
!!$ Will be used in the constructions of a hybrid bases for H2
!!$ H2+ built from Sturmians in the data.in file H2 A.O.
!!$ built from Sturmians in F5 file
    use MPI_module
    use grid_radial
    implicit none

    type(basis_sturmian_nr), intent(inout) :: self
    type(basis_sturmian_nr), intent(inout) :: basis_F5,basis_data
    integer, intent(in) :: basis_type

    type(sturmian_nr), pointer :: bi,bj
    integer :: i,j, ki,kj, li,lj, r1,r2, Nmax,Nmax_F5,Nmax_data
    real(dp) :: alphai,alphaj
    real(dp), pointer, dimension(:):: fi, fj
    
    Nmax_F5 = basis_size(basis_F5)
    Nmax_data = basis_size(basis_data)
!!$ Allocate a basis of size Nmax_F5 + Nmax_data 
    Nmax = Nmax_F5 + Nmax_data 
    call new_basis_nr(self,Nmax)

!!$ Copy basis functions. basis_F5 first.
    do i = 1, Nmax_F5
       call copy(self%b(i), basis_F5%b(i))
    end do
    do i = 1, Nmax_data
       j = i + Nmax_F5
       call copy(self%b(j), basis_data%b(i))
    end do

!$$ Overlaps of new basis. Doing from scratch. No symmetry in array.
    self%ortint(:,:) = 0d0
    self%ham1el(:,:) = 0d0

    do i=1, Nmax
       bi => self%b(i)
       ki = get_k(bi)
       fi => fpointer(bi)
       li =  get_ang_mom(bi)
       alphai = get_alpha(bi)
       
       do j = 1, Nmax
          bj => self%b(j)
          kj = get_k(bj)
          fj => fpointer(bj)
          lj = get_ang_mom(bj)
          alphaj = get_alpha(bj)
             
          if ( li /= lj ) cycle ! Selection rules
!$$ Obtaining analytic overlaps for functions with the same exponential fall offs 
          if ( alphai == alphaj ) then
             
             if ( ki == kj  ) then
                self%ortint(i,j) = 1d0
             else if (kj + 1  == ki ) then  ! Need to check this Mark
                self%ortint(i,j) =  -0.5d0*sqrt(1.0-dble(lj*(lj+1))/dble((lj+kj)*(lj+kj+1)))
             else if (ki + 1  == kj ) then  
                self%ortint(i,j) =  -0.5d0*sqrt(1.0-dble(li*(li+1))/dble((li+ki)*(li+ki+1)))
             end if

!$$ Overlaps obtained numerically for functions with different exponential fall offs            
          else
             r1 = max( get_minf(bi), get_minf(bj) )
             r2 = min( get_maxf(bi), get_maxf(bj) )
             
             self%ortint(i,j) = SUM( fi(r1:r2) * fj(r1:r2) * grid%weight(r1:r2) )
          end if
             
       end do
    end do

  end subroutine combine_basis_nr
!
!
! create basis of n Sturmian functions and allocate space for each radial part of the function in the basis
  subroutine construct_nr(self,n,l,al)
    use grid_radial
    use MPI_module

    implicit none
    type(basis_sturmian_nr), intent(inout):: self
    integer, intent(in):: n
    integer, intent(in):: l
    real(dp), intent(in):: al

    integer:: i, k, i1, i2, la
    real(dp):: f8(grid%nr,n), x, tmp
    real(dp), dimension(grid%nr):: temp   
    real(dp):: lambda
    real(dp), pointer, dimension(:):: weight, gridr
    integer:: ortog, m
    
    m = 0   !   NOTE to be corrected...

!
    weight => grid%weight
!     weight => grid%bool
    gridr => grid%gridr

! create array of n one-electron functions for given l
    call new_basis(self,n)
    lambda = 2d0 * al
    ortog = 0
    if(ortog .eq. 1) then
       call  lagpol8(dble(2*l+2),lambda,f8,grid%nr,n,gridr,grid%nr)
    else
       call  lagpol8(dble(2*l+1),dble(2.0d0*al),f8,grid%nr,n,gridr,grid%nr)
    endif

    do k=1,n
       if(ortog .eq. 1) then
          tmp = k
          do i=1,2*l+1
             tmp = tmp*(k+i)
          enddo
          tmp = sqrt(lambda/tmp)
       else
          tmp =  2.0d0*(k+l)
          do i=0,2*l
             tmp = tmp*(k+i)
          enddo
          tmp = sqrt(tmp)
          tmp = sqrt(2d0*al)/tmp
       endif

       do i=1,grid%nr
          x = gridr(i) * 2d0*al
          temp(i) = tmp*(x**(l+1))*exp(-x/2.0D0)*f8(i,k)
       end do


       call minmaxi(temp,grid%nr,i1,i2)
       call init_function(self%b(k),l,k,i1,i2,temp,grid%nr,al)
!      if (myid == 0) print*, 'k=',k,', size=',SIZE(self%b(k)%f)
    enddo

!!$ make overlap matrix
    la = l
    do k=1,n
       self%ortint(k,k) = 1d0
    enddo
    do k=1,n-1
          self%ortint(k,k+1) =  -0.5d0*sqrt(1.0-dble(la*(la+1))/dble((la+k)*(la+k+1)))
          self%ortint(k+1,k) = self%ortint(k,k+1) 
    enddo

    if (myid == 0) write(*,'("Created nonrel. Sturmian basis with n=",i5,"  functions for angular momentum l=",i5)') n, l
  end subroutine construct_nr
!
!-------------------------------------------------------------------------------------------------!
! Make a radial basis of all Sturmians.
! Functions are in the form N * exp(-x/2) * x^a/2 * L_n^k(x). If a = k the basis is orthogonal.
!
! INPUT: self (type basis_sturmian_nr)
!
! OUTPUT: self
!
!
  subroutine construct_all_nr(self, dataARG)

    use grid_radial
    use input_data
    use MPI_module

    implicit none

    type(basis_sturmian_nr), intent(inout):: self
    type(input), intent(in):: dataARG  

    integer:: i, j, k, i1, i2, kn, nall, l, labot, latop, n, la, nsp, ni, nj, li, lj, mi, mj
    real(dp):: f8(grid%nr,MAXVAL(dataARG%nps(dataARG%labot:dataARG%latop))), x, tmp
    real(dp), dimension(grid%nr):: temp   
    real(dp):: lambda, lagarg
    real(dp), pointer, dimension(:):: weight, gridr
    integer:: ortog, m

    ortog = dataARG%calculation_type   ! Orthogonality of the basis functions. 0: nonorthogonal (2l+1 type); 1: orthogonal (2l+2 type).

    weight => grid%weight
    gridr => grid%gridr

    ! create array of nall one-electron functions for given all l.
    labot = dataARG%labot
    latop = dataARG%latop

    nall = SUM(dataARG%nps(labot:latop))
    call new_basis(self,nall)

    kn = 0   ! Counter of functions calculated.

    do l=labot,latop
       
       n = dataARG%nps(l)
       lambda = 2d0 * dataARG%alpha(l)
       
       ! Calculate the Laguerre polynomials up to order n on the given r-grid.
       f8(:,:) = 0d0
       if(ortog .eq. 0) then   ! Nonorthogonal spherical.
          lagarg = dble(2*l+1)
       elseif(ortog .eq. 1) then   ! Orthogonal spherical.
          lagarg = dble(2*l+2)
       endif
       call lagpol8( lagarg, lambda, f8, grid%nr, n, gridr, grid%nr )
       
       ! Determine the normalisation constant.
       do k=1,n
          tmp = 1d0   ! Orthogonal case.
          if(ortog .eq. 0) tmp = 2d0*(k+l)   ! Nonorthogonal case.

          do i=0, nint(lagarg)-1   ! Factorial quotient (k-1+lagarg)!/(k-1)!
             tmp = tmp*(k+i)
          enddo
          tmp = sqrt( lambda / tmp )

          ! Multiply the Laguerre polynomials by the normalisation, exponential and polynomial factors.
          do i=1,grid%nr
             x = gridr(i) * lambda
             temp(i) = tmp * x**(l+1d0) * exp(-x/2.0d0) * f8(i,k)
          end do
          call minmaxi(temp,grid%nr,i1,i2)
          
          kn = kn + 1
          
          ! Initialise each Sturmian function (type sturmian_nr) within the basis.
          call init_function(self%b(kn),l,k,i1,i2,temp,grid%nr,lambda/2d0)
!         if (myid == 0) print*, 'k=',k,', size=',SIZE(self%b(k)%f)

       enddo ! k -- each loop creates a single basis function.
       if (myid == 0) write(*,'(2(A,I2))') 'Created radial basis with n=',n, ' Sturmian functions for angular momentum=',l
    enddo ! l -- one loop for each labot<=l<=latop.

    ! Check if all requested functions have been created.
    if(nall .ne. kn) then
       if (myid == 0) print*,'sturmian.f90: nall != kn:', nall,  kn
       stop
    endif

    self%n = nall

    ! Overlap matrix.
    do i=1,nall   ! nallxnall identity matrix.
       self%ortint(i,i) = 1d0
    enddo
    if(ortog .eq. 0) then   ! Nonorthogonal--off-diagonal overlap matrix elements.
       do nsp=1,nall-1
          la = get_ang_mom(self%b(nsp))
          i = get_k(self%b(nsp))
          if( la .eq.  get_ang_mom(self%b(nsp+1))) then
             self%ortint(nsp,nsp+1) = -0.5d0*sqrt(1.0-dble(la*(la+1))/dble((la+i)*(la+i+1)))
             self%ortint(nsp+1,nsp) = self%ortint(nsp,nsp+1) 
          endif
       enddo
    endif


    if (myid == 0) write(*,'(3(A,I3))') 'Created radial basis with n=',nall, ' Sturmian functions for angular momentum from ',labot, ' to ',latop       


  end subroutine construct_all_nr
!
!
  subroutine construct_all_alpha_nr(self, dataARG, alpha_nl)
!!$ This routine is exactle the same as the one above, however, it 
!!$ constructs Laguerre functions with varying alpha designated by the 
!!$ array alpha_nl(self%nps(labot:latop),labot:latop)
    use grid_radial
    use input_data
    use MPI_module

    implicit none

    type(basis_sturmian_nr), intent(inout):: self
    type(input), intent(in):: dataARG
    real(dp), intent(in):: alpha_nl(MAXVAL(dataARG%nps(dataARG%labot:dataARG%latop)),dataARG%labot:dataARG%latop) 

    integer:: i, j, k, i1, i2, kn, nall, l, labot, latop, n, la, nsp, ni, nj, li, lj, mi, mj
    real(dp):: f8(grid%nr,MAXVAL(dataARG%nps(dataARG%labot:dataARG%latop))), x, tmp
    real(dp), dimension(grid%nr):: temp
    real(dp):: lambda, lagarg
    real(dp), pointer, dimension(:):: weight, gridr
    integer:: ortog, m
    type(sturmian_nr), pointer :: bi, bj
    integer :: ki, kj 
    real(dp):: alphai, alphaj
    real(dp), pointer, dimension(:):: fi, fj

    ortog = dataARG%calculation_type   ! Orthogonality of the basis functions. 0: nonorthogonal (2l+1 type); 1: orthogonal (2l+2 type).

    weight => grid%weight
    gridr => grid%gridr

    ! create array of nall one-electron functions for given all l.
    labot = dataARG%labot
    latop = dataARG%latop

    nall = SUM(dataARG%nps(labot:latop))
    call new_basis(self,nall)

    kn = 0   ! Counter of functions calculated.

    do l=labot,latop

       n = dataARG%nps(l)

       ! Calculate the Laguerre polynomials up to order n on the given r-grid.
       f8(:,:) = 0d0
       if(ortog .eq. 0) then   ! Nonorthogonal spherical.
          lagarg = dble(2*l+1)
       elseif(ortog .eq. 1) then   ! Orthogonal spherical.
          lagarg = dble(2*l+2)
       endif

       ! Determine the normalisation constant.
       do k=1,n
          lambda = 2d0 * dataARG%alpha_nl(k,l)
          call lagpol8( lagarg, lambda, f8, grid%nr, n, gridr, grid%nr )

          if(ortog .eq. 0) then   ! Nonorthogonal--extra factor of 1/(2k+2l).
             tmp = 2.0d0*(k+l)
          elseif(ortog .eq. 1) then   ! Orthogonal.
             tmp = 1.0d0
          endif
          do i=0, nint(lagarg)-1   ! Factorial quotient (k-1+lagarg)!/(k-1)!
             tmp = tmp*(k+i)
          enddo
          tmp = sqrt( lambda / tmp )

          ! Multiply the Laguerre polynomials by the normalisation, exponential
          ! and polynomial factors.
          if(ortog .eq. 0) then   ! Nonorthogonal--correction to the power of x.
             lagarg = lagarg + 1.0d0   ! 2l+1 -> 2(l+1)
          endif
          do i=1,grid%nr
             x = gridr(i) * lambda
             temp(i) = tmp * x**(lagarg/2.0d0) * exp(-x/2.0d0) * f8(i,k)
          end do
          if( ortog .eq. 0 ) then   ! Uncorrection correction.
             lagarg = lagarg - 1.0d0   ! 2(l+1) -> 2l+1
          endif

          call minmaxi(temp,grid%nr,i1,i2)

          kn = kn + 1

          ! Initialise each Sturmian function (type sturmian_nr) within the
          ! basis.
          call init_function(self%b(kn),l,k,i1,i2,temp,grid%nr,lambda/2d0)
!         if (myid == 0) print*, 'k=',k,', size=',SIZE(self%b(k)%f)

       enddo ! k -- each loop creates a single basis function.
       if (myid == 0) write(*,'(2(A,I2))') 'Created radial basis with n=',n, ' Sturmian functions for angular momentum=',l
    enddo ! l -- one loop for each labot<=l<=latop.

    ! Check if all requested functions have been created.
    if(nall .ne. kn) then
       if (myid == 0) print*,'sturmian.f90: nall != kn:', nall,  kn
       stop
    endif

    self%n = nall

    ! Overlap matrix.
    do i=1,nall   ! nallxnall identity matrix.
       self%ortint(i,i) = 1d0
    enddo

    do i = 1, basis_size(self) 
       bi => self%b(i)
       ki = get_k(bi)
       fi => fpointer(bi)
       li =  get_ang_mom(bi)
       alphai = get_alpha(bi)
       
       do j = 1, basis_size(self)
          bj => self%b(j)
          kj = get_k(bj)
          fj => fpointer(bj)
          lj =  get_ang_mom(bj)
          alphaj = get_alpha(bj)

          if ( li /= lj ) cycle ! Selection rules
!$$ Obtaining analytic overlaps for functions with the same exponential fall-offs 
          if (  alphai == alphaj ) then
             if ( ki == kj  ) then
                self%ortint(i,j) = 1d0
             else if (kj + 1  == ki ) then  ! Need to check this Mark
                self%ortint(i,j) = -0.5d0*sqrt(1.0-dble(lj*(lj+1))/dble((lj+kj)*(lj+kj+1)))
             else if (ki + 1  == kj ) then
                self%ortint(i,j) = -0.5d0*sqrt(1.0-dble(li*(li+1))/dble((li+ki)*(li+ki+1)))
             end if
!$$ Overlaps obtained numerically for functions with different exponential fall-offs            
          else
             i1 = max( get_minf(bi), get_minf(bj) )
             i2 = min( get_maxf(bi), get_maxf(bj) )
             self%ortint(i,j) = SUM( fi(i1:i2)*fj(i1:i2) * grid%weight(i1:i2) )
          end if
        end do ! j
    end do ! i

    ! Code for verifying the overlap matrix if one wishes to do so.
    if( 2 + 2 == 5 ) then
       do n=1,nall
          do i=1,n
             call ovlp_nr( self%b(n), self%b(i), tmp )
             if (myid == 0) write(*,'(2I5, 2F20.10)') n, i, tmp, self%ortint(n,i)
          enddo
          if (myid == 0) print*
       enddo
    endif

    if (myid == 0) write(*,'(3(A,I3))') 'Created radial basis with n=',nall, ' Sturmian functions for angular momentum from ',labot, ' to ',latop


  end subroutine construct_all_alpha_nr
!
!-------------------------------------------------------------------------------------------------!
!
  subroutine construct_spheroidal_basis( self, basis, Vi2,Vvec )
    use grid_radial
    use input_data
    use MPI_module
    implicit none

    type(basis_sturmian_nr), intent(inout) :: self
    type(basis_input), intent(in) :: basis
    integer, optional :: Vi2
    real(dp), dimension(:), optional :: Vvec

    integer :: nTot, lMin,lMax, mMin,mMax, k,l,m, ind, i,j
    real(dp) :: R,alpha
    real(dp), dimension(:), pointer :: rhoVec,weightVec, fnVec
    real(dp), dimension(1:grid%nr) :: leftVec,rightVec


    lMin = basis%labot ; lMax = basis%latop
    mMin = basis%mabot ; mMax = basis%matop

    nTot = 0
    do m = mMin,mMax
!!$       l = lMax - m + 1   ! Number of l's possible for this m.
!!$       k = basis%nps(m)   ! Number of k's chosen for this m.
!!$       nTot = nTot + k*l
       k = sum( basis%nps(m:basis%latop) )
       nTot = nTot + k
    end do
    call new_basis( self, nTot )

    ind = 0   ! The subroutine will start loading functions into ind+1.
    do m = mMin,mMax
       alpha = basis%alpha(m)
       call construct_spheroidal_m( self,ind, basis%nps,lMax,m, alpha )
    end do
    if (myid==0) write(*,'(4(A,I4))') 'Created a spheroidal basis of ', ind, ' (requested ', nTot, ') functions with m =', mMin, ' -->', mMax, ' and m <= l <=', lMax
    if (ind /= nTot) stop 'sturmian.f90/construct_spheroidal_basis() : ind /= nTot!'

  end subroutine construct_spheroidal_basis

  subroutine construct_spheroidal_m( self,ind, kVec,lMax,m, alpha )
    use grid_radial
    implicit none

    type(basis_sturmian_nr), intent(inout):: self
    integer, intent(inout) :: ind
    integer, dimension(0:lMax) :: kVec
    integer, intent(in):: lMax, m
    real(dp), intent(in):: alpha

    integer :: nPts, k,kMax, l, i, i1,i2, j1,j2
    real(dp) :: lambda, rk,rm, norma
    real(dp), dimension(1:grid%nr) :: xVec, fVec,gVec
    real(dp), dimension(:), pointer :: rhoVec
    real(dp), dimension(1:grid%nr,1:maxval(kVec(:))) :: f8


    ! Set up and calculate the associated Laguerre functions.
    kMax = maxval( kVec(:) )
    rm = dble(m)
    lambda = 2d0 * alpha
    nPts = grid%nr
    rhoVec => grid%gridr
    xVec(:) = exp(-lambda*rhoVec(:)/2d0) * (lambda*rhoVec(:))**(dble(m)/2d0)
    call lagpol8( dble(m),lambda, f8,nPts,kMax, rhoVec(:),nPts )

    ! Calculate radial functions and store them in the basis.
    do l = m, lMax
       do k = 1, kVec(l)
!!$       do k = 1, kVec(m)
          rk = dble(k); norma = lambda
          do i = k, k-1+m
             norma = norma / dble(i)
          end do
          norma = dsqrt(norma)   ! sqrt( lambda * (k-1)!/(k-1+m)! )

          ! Derivative from recursion relation requires previous function.
          gVec(:) = 0d0
!!$          if (k > 1) gVec(i1:i2) = -dsqrt((rk-1d0)*(rk-1d0+rm)) / rhoVec(i1:i2) * fVec(i1:i2)
          if (k > 1) gVec(i1:i2) = -dsqrt((rk-1d0)*(rk-1d0+rm)) * fVec(i1:i2)

          fVec(:) = norma * xVec(:) * f8(:,k)
          gVec(:) = (gVec(:) + (rk-1d0+rm/2d0)*fVec(:)) /rhoVec(:)
          gVec(:) = gVec(:) - alpha*fVec(:)
          call minmaxi(fVec,nPts, i1,i2)
          call minmaxi(gVec,nPts, j1,j2)
          i1 = min(i1,j1); i2 = max(i2,j2)             

          ind = ind+1
          call init_function( self%b(ind), k,l,m, alpha, i1,i2,fVec(1:i2),gVec(1:i2) )
       end do
    enddo


  end subroutine construct_spheroidal_m
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!
  subroutine destruct_nr_bf(self)
    implicit none
    type(sturmian_nr), intent(inout):: self
    integer:: stat

    deallocate(self%f)
  end subroutine destruct_nr_bf
!
!
  subroutine destruct_nr(self)
    implicit none
    type(basis_sturmian_nr), intent(inout):: self
    integer:: n
    integer:: i
    integer:: stat

    n= self%n
    do i=1,n
!       print*,'dealocate: i=', i, ', size=', SIZE(self%b(i)%f)
       if (associated(self%b(i)%f)) deallocate(self%b(i)%f)
    enddo
    deallocate(self%b, STAT=stat)

    if (associated(self%ortint)) deallocate(self%ortint)
    if (associated(self%ham1el)) deallocate(self%ham1el)

  end subroutine destruct_nr
!
!
  function basis_size_nr(self)
    implicit none
    integer:: basis_size_nr
    type(basis_sturmian_nr), intent(in):: self
    
    basis_size_nr = self%n
  end function basis_size_nr
!
  function get_num_m0(self)
    implicit none
    type(basis_sturmian_nr), intent(in) :: self
    integer :: get_num_m0
    integer :: i
    
    get_num_m0 = 0
    do i = 1, self%n
       if (self%b(i)%m == 0) get_num_m0 = get_num_m0 + 1
    end do
  end function get_num_m0
!
!
  function value_nr(self,i)
    implicit none
    real(dp):: value_nr
    type(sturmian_nr), intent(in):: self
    integer, intent(in):: i  !

    if ( i == 0 ) then
       value_nr = self%f(self%minf)
    elseif( i<self%minf .or. self%maxf<i) then
       value_nr = 0d0
    else
       value_nr = self%f(i)
    endif
    
  end function value_nr
!
!
  function get_k_nr(self)
    implicit none
    integer:: get_k_nr
    type(sturmian_nr), intent(in):: self
    
    get_k_nr = self%k
  end function get_k_nr
!
  subroutine set_k_nr(self,k)
    implicit none
    type(sturmian_nr):: self
    integer, intent(in):: k

    self%k = k
  end subroutine set_k_nr
!
!
  function get_ang_mom_nr(self)
    implicit none
    integer:: get_ang_mom_nr
    type(sturmian_nr), intent(in):: self
    
    get_ang_mom_nr = self%l
  end function get_ang_mom_nr
!
  subroutine set_ang_mom_nr(self,L)
    implicit none
    type(sturmian_nr):: self
    integer, intent(in):: L
    
    self%l = L
  end subroutine set_ang_mom_nr
!
!
  function get_ang_mom_proj_nr(self)
    implicit none
    integer:: get_ang_mom_proj_nr
    type(sturmian_nr), intent(in):: self
    
    get_ang_mom_proj_nr = self%m
  end function get_ang_mom_proj_nr
!
  subroutine set_ang_mom_proj_nr(self,M)
    implicit none
    type(sturmian_nr):: self
    integer, intent(in):: M
    
    self%m = M
  end subroutine set_ang_mom_proj_nr
!
!
function   get_min_L(self)
  implicit none
  integer:: get_min_L
  type(basis_sturmian_nr), intent(in):: self
  integer:: n, l, ltmp, Nmax
  
  Nmax = basis_size(self)
  if(Nmax .le. 0) then
     print*,"sturmian.f90:  Nmax <=0 in get_min_L(self)"
     stop
  endif

  l = get_ang_mom(self%b(1))
  do n=2,Nmax
     ltmp = get_ang_mom(self%b(n))
     l = min(l,ltmp)
  enddo

  get_min_L = l
  
  return
end function get_min_L
!
function   get_max_L(self)
  implicit none
  integer:: get_max_L
  type(basis_sturmian_nr), intent(in):: self
  integer:: n, l, ltmp, Nmax
  
  Nmax = basis_size(self)
  if(Nmax .le. 0) then
     print*,"sturmian.f90:  Nmax <=0 in get_max_L(self)"
     stop
  endif

  l = get_ang_mom(self%b(1))
  do n=2,Nmax
     ltmp = get_ang_mom(self%b(n))
     l = max(l,ltmp)
  enddo

  get_max_l = l
  
  return
end function get_max_L
!
!  Finds how many orbitals with given L are in the basis
function   get_max_nL(self,l)
  implicit none
  integer:: get_max_nL
  type(basis_sturmian_nr), intent(in):: self
  integer, intent(in):: l
  integer:: n, ltmp, Nmax, numk
  
  Nmax = basis_size(self)
  if(Nmax .le. 0) then
     print*,"sturmian.f90:  Nmax <=0 in get_maxnL(self,l)"
     stop
  endif

  numk = 0
  do n=1,Nmax
     ltmp = get_ang_mom(self%b(n))
     if(l .eq. ltmp) then
        numk = numk + 1
     endif
  enddo

  get_max_nL = numk
  
  return
end function get_max_nL
!
!  Finds how many orbitals with given L are in the basis up to a basis function M
function   get_nL_uptoM(self,l,M)
  implicit none
  integer:: get_nL_uptoM
  type(basis_sturmian_nr), intent(in):: self
  integer, intent(in):: l, M
  integer:: n, ltmp, Nmax, numk
  
  Nmax = basis_size(self)
  if(Nmax .le. 0) then
     print*,"sturmian.f90:  Nmax <=0 in get_nL_uptoM(self)"
     stop
  endif

  if(M .le. 0) then
     print*,"sturmian.f90:  M <=0 in get_min_kappa_uptoM(self,kappa,M)"
     stop
  endif

  if(M .gt. Nmax) then
     print*,"sturmian.f90:  M <= Nmax in get_nL_upto(self,kappa,M)"
     stop
  endif


  numk = 0
  do n=1,M
     ltmp = get_ang_mom(self%b(n))
     if(l .eq. ltmp) then
        numk = numk + 1
     endif
  enddo

  get_nL_uptoM = numk
  
  return
end function get_nL_uptoM
!
!  Finds the maximum value of orbitals for given L over all possible L are in the basis
function   get_maxall_nL(self)
  implicit none
  integer:: get_maxall_nL
  type(basis_sturmian_nr), intent(in):: self
  integer:: l, l_min, l_max, numk, numktmp
  
  l_max = get_max_L(self)
  l_min = get_min_L(self)

  numk = 0
  do l=l_min,l_max
     numktmp =  get_max_nl(self,l)
     numk = max(numk,numktmp)
  enddo

  get_maxall_nL = numk
  
  return
end function get_maxall_nL
!
!
  function fpointer_nr(self)
    implicit none
    real(dp), pointer, dimension(:):: fpointer_nr
    type(sturmian_nr), intent(in):: self
    
    fpointer_nr => self%f
  end function fpointer_nr
!
  function gpointer_nr(self)
    implicit none
    real(dp), pointer, dimension(:):: gpointer_nr
    type(sturmian_nr), intent(in):: self
    
    gpointer_nr => self%g
  end function gpointer_nr
!
!
  function get_minf_nr(self)
    implicit none
    integer:: get_minf_nr
    type(sturmian_nr), intent(in):: self
    
    get_minf_nr = self%minf
  end function get_minf_nr
!
!
  function get_maxf_nr(self)
    implicit none
    integer:: get_maxf_nr
    type(sturmian_nr), intent(in):: self
    
    get_maxf_nr = self%maxf
  end function get_maxf_nr
!
!
  function get_alpha(self)
    implicit none
    real(dp):: get_alpha
    type(sturmian_nr), intent(in):: self
    
    get_alpha = self%alpha
  end function get_alpha
!
!
 function get_energy_nr(self)
    implicit none
    real(dp):: get_energy_nr
    type(sturmian_nr), intent(in):: self
    
    get_energy_nr = self%en
  end function get_energy_nr
!
  subroutine  print_energy_nr(self)
    use input_data
    implicit none
    type(basis_sturmian_nr), intent(in):: self
    integer:: n
    real(dp):: ioniz_en, exit_en 

     write(*,'("    N     J     exitation en.   ionization en.")')
    do n=1,self%N
       exit_en = (self%b(n)%en-self%b(1)%en)*data_in%eV  !  27.2116
       ioniz_en = (self%b(n)%en)*data_in%eV  !  27.2116  ! in eV
       write(*,'(i5,I6,2F17.5)') n, get_ang_mom(self%b(n)), exit_en, ioniz_en
    enddo
    print*

  end subroutine print_energy_nr
!
!
  subroutine  print_wf_nr(self,N,Rmax,filename)
    use grid_radial

    implicit none
    type(basis_sturmian_nr), intent(in):: self  ! basis
    integer, intent(in):: N                     ! number of functions to be printed
    real(dp), intent(in):: Rmax                   ! max value of R
    character(LEN=40), intent(in):: filename    ! name of file where to print functions
    integer:: i,m

    open(150,file=filename,status='REPLACE')
    
    write(150,'("# n:",12X,100(I10,5X))') (get_k(self%b(m)), m=1,N)
    write(150,'("# kappa:",8X,100(I10,5X))') (get_ang_mom(self%b(m)), m=1,N)
    write(150,'("# energy:",8X,100(10X,ES15.6,5X))') (get_energy(self%b(m)), m=1,N)
    do i=1,grid%nr
       if(grid%gridr(i) .le. grid%rmax .and. grid%gridr(i) .le. Rmax) then
          write(150,'(F15.5,100(1X,E14.5))') grid%gridr(i), (value(self%b(m),i), m=1,N)
       endif
    enddo
    close(150)

  end subroutine print_wf_nr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shellsort algorithm
  subroutine sort_by_energy_nr(self)
    type(basis_sturmian_nr), intent(inout):: self
    integer:: gap, i, j, N
    type(sturmian_nr):: Tmp

    N = self%n
    gap = N/2
    do 
       if(gap .le. 0) exit
       do i=gap+1,N
          call copy(Tmp,self%b(i))
          do j=i,gap+1,-gap
             if(Tmp%en .lt. self%b(j-gap)%en) then
                call copy(self%b(j),self%b(j-gap))
             else
                exit
             endif
             call copy(self%b(j-gap),Tmp)
          enddo
       enddo
       if ( gap .eq. 2) then
          gap = 1
       else
          gap = nint( real(gap)/2.2d0 )
       endif       
    enddo

  end subroutine sort_by_energy_nr
!
!******************************************************************
!
   subroutine ovlp3_nr(pn,pnp,v,i1,i2,result)
      use grid_radial 
      implicit none

      type(sturmian_nr), intent(in):: pn,pnp !  one-electron states
      real(dp), dimension(grid%nr), intent(in):: v
      integer, intent(in)::  i1, i2
      real(dp), intent(out):: result
!
      real(dp), dimension(grid%nr):: temp, fun
      real(dp), pointer, dimension(:):: f, fp
      integer:: minf,maxf, minfp,maxfp
      integer::  minfun, maxfun
!
      result = 0d0

      if(get_ang_mom(pn) .ne. get_ang_mom(pnp)) return
      if(get_ang_mom_proj(pn) .ne. get_ang_mom_proj(pnp)) return

      f  => fpointer(pn)
      fp => fpointer(pnp)
      minf = get_minf(pn)
      maxf = get_maxf(pn)
      minfp = get_minf(pnp)
      maxfp = get_maxf(pnp)


      result = 0.0d0
      minfun = max(max(minf,minfp),i1)
      maxfun = min(min(maxf,maxfp),i2)
         
      fun(minfun:maxfun) = f(minfun:maxfun)*fp(minfun:maxfun)*grid%weight(minfun:maxfun)*v(minfun:maxfun)
      result = SUM(fun(minfun:maxfun))         

    end subroutine ovlp3_nr
!
   subroutine ovlp_nr(pn,pnp,result)
      use grid_radial 
      implicit none

      type(sturmian_nr), intent(in):: pn,pnp !  one-electron states
      real(dp), intent(out):: result
!
      real(dp), dimension(grid%nr):: temp, fun
      real(dp), pointer, dimension(:):: f, fp
      integer:: minf,maxf, minfp,maxfp
      integer::  minfun, maxfun
!

      result = 0.0d0
      if(get_ang_mom(pn) .ne. get_ang_mom(pnp)) return

      f  => fpointer(pn)
      fp => fpointer(pnp)
      minf = get_minf(pn)
      maxf = get_maxf(pn)
      minfp = get_minf(pnp)
      maxfp = get_maxf(pnp)


      minfun = max(minf,minfp)
      maxfun = min(maxf,maxfp)
         
      fun(minfun:maxfun) = f(minfun:maxfun)*fp(minfun:maxfun)*grid%weight(minfun:maxfun)
      result = SUM(fun(minfun:maxfun))         

    end subroutine ovlp_nr
!
!
   subroutine me_1el_nr(lam,pn,pnp,v,i1,i2,result)
      use grid_radial 
      implicit none

      integer, intent(in):: lam
      type(sturmian_nr), intent(in):: pn,pnp !  one-electron states
      real(dp), dimension(grid%nr), intent(in):: v   ! v is a scalar
      integer, intent(in)::  i1, i2
      real(dp), intent(out):: result
!
      real(dp), dimension(grid%nr):: temp, fun
      real(dp), pointer, dimension(:):: f, fp
      integer:: minf,maxf, minfp,maxfp
      integer:: l, lp
      real(dp):: j, jp, rlam
      integer:: k
      real(dp):: CLAM
      real(dp):: tmp
      integer::  minfun, maxfun, i1old, i2old
!


      if(get_ang_mom(pn) .ne. get_ang_mom(pnp)) return
      if(get_ang_mom_proj(pn) .ne. get_ang_mom_proj(pnp)) return

      rlam = lam

      l  =  get_ang_mom(pn)
      lp =  get_ang_mom(pnp)
      f  => fpointer(pn)
      fp => fpointer(pnp)
      minf = get_minf(pn)
      maxf = get_maxf(pn)
      minfp = get_minf(pnp)
      maxfp = get_maxf(pnp)


      result = 0.0d0
      minfun = max(minf,minfp)
      maxfun = min(maxf,maxfp)
      tmp = sqrt(2*lp +1d0)*CLAM(dble(l),rlam,dble(lp))
         
      if(tmp .ne. 0d0) then

         fun(minfun:maxfun) = f(minfun:maxfun)*fp(minfun:maxfun)*grid%weight(minfun:maxfun)*v(minfun:maxfun)
         result = tmp * SUM(fun(minfun:maxfun))         
      endif

    end subroutine me_1el_nr

!
!!$----------------------------------------------------------------------------------------

subroutine gsort_nr(self)
  use  grid_radial

  implicit none

  type(basis_sturmian_nr), intent(inout):: self

  real(dp), dimension(grid%nr):: v
  real(dp), dimension(self%n,self%n):: ortint
  real(dp), pointer, dimension(:):: f1, f2
  integer:: nspm, n, m, i1, i2, i, nr
  integer:: l1, l2, m1, m2
  real(dp):: tmp, sum1
  real(dp)::  sum2
  integer:: j,k
  
  write(*,'(" G - S orthogonalization")') 
  
  nr = grid%nr
  nspm= self%n
  v(:) = 0.0d0

  do n=1,nspm
     f1 => fpointer(self%b(n))
     l1 = get_ang_mom(self%b(n))
     m1 = get_ang_mom_proj(self%b(n))
     i1 = get_minf(self%b(n))
     i2 = get_maxf(self%b(n))
     do m=1,n
        f2 => fpointer(self%b(m))         
        l2 = get_ang_mom(self%b(m))
        m2 = get_ang_mom_proj(self%b(m))
        i1 = max(i1,get_minf(self%b(m)))
        i2 = min(i2,get_maxf(self%b(m)))
        tmp = 0.0d0
        if(l1 .eq. l2 .and. m1 .eq. m2) then
           tmp = SUM( (f1(i1:i2)* f2(i1:i2) )*grid%weight(i1:i2) )
        endif
        ortint(n,m) = tmp
        ortint(m,n) = tmp
     enddo
  enddo

!!$     form overlap array <u_j|v_k>, u_j - old nonorthogonal set,
!!$     v_k - new set but not yet normalized
!!$     Only elements with  j >= k  required.
  do j=1,nspm
     do k=1,j
        sum1 = 0d0
        do n=1,k-1
           sum1 = sum1 + ortint(k,n)*ortint(j,n)/ortint(n,n)
        end do
        ortint(j,k) = ortint(j,k) - sum1
!!$            write(20,'("j,k =",2I3,", <j|k> =",F10.5)')  j, k, real(ortint(j,k))
     end do
  end do
      
      
!!$    form new orthonomal set vb_k, f_k = vb_k
  do k=1,nspm
     f1 => fpointer(self%b(k))
     v(:) = 0.0d0
     do i=1,nr
        sum2 = 0d0
        do n=1,k-1
           f2 => fpointer(self%b(n))
           i1 = get_minf(self%b(n))
           i2 = get_maxf(self%b(n))
           if(i .ge. i1 .and. i .le. i2) then
              sum2 = sum2 + dble(f2(i))*ortint(k,n)/dsqrt(ortint(n,n))
           endif
        end do
        v(i) = (value(self%b(k),i) - sum2)/dsqrt(ortint(k,k))
     end do
     call minmaxi2(v,nr,i1,i2)
     self%b(k)%maxf = i2
     self%b(k)%minf = i1
     deallocate(self%b(k)%f)
     allocate(self%b(k)%f(1:i2))

     tmp = sqrt(SUM( (v(i1:i2)*v(i1:i2))*grid%weight(i1:i2)))
!     print*,k,', tmp=', tmp, ortint(k,k)
     v(i1:i2) = v(i1:i2)/tmp

     self%b(k)%f(1:i2) = v(1:i2)
  end do

  return
end subroutine gsort_nr
!!$------------------------------------------------------------------------------------
! This subroutine make use nonrelativistic Sturmians
!!   V is a scalar
!!   the basis  self  is built from functions with the same angular momentum la  and its projection 
  subroutine construct_wflocalpot_nr(self,Nwf,V,nd,la,al)
    use grid_radial
    use MPI_module    

    implicit none
    
    type(basis_sturmian_nr), intent(out):: self
    integer, intent(in):: Nwf
    real(dp), dimension(grid%nr), intent(in):: V  !!!  V is a scalar
    integer, intent(in):: nd, la 
    real(dp), intent(in):: al


    type(basis_sturmian_nr):: bst
    real(dp), pointer, dimension(:):: weight, gridr
    real(dp), dimension(:,:), allocatable:: H, b, CI
    real(dp), dimension(:), allocatable:: w
    integer:: i, j, m, n, nr, nrj,  imax, imin, im1,im2, N_oneel, kk
    integer:: matz, ioerr
    real(dp), pointer, dimension(:):: fi, fj, fm
    integer:: minfi, maxfi, minfj, maxfj, minf, maxf
    real(dp):: tmpC
    real(dp), dimension(grid%nr):: temp
    real(dp):: energy, overlap, tmp, tmp1

    if(Nwf .gt. nd) then
       if (myid == 0) print*,'sturmian.f90: construct_wflocalpot_nr(): Nwf>nd:',Nwf , nd
    endif

    weight => grid%weight
    gridr => grid%gridr

    call new_basis(self,Nwf)

!!$ For given (la,nd,al) construct Sturmian basis.
    call construct(bst,nd, la, al)

!!$ Temporary arrays
    allocate(H(nd,nd))
    allocate(b(nd,nd))
    allocate(CI(nd,nd))
    allocate(w(nd))

!!$ Calculate H matrix
    b = 0d0
    H = 0d0
!!$ Overlap matrix
    b(1:nd,1:nd) = bst%ortint(1:nd,1:nd)

!!$ get kinetic energy: (using special properties of Lag func)
    H(1:nd,1:nd) = -al*al*b(1:nd,1:nd)/2d0
    do i = 1,nd
       H(i,i) = H(i,i) + al*al
    end do


!!$ Hamiltonian  matrix
    do i = 1,nd
       fi => fpointer(bst%b(i))
       minfi = get_minf(bst%b(i))
       maxfi = get_maxf(bst%b(i))
       do j = 1, i 
          fj => fpointer(bst%b(j))
          minfj = get_minf(bst%b(j))
          maxfj = get_maxf(bst%b(j))
          minf = max(minfi,minfj)
          maxf = min(maxfi,maxfj)
          tmp = SUM(weight(minf:maxf)*fi(minf:maxf)*fj(minf:maxf)*V(minf:maxf))
          H(i,j) = H(i,j)  + tmp
          H(j,i) = H(i,j) 
       end do
    end do

!!$      if (myid == 0) write(*,'(1P,5E15.5)') b
!!$      if (myid == 0) print*
!!$      if (myid == 0) print*
!!$      if (myid == 0) write(*,'(1P,5E15.5)') H


    matz=2
    call  rsg(nd,nd,H,b,w,matz,CI,ioerr)
    if(ioerr .ne. 0) then
       if (myid == 0) write(*,'("ioerr =",I3)') ioerr
       stop 'construct_wflocalpot: ERROR IN DIAGONALIZATION'
    endif

    if (myid == 0) write(*,'("Diagonalizing with L, nd,alpha:",2I5,F15.5)')la,nd,al
    do j=1, nd
       if (myid == 0)write(*,'(I5,1P,E15.5)') j,w(j)
    enddo
    if (myid == 0) print*

    N_oneel = 0
    do j=1,Nwf
       temp(:) = 0d0
       imin = grid%nr ; imax = 1
       do m=1,nd
          fm => fpointer(bst%b(m))
          tmpC =  CI(m,j)
          imax = get_maxf(bst%b(m))
          imin = get_minf(bst%b(m))
          do i=imin,imax
             temp(i) = temp(i) + tmpC*fm(i)
          enddo
       enddo
       call minmaxi(temp(1:grid%nr),grid%nr,im1,im2)
!!$ Fix sign of the one-electron functions
       if( sign(1d0,temp(im1+1)) .lt. 0)  then
          temp = -temp
       end if
       im1 = min(im1,imin)
       im2 = max(im2,imax)
       N_oneel = N_oneel + 1
       if(N_oneel .gt. Nwf) then
          if (myid == 0) print*,'N_oneel=',N_oneel, ', Nwf=', Nwf, ' nd=', nd
          stop 'N_oneel .gt. Nwf'
       endif
       ! record one-electron function
       energy = w(j)
       kk = j
       !              print*, '*** n,kk:', N_oneel, kk
       call init_function(self%b(N_oneel),energy,la,kk,im1,im2,temp,grid%nr) 
       
    end do
    
    deallocate(H)
    deallocate(b)
    deallocate(CI)
    deallocate(w)
    call destruct(bst)
 
    if(N_oneel .ne. Nwf) then
       if (myid == 0) print*,'N_oneel=',N_oneel, ', Nwf=', Nwf
       stop 'They should be equal!!!'
    endif
    
    call sort_by_energy(self)
!!$    call print_energy(self)
 
  end subroutine construct_wflocalpot_nr

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine expand_basis_m(bst)
  !
  implicit none

  type(basis_sturmian_nr), intent(inout) :: bst

  type(basis_sturmian_nr) :: bst_tmp
  type(sturmian_nr), pointer :: sturmOld,sturmNew
  integer :: nOld,nNon0,nNew, i,j,m


  nOld = basis_size(bst)   ! Number of functions in the original basis.
  nNon0 = nOld - get_num_m0(bst)   ! Number of these which have non-zero m.
  nNew = nOld + nNon0   ! Number of functions in the new, expanded basis.
  call new_basis(bst_tmp, nNew)

! THIS IS THE OLD METHOD FOR EXPANDING THE BASIS. IT RESULTS IN:
! m = 0 functions (l = 0, 1, 2, ...)
! m =+1 functions (l = 1, 2, ...)
! m =+2 functions (l = 2, ...)
! ...
! m =-1 functions (l = 1, 2, ...)
! m =-2 functions (l = 2, ...)
! ...
! THIS CAUSED DEGENERATE (+/-M) STATES TO BE REPRESENTED ASYMMETRICALLY.
!!$  ! First expand the underlying basis.
!!$  do i = 1, nOld
!!$     sturmOld => bst%b(i)
!!$     sturmNew => bst_tmp%b(i)
!!$     call copy(sturmNew, sturmOld)
!!$
!!$     m = get_ang_mom_proj(sturmNew)
!!$     if (m > 0) then
!!$        sturmNew => bst_tmp%b(i+nNon0)
!!$        call copy(sturmNew, sturmOld)
!!$        call set_ang_mom_proj(sturmNew, -m)
!!$     end if
!!$  end do

! THIS IS THE NEW METHOD FOR EXPANDING THE BASIS. IT RESULTS IN:
! m = 0 functions (l = 0, 1, 2, ...)
! m =+1 functions (l = 1, 2, ...)
! m =-1 functions (l = 1, 2, ...)
! m =+2 functions (l = 2, ...)
! m =-2 functions (l = 2, ...)
! ...
! WHILE NOT EXACTLY THE SAME AS SPHERICAL IT SHOULD REPRESENT SYMMETRICALLY.
  j = 0   ! Counter for the position in the new basis.
  do i = 1, nOld
     j = j+1   ! Basis function with original angular momentum projection = m.
     sturmOld => bst%b(i)
     sturmNew => bst_tmp%b(j)
     call copy(sturmNew, sturmOld)

     m = get_ang_mom_proj(sturmNew)
     if (m > 0) then   ! Don't need |m| because sturmians are created with m>0.
        j = j+1   ! New basis function with angular momentum projection = -m.
        sturmNew => bst_tmp%b(j)
        call copy(sturmNew, sturmOld)
        call set_ang_mom_proj(sturmNew, -m)
     end if
  end do
  if (j /= nNew) stop 'ERROR: Spheroidal basis was expanded incorrectly.'
  
  ! Then replace the original basis with the expanded basis.
  call destruct(bst)
  call new_basis(bst, nNew)
  do i = 1, nNew
     sturmOld => bst_tmp%b(i)
     sturmNew => bst%b(i)
     call copy(sturmNew, sturmOld)
  end do


end subroutine expand_basis_m

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine write_basis(bst, filename)
  implicit none

  ! Write a maximum of this many non-zero overlap & Hamiltonian matrix elements.
  integer, parameter :: countMax = 20   ! Can get messy for a big basis.

  type(basis_sturmian_nr), intent(in) :: bst
  character(len=*), intent(in) :: filename

  type(sturmian_nr), pointer :: sturm
  integer :: i,j,n, i1,i2, count
  real(dp) :: g1,g2, overlap,Ham_1e


  write(*,'(2A)') 'Writing sturmian basis to ', filename
  open(11, file=filename)

  n = basis_size(bst); count = 1
  do i = 1, n
     sturm => bst%b(i)
     i1 = get_minf(sturm); i2 = get_maxf(sturm)

     ! Write information about the sturmian to the file.
     if (count > 0) write(11,'(2A8,2A4,A8,2A6,4A16)') 'index', 'k','l','m', 'alpha', 'minf','maxf', 'first(f)','last(f)', 'first(g)?','last(g)?'
     write(11,'(2I8,2I4,F8.2,2I6,2F16.6)',advance='no') i, get_k(sturm),get_ang_mom(sturm),get_ang_mom_proj(sturm), get_alpha(sturm), i1,i2, value(sturm,i1),value(sturm,i2)
     if (associated(sturm%g)) then
        g1 = sturm%g(i1); g2 = sturm%g(i2)
        if ( g1/=0d0 .or. g2/=0d0 ) write(11,'(2F16.6)',advance='no') g1,g2
     end if
     write(11,*)

     ! Write up to parameter countMax non-zero matrix elements to file.
     j = i-1; count = 0
     do while ( count<countMax .and. j<n )
        j = j + 1
        sturm => bst%b(j)
        overlap = bst%ortint(i,j); Ham_1e = bst%ham1el(i,j)
        if ( overlap/=0d0 .or. Ham_1e/=0d0 ) then
           count = count + 1
           if (count == 1) write(11,'(A12,A8,2A4,2A12)') "index'", "k'","l'","m'", 'overlap','Ham_1e'
           write(11,'(A4,2I8,2I4,2F12.6)') '--->', j, get_k(sturm),get_ang_mom(sturm),get_ang_mom_proj(sturm), overlap, Ham_1e
        end if
     end do
     if (count > 0) write(11,*)

  end do

  close(11)

end subroutine write_basis


subroutine write_sturmians()
  use input_data
  implicit none

  type(sturmian_nr), pointer :: fn
  integer :: coord, lMax,mMax, i, l,m, matLen,matWid
  integer, dimension(:), allocatable :: mVec
  integer, dimension(:,:), allocatable :: lMat
  real(dp) :: R
  real(dp), dimension(:,:), allocatable :: coeffMat
  

  coord = data_in%calculation_type
  R = data_in%Rd
  lmax = data_in%latop
  mmax = data_in%Mt_max
  matLen = 0   ! Number of (l,m) combinations to write to file.
  matWid = 1   ! Size of expansion for spheroidal harmonics (N/A for sturmians).

  do l = 0, lmax
     do m = 0, min(l,mmax)
        matLen = matLen+1
     end do
  end do
  allocate( lMat(matLen,matWid), mVec(matLen), coeffMat(matLen,matWid) )

  i = 0
  do l = 0, lmax
     do m = 0, min(l,mmax)
        i = i+1
        lMat(i,1) = l
        mVec(i) = m
        coeffMat(i,1) = 1d0
     end do
  end do

 ! call write_harmonics( coord,R, lMax,mMax, matLen,matWid, lMat,mVec,coeffMat )

end subroutine write_sturmians


subroutine write_harmonics( coord,Rd, lMax,mMax, matLen,matWid, lmat,mVec,coeffMat )
  use Associated_Legendre_Functions
  use Special_Functions
  implicit none

  integer, intent(in) :: coord, lMax,mMax, matLen,matWid
  integer, dimension(matLen), intent(in) :: mVec
  integer, dimension(matLen,matWid), intent(in) :: lMat
  real(dp), intent(in) :: Rd
  real(dp), dimension(matLen,matWid), intent(in) :: coeffMat

  integer :: nazi,npol, iazi,ipol, i,j, fortNo, i1,i2
  real(dp) :: rad,pol,azi, coeff, rhoi,etai, re,im
  real(dp), dimension(:), allocatable :: thetaVec,phiVec, etaVec, fnVec


  nazi = 500
  npol = 500

  n_points = nazi
  l_max = lMax
  m_max = mMax
  directive = 'regular'
  Leg%D%A%Dir = recur
!  normalize = .false.
  write(*,3) l_max, m_max, n_points, Directive, Control, normalize, eps, recur, test_wron

  allocate( thetaVec(nazi),phiVec(npol), fnVec(nazi), Factor(0:l_max+m_max), x(n_points), Leg%R_LM%F(0:l_max,0:m_max),Leg%I_LM%F(0:l_max,0:m_max), Leg%R_LM%P(0:l_max,0:m_max,n_points),Leg%I_LM%Q(0:l_max,0:m_max,n_points) )

  if ( coord == 0 ) then
     thetaVec(:) = (/( pi*dble(iazi-1)/dble(nazi-1), iazi=1,nazi )/)
     x(:) = -cos( thetaVec(:) )
  elseif ( coord == 2 ) then
     allocate( etaVec(nazi) )
     etaVec(:) = (/( 1d0 - 2d0*dble(iazi-1)/dble(nazi-1), iazi=1,nazi )/)
     x(:) = etaVec(:)
  end if
  phiVec(:) = (/( 2d0*pi*dble(ipol-1)/dble(npol-1), ipol=1,npol )/)

  call Factorials
  call Legendre( R_LM=Leg%R_LM )
  Leg%R_LM%P = Leg%R_LM%P / dsqrt(2d0*pi)


  do i = 1, matLen
     fnVec(:) = 0d0
     m = mVec(i)
     if ( m < 0 ) cycle

     fortNo = 500+i
     do j = 1, matWid
        coeff = coeffMat(i,j)
        if ( abs(coeff) < 1E-9 ) cycle

        l = lMat(i,j)
!        fortNo = 900 + 10*l + m
        fnVec(:) = fnVec(:) + coeff * Leg%R_LM%P(l,m,:)
     end do

     open(fortNo)
     do iazi = 1, nazi

        rad = 0d0; azi = 0d0
        if ( coord == 0 ) then
           rad = fnVec(iazi)
           azi = thetaVec(iazi) - pi/2d0
        elseif ( coord==2 .or. coord==3 ) then
           rhoi = abs( fnVec(iazi) )
           etai = etaVec(iazi)
           rad = dsqrt( rhoi*(rhoi+Rd) + Rd*Rd/4d0*etai*etai )
           azi = (rhoi+Rd/2d0)*etai / rad
           
           if ( abs(azi) < 1d0 ) then
              azi = acos( (rhoi+Rd/2d0)*etai / rad ) - pi/2d0
           elseif ( abs(azi) < 1.01d0 ) then
              azi = -pi/2d0 * sign(1d0,azi)
           else
              stop 'sturmian.f90/write_harmonics() --- Spheroidal problems!'
           end if
        end if
        
        do ipol = 1, npol
           pol = phiVec(ipol)
           re = cos(dble(m)*pol)
           im = sin(dble(m)*pol)
           write(fortNo,'(5F20.10)') pol, azi, rad*rad, rad*re, rad*im
        end do
        write(fortNo,*)
     end do
     close(fortNo)
  end do

  deallocate( thetaVec,phiVec, Factor, x, Leg%R_LM%F,Leg%I_LM%F, Leg%R_LM%P,Leg%I_LM%Q )
  if (coord == 2) deallocate( etaVec )


stop!!
  return
3 Format(/,10x,'Maximum L                             = ',i6,1x, &
       'Maximum M                      = ',i6, &
       /,10x,'Number of Points                      = ',i6,1x, &
       'Type Functions Computed        = ',a10,1x, &
       /,10x,'Type Calculation               = ',a24,1x, &
       'Normalize Regular Functions on (-1,1) = ',l1,1x, &
       /,10x,'Continued Fraction Convergence = ',1pe15.8,1x, &
       'Backward Recurrence Method            = ',a24,1x, &
       /,10x, 'Test Wronskian                 = ',l1) 
end subroutine write_harmonics


end module sturmian_class

!----------------------------------------------------------------!
!
subroutine Ylam_nr(lam,pf1,pf2,maxfm,temp,i1,i2)
  use grid_radial 
  use sturmian_class
  implicit none
  
  integer, intent(in):: lam
  type(sturmian_nr), intent(in):: pf1,pf2 !  one-electron states
  integer, intent(in):: maxfm
  real(dp), dimension(grid%nr), intent(out):: temp
  integer, intent(out):: i1, i2
!
  real(dp), pointer, dimension(:):: f1, f2
  integer:: m1, m2
  real(dp), dimension(grid%nr):: fun

  m1 = max(get_minf(pf1),get_minf(pf2))
  m2 = min(get_maxf(pf1),get_maxf(pf2))

  f1 => fpointer(pf1)
  f2 => fpointer(pf2)
  
  fun(m1:m2) = f1(m1:m2)*f2(m1:m2) *grid%weight(m1:m2)

  call form(lam,fun,m1,m2,maxfm,temp,i1,i2)


end subroutine Ylam_nr





!===================================================================
!                                     m
!   Laguerre's polinomials  f(i,n) = L   (dlambda * grid(i))
!                                     n-1
!===================================================================
subroutine lagpol8 (m, dlambda, f, nload, nmax, grid, nr)
  implicit none
  integer:: nload, nmax, nr
  real(dp):: m, dlambda, f(nload, nmax), grid(nr) 
  real(dp):: pl(nmax)
  real(dp):: L2, r, x, pl1, pl2, pl3
  integer:: i, n
!
! INPUT:
! -----
!  m - parameter of Laguerre polinomial.
!  nmax  - the max number of n.
!  nload - dimension of f.
!  grid  - "R" - grid of "NR" points.
!  pl(nmax) - work space for this program.
!
! OUTPUT:
! ------
!  f(i,n) - Laguerre polinomials.
!
  L2 = m - 2
!
!     Loop by  r-grid
!
  do i = 1, nr
     r = grid(i)
     x = dlambda * r
!
!        Define Laguerre's polinomials
!        and store them in array 'f'
!
     pl1    = 1.0d0
     pl(1)  = pl1
     f(i,1) = pl1
!
     pl2    = dble(L2 + 3) - x
     if (nmax.gt.1) then
        pl(2)  = pl2
        f(i,2) = pl2
     endif
!
     do n = 3, nmax
        pl3 = ((dble(2*n-1+L2)-x)*pl2 - dble(n+L2)*pl1) / dble(n-1)
        pl(n) = pl3
        f(i,n) = pl3
        pl1 = pl2
        pl2 = pl3
     end do
  end do
  return
end subroutine lagpol8
!
!==================================================================
!
subroutine minmaxi(f,nr,i1,i2)
  use grid_radial
  implicit none
  integer, intent(in):: nr
  real(dp), intent(in), dimension(nr)::  f
  integer, intent(out):: i1, i2
  
  integer:: i

  i=1
  do while (i.lt.nr.and.abs(f(i)).lt.grid%regcut)
     i=i+1
  end do
  i1=i
  i=nr
  do while (i.gt.1.and.abs(f(i)).lt.grid%expcut)
     i=i-1
  end do
  i2=i
  if (i1 > i2) i1 = i2   ! If entire function is 0, set i1=i2=1.
  return
end subroutine minmaxi
!
subroutine minmaxi2(f,nr,i1,i2)
  use grid_radial
  implicit none
  integer, intent(in):: nr
  real(dp), intent(in), dimension(nr,2)::  f
  integer, intent(out):: i1, i2
  
  integer:: i, i1s,i1l,i2s,i2l

  i=1
  do while (i.lt.nr.and.abs(f(i,1)).lt.grid%regcut)
     i=i+1
  end do
  i1l=i
  i=nr
  do while (i.gt.1.and.abs(f(i,1)).lt.grid%expcut)
     i=i-1
  end do
  i2l=i

  i=1
  do while (i.lt.nr.and.abs(f(i,2)).lt.grid%regcut)
     i=i+1
  end do
  i1s=i
  i=nr
  do while (i.gt.1.and.abs(f(i,2)).lt.grid%expcut)
     i=i-1
  end do
  i2s=i

  i1 = min(i1s,i1l)
  i2 = max(i2s,i2l)

  return
end subroutine minmaxi2
!
!==================================================================
!
!!$function radial_derivative( inVec, i1,i2 )
!!$  !
!!$  use grid_radial
!!$  implicit none
!!$
!!$  integer, intent(in) :: i1,i2
!!$  real(dp), dimension(i1:i2), intent(in) :: inVec
!!$  real(dp), dimension(i1:i2) :: radial_derivative
!!$
!!$  real(dp), dimension(:), pointer :: rVec
!!$
!!$
!!$  radial_derivative(:) = 0d0
!!$  rVec => grid%gridr
!!$
!!$  radial_derivative(i1:i2-1) = (inVec(i1+1:i2) - inVec(i1:i2-1))/(rVec(i1+1:i2) - rVec(i1:i2-1))  
!!$
!!$end function radial_derivative
!
!==================================================================



