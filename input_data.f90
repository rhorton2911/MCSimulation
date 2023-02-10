module input_data

  use numbers
  public:: readin


  type input
         character(len=10):: target          ! target label
         real(dp):: energy             ! incident energy in au
         real(dp):: energyeV        ! incident energy in eV
         real(dp):: density         ! Density of the medium 
	       integer:: totalSims        ! total number of simulations
         integer:: debugOp          ! Debug option, runs tests on cross sections, printing results to files.
         integer:: posmode          !Positron mode, use positron as incident particle
				 integer:: trackSecEl            !Debugging option to switch of secondary electron tracking (0 = switch off, 1 = leave on)
				 integer:: numStatesIn
				 character(len=8), dimension(100):: statesToUse !List of H2 electronic states (including ionisation) to use.  
         integer:: totalVariations  !total number of time to vary input parameters for unceratainty propagation
         integer:: numToWrite     !Number of simulations to write to file for visualization.
         integer:: benchmarkOp    !run using benchmark cross sections (0 = don't run, 1 = use Garvey and Green CS) 
	       character(len=40):: ionop     ! ionisation option
         integer:: sdcsScaleOp    !OPtion of whether or not to scale SDCS
         integer:: stateIonop     !Choice to symmetrise +ve E pseudostate cross sections
         integer:: angleop        !scattering angle option. 1 = average, 2= distribution
         integer:: dicsop         !Dissociative Ionisation op. 1=CCC df, 2= DICS cross sections
         integer:: enlossop       !Inelastic scattering energy loss option (1=elastic, 2=inelastic)
         integer:: momOp           !momOp (0=formula,1=momentum transfer cs)
         integer::distortSdcsOp   !option to distort sdcs to investigate sensitivity to SDCS values
         integer:: vcsop          !Option of whether or not to use vibrational cross sections
	       integer:: groundVibOp    !Option of whether or not to restrict recording of vibrational excitations by primary to ground state only (1=yes, 0=no)
         integer:: dissOn         !Option for whether or not to allow dissociative excitation (1=yes, 0=no)
         !integer:: numVib
         integer:: radop          !Path length distribution option(1=exponential,0=constant mean)
         integer:: elScale        !elastic scattering CS scaling option for sensitivity analysis
         integer:: inelScale      !inelastic scattering CS scaling option for sensitivity analysis
         integer:: ionScale       !ionisation CS scaling option for sensitivity analysis
         integer:: Ntype              ! number of data types (positron, electron, etc...) 
         integer:: tcstype            ! type of totalcs files, 1 for Molecules,...
     character(len=50):: DATApath      ! path to the folder with the data on given type of particle
     character(len=40):: DATApathVcs   ! path to folder with vibrational cross sections for e-H2
     character(len=40):: DATApathdcsinel  !Path to inelastic dcs files 
     character(len=40):: DATApathDf
     character(len=50):: DATApathVcsPseudo  !Path to vibrational pseudostates
     integer:: Ntcs               ! number of totalcs files
     integer:: Naddics            ! number of additional ics files
     integer:: Ndf                ! number of diss. fraction files
     integer:: Nsdcs              ! number of SDCS files
     integer:: Ndcs               ! number of DCS files
     integer:: Ndcs1              ! number of dcs files of format 1
     integer:: Ndcs2              ! number of dcs file with second format
     integer:: Ndcsinel           ! number of files containing inelastic scattering dcs
     integer:: Ndcsinelcor        ! number of corrected inelastic dcs files
     integer:: NvcsEn             ! number of files explicitly containing vibrational state energies
     integer:: Nvcs               ! number of VCS files
     integer:: NvcsNew            ! number of vcs files in new data set
     integer:: NDissFrac          ! number of dissociation fractions for vibrationally resolved states
     integer:: NVcsPseudo
     !was length 40
     character(len=80), pointer, dimension(:) :: filename_tcs => NULL(), filename_addics => NULL(), filename_df => NULL(), filename_sdcs => NULL(), &
                                            filename_dcs => NULL(), filename_dcsinel => NULL(), filename_dcsinelcor => NULL(), filename_sdcsNew => NULL(), &
                                            filename_vcsEn => NULL(), filename_vcs => NULL(), filename_vcsNew => NULL(), filename_dfVib => NULL(),  &
                                            filename_VcsPseudo => NULL()
     character(len=80):: filename_X1SgStates
     real(dp):: eV = 27.21138505_dp    !SI to Hartees conversion factor
     real(dp):: bohrRadius = 5.29e-11_dp  !m^2
     real(dp):: eVToSi = 1.60e-19_dp
	 logical:: benchmark		  ! True if testing Ps formation Benchmark

  end type input
!****
  type(input):: data_in     ! contains all input data
!****

contains

  subroutine readin( self, nfile,iwrite )
      implicit none
      type(input):: self
      integer, intent(in):: nfile
      integer, intent(in):: iwrite  ! if =1 write to screen the read data, =0 no writing to screen
      real(dp):: eV, en_eV
      logical:: ex
      integer:: iostat_data, itcstype, iaddtcstype, idftype, isdcstype, ii
      character(len=20):: mode
      character(len=60)::tempDataPath

      eV = self%eV  !  27.2116
      
      if(nfile .eq. 10) then
         inquire(file='data.in',exist=ex)
         if(.not. ex) then
            print*,'File data.in does not exists. STOP'
            STOP
         end if
         open(nfile,file='data.in',iostat=iostat_data,action='read')
         if(iostat_data.ne.0) then
            print*, '********  File  canot open file data.in'
            stop
         end if
      endif
      
      read(nfile,*) mode
      !if(iwrite .eq. 1) write(*,*) 'PsFormationBenchmark: ', self%benchmark 
      if(mode .eq. 'PsBenchmarkTest') then
         self%benchmark = .true.
      else
         self%benchmark = .false.
      end if
      if(iwrite .eq. 1) write(*,*) 'mode: ', mode
    
      read(nfile,'(A10)') self%target
      if(iwrite .eq. 1) write(*,*) 'target: ', self%target

      en_eV = 0d0
      read(nfile,*) en_eV
    	self%energyeV = en_eV
      self%energy = en_eV / eV   ! convert to atomic units
      if(iwrite .eq. 1) write(*,*) 'energy: ', en_eV , ' eV'

     	read(nfile,*) self%totalSims
     	if(iwrite .eq. 1) write(*,*) 'totalSims: ', self%totalSims
   
      read(nfile,*) self%debugOp
      if(iwrite .eq. 1) write(*,*) 'debugOp: ', self%debugOp

      read(nfile,*) self%posmode
      if(iwrite .eq. 1) write(*,*) 'posmode: ', self%posmode

      read(nfile,*) self%trackSecEl
      if(iwrite .eq. 1) write(*,*) 'trackSecE: ', self%trackSecEl

      if (self%posmode .eq. 1) then
         self%trackSecEl = 0
      end if

      read(nfile,*) self%numStatesIn
      if(iwrite .eq. 1) write(*,*) 'numStatesIn: ', self%numStatesIn
				
			if (self%numStatesIn .ne. 0) then
         read(nfile,*) (self%statesToUse(ii), ii=1, self%numStatesIn) 
         if(iwrite .eq. 1) write(*,*) 'statesToUse: ',  (self%statesToUse(ii), ii=1, self%numStatesIn)
			else if (self%numStatesIn .eq. 0) then
				 !Default simulation, use all states.
				 self%statesToUse(:) = ""
				 read(nfile,* ) 
			end if

    	read(nfile,*) self%totalVariations
    	if(iwrite .eq. 1) write(*,*) 'totalVariations: ', self%totalVariations
    
    	read(nfile,*) self%numToWrite
    	if(iwrite .eq. 1) write(*,*) 'numToWrite: ', self%numToWrite
    
    	read(nfile,*) self%benchmarkOp
    	if(iwrite .eq. 1) write(*,*) 'benchmarkOp: ', self%benchmarkOp
  
      read(nfile,*) self%density	 
    	if(iwrite .eq. 1) write(*,*) 'density: ', self%density

      read(nfile,*) self%angleop
      if(iwrite .eq. 1) write(*,*) 'angleop: ', self%angleop

      read(nfile,*) self%dicsop
      if(iwrite .eq. 1) write(*,*) 'dicsop: ', self%dicsop

      read(nfile,*) self%enlossop
      if(iwrite .eq. 1) write(*,*) 'enlossop: ', self%enlossop

      read(nfile,*) self%momOp
      if(iwrite .eq. 1) write(*,*) 'momOp: ', self%momOp
	
    	read(nfile,*) self%ionop
    	if(iwrite .eq. 1) write(*,*) 'ionop: ', self%ionop

      read(nfile,*) self%sdcsScaleOp
      if(iwrite .eq. 1) write(*,*) 'sdcsScaleOp: ', self%sdcsScaleOp

      read(nfile,*) self%stateIonop 
    	if(iwrite .eq. 1) write(*,*) 'stateIonop: ', self%stateIonop

        read(nfile,*) self%distortSdcsOp
	if(iwrite .eq. 1) write(*,*) 'distortSdcsOp: ', self%distortSdcsOp

        read(nfile,*) self%vcsop
	if(iwrite .eq. 1) write(*,*) 'vcsop: ', self%vcsop

        read(nfile,*) self%groundVibOp
	if(iwrite .eq. 1) write(*,*) 'groundVibOp: ', self%groundVibOp

        read(nfile,*) self%dissOn
	if(iwrite .eq. 1) write(*,*) 'dissOn: ', self%dissOn
 
        !read(nfile,*) self%numVib
        !if(iwrite .eq. 1) write(*,*) 'numVib: ', self%numVib

        read(nfile,*) self%radop
        if(iwrite .eq. 1) write(*,*) 'radop: ', self%radop

        read(nfile,*) self%elScale
        if(iwrite .eq. 1) write(*,*) 'elScale: ', self%elScale

        read(nfile,*) self%inelScale
        if(iwrite .eq. 1) write(*,*) 'inelScale: ', self%inelScale

        read(nfile,*) self%ionScale
        if(iwrite .eq. 1) write(*,*) 'ionScale: ', self%ionScale

    read(nfile,*) self%Ntype     
    if(iwrite .eq. 1) write(*,*) 'Ntype: ', self%Ntype

    read(nfile,*) self%tcstype
    read(nfile,'(a)') self%DATApath
    if(iwrite .eq. 1) write(*,*) 'tcstype, DATApath: ', self%tcstype, self%DATApath

    read(nfile,*) self%DATApathdcsinel
    if(iwrite .eq. 1) write(*,*) 'DATApathdcsinel: ', self%DATApathdcsinel
 
    read(nfile,*) self%DATApathVcs
    if(iwrite .eq. 1) write(*,*) 'DATApathVcs: ', self%DATApathVcs

    read(nfile,*) self%DATApathDf
    if(iwrite .eq. 1) write(*,*) 'DATApathDf: ', self%DATApathDf

    read(nfile,*) self%DATApathVcsPseudo
    if(iwrite .eq. 1) write(*,*) 'DATApathVcsPseudo: ', self%DATApathVcsPseudo

    read(nfile,*) self%Ntcs, itcstype
    if(iwrite .eq. 1) write(*,*) 'Ntcs: ', self%Ntcs
    ! create array for totalcs file names    
    call read_files_block(self,nfile,iwrite,self%Ntcs,self%filename_tcs)

    read(nfile,*) self%Naddics, iaddtcstype
    if(iwrite .eq. 1) write(*,*) 'Naddics: ', self%Naddics
    ! create array for additional ics file names
    call read_files_block(self,nfile,iwrite,self%Naddics,self%filename_addics)

    read(nfile,*) self%Ndf, idftype
    if(iwrite .eq. 1) write(*,*) 'Ndf: ', self%Ndf
    ! create array for additional diss. fraction file names
    call read_files_block(self,nfile,iwrite,self%Ndf,self%filename_df)

    read(nfile,*) self%Nsdcs, isdcstype
    if(iwrite .eq. 1) write(*,*) 'Nsdcs: ', self%Nsdcs
    ! create array for additional diss. fraction file names
    call read_files_block(self,nfile,iwrite,self%Nsdcs,self%filename_sdcs)

    read(nfile,*) self%Ndcs, self%Ndcs1, self%Ndcs2
    if(iwrite .eq. 1) write(*,*) 'Ndcs, type1, type2: ', self%Ndcs, self%Ndcs1, self%Ndcs2
    ! create array for elastic dcs  file names
    call read_files_block(self,nfile,iwrite,self%Ndcs,self%filename_dcs)

    read(nfile,*) self%Ndcsinel
    if(iwrite .eq. 1) write(*,*) 'Ndcsinel: ', self%Ndcsinel
    ! create array for inelastic dcs  file names
    tempDataPath = TRIM(self%DATApath)
    self%DATApath = TRIM(self%DATApath)//'/'//TRIM(self%DATApathdcsinel)
    call read_files_block(self,nfile,iwrite,self%Ndcsinel,self%filename_dcsinel)
    self%DATApathdcsinel= TRIM(self%DATApath)
    self%DATApath = TRIM(tempDataPath)

    read(nfile,*) self%Ndcsinelcor
    if(iwrite .eq. 1) write(*,*) "Ndcsinelcor: ", self%Ndcsinelcor
    call read_files_block(self,nfile,iwrite,self%Ndcsinelcor,self%filename_dcsinelcor)
 
    read(nfile,*) self%NvcsEn
    if(iwrite .eq. 1) write(*,*) "NvcsEN: ", self%NvcsEN
    call read_files_block(self,nfile,iwrite,self%NvcsEn,self%filename_vcsEn)

    read(nfile,*) self%Nvcs
    if(iwrite .eq. 1) write(*,*) 'Nvcs: ', self%Nvcs
    !Create array for vibrational cs file names
    call read_files_block(self,nfile,iwrite,self%Nvcs,self%filename_vcs)

    read(nfile,*) self%NvcsNew
    if(iwrite .eq. 1) write(*,*) 'NvcsNew: ', self%NvcsNew
    !Vcs files in subdirectory of usual (DATA_1), due to large number
    tempDataPath = TRIM(self%DATApath)
    self%DATApath = TRIM(self%DATApath)//'/'//TRIM(self%DATApathVcs)
    call read_files_block(self,nfile,iwrite,self%NvcsNew,self%filename_vcsNew)
    self%DATApathVcs= TRIM(self%DATApath)
    self%DATApath = TRIM(tempDataPath)

    read(nfile,*) self%NDissFrac
    if(iwrite .eq. 1) write(*,*) 'NDissFrac: ', self%NDissFrac
    !Vcs files in subdirectory of usual (DATA_1), due to large number
    tempDataPath = self%DATApath
    self%DATApath = TRIM(self%DATApath)//'/'//TRIM(self%DATApathDf)
    call read_files_block(self,nfile,iwrite,self%NDissFrac,self%filename_dfVib)
    self%DATApathDf= TRIM(self%DATApath)
    self%DATApath = TRIM(tempDataPath)

    read(nfile,*) self%NVcsPseudo
    if(iwrite .eq. 1) write(*,*) 'NVcsPseudo: ', self%NVcsPseudo
    !Pseudostate files in subdirectory of usual (DATA_1), due to large number
    tempDataPath = self%DATApath
    self%DATApath = TRIM(self%DATApath)//'/'//TRIM(self%DATApathVcsPseudo)
    call read_files_block(self,nfile,iwrite,self%NVcsPseudo,self%filename_VcsPseudo)
    self%DATApathVcsPseudo= TRIM(self%DATApath)
    self%DATApath = TRIM(tempDataPath) 

    close(nfile)

    print*,'finish reading data_in'
    print*

  end subroutine readin
!
!
  subroutine read_files_block(self,nopenfile,iwrite,Nfiles,filename_ar)
    
    implicit none
    
    type(input):: self
    integer, intent(in):: nopenfile,iwrite,Nfiles
    character(len=80), pointer, dimension(:) :: filename_ar  !=> NULL()
    integer:: i

    if(associated(filename_ar)) then
       deallocate(filename_ar)
    endif
    allocate(filename_ar(Nfiles))
    do i=1,Nfiles
       read(nopenfile,'(A)') filename_ar(i)
       filename_ar(i) = TRIM(self%DATApath)//'/'//TRIM(filename_ar(i))
    enddo
    if(iwrite .eq. 1) then
       do i=1,Nfiles
          write(*,*) i, filename_ar(i)
       enddo
    endif

  end subroutine read_files_block


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: destruct_input
  !Purpose: deallocates any allocated memory in the imported 'input' type.
  !Date last modified: 06/02/2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine destruct_input(self)
      implicit none
      type(input)::self
    
      !Deallocates memory
      if (associated(self%filename_tcs)) then
         deallocate(self%filename_tcs)
      endif

      if (associated(self%filename_addics)) then
         deallocate(self%filename_addics)
      endif

      if (associated(self%filename_df)) then
         deallocate(self%filename_df)
      endif

      if (associated(self%filename_sdcs)) then
         deallocate(self%filename_sdcs)
      endif

      if (associated(self%filename_dcs)) then
         deallocate(self%filename_dcs)
      endif

      if (associated(self%filename_sdcsNew)) then
         deallocate(self%filename_dcs)
      endif

      if (associated(self%filename_vcsEn)) then
         deallocate(self%filename_vcsEn)
      endif

      if (associated(self%filename_vcs)) then
         deallocate(self%filename_vcs)
      endif

      if (associated(self%filename_vcsNew)) then
         deallocate(self%filename_vcsNew)
      endif

  end subroutine destruct_input
  
end module input_data
