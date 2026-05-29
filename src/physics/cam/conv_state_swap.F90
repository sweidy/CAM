module conv_state_swap
!=====================================================================
!
! Purpose: Implement conv_state_swap: subtract out bias in U,V,T,Q before
! CRM in SPCAM based on files 
!
! Author: Sarah Weidman
!
! Description:
!         
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size
  use phys_grid   ,   only:scatter_field_to_chunk, get_ncols_p, get_lat_all_p,get_lon_all_p, &
        get_rlat_all_p,get_lat_p,get_lon_p
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc
  use cam_logfile ,   only:iulog
  use physics_types,    only: physics_state, physics_tend, physics_ptend, physics_update, &
        physics_ptend_init
  use ppgrid,           only: begchunk, endchunk, pcols, pver, pverp, psubcols
  use phys_control, only: phys_getopts
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: conv_state_swap_readnl ! read namelist definitions
  public:: conv_state_swap_init ! initialize all the variables
  public:: update_conv_state_swap_profile ! determine whether to read new profile or not
  public:: conv_state_swap_in ! add in bias to state
  public:: conv_state_swap_out ! subtract bias to state
  public:: read_netcdf_conv_state_swap ! read convective state file
  public:: ConvStateSwap_Model ! for if statements outside of module
  ! 0 for tendency only from file, 1 for tendency from forcing only, 2 from difference from forcing to file
  public:: ConvStateSwap_use_forcing_opt
  public:: ConvStateSwap_tau

  ! Conv state swap parameters
  logical          :: ConvStateSwap_Model       =.false.
  integer          :: ConvStateSwap_use_forcing_opt
  character(len=cl):: ConvStateSwap_Path
  character(len=cs):: ConvStateSwap_File,ConvStateSwap_File_Template
  real(r8)         :: ConvStateSwap_tau
  integer          :: ConvStateSwap_Next_Year,ConvStateSwap_Next_Month
  integer          :: ConvStateSwap_Next_Day ,ConvStateSwap_Next_Sec
  integer          :: ConvStateSwap_Step
  real(r8)         :: ConvStateSwap_Ucoef
  real(r8)         :: ConvStateSwap_Vcoef
  real(r8)         :: ConvStateSwap_Qcoef
  real(r8)         :: ConvStateSwap_Tcoef
  integer          :: ConvStateSwap_Beg_Year ,ConvStateSwap_Beg_Month
  integer          :: ConvStateSwap_Beg_Day  ,ConvStateSwap_Beg_Sec


  logical :: ConvStateSwap_File_Present
  logical :: ConvStateSwap_Initialized =.false.
  integer ConvStateSwap_nlon,ConvStateSwap_nlat,ConvStateSwap_nlev

  ! ConvStateSwap observation arrays
  real(r8),allocatable::Ufield3d (:,:,:) !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Vfield3d (:,:,:) !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Tfield3d (:,:,:) !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Qfield3d (:,:,:) !(pcols,pver,begchunk:endchunk)

contains
  !================================================================
  subroutine conv_state_swap_readnl(nlfile)
   ! 
   ! conv_state_swap_READNL: Initialize default values controlling the conv_state_swap 
   !                 process. Then read namelist values to override 
   !                 them.
   !===============================================================
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn

   namelist /conv_state_swap_nl/ ConvStateSwap_Model,ConvStateSwap_Path,                       &
                         ConvStateSwap_File_Template, ConvStateSwap_tau,                       &
                         ConvStateSwap_Step,ConvStateSwap_use_forcing_opt,                     &
                         ConvStateSwap_Ucoef, ConvStateSwap_Vcoef, ConvStateSwap_Qcoef,        &
                         ConvStateSwap_Tcoef, ConvStateSwap_Beg_Year, ConvStateSwap_Beg_Month, &
                         ConvStateSwap_Beg_Day

   ! Set Default Namelist values
   !-----------------------------
   ConvStateSwap_Model         = .false.
   ConvStateSwap_use_forcing_opt = 0
   ConvStateSwap_Path          = '/n/holylfs06/LABS/kuang_lab/Lab/sweidman/IC_CESM2/'
   ConvStateSwap_File_Template = 'spcam_replay_uvtq.%m-%d-%s.nc'
   ConvStateSwap_tau           = 21600._r8 ! 30 minute forcing timescale
   ConvStateSwap_Step          = 21600._r8 ! read every 6 hrs
   ConvStateSwap_Ucoef         = 1._r8
   ConvStateSwap_Vcoef         = 1._r8
   ConvStateSwap_Qcoef         = 1._r8
   ConvStateSwap_Tcoef         = 1._r8
   ConvStateSwap_Beg_Year      = 1980
   ConvStateSwap_Beg_Month     = 1
   ConvStateSwap_Beg_Day       = 1
   ConvStateSwap_Beg_Sec       = 0

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'conv_state_swap_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,conv_state_swap_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('conv_state_swap_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(ConvStateSwap_Path         ,len(ConvStateSwap_Path)         ,mpichar,0,mpicom)
   call mpibcast(ConvStateSwap_File_Template,len(ConvStateSwap_File_Template),mpichar,0,mpicom)
   call mpibcast(ConvStateSwap_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(ConvStateSwap_use_forcing_opt        , 1, mpiint, 0, mpicom)
   call mpibcast(ConvStateSwap_tau          , 1, mpir8 , 0, mpicom)
   call mpibcast(ConvStateSwap_Step          , 1, mpir8 , 0, mpicom)
   call mpibcast(ConvStateSwap_Ucoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(ConvStateSwap_Vcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(ConvStateSwap_Tcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(ConvStateSwap_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(ConvStateSwap_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(ConvStateSwap_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(ConvStateSwap_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(ConvStateSwap_Beg_Sec      , 1, mpiint, 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! conv_state_swap_readnl

  !================================================================

  !================================================================
  subroutine conv_state_swap_init
    ! 
    ! conv_state_swap_INIT: Allocate space and initialize corrector values
    !===============================================================
    use ppgrid        ,only: pver,pcols,begchunk,endchunk
    use error_messages,only: alloc_err
    use dycore        ,only: dycore_is
    use dyn_grid      ,only: get_horiz_grid_dim_d
    use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
    use cam_history   ,only: addfld,add_default
    use shr_const_mod ,only: SHR_CONST_PI
    use filenames     ,only: interpret_filename_spec
 
    ! Local values
    !----------------
    integer  Year,Month,Day,Sec, Beg_Sec
    integer  YMD1,YMD
    logical  After_Beg
    integer  istat,lchnk,ncol,icol,ilev
    integer  hdim1_d,hdim2_d
    integer  dtime
    real(r8) rlat,rlon

    ! Get the time step size
    !------------------------
    dtime = get_step_size()
    Beg_sec = 0 ! run must start at midnight
 
    ! Allocate Space for corrector data arrays
    !-----------------------------------------
    allocate(Ufield3d(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Ufield3d',pcols*pver*((endchunk-begchunk)+1))
    allocate(Vfield3d(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Vfield3d',pcols*pver*((endchunk-begchunk)+1))
    allocate(Tfield3d(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Tfield3d',pcols*pver*((endchunk-begchunk)+1))
    allocate(Qfield3d(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'corrector_init','Qfield3d',pcols*pver*((endchunk-begchunk)+1))


    ! Values initialized only by masterproc
    !-----------------------------------------
    if(masterproc) then

      call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
      ConvStateSwap_nlon=hdim1_d
      ConvStateSwap_nlat=hdim2_d
      ConvStateSwap_nlev=pver

      ! Check the time relative to the corrector window
      !------------------------------------------------
      call get_curr_date(Year,Month,Day,Sec)
      YMD=(Year*10000) + (Month*100) + Day
      YMD1=(ConvStateSwap_Beg_Year*10000) + (ConvStateSwap_Beg_Month*100) + ConvStateSwap_Beg_Day
      call timemgr_time_ge(YMD1,ConvStateSwap_Beg_Sec,         &
                            YMD ,Sec          ,After_Beg)

      ! Set Time indicies so that the next call to 
      ! timestep_init will initialize the data arrays.
      !--------------------------------------------
      if (After_Beg) then
        ConvStateSwap_Next_Year =Year
        ConvStateSwap_Next_Month=Month
        ConvStateSwap_Next_Day  =Day
        ConvStateSwap_Next_Sec  =(Sec/ConvStateSwap_Step)*ConvStateSwap_Step
      else
        ConvStateSwap_Next_Year =ConvStateSwap_Beg_Year
        ConvStateSwap_Next_Month=ConvStateSwap_Beg_Month
        ConvStateSwap_Next_Day  =ConvStateSwap_Beg_Day
        ConvStateSwap_Next_Sec  =ConvStateSwap_Beg_Sec
      endif

      ConvStateSwap_File_Present=.false.

      ! Initialization is done, 
      !--------------------------
      ConvStateSwap_Initialized=.true.

      ! Check that this is a valid DYCORE model
      !------------------------------------------
      if(.not.dycore_is('LR')) then
        call endrun('corrector IS CURRENTLY ONLY CONFIGURED FOR FV')
      endif

      ! Informational Output
      !---------------------------
      write(iulog,*) ' '
      write(iulog,*) '---------------------------------------------------------'
      write(iulog,*) '  MODEL ConvStateSwap INITIALIZED WITH THE FOLLOWING SETTINGS: '
      write(iulog,*) '---------------------------------------------------------'
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Model=',ConvStateSwap_Model
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_use_forcing_opt=',ConvStateSwap_use_forcing_opt
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Path=',ConvStateSwap_Path
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_File_Template =',ConvStateSwap_File_Template
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_tau  =',ConvStateSwap_tau
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Step  =',ConvStateSwap_Step
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Ucoef  =',ConvStateSwap_Ucoef
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Vcoef  =',ConvStateSwap_Vcoef
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Qcoef  =',ConvStateSwap_Qcoef
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Tcoef  =',ConvStateSwap_Tcoef
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Beg_Year  =',ConvStateSwap_Beg_Year
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Beg_Month  =',ConvStateSwap_Beg_Month
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Beg_Day  =',ConvStateSwap_Beg_Day
      write(iulog,*) 'ConvStateSwap: ConvStateSwap_Beg_Sec  =',ConvStateSwap_Beg_Sec
end if ! masterproc

#ifdef SPMD
    call mpibcast(ConvStateSwap_Step          ,            1, mpir8 , 0, mpicom)
    call mpibcast(ConvStateSwap_Next_Year     ,            1, mpiint, 0, mpicom)
    call mpibcast(ConvStateSwap_Next_Month    ,            1, mpiint, 0, mpicom)
    call mpibcast(ConvStateSwap_Next_Day      ,            1, mpiint, 0, mpicom)
    call mpibcast(ConvStateSwap_Next_Sec      ,            1, mpiint, 0, mpicom)
    call mpibcast(ConvStateSwap_Model         ,            1, mpilog, 0, mpicom)
    call mpibcast(ConvStateSwap_use_forcing_opt         ,  1, mpiint, 0, mpicom)
    call mpibcast(ConvStateSwap_Initialized   ,            1, mpilog, 0, mpicom)
    call mpibcast(ConvStateSwap_nlev          ,            1, mpiint, 0, mpicom)
    call mpibcast(ConvStateSwap_nlon          ,            1, mpiint, 0, mpicom)
    call mpibcast(ConvStateSwap_nlat          ,            1, mpiint, 0, mpicom)
#endif

    ! Initialize the analysis filename at the NEXT time for startup.
      !---------------------------------------------------------------

    if ((ConvStateSwap_use_forcing_opt.eq.0).or.(ConvStateSwap_use_forcing_opt.eq.2)) then
      ConvStateSwap_File=interpret_filename_spec(ConvStateSwap_File_Template      , &
                                        yr_spec=Year , &
                                        mon_spec=Month, &
                                        day_spec=Day  , &
                                        sec_spec=Beg_Sec    )

      if(masterproc) then
      write(iulog,*) 'ConvStateSwap: Reading forcing:',trim(ConvStateSwap_Path)//trim(ConvStateSwap_File)
      endif

      call read_netcdf_conv_state_swap (trim(ConvStateSwap_Path)//trim(ConvStateSwap_File))
    endif

    ! Load zeros into arrays
    !------------------------------------------------------
    do lchnk=begchunk,endchunk
      Ufield3d(:pcols,:pver,lchnk)=0._r8
      Vfield3d(:pcols,:pver,lchnk)=0._r8
      Qfield3d(:pcols,:pver,lchnk)=0._r8
      Tfield3d(:pcols,:pver,lchnk)=0._r8
    end do

    ! End Routine
    !------------

    return
  end subroutine



  subroutine read_netcdf_conv_state_swap(analysis_file)
    use netcdf
    use ppgrid,    only: pver,pcols,begchunk,endchunk
    use dyn_grid      ,only: get_horiz_grid_dim_d
    use error_messages,only: alloc_err
    use cam_abortutils, only:endrun
    !implicit none
 
    ! Arguments
    !-------------
    character(len=*),intent(in):: analysis_file
 
    ! Local values
    !-------------
    integer lev
    integer nlon,nlat,plev,istat
    integer ncid,varid
    integer ilat,ilon,ilev
    real(r8) Xanal(ConvStateSwap_nlon,ConvStateSwap_nlat,pver) 
    real(r8) Lat_anal(ConvStateSwap_nlat)
    real(r8) Lon_anal(ConvStateSwap_nlon)
    real(r8) Xtrans(ConvStateSwap_nlon,pver,ConvStateSwap_nlat)
    !integer  hdim1_d,hdim2_d
  
    if(masterproc) then
    ! Open the given file
      !-----------------------
    istat=nf90_open(trim(analysis_file),NF90_NOWRITE,ncid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*)'NF90_OPEN: failed for file ',trim(analysis_file)
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    ! Read in Dimensions
      !--------------------
    istat=nf90_inq_dimid(ncid,'lon',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_inquire_dimension(ncid,varid,len=nlon)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    istat=nf90_inq_dimid(ncid,'lat',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_inquire_dimension(ncid,varid,len=nlat)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    istat=nf90_inq_dimid(ncid,'lev',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_inquire_dimension(ncid,varid,len=plev)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    istat=nf90_inq_varid(ncid,'lon',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Lon_anal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    istat=nf90_inq_varid(ncid,'lat',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Lat_anal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
 
    if((ConvStateSwap_nlon.ne.nlon).or.(ConvStateSwap_nlat.ne.nlat).or.(plev.ne.pver)) then
     write(iulog,*) 'ERROR: ConvStateSwap_update_analyses_fv: nlon=',nlon,' ConvStateSwap_nlon=',ConvStateSwap_nlon
     write(iulog,*) 'ERROR: ConvStateSwap_update_analyses_fv: nlat=',nlat,' ConvStateSwap_nlat=',ConvStateSwap_nlat
     write(iulog,*) 'ERROR: ConvStateSwap_update_analyses_fv: plev=',plev,' pver=',pver
     call endrun('ConvStateSwap_update_analyses_fv: analyses dimension mismatch')
    endif
 
    ! Read in, transpose lat/lev indices, 
      ! and scatter data arrays
      !----------------------------------
    istat=nf90_inq_varid(ncid,'U',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Xanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    do ilat=1,nlat
    do ilev=1,plev
    do ilon=1,nlon
      Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    end do
    end do
    end do
  endif ! (masterproc) then
  call scatter_field_to_chunk(1,ConvStateSwap_nlev,1,ConvStateSwap_nlon,Xtrans,   &
                              Ufield3d(1,1,begchunk))
 
  if(masterproc) then
    istat=nf90_inq_varid(ncid,'V',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Xanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    do ilat=1,nlat
    do ilev=1,plev
    do ilon=1,nlon
      Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    end do
    end do
    end do
  endif ! (masterproc) then
  call scatter_field_to_chunk(1,ConvStateSwap_nlev,1,ConvStateSwap_nlon,Xtrans,   &
                              Vfield3d(1,1,begchunk))
 
  if(masterproc) then
    istat=nf90_inq_varid(ncid,'T',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Xanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    do ilat=1,nlat
    do ilev=1,plev
    do ilon=1,nlon
      Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    end do
    end do
    end do
  endif ! (masterproc) then
  call scatter_field_to_chunk(1,ConvStateSwap_nlev,1,ConvStateSwap_nlon,Xtrans,   &
                              Tfield3d(1,1,begchunk))
 
  if(masterproc) then
    istat=nf90_inq_varid(ncid,'Q',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    istat=nf90_get_var(ncid,varid,Xanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    do ilat=1,nlat
    do ilev=1,plev
    do ilon=1,nlon
      Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
    end do
    end do
    end do
 
  ! Close analysis file
  istat=nf90_close(ncid)
  if(istat.ne.NF90_NOERR) then
    write(iulog,*) nf90_strerror(istat)
    call endrun ('UPDATE_ANALYSES_EUL')
  endif
 
  endif ! (masterproc) then
  call scatter_field_to_chunk(1,ConvStateSwap_nlev,1,ConvStateSwap_nlon,Xtrans,   &
                              Qfield3d(1,1,begchunk))
 
    ! End Routine
    !------------
    return
 
   end subroutine read_netcdf_conv_state_swap


   subroutine update_conv_state_swap_profile (state)

  !----------------------------------------------------------------------- 
  ! Purpose: 
  !-----------------------------------------------------------------------
   use physics_buffer, only : pbuf_get_index, dtype_r8
   use physconst,    only: cpair
   use dyn_grid,     only: get_horiz_grid_dim_d
   use time_manager, only: get_nstep, get_curr_date
   use filenames     ,only: interpret_filename_spec
   use constituents,       only: cnst_get_ind
   use corrector,       only: Force_Model,Force_ON,corrector_timestep_tend
   use running_mean,    only: Running_mean_Model,Running_mean_ON,running_mean_timestep_tend,Running_nudge_U,Running_nudge_V,Running_nudge_T,Running_nudge_Q
   use physics_types,   only: physics_ptend_reset
   
   integer, save :: nstep_count

  ! Arguments
   type(physics_state), intent(inout) :: state(begchunk:endchunk)

  ! Local workspace
   type(physics_ptend)   :: ptend(begchunk:endchunk)                  ! indivdual parameterization tendencies

   integer :: i, j,n ,ilat                 ! longitude, latitude,field, and global column indices
   integer :: c, ncols, k, istep
   real(r8) ::rlat(pcols),conv_forcingtime!,tmprand(128,64)
   integer :: ilat_all(pcols)
   integer :: ndays, Day, Month, Year, ncsec, dtime, Sec
   integer YMD1,YMD2,YMD
   logical After_Beg
   integer :: modstep6hr
   logical :: Update_ConvState
   logical :: fileexists
   logical,save :: conv_started = .false. 
   integer :: ierr,csize     
   integer :: indw                      
  !-----
      
   real(r8), allocatable :: tmpfield_conv(:)

   istep=get_nstep()
   nstep_count=istep
   call cnst_get_ind('Q',indw)

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(ConvStateSwap_Beg_Year*10000) + (ConvStateSwap_Beg_Month*100) + ConvStateSwap_Beg_Day
   call timemgr_time_ge(YMD1,ConvStateSwap_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   ! allocate state forcing
   if (conv_started .ne. .TRUE.) then
      conv_started=.TRUE.
      
#if ( defined SPMD )
          do c = begchunk, endchunk
              ncols = get_ncols_p(c)
              do i = 1, ncols
                  do k=1,pver
                      state(c)%qconvforce(i,k) = 0.0
                      state(c)%uconvforce(i,k) = 0.0
                      state(c)%vconvforce(i,k) = 0.0
                      state(c)%sconvforce(i,k) = 0.0
                  end do
              end do
          end do
#endif
    endif

    if (After_Beg) then

    if (ConvStateSwap_use_forcing_opt.eq.0) then

      ! determine whether time to update convection state file
      call get_curr_date(Year,Month,Day,ncsec)
      YMD=(Year*10000) + (Month*100) + Day
      YMD1=(ConvStateSwap_Next_Year*10000) + (ConvStateSwap_Next_Month*100) + ConvStateSwap_Next_Day
      call timemgr_time_ge(YMD1,ConvStateSwap_Next_Sec,            &
                          YMD ,ncsec           ,Update_ConvState)

      ! if time to read
      if (Update_ConvState) then

        if (masterproc) print*, 'time to update convection bias file'
        if (masterproc) write(iulog,*) "ncsec, next sec", ncsec,ConvStateSwap_Next_Sec 

        fileexists=.FALSE.
      
        ! check if file exists and if so, read in
        ConvStateSwap_File=interpret_filename_spec(ConvStateSwap_File_Template      , &
              yr_spec=ConvStateSwap_Next_Year , &
              mon_spec=ConvStateSwap_Next_Month, &
              day_spec=ConvStateSwap_Next_Day  , &
              sec_spec=ConvStateSwap_Next_Sec    )
        
        INQUIRE(FILE=trim(ConvStateSwap_Path)//trim(ConvStateSwap_File), EXIST=fileexists)
        
        if (.not. fileexists) print*, 'file missing', trim(ConvStateSwap_Path)//trim(ConvStateSwap_File)

        if (fileexists) then
          if(masterproc) then
            write(iulog,*) 'ConvStateSwap: Reading convective state file:',trim(ConvStateSwap_Path)//trim(ConvStateSwap_File)
          endif

          call read_netcdf_conv_state_swap (trim(ConvStateSwap_Path)//trim(ConvStateSwap_File))

        end if

        ! increment time for next read
        call timemgr_time_inc(YMD1,ConvStateSwap_Next_Sec,              &
        YMD2,ConvStateSwap_Next_Sec,ConvStateSwap_Step,0,0)
        ConvStateSwap_Next_Year =(YMD2/10000)
        YMD2            = YMD2-(ConvStateSwap_Next_Year*10000)
        ConvStateSwap_Next_Month=(YMD2/100)
        ConvStateSwap_Next_Day  = YMD2-(ConvStateSwap_Next_Month*100)
      
        do c = begchunk, endchunk
          ncols = get_ncols_p(c)
          do i = 1, ncols
              do k=1,pver
                  state(c)%qconvforce(i,k) = (Qfield3d(i,k,c))/ConvStateSwap_tau*ConvStateSwap_Qcoef
                  state(c)%uconvforce(i,k) = (Ufield3d(i,k,c))/ConvStateSwap_tau*ConvStateSwap_Ucoef
                  state(c)%vconvforce(i,k) = (Vfield3d(i,k,c))/ConvStateSwap_tau*ConvStateSwap_Vcoef
                  state(c)%sconvforce(i,k) = (Tfield3d(i,k,c))/ConvStateSwap_tau*ConvStateSwap_Tcoef
              end do
          end do
        end do
      endif

    elseif(ConvStateSwap_use_forcing_opt.eq.1) then

      if((Running_mean_Model).and.(Running_mean_ON)) then
        
        ! subtract state added from nudging
        ! should be negative of whatever is being applied from nudging. conv_swap_in will add this as a positive tendency
#if ( defined SPMD )
        do c = begchunk, endchunk
          call running_mean_timestep_tend(state(c),ptend(c))
          ncols = get_ncols_p(c)

          state(c)%uconvforce(:ncols,:pver) = -ptend(c)%u(:ncols,:pver)*21600._r8/ConvStateSwap_tau*ConvStateSwap_Ucoef
          state(c)%vconvforce(:ncols,:pver) = -ptend(c)%v(:ncols,:pver)*21600._r8/ConvStateSwap_tau*ConvStateSwap_Vcoef
          state(c)%sconvforce(:ncols,:pver) = -ptend(c)%s(:ncols,:pver)*21600._r8/ConvStateSwap_tau/cpair*ConvStateSwap_Tcoef
          state(c)%qconvforce(:ncols,:pver) = -ptend(c)%q(:ncols,:pver,indw)*21600._r8/ConvStateSwap_tau*ConvStateSwap_Qcoef

          call physics_ptend_reset(ptend(c))
        enddo
#endif

      else if((Force_Model).and.(Force_ON)) then
        
        ! subtract state added from corrector
#if ( defined SPMD )
        do c = begchunk, endchunk
          call corrector_timestep_tend(state(c),ptend(c))
          ncols = get_ncols_p(c)

          state(c)%uconvforce(:ncols,:pver) = -ptend(c)%u(:ncols,:pver)*21600._r8/ConvStateSwap_tau*ConvStateSwap_Ucoef
          state(c)%vconvforce(:ncols,:pver) = -ptend(c)%v(:ncols,:pver)*21600._r8/ConvStateSwap_tau*ConvStateSwap_Vcoef
          state(c)%sconvforce(:ncols,:pver) = -ptend(c)%s(:ncols,:pver)*21600._r8/ConvStateSwap_tau/cpair*ConvStateSwap_Tcoef
          state(c)%qconvforce(:ncols,:pver) = -ptend(c)%q(:ncols,:pver,indw)*21600._r8/ConvStateSwap_tau*ConvStateSwap_Qcoef

          call physics_ptend_reset(ptend(c))
        enddo
#endif

      else
        if(masterproc) then 
          write(iulog,*) "WARNING: running mean and correction not applied, need different opt "
        endif 
      endif 

    elseif(ConvStateSwap_use_forcing_opt.eq.2) then ! use file and tendency

     ! determine whether time to update convection state file
      call get_curr_date(Year,Month,Day,ncsec)
      YMD=(Year*10000) + (Month*100) + Day
      YMD1=(ConvStateSwap_Next_Year*10000) + (ConvStateSwap_Next_Month*100) + ConvStateSwap_Next_Day
      call timemgr_time_ge(YMD1,ConvStateSwap_Next_Sec,            &
                          YMD ,ncsec           ,Update_ConvState)

      ! if time to read
      if (Update_ConvState) then

        if (masterproc) write(iulog,*) 'time to update convection bias file'
        if (masterproc) write(iulog,*) "ncsec, next sec", ncsec,ConvStateSwap_Next_Sec 

        fileexists=.FALSE.
      
        ! check if file exists and if so, read in
        ConvStateSwap_File=interpret_filename_spec(ConvStateSwap_File_Template      , &
              yr_spec=ConvStateSwap_Next_Year , &
              mon_spec=ConvStateSwap_Next_Month, &
              day_spec=ConvStateSwap_Next_Day  , &
              sec_spec=ConvStateSwap_Next_Sec    )
        
        INQUIRE(FILE=trim(ConvStateSwap_Path)//trim(ConvStateSwap_File), EXIST=fileexists)
        
        if (.not. fileexists) print*, 'file missing', trim(ConvStateSwap_Path)//trim(ConvStateSwap_File)

        if (fileexists) then
          if(masterproc) then
            write(iulog,*) 'ConvStateSwap: Reading convective state file:',trim(ConvStateSwap_Path)//trim(ConvStateSwap_File)
          endif

          call read_netcdf_conv_state_swap (trim(ConvStateSwap_Path)//trim(ConvStateSwap_File))

          if(masterproc) then 
            write(iulog,*) "reading file, Qfield3d(i,k,c) ", Qfield3d(1,20,begchunk)
          endif
        end if

        ! increment time for next read
        call timemgr_time_inc(YMD1,ConvStateSwap_Next_Sec,              &
        YMD2,ConvStateSwap_Next_Sec,ConvStateSwap_Step,0,0)
        ConvStateSwap_Next_Year =(YMD2/10000)
        YMD2            = YMD2-(ConvStateSwap_Next_Year*10000)
        ConvStateSwap_Next_Month=(YMD2/100)
        ConvStateSwap_Next_Day  = YMD2-(ConvStateSwap_Next_Month*100)
      endif  ! update file 

        if((Running_mean_Model).and.(Running_mean_ON)) then
#if ( defined SPMD )        
        do c = begchunk, endchunk
          ncols = get_ncols_p(c)

          do i = 1, ncols
              do k=1,pver

              state(c)%uconvforce(i,k) = (Ufield3d(i,k,c)-Running_nudge_U(i,k,c))/ConvStateSwap_tau*ConvStateSwap_Ucoef
              state(c)%vconvforce(i,k) = (Vfield3d(i,k,c)-Running_nudge_V(i,k,c))/ConvStateSwap_tau*ConvStateSwap_Vcoef
              state(c)%sconvforce(i,k) = (Tfield3d(i,k,c)-Running_nudge_T(i,k,c))/ConvStateSwap_tau*ConvStateSwap_Tcoef
              state(c)%qconvforce(i,k) = (Qfield3d(i,k,c)-Running_nudge_Q(i,k,c))/ConvStateSwap_tau*ConvStateSwap_Qcoef

              enddo
          enddo
        enddo
#endif
      else if((Force_Model).and.(Force_ON)) then

        ! subtract state added from corrector
        do c = begchunk, endchunk
          call corrector_timestep_tend(state(c),ptend(c))

          if(masterproc) then 
            write(iulog,*) "in Force Model c",c
            write(iulog,*) "Ufield3d(i,k,c) ",Ufield3d(1,20,c)
            write(iulog,*) "ptend(c)%u(:ncols,:pver) ",ptend(c)%u(1,20)
            write(iulog,*) "state(c)%uconvforce(:ncols,:pver) ",state(c)%uconvforce(1,20)
          endif 
          ncols = get_ncols_p(c)

          state(c)%uconvforce(:ncols,:pver) = (Ufield3d(:ncols,:pver,c)-ptend(c)%u(:ncols,:pver)*21600._r8)/ConvStateSwap_tau*ConvStateSwap_Ucoef
          state(c)%vconvforce(:ncols,:pver) = (Vfield3d(:ncols,:pver,c)-ptend(c)%v(:ncols,:pver)*21600._r8)/ConvStateSwap_tau*ConvStateSwap_Vcoef
          state(c)%sconvforce(:ncols,:pver) = (Tfield3d(:ncols,:pver,c)-ptend(c)%s(:ncols,:pver)*21600._r8/cpair)/ConvStateSwap_tau*ConvStateSwap_Tcoef
          state(c)%qconvforce(:ncols,:pver) = (Qfield3d(:ncols,:pver,c)-ptend(c)%q(:ncols,:pver,indw)*21600._r8)/ConvStateSwap_tau*ConvStateSwap_Qcoef

          call physics_ptend_reset(ptend(c))
        enddo

      else
        if(masterproc) then 
          write(iulog,*) "WARNING: running mean not applied, need different opt "
        endif 
      endif  ! if running_mean (updates with climo)

      
    endif ! tend_opt    
    endif ! After_Beg   
  end subroutine update_conv_state_swap_profile


   subroutine conv_state_swap_in (ztodt, state,tend)
    !----------------------------------------------------------------------- 
      ! Purpose: 
      !-----------------------------------------------------------------------
    use physconst,    only: cpair
    use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
    use time_manager, only: get_nstep
    use constituents,     only: cnst_get_ind, pcnst
    use check_energy,    only: check_energy_chng
    
 ! Arguments
    real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
    type(physics_state), intent(inout) :: state            ! physics state(c) 
    type(physics_tend), intent(inout) :: tend
 
 ! Local workspace
    type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies
    integer :: i, j         ! longitude, latitude,field, and global column indices
    integer :: ncols, k, indw, istep
    real(r8) :: zero(pcols)                    ! array of zeros
    logical  :: lq(pcnst)
                       
 !-----

    zero = 0._r8
    istep=get_nstep()
 
    ! if (masterproc) then
    ! print*, "starting swap to convection state"
    ! !print*, state%sconvforce(1,5)*cpair ! lower number is higher level
    ! print*, state%s(1,5)
    ! !print*, state%uconvforce(1,5)
    ! !print*, state%u(1,5)
    ! endif
    
#if ( defined SPMD )
    call cnst_get_ind('Q',indw)
    lq(:)   =.false.
    lq(indw)=.true.
    call physics_ptend_init(ptend, state%psetcols, "none", ls=.true.,  lu=.true., lv=.true., lq=lq)     !  print *, "=========doloop ====>  ncols= ", ncols

    ncols  = state%ncol
    do i = 1, ncols
    do k=1,pver   
    
        ptend%q(i,k,indw) = ptend%q(i,k,indw) + state%qconvforce(i,k) !/forcingtime
        ptend%u(i,k) = ptend%u(i,k) + state%uconvforce(i,k) !/forcingtime
        ptend%v(i,k) = ptend%v(i,k) + state%vconvforce(i,k) !/forcingtime
        ptend%s(i,k) = ptend%s(i,k) + state%sconvforce(i,k)*cpair !/forcingtime
    
    end do
    end do
 
    !apply tendencies to model
    call physics_update (state, ptend, ztodt, tend) ! this calls ptend deallocate
    call check_energy_chng(state, tend, "convstateswap", istep, ztodt, zero, zero, zero, zero)

    ! if (masterproc) then
    !   print*, "finished swap to convection state"
    !   !print*, state%sconvforce(1,5)*cpair
    !   print*, state%s(1,5)
    !   !print*, state%uconvforce(1,5)
    !   !print*, state%u(1,5)
    ! endif
#endif
 
   end subroutine conv_state_swap_in

   subroutine conv_state_swap_out (ztodt, state,tend)

    !----------------------------------------------------------------------- 
      ! Purpose: 
      !-----------------------------------------------------------------------
    use physconst,    only: cpair
    use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
    use time_manager, only: get_nstep
    use constituents,     only: cnst_get_ind, pcnst
    use check_energy,    only: check_energy_chng
    
 ! Arguments
    real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
    type(physics_state), intent(inout) :: state            ! physics state(c) 
    type(physics_tend), intent(inout) :: tend
 
 ! Local workspace
    type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies
    integer :: i, j         ! longitude, latitude,field, and global column indices
    integer :: ncols, k, indw, istep
    real(r8) :: zero(pcols)                    ! array of zeros
    logical  :: lq(pcnst)
                       
 !-----

    zero = 0._r8
    istep=get_nstep()
 
    ! if (masterproc) then
    ! print*, "starting swap from convection state"
    ! !print*, state%sconvforce(1,5)*cpair ! lower number is higher level
    ! print*, state%s(1,5)
    ! !print*, state%uconvforce(1,5)
    ! !print*, state%u(1,5)
    ! endif
    
#if ( defined SPMD )
    call cnst_get_ind('Q',indw)
    lq(:)   =.false.
    lq(indw)=.true.
    call physics_ptend_init(ptend, state%psetcols, "none", ls=.true.,  lu=.true., lv=.true., lq=lq)     !  print *, "=========doloop ====>  ncols= ", ncols
    ncols  = state%ncol
    
    do i = 1, ncols
    do k=1,pver   
    
        ptend%q(i,k,indw) = ptend%q(i,k,indw) - state%qconvforce(i,k) !/forcingtime
        ptend%u(i,k) = ptend%u(i,k) - state%uconvforce(i,k) !/forcingtime
        ptend%v(i,k) = ptend%v(i,k) - state%vconvforce(i,k) !/forcingtime
        ptend%s(i,k) = ptend%s(i,k) - state%sconvforce(i,k)*cpair !/forcingtime
    
    end do
    end do

 
    !apply tendencies to model
    call physics_update (state, ptend, ztodt, tend) ! this calls ptend deallocate
    call check_energy_chng(state, tend, "convstateswap", istep, ztodt, zero, zero, zero, zero)

    ! if (masterproc) then
    !   print*, "finished swap from convection state"
    !   !print*, state%sconvforce(1,5)*cpair
    !   print*, state%s(1,5)
    !   !print*, state%uconvforce(1,5)
    !   !print*, state%u(1,5)
    ! endif
#endif
    
   end subroutine conv_state_swap_out

end module conv_state_swap 
