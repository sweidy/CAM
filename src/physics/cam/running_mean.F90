module running_mean
!=====================================================================
!
! Purpose: Implement Nudging of the model state of U,V,T,Q, and/or PS
!          toward specified values from analyses. 
!
! Author: Patrick Callaghan
!
! Description:
!    
!    This module assumes that the user has {U,V,T,Q,PS} values from analyses 
!    which have been preprocessed onto the current model grid and adjusted 
!    for differences in topography. It is also assumed that these resulting 
!    values and are stored in individual files which are indexed with respect 
!    to year, month, day, and second of the day. When the model is inbetween 
!    the given begining and ending times, a relaxation forcing is added to 
!    nudge the model toward the analyses values determined from the forcing 
!    option specified. After the model passes the ending analyses time, the 
!    forcing discontinues.
!
!    Some analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the nudging 
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap. 
!
!    Currently, the nudging module is set up to accomodate nudging of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented. 
!
!    The nudging of the model toward the analyses data is controlled by 
!    the 'nudging_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which nudging is applied, the strength of the nudging
!    tendencies, and its spatial distribution. 
!
!    FORCING:
!    --------
!    Nudging tendencies are applied as a relaxation force between the current
!    model state values and target state values derived from the avalilable
!    analyses. The form of the target values is selected by the 'Running_mean_Force_Opt'
!    option, the timescale of the forcing is determined from the given 
!    'Running_mean_TimeScale_Opt', and the nudging strength Alpha=[0.,1.] for each 
!    variable is specified by the 'Running_mean_Xcoef' values. Where X={U,V,T,Q,PS}
!
!           F_Running_mean = Alpha*((Target-Model(t_curr))/TimeScale
!
!
!    WINDOWING:
!    ----------
!    The region of applied nudging can be limited using Horizontal/Vertical 
!    window functions that are constructed using a parameterization of the 
!    Heaviside step function. 
!
!    The Heaviside window function is the product of separate horizonal and vertical 
!    windows that are controled via 12 parameters:
!
!        Running_mean_Hwin_lat0:     Specify the horizontal center of the window in degrees. 
!        Running_mean_Hwin_lon0:     The longitude must be in the range [0,360] and the 
!                             latitude should be [-90,+90].
!        Running_mean_Hwin_latWidth: Specify the lat and lon widths of the window as positive 
!        Running_mean_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999) 
!                             renders the window a constant in that direction.
!        Running_mean_Hwin_latDelta: Controls the sharpness of the window transition with a 
!        Running_mean_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step 
!                             function while a large value yeilds a smoother transition.
!        Running_mean_Hwin_Invert  : A logical flag used to invert the horizontal window function 
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Running_mean_Vwin_Lindex:   In the vertical, the window is specified in terms of model 
!        Running_mean_Vwin_Ldelta:   level indcies. The High and Low transition levels should 
!        Running_mean_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also 
!        Running_mean_Vwin_Hdelta:   specified in terms of model indices. For a window function 
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition 
!                             lengths should be set to 0.001 
!        Running_mean_Vwin_Invert  : A logical flag used to invert the vertical window function 
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent 
!                 of the vertical (30 levels):
!                        Running_mean_Hwin_lat0     = 0.         Running_mean_Vwin_Lindex = 0.
!                        Running_mean_Hwin_latWidth = 30.        Running_mean_Vwin_Ldelta = 0.001
!                        Running_mean_Hwin_latDelta = 5.0        Running_mean_Vwin_Hindex = 31.
!                        Running_mean_Hwin_lon0     = 180.       Running_mean_Vwin_Hdelta = 0.001 
!                        Running_mean_Hwin_lonWidth = 999.       Running_mean_Vwin_Invert = .false.
!                        Running_mean_Hwin_lonDelta = 1.0
!                        Running_mean_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Running_mean_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before 
!    running the model. Lookat_NudgeWindow.ncl is a script avalable in the tools directory 
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be 
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!    Forcing contributions are available for history file output by 
!    the names:    {'Running_mean_U','Running_mean_V','Running_mean_T',and 'Running_mean_Q'}
!    The target values that the model state is nudged toward are available for history 
!    file output via the variables:  {'Target_U','Target_V','Target_T',and 'Target_Q'}
!
!    &nudging_nl
!      Running_mean_Model         - LOGICAL toggle to activate nudging.
!                              TRUE  -> Nudging is on.
!                              FALSE -> Nudging is off.                            [DEFAULT]
!
!      Running_mean_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/nudging/ERAI-Data/')
!
!      Running_mean_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. '%y/ERAI_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc')
!
!      Running_mean_Times_Per_Day - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.
!                              8 --> 3 hourly.
!
!      Running_mean_Model_times_Per_Day - INT Number of times to update the model state (used for nudging) 
!                                each day. The value is restricted to be longer than the 
!                                current model timestep and shorter than the analyses 
!                                timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      Running_mean_Beg_Year      - INT nudging begining year.  [1979- ]
!      Running_mean_Beg_Month     - INT nudging begining month. [1-12]
!      Running_mean_Beg_Day       - INT nudging begining day.   [1-31]
!      Running_mean_End_Year      - INT nudging ending year.    [1979-]
!      Running_mean_End_Month     - INT nudging ending month.   [1-12]
!      Running_mean_End_Day       - INT nudging ending day.     [1-31]
!
!      Running_mean_Force_Opt     - INT Index to select the nudging Target for a relaxation 
!                                forcing of the form: 
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)                 [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      Running_mean_TimeScale_Opt - INT Index to select the timescale for nudging.
!                                where (t'==Analysis times ; t==Model Times) 
!
!                              0 -->  TimeScale = 1/Tdlt_Anal                      [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!      Running_mean_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Running_mean_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Running_mean_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Running_mean_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Running_mean_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Nudging of this variable)
!                                         1 == CONSTANT (Spatially Uniform Nudging)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Running_mean_Ucoef         - REAL fractional nudging coeffcient for U. 
!      Running_mean_Vcoef         - REAL fractional nudging coeffcient for V. 
!      Running_mean_Tcoef         - REAL fractional nudging coeffcient for T. 
!      Running_mean_Qcoef         - REAL fractional nudging coeffcient for Q. 
!      Running_mean_PScoef        - REAL fractional nudging coeffcient for PS. 
!
!                                 The strength of the nudging is specified as a fractional 
!                                 coeffcient between [0,1].
!           
!      Running_mean_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Running_mean_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Running_mean_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Running_mean_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Running_mean_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Running_mean_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Running_mean_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      Running_mean_Vwin_Lindex   - REAL LO model index of transition
!      Running_mean_Vwin_Hindex   - REAL HI model index of transition
!      Running_mean_Vwin_Ldelta   - REAL LO transition length 
!      Running_mean_Vwin_Hdelta   - REAL HI transition length 
!      Running_mean_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!    /
!
!================
!
! TO DO:
! -----------
!    ** Implement Ps Nudging????
!          
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,get_curr_date,get_step_size,get_nstep
  use phys_grid   ,   only:scatter_field_to_chunk, gather_chunk_to_field
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default 
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: Running_mean_Model,Running_mean_ON
  public:: running_mean_readnl
  public:: running_mean_init
  public:: running_mean_timestep_init
  public:: running_mean_timestep_tend
  private::running_mean_update_model_fv
  private::running_mean_update_analyses_fv
  private::interpret_filename_climo
  private::running_mean_set_profile
  

  ! running_mean Parameters
  !--------------------
  logical          :: Running_mean_Model       =.false.
  logical          :: Running_mean_ON          =.false.
  logical          :: Running_mean_Initialized =.false.
  character(len=cl):: Target_Path
  character(len=cs):: Target_File,Target_File_Template
  character(len=cl):: Running_mean_Path
  character(len=cs):: Running_mean_File,Running_mean_File_Template
  integer          :: Running_mean_Force_Opt
  integer          :: Running_mean_TimeScale_Opt
  integer          :: Running_mean_TSmode
  integer          :: Running_mean_Times_Per_Day
  integer          :: Running_mean_Model_times_Per_Day
  real(r8)         :: Running_mean_Ucoef,Running_mean_Vcoef
  integer          :: Running_mean_Uprof,Running_mean_Vprof
  real(r8)         :: Running_mean_Qcoef,Running_mean_Tcoef
  integer          :: Running_mean_Qprof,Running_mean_Tprof
  integer          :: Running_mean_Beg_Year ,Running_mean_Beg_Month
  integer          :: Running_mean_Beg_Day  ,Running_mean_Beg_Sec
  integer          :: Running_mean_End_Year ,Running_mean_End_Month
  integer          :: Running_mean_End_Day  ,Running_mean_End_Sec
  integer          :: Running_mean_Curr_Year,Running_mean_Curr_Month
  integer          :: Running_mean_Curr_Day ,Running_mean_Curr_Sec
  integer          :: Running_mean_Next_Year,Running_mean_Next_Month
  integer          :: Running_mean_Next_Day ,Running_mean_Next_Sec
  integer          :: Running_mean_Step
  integer          :: Model_Curr_Year,Model_Curr_Month
  integer          :: Model_Curr_Day ,Model_Curr_Sec
  integer          :: Model_Next_Year,Model_Next_Month
  integer          :: Model_Next_Day ,Model_Next_Sec
  integer          :: Model_Step
  real(r8)         :: Running_mean_Hwin_lat0
  real(r8)         :: Running_mean_Hwin_latWidth
  real(r8)         :: Running_mean_Hwin_latDelta
  real(r8)         :: Running_mean_Hwin_lon0
  real(r8)         :: Running_mean_Hwin_lonWidth
  real(r8)         :: Running_mean_Hwin_lonDelta
  logical          :: Running_mean_Hwin_Invert = .false.
  real(r8)         :: Running_mean_Hwin_lo
  real(r8)         :: Running_mean_Hwin_hi
  real(r8)         :: Running_mean_Vwin_Hindex
  real(r8)         :: Running_mean_Vwin_Hdelta
  real(r8)         :: Running_mean_Vwin_Lindex
  real(r8)         :: Running_mean_Vwin_Ldelta
  logical          :: Running_mean_Vwin_Invert =.false.
  real(r8)         :: Running_mean_Vwin_lo
  real(r8)         :: Running_mean_Vwin_hi
  real(r8)         :: Running_mean_Hwin_latWidthH
  real(r8)         :: Running_mean_Hwin_lonWidthH
  real(r8)         :: Running_mean_Hwin_max
  real(r8)         :: Running_mean_Hwin_min
  integer          :: Running_mean_win_size
  character(len=10) :: Running_mean_weight_type

  ! running_mean State Arrays
  !-----------------------
  integer Running_mean_nlon,Running_mean_nlat,Running_mean_ncol,Running_mean_nlev
  real(r8),allocatable::Target_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_U     (:,:,:,:)  !(pcols,pver,begchunk:endchunk,winsize)
  real(r8),allocatable:: Running_mean_V     (:,:,:,:)  !(pcols,pver,begchunk:endchunk,winsize)
  real(r8),allocatable:: Running_mean_T     (:,:,:,:)  !(pcols,pver,begchunk:endchunk,winsize)
  real(r8),allocatable:: Running_mean_Q     (:,:,:,:)  !(pcols,pver,begchunk:endchunk,winsize)
  real(r8),allocatable:: Running_nudge_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_nudge_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  !real(r8),allocatable:: Running_nudge_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_nudge_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_nudge_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_Utau  (:,:,:)  !(pcols,pver,begchunk:endchunk) 
  real(r8),allocatable:: Running_mean_Vtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_Stau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_Qtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  integer,allocatable :: Running_mean_nstep (:)      ! accumulated number of steps used for running mean to date (time))

  ! running_mean Observation Arrays
  !-----------------------------
  logical :: Target_File_Present
  logical :: Running_mean_File_Present

contains
  !================================================================
  subroutine running_mean_readnl(nlfile)
   ! 
   ! running_mean_READNL: Initialize default values controlling the running_mean 
   !                 process. Then read namelist values to override 
   !                 them.
   !===============================================================
   use ppgrid        ,only: pver
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

   namelist /running_mean_nl/ Running_mean_Model,Target_Path,                       &
                         Target_File_Template,Running_mean_Force_Opt,          &
                         Running_mean_Path,Running_mean_File_Template, &
                         Running_mean_TimeScale_Opt,                          &
                         Running_mean_Times_Per_Day,Running_mean_Model_times_Per_Day,      &
                         Running_mean_Ucoef ,Running_mean_Uprof,                     &
                         Running_mean_Vcoef ,Running_mean_Vprof,                     &
                         Running_mean_Qcoef ,Running_mean_Qprof,                     &
                         Running_mean_Tcoef ,Running_mean_Tprof,                     &
                         Running_mean_Beg_Year,Running_mean_Beg_Month,Running_mean_Beg_Day, &
                         Running_mean_End_Year,Running_mean_End_Month,Running_mean_End_Day, &
                         Running_mean_Hwin_lat0,Running_mean_Hwin_lon0,              &
                         Running_mean_Hwin_latWidth,Running_mean_Hwin_lonWidth,      &
                         Running_mean_Hwin_latDelta,Running_mean_Hwin_lonDelta,      &
                         Running_mean_Hwin_Invert,                            &
                         Running_mean_Vwin_Lindex,Running_mean_Vwin_Hindex,          &
                         Running_mean_Vwin_Ldelta,Running_mean_Vwin_Hdelta,          &
                         Running_mean_Vwin_Invert,                            &
                         Running_mean_win_size, Running_mean_weight_type 

   ! running_mean is NOT initialized yet, For now
   ! running_mean will always begin/end at midnight.
   !--------------------------------------------
   Running_mean_Initialized =.false.
   Running_mean_ON          =.false.
   Running_mean_Beg_Sec=0
   Running_mean_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Running_mean_Model         = .false.
   Target_Path          = '/n/holylfs04/LABS/kuang_lab/Lab/sweidman/MERRA2_OG/MERRA2_f19/'
   Target_File_Template = 'MERRA2_%m%d_%h.nc'
   Running_mean_Force_Opt     = 0
   Running_mean_Path          = '/n/home04/sweidman/holylfs04/IC_CESM2/'
   Running_mean_File_Template = 'cam_running_mean.%m-%d-%s.nc'
   Running_mean_TimeScale_Opt = 0
   Running_mean_TSmode        = 0
   Running_mean_Times_Per_Day = 4
   Running_mean_Model_times_Per_Day = 4
   Running_mean_Ucoef         = 1._r8
   Running_mean_Vcoef         = 1._r8
   Running_mean_Qcoef         = 1._r8
   Running_mean_Tcoef         = 1._r8
   Running_mean_Uprof         = 1
   Running_mean_Vprof         = 1
   Running_mean_Qprof         = 1
   Running_mean_Tprof         = 1
   Running_mean_Beg_Year      = 1980
   Running_mean_Beg_Month     = 1
   Running_mean_Beg_Day       = 1
   Running_mean_End_Year      = 2019
   Running_mean_End_Month     = 12
   Running_mean_End_Day       = 31
   Running_mean_Hwin_lat0     = 0._r8
   Running_mean_Hwin_latWidth = 9999._r8
   Running_mean_Hwin_latDelta = 1.0_r8
   Running_mean_Hwin_lon0     = 180._r8
   Running_mean_Hwin_lonWidth = 9999._r8
   Running_mean_Hwin_lonDelta = 1.0_r8
   Running_mean_Hwin_Invert   = .false.
   Running_mean_Hwin_lo       = 0.0_r8
   Running_mean_Hwin_hi       = 1.0_r8
   Running_mean_Vwin_Hindex   = float(pver+1)
   Running_mean_Vwin_Hdelta   = 0.001_r8
   Running_mean_Vwin_Lindex   = 0.0_r8
   Running_mean_Vwin_Ldelta   = 0.001_r8
   Running_mean_Vwin_Invert   = .false.
   Running_mean_Vwin_lo       = 0.0_r8
   Running_mean_Vwin_hi       = 1.0_r8
   Running_mean_win_size           = 15
   Running_mean_weight_type        = 'uniform'


   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'running_mean_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,running_mean_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('running_mean_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Running_mean_Hwin_Invert) then
     Running_mean_Hwin_lo = 1.0_r8
     Running_mean_Hwin_hi = 0.0_r8
   else
     Running_mean_Hwin_lo = 0.0_r8
     Running_mean_Hwin_hi = 1.0_r8
   endif

   if(Running_mean_Vwin_Invert) then
     Running_mean_Vwin_lo = 1.0_r8
     Running_mean_Vwin_hi = 0.0_r8
   else
     Running_mean_Vwin_lo = 0.0_r8
     Running_mean_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values 
   !----------------------------------
   if((Running_mean_Hwin_lat0.lt.-90._r8).or.(Running_mean_Hwin_lat0.gt.+90._r8)) then
     write(iulog,*) 'running_mean: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'running_mean:  Running_mean_Hwin_lat0=',Running_mean_Hwin_lat0
     call endrun('running_mean_readnl:: ERROR in namelist')
   endif

   if((Running_mean_Hwin_lon0.lt.0._r8).or.(Running_mean_Hwin_lon0.ge.360._r8)) then
     write(iulog,*) 'running_mean: Window lon0 must be in [0,+360)'
     write(iulog,*) 'running_mean:  Running_mean_Hwin_lon0=',Running_mean_Hwin_lon0
     call endrun('running_mean_readnl:: ERROR in namelist')
   endif

   if((Running_mean_Vwin_Lindex.gt.Running_mean_Vwin_Hindex)                         .or. &
      (Running_mean_Vwin_Hindex.gt.float(pver+1)).or.(Running_mean_Vwin_Hindex.lt.0._r8).or. &
      (Running_mean_Vwin_Lindex.gt.float(pver+1)).or.(Running_mean_Vwin_Lindex.lt.0._r8)   ) then
     write(iulog,*) 'running_mean: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'running_mean: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'running_mean: Lindex must be LE than Hindex'
     write(iulog,*) 'running_mean:  Running_mean_Vwin_Lindex=',Running_mean_Vwin_Lindex
     write(iulog,*) 'running_mean:  Running_mean_Vwin_Hindex=',Running_mean_Vwin_Hindex
     call endrun('running_mean_readnl:: ERROR in namelist')
   endif

   if((Running_mean_Hwin_latDelta.le.0._r8).or.(Running_mean_Hwin_lonDelta.le.0._r8).or. &
      (Running_mean_Vwin_Hdelta  .le.0._r8).or.(Running_mean_Vwin_Ldelta  .le.0._r8)    ) then
     write(iulog,*) 'running_mean: Window Deltas must be positive'
     write(iulog,*) 'running_mean:  Running_mean_Hwin_latDelta=',Running_mean_Hwin_latDelta
     write(iulog,*) 'running_mean:  Running_mean_Hwin_lonDelta=',Running_mean_Hwin_lonDelta
     write(iulog,*) 'running_mean:  Running_mean_Vwin_Hdelta=',Running_mean_Vwin_Hdelta
     write(iulog,*) 'running_mean:  Running_mean_Vwin_Ldelta=',Running_mean_Vwin_Ldelta
     call endrun('running_mean_readnl:: ERROR in namelist')

   endif

   if((Running_mean_Hwin_latWidth.le.0._r8).or.(Running_mean_Hwin_lonWidth.le.0._r8)) then
     write(iulog,*) 'running_mean: Window widths must be positive'
     write(iulog,*) 'running_mean:  Running_mean_Hwin_latWidth=',Running_mean_Hwin_latWidth
     write(iulog,*) 'running_mean:  Running_mean_Hwin_lonWidth=',Running_mean_Hwin_lonWidth
     call endrun('running_mean_readnl:: ERROR in namelist')
   endif

   if ( (trim(Running_mean_weight_type) /= "uniform") .and. &
      (trim(Running_mean_weight_type) /= "gaussian") .and. &
      (trim(Running_mean_weight_type) /= "triangular") ) then

     write(iulog,*) "running_mean: Invalid Running_mean_weight_type specified"
     write(iulog,*) "running_mean: Running_mean_weight_type = ", trim(Running_mean_weight_type)
     write(iulog,*) "running_mean: Must be one of: 'uniform', 'gaussian', 'triangular'"
     call endrun("running_mean_readnl:: ERROR in namelist (invalid Running_mean_weight_type)")
   end if

   if (Running_mean_win_size < 1 .or. mod(Running_mean_win_size,2) /= 1) then
     write(iulog,*) "running_mean: Invalid Running_mean_win_size specified"
     write(iulog,*) "running_mean: Running_mean_win_size must be positive and odd"
     call endrun("running_mean_readnl:: ERROR in namelist (invalid Running_mean_win_size)")
   end if

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Target_Path         ,len(Target_Path)         ,mpichar,0,mpicom)
   call mpibcast(Target_File_Template,len(Target_File_Template),mpichar,0,mpicom)
   call mpibcast(Running_mean_Path         ,len(Running_mean_Path)         ,mpichar,0,mpicom)
   call mpibcast(Running_mean_File_Template,len(Running_mean_File_Template),mpichar,0,mpicom)
   call mpibcast(Running_mean_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_Force_Opt    , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_TimeScale_Opt, 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_TSmode       , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Model_times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Ucoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Vcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Tcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Uprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Vprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Tprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Qprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Beg_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_End_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_End_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_End_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_End_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Vwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_win_size,   1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_weight_type,   1, mpilog, 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! running_mean_readnl
  !================================================================


  !================================================================
  subroutine running_mean_init
   ! 
   ! running_mean_INIT: Allocate space and initialize running_mean values
   !===============================================================
   use ppgrid        ,only: pver,pcols,begchunk,endchunk
   use error_messages,only: alloc_err
   use dycore        ,only: dycore_is
   use dyn_grid      ,only: get_horiz_grid_dim_d
   use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
   use cam_history   ,only: addfld
   use shr_const_mod ,only: SHR_CONST_PI
   use filenames     ,only: interpret_filename_spec

   ! Local values
   !----------------
   integer  Year,Month,Day,Sec
   integer  YMD1,YMD
   logical  After_Beg,Before_End
   integer  istat,lchnk,ncol,icol,ilev
   integer  hdim1_d,hdim2_d
   integer  dtime
   real(r8) rlat,rlon
   real(r8) Wprof(pver)
   real(r8) lonp,lon0,lonn,latp,lat0,latn
   real(r8) Val1_p,Val2_p,Val3_p,Val4_p
   real(r8) Val1_0,Val2_0,Val3_0,Val4_0
   real(r8) Val1_n,Val2_n,Val3_n,Val4_n
   integer  nn
   integer  modstep

   ! Get the time step size
   !------------------------
   dtime = get_step_size()

   ! Allocate Space for running_mean data arrays
   !-----------------------------------------
   allocate(Target_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Target_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Target_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Target_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Target_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Target_Q',pcols*pver*((endchunk-begchunk)+1))

   allocate(Model_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Model_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Model_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Model_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Model_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Model_Q',pcols*pver*((endchunk-begchunk)+1))

   allocate(Running_nudge_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_nudge_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_nudge_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_nudge_V',pcols*pver*((endchunk-begchunk)+1))
   !allocate(Running_nudge_T(pcols,pver,begchunk:endchunk),stat=istat)
   !call alloc_err(istat,'running_mean_init','Running_nudge_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_nudge_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_nudge_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_nudge_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_nudge_Q',pcols*pver*((endchunk-begchunk)+1))

   allocate(Running_mean_U(pcols,pver,begchunk:endchunk,Running_mean_win_size),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_U',pcols*pver*((endchunk-begchunk)+1)*Running_mean_win_size)
   allocate(Running_mean_V(pcols,pver,begchunk:endchunk,Running_mean_win_size),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_V',pcols*pver*((endchunk-begchunk)+1)*Running_mean_win_size)
   allocate(Running_mean_T(pcols,pver,begchunk:endchunk,Running_mean_win_size),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_T',pcols*pver*((endchunk-begchunk)+1)*Running_mean_win_size)
   allocate(Running_mean_Q(pcols,pver,begchunk:endchunk,Running_mean_win_size),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Q',pcols*pver*((endchunk-begchunk)+1)*Running_mean_win_size)

   allocate(Running_mean_nstep(Running_mean_win_size),stat=istat)
   call alloc_err(istat, 'running_mean_init', 'Running_mean_nstep',Running_mean_win_size)

   ! Allocate Space for spatial dependence of 
   ! running_mean Coefs and running_mean Forcing.
   !-------------------------------------------
   allocate(Running_mean_Utau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Utau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_Vtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Vtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_Stau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Stau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Qtau',pcols*pver*((endchunk-begchunk)+1))
   
   allocate(Running_mean_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Ustep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Vstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Sstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'running_mean_init','Running_mean_Qstep',pcols*pver*((endchunk-begchunk)+1))
  

   ! Register output fields with the cam history module
   !-----------------------------------------------------
   call addfld( 'Running_nudge_U',(/ 'lev' /),'A','m/s/s'  ,'U running_mean nudging Tendency')
   call addfld( 'Running_nudge_V',(/ 'lev' /),'A','m/s/s'  ,'V running_mean nudging Tendency')
   call addfld( 'Running_nudge_T',(/ 'lev' /),'A','K/s'    ,'T running_mean nudging Tendency')
   call addfld( 'Running_nudge_Q',(/ 'lev' /),'A','kg/kg/s','Q running_mean nudging Tendency')
   call addfld('Target_U',(/ 'lev' /),'A','m/s'    ,'U running_mean Target'  )
   call addfld('Target_V',(/ 'lev' /),'A','m/s'    ,'V running_mean Target'  )
   call addfld('Target_T',(/ 'lev' /),'A','K'      ,'T running_mean Target'  )
   call addfld('Target_Q',(/ 'lev' /),'A','kg/kg'  ,'Q running_mean Target  ')

   !-----------------------------------------
   ! Values initialized only by masterproc
   !-----------------------------------------
   if(masterproc) then

     ! Set the Stepping intervals for Model and running_mean values
     ! Ensure that the Model_Step is not smaller then one timestep
     !  and not larger then the Running_mean_Step.
     !--------------------------------------------------------
     Model_Step=86400/Running_mean_Model_times_Per_Day
     Running_mean_Step=86400/Running_mean_Times_Per_Day
     if(Model_Step.lt.dtime) then
       write(iulog,*) ' '
       write(iulog,*) 'running_mean: Model_Step cannot be less than a model timestep'
       write(iulog,*) 'running_mean:  Setting Model_Step=dtime , dtime=',dtime
       write(iulog,*) ' '
       Model_Step=dtime
     endif
     if(Model_Step.gt.Running_mean_Step) then
       write(iulog,*) ' '
       write(iulog,*) 'running_mean: Model_Step cannot be more than Running_mean_Step'
       write(iulog,*) 'running_mean:  Setting Model_Step=Running_mean_Step, Running_mean_Step=',Running_mean_Step
       write(iulog,*) ' '
       Model_Step=Running_mean_Step
     endif

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Running_mean_nlon=hdim1_d
     Running_mean_nlat=hdim2_d
     Running_mean_ncol=hdim1_d*hdim2_d
     Running_mean_nlev=pver

     ! Check the time relative to the running_mean window
     !------------------------------------------------
     call get_curr_date(Year,Month,Day,Sec)
     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Running_mean_Beg_Year*10000) + (Running_mean_Beg_Month*100) + Running_mean_Beg_Day
     call timemgr_time_ge(YMD1,Running_mean_Beg_Sec,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Running_mean_End_Year*10000) + (Running_mean_End_Month*100) + Running_mean_End_Day
     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Running_mean_End_Sec,Before_End)
  
     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to 
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Year
       Model_Next_Month=Month
       Model_Next_Day  =Day
       Model_Next_Sec  =(Sec/Model_Step)*Model_Step
       Running_mean_Next_Year =Year
       Running_mean_Next_Month=Month
       Running_mean_Next_Day  =Day
       Running_mean_Next_Sec  =(Sec/Running_mean_Step)*Running_mean_Step
     elseif(.not.After_Beg) then
       ! Set Time indicies to running_mean start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Running_mean_Beg_Year
       Model_Next_Month=Running_mean_Beg_Month
       Model_Next_Day  =Running_mean_Beg_Day
       Model_Next_Sec  =Running_mean_Beg_Sec
       Running_mean_Next_Year =Running_mean_Beg_Year
       Running_mean_Next_Month=Running_mean_Beg_Month
       Running_mean_Next_Day  =Running_mean_Beg_Day
       Running_mean_Next_Sec  =Running_mean_Beg_Sec
     elseif(.not.Before_End) then
       ! running_mean will never occur, so switch it off
       !--------------------------------------------
       Running_mean_Model=.false.
       Running_mean_ON   =.false.
       write(iulog,*) ' '
       write(iulog,*) 'running_mean: WARNING - running_mean has been requested by it will'
       write(iulog,*) 'running_mean:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function  
     !----------------------------------------
     lonp= 180._r8
     lon0=   0._r8
     lonn=-180._r8
     latp=  90._r8-Running_mean_Hwin_lat0
     lat0=   0._r8
     latn= -90._r8-Running_mean_Hwin_lat0
    
     Running_mean_Hwin_lonWidthH=Running_mean_Hwin_lonWidth/2._r8
     Running_mean_Hwin_latWidthH=Running_mean_Hwin_latWidth/2._r8

     Val1_p=(1._r8+tanh((Running_mean_Hwin_lonWidthH+lonp)/Running_mean_Hwin_lonDelta))/2._r8
     Val2_p=(1._r8+tanh((Running_mean_Hwin_lonWidthH-lonp)/Running_mean_Hwin_lonDelta))/2._r8
     Val3_p=(1._r8+tanh((Running_mean_Hwin_latWidthH+latp)/Running_mean_Hwin_latDelta))/2._r8
     Val4_p=(1._r8+tanh((Running_mean_Hwin_latWidthH-latp)/Running_mean_Hwin_latDelta))/2_r8
     Val1_0=(1._r8+tanh((Running_mean_Hwin_lonWidthH+lon0)/Running_mean_Hwin_lonDelta))/2._r8
     Val2_0=(1._r8+tanh((Running_mean_Hwin_lonWidthH-lon0)/Running_mean_Hwin_lonDelta))/2._r8
     Val3_0=(1._r8+tanh((Running_mean_Hwin_latWidthH+lat0)/Running_mean_Hwin_latDelta))/2._r8
     Val4_0=(1._r8+tanh((Running_mean_Hwin_latWidthH-lat0)/Running_mean_Hwin_latDelta))/2._r8

     Val1_n=(1._r8+tanh((Running_mean_Hwin_lonWidthH+lonn)/Running_mean_Hwin_lonDelta))/2._r8
     Val2_n=(1._r8+tanh((Running_mean_Hwin_lonWidthH-lonn)/Running_mean_Hwin_lonDelta))/2._r8
     Val3_n=(1._r8+tanh((Running_mean_Hwin_latWidthH+latn)/Running_mean_Hwin_latDelta))/2._r8
     Val4_n=(1._r8+tanh((Running_mean_Hwin_latWidthH-latn)/Running_mean_Hwin_latDelta))/2._r8

     Running_mean_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Running_mean_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     Target_File_Present=.false.
     Running_mean_File_Present=.false.

     ! Initialization is done, 
     !--------------------------
     Running_mean_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if((.not.dycore_is('UNSTRUCTURED')).and. &
        (.not.dycore_is('EUL')         ).and. &
        (.not.dycore_is('LR')          )      ) then
       call endrun('running_mean IS CURRENTLY ONLY CONFIGURED FOR CAM-SE, FV, or EUL')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL running_mean INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) 'running_mean: Running_mean_Model=',Running_mean_Model
     write(iulog,*) 'running_mean: Target_Path=',Target_Path
     write(iulog,*) 'running_mean: Target_File_Template =',Target_File_Template
     write(iulog,*) 'running_mean: Running_mean_Path=',Running_mean_Path
     write(iulog,*) 'running_mean: Running_mean_File_Template =',Running_mean_File_Template
     write(iulog,*) 'running_mean: Running_mean_Force_Opt=',Running_mean_Force_Opt    
     write(iulog,*) 'running_mean: Running_mean_TimeScale_Opt=',Running_mean_TimeScale_Opt    
     write(iulog,*) 'running_mean: Running_mean_TSmode=',Running_mean_TSmode
     write(iulog,*) 'running_mean: Running_mean_Times_Per_Day=',Running_mean_Times_Per_Day
     write(iulog,*) 'running_mean: Running_mean_Model_times_Per_Day=',Running_mean_Model_times_Per_Day
     write(iulog,*) 'running_mean: Running_mean_Step=',Running_mean_Step
     write(iulog,*) 'running_mean: Model_Step=',Model_Step
     write(iulog,*) 'running_mean: Running_mean_Ucoef  =',Running_mean_Ucoef
     write(iulog,*) 'running_mean: Running_mean_Vcoef  =',Running_mean_Vcoef
     write(iulog,*) 'running_mean: Running_mean_Qcoef  =',Running_mean_Qcoef
     write(iulog,*) 'running_mean: Running_mean_Tcoef  =',Running_mean_Tcoef
     write(iulog,*) 'running_mean: Running_mean_Uprof  =',Running_mean_Uprof
     write(iulog,*) 'running_mean: Running_mean_Vprof  =',Running_mean_Vprof
     write(iulog,*) 'running_mean: Running_mean_Qprof  =',Running_mean_Qprof
     write(iulog,*) 'running_mean: Running_mean_Tprof  =',Running_mean_Tprof
     write(iulog,*) 'running_mean: Running_mean_Beg_Year =',Running_mean_Beg_Year
     write(iulog,*) 'running_mean: Running_mean_Beg_Month=',Running_mean_Beg_Month
     write(iulog,*) 'running_mean: Running_mean_Beg_Day  =',Running_mean_Beg_Day
     write(iulog,*) 'running_mean: Running_mean_End_Year =',Running_mean_End_Year
     write(iulog,*) 'running_mean: Running_mean_End_Month=',Running_mean_End_Month
     write(iulog,*) 'running_mean: Running_mean_End_Day  =',Running_mean_End_Day
     write(iulog,*) 'running_mean: Running_mean_Hwin_lat0     =',Running_mean_Hwin_lat0
     write(iulog,*) 'running_mean: Running_mean_Hwin_latWidth =',Running_mean_Hwin_latWidth
     write(iulog,*) 'running_mean: Running_mean_Hwin_latDelta =',Running_mean_Hwin_latDelta
     write(iulog,*) 'running_mean: Running_mean_Hwin_lon0     =',Running_mean_Hwin_lon0
     write(iulog,*) 'running_mean: Running_mean_Hwin_lonWidth =',Running_mean_Hwin_lonWidth
     write(iulog,*) 'running_mean: Running_mean_Hwin_lonDelta =',Running_mean_Hwin_lonDelta
     write(iulog,*) 'running_mean: Running_mean_Hwin_Invert   =',Running_mean_Hwin_Invert  
     write(iulog,*) 'running_mean: Running_mean_Hwin_lo       =',Running_mean_Hwin_lo
     write(iulog,*) 'running_mean: Running_mean_Hwin_hi       =',Running_mean_Hwin_hi
     write(iulog,*) 'running_mean: Running_mean_Vwin_Hindex   =',Running_mean_Vwin_Hindex
     write(iulog,*) 'running_mean: Running_mean_Vwin_Hdelta   =',Running_mean_Vwin_Hdelta
     write(iulog,*) 'running_mean: Running_mean_Vwin_Lindex   =',Running_mean_Vwin_Lindex
     write(iulog,*) 'running_mean: Running_mean_Vwin_Ldelta   =',Running_mean_Vwin_Ldelta
     write(iulog,*) 'running_mean: Running_mean_Vwin_Invert   =',Running_mean_Vwin_Invert  
     write(iulog,*) 'running_mean: Running_mean_Vwin_lo       =',Running_mean_Vwin_lo
     write(iulog,*) 'running_mean: Running_mean_Vwin_hi       =',Running_mean_Vwin_hi
     write(iulog,*) 'running_mean: Running_mean_Hwin_latWidthH=',Running_mean_Hwin_latWidthH
     write(iulog,*) 'running_mean: Running_mean_Hwin_lonWidthH=',Running_mean_Hwin_lonWidthH
     write(iulog,*) 'running_mean: Running_mean_Hwin_max      =',Running_mean_Hwin_max
     write(iulog,*) 'running_mean: Running_mean_Hwin_min      =',Running_mean_Hwin_min
     write(iulog,*) 'running_mean: Running_mean_Initialized   =',Running_mean_Initialized
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Model_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Model         ,            1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_ON            ,            1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_Initialized   ,            1, mpilog, 0, mpicom)
   call mpibcast(Running_mean_ncol          ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_nlev          ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_nlon          ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_nlat          ,            1, mpiint, 0, mpicom)
   call mpibcast(Running_mean_Hwin_max      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_min      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_lonWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Running_mean_Hwin_latWidthH,            1, mpir8 , 0, mpicom)
#endif


!!DIAG
   if(masterproc) then
     write(iulog,*) 'running_mean: running_mean_init() OBS arrays allocated and initialized'
     write(iulog,*) 'running_mean: running_mean_init() SIZE#',(9*pcols*pver*((endchunk-begchunk)+1))
     write(iulog,*) 'running_mean: running_mean_init() MB:',float(8*9*pcols*pver*((endchunk-begchunk)+1))/(1024._r8*1024._r8)
     write(iulog,*) 'running_mean: running_mean_init() pcols=',pcols,' pver=',pver
     write(iulog,*) 'running_mean: running_mean_init() begchunk:',begchunk,' endchunk=',endchunk
     write(iulog,*) 'running_mean: running_mean_init() Running_mean_File_Present=',Running_mean_File_Present
     write(iulog,*) 'running_mean: running_mean_init() Target_File_Present=',Target_File_Present
   endif
!!DIAG

   ! Initialize the analysis filename at the NEXT time for startup.
   ! TODO: this is not the right file to read -- use the previous timestep
   !---------------------------------------------------------------
   
   modstep=int(Running_mean_Next_Sec / 10800)

   Target_File=interpret_filename_climo(Target_File_Template      , &
          mon_spec=Running_mean_Next_Month, &
          day_spec=Running_mean_Next_Day  , &
          hr_spec=modstep, &
          sec_spec=Running_mean_Next_Sec    )
   if(masterproc) then
    write(iulog,*) 'running_mean: Reading analyses:',trim(Target_Path)//trim(Target_File)
   endif

   !----------------------------------------------------------

    call running_mean_update_analyses_fv (trim(Target_Path)//trim(Target_File))

   ! Initialize running_mean Coeffcient profiles in local arrays
   ! Load zeros into running_mean arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call running_mean_set_profile(rlat,rlon,Running_mean_Uprof,Wprof,pver)
       Running_mean_Utau(icol,:,lchnk)=Wprof(:)
       call running_mean_set_profile(rlat,rlon,Running_mean_Vprof,Wprof,pver)
       Running_mean_Vtau(icol,:,lchnk)=Wprof(:)
       call running_mean_set_profile(rlat,rlon,Running_mean_Tprof,Wprof,pver)
       Running_mean_Stau(icol,:,lchnk)=Wprof(:)
       call running_mean_set_profile(rlat,rlon,Running_mean_Qprof,Wprof,pver)
       Running_mean_Qtau(icol,:,lchnk)=Wprof(:)

     end do
     Running_mean_Utau(:ncol,:pver,lchnk) =                             &
     Running_mean_Utau(:ncol,:pver,lchnk) * Running_mean_Ucoef/float(Running_mean_Step)
     Running_mean_Vtau(:ncol,:pver,lchnk) =                             &
     Running_mean_Vtau(:ncol,:pver,lchnk) * Running_mean_Vcoef/float(Running_mean_Step)
     Running_mean_Stau(:ncol,:pver,lchnk) =                             &
     Running_mean_Stau(:ncol,:pver,lchnk) * Running_mean_Tcoef/float(Running_mean_Step)
     Running_mean_Qtau(:ncol,:pver,lchnk) =                             &
     Running_mean_Qtau(:ncol,:pver,lchnk) * Running_mean_Qcoef/float(Running_mean_Step)

     Running_mean_U(:pcols,:pver,lchnk)=0._r8
     Running_mean_V(:pcols,:pver,lchnk)=0._r8
     Running_mean_T(:pcols,:pver,lchnk)=0._r8
     Running_mean_S(:pcols,:pver,lchnk)=0._r8
     Running_mean_Q(:pcols,:pver,lchnk)=0._r8

     Running_mean_Ustep(:pcols,:pver,lchnk)=0._r8
     Running_mean_Vstep(:pcols,:pver,lchnk)=0._r8
     Running_mean_Sstep(:pcols,:pver,lchnk)=0._r8
     Running_mean_Qstep(:pcols,:pver,lchnk)=0._r8
     Target_U(:pcols,:pver,lchnk)=0._r8
     Target_V(:pcols,:pver,lchnk)=0._r8
     Target_T(:pcols,:pver,lchnk)=0._r8
     Target_S(:pcols,:pver,lchnk)=0._r8
     Target_Q(:pcols,:pver,lchnk)=0._r8
   end do

   Running_mean_nstep(:)=0

   ! End Routine
   !------------
   return
  end subroutine ! running_mean_init
  !================================================================


  !================================================================
  subroutine running_mean_timestep_init(phys_state)
   ! 
   ! running_mean_TIMESTEP_INIT: 
   !                 Check the current time and update Model/running_mean 
   !                 arrays when necessary. Toggle the running_mean flag
   !                 when the time is withing the running_mean window.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use dycore       ,only: dycore_is
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use filenames    ,only: interpret_filename_spec
   use ESMF

   ! Arguments
   !-----------
   type(physics_state),intent(in):: phys_state(begchunk:endchunk)

   ! Local values
   !----------------
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD
   logical Update_Model,Update_Running_mean,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,indw

   type(ESMF_Time)         Date1,Date2
   type(ESMF_TimeInterval) DateDiff
   integer                 DeltaT
   real(r8)                Tscale
   real(r8)                Tfrac
   integer                 rc
   integer                 nn
   integer                 kk
   real(r8)                Sbar,Qbar,Wsum
   integer                 modstep, nstep

   real(r8)                wrk  ! nudging timescale, adjusted by running_mean_nstep
   integer                 iw, nstep_idx, half

   ! Check if running_mean is initialized
   !---------------------------------
   if(.not.Running_mean_Initialized) then
     call endrun('running_mean_timestep_init:: running_mean NOT Initialized')
   endif

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(Running_mean_Beg_Year*10000) + (Running_mean_Beg_Month*100) + Running_mean_Beg_Day
   call timemgr_time_ge(YMD1,Running_mean_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Running_mean_End_Year*10000) + (Running_mean_End_Month*100) + Running_mean_End_Day
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Running_mean_End_Sec,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_Next_Year*10000) + (Model_Next_Month*100) + Model_Next_Day
   call timemgr_time_ge(YMD1,Model_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)

   if((Before_End).and.(Update_Model)) then
     ! Increment the Model times by the current interval
     !---------------------------------------------------
     Model_Curr_Year =Model_Next_Year
     Model_Curr_Month=Model_Next_Month
     Model_Curr_Day  =Model_Next_Day
     Model_Curr_Sec  =Model_Next_Sec
     YMD1=(Model_Curr_Year*10000) + (Model_Curr_Month*100) + Model_Curr_Day
     call timemgr_time_inc(YMD1,Model_Curr_Sec,              &
                           YMD2,Model_Next_Sec,Model_Step,0,0)

     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model 
     ! time to a Model_Step after the current time.
     !--------------------------------------------------------------
     call timemgr_time_ge(YMD2,Model_Next_Sec,            &
                          YMD ,Sec           ,Sync_Error)
     if(Sync_Error) then
       Model_Curr_Year =Year
       Model_Curr_Month=Month
       Model_Curr_Day  =Day
       Model_Curr_Sec  =Sec
       call timemgr_time_inc(YMD ,Model_Curr_Sec,              &
                             YMD2,Model_Next_Sec,Model_Step,0,0)
       write(iulog,*) 'running_mean: WARNING - Model_Time Sync ERROR... CORRECTED'
     endif
     Model_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Model_Next_Year*10000)
     Model_Next_Month=(YMD2/100)
     Model_Next_Day  = YMD2-(Model_Next_Month*100)

     ! Load values at Current into the Model arrays
     !-----------------------------------------------
     
     call cnst_get_ind('Q',indw)
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Model_U(:ncol,:pver,lchnk)=phys_state(lchnk)%u(:ncol,:pver)
       Model_V(:ncol,:pver,lchnk)=phys_state(lchnk)%v(:ncol,:pver)
       Model_T(:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
       Model_Q(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,indw)
     end do

      ! DSE tendencies from Temperature only
      !---------------------------------------
      do lchnk=begchunk,endchunk
        ncol=phys_state(lchnk)%ncol
        Model_S(:ncol,:pver,lchnk)=cpair*Model_T(:ncol,:pver,lchnk)
      end do

      Running_mean_File=interpret_filename_spec(Running_mean_File_Template      , &
                                        sec_spec=Model_Curr_Sec    )
      INQUIRE(FILE=trim(Running_mean_Path)//trim(Running_mean_File), EXIST=Running_mean_File_Present)
      
      if (.not. Running_mean_File_Present) print*, 'running mean file missing', Running_mean_File
     
     if(masterproc) then
      write(iulog,*) 'running_mean: Reading analyses:',trim(Running_mean_Path)//trim(Running_mean_File)
     endif

     !----------------------------------------------------------
    call running_mean_update_model_fv (trim(Running_mean_Path)//trim(Running_mean_File), Model_Curr_Month, Running_mean_Curr_Day)

     ! Load Dry Static Energy values for Model
    ! DSE tendencies from Temperature only
    !---------------------------------------
    
    ! only use centered value as nudging, but calculate weight for each
    do iw = 1, Running_mean_win_size
      nstep_idx = Running_mean_nstep(iw)

      if (nstep_idx > 96000) then
          wrk = 1._r8/96000._r8
      elseif (nstep_idx < 1._r8) then ! in case of nstep = 0
          wrk = 1._r8
      else 
          wrk = 1._r8 / (1._r8 + nstep_idx-1)
      endif

      if (masterproc) then
          write(iulog,*) "iw; wrk; nstep", iw, wrk, nstep_idx
      end if

    do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Running_mean_U(:ncol,:pver,lchnk,iw)=Running_mean_U(:ncol,:pver,lchnk,iw)*(1-wrk) + Model_U(:ncol,:pver,lchnk)*wrk
         Running_mean_V(:ncol,:pver,lchnk,iw)=Running_mean_V(:ncol,:pver,lchnk,iw)*(1-wrk) + Model_V(:ncol,:pver,lchnk)*wrk
         Running_mean_T(:ncol,:pver,lchnk,iw)=Running_mean_T(:ncol,:pver,lchnk,iw)*(1-wrk) + Model_S(:ncol,:pver,lchnk)*wrk
         Running_mean_Q(:ncol,:pver,lchnk,iw)=Running_mean_Q(:ncol,:pver,lchnk,iw)*(1-wrk) + Model_Q(:ncol,:pver,lchnk)*wrk
    end do
    end do ! do iw

    half = (Running_mean_win_size - 1)/2 ! center index
    if(masterproc) then
      write(iulog,*) 'half index for nudging: ', half
    endif
    ! save nudging tendency as centered running mean
    do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Running_nudge_U(:ncol,:pver,lchnk)=Running_mean_U(:ncol,:pver,lchnk,half)
       Running_nudge_V(:ncol,:pver,lchnk)=Running_mean_V(:ncol,:pver,lchnk,half)
       Running_nudge_S(:ncol,:pver,lchnk)=Running_mean_T(:ncol,:pver,lchnk,half)*cpair
       Running_nudge_Q(:ncol,:pver,lchnk)=Running_mean_Q(:ncol,:pver,lchnk,half)
     end do


    if (.not. Running_mean_File_Present) print*, 'running mean file missing', Running_mean_File
    if(masterproc) then
      write(iulog,*) 'running_mean: Writing to file:',trim(Running_mean_Path)//trim(Running_mean_File)
    endif

    call running_mean_write_model_fv(trim(Running_mean_Path)//trim(Running_mean_File), Model_Curr_Month, Running_mean_Curr_Day) 

   endif ! ((Before_End).and.(Update_Model)) then

   !----------------------------------------------------------------
   ! When past the NEXT time, Update running_mean Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Running_mean_Next_Year*10000) + (Running_mean_Next_Month*100) + Running_mean_Next_Day
   call timemgr_time_ge(YMD1,Running_mean_Next_Sec,            &
                        YMD ,Sec           ,Update_Running_mean)

   if((Before_End).and.(Update_Running_mean)) then
     ! Increment the Running_mean times by the current interval
     !---------------------------------------------------
     Running_mean_Curr_Year =Running_mean_Next_Year
     Running_mean_Curr_Month=Running_mean_Next_Month
     Running_mean_Curr_Day  =Running_mean_Next_Day
     Running_mean_Curr_Sec  =Running_mean_Next_Sec
     YMD1=(Running_mean_Curr_Year*10000) + (Running_mean_Curr_Month*100) + Running_mean_Curr_Day
     call timemgr_time_inc(YMD1,Running_mean_Curr_Sec,              &
                           YMD2,Running_mean_Next_Sec,Running_mean_Step,0,0)
     Running_mean_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Running_mean_Next_Year*10000)
     Running_mean_Next_Month=(YMD2/100)
     Running_mean_Next_Day  = YMD2-(Running_mean_Next_Month*100)

     ! Set the analysis filename at the NEXT time. (MERRA)
     !---------------------------------------------------------------
     modstep=int(Running_mean_Next_Sec / 10800)
     Target_File=interpret_filename_climo(Target_File_Template      , &
          mon_spec=Running_mean_Next_Month, &
          day_spec=Running_mean_Next_Day  , &
          hr_spec=modstep, &
          sec_spec=Running_mean_Next_Sec    )

      if(masterproc) then
        write(iulog,*) trim(Target_Path)//trim(Target_File)
      endif
      
      INQUIRE(FILE=trim(Target_Path)//trim(Target_File), EXIST=Target_File_Present)
      if (.not. Target_File_Present) print*, 'running_mean target file missing', Target_File

     !----------------------------------------------------------
    call running_mean_update_analyses_fv (trim(Target_Path)//trim(Target_File))
   endif ! ((Before_End).and.(Update_Running_mean)) then

   !----------------------------------------------------------------
   ! Toggle running_mean flag when the time interval is between 
   ! beginning and ending times, and all of the analyses files exist.
   !----------------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
       Running_mean_ON=Target_File_Present
   else
     Running_mean_ON=.false.
   endif

   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if((Before_End).and.((Update_Running_mean).or.(Update_Model))) then

     ! Now Load the Target values for running_mean tendencies
     !---------------------------------------------------
     ! Target is OBS data at NEXT time
     !----------------------------------
     ! Now load Dry Static Energy values for Target
       ! DSE tendencies from Temperature only
       !---------------------------------------
      do lchnk=begchunk,endchunk
        ncol=phys_state(lchnk)%ncol
        Target_S(:ncol,:pver,lchnk)=cpair*Target_T(:ncol,:pver,lchnk)
      end do

     ! Set Tscale for the specified Forcing Option 
     !-----------------------------------------------
     if(Running_mean_TimeScale_Opt.eq.0) then
       Tscale=1._r8
     elseif(Running_mean_TimeScale_Opt.eq.1) then
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Running_mean_Next_Year,MM=Running_mean_Next_Month, &
                               DD=Running_mean_Next_Day , S=Running_mean_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tscale=float(Running_mean_Step)/float(DeltaT)
     else
       write(iulog,*) 'running_mean: Unknown Running_mean_TimeScale_Opt=',Running_mean_TimeScale_Opt
       call endrun('running_mean_timestep_init:: ERROR unknown running_mean_TimeScale_Opt')
     endif

     ! Update the running_mean tendencies with center idx
     !--------------------------------
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Running_mean_Ustep(:ncol,:pver,lchnk)=(  Target_U(:ncol,:pver,lchnk)      &
                                         -Running_nudge_U(:ncol,:pver,lchnk))     &
                                      *Tscale*Running_mean_Utau(:ncol,:pver,lchnk)
       Running_mean_Vstep(:ncol,:pver,lchnk)=(  Target_V(:ncol,:pver,lchnk)      &
                                         -Running_nudge_V(:ncol,:pver,lchnk))     &
                                      *Tscale*Running_mean_Vtau(:ncol,:pver,lchnk)
       Running_mean_Sstep(:ncol,:pver,lchnk)=(  Target_S(:ncol,:pver,lchnk)      &
                                         -Running_nudge_S(:ncol,:pver,lchnk))     &
                                      *Tscale*Running_mean_Stau(:ncol,:pver,lchnk)
       Running_mean_Qstep(:ncol,:pver,lchnk)=(  Target_Q(:ncol,:pver,lchnk)      &
                                         -Running_nudge_Q(:ncol,:pver,lchnk))     &
                                      *Tscale*Running_mean_Qtau(:ncol,:pver,lchnk)
     end do

     if (masterproc) then
        write(iulog,*) 'day, sec', Running_mean_Curr_Day, Running_mean_Curr_Sec
        write(iulog,*) 'Running_mean_Utau(1,20,1) = ', Running_mean_Utau(1,20,begchunk)
        write(iulog,*) 'Target_U(1,20,1) = ', Target_U(1,20,begchunk)
        write(iulog,*) 'Running_mean_U(1,20,1) = ', Running_mean_U(1,20,begchunk) 
        write(iulog,*) 'Running_mean_Ustep(1,20,1) = ', Running_mean_Ustep(1,20,begchunk)
     end if

     !******************
     ! DIAG
     !******************
!    if(masterproc) then
!      write(iulog,*) 'PFC: Target_T(1,:pver,begchunk)=',Target_T(1,:pver,begchunk)  
!      write(iulog,*) 'PFC:  Model_T(1,:pver,begchunk)=',Model_T(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Target_S(1,:pver,begchunk)=',Target_S(1,:pver,begchunk)  
!      write(iulog,*) 'PFC:  Model_S(1,:pver,begchunk)=',Model_S(1,:pver,begchunk)
!      write(iulog,*) 'PFC:      Target_PS(1,begchunk)=',Target_PS(1,begchunk)  
!      write(iulog,*) 'PFC:       Model_PS(1,begchunk)=',Model_PS(1,begchunk)
!      write(iulog,*) 'PFC: Running_mean_Sstep(1,:pver,begchunk)=',Running_mean_Sstep(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Running_mean_Xstep arrays updated:'
!    endif
   endif ! ((Before_End).and.((Update_Running_mean).or.(Update_Model))) then

   ! End Routine
   !------------
   return
  end subroutine ! running_mean_timestep_init
  !================================================================


  !================================================================
  subroutine running_mean_timestep_tend(phys_state,phys_tend)
   ! 
   ! running_mean_TIMESTEP_TEND: 
   !                If running_mean is ON, return the running_mean contributions 
   !                to forcing using the current contents of the Running_mean 
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'running_mean',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Running_mean_ON) then
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
     phys_tend%u(:ncol,:pver)     =Running_mean_Ustep(:ncol,:pver,lchnk)
     phys_tend%v(:ncol,:pver)     =Running_mean_Vstep(:ncol,:pver,lchnk)
     phys_tend%s(:ncol,:pver)     =Running_mean_Sstep(:ncol,:pver,lchnk)
     phys_tend%q(:ncol,:pver,indw)=Running_mean_Qstep(:ncol,:pver,lchnk)

     call outfld( 'Running_mean_U',phys_tend%u                ,pcols,lchnk)
     call outfld( 'Running_mean_V',phys_tend%v                ,pcols,lchnk)
     call outfld( 'Running_mean_T',phys_tend%s/cpair          ,pcols,lchnk)
     call outfld( 'Running_mean_Q',phys_tend%q(1,1,indw)      ,pcols,lchnk)

   endif

   ! End Routine
   !------------
   return
  end subroutine ! running_mean_timestep_tend
  !================================================================


  !================================================================
  subroutine running_mean_update_model_fv(running_mean_file, target_month, target_day)
   ! 
   ! running_mean_UPDATE_ANALYSES_FV: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: running_mean_file
   integer,intent(in):: target_month, target_day  ! adding for centered mean ++SW

   ! Local values
   !-------------
   ! adding time variables ++SW
   integer nlon,nlat,plev,istat,ntime
   integer ncid,varid
   integer ilat,ilon,ilev, iw
   !real(r8) Xmean(Running_mean_nlev,Running_mean_nlat,Running_mean_nlon)
   real(r8) Lat_anal(Running_mean_nlat)
   real(r8) Lon_anal(Running_mean_nlon)
   real(r8) Xtrans(Running_mean_nlon,Running_mean_nlev,Running_mean_nlat)

   ! adding for taking weighted mean ++SW
   real(r8), allocatable :: Time_anal(:)                      ! time dimension (e.g., days since ref)
   real(r8), allocatable :: Uslab(:,:,:,:)                      ! lon x lat x lev for one time
   !real(r8), allocatable :: w(:)                              ! window weights
   integer, allocatable  :: t_indices(:)                      ! actual time indices used

   ! window config 
   integer, dimension(Running_mean_win_size) :: win_offsets 
   !real(r8) :: sigma_days, wsum
   integer :: half, it_center   ! for creating centered window indices
   integer :: itime
   integer, dimension(4) :: start, count ! for reading time dimension
   integer, dimension(1) :: start_t, count_t ! for reading nstep_nudge

   ! calculate doy
   integer, dimension(12) :: cum = (/ 0,31,59,90,120,151,181,212,243,273,304,334 /)
   integer :: ndoys = 365
   integer :: doy, nstep_nudge

   integer :: ndims, dimids(4), dimid, dimlen, k
   character(len=NF90_MAX_NAME) :: dimname

   ! Just read in one file but will select times inside it
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(running_mean_file),EXIST=Running_mean_File_Present)
     write(iulog,*)'running mean nudge: Running_mean_File_Present=',Running_mean_File_Present
   endif
#ifdef SPMD
   call mpibcast(Running_mean_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Running_mean_File_Present) return

   ! masterproc does all of the work here
   !-----------------------------------------
   if(masterproc) then
   
     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(running_mean_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(running_mean_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV lon')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV lon')
     endif

     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV lat')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV lat')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV lev')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV lev')
     endif

     ! added read in time dimension ++SW
     istat=nf90_inq_dimid(ncid,'time',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV time')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=ntime)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV time')
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
     
     if((Running_mean_nlon.ne.nlon).or.(Running_mean_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: running_mean_update_analyses_fv: nlon=',nlon,' Running_mean_nlon=',Running_mean_nlon
      write(iulog,*) 'ERROR: running_mean_update_analyses_fv: nlat=',nlat,' Running_mean_nlat=',Running_mean_nlat
      write(iulog,*) 'ERROR: running_mean_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('running_mean_update_analyses_fv: analyses dimension mismatch')
     endif

     ! allocate extra time variables
     allocate(Time_anal(ntime))
     allocate(Uslab(Running_mean_nlev,Running_mean_nlat,Running_mean_nlon))
     !allocate(w(Running_mean_win_size))
     allocate(t_indices(Running_mean_win_size))

     istat=nf90_inq_varid(ncid,'time',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Time_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

    ! calculate doy for time index of file
    doy = cum(target_month) + target_day
    it_center = doy+1 !modulo(doy-1, ndoys) + 1
    write(iulog,*) 'calculated it_center ', it_center

    ! Define a centered window
    half = (Running_mean_win_size - 1)/2
    win_offsets = [(iw, iw=-half, half)] 

    ! Map offsets to legal indices with wrap-around (use modulo year logic). 
    do iw = 1, Running_mean_win_size
      t_indices(iw) = modulo(it_center - 1 + win_offsets(iw), ntime) + 1  ! Fortran 1-based, modulo wrap
    end do
    write(iulog,*) 'calculated t_indices ', t_indices

    ! ! create weights. needs to be one of uniform, gaussian, or triangular
    ! select case (trim(Running_mean_weight_type))
    ! case ("uniform")
    !   do iw = 1, Running_mean_win_size
    !     w(iw) = 1.0_r8 
    !   end do
    ! case ("gaussian")
    !   ! sigma relative to window half-width; ~ 0.5 works well, using as default
    !   sigma_days = max(1.0e-6, 0.5*real(max(1,half),kind=r8))
    !   do iw = 1, Running_mean_win_size
    !     w(iw) = exp( -0.5 * ( real(win_offsets(iw), r8) / sigma_days )**2 )
    !   end do
    ! case ("triangular")
    !   ! Triangular: weight drops linearly with |offset|, peak at center. Ensure non-negative.
    !   do iw = 1, Running_mean_win_size
    !     w(iw) = real(half + 1 - abs(win_offsets(iw)), r8)   
    !   end do
    ! end select
    ! wsum = sum(w);  if (wsum <= 0.0_r8) then
    !   call endrun('UPDATE_ANALYSES_FV: zero/neg window weight sum')
    ! end if
    ! w = w / wsum   ! normalize

    ! write(iulog,*) 'calculated weights ', w

    ! start reading in nstep
    istat=nf90_inq_varid(ncid,'nstep_nudge',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif

    start_t = t_indices(1)
    count_t = Running_mean_win_size
    
    istat = nf90_get_var(ncid, varid, nstep_nudge, start=start_t, count=count_t)
    if (istat /= NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun('UPDATE_ANALYSES_FV(nstep slab read)')
    end if
    Running_mean_nstep(:) = nstep_nudge(:)
    write(iulog,*) 'running_mean_nstep from file: ',Running_mean_nstep

    ! do iw = 1, Running_mean_win_size
      
    !   itime = t_indices(iw)

    !   start_t = itime
    !   count_t = 1
    !   ! TODO: can't do a start count situation for scalars (but how to read correct timestep)
    !   istat = nf90_get_var(ncid, varid, nstep_nudge, start=start_t)
    !   if (istat /= NF90_NOERR) then
    !     write(iulog,*) nf90_strerror(istat)
    !     call endrun('UPDATE_ANALYSES_FV(nstep slab read)')
    !   end if

    !   Running_mean_nstep(iw) = nstep_nudge
    ! end do
    endif ! (masterproc) then

! TODO: send to other processors?
!#ifdef SPMD
!    call mpibcast(Running_mean_nstep, Running_mean_win_size, mpiint, 0, mpicom)
!#endif

    if(masterproc) then
    ! start reading in U
    istat=nf90_inq_varid(ncid,'U',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_FV')
    endif
    
    ! Accumulate weighted mean by reading one time slab at a time
    ! TODO: change to not taking a mean here. 
    ! need to have start / count be a vector for full window
    ! xanal / xtrans need to have more dimensions to send to chunks
    ! need to figure out scatter to chunks situation and if that's allowed with time
    ! Xmean = 0.0_r8

    ! U saved in lev, lat, lon, time order
    start = (/ 1, 1, 1, t_indices(1) /)
    count = (/ plev, nlat, nlon, Running_mean_win_size/)

    istat = nf90_get_var(ncid, varid, Uslab, start=start, count=count)
    if (istat /= NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      write(iulog,*) 'itime, start, count ', itime, start, count
      call endrun('UPDATE_ANALYSES_FV(U slab read)')
    end if

    ! do iw = 1, Running_mean_win_size
    !   itime = t_indices(iw)
    !   Xmean = Xmean + w(iw) * Uslab
    ! end do
    ! Xmean(lev,lat,lon) is the weighted climatology for the date window.

    do ilat = 1, nlat
      do ilev = 1, plev
        do ilon = 1, nlon
          Xtrans(ilon, ilev, ilat) = Uslab(ilev, ilat, ilon, iw)
        end do
      end do
    end do
    endif ! (masterproc) then
    call scatter_field_to_chunk(1,Running_mean_nlev,1,Running_mean_nlon,Xtrans,   &
                               Running_mean_U(1,1,begchunk))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     ! Accumulate weighted mean by reading one time slab at a time
    Xmean = 0.0_r8
    do iw = 1, Running_mean_win_size
      
      itime = t_indices(iw)

      ! Assuming U(time,lon,lat,lev) dimension order
      start = (/ 1, 1, 1, itime /)
      count = (/ plev, nlat, nlon, 1/)
      
      istat = nf90_get_var(ncid, varid, Uslab, start=start, count=count)
      if (istat /= NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV(V slab read)')
      end if

      Xmean = Xmean + w(iw) * Uslab
    end do

    ! Xmean(lon,lat,lev) is the weighted climatology for the date window.
    do ilat = 1, nlat
      do ilev = 1, plev
        do ilon = 1, nlon
          Xtrans(ilon, ilev, ilat) = Xmean(ilev, ilat, ilon)
        end do
      end do
    end do
    endif ! (masterproc) then
   call scatter_field_to_chunk(1,Running_mean_nlev,1,Running_mean_nlon,Xtrans,   &
                               Running_mean_V(1,1,begchunk))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     ! Accumulate weighted mean by reading one time slab at a time
    Xmean = 0.0_r8
    do iw = 1, Running_mean_win_size
      
      itime = t_indices(iw)

      ! Assuming U(time,lon,lat,lev) dimension order
      start = (/ 1, 1, 1, itime /)
      count = (/ plev, nlat, nlon, 1/)
      
      istat = nf90_get_var(ncid, varid, Uslab, start=start, count=count)
      if (istat /= NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV(T slab read)')
      end if

      Xmean = Xmean + w(iw) * Uslab
    end do

    ! Xmean(lon,lat,lev) is the weighted climatology for the date window.
    do ilat = 1, nlat
      do ilev = 1, plev
        do ilon = 1, nlon
          Xtrans(ilon, ilev, ilat) = Xmean(ilev, ilat, ilon)
        end do
      end do
    end do
    endif ! (masterproc) then
   call scatter_field_to_chunk(1,Running_mean_nlev,1,Running_mean_nlon,Xtrans,   &
                              Running_mean_T(1,1,begchunk))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     ! Accumulate weighted mean by reading one time slab at a time
    Xmean = 0.0_r8
    do iw = 1, Running_mean_win_size
      
      itime = t_indices(iw)

      ! Assuming U(time,lon,lat,lev) dimension order
      start = (/ 1, 1, 1, itime /)
      count = (/ plev, nlat, nlon, 1/)
      
      istat = nf90_get_var(ncid, varid, Uslab, start=start, count=count)
      if (istat /= NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV(Q slab read)')
      end if

      Xmean = Xmean + w(iw) * Uslab
    end do

    ! Xmean(lon,lat,lev) is the weighted climatology for the date window.
    do ilat = 1, nlat
      do ilev = 1, plev
        do ilon = 1, nlon
          Xtrans(ilon, ilev, ilat) = Xmean(ilev, ilat, ilon)
        end do
      end do
    end do

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Running_mean_nlev,1,Running_mean_nlon,Xtrans,   &
                               Running_mean_Q(1,1,begchunk))

   if (allocated(Time_anal)) deallocate(Time_anal)
   if (allocated(Uslab))     deallocate(Uslab)
   if (allocated(w))         deallocate(w)
   if (allocated(t_indices)) deallocate(t_indices)


   ! End Routine
   !------------
   return
  end subroutine ! running_mean_update_model_fv
  !================================================================

  !================================================================
  subroutine running_mean_update_analyses_fv(anal_file)
   ! 
   ! running_mean_UPDATE_ANALYSES_FV: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Running_mean_nlon,Running_mean_nlat,Running_mean_nlev)
   real(r8) Lat_anal(Running_mean_nlat)
   real(r8) Lon_anal(Running_mean_nlon)
   real(r8) Xtrans(Running_mean_nlon,Running_mean_nlev,Running_mean_nlat)
   integer  nn

   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Target_File_Present)
     write(iulog,*)'running_mean: Target_File_Present=',Target_File_Present
   endif
#ifdef SPMD
   call mpibcast(Target_File_Present, 1, mpilog, 0, mpicom)

#endif
   if(.not.Target_File_Present) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then
   
     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
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

     if((Running_mean_nlon.ne.nlon).or.(Running_mean_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: running_mean_update_analyses_fv: nlon=',nlon,' Running_mean_nlon=',Running_mean_nlon
      write(iulog,*) 'ERROR: running_mean_update_analyses_fv: nlat=',nlat,' Running_mean_nlat=',Running_mean_nlat
      write(iulog,*) 'ERROR: running_mean_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('running_mean_update_analyses_fv: analyses dimension mismatch')
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
   call scatter_field_to_chunk(1,Running_mean_nlev,1,Running_mean_nlon,Xtrans,   &
                               Target_U(1,1,begchunk))

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
   call scatter_field_to_chunk(1,Running_mean_nlev,1,Running_mean_nlon,Xtrans,   &
                               Target_V(1,1,begchunk))

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
   call scatter_field_to_chunk(1,Running_mean_nlev,1,Running_mean_nlon,Xtrans,   &
                               Target_T(1,1,begchunk))

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

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Running_mean_nlev,1,Running_mean_nlon,Xtrans,   &
                               Target_Q(1,1,begchunk))

   ! End Routine
   !------------
   return
  end subroutine ! running_mean_update_analyses_fv
  !================================================================

    !================================================================
  subroutine running_mean_write_model_fv(running_mean_file, target_month, target_day)
   ! 
   ! running_mean_UPDATE_ANALYSES_FV: 
   !                 Open the given analyses data file, write out in 
   !                 U,V,T,Q, and PS values and after gathering from chunks
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: running_mean_file
   integer, intent(in)        :: target_month, target_day

   ! Local values
   !-------------
   ! adding time variables ++SW
   integer nlon,nlat,plev,istat,ntime
   integer ncid,varid,varid_t
   integer ilat,ilon,ilev, iw
   real(r8) Xtrans(Running_mean_nlon,Running_mean_nlev,Running_mean_nlat)

   ! adding for taking weighted mean ++SW
   real, allocatable :: Time_anal(:)                      ! time dimension (e.g., days since ref)
   real, allocatable :: Uslab(:,:,:)                      ! lon x lat x lev for one timestep
   real, allocatable :: w(:)                              ! window weights
   integer, allocatable :: t_indices(:)                   ! actual time indices used

   ! window config (11-day window)
   integer, dimension(Running_mean_win_size) :: win_offsets 
   real :: sigma_days, wsum
   integer :: half, it_center   ! for creating centered window indices
   integer :: itime
   integer, dimension(4) :: start, count

   ! calculate doy
   integer, dimension(12) :: cum = (/ 0,31,59,90,120,151,181,212,243,273,304,334 /)
   integer :: ndoys = 365
   integer :: doy
   integer :: ndims, dimids(4), dimlen, k, dimid
   character(len=NF90_MAX_NAME) :: dimname

   ! Just read in one file but will select times inside it
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(running_mean_file),EXIST=Running_mean_File_Present)
     write(iulog,*)'running mean nudge: Running_mean_File_Present=',Running_mean_File_Present
   endif
#ifdef SPMD
   call mpibcast(Running_mean_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Running_mean_File_Present) return

   ! masterproc does all of the work here
   !-----------------------------------------
   if(masterproc) then

    ! Open the given file
     !-----------------------
     istat=nf90_open(trim(running_mean_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(running_mean_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
   
   nlon = Running_mean_nlon
   nlat = Running_mean_nlat
   plev = pver

   ! added read in time dimension ++SW
     istat=nf90_inq_dimid(ncid,'time',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=ntime)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! allocate extra time variables
     allocate(Time_anal(ntime))
     allocate(Uslab(plev,Running_mean_nlat,Running_mean_nlon))
     allocate(w(Running_mean_win_size))
     allocate(t_indices(Running_mean_win_size))

     istat=nf90_inq_varid(ncid,'time',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Time_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

    ! calculate doy for time index of file
    doy = cum(target_month) + target_day
    it_center = doy+1 !modulo(doy-1, ndoys) + 1
    write(iulog,*) 'calculated it_center ', it_center

    ! Define a centered window
    half = (Running_mean_win_size - 1)/2
    win_offsets = [(iw, iw=-half, half)] 

    ! Map offsets to legal indices with wrap-around (use modulo year logic). 
    do iw = 1, Running_mean_win_size
      t_indices(iw) = modulo(it_center - 1 + win_offsets(iw), ntime) + 1  ! Fortran 1-based, modulo wrap
    end do
    write(iulog,*) 'calculated t_indices ', t_indices

    ! create weights. needs to be one of uniform, gaussian, or triangular
    select case (trim(Running_mean_weight_type))
    case ("uniform")
      do iw = 1, Running_mean_win_size
        w(iw) = 1.0_r8 
      end do
    case ("gaussian")
      ! sigma relative to window half-width; ~ 0.5 works well, using as default
      sigma_days = max(1.0e-6, 0.5*real(max(1,half),kind=r8))
      do iw = 1, Running_mean_win_size
        w(iw) = exp( -0.5 * ( real(win_offsets(iw), r8) / sigma_days )**2 )
      end do
    case ("triangular")
      ! Triangular: weight drops linearly with |offset|, peak at center. Ensure non-negative.
      do iw = 1, Running_mean_win_size
        w(iw) = real(half + 1 - abs(win_offsets(iw)), r8)   
      end do
    end select
    wsum = sum(w);  if (wsum <= 0.0_r8) then
      call endrun('UPDATE_ANALYSES_FV: zero/neg window weight sum')
    end if
    w = w / wsum   ! normalize
    endif ! masterproc
    
    ! Zeyuan Hu 12/23/2024: gather global state variables
    ! Modified by Sarah Weidman 9/2025
    !---------------------------------------------------

    call gather_chunk_to_field(1,Running_mean_nlev,1,Running_mean_nlon,Running_mean_U,Xtrans)
    
    if (masterproc) then
      ! open file and find U variable
      istat = nf90_open(trim(running_mean_file), NF90_WRITE, ncid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif
      istat = nf90_inq_varid(ncid, 'U', varid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif

      istat = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif
      do k=1, ndims
        dimid = dimids(k)
        istat = nf90_inquire_dimension(ncid, dimid, name=dimname, len=dimlen)
        write(iulog,*) 'U dim', k, ':', trim(dimname), ' len=', dimlen
      end do

      do iw = 1, Running_mean_win_size
        do ilat=1,nlat
        do ilev=1,plev
        do ilon=1,nlon
          ! weight the file by time within window
          Uslab(ilev,ilat,ilon)=Xtrans(ilon,ilev,ilat)*w(iw)
        end do
        end do
        end do
        ! write timestep to file
        itime = t_indices(iw)
        start = (/ 1, 1, 1, itime /)
        count = (/ plev, nlat, nlon, 1/)
        istat = nf90_put_var(ncid, varid, Uslab, start=start, count=count)
        if(istat.ne.NF90_NOERR) then
          write(iulog,*) nf90_strerror(istat)
          call endrun ('UPDATE_ANALYSES_FV')
        endif

      end do
      istat = nf90_sync(ncid)
      if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV')
      endif

      istat = nf90_close(ncid)
      if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV')
      endif
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Running_mean_nlev,1,Running_mean_nlon,Running_mean_V,Xtrans)
    
    if (masterproc) then
      ! open file and find V variable
      istat = nf90_open(trim(running_mean_file), NF90_WRITE, ncid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif
      istat = nf90_inq_varid(ncid, 'V', varid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif

      do iw = 1, Running_mean_win_size
        do ilat=1,nlat
        do ilev=1,plev
        do ilon=1,nlon
          ! weight the file by time within window
          Uslab(ilev,ilat,ilon)=Xtrans(ilon,ilev,ilat)*w(iw)
        end do
        end do
        end do
        ! write timestep to file
        itime = t_indices(iw)
        start = (/ 1, 1, 1, itime /)
        count = (/ plev, nlat, nlon, 1/)
        istat = nf90_put_var(ncid, varid, Uslab, start=start, count=count)
        if(istat.ne.NF90_NOERR) then
          write(iulog,*) nf90_strerror(istat)
          call endrun ('UPDATE_ANALYSES_FV')
        endif
      end do
      istat = nf90_sync(ncid)
      if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV')
      endif

      istat = nf90_close(ncid)
      if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV')
      endif
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Running_mean_nlev,1,Running_mean_nlon,Running_mean_T,Xtrans)
    
    if (masterproc) then
      ! open file and find T variable
      istat = nf90_open(trim(running_mean_file), NF90_WRITE, ncid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif
      istat = nf90_inq_varid(ncid, 'T', varid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif

      do iw = 1, Running_mean_win_size
        do ilat=1,nlat
        do ilev=1,plev
        do ilon=1,nlon
          ! weight the file by time within window
          Uslab(ilev,ilat,ilon)=Xtrans(ilon,ilev,ilat)*w(iw)
        end do
        end do
        end do
        ! write timestep to file
        itime = t_indices(iw)
        start = (/ 1, 1, 1, itime /)
        count = (/ plev, nlat, nlon, 1/)
        istat = nf90_put_var(ncid, varid, Uslab, start=start, count=count)
        if(istat.ne.NF90_NOERR) then
          write(iulog,*) nf90_strerror(istat)
          call endrun ('UPDATE_ANALYSES_FV')
        endif
      end do
      istat = nf90_sync(ncid)
      if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV')
      endif

      istat = nf90_close(ncid)
      if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV')
      endif
    endif ! (masterproc) then

    call gather_chunk_to_field(1,Running_mean_nlev,1,Running_mean_nlon,Running_mean_Q,Xtrans)
    
    if (masterproc) then
      ! open file and find Q variable
      istat = nf90_open(trim(running_mean_file), NF90_WRITE, ncid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif
      istat = nf90_inq_varid(ncid, 'Q', varid)
      if(istat.ne.NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun ('UPDATE_ANALYSES_FV')
      endif

      do iw = 1, Running_mean_win_size
        do ilat=1,nlat
        do ilev=1,plev
        do ilon=1,nlon
          ! weight the file by time within window
          Uslab(ilev,ilat,ilon)=Xtrans(ilon,ilev,ilat)*w(iw)
        end do
        end do
        end do
        ! write timestep to file
        itime = t_indices(iw)
        start = (/ 1, 1, 1, itime /)
        count = (/ plev, nlat, nlon, 1/)
        istat = nf90_put_var(ncid, varid, Uslab, start=start, count=count)
        if(istat.ne.NF90_NOERR) then
          write(iulog,*) nf90_strerror(istat)
          call endrun ('UPDATE_ANALYSES_FV')
        endif
      end do
      istat = nf90_sync(ncid)
      if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV')
      endif

      istat = nf90_close(ncid)
      if (istat .ne. NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat)
        call endrun('UPDATE_ANALYSES_FV')
      endif
    endif ! (masterproc) then

   if (allocated(Time_anal)) deallocate(Time_anal)
   if (allocated(Uslab))     deallocate(Uslab)
   if (allocated(w))         deallocate(w)
   if (allocated(t_indices)) deallocate(t_indices)

   ! End Routine
   !------------
   return
  end subroutine ! running_mean_write_model_fv
  
  !================================================================

   character(len=cl) function interpret_filename_climo( filename_spec, case, &
   mon_spec, day_spec, hr_spec, sec_spec )

! Create a filename from a filename specifier. The 
! filename specifyer includes codes for setting things such as the
! month, day, seconds in day, caseid, and tape number. 
!
! Interpret filename specifyer string with: 
!
!      %c for case, 
!      %m for month
!      %d for day
!      %h for modstep
!      %% for the "%" character
!
! If the filename specifyer has spaces " ", they will be trimmed out
! of the resulting filename.

   ! arguments
   character(len=*), intent(in)           :: filename_spec   ! Filename specifier to use
   character(len=*), intent(in), optional :: case            ! Optional casename
   integer         , intent(in), optional :: mon_spec        ! Simulation month
   integer         , intent(in), optional :: day_spec        ! Simulation day
   integer         , intent(in), optional :: hr_spec         ! Modstep
   integer         , intent(in), optional :: sec_spec        ! Simulation seconds of day

   ! Local variables
   integer :: month ! Simulation month
   integer :: day   ! Simulation day
   integer :: ncsec   ! Seconds into current simulation day
   integer :: modstep ! modstep into current simulation day
   character(len=cl) :: string    ! Temporary character string 
   character(len=cl) :: format    ! Format character string 
   integer :: i, n  ! Loop variables
   logical :: done
   !-----------------------------------------------------------------------------


   if ( len_trim(filename_spec) == 0 )then
      call endrun ('INTERPRET_FILENAME_CLIMO: filename specifier is empty')
   end if
   if ( index(trim(filename_spec)," ") /= 0 )then
      call endrun ('INTERPRET_FILENAME_CLIMO: filename specifier can not contain a space:'//trim(filename_spec))
   end if
   !
   ! Determine month, day and sec to put in filename
   !
   if (present(mon_spec) .and. present(day_spec) .and. present(hr_spec) .and. present(sec_spec)) then
      month = mon_spec
      day   = day_spec
      modstep = hr_spec
      ncsec = sec_spec
   end if
   !
   ! Go through each character in the filename specifyer and interpret if special string
   !
   i = 1
   interpret_filename_climo = ''
   do while ( i <= len_trim(filename_spec) )
      !
      ! If following is an expansion string
      !
      if ( filename_spec(i:i) == "%" )then
         i = i + 1
         select case( filename_spec(i:i) )
         case( 'm' )   ! month
            write(string,'(i2.2)') month
         case( 'd' )   ! day
            write(string,'(i2.2)') day
         case( 'h' )   ! 3-hour period
            write(string,'(i1.1)') modstep
         case( 's' )   ! second
            write(string,'(i5.5)') ncsec
         case( '%' )   ! percent character
            string = "%"
         case default
            call endrun ('INTERPRET_FILENAME_CLIMO: Invalid expansion character: '//filename_spec(i:i))
         end select
         !
         ! Otherwise take normal text up to the next "%" character
         !
      else
         n = index( filename_spec(i:), "%" )
         if ( n == 0 ) n = len_trim( filename_spec(i:) ) + 1
         if ( n == 0 ) exit 
         string = filename_spec(i:n+i-2)
         i = n + i - 2
      end if
      if ( len_trim(interpret_filename_climo) == 0 )then
        interpret_filename_climo = trim(string)
      else
         if ( (len_trim(interpret_filename_climo)+len_trim(string)) >= cl )then
            call endrun ('INTERPRET_FILENAME_CLIMO: Resultant filename too long')
         end if
         interpret_filename_climo = trim(interpret_filename_climo) // trim(string)
      end if
      i = i + 1

   end do
   if ( len_trim(interpret_filename_climo) == 0 )then
      call endrun ('INTERPRET_FILENAME_CLIMO: Resulting filename is empty')
   end if

end function interpret_filename_climo

  !================================================================
  subroutine running_mean_set_profile(rlat,rlon,Running_mean_prof,Wprof,nlev)
   ! 
   ! running_mean_SET_PROFILE: for the given lat,lon, and running_mean_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on running_mean strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Running_mean_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   ! Local values
   !----------------
   integer  ilev
   real(r8) Hcoef,latx,lonx,Vmax,Vmin
   real(r8) lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(Running_mean_prof.eq.0) then
     ! No running_mean
     !-------------
     Wprof(:)=0.0_r8
   elseif(Running_mean_prof.eq.1) then
     ! Uniform running_mean
     !-----------------
     Wprof(:)=1.0_r8
   elseif(Running_mean_prof.eq.2) then
     ! Localized running_mean with specified Heaviside window function
     !------------------------------------------------------------
     if(Running_mean_Hwin_max.le.Running_mean_Hwin_min) then
       ! For a constant Horizontal window function, 
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Running_mean_Hwin_lo,Running_mean_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Running_mean_Hwin_lat0
       lonx=rlon-Running_mean_Hwin_lon0
       if(lonx.gt. 180._r8) lonx=lonx-360._r8
       if(lonx.le.-180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Running_mean_Hwin_lonWidthH+lonx)/Running_mean_Hwin_lonDelta
       lon_hi=(Running_mean_Hwin_lonWidthH-lonx)/Running_mean_Hwin_lonDelta
       lat_lo=(Running_mean_Hwin_latWidthH+latx)/Running_mean_Hwin_latDelta
       lat_hi=(Running_mean_Hwin_latWidthH-latx)/Running_mean_Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Running_mean_Hwin_min)/(Running_mean_Hwin_max-Running_mean_Hwin_min)
       Hcoef=(1._r8-Hcoef)*Running_mean_Hwin_lo + Hcoef*Running_mean_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Running_mean_Vwin_Lindex)/Running_mean_Vwin_Ldelta
       lev_hi=(Running_mean_Vwin_Hindex-float(ilev))/Running_mean_Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do 

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax.le.Vmin).or.((Running_mean_Vwin_Hindex.ge.(nlev+1)).and. &
                           (Running_mean_Vwin_Lindex.le. 0      )     )) then
       ! For a constant Vertical window function, 
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Running_mean_Vwin_lo,Running_mean_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Running_mean_Vwin_lo + Wprof(:)*(Running_mean_Vwin_hi-Running_mean_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile 
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('running_mean_set_profile:: Unknown Running_mean_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! running_mean_set_profile
  !================================================================

end module running_mean
