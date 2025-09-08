module nudging
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
!    analyses. The form of the target values is selected by the 'Nudge_Force_Opt'
!    option, the timescale of the forcing is determined from the given 
!    'Nudge_TimeScale_Opt', and the nudging strength Alpha=[0.,1.] for each 
!    variable is specified by the 'Nudge_Xcoef' values. Where X={U,V,T,Q,PS}
!
!           F_nudge = Alpha*((Target-Model(t_curr))/TimeScale
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
!        Nudge_Hwin_lat0:     Specify the horizontal center of the window in degrees. 
!        Nudge_Hwin_lon0:     The longitude must be in the range [0,360] and the 
!                             latitude should be [-90,+90].
!        Nudge_Hwin_latWidth: Specify the lat and lon widths of the window as positive 
!        Nudge_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999) 
!                             renders the window a constant in that direction.
!        Nudge_Hwin_latDelta: Controls the sharpness of the window transition with a 
!        Nudge_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step 
!                             function while a large value yeilds a smoother transition.
!        Nudge_Hwin_Invert  : A logical flag used to invert the horizontal window function 
!                             to get its compliment.(e.g. to nudge outside a given window).
!
!        Nudge_Vwin_Lindex:   In the vertical, the window is specified in terms of model 
!        Nudge_Vwin_Ldelta:   level indcies. The High and Low transition levels should 
!        Nudge_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also 
!        Nudge_Vwin_Hdelta:   specified in terms of model indices. For a window function 
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition 
!                             lengths should be set to 0.001 
!        Nudge_Vwin_Invert  : A logical flag used to invert the vertical window function 
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent 
!                 of the vertical (30 levels):
!                        Nudge_Hwin_lat0     = 0.         Nudge_Vwin_Lindex = 0.
!                        Nudge_Hwin_latWidth = 30.        Nudge_Vwin_Ldelta = 0.001
!                        Nudge_Hwin_latDelta = 5.0        Nudge_Vwin_Hindex = 31.
!                        Nudge_Hwin_lon0     = 180.       Nudge_Vwin_Hdelta = 0.001 
!                        Nudge_Hwin_lonWidth = 999.       Nudge_Vwin_Invert = .false.
!                        Nudge_Hwin_lonDelta = 1.0
!                        Nudge_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply nudging at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Nudge_Hwin_Invert = .true.
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
!    the names:    {'Nudge_U','Nudge_V','Nudge_T',and 'Nudge_Q'}
!    The target values that the model state is nudged toward are available for history 
!    file output via the variables:  {'Target_U','Target_V','Target_T',and 'Target_Q'}
!
!    &nudging_nl
!      Nudge_Model         - LOGICAL toggle to activate nudging.
!                              TRUE  -> Nudging is on.
!                              FALSE -> Nudging is off.                            [DEFAULT]
!
!      Nudge_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/nudging/ERAI-Data/')
!
!      Nudge_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. '%y/ERAI_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc')
!
!      Nudge_Times_Per_Day - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.
!                              8 --> 3 hourly.
!
!      Model_Times_Per_Day - INT Number of times to update the model state (used for nudging) 
!                                each day. The value is restricted to be longer than the 
!                                current model timestep and shorter than the analyses 
!                                timestep. As this number is increased, the nudging
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      Nudge_Beg_Year      - INT nudging begining year.  [1979- ]
!      Nudge_Beg_Month     - INT nudging begining month. [1-12]
!      Nudge_Beg_Day       - INT nudging begining day.   [1-31]
!      Nudge_End_Year      - INT nudging ending year.    [1979-]
!      Nudge_End_Month     - INT nudging ending month.   [1-12]
!      Nudge_End_Day       - INT nudging ending day.     [1-31]
!
!      Nudge_Force_Opt     - INT Index to select the nudging Target for a relaxation 
!                                forcing of the form: 
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)                 [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      Nudge_TimeScale_Opt - INT Index to select the timescale for nudging.
!                                where (t'==Analysis times ; t==Model Times) 
!
!                              0 -->  TimeScale = 1/Tdlt_Anal                      [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!      Nudge_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Nudge_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Nudge_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Nudge_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Nudge_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Nudging of this variable)
!                                         1 == CONSTANT (Spatially Uniform Nudging)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Nudge_Ucoef         - REAL fractional nudging coeffcient for U. 
!      Nudge_Vcoef         - REAL fractional nudging coeffcient for V. 
!      Nudge_Tcoef         - REAL fractional nudging coeffcient for T. 
!      Nudge_Qcoef         - REAL fractional nudging coeffcient for Q. 
!      Nudge_PScoef        - REAL fractional nudging coeffcient for PS. 
!
!                                 The strength of the nudging is specified as a fractional 
!                                 coeffcient between [0,1].
!           
!      Nudge_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Nudge_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Nudge_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Nudge_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Nudge_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Nudge_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Nudge_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      Nudge_Vwin_Lindex   - REAL LO model index of transition
!      Nudge_Vwin_Hindex   - REAL HI model index of transition
!      Nudge_Vwin_Ldelta   - REAL LO transition length 
!      Nudge_Vwin_Hdelta   - REAL HI transition length 
!      Nudge_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
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
  use phys_grid   ,   only:scatter_field_to_chunk
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

  public:: Nudge_Model,Nudge_ON
  public:: nudging_readnl
  public:: nudging_init
  public:: nudging_timestep_init
  public:: nudging_timestep_tend
  private::nudging_update_model_fv
  private::nudging_update_analyses_fv
  private::find_nearest_time_index
  private::interpret_filename_climo
  private::nudging_set_PSprofile
  private::nudging_set_profile
  private::calc_DryStaticEnergy
  

  ! Nudging Parameters
  !--------------------
  logical          :: Nudge_Model       =.false.
  logical          :: Nudge_ON          =.false.
  logical          :: Nudge_Initialized =.false.
  character(len=cl):: Nudge_Path
  character(len=cs):: Nudge_File,Nudge_File_Template
  character(len=cl):: Running_mean_Path
  character(len=cs):: Running_mean_File,Running_mean_File_Template
  integer          :: Nudge_Force_Opt
  integer          :: Nudge_TimeScale_Opt
  integer          :: Nudge_TSmode
  integer          :: Nudge_Times_Per_Day
  integer          :: Model_Times_Per_Day
  real(r8)         :: Nudge_Ucoef,Nudge_Vcoef
  integer          :: Nudge_Uprof,Nudge_Vprof
  real(r8)         :: Nudge_Qcoef,Nudge_Tcoef
  integer          :: Nudge_Qprof,Nudge_Tprof
  real(r8)         :: Nudge_PScoef
  integer          :: Nudge_PSprof
  integer          :: Nudge_Beg_Year ,Nudge_Beg_Month
  integer          :: Nudge_Beg_Day  ,Nudge_Beg_Sec
  integer          :: Nudge_End_Year ,Nudge_End_Month
  integer          :: Nudge_End_Day  ,Nudge_End_Sec
  integer          :: Nudge_Curr_Year,Nudge_Curr_Month
  integer          :: Nudge_Curr_Day ,Nudge_Curr_Sec
  integer          :: Nudge_Next_Year,Nudge_Next_Month
  integer          :: Nudge_Next_Day ,Nudge_Next_Sec
  integer          :: Nudge_Step
  integer          :: Model_Curr_Year,Model_Curr_Month
  integer          :: Model_Curr_Day ,Model_Curr_Sec
  integer          :: Model_Next_Year,Model_Next_Month
  integer          :: Model_Next_Day ,Model_Next_Sec
  integer          :: Model_Step
  real(r8)         :: Nudge_Hwin_lat0
  real(r8)         :: Nudge_Hwin_latWidth
  real(r8)         :: Nudge_Hwin_latDelta
  real(r8)         :: Nudge_Hwin_lon0
  real(r8)         :: Nudge_Hwin_lonWidth
  real(r8)         :: Nudge_Hwin_lonDelta
  logical          :: Nudge_Hwin_Invert = .false.
  real(r8)         :: Nudge_Hwin_lo
  real(r8)         :: Nudge_Hwin_hi
  real(r8)         :: Nudge_Vwin_Hindex
  real(r8)         :: Nudge_Vwin_Hdelta
  real(r8)         :: Nudge_Vwin_Lindex
  real(r8)         :: Nudge_Vwin_Ldelta
  logical          :: Nudge_Vwin_Invert =.false.
  real(r8)         :: Nudge_Vwin_lo
  real(r8)         :: Nudge_Vwin_hi
  real(r8)         :: Nudge_Hwin_latWidthH
  real(r8)         :: Nudge_Hwin_lonWidthH
  real(r8)         :: Nudge_Hwin_max
  real(r8)         :: Nudge_Hwin_min
  integer          :: win_size
  LOGICAL          :: use_gaussian = .true. 

  ! Nudging State Arrays
  !-----------------------
  integer Nudge_nlon,Nudge_nlat,Nudge_ncol,Nudge_nlev
  real(r8),allocatable::Target_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Model_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Running_mean_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Utau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Vtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Stau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Qtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_PStau (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Nudge_PSstep(:,:)    !(pcols,begchunk:endchunk)

  ! Nudging Observation Arrays
  !-----------------------------
  logical :: Nudge_File_Present

contains
  !================================================================
  subroutine nudging_readnl(nlfile)
   ! 
   ! NUDGING_READNL: Initialize default values controlling the Nudging 
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

   namelist /nudging_nl/ Nudge_Model,Nudge_Path,                       &
                         Nudge_File_Template,Nudge_Force_Opt,          &
                         Running_mean_Path,Running_mean_File_Template, &
                         Nudge_TimeScale_Opt,                          &
                         Nudge_Times_Per_Day,Model_Times_Per_Day,      &
                         Nudge_Ucoef ,Nudge_Uprof,                     &
                         Nudge_Vcoef ,Nudge_Vprof,                     &
                         Nudge_Qcoef ,Nudge_Qprof,                     &
                         Nudge_Tcoef ,Nudge_Tprof,                     &
                         Nudge_PScoef,Nudge_PSprof,                    &
                         Nudge_Beg_Year,Nudge_Beg_Month,Nudge_Beg_Day, &
                         Nudge_End_Year,Nudge_End_Month,Nudge_End_Day, &
                         Nudge_Hwin_lat0,Nudge_Hwin_lon0,              &
                         Nudge_Hwin_latWidth,Nudge_Hwin_lonWidth,      &
                         Nudge_Hwin_latDelta,Nudge_Hwin_lonDelta,      &
                         Nudge_Hwin_Invert,                            &
                         Nudge_Vwin_Lindex,Nudge_Vwin_Hindex,          &
                         Nudge_Vwin_Ldelta,Nudge_Vwin_Hdelta,          &
                         Nudge_Vwin_Invert,                            &
                         win_size, use_gaussian 

   ! Nudging is NOT initialized yet, For now
   ! Nudging will always begin/end at midnight.
   !--------------------------------------------
   Nudge_Initialized =.false.
   Nudge_ON          =.false.
   Nudge_Beg_Sec=0
   Nudge_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Nudge_Model         = .false.
   Nudge_Path          = '/n/holylfs04/LABS/kuang_lab/Lab/sweidman/MERRA2_OG/MERRA2_f19/'
   Nudge_File_Template = 'MERRA2_%m%d_%h.nc'
   Nudge_Force_Opt     = 0
   Nudge_Path          = '/n/home04/sweidman/holylfs04/IC_CESM2/'
   Nudge_File_Template = 'cam_running_mean.%m-%d-%s.nc'
   Nudge_TimeScale_Opt = 0
   Nudge_TSmode        = 0
   Nudge_Times_Per_Day = 4
   Model_Times_Per_Day = 4
   Nudge_Ucoef         = 1._r8
   Nudge_Vcoef         = 1._r8
   Nudge_Qcoef         = 1._r8
   Nudge_Tcoef         = 1._r8
   Nudge_PScoef        = 0._r8
   Nudge_Uprof         = 1
   Nudge_Vprof         = 1
   Nudge_Qprof         = 1
   Nudge_Tprof         = 1
   Nudge_PSprof        = 0
   Nudge_Beg_Year      = 1980
   Nudge_Beg_Month     = 1
   Nudge_Beg_Day       = 1
   Nudge_End_Year      = 2019
   Nudge_End_Month     = 12
   Nudge_End_Day       = 31
   Nudge_Hwin_lat0     = 0._r8
   Nudge_Hwin_latWidth = 9999._r8
   Nudge_Hwin_latDelta = 1.0_r8
   Nudge_Hwin_lon0     = 180._r8
   Nudge_Hwin_lonWidth = 9999._r8
   Nudge_Hwin_lonDelta = 1.0_r8
   Nudge_Hwin_Invert   = .false.
   Nudge_Hwin_lo       = 0.0_r8
   Nudge_Hwin_hi       = 1.0_r8
   Nudge_Vwin_Hindex   = float(pver+1)
   Nudge_Vwin_Hdelta   = 0.001_r8
   Nudge_Vwin_Lindex   = 0.0_r8
   Nudge_Vwin_Ldelta   = 0.001_r8
   Nudge_Vwin_Invert   = .false.
   Nudge_Vwin_lo       = 0.0_r8
   Nudge_Vwin_hi       = 1.0_r8
   win_size            = 11
   use_gaussian        = .true. 


   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'nudging_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,nudging_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('nudging_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Nudge_Hwin_Invert) then
     Nudge_Hwin_lo = 1.0_r8
     Nudge_Hwin_hi = 0.0_r8
   else
     Nudge_Hwin_lo = 0.0_r8
     Nudge_Hwin_hi = 1.0_r8
   endif

   if(Nudge_Vwin_Invert) then
     Nudge_Vwin_lo = 1.0_r8
     Nudge_Vwin_hi = 0.0_r8
   else
     Nudge_Vwin_lo = 0.0_r8
     Nudge_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values 
   !----------------------------------
   if((Nudge_Hwin_lat0.lt.-90._r8).or.(Nudge_Hwin_lat0.gt.+90._r8)) then
     write(iulog,*) 'NUDGING: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lat0=',Nudge_Hwin_lat0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_lon0.lt.0._r8).or.(Nudge_Hwin_lon0.ge.360._r8)) then
     write(iulog,*) 'NUDGING: Window lon0 must be in [0,+360)'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lon0=',Nudge_Hwin_lon0
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Vwin_Lindex.gt.Nudge_Vwin_Hindex)                         .or. &
      (Nudge_Vwin_Hindex.gt.float(pver+1)).or.(Nudge_Vwin_Hindex.lt.0._r8).or. &
      (Nudge_Vwin_Lindex.gt.float(pver+1)).or.(Nudge_Vwin_Lindex.lt.0._r8)   ) then
     write(iulog,*) 'NUDGING: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'NUDGING: Lindex must be LE than Hindex'
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Lindex=',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hindex=',Nudge_Vwin_Hindex
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   if((Nudge_Hwin_latDelta.le.0._r8).or.(Nudge_Hwin_lonDelta.le.0._r8).or. &
      (Nudge_Vwin_Hdelta  .le.0._r8).or.(Nudge_Vwin_Ldelta  .le.0._r8)    ) then
     write(iulog,*) 'NUDGING: Window Deltas must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latDelta=',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonDelta=',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Hdelta=',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING:  Nudge_Vwin_Ldelta=',Nudge_Vwin_Ldelta
     call endrun('nudging_readnl:: ERROR in namelist')

   endif

   if((Nudge_Hwin_latWidth.le.0._r8).or.(Nudge_Hwin_lonWidth.le.0._r8)) then
     write(iulog,*) 'NUDGING: Window widths must be positive'
     write(iulog,*) 'NUDGING:  Nudge_Hwin_latWidth=',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING:  Nudge_Hwin_lonWidth=',Nudge_Hwin_lonWidth
     call endrun('nudging_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Nudge_Path         ,len(Nudge_Path)         ,mpichar,0,mpicom)
   call mpibcast(Nudge_File_Template,len(Nudge_File_Template),mpichar,0,mpicom)
   call mpibcast(Running_mean_Path         ,len(Running_mean_Path)         ,mpichar,0,mpicom)
   call mpibcast(Running_mean_File_Template,len(Running_mean_File_Template),mpichar,0,mpicom)
   call mpibcast(Nudge_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Force_Opt    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_TimeScale_Opt, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_TSmode       , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Model_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Ucoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Tcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_PScoef       , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Uprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Vprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Tprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Qprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_PSprof       , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Beg_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_End_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Vwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(win_size,   1, mpilog, 0, mpicom)
   call mpibcast(use_gaussian,   1, mpilog, 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_readnl
  !================================================================


  !================================================================
  subroutine nudging_init
   ! 
   ! NUDGING_INIT: Allocate space and initialize Nudging values
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
   integer               nn

   ! Get the time step size
   !------------------------
   dtime = get_step_size()

   ! Allocate Space for Nudging data arrays
   !-----------------------------------------
   allocate(Target_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Target_PS',pcols*((endchunk-begchunk)+1))

   allocate(Model_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Model_PS',pcols*((endchunk-begchunk)+1))

   allocate(Running_mean_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Running_mean_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Running_mean_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Running_mean_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Running_mean_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Running_mean_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Running_mean_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Running_mean_PS',pcols*((endchunk-begchunk)+1))

   ! Allocate Space for spatial dependence of 
   ! Nudging Coefs and Nudging Forcing.
   !-------------------------------------------
   allocate(Nudge_Utau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Utau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Vtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Vtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Stau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Stau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Qtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_PStau(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_PStau',pcols*((endchunk-begchunk)+1))

   allocate(Nudge_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Ustep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Vstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Sstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_Qstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Nudge_PSstep(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'nudging_init','Nudge_PSstep',pcols*((endchunk-begchunk)+1))


   ! Register output fields with the cam history module
   !-----------------------------------------------------
   call addfld( 'Nudge_U',(/ 'lev' /),'A','m/s/s'  ,'U Nudging Tendency')
   call addfld( 'Nudge_V',(/ 'lev' /),'A','m/s/s'  ,'V Nudging Tendency')
   call addfld( 'Nudge_T',(/ 'lev' /),'A','K/s'    ,'T Nudging Tendency')
   call addfld( 'Nudge_Q',(/ 'lev' /),'A','kg/kg/s','Q Nudging Tendency')
   call addfld('Target_U',(/ 'lev' /),'A','m/s'    ,'U Nudging Target'  )
   call addfld('Target_V',(/ 'lev' /),'A','m/s'    ,'V Nudging Target'  )
   call addfld('Target_T',(/ 'lev' /),'A','K'      ,'T Nudging Target'  )
   call addfld('Target_Q',(/ 'lev' /),'A','kg/kg'  ,'Q Nudging Target  ')

   !-----------------------------------------
   ! Values initialized only by masterproc
   !-----------------------------------------
   if(masterproc) then

     ! Set the Stepping intervals for Model and Nudging values
     ! Ensure that the Model_Step is not smaller then one timestep
     !  and not larger then the Nudge_Step.
     !--------------------------------------------------------
     Model_Step=86400/Model_Times_Per_Day
     Nudge_Step=86400/Nudge_Times_Per_Day
     if(Model_Step.lt.dtime) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: Model_Step cannot be less than a model timestep'
       write(iulog,*) 'NUDGING:  Setting Model_Step=dtime , dtime=',dtime
       write(iulog,*) ' '
       Model_Step=dtime
     endif
     if(Model_Step.gt.Nudge_Step) then
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: Model_Step cannot be more than Nudge_Step'
       write(iulog,*) 'NUDGING:  Setting Model_Step=Nudge_Step, Nudge_Step=',Nudge_Step
       write(iulog,*) ' '
       Model_Step=Nudge_Step
     endif

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Nudge_nlon=hdim1_d
     Nudge_nlat=hdim2_d
     Nudge_ncol=hdim1_d*hdim2_d
     Nudge_nlev=pver

     ! Check the time relative to the nudging window
     !------------------------------------------------
     call get_curr_date(Year,Month,Day,Sec)
     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day
     call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day
     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Nudge_End_Sec,Before_End)
  
     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to 
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Year
       Model_Next_Month=Month
       Model_Next_Day  =Day
       Model_Next_Sec  =(Sec/Model_Step)*Model_Step
       Nudge_Next_Year =Year
       Nudge_Next_Month=Month
       Nudge_Next_Day  =Day
       Nudge_Next_Sec  =(Sec/Nudge_Step)*Nudge_Step
     elseif(.not.After_Beg) then
       ! Set Time indicies to Nudging start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Nudge_Beg_Year
       Model_Next_Month=Nudge_Beg_Month
       Model_Next_Day  =Nudge_Beg_Day
       Model_Next_Sec  =Nudge_Beg_Sec
       Nudge_Next_Year =Nudge_Beg_Year
       Nudge_Next_Month=Nudge_Beg_Month
       Nudge_Next_Day  =Nudge_Beg_Day
       Nudge_Next_Sec  =Nudge_Beg_Sec
     elseif(.not.Before_End) then
       ! Nudging will never occur, so switch it off
       !--------------------------------------------
       Nudge_Model=.false.
       Nudge_ON   =.false.
       write(iulog,*) ' '
       write(iulog,*) 'NUDGING: WARNING - Nudging has been requested by it will'
       write(iulog,*) 'NUDGING:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function  
     !----------------------------------------
     lonp= 180._r8
     lon0=   0._r8
     lonn=-180._r8
     latp=  90._r8-Nudge_Hwin_lat0
     lat0=   0._r8
     latn= -90._r8-Nudge_Hwin_lat0
    
     Nudge_Hwin_lonWidthH=Nudge_Hwin_lonWidth/2._r8
     Nudge_Hwin_latWidthH=Nudge_Hwin_latWidth/2._r8

     Val1_p=(1._r8+tanh((Nudge_Hwin_lonWidthH+lonp)/Nudge_Hwin_lonDelta))/2._r8
     Val2_p=(1._r8+tanh((Nudge_Hwin_lonWidthH-lonp)/Nudge_Hwin_lonDelta))/2._r8
     Val3_p=(1._r8+tanh((Nudge_Hwin_latWidthH+latp)/Nudge_Hwin_latDelta))/2._r8
     Val4_p=(1._r8+tanh((Nudge_Hwin_latWidthH-latp)/Nudge_Hwin_latDelta))/2_r8
     Val1_0=(1._r8+tanh((Nudge_Hwin_lonWidthH+lon0)/Nudge_Hwin_lonDelta))/2._r8
     Val2_0=(1._r8+tanh((Nudge_Hwin_lonWidthH-lon0)/Nudge_Hwin_lonDelta))/2._r8
     Val3_0=(1._r8+tanh((Nudge_Hwin_latWidthH+lat0)/Nudge_Hwin_latDelta))/2._r8
     Val4_0=(1._r8+tanh((Nudge_Hwin_latWidthH-lat0)/Nudge_Hwin_latDelta))/2._r8

     Val1_n=(1._r8+tanh((Nudge_Hwin_lonWidthH+lonn)/Nudge_Hwin_lonDelta))/2._r8
     Val2_n=(1._r8+tanh((Nudge_Hwin_lonWidthH-lonn)/Nudge_Hwin_lonDelta))/2._r8
     Val3_n=(1._r8+tanh((Nudge_Hwin_latWidthH+latn)/Nudge_Hwin_latDelta))/2._r8
     Val4_n=(1._r8+tanh((Nudge_Hwin_latWidthH-latn)/Nudge_Hwin_latDelta))/2._r8

     Nudge_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Nudge_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     Nudge_File_Present=.false.

     ! Initialization is done, 
     !--------------------------
     Nudge_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if((.not.dycore_is('UNSTRUCTURED')).and. &
        (.not.dycore_is('EUL')         ).and. &
        (.not.dycore_is('LR')          )      ) then
       call endrun('NUDGING IS CURRENTLY ONLY CONFIGURED FOR CAM-SE, FV, or EUL')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL NUDGING INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) 'NUDGING: Nudge_Model=',Nudge_Model
     write(iulog,*) 'NUDGING: Nudge_Path=',Nudge_Path
     write(iulog,*) 'NUDGING: Nudge_File_Template =',Nudge_File_Template
     write(iulog,*) 'NUDGING: Nudge_Force_Opt=',Nudge_Force_Opt    
     write(iulog,*) 'NUDGING: Nudge_TimeScale_Opt=',Nudge_TimeScale_Opt    
     write(iulog,*) 'NUDGING: Nudge_TSmode=',Nudge_TSmode
     write(iulog,*) 'NUDGING: Nudge_Times_Per_Day=',Nudge_Times_Per_Day
     write(iulog,*) 'NUDGING: Model_Times_Per_Day=',Model_Times_Per_Day
     write(iulog,*) 'NUDGING: Nudge_Step=',Nudge_Step
     write(iulog,*) 'NUDGING: Model_Step=',Model_Step
     write(iulog,*) 'NUDGING: Nudge_Ucoef  =',Nudge_Ucoef
     write(iulog,*) 'NUDGING: Nudge_Vcoef  =',Nudge_Vcoef
     write(iulog,*) 'NUDGING: Nudge_Qcoef  =',Nudge_Qcoef
     write(iulog,*) 'NUDGING: Nudge_Tcoef  =',Nudge_Tcoef
     write(iulog,*) 'NUDGING: Nudge_PScoef =',Nudge_PScoef
     write(iulog,*) 'NUDGING: Nudge_Uprof  =',Nudge_Uprof
     write(iulog,*) 'NUDGING: Nudge_Vprof  =',Nudge_Vprof
     write(iulog,*) 'NUDGING: Nudge_Qprof  =',Nudge_Qprof
     write(iulog,*) 'NUDGING: Nudge_Tprof  =',Nudge_Tprof
     write(iulog,*) 'NUDGING: Nudge_PSprof =',Nudge_PSprof
     write(iulog,*) 'NUDGING: Nudge_Beg_Year =',Nudge_Beg_Year
     write(iulog,*) 'NUDGING: Nudge_Beg_Month=',Nudge_Beg_Month
     write(iulog,*) 'NUDGING: Nudge_Beg_Day  =',Nudge_Beg_Day
     write(iulog,*) 'NUDGING: Nudge_End_Year =',Nudge_End_Year
     write(iulog,*) 'NUDGING: Nudge_End_Month=',Nudge_End_Month
     write(iulog,*) 'NUDGING: Nudge_End_Day  =',Nudge_End_Day
     write(iulog,*) 'NUDGING: Nudge_Hwin_lat0     =',Nudge_Hwin_lat0
     write(iulog,*) 'NUDGING: Nudge_Hwin_latWidth =',Nudge_Hwin_latWidth
     write(iulog,*) 'NUDGING: Nudge_Hwin_latDelta =',Nudge_Hwin_latDelta
     write(iulog,*) 'NUDGING: Nudge_Hwin_lon0     =',Nudge_Hwin_lon0
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidth =',Nudge_Hwin_lonWidth
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonDelta =',Nudge_Hwin_lonDelta
     write(iulog,*) 'NUDGING: Nudge_Hwin_Invert   =',Nudge_Hwin_Invert  
     write(iulog,*) 'NUDGING: Nudge_Hwin_lo       =',Nudge_Hwin_lo
     write(iulog,*) 'NUDGING: Nudge_Hwin_hi       =',Nudge_Hwin_hi
     write(iulog,*) 'NUDGING: Nudge_Vwin_Hindex   =',Nudge_Vwin_Hindex
     write(iulog,*) 'NUDGING: Nudge_Vwin_Hdelta   =',Nudge_Vwin_Hdelta
     write(iulog,*) 'NUDGING: Nudge_Vwin_Lindex   =',Nudge_Vwin_Lindex
     write(iulog,*) 'NUDGING: Nudge_Vwin_Ldelta   =',Nudge_Vwin_Ldelta
     write(iulog,*) 'NUDGING: Nudge_Vwin_Invert   =',Nudge_Vwin_Invert  
     write(iulog,*) 'NUDGING: Nudge_Vwin_lo       =',Nudge_Vwin_lo
     write(iulog,*) 'NUDGING: Nudge_Vwin_hi       =',Nudge_Vwin_hi
     write(iulog,*) 'NUDGING: Nudge_Hwin_latWidthH=',Nudge_Hwin_latWidthH
     write(iulog,*) 'NUDGING: Nudge_Hwin_lonWidthH=',Nudge_Hwin_lonWidthH
     write(iulog,*) 'NUDGING: Nudge_Hwin_max      =',Nudge_Hwin_max
     write(iulog,*) 'NUDGING: Nudge_Hwin_min      =',Nudge_Hwin_min
     write(iulog,*) 'NUDGING: Nudge_Initialized   =',Nudge_Initialized
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Model_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Model         ,            1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ON            ,            1, mpilog, 0, mpicom)
   call mpibcast(Nudge_Initialized   ,            1, mpilog, 0, mpicom)
   call mpibcast(Nudge_ncol          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlev          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlon          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_nlat          ,            1, mpiint, 0, mpicom)
   call mpibcast(Nudge_Hwin_max      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_min      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_lonWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Nudge_Hwin_latWidthH,            1, mpir8 , 0, mpicom)
#endif


!!DIAG
   if(masterproc) then
     write(iulog,*) 'NUDGING: nudging_init() OBS arrays allocated and initialized'
     write(iulog,*) 'NUDGING: nudging_init() SIZE#',(9*pcols*pver*((endchunk-begchunk)+1))
     write(iulog,*) 'NUDGING: nudging_init() MB:',float(8*9*pcols*pver*((endchunk-begchunk)+1))/(1024._r8*1024._r8)
     write(iulog,*) 'NUDGING: nudging_init() pcols=',pcols,' pver=',pver
     write(iulog,*) 'NUDGING: nudging_init() begchunk:',begchunk,' endchunk=',endchunk
     write(iulog,*) 'NUDGING: nudging_init() Nudge_File_Present=',Nudge_File_Present
   endif
!!DIAG

   ! TODO: initialize the running mean file: should it be all 0s or something else?

   ! Initialize the analysis filename at the NEXT time for startup.
   ! TODO: this is not the right file to read -- use the previous timestep
   !---------------------------------------------------------------
   Nudge_File=interpret_filename_spec(Nudge_File_Template      , &
                                       yr_spec=Nudge_Next_Year , &
                                      mon_spec=Nudge_Next_Month, &
                                      day_spec=Nudge_Next_Day  , &
                                      sec_spec=Nudge_Next_Sec    )
   if(masterproc) then
    write(iulog,*) 'NUDGING: Reading analyses:',trim(Nudge_Path)//trim(Nudge_File)
   endif

   !----------------------------------------------------------

    call nudging_update_analyses_fv (trim(Nudge_Path)//trim(Nudge_File))

   ! Initialize Nudging Coeffcient profiles in local arrays
   ! Load zeros into nudging arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call nudging_set_profile(rlat,rlon,Nudge_Uprof,Wprof,pver)
       Nudge_Utau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Vprof,Wprof,pver)
       Nudge_Vtau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Tprof,Wprof,pver)
       Nudge_Stau(icol,:,lchnk)=Wprof(:)
       call nudging_set_profile(rlat,rlon,Nudge_Qprof,Wprof,pver)
       Nudge_Qtau(icol,:,lchnk)=Wprof(:)

       Nudge_PStau(icol,lchnk)=nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
     end do
     Nudge_Utau(:ncol,:pver,lchnk) =                             &
     Nudge_Utau(:ncol,:pver,lchnk) * Nudge_Ucoef/float(Nudge_Step)
     Nudge_Vtau(:ncol,:pver,lchnk) =                             &
     Nudge_Vtau(:ncol,:pver,lchnk) * Nudge_Vcoef/float(Nudge_Step)
     Nudge_Stau(:ncol,:pver,lchnk) =                             &
     Nudge_Stau(:ncol,:pver,lchnk) * Nudge_Tcoef/float(Nudge_Step)
     Nudge_Qtau(:ncol,:pver,lchnk) =                             &
     Nudge_Qtau(:ncol,:pver,lchnk) * Nudge_Qcoef/float(Nudge_Step)
     Nudge_PStau(:ncol,lchnk)=                             &
     Nudge_PStau(:ncol,lchnk)* Nudge_PScoef/float(Nudge_Step)

     Nudge_Ustep(:pcols,:pver,lchnk)=0._r8
     Nudge_Vstep(:pcols,:pver,lchnk)=0._r8
     Nudge_Sstep(:pcols,:pver,lchnk)=0._r8
     Nudge_Qstep(:pcols,:pver,lchnk)=0._r8
     Nudge_PSstep(:pcols,lchnk)=0._r8
     Target_U(:pcols,:pver,lchnk)=0._r8
     Target_V(:pcols,:pver,lchnk)=0._r8
     Target_T(:pcols,:pver,lchnk)=0._r8
     Target_S(:pcols,:pver,lchnk)=0._r8
     Target_Q(:pcols,:pver,lchnk)=0._r8
     Target_PS(:pcols,lchnk)=0._r8
   end do

   ! End Routine
   !------------
   return
  end subroutine ! nudging_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_init(phys_state)
   ! 
   ! NUDGING_TIMESTEP_INIT: 
   !                 Check the current time and update Model/Nudging 
   !                 arrays when necessary. Toggle the Nudging flag
   !                 when the time is withing the nudging window.
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
   logical Update_Model,Update_Nudge,Sync_Error
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
   integer                 dtime
   integer                 modstep

   real(r8)                wrk
   ! Check if Nudging is initialized
   !---------------------------------
   if(.not.Nudge_Initialized) then
     call endrun('nudging_timestep_init:: Nudging NOT Initialized')
   endif

   ! Get time step size
   !--------------------
   dtime = get_step_size()

   ! TODO: set wrk to a reasonable value, change over time (I think)
   if (nstep_count > 96000) then
       wrk = 1._r8/96000._r8
   else 
       wrk = 1._r8 / (1._r8 + nstep_count-1)
   endif

   if (masterproc) then
     write(iulog,*) "wrk XXXXXXXXXXXXX", wrk
   end if

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(Nudge_Beg_Year*10000) + (Nudge_Beg_Month*100) + Nudge_Beg_Day
   call timemgr_time_ge(YMD1,Nudge_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Nudge_End_Year*10000) + (Nudge_End_Month*100) + Nudge_End_Day
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Nudge_End_Sec,Before_End)

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
       write(iulog,*) 'NUDGING: WARNING - Model_Time Sync ERROR... CORRECTED'
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
       Model_PS(:ncol,lchnk)=phys_state(lchnk)%ps(:ncol)
     end do

     ! Load Dry Static Energy values for Model
     !-----------------------------------------
     if(Nudge_TSmode.eq.0) then
       ! DSE tendencies from Temperature only
       !---------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Model_S(:ncol,:pver,lchnk)=cpair*Model_T(:ncol,:pver,lchnk)
       end do
     elseif(Nudge_TSmode.eq.1) then
       ! Caluculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
       !------------------------------------------------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         call calc_DryStaticEnergy(Model_T(:,:,lchnk)  , Model_Q(:,:,lchnk), &
                                 phys_state(lchnk)%phis,  Model_PS(:,lchnk), &
                                                  Model_S(:,:,lchnk), ncol)
       end do
     endif 

     ! TODO: This is where we need to load in the running mean from current timestep
      Running_mean_File=interpret_filename_spec(Running_mean_File_Template      , &
                                        sec_spec=Model_Curr_Sec    )
      INQUIRE(FILE=trim(Running_mean_Path)//trim(Running_mean_File), EXIST=fileexists)
      
      if (.not. fileexists) print*, 'running mean file missing', Running_mean_File
     
     if(masterproc) then
      write(iulog,*) 'NUDGING: Reading analyses:',trim(Running_mean_Path)//trim(Running_mean_File)
     endif

     !----------------------------------------------------------
    call nudging_update_model_fv (trim(Running_mean_Path)//trim(Running_mean_File), Model_Curr_Month, Nudge_Curr_Day)

     ! Load Dry Static Energy values for Model
     !-----------------------------------------
     if(Nudge_TSmode.eq.0) then
       ! DSE tendencies from Temperature only
       !---------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Running_mean_S(:ncol,:pver,lchnk)=cpair*Running_mean_T(:ncol,:pver,lchnk)
       end do
     elseif(Nudge_TSmode.eq.1) then
       ! Caluculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
       !------------------------------------------------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         call calc_DryStaticEnergy(Running_mean_T(:,:,lchnk)  , Running_mean_Q(:,:,lchnk), &
                                 phys_state(lchnk)%phis,  Running_mean_PS(:,lchnk), &
                                                  Running_mean_S(:,:,lchnk), ncol)
       end do
     endif 

    do lchnk=begchunk,endchunk
        ! TODO: this is where I update the running mean 
         ncol=phys_state(lchnk)%ncol
         Running_mean_U(:ncol,:pver,lchnk)=Running_mean_U(:ncol,:pver,lchnk)*wrk + Model_U(:ncol,:pver,lchnk)*(1-wrk)
         Running_mean_V(:ncol,:pver,lchnk)=Running_mean_V(:ncol,:pver,lchnk)*wrk + Model_V(:ncol,:pver,lchnk)*(1-wrk)
         Running_mean_S(:ncol,:pver,lchnk)=Running_mean_S(:ncol,:pver,lchnk)*wrk + Model_S(:ncol,:pver,lchnk)*(1-wrk)
         Running_mean_Q(:ncol,:pver,lchnk)=Running_mean_Q(:ncol,:pver,lchnk)*wrk + Model_Q(:ncol,:pver,lchnk)*(1-wrk)
         Running_mean_PS(:ncol     ,lchnk)=Running_mean_PS(:ncol     ,lchnk)*wrk + Model_PS(:ncol     ,lchnk)*(1-wrk)
       end do

   endif ! ((Before_End).and.(Update_Model)) then



   !----------------------------------------------------------------
   ! When past the NEXT time, Update Nudging Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Nudge_Next_Year*10000) + (Nudge_Next_Month*100) + Nudge_Next_Day
   call timemgr_time_ge(YMD1,Nudge_Next_Sec,            &
                        YMD ,Sec           ,Update_Nudge)

   if((Before_End).and.(Update_Nudge)) then
     ! Increment the Nudge times by the current interval
     !---------------------------------------------------
     Nudge_Curr_Year =Nudge_Next_Year
     Nudge_Curr_Month=Nudge_Next_Month
     Nudge_Curr_Day  =Nudge_Next_Day
     Nudge_Curr_Sec  =Nudge_Next_Sec
     YMD1=(Nudge_Curr_Year*10000) + (Nudge_Curr_Month*100) + Nudge_Curr_Day
     call timemgr_time_inc(YMD1,Nudge_Curr_Sec,              &
                           YMD2,Nudge_Next_Sec,Nudge_Step,0,0)
     Nudge_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Nudge_Next_Year*10000)
     Nudge_Next_Month=(YMD2/100)
     Nudge_Next_Day  = YMD2-(Nudge_Next_Month*100)

     ! Set the analysis filename at the NEXT time. (MERRA)
     !---------------------------------------------------------------
     modstep=int(Nudge_Next_Sec / 10800)
     filename=interpret_filename_climo(Nudge_File_Template      , &
          mon_spec=Nudge_Next_Month, &
          day_spec=Nudge_Next_Day  , &
          hr_spec=modstep, &
          sec_spec=Nudge_Next_Sec    )

      if(masterproc) then
        write(iulog,*) trim(Nudge_Path)//trim(filename)
      endif
      
      INQUIRE(FILE=trim(Nudge_Path)//trim(filename), EXIST=fileexists)
      if (.not. fileexists) print*, 'nudging file missing', filename

     !----------------------------------------------------------
    call nudging_update_analyses_fv (trim(Nudge_Path)//trim(Nudge_File))
   endif ! ((Before_End).and.(Update_Nudge)) then

   !----------------------------------------------------------------
   ! Toggle Nudging flag when the time interval is between 
   ! beginning and ending times, and all of the analyses files exist.
   !----------------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
       Nudge_ON=Nudge_File_Present
   else
     Nudge_ON=.false.
   endif

   !-------------------------------------------------------
   ! HERE Implement time dependence of Nudging Coefs HERE
   !-------------------------------------------------------
   ! TODO: calculate nudging coefficient

   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if((Before_End).and.((Update_Nudge).or.(Update_Model))) then

     ! Now Load the Target values for nudging tendencies
     !---------------------------------------------------
     ! Target is OBS data at NEXT time
     !----------------------------------
     ! Now load Dry Static Energy values for Target
     !---------------------------------------------
     if(Nudge_TSmode.eq.0) then
       ! DSE tendencies from Temperature only
       !---------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_S(:ncol,:pver,lchnk)=cpair*Target_T(:ncol,:pver,lchnk)
       end do
     elseif(Nudge_TSmode.eq.1) then
       ! Caluculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
       !------------------------------------------------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         call calc_DryStaticEnergy(Target_T(:,:,lchnk), Target_Q(:,:,lchnk), &
                                 phys_state(lchnk)%phis, Target_PS(:,lchnk), &
                                                  Target_S(:,:,lchnk), ncol)
       end do
     endif

     ! Set Tscale for the specified Forcing Option 
     !-----------------------------------------------
     if(Nudge_TimeScale_Opt.eq.0) then
       Tscale=1._r8
     elseif(Nudge_TimeScale_Opt.eq.1) then
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Nudge_Next_Year,MM=Nudge_Next_Month, &
                               DD=Nudge_Next_Day , S=Nudge_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tscale=float(Nudge_Step)/float(DeltaT)
     else
       write(iulog,*) 'NUDGING: Unknown Nudge_TimeScale_Opt=',Nudge_TimeScale_Opt
       call endrun('nudging_timestep_init:: ERROR unknown Nudging_TimeScale_Opt')
     endif

     ! Update the nudging tendencies
     !--------------------------------
     ! TODO: change to update Model with running mean and current timestep state
     ! ptend = ptend + q_nudge*tau
     ! save q_nudge also to file again, updating nearest +/- 5 day neighbors (?)
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Nudge_Ustep(:ncol,:pver,lchnk)=(  Target_U(:ncol,:pver,lchnk)      &
                                         -Running_mean_U(:ncol,:pver,lchnk))     &
                                      *Tscale*Nudge_Utau(:ncol,:pver,lchnk)
       Nudge_Vstep(:ncol,:pver,lchnk)=(  Target_V(:ncol,:pver,lchnk)      &
                                         -Running_mean_V(:ncol,:pver,lchnk))     &
                                      *Tscale*Nudge_Vtau(:ncol,:pver,lchnk)
       Nudge_Sstep(:ncol,:pver,lchnk)=(  Target_S(:ncol,:pver,lchnk)      &
                                         -Running_mean_S(:ncol,:pver,lchnk))     &
                                      *Tscale*Nudge_Stau(:ncol,:pver,lchnk)
       Nudge_Qstep(:ncol,:pver,lchnk)=(  Target_Q(:ncol,:pver,lchnk)      &
                                         -Running_mean_Q(:ncol,:pver,lchnk))     &
                                      *Tscale*Nudge_Qtau(:ncol,:pver,lchnk)
       Nudge_PSstep(:ncol,     lchnk)=(  Target_PS(:ncol,lchnk)      &
                                         -Running_mean_PS(:ncol,lchnk))     &
                                      *Tscale*Nudge_PStau(:ncol,lchnk)
     end do

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
!      write(iulog,*) 'PFC: Nudge_Sstep(1,:pver,begchunk)=',Nudge_Sstep(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Nudge_Xstep arrays updated:'
!    endif
   endif ! ((Before_End).and.((Update_Nudge).or.(Update_Model))) then

   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_init
  !================================================================


  !================================================================
  subroutine nudging_timestep_tend(phys_state,phys_tend)
   ! 
   ! NUDGING_TIMESTEP_TEND: 
   !                If Nudging is ON, return the Nudging contributions 
   !                to forcing using the current contents of the Nudge 
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
   call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Nudge_ON) then
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
     phys_tend%u(:ncol,:pver)     =Nudge_Ustep(:ncol,:pver,lchnk)
     phys_tend%v(:ncol,:pver)     =Nudge_Vstep(:ncol,:pver,lchnk)
     phys_tend%s(:ncol,:pver)     =Nudge_Sstep(:ncol,:pver,lchnk)
     phys_tend%q(:ncol,:pver,indw)=Nudge_Qstep(:ncol,:pver,lchnk)

     call outfld( 'Nudge_U',phys_tend%u                ,pcols,lchnk)
     call outfld( 'Nudge_V',phys_tend%v                ,pcols,lchnk)
     call outfld( 'Nudge_T',phys_tend%s/cpair          ,pcols,lchnk)
     call outfld( 'Nudge_Q',phys_tend%q(1,1,indw)      ,pcols,lchnk)

     ! TODO: save Running_mean_x back into file 

   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_tend
  !================================================================


  !================================================================
  subroutine nudging_update_model_fv(running_mean_file, target_month, target_day)
   ! 
   ! NUDGING_UPDATE_ANALYSES_FV: 
   !                 Open the given analyses data file, read in 
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: running_mean_file

   ! Local values
   !-------------
   ! adding time variables ++SW
   integer lev
   integer nlon,nlat,plev,istat,ntime
   integer ncid,varid,varid_t
   integer ilat,ilon,ilev,it, iw
   integer Nudge_ntime ! unsure if I need to set this value earlier ++SW
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
   integer  nn,Nindex

   integer,intent(in):: target_month, target_day  ! adding for centered mean ++SW

   ! adding for taking weighted mean ++SW
   real, allocatable :: Time_anal(:)                      ! time dimension (e.g., days since ref)
   real, allocatable :: Uslab(:,:,:)                      ! lon x lat x lev for one time
   real, allocatable :: w(:)                              ! window weights
   integer, allocatable :: t_indices(:)                   ! actual time indices used

   ! window config (11-day window)
   integer, parameter :: win_size = 11 ! TODO: make this a namelist option
   integer, dimension(win_size) :: win_offsets 
   logical, parameter :: use_gaussian = .true.     ! TODO: make this a namelist option
   real :: sigma_days, wsum
   integer :: halfm, it_center   ! for creating centered window indices

   ! Just read in one file but will select times inside it
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(running_mean_file),EXIST=Nudge_File_Present)
     write(iulog,*)'running mean nudge: Nudge_File_Present=',Nudge_File_Present
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, 1, mpilog, 0, mpicom)
#endif
   if(.not.Nudge_File_Present) return

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

     ! added read in time dimension ++SW
     istat=nf90_inq_varid(ncid,'time',varid_t)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     
     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_fv: analyses dimension mismatch')
     endif

     ! allocate extra time variables
     allocate(Time_anal(ntime))
     allocate(Uslab(Nudge_nlon,Nudge_nlat,Nudge_nlev)) ! using Xanal dimensions (lat,lev flipped)
     allocate(w(win_size))
     allocate(t_indices(win_size))

     istat=nf90_get_var(ncid,varid_t,Time_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! TODO: added chatgpt code below here ++SW

    ! Convert requested Y-M-D to the same numeric "days since ..." as time_vals
    ! (We avoid parsing units by matching to the nearest available entry via a helper.)
    ! If you want exact conversion from units, extend the helper accordingly.
    call find_nearest_time_index(time_vals, target_year, target_month, target_day, it_center)

    ! Define a centered window
    half = (win_size - 1)/2
    win_offsets = [(i, i=-half, half)] 

    ! Map offsets to legal indices with wrap-around (use modulo year logic). 
    do iw = 1, win_size
      t_indices(iw) = modulo(it_center - 1 + win_offsets(iw), ntime) + 1  ! Fortran 1-based, modulo wrap
    end do

    ! Build weights (choose triangular or Gaussian)
    if (.not. use_gaussian) then
      ! Triangular: weight drops linearly with |offset|, peak at center. Ensure non-negative.
      do iw = 1, win_size
        w(iw) = real(half + 1 - abs(win_offsets(iw)))   ! e.g., 3,2,1,2,3 for half=2
      end do
    else
      ! sigma relative to window half-width; sigma_scale ~ 0.5 works well
      sigma_days = max(1.0e-6, sigma_scale*real(max(1,half)))
      do iw = 1, win_size
        w(iw) = exp( -0.5 * ( real(win_offsets(iw)) / sigma_days )**2 )
      end do
    end if
    wsum = sum(w);  if (wsum <= 0.0) then
      call endrun('UPDATE_ANALYSES_FV: zero/neg window weight sum')
    end if
    w = w / wsum   ! normalize

    ! Accumulate weighted mean by reading one time slab at a time
    do iw = 1, win_size
      integer :: itime
      integer, dimension(4) :: start, count
      itime = t_indices(iw)

      ! Assuming U(lon,lat,lev,time) dimension order
      start = (/ 1, 1, 1, itime /)
      count = (/ nlon, nlat, nlev, 1 /)

      istat = nf90_get_var(ncid, varid_U, Uslab, start=start, count=count)
      if (istat /= NF90_NOERR) then
        write(iulog,*) nf90_strerror(istat); call endrun('UPDATE_ANALYSES_FV(U slab read)')
      end if

      Xmean = Xmean + w(iw) * Uslab
    end do

    ! Xmean(lon,lat,lev) is the weighted climatology for the date window.
    do ilat = 1, nlat
      do ilev = 1, nlev
        do ilon = 1, nlon
          Xtrans(ilon, ilev, ilat) = Xmean(ilon, ilat, ilev)
        end do
      end do
    end do

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
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
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
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
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
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
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
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Target_Q(1,1,begchunk))


     

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_model_fv
  !================================================================

  !================================================================
  subroutine nudging_update_analyses_fv(anal_file)
   ! 
   ! NUDGING_UPDATE_ANALYSES_FV: 
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
   real(r8) Xanal(Nudge_nlon,Nudge_nlat,Nudge_nlev)
   real(r8) PSanal(Nudge_nlon,Nudge_nlat)
   real(r8) Lat_anal(Nudge_nlat)
   real(r8) Lon_anal(Nudge_nlon)
   real(r8) Xtrans(Nudge_nlon,Nudge_nlev,Nudge_nlat)
   integer  nn

   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     inquire(FILE=trim(anal_file),EXIST=Nudge_File_Present)
     write(iulog,*)'NUDGING: Nudge_File_Present=',Nudge_File_Present
   endif
#ifdef SPMD
   call mpibcast(Nudge_File_Present, 1, mpilog, 0, mpicom)

#endif
   if(.not.Nudge_File_Present) return

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

     if((Nudge_nlon.ne.nlon).or.(Nudge_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlon=',nlon,' Nudge_nlon=',Nudge_nlon
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: nlat=',nlat,' Nudge_nlat=',Nudge_nlat
      write(iulog,*) 'ERROR: nudging_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('nudging_update_analyses_fv: analyses dimension mismatch')
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
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
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
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
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
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
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
   call scatter_field_to_chunk(1,Nudge_nlev,1,Nudge_nlon,Xtrans,   &
                               Target_Q(1,1,begchunk))

   ! End Routine
   !------------
   return
  end subroutine ! nudging_update_analyses_fv
  !================================================================

  subroutine find_nearest_time_index(time_vals, mm, dd, idx)
    real,    intent(in)  :: time_vals(:)
    integer, intent(in)  :: mm, dd
    integer, dimension(12) :: cum = (/ 0,31,59,90,120,151,181,212,243,273,304,334 /)
    integer, intent(out) :: idx
    ! Strategy: if the file is daily and covers a single year, we can map to DOY and
    ! then pick that DOY (with clamp or wrap). Otherwise, pick the time entry whose
    ! calendar date matches best by comparing to the list of dates you derive from
    ! an origin. To keep this simple and robust for daily files, we do DOY.

    integer :: doy, n
    doy = cum(mm) + dd
    n = size(time_vals)

    ! probably have to shift for 3-hourly data (though maybe I should save as hourly and with separate variables)
    if (n >= 365) then
      idx = modulo(doy-1, n) + 1   !  wrap
    else
      ! Fallback: nearest by numeric distance to the median time of requested DOY
      ! (This path is rarely hit if you truly have daily data.)
      real :: target
      ! crude: assume time_vals are “days since” and Jan 1 is near integer boundary
      target = real(doy - 1)
      idx = 1 + minloc( abs(time_vals - target), dim=1 ) - 1
    end if

    ! TODO: make sure this works for 3-hourly data
  end subroutine find_nearest_time_index

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
      call endrun ('INTERPRET_FILENAME_SPEC: filename specifier is empty')
   end if
   if ( index(trim(filename_spec)," ") /= 0 )then
      call endrun ('INTERPRET_FILENAME_SPEC: filename specifier can not contain a space:'//trim(filename_spec))
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
            call endrun ('INTERPRET_FILENAME_SPEC: Invalid expansion character: '//filename_spec(i:i))
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
            call endrun ('INTERPRET_FILENAME_SPEC: Resultant filename too long')
         end if
         interpret_filename_climo = trim(interpret_filename_climo) // trim(string)
      end if
      i = i + 1

   end do
   if ( len_trim(interpret_filename_climo) == 0 )then
      call endrun ('INTERPRET_FILENAME_SPEC: Resulting filename is empty')
   end if

end function interpret_filename_climo

  !================================================================
  subroutine nudging_set_profile(rlat,rlon,Nudge_prof,Wprof,nlev)
   ! 
   ! NUDGING_SET_PROFILE: for the given lat,lon, and Nudging_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Nudge_prof
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
   if(Nudge_prof.eq.0) then
     ! No Nudging
     !-------------
     Wprof(:)=0.0_r8
   elseif(Nudge_prof.eq.1) then
     ! Uniform Nudging
     !-----------------
     Wprof(:)=1.0_r8
   elseif(Nudge_prof.eq.2) then
     ! Localized Nudging with specified Heaviside window function
     !------------------------------------------------------------
     if(Nudge_Hwin_max.le.Nudge_Hwin_min) then
       ! For a constant Horizontal window function, 
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Nudge_Hwin_lo,Nudge_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Nudge_Hwin_lat0
       lonx=rlon-Nudge_Hwin_lon0
       if(lonx.gt. 180._r8) lonx=lonx-360._r8
       if(lonx.le.-180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Nudge_Hwin_lonWidthH+lonx)/Nudge_Hwin_lonDelta
       lon_hi=(Nudge_Hwin_lonWidthH-lonx)/Nudge_Hwin_lonDelta
       lat_lo=(Nudge_Hwin_latWidthH+latx)/Nudge_Hwin_latDelta
       lat_hi=(Nudge_Hwin_latWidthH-latx)/Nudge_Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Nudge_Hwin_min)/(Nudge_Hwin_max-Nudge_Hwin_min)
       Hcoef=(1._r8-Hcoef)*Nudge_Hwin_lo + Hcoef*Nudge_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Nudge_Vwin_Lindex)/Nudge_Vwin_Ldelta
       lev_hi=(Nudge_Vwin_Hindex-float(ilev))/Nudge_Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do 

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax.le.Vmin).or.((Nudge_Vwin_Hindex.ge.(nlev+1)).and. &
                           (Nudge_Vwin_Lindex.le. 0      )     )) then
       ! For a constant Vertical window function, 
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Nudge_Vwin_lo,Nudge_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Nudge_Vwin_lo + Wprof(:)*(Nudge_Vwin_hi-Nudge_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile 
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('nudging_set_profile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! nudging_set_profile
  !================================================================


  !================================================================
  real(r8) function nudging_set_PSprofile(rlat,rlon,Nudge_PSprof)
   ! 
   ! NUDGING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on nudging strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Nudge_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Nudge_PSprof.eq.0) then
     ! No Nudging
     !-------------
     nudging_set_PSprofile=0.0_r8
   elseif(Nudge_PSprof.eq.1) then
     ! Uniform Nudging
     !-----------------
     nudging_set_PSprofile=1.0_r8
   else
     call endrun('nudging_set_PSprofile:: Unknown Nudge_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! nudging_set_PSprofile
  !================================================================


  !================================================================
  subroutine calc_DryStaticEnergy(t, q, phis, ps, dse, ncol)
   ! 
   ! calc_DryStaticEnergy: Given the temperature, specific humidity, surface pressure,
   !                       and surface geopotential for a chunk containing 'ncol' columns, 
   !                       calculate and return the corresponding dry static energy values.
   !--------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pver, pverp
   use dycore,       only: dycore_is
   use hycoef,       only: hyai, hybi, ps0, hyam, hybm
   use physconst,    only: zvir, gravit, cpair, rair
   !
   ! Input/Output arguments
   !-----------------------
   integer , intent(in) :: ncol      ! Number of columns in chunk
   real(r8), intent(in) :: t(:,:)    ! (pcols,pver) - temperature
   real(r8), intent(in) :: q(:,:)    ! (pcols,pver) - specific humidity
   real(r8), intent(in) :: ps(:)     ! (pcols)      - surface pressure 
   real(r8), intent(in) :: phis(:)   ! (pcols)      - surface geopotential
   real(r8), intent(out):: dse(:,:)  ! (pcols,pver)  - dry static energy
   !
   ! Local variables
   !------------------
   logical  :: fvdyn                 ! finite volume dynamics
   integer  :: ii,kk                 ! Lon, level, level indices
   real(r8) :: tvfac                 ! Virtual temperature factor
   real(r8) :: hkk(ncol)             ! diagonal element of hydrostatic matrix
   real(r8) :: hkl(ncol)             ! off-diagonal element
   real(r8) :: pint(ncol,pverp)      ! Interface pressures
   real(r8) :: pmid(ncol,pver )      ! Midpoint pressures
   real(r8) :: zi(ncol,pverp)        ! Height above surface at interfaces
   real(r8) :: zm(ncol,pver )        ! Geopotential height at mid level

   ! Set dynamics flag
   !-------------------
   fvdyn = dycore_is ('LR')

   ! Load Pressure values and midpoint pressures 
   !----------------------------------------------
   do kk=1,pverp
     do ii=1,ncol
       pint(ii,kk)=(hyai(kk)*ps0)+(hybi(kk)*ps(ii))
     end do
   end do
   do kk=1,pver
     do ii=1,ncol
       pmid(ii,kk)=(hyam(kk)*ps0)+(hybm(kk)*ps(ii))
     end do
   end do

   ! The surface height is zero by definition.
   !-------------------------------------------
   do ii = 1,ncol
     zi(ii,pverp) = 0.0_r8
   end do

   ! Compute the dry static energy, zi, zm from bottom up
   ! Note, zi(i,k) is the interface above zm(i,k)
   !---------------------------------------------------------
   do kk=pver,1,-1

     ! First set hydrostatic elements consistent with dynamics
     !--------------------------------------------------------
     if(fvdyn) then
       do ii=1,ncol
         hkl(ii)=log(pint(ii,kk+1))-log(pint(ii,kk))
         hkk(ii)=1._r8-(hkl(ii)*pint(ii,kk)/(pint(ii,kk+1)-pint(ii,kk)))
       end do
     else
       do ii=1,ncol
         hkl(ii)=(pint(ii,kk+1)-pint(ii,kk))/pmid(ii,kk)
         hkk(ii)=0.5_r8*hkl(ii)
       end do
     endif

     ! Now compute zm, zi, and dse  (WACCM-X vars rairv/zairv/cpairv not used!)
     !------------------------------------------------------------------------
     do ii=1,ncol
       tvfac=t(ii,kk)*rair*(1._r8+(zvir*q(ii,kk)))/gravit
       zm (ii,kk)=zi(ii,kk+1) + (tvfac*hkk(ii))
       zi (ii,kk)=zi(ii,kk+1) + (tvfac*hkl(ii))
       dse(ii,kk)=(t(ii,kk)*cpair) + phis(ii) + (gravit*zm(ii,kk))
     end do

   end do ! kk=pver,1,-1

   ! End Routine
   !-----------
   return
  end subroutine calc_DryStaticEnergy
  !================================================================

end module nudging
