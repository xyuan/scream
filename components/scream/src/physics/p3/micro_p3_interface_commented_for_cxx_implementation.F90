module micro_p3_interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Interface between E3SM and P3 microphysics
  !!
  !! Author: Peter Caldwell
  !!
  !! Last updated: 2018-09-12
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use shr_kind_mod,   only: rtype=>shr_kind_r8

!comment: I think Kai added handle_errmsg. It would be better to 
!use standard E3SM libraries if possible.

  use micro_p3_utils, only: p3_qc_autocon_expon, p3_qc_accret_expon
       
  implicit none
  save

  public :: micro_p3_init, micro_p3_register, micro_p3_tend, &
            micro_p3_init_cnst, micro_p3_implements_cnst &
            ,micro_p3_readnl

  character(len=16), parameter :: unset_str = 'UNSET'

  private

  !Define indices for state%q constituents at module level so
  !defining them in micro_p3_register makes them permanently 
  !available.
  CHARACTER(len=16) :: precip_frac_method = 'max_overlap'  ! AaronDonahue, Hard-coded for now, should be fixed in the future

  integer, public ::    &
       ixcldliq = -1,   & ! cloud liquid amount index
       ixcldice = -1,      & ! ice index
       ixnumliq = -1,   & ! cloud liquid number index
       ixnumice = -1,   & ! cloud ice number index
       ixrain   = -1,   & ! rain index
       ixnumrain= -1,   & ! rain number index
       ixcldrim = -1,      & ! rime index ??
       ixrimvol  = -1,  & ! rime volume index ??
       ixqm  = -1      ! ?? index ??

!! pbuf 
   integer :: &
      cldo_idx,           &
      qme_idx,            &
      precip_total_tend_idx,          &
      nevapr_idx,         &
      dei_idx,            &
      rate1_cw2pr_st_idx, &
      mu_idx,             &
      lambdac_idx,        &
      rei_idx,            &
      rel_idx,            &
      ls_flxprc_idx,      &
      ls_flxsnw_idx,      &
      ls_reffrain_idx,    &
      ls_reffsnow_idx,    &
      cv_reffliq_idx,     &
      cv_reffice_idx,     &
      qr_evap_tend_idx,      &
      cmeliq_idx,         &
      relvar_idx,         &
      accre_enhan_idx     

! Physics buffer indices for fields registered by other modules
   integer :: &
      ast_idx = -1            

   integer :: &
      ni_activated_idx = -1,           &
      npccn_idx = -1,          &
      prec_str_idx = -1,       &
      prec_pcw_idx = -1,       &
      prec_sed_idx = -1,       &
      snow_str_idx = -1,       &
      snow_pcw_idx = -1,       &
      snow_sed_idx = -1

   real(rtype) :: &
      micro_mg_accre_enhan_fac = huge(1.0_rtype), & !Accretion enhancement factor from namelist
      prc_coef1_in             = huge(1.0_rtype), &
      prc_exp_in               = huge(1.0_rtype), &
      prc_exp1_in              = huge(1.0_rtype)

   integer :: ncnst

   character(len=8), parameter :: &      ! Constituent names
      cnst_names(8) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE', &
                      'RAINQM', 'CLDRIM','NUMRAI','BVRIM '/)

   character(len=128) :: micro_p3_lookup_dir     = unset_str ! location of p3 input files
   character(len=16)  :: micro_p3_tableversion   = unset_str ! P3 table version
   logical            :: micro_aerosolactivation = .false.   ! Use aerosol activation
   logical            :: micro_subgrid_cloud     = .false.   ! Use subgrid cloudiness
   logical            :: micro_tend_output       = .false.   ! Default microphysics tendencies to output file
   contains
!===============================================================================
subroutine micro_p3_readnl(nlfile)


  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'micro_p3_cam_readnl'

  namelist /micro_nl/ &
       micro_p3_tableversion, micro_p3_lookup_dir, micro_aerosolactivation, micro_subgrid_cloud, &
       micro_tend_output, p3_qc_autocon_expon, p3_qc_accret_expon

  !-----------------------------------------------------------------------------


#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(micro_p3_tableversion,   len(micro_p3_tableversion), mpichar, 0, mpicom)
  call mpibcast(micro_p3_lookup_dir,     len(micro_p3_lookup_dir),   mpichar, 0, mpicom)
  call mpibcast(micro_aerosolactivation, 1,                          mpilog,  0, mpicom)
  call mpibcast(micro_subgrid_cloud,     1,                          mpilog,  0, mpicom)
  call mpibcast(micro_tend_output,       1,                          mpilog,  0, mpicom)
  call mpibcast(p3_qc_autocon_expon,      1,                          mpir8,   0, mpicom)
  call mpibcast(p3_qc_accret_expon,       1,                          mpir8,   0, mpicom)

#endif

  ! Check to make sure p3 table version is valid
  select case (trim(micro_p3_tableversion))
    case ('4')
       ! Version 4 is valid
    case default
       print *, micro_p3_tableversion
       call bad_version_endrun()
  end select


contains
 
  subroutine bad_version_endrun
    ! Endrun wrapper with a more useful error message.
    character(len=128) :: errstring
    write(errstring,*) "Invalid version number specified for P3 microphysics: ", &
         micro_p3_tableversion
    call endrun(errstring)
  end subroutine bad_version_endrun

end subroutine micro_p3_readnl
  !================================================================================================

  subroutine micro_p3_register()

  logical :: prog_modal_aero ! prognostic aerosols


  call phys_getopts( prog_modal_aero_out   = prog_modal_aero )

   ncnst = 0
    ! Register Microphysics Constituents 
    ! (i.e. members of state%q) and save indices.
    ! TODO make sure the cnst_names match what we think they are here.
    !================
   call cnst_add(cnst_names(1), mwdry, cpair, 0._rtype, ixcldliq, &
         longname='Grid box averaged cloud liquid amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(2), mwdry, cpair, 0._rtype, ixcldice, &
         longname='Grid box averaged cloud ice amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(3), mwh2o, cpair, 0._rtype, ixnumliq, &
         longname='Grid box averaged cloud liquid number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(4), mwh2o, cpair, 0._rtype, ixnumice, &
         longname='Grid box averaged cloud ice number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(5), mwh2o, cpair, 0._rtype, ixrain, &
         longname='Grid box averaged rain amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(6), mwh2o, cpair, 0._rtype, ixcldrim, &
         longname='Grid box averaged riming amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(7), mwh2o, cpair, 0._rtype, ixnumrain, &
         longname='Grid box averaged rain number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(8), mwh2o, cpair, 0._rtype, ixrimvol, &
         longname='Grid box averaged riming volume', &
         is_convtran1=.true.)
   ncnst = ncnst + 1

  end subroutine micro_p3_register

  !================================================================================================
  function micro_p3_implements_cnst(name)

    ! Return true if specified constituent is implemented by the
    ! microphysics package

    character(len=*), intent(in) :: name        ! constituent name
    logical :: micro_p3_implements_cnst    ! return value

    micro_p3_implements_cnst = any(name == cnst_names)

  end function micro_p3_implements_cnst


  !================================================================================================

  subroutine micro_p3_init_cnst(name, q)

    ! Initialize the microphysics constituents, if they are
    ! not read from the initial file.

    character(len=*), intent(in) :: name     ! constituent name
    real(rtype), intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)

    if (micro_p3_implements_cnst(name)) q = 0.0_rtype

  end subroutine micro_p3_init_cnst

  !================================================================================================

  subroutine micro_p3_init()
    use micro_p3,       only: p3_init
    use micro_p3_utils, only: micro_p3_utils_init

    integer        :: m, mm
    integer        :: ierr
    logical :: history_amwg         ! output the variables used by the AMWG diag package
    logical :: history_verbose      ! produce verbose history output
    logical :: history_budget       ! Output tendencies and state variables for CAM4
    integer :: budget_histfile      ! output history file number for budget fields
                                   ! temperature, water vapor, cloud ice and cloud

    call micro_p3_utils_init(cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice, &
             cpliq,tmelt,pi,iulog,masterproc)



    ! INITIALIZE OUTPUT
    !==============
    do m = 1, ncnst
       call cnst_get_ind(cnst_names(m), mm)
       if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixcldrim /)) ) then
          ! mass mixing ratios
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', 'kg/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain /)) ) then
          ! number concentrations
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', '1/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else if ( mm == ixrimvol ) then
          ! number concentrations
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'm3/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', 'm3/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else
          call endrun( "micro_p3_acme_init: &
               &Could not call addfld for constituent with unknown units.")
       endif
    end do
      end if
   end if

  end subroutine micro_p3_init

  !================================================================================================
    subroutine get_cloud_fraction(its,ite,kts,kte,ast,qc,qr,qi,method, &
                  cld_frac_i,cld_frac_l,cld_frac_r)
      
       use micro_p3_utils, only: mincld, qsmall

       integer,intent(in)                                 :: its,ite,kts,kte
       real(rtype),dimension(its:ite,kts:kte),intent(in)  :: ast, qc, qr, qi
       character(len=16),intent(in)                       :: method
       real(rtype),dimension(its:ite,kts:kte),intent(out) :: cld_frac_i, cld_frac_l, cld_frac_r
       real(rtype),dimension(its:ite,kts:kte)             :: cldm

       integer  :: i,k
       integer  :: ktop, kbot, kdir

       call t_startf('micro_p3_get_cloud_fraction')
       ktop = kts        !k of top level
       kbot = kte        !k of bottom level
       kdir = -1         !(k: 1=top, nk=bottom)

       cldm(:,:)  = mincld
       cld_frac_i(:,:) = mincld
       cld_frac_l(:,:) = mincld
       do k = kbot,ktop,kdir
          do i=its,ite
             cldm(i,k)  = max(ast(i,k), mincld)
             cld_frac_i(i,k) = max(ast(i,k), mincld)
             cld_frac_l(i,k) = max(ast(i,k), mincld)
             cld_frac_r(i,k) = cldm(i,k)
          end do
       end do

       !! 
       !! precipitation fraction 
       !! 
       IF (trim(method) == 'in_cloud') THEN
          DO k = ktop-kdir,kbot,-kdir 
             DO i=its,ite
                ! in_cloud means that precip_frac (cld_frac_r) = cloud (cldm) frac when cloud mass
                ! is present. Below cloud, precip frac is equal to the cloud
                ! fraction from the last layer that had cloud. Since presence or
                ! absence of cloud is defined as mass > qsmall, sub-cloud precip
                ! frac for the in_cloud method tends to be very small and is
                ! very sensitive to tiny changes in condensate near cloud base.
                IF (qc(i,k) .lt. qsmall .and. qi(i,k) .lt. qsmall) THEN
                   ! max(cld_frac_r above and cld_frac_r for this layer) is taken here
                   ! because code is incapable of handling cld_frac_r<cldm for a
                   ! given grid cell
                   cld_frac_r(i,k) = max(cld_frac_r(i,k+kdir),cld_frac_r(i,k))
                END IF
             END DO !i
          END DO !k
       ELSE IF (trim(method) == 'max_overlap') THEN
       ! max overlap is the max cloud fraction in all layers above which are
       ! connected to this one by a continuous band of precip mass. If
       ! there's no precip mass falling into a cell, it's precip frac is equal
       ! to the cloud frac, which is probably ~zero.

       ! IF rain or ice mix ratios are smaller than threshold,
       ! then leave cld_frac_r as cloud fraction at current level
          DO k = ktop-kdir,kbot,-kdir 
             DO i=its,ite
                IF (qr(i,k+kdir) .ge. qsmall .or. qi(i,k+kdir) .ge. qsmall) THEN
                   cld_frac_r(i,k) = max(cld_frac_r(i,k+kdir),cld_frac_r(i,k))
                END IF
             END DO ! i
          END DO ! k
       END IF


       call t_stopf('micro_p3_get_cloud_fraction')
       return
    end subroutine get_cloud_fraction

  !================================================================================================
  subroutine micro_p3_tend(dtime, numcols)

    use micro_p3,       only: p3_main
    use micro_p3_utils, only: avg_diameter, &
                              rho_h2o, &
                              rho_h2os, &
                              qsmall, &
                              mincld, & 
                              inv_cp 

    !INPUT/OUTPUT VARIABLES
    real(rtype),                 intent(in)    :: dtime
    logical :: lq(pcnst)   !list of what constituents to update
    
    integer :: numcols
    !INTERNAL VARIABLES
    real(rtype) :: dz(pcols,pver)        !geometric layer thickness              m
    real(rtype) :: cldliq(pcols,pver)     !cloud liquid water mixing ratio        kg/kg
    real(rtype) :: numliq(pcols,pver)     !cloud liquid water drop concentraiton  #/kg
    real(rtype) :: rain(pcols,pver)       !rain water mixing ratio                kg/kg
    real(rtype) :: numrain(pcols,pver)    !rain water number concentration        #/kg
    real(rtype) :: qv(pcols,pver)         !water vapor mixing ratio               kg/kg
    real(rtype) :: ice(pcols,pver)        !total ice water mixing ratio           kg/kg
    real(rtype) :: qm(pcols,pver)      !rime ice mixing ratio                  kg/kg
    real(rtype) :: numice(pcols,pver)     !total ice crystal number concentration #/kg
    real(rtype) :: rimvol(pcols,pver)     !rime volume mixing ratio               m3/kg
    real(rtype) :: temp(pcols,pver)       !temperature copy needed for tendency   K
    real(rtype) :: th(pcols,pver)         !potential temperature                  K
    real(rtype) :: precip_liq_surf(pcols)         !precipitation rate, liquid             m s-1
    real(rtype) :: precip_ice_surf(pcols)         !precipitation rate, solid              m s-1

    real(rtype) :: rho_qi(pcols,pver)  !bulk density of ice                    kg m-1
    real(rtype) :: pres(pcols,pver)       !pressure at midlevel                   hPa
    real(rtype) :: cmeiout(pcols,pver)
    real(rtype) :: precip_liq_flux(pcols,pver+1)     !grid-box average rain flux (kg m^-2s^-1) pverp
    real(rtype) :: precip_ice_flux(pcols,pver+1)     !grid-box average ice/snow flux (kg m^-2s^-1) pverp
    real(rtype) :: exner(pcols,pver)      !exner formula for converting between potential and normal temp
    real(rtype) :: cld_frac_r(pcols,pver)      !rain cloud fraction
    real(rtype) :: cld_frac_l(pcols,pver)      !liquid cloud fraction
    real(rtype) :: cld_frac_i(pcols,pver)      !ice cloud fraction
    real(rtype) :: tend_out(pcols,pver,49) !microphysical tendencies
    real(rtype), dimension(pcols,pver) :: liq_ice_exchange ! sum of liq-ice phase change tendenices
    real(rtype), dimension(pcols,pver) :: vap_liq_exchange ! sum of vap-liq phase change tendenices
    real(rtype), dimension(pcols,pver) :: vap_ice_exchange ! sum of vap-ice phase change tendenices
    real(rtype) :: dummy_out(pcols,pver)    ! dummy_output variable for p3_main to replace unused variables.

    ! PBUF Variables
    real(rtype), pointer :: ast(:,:)      ! Relative humidity cloud fraction
    real(rtype), pointer :: ni_activated(:,:)     ! ice nucleation number
    real(rtype), pointer :: npccn(:,:)    ! liquid activation number tendency
    real(rtype), pointer :: cmeliq(:,:)
    !!
    real(rtype), pointer :: prec_str(:)    ! [Total] Sfc flux of precip from stratiform [ m/s ]
    real(rtype), pointer :: prec_sed(:)    ! Surface flux of total cloud water from sedimentation
    real(rtype), pointer :: prec_pcw(:)    ! Sfc flux of precip from microphysics [ m/s ]
    real(rtype), pointer :: snow_str(:)    ! [Total] Sfc flux of snow from stratiform   [ m/s ]
    real(rtype), pointer :: snow_pcw(:)    ! Sfc flux of snow from microphysics [ m/s ]
    real(rtype), pointer :: snow_sed(:)    ! Surface flux of cloud ice from sedimentation
    real(rtype), pointer :: relvar(:,:)    ! cloud liquid relative variance [-]
    real(rtype), pointer :: cldo(:,:)      ! Old cloud fraction
    real(rtype), pointer :: qr_evap_tend(:,:) ! precipitation evaporation rate 
    !! wetdep 
    real(rtype), pointer :: qme(:,:)
    real(rtype), pointer :: precip_total_tend(:,:)        ! Total precipitation (rain + snow)
    real(rtype), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
    !! COSP simulator
    real(rtype), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
    real(rtype), pointer :: rei(:,:)          ! Ice effective drop size (microns)
    real(rtype), pointer :: flxprc(:,:)     ! P3 grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
    real(rtype), pointer :: flxsnw(:,:)     ! P3 grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
    real(rtype), pointer :: reffrain(:,:)   ! P3 diagnostic rain effective radius (um)
    real(rtype), pointer :: reffsnow(:,:)   ! P3 diagnostic snow effective radius (um)
    real(rtype), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
    real(rtype), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)
    !! radiation 
    real(rtype), pointer :: dei(:,:)          ! Ice effective diameter (um)
    real(rtype), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
    real(rtype), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
    ! DONE PBUF
    ! For recording inputs/outputs to p3_main
    real(rtype) :: p3_main_inputs(pcols,pver+1,17) ! Record of inputs for p3_main
    real(rtype) :: p3_main_outputs(pcols,pver+1,31) ! Record of outputs for p3_main

    ! Derived Variables
    real(rtype) :: icimrst(pcols,pver) ! stratus ice mixing ratio - on grid
    real(rtype) :: icwmrst(pcols,pver) ! stratus water mixing ratio - on grid
    real(rtype) :: rho(pcols,pver)
    real(rtype) :: drout2(pcols,pver)
    real(rtype) :: reff_rain(pcols,pver)
    real(rtype) :: col_location(pcols,3),tmp_loc(pcols)  ! Array of column lon (index 1) and lat (index 2)
    integer     :: tmpi_loc(pcols) ! Global column index temp array

    ! Variables used for microphysics output
    real(rtype) :: aqrain(pcols,pver)
    real(rtype) :: anrain(pcols,pver)
    real(rtype) :: nfice(pcols,pver)
    real(rtype) :: efcout(pcols,pver)      
    real(rtype) :: efiout(pcols,pver)      
    real(rtype) :: ncout(pcols,pver)      
    real(rtype) :: niout(pcols,pver)      
    real(rtype) :: freqr(pcols,pver)      
    real(rtype) :: freql(pcols,pver)      
    real(rtype) :: freqi(pcols,pver)      
    real(rtype) :: cdnumc(pcols)      
    real(rtype) :: icinc(pcols,pver) 
    real(rtype) :: icwnc(pcols,pver) 

 
    integer :: it                      !timestep counter                       -
    integer :: its, ite                !horizontal bounds (column start,finish)
    integer :: kts                     !closest level to TOM                   -
    integer :: kte                     !near surface level                     -

    logical :: do_predict_nc           !prognostic droplet concentration or not?
    logical :: do_subgrid_clouds       !use subgrid cloudiness in tendency calculations?
    integer :: icol, ncol, k
    integer :: psetcols, lchnk
    integer :: itim_old

    ! For rrtmg optics. specified distribution.
    real(rtype), parameter :: dcon   = 25.e-6_rtype         ! Convective size distribution effective radius (um)
    real(rtype), parameter :: mucon  = 5.3_rtype            ! Convective size distribution shape parameter
    real(rtype), parameter :: deicon = 50._rtype            ! Convective ice effective diameter (um)

    call t_startf('micro_p3_tend_init')
 
    !psetcols = state%psetcols
    !lchnk = state%lchnk

    !+++ Aaron Donahue
    itim_old = pbuf_old_tim_idx()

    !============================ 
    ! All external PBUF variables:
    ! INPUTS
    ! AaronDonahue:  All of these inputs will need be passed from AD to P3
    ! through interface.
    call pbuf_get_field(pbuf, ast_idx,         ast, start=(/1,1,itim_old/), kount=(/psetcols,pver,1/))
    call pbuf_get_field(pbuf, ni_activated_idx,        ni_activated                                                  ) 
    call pbuf_get_field(pbuf, npccn_idx,       npccn                                                 )
    ! AaronDonahue - will need to be passed as num_cols.
    !ncol = state%ncol
    !==============
    ! Some pre-microphysics INITIALIZATION
    !==============
    ! HANDLE AEROSOL ACTIVATION
    !==============
    ! AaronDonahue - do_predict_nc will be an explcit input.
    do_predict_nc = micro_aerosolactivation 

    ! COMPUTE GEOMETRIC THICKNESS OF GRID & CONVERT T TO POTENTIAL TEMPERATURE
    !==============
    ! AaronDonahue - We will caclulate these values in the interface, just like
    ! done here.  Note, that it is possible that "th" may be a state variable in
    ! scream, if so, we won't have to calculate it.  For now, go ahead.
    exner(:ncol,:pver) = 1._rtype/((state%pmid(:ncol,:pver)*1.e-5_rtype)**(rair*inv_cp))
    do icol = 1,ncol
       do k = 1,pver
! Note: dz is calculated in the opposite direction that pdel is calculated,
! thus when considering any dp/dz calculation we must also change the sign.
          dz(icol,k) = state%zi(icol,k) - state%zi(icol,k+1)
          th(icol,k)  = state%t(icol,k)*exner(icol,k) !/(state%pmid(icol,k)*1.e-5)**(rd*inv_cp) 
       end do
    end do

    ! ASSIGN TOP AND BOTTOM INDICES FOR GRID
    !==============
    !kts is closest level to top of model. Instead of 1 (top-of-model), 
    !we use the previously-defined trop_cloud_top_lev to reduce the number of 
    !levels we need to calculate and to avoid upper-atmos regions where this
    !micro-physics is inappropriate. kte is the near-surface level = pver.

    ! AaronDonahue - In scream these will just be kts=1 and kte=num_lev
    kts=1
    kte=num_lev 

    ! HANDLE TIMESTEP COUNTER, GET LON/LAT values (used by error warning code)
    !==============
    !p3 wants to know the timestep number (which it calls "it") because it
    !handles things differently on the first step, where it doesn't have values
    !yet. E3SM has a handy function for deciding if this is the first step, so 
    !we hack "it" with "is_first_step()" for now. Eventually, we should replace
    !"it" with a logical.
    ! AaronDonahue - "it" will be an explicit input from AD.
    it = get_nstep()
!ASD - Ignore For Now    tmp_loc =-999.0_rtype
!ASD - Ignore For Now    call get_rlon_all_p(lchnk,ncol,tmp_loc)
!ASD - Ignore For Now    col_location(:ncol,2) = tmp_loc(:ncol)*180.0_rtype/pi
!ASD - Ignore For Now    call get_rlat_all_p(lchnk,ncol,tmp_loc)
!ASD - Ignore For Now    col_location(:ncol,3) = tmp_loc(:ncol)*180.0_rtype/pi
!ASD - Ignore For Now    call get_gcol_all_p(lchnk,ncol,tmpi_loc)
!ASD - Ignore For Now    col_location(:ncol,1) = real(tmpi_loc(:ncol))

    ! MAKE LOCAL COPIES OF VARS MODIFIED BY P3
    !==============
    !local copies are needed because state is passed into this routine as intent=in
    !while P3 seeks to modify state variables in-place. Also, we need a copy of 
    !old values in order to back out ptend values later. Traditionally, a local copy 
    !is created by copying the whole state. It is much cheaper to just copy the 
    !variables we need. 
    ! AaronDonahue - the following just parses the "q" array into the individual
    ! components that P3 cares about.  We will want to do the same.  Note
    ! instead of having ixXXXX we will need to come up with another way to make
    ! sure we have the right index for each entry in the array of "q".
    cldliq  = state%q(:,:,ixcldliq)
    numliq  = state%q(:,:,ixnumliq)
    rain    = state%q(:,:,ixrain)
    numrain = state%q(:,:,ixnumrain)
    qv      = state%q(:,:,1)
    ice     = state%q(:,:,ixcldice)
    qm   = state%q(:,:,ixcldrim) !Aaron, changed ixqm to ixcldrim to match Kai's code
    numice  = state%q(:,:,ixnumice)
    rimvol  = state%q(:,:,ixrimvol)
    ! AaronDonahue - these four will be explicit inputs from AD.
    its     = 1
    ite     = state%ncol
    kts     = 1
    kte     = pver
    ! AaronDonahue - this is just a renaming for asthetics.  We can just pass
    ! the pressure directly to p3_main from the input passed from the AD.
    pres    = state%pmid(:,:)
    ! Initialize the raidation dependent variables.
    ! AaronDonahue - We will need these to be initialized inside our 
    mu      = 0.0_rtype !mucon
    lambdac = 0.0_rtype !(mucon + 1._rtype)/dcon
!ASD - Ignore For Now    dei     = 50.0_rtype !deicon
    ! Determine the cloud fraction and precip cover
    ! AaronDonahue - We will also want to create a subroutine in our interface
    ! that calculates the cloud_fraction like this.  We can set
    ! "do_subgrid_clouds=true" for now and maybe change if there is a need
    ! later.
    ! Note: we may want to ping the whole group about whether or not the cloud
    ! fraction routine really belongs as a universal function rather than have
    ! it's home in P3.
    cld_frac_i(:,:) = 1.0_rtype
    cld_frac_l(:,:) = 1.0_rtype
    cld_frac_r(:,:) = 1.0_rtype
    do_subgrid_clouds = micro_subgrid_cloud
    if (do_subgrid_clouds) &
        call get_cloud_fraction(its,ite,kts,kte,ast(its:ite,kts:kte),cldliq(its:ite,kts:kte), &
                rain(its:ite,kts:kte),ice(its:ite,kts:kte),precip_frac_method, &
                cld_frac_i(its:ite,kts:kte),cld_frac_l(its:ite,kts:kte),cld_frac_r(its:ite,kts:kte))
    ! AaronDonahue - You can ignore all t_stopf and t_startf calls (in case I
    ! miss any.)
!ASD - Ignore For Now    call t_stopf('micro_p3_tend_init')
!ASD - Ignore For Now
!ASD - Ignore For Now    p3_main_inputs(:,:,:) = -999._rtype
!ASD - Ignore For Now    do k = 1,pver
!ASD - Ignore For Now      p3_main_inputs(1,k,1)  = ast(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,2)  = ni_activated(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,3)  = npccn(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,4)  = state%pmid(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,5)  = state%zi(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,6)  = state%T(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,7)  = qv(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,8)  = cldliq(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,9)  = ice(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,10) = numliq(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,11) = numice(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,12) = rain(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,13) = numrain(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,14) = qm(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,15) = rimvol(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,16) = state%pdel(1,k)
!ASD - Ignore For Now      p3_main_inputs(1,k,17) = relvar(1,k)
!ASD - Ignore For Now    end do
!ASD - Ignore For Now    p3_main_inputs(1,pver+1,5) = state%zi(1,pver+1)

    ! CALL P3
    !==============
    ! TODO: get proper value for 'it' from time module

    ! AaronDonahue - These initializations might not be needed, but always good
    ! to initialize.  Note: dummy_out is only a F90 thing for in/out's that we
    ! removed in the C++ version.  You shouldn't have anything "dummy" in the
    ! C++ code. 
!ASD - Ignore For Now    dummy_out(:,:) = 0.0_rtype
    precip_liq_surf = 0.0_rtype
    precip_ice_surf = 0.0_rtype
    prec_pcw = 0.0_rtype
    snow_pcw = 0.0_rtype
    vap_liq_exchange = 0.0_rtype

!ASD - Ignore For Now    call t_startf('micro_p3_tend_loop')
    ! AaronDonahue - of course, the C++ interface will call p3_main in a
    ! different way using structures.
    call p3_main( &
         cldliq(its:ite,kts:kte),     & ! INOUT  cloud, mass mixing ratio         kg kg-1
         numliq(its:ite,kts:kte),     & ! INOUT  cloud, number mixing ratio       #  kg-1
         rain(its:ite,kts:kte),       & ! INOUT  rain, mass mixing ratio          kg kg-1
         numrain(its:ite,kts:kte),    & ! INOUT  rain, number mixing ratio        #  kg-1
         th(its:ite,kts:kte),         & ! INOUT  potential temperature            K
         qv(its:ite,kts:kte),         & ! INOUT  water vapor mixing ratio         kg kg-1
         dtime,                       & ! IN     model time step                  s
         ice(its:ite,kts:kte),        & ! INOUT  ice, total mass mixing ratio     kg kg-1
         qm(its:ite,kts:kte),      & ! INOUT  ice, rime mass mixing ratio      kg kg-1
         numice(its:ite,kts:kte),     & ! INOUT  ice, total number mixing ratio   #  kg-1
         rimvol(its:ite,kts:kte),     & ! INOUT  ice, rime volume mixing ratio    m3 kg-1
         pres(its:ite,kts:kte),       & ! IN     pressure at cell midpoints       Pa
         dz(its:ite,kts:kte),        & ! IN     vertical grid spacing            m
         npccn(its:ite,kts:kte),      & ! IN ccn activation number tendency kg-1 s-1
         ni_activated(its:ite,kts:kte),       & ! IN activated ice nuclei concentration kg-1
         relvar(its:ite,kts:kte),     & ! IN cloud liquid relative variance
         it,                          & ! IN     time step counter NOTE: starts at 1 for first time step
         precip_liq_surf(its:ite),            & ! OUT    surface liquid precip rate       m s-1
         precip_ice_surf(its:ite),            & ! OUT    surface frozen precip rate       m s-1
         its,                         & ! IN     horizontal index lower bound     -
         ite,                         & ! IN     horizontal index upper bound     -
         kts,                         & ! IN     vertical index lower bound       -
         kte,                         & ! IN     vertical index upper bound       -
         rel(its:ite,kts:kte),        & ! OUT    effective radius, cloud          m
         rei(its:ite,kts:kte),        & ! OUT    effective radius, ice            m
         rho_qi(its:ite,kts:kte),  & ! OUT    bulk density of ice              kg m-3
         do_predict_nc,               & ! IN     .true.=prognostic Nc, .false.=specified Nc
         ! AaronDonahue new stuff
         state%pdel(its:ite,kts:kte), & ! IN pressure level thickness for computing total mass
         exner(its:ite,kts:kte),      & ! IN exner values
         cmeiout(its:ite,kts:kte),    & ! OUT Deposition/sublimation rate of cloud ice 
         precip_total_tend(its:ite,kts:kte),      & ! OUT Total precipitation (rain + snow)
         nevapr(its:ite,kts:kte),     & ! OUT evaporation of total precipitation (rain + snow)
         qr_evap_tend(its:ite,kts:kte),  & ! OUT rain evaporation
         precip_liq_flux(its:ite,kts:kte+1),     & ! OUT grid-box average rain flux (kg m^-2s^-1) pverp 
         precip_ice_flux(its:ite,kts:kte+1),     & ! OUT grid-box average ice/snow flux (kgm^-2 s^-1) pverp
         cld_frac_r(its:ite,kts:kte),      & ! IN rain cloud fraction
         cld_frac_l(its:ite,kts:kte),      & ! IN liquid cloud fraction
         cld_frac_i(its:ite,kts:kte),      & ! IN ice cloud fraction
         tend_out(its:ite,kts:kte,:), & ! OUT p3 microphysics tendencies
         mu(its:ite,kts:kte),         & ! OUT Size distribution shape parameter for radiation
         lambdac(its:ite,kts:kte),    & ! OUT Size distribution slope parameter for radiation
         liq_ice_exchange(its:ite,kts:kte),& ! OUT sum of liq-ice phase change tendenices   
         vap_liq_exchange(its:ite,kts:kte),& ! OUT sun of vap-liq phase change tendencies
         vap_ice_exchange(its:ite,kts:kte),& ! OUT sum of vap-ice phase change tendencies
         col_location(its:ite,:3)          & ! IN column locations
         )

!ASD - Ignore For Now    p3_main_outputs(:,:,:) = -999._rtype
!ASD - Ignore For Now    do k = 1,pver
!ASD - Ignore For Now      p3_main_outputs(1,k, 1) = cldliq(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k, 2) = numliq(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k, 3) = rain(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k, 4) = numrain(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k, 5) = th(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k, 6) = qv(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k, 7) = ice(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k, 8) = qm(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k, 9) = numice(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,10) = rimvol(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,14) = rel(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,15) = rei(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,18) = rho_qi(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,19) = cmeiout(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,20) = precip_total_tend(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,21) = nevapr(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,22) = qr_evap_tend(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,23) = precip_liq_flux(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,24) = precip_ice_flux(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,27) = mu(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,28) = lambdac(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,29) = liq_ice_exchange(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,30) = vap_liq_exchange(1,k)
!ASD - Ignore For Now      p3_main_outputs(1,k,31) = vap_ice_exchange(1,k)
!ASD - Ignore For Now    end do
!ASD - Ignore For Now    p3_main_outputs(1,1,11) = precip_liq_surf(1)
!ASD - Ignore For Now    p3_main_outputs(1,1,12) = precip_ice_surf(1)
!ASD - Ignore For Now    p3_main_outputs(1,pver+1,23) = precip_liq_flux(1,pver+1)
!ASD - Ignore For Now    p3_main_outputs(1,pver+1,24) = precip_ice_flux(1,pver+1)
!ASD - Ignore For Now    call outfld('P3_input',  p3_main_inputs,  pcols, lchnk)
!ASD - Ignore For Now    call outfld('P3_output', p3_main_outputs, pcols, lchnk)

    !MASSAGE OUTPUT TO FIT E3SM EXPECTATIONS
    !============= 

    !TODO: figure out what else other E3SM parameterizations need from micro and make sure 
    !they are assigned here. The comments below are a step in that direction.


    !cloud_rad_props also uses snow radiative properties which aren't available from 
    !P3 (perhaps because ice phase in p3 includes *all* ice already?).

    !BACK OUT TENDENCIES FROM STATE CHANGES
    !=============
    ! AaronDonahue - In the F90 implementation we are having P3 calculate
    ! tendencies.  In the SCREAM version we will have P3 pass the updated state.
    ! So this can be ignored, except maybe the need to update the temperature
    ! state "T_atm" (here "temp") from the "th" value produced by P3.
    temp(:ncol,:pver) = th(:ncol,:pver)/exner(:ncol,:pver) 
!ASD - Ignore For Now    ptend%s(:ncol,:pver)           = cpair*( temp(:ncol,:pver) - state%t(:ncol,:pver) )/dtime 
!ASD - Ignore For Now    ptend%q(:ncol,:pver,1)         = ( max(0._rtype,qv(:ncol,:pver)     ) - state%q(:ncol,:pver,1)         )/dtime
!ASD - Ignore For Now    ptend%q(:ncol,:pver,ixcldliq)  = ( max(0._rtype,cldliq(:ncol,:pver) ) - state%q(:ncol,:pver,ixcldliq)  )/dtime
!ASD - Ignore For Now    ptend%q(:ncol,:pver,ixnumliq)  = ( max(0._rtype,numliq(:ncol,:pver) ) - state%q(:ncol,:pver,ixnumliq)  )/dtime
!ASD - Ignore For Now    ptend%q(:ncol,:pver,ixrain)    = ( max(0._rtype,rain(:ncol,:pver)   ) - state%q(:ncol,:pver,ixrain)    )/dtime
!ASD - Ignore For Now    ptend%q(:ncol,:pver,ixnumrain) = ( max(0._rtype,numrain(:ncol,:pver)) - state%q(:ncol,:pver,ixnumrain) )/dtime
!ASD - Ignore For Now    ptend%q(:ncol,:pver,ixcldice)  = ( max(0._rtype,ice(:ncol,:pver)    ) - state%q(:ncol,:pver,ixcldice)  )/dtime
!ASD - Ignore For Now    ptend%q(:ncol,:pver,ixnumice)  = ( max(0._rtype,numice(:ncol,:pver) ) - state%q(:ncol,:pver,ixnumice)  )/dtime
!ASD - Ignore For Now    ptend%q(:ncol,:pver,ixcldrim)  = ( max(0._rtype,qm(:ncol,:pver)  ) - state%q(:ncol,:pver,ixcldrim)  )/dtime
!ASD - Ignore For Now    ptend%q(:ncol,:pver,ixrimvol)  = ( max(0._rtype,rimvol(:ncol,:pver) ) - state%q(:ncol,:pver,ixrimvol)  )/dtime

!ASD - Ignore For Now    call t_stopf('micro_p3_tend_loop')
!ASD - Ignore For Now    call t_startf('micro_p3_tend_finish')
   ! Following MG interface as a template:

    ! AaronDonahue - The following are a set of values that P3 would calculate
    ! for use by other processes.  I think we should skip these for now and only
    ! add them back into our interface as needed.
    ! Note/Thought:  This post-processing of variables is meant to provide
    ! information for the likes of macrophysics, the energy checker, COSP and
    ! radiation.  I think it would be great for our interface if we set up
    ! subroutines for each.  That way if you are focused on radiation (for
    ! example) you can jump right to the radiation post-processing code.  I
    ! include this thought here so we don't lose it.
!ASD - Ignore For Now    ! Net micro_p3 condensation rate
!ASD - Ignore For Now    qme(:ncol,top_lev:pver) = cmeliq(:ncol,top_lev:pver) + cmeiout(:ncol,top_lev:pver)  ! cmeiout is output from p3 micro
!ASD - Ignore For Now    ! Add cmeliq to  vap_liq_exchange
!ASD - Ignore For Now    vap_liq_exchange(:ncol,top_lev:pver) = vap_liq_exchange(:ncol,top_lev:pver) + cmeliq(:ncol,top_lev:pver) 
!ASD - Ignore For Now
!ASD - Ignore For Now!====================== Export variables/Conservation START ======================!
!ASD - Ignore For Now     !For precip, accumulate only total precip in prec_pcw and snow_pcw variables.
!ASD - Ignore For Now    ! Other precip output variables are set to 0
!ASD - Ignore For Now    ! Do not subscript by ncol here, because in physpkg we divide the whole
!ASD - Ignore For Now    ! array and need to avoid an FPE due to uninitialized data.
!ASD - Ignore For Now    prec_pcw = precip_liq_surf + precip_ice_surf
!ASD - Ignore For Now    prec_sed = 0._rtype
!ASD - Ignore For Now    prec_str = prec_pcw + prec_sed
!ASD - Ignore For Now
!ASD - Ignore For Now    snow_pcw = precip_ice_surf
!ASD - Ignore For Now    snow_sed = 0._rtype
!ASD - Ignore For Now    snow_str = snow_pcw + snow_sed
!ASD - Ignore For Now!====================== Export variables/Conservation END ======================!
!ASD - Ignore For Now
!ASD - Ignore For Now!====================== Radiation Specific Outputs START ======================!
!ASD - Ignore For Now
!ASD - Ignore For Now   ! Calculate rho for size distribution
!ASD - Ignore For Now   ! parameter calculations and average it if needed
!ASD - Ignore For Now   
!ASD - Ignore For Now   rho(:ncol,top_lev:) = &
!ASD - Ignore For Now      state%pmid(:ncol,top_lev:) / (rair*temp(:ncol,top_lev:))
!ASD - Ignore For Now   ! ------------------------------------------------------------ !
!ASD - Ignore For Now   ! Compute in cloud ice and liquid mixing ratios                !
!ASD - Ignore For Now   ! Note that 'iclwp, iciwp' are used for radiation computation. !
!ASD - Ignore For Now   ! ------------------------------------------------------------ !
!ASD - Ignore For Now      
!ASD - Ignore For Now   icinc = 0._rtype
!ASD - Ignore For Now   icwnc = 0._rtype
!ASD - Ignore For Now      
!ASD - Ignore For Now   do k = top_lev, pver
!ASD - Ignore For Now      do icol = 1, ncol
!ASD - Ignore For Now         ! Limits for in-cloud mixing ratios consistent with P3 microphysics
!ASD - Ignore For Now         ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
!ASD - Ignore For Now         icimrst(icol,k)   = min( state%q(icol,k,ixcldice) / max(mincld,cld_frac_i(icol,k)),0.005_rtype )
!ASD - Ignore For Now         icwmrst(icol,k)   = min( state%q(icol,k,ixcldliq) / max(mincld,cld_frac_l(icol,k)),0.005_rtype )
!ASD - Ignore For Now         icinc(icol,k)     = state%q(icol,k,ixnumice) / max(mincld,cld_frac_i(icol,k)) * &
!ASD - Ignore For Now              state%pmid(icol,k) / (287.15_rtype*state%t(icol,k))
!ASD - Ignore For Now         icwnc(icol,k)     = state%q(icol,k,ixnumliq) / max(mincld,cld_frac_l(icol,k)) * &
!ASD - Ignore For Now              state%pmid(icol,k) / (287.15_rtype*state%t(icol,k))
!ASD - Ignore For Now      end do                    
!ASD - Ignore For Now   end do
!ASD - Ignore For Now
!ASD - Ignore For Now   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ASD - Ignore For Now   !!
!ASD - Ignore For Now   !! derived fields 
!ASD - Ignore For Now   !!
!ASD - Ignore For Now   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ASD - Ignore For Now    !cloud_rad_props needs ice effective diameter, which Kai calculates as below:
!ASD - Ignore For Now    !   dei = rei*diag_rhopo(i,k,iice)/rho_h2os*2._rtype
!ASD - Ignore For Now    !where rhopo is bulk ice density from table lookup (taken from table_val_ice_bulk_dens, here written as rho_qi) and rho_h2os=917.0 is a constant parameter.
!ASD - Ignore For Now   !! Effective radius for cloud liquid
!ASD - Ignore For Now   rel(:ncol,top_lev:) = rel(:ncol,top_lev:) * 1e6_rtype  ! Rescale rel to be in microns
!ASD - Ignore For Now   !! Effective radius for cloud ice
!ASD - Ignore For Now   rei(:ncol,top_lev:) = rei(:ncol,top_lev:) * 1e6_rtype  ! Rescale rei to be in microns
!ASD - Ignore For Now   !! Effective diameter for cloud ice
!ASD - Ignore For Now   dei(:ncol,top_lev:) = rei(:ncol,top_lev:) * 2._rtype
!ASD - Ignore For Now
!ASD - Ignore For Now   !!
!ASD - Ignore For Now   !! Limiters for low cloud fraction
!ASD - Ignore For Now   !!
!ASD - Ignore For Now   
!ASD - Ignore For Now   do k = top_lev, pver
!ASD - Ignore For Now      do icol = 1, ncol
!ASD - Ignore For Now         if ( ast(icol,k) < 1.e-4_rtype ) then
!ASD - Ignore For Now            mu(icol,k) = mucon
!ASD - Ignore For Now            lambdac(icol,k) = (mucon + 1._rtype)/dcon
!ASD - Ignore For Now            dei(icol,k) = deicon
!ASD - Ignore For Now         end if
!ASD - Ignore For Now      end do
!ASD - Ignore For Now   end do
!ASD - Ignore For Now
!ASD - Ignore For Now   !!
!ASD - Ignore For Now   !! New output fields
!ASD - Ignore For Now   !!
!ASD - Ignore For Now   efcout      = 0._rtype
!ASD - Ignore For Now   efiout      = 0._rtype
!ASD - Ignore For Now   ncout       = 0._rtype
!ASD - Ignore For Now   niout       = 0._rtype
!ASD - Ignore For Now   freql       = 0._rtype
!ASD - Ignore For Now   freqi       = 0._rtype
!ASD - Ignore For Now   cdnumc      = 0._rtype
!ASD - Ignore For Now   nfice       = 0._rtype
!ASD - Ignore For Now
!ASD - Ignore For Now   ! FICE
!ASD - Ignore For Now   do k = top_lev, pver
!ASD - Ignore For Now      do icol = 1, ncol
!ASD - Ignore For Now      if (ice(icol,k).gt.qsmall .and. (rain(icol,k)+ice(icol,k)+cldliq(icol,k)).gt.qsmall) then
!ASD - Ignore For Now         nfice(icol,k)=min(ice(icol,k)/(rain(icol,k)+ice(icol,k)+cldliq(icol,k)),1._rtype)
!ASD - Ignore For Now      else
!ASD - Ignore For Now         nfice(icol,k)=0._rtype
!ASD - Ignore For Now      end if
!ASD - Ignore For Now      end do
!ASD - Ignore For Now   end do
!ASD - Ignore For Now
!ASD - Ignore For Now   ! Column droplet concentration
!ASD - Ignore For Now   cdnumc(:ncol) = sum(numliq(:ncol,top_lev:pver) * &
!ASD - Ignore For Now        state%pdel(:ncol,top_lev:pver)/gravit, dim=2)
!ASD - Ignore For Now   do k = top_lev, pver
!ASD - Ignore For Now      do icol = 1, ncol
!ASD - Ignore For Now         if ( cld_frac_l(icol,k) > 0.01_rtype .and. icwmrst(icol,k) > 5.e-5_rtype ) then
!ASD - Ignore For Now            efcout(icol,k) = rel(icol,k) * cld_frac_l(icol,k)
!ASD - Ignore For Now            ncout(icol,k)  = icwnc(icol,k) * cld_frac_l(icol,k)
!ASD - Ignore For Now            freql(icol,k)  = cld_frac_l(icol,k)
!ASD - Ignore For Now         end if
!ASD - Ignore For Now         if ( cld_frac_i(icol,k) > 0.01_rtype .and. icimrst(icol,k) > 1.e-6_rtype ) then
!ASD - Ignore For Now            efiout(icol,k) = rei(icol,k) * cld_frac_i(icol,k)
!ASD - Ignore For Now            niout(icol,k)  = icinc(icol,k) * cld_frac_i(icol,k)
!ASD - Ignore For Now            freqi(icol,k)  = cld_frac_i(icol,k)
!ASD - Ignore For Now         end if
!ASD - Ignore For Now      end do
!ASD - Ignore For Now   end do
!ASD - Ignore For Now
!ASD - Ignore For Now   ! note: 1e-6 kgho2/kgair/s * 1000. pa / (9.81 m/s2) / 1000 kgh2o/m3 = 1e-7 m/s
!ASD - Ignore For Now   ! this is 1ppmv of h2o in 10hpa
!ASD - Ignore For Now   ! alternatively: 0.1 mm/day * 1.e-4 m/mm * 1/86400 day/s = 1.e-9
!ASD - Ignore For Now 
!ASD - Ignore For Now   !!
!ASD - Ignore For Now   !! Rain/Snow effective diameter
!ASD - Ignore For Now   !!
!ASD - Ignore For Now   drout2    = 0._rtype
!ASD - Ignore For Now   reff_rain = 0._rtype
!ASD - Ignore For Now   aqrain    = 0._rtype
!ASD - Ignore For Now   anrain    = 0._rtype
!ASD - Ignore For Now   freqr     = 0._rtype
!ASD - Ignore For Now   ! Prognostic precipitation
!ASD - Ignore For Now   where (rain(:ncol,top_lev:) >= 1.e-7_rtype)
!ASD - Ignore For Now      drout2(:ncol,top_lev:) = avg_diameter( &
!ASD - Ignore For Now           rain(:ncol,top_lev:), &
!ASD - Ignore For Now           numrain(:ncol,top_lev:) * rho(:ncol,top_lev:), &
!ASD - Ignore For Now           rho(:ncol,top_lev:), rho_h2o)
!ASD - Ignore For Now
!ASD - Ignore For Now      aqrain = rain * cld_frac_r
!ASD - Ignore For Now      anrain = numrain * cld_frac_r
!ASD - Ignore For Now      freqr = cld_frac_r
!ASD - Ignore For Now      reff_rain(:ncol,top_lev:) = drout2(:ncol,top_lev:) * &
!ASD - Ignore For Now           1.5_rtype * 1.e6_rtype
!ASD - Ignore For Now   end where
!ASD - Ignore For Now
!ASD - Ignore For Now!====================== COSP Specific Outputs START ======================!
!ASD - Ignore For Now! LS_FLXPRC, LS_FLXSNW, LS_REFFRAIN, LS_REFFSNOW, CV_REFFLIQ, CV_REFFICE
!ASD - Ignore For Now   !== Grid-box mean flux_large_scale_cloud at interfaces (kg/m2/s)
!ASD - Ignore For Now! flxprc and flxsnw are used in COSP to compute precipitation fractional
!ASD - Ignore For Now! area and derive precipitation (rain and snow) mixing ratios. Including iflx
!ASD - Ignore For Now! and cflx in precipitation fluxes would result in additional effects of cloud liquid and
!ASD - Ignore For Now! ice on cosp's smiluated lidar and radar reflectivity signal through the rain/snow
!ASD - Ignore For Now! portion of calculations that are handled separately from that of cloud liquid
!ASD - Ignore For Now! and ice. If included, it would not exactly amount to double counting the effect of
!ASD - Ignore For Now! cloud liquid and ice because the mixing ratio derived from iflx and cflx epected to be much smaller
!ASD - Ignore For Now! than the actual grid-mean cldliq and cldice, and rain or snow size distribution
!ASD - Ignore For Now! would be used to compute the lidar/radar signal strength.
!ASD - Ignore For Now! 
!ASD - Ignore For Now! Note that it would need to include iflx and cflx to make the values at surface
!ASD - Ignore For Now! interface consistent with large scale precipitation rates.
!ASD - Ignore For Now
!ASD - Ignore For Now    ! array must be zeroed beyond trop_cloud_top_pre otherwise undefined values will be used in cosp.
!ASD - Ignore For Now    flxprc(:ncol,1:top_lev) = 0.0_rtype ! Rain+Snow
!ASD - Ignore For Now    flxsnw(:ncol,1:top_lev) = 0.0_rtype ! Snow
!ASD - Ignore For Now 
!ASD - Ignore For Now    flxprc(:ncol,top_lev:pverp) = precip_liq_flux(:ncol,top_lev:pverp) + precip_ice_flux(:ncol,top_lev:pverp) ! need output from p3
!ASD - Ignore For Now    flxsnw(:ncol,top_lev:pverp) = precip_ice_flux(:ncol,top_lev:pverp) ! need output from p3
!ASD - Ignore For Now
!ASD - Ignore For Now    cvreffliq(:ncol,top_lev:pver) = 9.0_rtype
!ASD - Ignore For Now    cvreffice(:ncol,top_lev:pver) = 37.0_rtype
!ASD - Ignore For Now
!ASD - Ignore For Now    reffrain(:,:) = 0._rtype
!ASD - Ignore For Now    reffsnow(:,:) = 0._rtype
!ASD - Ignore For Now    reffrain(:ncol,top_lev:pver) = reff_rain(:ncol,top_lev:pver)
!ASD - Ignore For Now    reffsnow(:ncol,top_lev:pver) = 1000._rtype !! dummy value, the choice here impacts the COSP output variable: CFAD_DBZE94_CS.  TODO: Figure out if this is ok, change if needed.
!ASD - Ignore For Now
!ASD - Ignore For Now!====================== COSP Specific Outputs  END ======================!
!ASD - Ignore For Now
!ASD - Ignore For Now    !WRITE OUTPUT
!ASD - Ignore For Now    !=============
!ASD - Ignore For Now   call outfld('AQRAIN',      aqrain,      psetcols, lchnk)
!ASD - Ignore For Now   call outfld('ANRAIN',      anrain,      psetcols, lchnk)
!ASD - Ignore For Now   call outfld('AREL',        efcout,      pcols,    lchnk)
!ASD - Ignore For Now   call outfld('AREI',        efiout,      pcols,    lchnk) 
!ASD - Ignore For Now   call outfld('AWNC' ,       ncout,       pcols,    lchnk)
!ASD - Ignore For Now   call outfld('AWNI' ,       niout,       pcols,    lchnk)
!ASD - Ignore For Now   call outfld('FICE',        nfice,       psetcols, lchnk)
!ASD - Ignore For Now   call outfld('FREQL',       freql,       pcols,    lchnk)
!ASD - Ignore For Now   call outfld('FREQI',       freqi,       pcols,    lchnk)
!ASD - Ignore For Now   call outfld('FREQR',       freqr,       psetcols, lchnk)
!ASD - Ignore For Now   call outfld('CDNUMC',      cdnumc,      pcols,    lchnk)
!ASD - Ignore For Now
!ASD - Ignore For Now   call outfld('CLOUDFRAC_LIQ_MICRO',  cld_frac_l,      pcols, lchnk)
!ASD - Ignore For Now   call outfld('CLOUDFRAC_ICE_MICRO',  cld_frac_i,      pcols, lchnk)
!ASD - Ignore For Now   call outfld('CLOUDFRAC_RAIN_MICRO', cld_frac_r,      pcols, lchnk)
!ASD - Ignore For Now
!ASD - Ignore For Now   ! Write p3 tendencies as output 
!ASD - Ignore For Now   ! warm-phase process rates
!ASD - Ignore For Now   call outfld('P3_qrcon',  tend_out(:,:, 1), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qc2qr_accret_tend',  tend_out(:,:, 2), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qc2qr_autoconv_tend',  tend_out(:,:, 3), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nc_accret_tend',  tend_out(:,:, 4), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nc2nr_autoconv_tend', tend_out(:,:, 5), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nc_selfcollect_tend',  tend_out(:,:, 6), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nr_selfcollect_tend',  tend_out(:,:, 7), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nc_nuceat_tend',  tend_out(:,:, 8), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qccon',  tend_out(:,:, 9), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qcnuc',  tend_out(:,:,10), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qr2qv_evap_tend',  tend_out(:,:,11), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qcevp',  tend_out(:,:,12), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nr_evap_tend',  tend_out(:,:,13), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_ncautr', tend_out(:,:,14), pcols, lchnk) 
!ASD - Ignore For Now   ! ice-phase process rate
!ASD - Ignore For Now   call outfld('P3_qccol',  tend_out(:,:,15), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qwgrth', tend_out(:,:,16), pcols, lchnk) ! Not a tendency in itself, it is used to build qccol and qrcol
!ASD - Ignore For Now   call outfld('P3_qidep',  tend_out(:,:,17), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qrcol',  tend_out(:,:,18), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qinuc',  tend_out(:,:,19), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nc_collect_tend',  tend_out(:,:,20), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nr_collect_tend',  tend_out(:,:,21), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_ni_nucleat_tend',  tend_out(:,:,22), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qi2qv_sublim_tend',  tend_out(:,:,23), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qi2qr_melt_tend',  tend_out(:,:,24), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_ni2nr_melt_tend',  tend_out(:,:,25), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_ni_sublim_tend',  tend_out(:,:,26), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_ni_selfcollect_tend',  tend_out(:,:,27), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qc2qi_hetero_frz_tend', tend_out(:,:,28), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qr2qi_immers_frz_tend', tend_out(:,:,29), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nc2ni_immers_frz_tend', tend_out(:,:,30), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nr2ni_immers_frz_tend', tend_out(:,:,31), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_nr_ice_shed_tend', tend_out(:,:,32), pcols, lchnk) 
!ASD - Ignore For Now   call outfld('P3_qc2qr_ice_shed_tend',  tend_out(:,:,33), pcols, lchnk) 
!ASD - Ignore For Now!   call outfld('P3_qcmul',  tend_out(:,:,34), pcols, lchnk) ! Not actually used, so not actually recorded.  Commented out here for continuity of the array.
!ASD - Ignore For Now   call outfld('P3_ncshdc', tend_out(:,:,35), pcols, lchnk)
!ASD - Ignore For Now   ! sedimentation 
!ASD - Ignore For Now   call outfld('P3_sed_CLDLIQ',  tend_out(:,:,36), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_sed_NUMLIQ',  tend_out(:,:,37), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_sed_CLDRAIN', tend_out(:,:,38), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_sed_NUMRAIN', tend_out(:,:,39), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_sed_CLDICE',  tend_out(:,:,40), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_sed_NUMICE',  tend_out(:,:,41), pcols, lchnk)
!ASD - Ignore For Now   ! microphysics processes 
!ASD - Ignore For Now   call outfld('P3_mtend_CLDLIQ',  tend_out(:,:,42), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_mtend_NUMLIQ',  tend_out(:,:,43), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_mtend_CLDRAIN', tend_out(:,:,44), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_mtend_NUMRAIN', tend_out(:,:,45), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_mtend_CLDICE',  tend_out(:,:,46), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_mtend_NUMICE',  tend_out(:,:,47), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_mtend_Q',       tend_out(:,:,48), pcols, lchnk)
!ASD - Ignore For Now   call outfld('P3_mtend_TH',      tend_out(:,:,49), pcols, lchnk)
!ASD - Ignore For Now   ! Phase change tendencies 
!ASD - Ignore For Now   call outfld('vap_ice_exchange',      vap_ice_exchange,      pcols, lchnk)
!ASD - Ignore For Now   call outfld('vap_liq_exchange',      vap_liq_exchange,      pcols, lchnk)
!ASD - Ignore For Now   call outfld('liq_ice_exchange',      liq_ice_exchange,      pcols, lchnk)
!ASD - Ignore For Now
!ASD - Ignore For Now   call t_stopf('micro_p3_tend_finish')
  end subroutine micro_p3_tend

  !================================================================================================

end module micro_p3_interface
