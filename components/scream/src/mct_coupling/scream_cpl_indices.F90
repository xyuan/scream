module scream_cpl_indices

  use iso_c_binding, only: c_int, c_char

  implicit none
  private

  ! Focus only on the ones that scream imports/exports (subsets of x2a and a2x)
  integer, parameter, public :: num_required_cpl_imports = 29
  integer, parameter, public :: num_scream_imports       = 8
  integer, parameter, public :: num_required_exports     = 33
  integer, parameter, public :: num_optional_cpl_imports = 4
  integer, parameter, public :: num_optional_exports     = 2
  integer, parameter, public :: num_cpl_imports          = num_required_cpl_imports + num_optional_cpl_imports
  integer, parameter, public :: num_exports              = num_required_exports + num_optional_exports

  integer(kind=c_int), public, allocatable, target :: index_x2a(:)
  integer(kind=c_int), public, allocatable, target :: index_a2x(:)

  ! Names used by the component coupler for import/export fields
  character(len=32,kind=c_char), public, allocatable, target :: cpl_names_a2x(:)
  character(len=32,kind=c_char), public, allocatable, target :: cpl_names_x2a(:)

  ! Names used by scream for import/export fields
  character(len=32,kind=c_char), public, allocatable, target :: scr_names_a2x(:)
  character(len=32,kind=c_char), public, allocatable, target :: scr_names_x2a(:)

  ! Vector component of import/export fields. If not a vector field, set to -1.
  integer(kind=c_int), public, allocatable, target :: vec_comp_a2x(:)
  integer(kind=c_int), public, allocatable, target :: vec_comp_x2a(:)



  public :: scream_set_cpl_indices

contains

  subroutine scream_set_cpl_indices (x2a, a2x)
    use iso_c_binding,  only: C_NULL_CHAR
    use mct_mod,        only: mct_aVect, mct_avect_indexra

    !
    ! Input(s)
    !
    type(mct_avect), intent(in) :: x2a, a2x
    !
    ! Local(s)
    !
    integer :: i,idx

    allocate (index_x2a(num_cpl_imports))
    allocate (index_a2x(num_exports))

    allocate (cpl_names_x2a(num_cpl_imports))
    allocate (cpl_names_a2x(num_exports))

    allocate (scr_names_x2a(num_cpl_imports))
    allocate (scr_names_a2x(num_exports))

    allocate (vec_comp_x2a(num_cpl_imports))
    allocate (vec_comp_a2x(num_exports))

    ! Determine attribute vector indices

    ! List of cpl names of inputs that scream cares about

    !------------------------------------------------------------------------------------------
    !Following inputs are surface values so they are dimensioned (1:ncol) for each chunk.
    !"cam_in" derived type is populated using these inputs in atm_import_export.F90
    !using values from other model components
    !"cam_in" is then used by SCREAM model
    !------------------------------------------------------------------------------------------

    !Comments after the inputs below are organized as follows:
    !Long name [units] (cam_in member which captures the input) [List of parameterizations which are using this input currently]

    cpl_names_x2a(1)  = 'Sx_avsdr'  ! short wave direct albedo    [no units] (cam_in%asdir)[RRTMGP]
    cpl_names_x2a(2)  = 'Sx_anidr'  ! long wave direct albedo     [no units] (cam_in%aldir)[RRTMGP]
    cpl_names_x2a(3)  = 'Sx_avsdf'  ! short wave difuse albedo    [no units] (cam_in%asdif)[RRTMGP]
    cpl_names_x2a(4)  = 'Sx_anidf'  ! long wave difuse albedo     [no units] (cam_in%aldif)[RRTMGP]
    cpl_names_x2a(5)  = 'Sx_t'      ! Surface temperature         [K]        (cam_in%ts)   [check_energy/output- not used anywhere else]
    cpl_names_x2a(6)  = 'So_t'      ! sea surface temperature
    cpl_names_x2a(7)  = 'Sl_snowh'  ! Water equivalent snow depth [m]        (cam_in%snowhland) [SHOC]
    cpl_names_x2a(8)  = 'Si_snowh'  ! Snow depth over ice         [m]        (cam_in%snowhice)  [***UNUSED***]

    cpl_names_x2a(9)  = 'Sl_fv'     ! friction velocity
    cpl_names_x2a(10) = 'Sl_ram1'   ! aerodynamical resistance

    cpl_names_x2a(11) = 'Sx_tref'   ! Reference height temperature[K]        (cam_in%tref)      [***UNUSED***]
    cpl_names_x2a(12) = 'Sx_qref'   ! Reference height humidity   [kg/kg]    (cam_in%qref)      [***UNUSED***]
    cpl_names_x2a(13) = 'Sf_ifrac'  ! Fraction of sfc area covered by sea-ice [no units] (cam_in%icefrac) [RRTMGP]
    !NOTE: Sf_ofrac (or ocean frac) is being used by aqua_planet and old schemes like vertical_diffusion,
    !Park stratiform_tend and macrophysics
    cpl_names_x2a(14) = 'Sf_ofrac'  ! Fraction of sfc area covered by ocean   [no units] (cam_in%ocnfrac) [***UNUSED***]
    cpl_names_x2a(15) = 'Sf_lfrac'  ! Fraction of sfc area covered by land    [no units] (cam_in%landfrac)[SHOC/RRTMGP/ZM]
    cpl_names_x2a(16) = 'Sx_u10'    ! 10m wind speed              [m/s]      (cam_in%u10)       [***UNUSED***]
    cpl_names_x2a(17) = 'Faxx_taux' ! Surface stress in X         [N/m2] (cam_in%wsx)  [SHOC]
    cpl_names_x2a(18) = 'Faxx_tauy' ! Surface stress in Y         [N/m2] (cam_in%wsx)  [SHOC]
    cpl_names_x2a(19) = 'Faxx_lat'  ! Surface latent heat flux    [W/m2] (cam_in%lhf)  [energy fixer qqflx_fixer/qneg4]
    cpl_names_x2a(20) = 'Faxx_sen'  ! Surface sensible heat flux  [W/m2] (cam_in%shf)  [SHOC/check_energy_chng]
    cpl_names_x2a(21) = 'Faxx_lwup' ! long wave up radiation flux [W/m2] (cam_in%lwup) [RRTMGP]
    cpl_names_x2a(22) = 'Faxx_evap' ! Surface water vapor flux    [kg/kg](cam_in%cflx(:,1)) [SHOC/check_energy_chng]
    !NOTE:SHOC computes So_ustar (or ustar) internally
    cpl_names_x2a(23) = 'So_ustar'  ! Friction/shear velocity     [m/s]      (cam_in%ustar) [***UNUSED***]
    cpl_names_x2a(24) = 'So_re'     ! ???? (cam_in%re) [***UNUSED***]

    cpl_names_x2a(25) = 'So_ssq' ! surface saturation specific humidity in ocean
    cpl_names_x2a(26) = 'Fall_flxdst1' ! dust flux size bin 1
    cpl_names_x2a(27) = 'Fall_flxdst2' ! dust flux size bin 2
    cpl_names_x2a(28) = 'Fall_flxdst3' ! dust flux size bin 3
    cpl_names_x2a(29) = 'Fall_flxdst4' ! dust flux size bin 4

    ! The following fields are optional
    cpl_names_x2a(30) = 'Sl_soilw' ! volumetric soil water
    cpl_names_x2a(31) = 'Fall_fco2_lnd' ! co2 flux from land
    cpl_names_x2a(32) = 'Faoo_fco2_ocn' ! co2 flux from ocean
    cpl_names_x2a(33) = 'Faoo_fdms_ocn' ! dms flux from ocean

    ! Names used by scream for the input fields above.
    ! Unused fields are marked and skipped during surface coupling.
    scr_names_x2a(1)  = 'sfc_alb_dir_vis'
    scr_names_x2a(2)  = 'sfc_alb_dir_nir'
    scr_names_x2a(3)  = 'sfc_alb_dif_vis'
    scr_names_x2a(4)  = 'sfc_alb_dif_nir'
    scr_names_x2a(5)  = 'unused'
    scr_names_x2a(6)  = 'unused'
    scr_names_x2a(7)  = 'unused'
    scr_names_x2a(8)  = 'unused'
    scr_names_x2a(9)  = 'unused'
    scr_names_x2a(10) = 'unused'
    scr_names_x2a(11) = 'unused'
    scr_names_x2a(12) = 'unused'
    scr_names_x2a(13) = 'unused'
    scr_names_x2a(14) = 'unused'
    scr_names_x2a(15) = 'unused'
    scr_names_x2a(16) = 'unused'
    scr_names_x2a(17) = 'surf_mom_flux'
    scr_names_x2a(18) = 'surf_mom_flux'
    scr_names_x2a(19) = 'unused'
    scr_names_x2a(20) = 'surf_sens_flux'
    scr_names_x2a(21) = 'unused'
    scr_names_x2a(22) = 'surf_latent_flux'
    scr_names_x2a(23) = 'unused'
    scr_names_x2a(24) = 'unused'
    scr_names_x2a(25) = 'unused'
    scr_names_x2a(26) = 'unused'
    scr_names_x2a(27) = 'unused'
    scr_names_x2a(28) = 'unused'
    scr_names_x2a(29) = 'unused'
    scr_names_x2a(30) = 'unused'
    scr_names_x2a(31) = 'unused'
    scr_names_x2a(32) = 'unused'
    scr_names_x2a(33) = 'unused'

    ! Default import vector components to -1. Set surf_mom_flux components.
    do i=1,num_cpl_imports
      vec_comp_x2a(i) = -1
    enddo
    vec_comp_x2a(17) = 0
    vec_comp_x2a(18) = 1

    index_x2a(1) = mct_avect_indexra(x2a,'Sx_avsdr')
    index_x2a(2) = mct_avect_indexra(x2a,'Sx_anidr')
    index_x2a(3) = mct_avect_indexra(x2a,'Sx_avsdf')
    index_x2a(4) = mct_avect_indexra(x2a,'Sx_anidf')
    index_x2a(5) = mct_avect_indexra(x2a,'Sx_t')
    index_x2a(6) = mct_avect_indexra(x2a,'So_t')
    index_x2a(7) = mct_avect_indexra(x2a,'Sl_snowh')
    index_x2a(8) = mct_avect_indexra(x2a,'Si_snowh')
    index_x2a(9) = mct_avect_indexra(x2a,'Sl_fv')
    index_x2a(10) = mct_avect_indexra(x2a,'Sl_ram1')
    index_x2a(11) = mct_avect_indexra(x2a,'Sx_tref')
    index_x2a(12) = mct_avect_indexra(x2a,'Sx_qref')
    index_x2a(13) = mct_avect_indexra(x2a,'Sf_ifrac')
    index_x2a(14) = mct_avect_indexra(x2a,'Sf_ofrac')
    index_x2a(15) = mct_avect_indexra(x2a,'Sf_lfrac')
    index_x2a(16) = mct_avect_indexra(x2a,'Sx_u10')
    index_x2a(17) = mct_avect_indexra(x2a,'Faxx_taux')
    index_x2a(18) = mct_avect_indexra(x2a,'Faxx_tauy')
    index_x2a(19) = mct_avect_indexra(x2a,'Faxx_lat')
    index_x2a(20) = mct_avect_indexra(x2a,'Faxx_sen')
    index_x2a(21) = mct_avect_indexra(x2a,'Faxx_lwup')
    index_x2a(22) = mct_avect_indexra(x2a,'Faxx_evap')
    index_x2a(23) = mct_avect_indexra(x2a,'So_ustar')
    index_x2a(24) = mct_avect_indexra(x2a,'So_re')
    index_x2a(25) = mct_avect_indexra(x2a,'So_ssq')
    index_x2a(26) = mct_avect_indexra(x2a,'Fall_flxdst1')
    index_x2a(27) = mct_avect_indexra(x2a,'Fall_flxdst2')
    index_x2a(28) = mct_avect_indexra(x2a,'Fall_flxdst3')
    index_x2a(29) = mct_avect_indexra(x2a,'Fall_flxdst4')
    index_x2a(30) = mct_avect_indexra(x2a,'Sl_soilw',perrWith='quiet')
    index_x2a(31) = mct_avect_indexra(x2a,'Fall_fco2_lnd',perrWith='quiet')
    index_x2a(32) = mct_avect_indexra(x2a,'Faoo_fco2_ocn',perrWith='quiet')
    index_x2a(33) = mct_avect_indexra(x2a,'Faoo_fdms_ocn',perrWith='quiet')

    ! List of cpl names of outputs that scream needs to pass back to cpl

    !Following outputs are computed in control/camsrfexch.F90 using internal SCREAM variables

    !comments after the outputs are organized as follows:
    !Long name [units] (cam_out derived type member) [optional comment about how it is computed]

    cpl_names_a2x(1)  = 'Sa_z'        ! Geopotential height above surface at midpoints [m] (cam_out%zbot)
    cpl_names_a2x(2)  = 'Sa_u'        ! Zonal wind        [m/s]  (cam_out%ubot)
    cpl_names_a2x(3)  = 'Sa_v'        ! Meridional wind   [m/s]  (cam_out%vbot)
    cpl_names_a2x(4)  = 'Sa_tbot'     ! Lowest model level temperature [K] (cam_out%tbot)
    cpl_names_a2x(5)  = 'Sa_ptem'     ! Potential temperature          [K] (cam_out%thbot)[Computed from temperature and exner function]
    cpl_names_a2x(6)  = 'Sa_pbot'     ! midpoint pressure [Pa]   (cam_out%pbot)
    cpl_names_a2x(7)  = 'Sa_pslv'     ! sea level atm pressure
    cpl_names_a2x(8)  = 'Sa_shum'     ! Specific humidity [kg/kg](cam_out%qbot(i,1)[surface water vapor, i.e., state%q(1:ncol,pver,1)]
    cpl_names_a2x(9)  = 'Sa_dens'     ! Density           [kg/m3](cam_out%rho) [Computed as pbot/(rair*tbot)]
    cpl_names_a2x(10) = 'Faxa_swnet'  ! sw: net
    cpl_names_a2x(11) = 'Faxa_lwdn'   ! downward lw heat flux

    !-------------------------------------------------------------------------------------------------
    !Important notes regarding following 4 cpl_names_a2x variables (for cpl_names_a2x indexed 9 to 12):
    !
    !1. All the prec* variables have the units of m/s in the model but they are converted to mm/s when
    !they are assigned to the respective cam_out members in components/eam/src/cpl/atm_import_export.F90
    !
    !2. Convective precip variables (precc and precsc, definitions below) should be zero for SCREAM since
    !   convection scheme is turned off in SCREAM
    !   'precc'  is Convective precipitation rate (liq + ice)
    !   'precsc' is Convective snow rate (water equivalent)
    !
    !3. Large scale precip is carried in the following variables:
    !   'precl'  is Large-scale (stable) precipitation rate (liq + ice)
    !   'precsl' is Large-scale (stable) snow rate (water equivalent)
    !-------------------------------------------------------------------------------------------------

    !Faxa_rainc is (precc-precsc), therefore it is just the "liquid" part of the convective prec
    !cam_out variable corresponding to "Faxa_rainc" should be zero for SCREAM
    cpl_names_a2x(12) = 'Faxa_rainc'    ! Liquid convective precip  [mm/s] (cam_out%precc-cam_out%precsc) [Obtained from Deep conv.]

    !Faxa_rainl is precl-precsl, therefore it is just the "liquid" part of the large scale prec
    cpl_names_a2x(13) = 'Faxa_rainl'    ! Liquid large-scale precip [mm/s] (cam_out%precl-cam_out%precsl) [obtained from P3]

    !cam_out variable corresponding to "Faxa_snowc" should be zero for SCREAM
    cpl_names_a2x(14) = 'Faxa_snowc'    ! Convective snow rate      [mm/s] (cam_out%precsc) [Obtained from Deep Conv.]
    cpl_names_a2x(15) = 'Faxa_snowl'    ! Large-scale (stable) snow rate [mm/s] (cam_out%precsl) [Obtained from P3]
    cpl_names_x2a(16) = 'Faxa_swndr'    ! sw: nir direct  downward
    cpl_names_x2a(17) = 'Faxa_swvdr'    ! sw: vis direct  downward
    cpl_names_x2a(18) = 'Faxa_swndf'    ! sw: nir diffuse downward
    cpl_names_x2a(19) = 'Faxa_swvdf'    ! sw: vis diffuse downward
    cpl_names_x2a(20) = 'Faxa_bcphidry' ! flux: Black Carbon hydrophilic dry deposition
    cpl_names_x2a(21) = 'Faxa_bcphodry' ! flux: Black Carbon hydrophobic dry deposition
    cpl_names_x2a(22) = 'Faxa_bcphiwet' ! flux: Black Carbon hydrophilic wet deposition
    cpl_names_x2a(23) = 'Faxa_ocphidry' ! flux: Organic Carbon hydrophilic dry deposition
    cpl_names_x2a(24) = 'Faxa_ocphodry' ! flux: Organic Carbon hydrophobic dry deposition
    cpl_names_x2a(25) = 'Faxa_ocphiwet' ! flux: Organic Carbon hydrophilic dry deposition
    cpl_names_x2a(26) = 'Faxa_dstdry1'  ! flux: Size 1 dust -- dry deposition
    cpl_names_x2a(27) = 'Faxa_dstdry2'  ! flux: Size 2 dust -- dry deposition
    cpl_names_x2a(28) = 'Faxa_dstdry3'  ! flux: Size 3 dust -- dry deposition
    cpl_names_x2a(29) = 'Faxa_dstdry4'  ! flux: Size 4 dust -- dry deposition
    cpl_names_x2a(30) = 'Faxa_dstwet1'  ! flux: Size 1 dust -- wet deposition
    cpl_names_x2a(31) = 'Faxa_dstwet2'  ! flux: Size 2 dust -- wet deposition
    cpl_names_x2a(32) = 'Faxa_dstwet3'  ! flux: Size 3 dust -- wet deposition
    cpl_names_x2a(33) = 'Faxa_dstwet4'  ! flux: Size 4 dust -- wet deposition
    cpl_names_x2a(34) = 'Sa_co2prog'    ! Always 0.0_r8 as it is not computed by SCREAM (prognostic co2 is turned off)
    cpl_names_x2a(35) = 'Sa_co2diag'    ! bottom atm level diagnostic co2

    ! Names used by scream for the output fields above. Some field retain
    ! their cpl_name, which will be combinations of multiple scream fields
    ! (this logic is taken car of in SurfaceCoupling). Others will be set
    ! to 0 (set_zero).
    scr_names_a2x(1)  = 'z_mid'
    scr_names_a2x(2)  = 'horiz_winds'
    scr_names_a2x(3)  = 'horiz_winds'
    scr_names_a2x(4)  = 'T_mid'
    scr_names_a2x(5)  = 'Sa_ptem'
    scr_names_a2x(6)  = 'p_mid'
    scr_names_a2x(7)  = 'set_zero'
    scr_names_a2x(8)  = 'qv'
    scr_names_a2x(9)  = 'Sa_dens'
    scr_names_a2x(10) = 'set_zero'
    scr_names_a2x(11) = 'precip_liq_surf'
    scr_names_a2x(12) = 'set_zero'
    scr_names_a2x(13) = 'set_zero'
    scr_names_a2x(14) = 'set_zero'
    scr_names_a2x(15) = 'set_zero'
    scr_names_a2x(16) = 'set_zero'
    scr_names_a2x(17) = 'set_zero'
    scr_names_a2x(18) = 'set_zero'
    scr_names_a2x(19) = 'set_zero'
    scr_names_a2x(20) = 'set_zero'
    scr_names_a2x(21) = 'set_zero'
    scr_names_a2x(22) = 'set_zero'
    scr_names_a2x(23) = 'set_zero'
    scr_names_a2x(24) = 'set_zero'
    scr_names_a2x(25) = 'set_zero'
    scr_names_a2x(26) = 'set_zero'
    scr_names_a2x(27) = 'set_zero'
    scr_names_a2x(28) = 'set_zero'
    scr_names_a2x(29) = 'set_zero'
    scr_names_a2x(30) = 'set_zero'
    scr_names_a2x(31) = 'set_zero'
    scr_names_a2x(32) = 'set_zero'
    scr_names_a2x(33) = 'set_zero'
    scr_names_a2x(34) = 'set_zero'
    scr_names_a2x(35) = 'set_zero'

    ! Set index
    index_a2x(1) = mct_avect_indexra(a2x,'Sa_z')
    index_a2x(2) = mct_avect_indexra(a2x,'Sa_u')
    index_a2x(3) = mct_avect_indexra(a2x,'Sa_v')
    index_a2x(4) = mct_avect_indexra(a2x,'Sa_tbot')
    index_a2x(5) = mct_avect_indexra(a2x,'Sa_ptem')
    index_a2x(6) = mct_avect_indexra(a2x,'Sa_pbot')
    index_a2x(7) = mct_avect_indexra(a2x,'Sa_pslv')
    index_a2x(8) = mct_avect_indexra(a2x,'Sa_shum')
    index_a2x(9) = mct_avect_indexra(a2x,'Sa_dens')
    index_a2x(10) = mct_avect_indexra(a2x,'Faxa_swnet')
    index_a2x(11) = mct_avect_indexra(a2x,'Faxa_lwdn')
    index_a2x(12) = mct_avect_indexra(a2x,'Faxa_rainc')
    index_a2x(13) = mct_avect_indexra(a2x,'Faxa_rainl')
    index_a2x(14) = mct_avect_indexra(a2x,'Faxa_snowc')
    index_a2x(15) = mct_avect_indexra(a2x,'Faxa_snowl')
    index_a2x(16) = mct_avect_indexra(a2x,'Faxa_swndr')
    index_a2x(17) = mct_avect_indexra(a2x,'Faxa_swvdr')
    index_a2x(18) = mct_avect_indexra(a2x,'Faxa_swndf')
    index_a2x(19) = mct_avect_indexra(a2x,'Faxa_swvdf')
    index_a2x(20) = mct_avect_indexra(a2x,'Faxa_bcphidry')
    index_a2x(21) = mct_avect_indexra(a2x,'Faxa_bcphodry')
    index_a2x(22) = mct_avect_indexra(a2x,'Faxa_bcphiwet')
    index_a2x(23) = mct_avect_indexra(a2x,'Faxa_ocphidry')
    index_a2x(24) = mct_avect_indexra(a2x,'Faxa_ocphodry')
    index_a2x(25) = mct_avect_indexra(a2x,'Faxa_ocphiwet')
    index_a2x(26) = mct_avect_indexra(a2x,'Faxa_dstdry1')
    index_a2x(27) = mct_avect_indexra(a2x,'Faxa_dstdry2')
    index_a2x(28) = mct_avect_indexra(a2x,'Faxa_dstdry3')
    index_a2x(29) = mct_avect_indexra(a2x,'Faxa_dstdry4')
    index_a2x(30) = mct_avect_indexra(a2x,'Faxa_dstwet1')
    index_a2x(31) = mct_avect_indexra(a2x,'Faxa_dstwet2')
    index_a2x(32) = mct_avect_indexra(a2x,'Faxa_dstwet3')
    index_a2x(33) = mct_avect_indexra(a2x,'Faxa_dstwet4')
    index_a2x(34) = mct_avect_indexra(a2x,'Sa_co2prog',perrWith='quiet')
    index_a2x(35) = mct_avect_indexra(a2x,'Sa_co2diag',perrWith='quiet')

    ! Default export vector components to -1. Set horiz_winds components.
    do i=1,num_exports
      vec_comp_a2x(i) = -1
    enddo
    vec_comp_a2x(2) = 0
    vec_comp_a2x(3) = 1

    do i=1,num_cpl_imports
      scr_names_x2a(i) = TRIM(scr_names_x2a(i)) // C_NULL_CHAR
    enddo

    do i=1,num_exports
      scr_names_a2x(i) = TRIM(scr_names_a2x(i)) // C_NULL_CHAR
    enddo

    ! We no longer need the cpl names
    deallocate(cpl_names_a2x)
    deallocate(cpl_names_x2a)

  end subroutine scream_set_cpl_indices

end module scream_cpl_indices
