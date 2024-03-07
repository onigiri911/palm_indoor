!> @file urban_surface_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 2015-2020 Czech Technical University in Prague
! Copyright 2015-2020 Institute of Computer Science of the Czech Academy of Sciences, Prague
! Copyright 1997-2020 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
!
! Current revisions:
! -----------------
!
!
! Former revisions:
! -----------------
! $Id: urban_surface_mod.f90 4783 2020-11-13 13:58:45Z raasch $
! - Default building parameters updated and extended (responsible: S. Rissmann)
! - First and second wall layer initialized with individual building properties rather than with
!   with the same properties (heat capacity and conductivity)
!
! 4780 2020-11-10 11:17:10Z suehring
! Enable 3D data output also with 64-bit precision
!
! 4750 2020-10-16 14:27:48Z suehring
! - bugfix in openmp directive
! - make t_green_h and t_green_v public (required in indoor model)
! 
! 4747 2020-10-16 09:19:57Z pavelkrc
! Fix window absorptivity calculation (correctly account for 2-sided reflection)
! 
! 4738 2020-10-14 08:05:07Z maronga
! Updating building data base (on behalf of Sascha Rissmann)
!
! 4733 2020-10-09 12:24:16Z maronga
! Bugfix in calculation of absored radiation by windows (reflected radiation was not taken into
! account)
!
! 4721 2020-10-02 10:21:52Z suehring
! Remove unused variables from USE statement
!
! 4717 2020-09-30 22:27:40Z pavelkrc
! Fixes and optimizations of OpenMP parallelization, formatting of OpenMP
! directives (J. Resler)
!
! 4713 2020-09-29 12:02:05Z pavelkrc
! - Do not change original fractions in USM energy balance
! - Correct OpenMP parallelization
! Author: J. Resler
!
! 4701 2020-09-27 11:02:15Z maronga
! Corrected parameter 33 for building_type 2 (ground floor window emissivity
!
! 4698 2020-09-25 08:37:55Z maronga
! Bugfix in calculation of iwghf_eb and iwghf_eb_window (introduced in previous revisions)
!
! 4694 2020-09-23 15:09:19Z pavelkrc
! Fix writing and reading of surface data to/from MPI restart file
!
! 4693 2020-09-22 19:47:04Z maronga
! Bugfix for last commit
!
! 4692 2020-09-22 17:17:52Z maronga
! Bugfix for previous revision
!
! 4687 2020-09-21 19:40:16Z
! Optimized code structure for treatment of inner wall and window heat flux
!
! 4671 2020-09-09 20:27:58Z pavelkrc
! Radiative transfer model RTM version 4.1
! - Implementation of downward facing USM and LSM surfaces
! - Restructuralization EB call
! - Improved debug logging
! - Removal of deprecated CSV inputs
! - Bugfixes
! Author: J. Resler (Institute of Computer Science, Prague)
!
! 4669 2020-09-09 13:43:47Z pavelkrc
! Fix calculation of force_radiation_call
!
! 4668 2020-09-09 13:00:16Z pavelkrc
! Limit vertical r_a similarly to horizontal
!
! 4660 2020-09-01 14:49:39Z pavelkrc
! Fix of wrong limitation of calculation of ueff (avoid r_a=inf) (J. Resler)
!
! 4653 2020-08-27 08:54:43Z pavelkrc
! Remove separate direction indices for urban and land faces (part of RTM v4.0)
!
! 4630 2020-07-30 14:54:34Z suehring
! - Bugfix in resistance calculation - avoid potential divisions by zero
! - Minor formatting adjustment
!
! 4602 2020-07-14 14:49:45Z suehring
! Add missing initialization of albedo type with values given from static input file
! 
! 4581 2020-06-29 08:49:58Z suehring
! Missing initialization in case of cyclic_fill runs
! 
! 4535 2020-05-15 12:07:23Z raasch
! bugfix for restart data format query
! 
! 4517 2020-05-03 14:29:30Z raasch
! added restart with MPI-IO for reading local arrays
! 
! 4510 2020-04-29 14:19:18Z raasch
! Further re-formatting to follow the PALM coding standard
!
! 4509 2020-04-26 15:57:55Z raasch
! File re-formatted to follow the PALM coding standard
!
! 4500 2020-04-17 10:12:45Z suehring
! Allocate array for wall heat flux, which is further used to aggregate tile
! fractions in the surface output
!
! 4495 2020-04-13 20:11:20Z raasch
! Restart data handling with MPI-IO added
!
! 4493 2020-04-10 09:49:43Z pavelkrc
! J.Resler, 2020/03/19
! - Remove reading of deprecated input parameters c_surface and lambda_surf
! - And calculate them from parameters of the outer wall/roof layer
!
! 4481 2020-03-31 18:55:54Z maronga
! Use statement for exchange horiz added
!
! 4442 2020-03-04 19:21:13Z suehring
! Change order of dimension in surface arrays %frac, %emissivity and %albedo to allow for better
! vectorization in the radiation interactions.
!
! 4441 2020-03-04 19:20:35Z suehring
! Removed wall_flags_static_0 from USE statements as it's not used within the module
!
! 4329 2019-12-10 15:46:36Z motisi
! Renamed wall_flags_0 to wall_flags_static_0
!
! 4309 2019-11-26 18:49:59Z suehring
! - Bugfix, include m_liq into restarts
! - Remove unused arrays for liquid water and saturation moisture at vertical walls
!
! 4305 2019-11-25 11:15:40Z suehring
! Revision of some indoor-model parameters
!
! 4259 2019-10-09 10:05:22Z suehring
! Instead of terminate the job in case the relative wall fractions do not sum-up to one, give only
! an informative message and normalize the fractions.
!
! 4258 2019-10-07 13:29:08Z suehring
! - Add checks to ensure that relative fractions of walls, windowns and green surfaces sum-up to one.
! - Revise message calls dealing with local checks.
!
! 4245 2019-09-30 08:40:37Z pavelkrc
! Initialize explicit per-surface parameters from building_surface_pars
!
! 4238 2019-09-25 16:06:01Z suehring
! Indoor-model parameters for some building types adjusted in order to avoid unrealistically high
! indoor temperatures (S. Rissmann)
!
! 4230 2019-09-11 13:58:14Z suehring
! Bugfix, initialize canopy resistance. Even if no green fraction is set, r_canopy must be
! initialized for output purposes.
!
! 4227 2019-09-10 18:04:34Z gronemeier
! Implement new palm_date_time_mod
!
! 4214 2019-09-02 15:57:02Z suehring
! Bugfix, missing initialization and clearing of soil-moisture tendency (J.Resler)
!
! 4182 2019-08-22 15:20:23Z scharf
! Corrected 'Former revisions' section
!
! 4168 2019-08-16 13:50:17Z suehring
! Replace function get_topography_top_index by topo_top_ind
!
! 4148 2019-08-08 11:26:00Z suehring
! - Add anthropogenic heat output factors for heating and cooling to building data base
! - Move definition of building_pars to usm_init_arrays since it is already required in the indoor
!   model
!
! 4127 2019-07-30 14:47:10Z suehring
! Do not add anthopogenic energy during wall/soil spin-up (merge from branch resler)
!
! 4077 2019-07-09 13:27:11Z gronemeier
! Set roughness length z0 and z0h/q at ground-floor level to same value as those above ground-floor
! level
!
! 4051 2019-06-24 13:58:30Z suehring
! Remove work-around for green surface fraction on buildings (do not set it zero)
!
! 4050 2019-06-24 13:57:27Z suehring
! In order to avoid confusion with global control parameter, rename the USM-internal flag spinup
! into during_spinup.
!
! 3987 2019-05-22 09:52:13Z kanani
! Introduce alternative switch for debug output during timestepping
!
! 3943 2019-05-02 09:50:41Z maronga
! Removed qsws_eb. Bugfix in calculation of qsws.
!
! 3933 2019-04-25 12:33:20Z kanani
! Remove allocation of pt_2m, this is done in surface_mod now (surfaces%pt_2m)
!
! 3921 2019-04-18 14:21:10Z suehring
! Undo accidentally commented initialization
!
! 3918 2019-04-18 13:33:11Z suehring
! Set green fraction to zero also at vertical surfaces
!
! 3914 2019-04-17 16:02:02Z suehring
! In order to obtain correct surface temperature during spinup set window fraction to zero
! (only during spinup) instead of just disabling time-integration of window-surface temperature.
!
! 3901 2019-04-16 16:17:02Z suehring
! Workaround - set green fraction to zero ( green-heat model crashes ).
!
! 3896 2019-04-15 10:10:17Z suehring
!
!
! 3896 2019-04-15 10:10:17Z suehring
! Bugfix, wrong index used for accessing building_pars from PIDS
!
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction of additional debug
! messages
!
! 3882 2019-04-10 11:08:06Z suehring
! Avoid different type kinds
! Move definition of building-surface properties from declaration block to an extra routine
!
! 3881 2019-04-10 09:31:22Z suehring
! Revise determination of local ground-floor level height.
! Make level 3 initalization conform with Palm-input-data standard
! Move output of albedo and emissivity to radiation module
!
! 3832 2019-03-28 13:16:58Z raasch
! Instrumented with openmp directives
!
! 3824 2019-03-27 15:56:16Z pavelkrc
! Remove unused imports
!
!
! 3814 2019-03-26 08:40:31Z pavelkrc
! Unused subroutine commented out
!
! 3769 2019-02-28 10:16:49Z moh.hefny
! Removed unused variables
!
! 3767 2019-02-27 08:18:02Z raasch
! Unused variables removed from rrd-subroutines parameter list
!
! 3748 2019-02-18 10:38:31Z suehring
! Revise conversion of waste-heat flux (do not divide by air density, will be done in diffusion_s)
!
! 3745 2019-02-15 18:57:56Z suehring
! - Remove internal flag indoor_model (is a global control parameter)
! - Add waste heat from buildings to the kinmatic heat flux
! - Consider waste heat in restart data
! - Remove unused USE statements
!
! 3744 2019-02-15 18:38:58Z suehring
! Fixed surface heat capacity in the building parameters convert the file back to unix format
!
! 3730 2019-02-11 11:26:47Z moh.hefny
! Formatting and clean-up (rvtils)
!
! 3710 2019-01-30 18:11:19Z suehring
! Check if building type is set within a valid range.
!
! 3705 2019-01-29 19:56:39Z suehring
! Make nzb_wall public, required for virtual-measurements
!
! 3704 2019-01-29 19:51:41Z suehring
! Some interface calls moved to module_interface + cleanup
!
! 3655 2019-01-07 16:51:22Z knoop
! Implementation of the PALM module interface
!
! 2007 2016-08-24 15:47:17Z kanani
! Initial revision
!
!
! Description:
! ------------
! 2016/6/9 - Initial version of the USM (Urban Surface Model)
!            authors: Jaroslav Resler, Pavel Krc (Czech Technical University in Prague and Institute
!            of Computer Science of the Czech Academy of Sciences, Prague)
!            with contributions: Michal Belda, Nina Benesova, Ondrej Vlcek
!            partly inspired by PALM LSM (B. Maronga)
!            parameterizations of Ra checked with TUF3D (E. S. Krayenhoff)
!> Module for Urban Surface Model (USM)
!> The module includes:
!>    1. Radiation model with direct/diffuse radiation, shading, reflections and integration with
!>       plant canopy
!>    2. Wall and wall surface model
!>    3. Surface layer energy balance
!>    4. Anthropogenic heat (only from transportation so far)
!>    5. Necessary auxiliary subroutines (reading inputs, writing outputs, restart simulations, ...)
!> It also makes use of standard radiation and integrates it into urban surface model.
!>
!> Further work:
!> -------------
!> @todo Revise sorting of building_pars
!> @todo Revise initialization when building_pars / building_surface_pars are provided - 
!>       intialization is not consistent to building_pars
!> @todo Revise flux conversion in energy-balance solver
!> @todo Check divisions in wtend (etc.) calculations for possible division by zero, e.g. in case
!> fraq(0,m) + fraq(1,m) = 0?!
!> @todo Use unit 90 for OPEN/CLOSE of input files (FK)
!--------------------------------------------------------------------------------------------------!
 MODULE urban_surface_mod

    USE arrays_3d,                                                                                 &
        ONLY:  exner,                                                                              &
               hyp,                                                                                &
               hyrho,                                                                              &
               p,                                                                                  &
               prr,                                                                                &
               pt,                                                                                 &
               q,                                                                                  &
               ql,                                                                                 &
               tend,                                                                               &
               u,                                                                                  &
               v,                                                                                  &
               vpt,                                                                                &
               w,                                                                                  &
               zu

    USE calc_mean_profile_mod,                                                                     &
        ONLY:  calc_mean_profile

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  c_p,                                                                                &
               degc_to_k,                                                                          &
               g,                                                                                  &
               kappa,                                                                              &
               l_v,                                                                                &
               magnus_tl,                                                                       &
               pi,                                                                                 &
               r_d,                                                                                &
               rho_l,                                                                              &
               sigma_sb

    USE control_parameters,                                                                        &
        ONLY:  average_count_3d,                                                                   &
               coupling_char,                                                                      &
               coupling_start_time,                                                                &
               debug_output,                                                                       &
               debug_output_timestep,                                                              &
               debug_string,                                                                       &
               dt_do3d,                                                                            &
               dt_3d,                                                                              &
               dz,                                                                                 &
               end_time,                                                                           &
               humidity,                                                                           &
               indoor_model,                                                                       &
               initializing_actions,                                                               &
               intermediate_timestep_count,                                                        &
               intermediate_timestep_count_max,                                                    &
               io_blocks,                                                                          &
               io_group,                                                                           &
               large_scale_forcing,                                                                &
               lsf_surf,                                                                           &
               message_string,                                                                     &
               pt_surface,                                                                         &
               restart_data_format_output,                                                         &
               surface_pressure,                                                                   &
               time_since_reference_point,                                                         &
               timestep_scheme,                                                                    &
               topography,                                                                         &
               tsc,                                                                                &
               urban_surface,                                                                      &
               varnamelength


    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bulk_cloud_model,                                                                   &
               precipitation

    USE cpulog,                                                                                    &
        ONLY:  cpu_log,                                                                            &
               log_point,                                                                          &
               log_point_s

    USE grid_variables,                                                                            &
        ONLY:  ddx,                                                                                &
               ddx2,                                                                               &
               ddy,                                                                                &
               ddy2,                                                                               &
               dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nnx,                                                                                &
               nny,                                                                                &
               nnz,                                                                                &
               nx,                                                                                 &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               ny,                                                                                 &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_top_ind

    USE, INTRINSIC :: iso_c_binding

    USE kinds

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time,                                                                      &
               seconds_per_hour

    USE pegrid

    USE radiation_model_mod,                                                                       &
        ONLY:  albedo_type,                                                                        &
               dirname,                                                                            &
               diridx,                                                                             &
               dirint,                                                                             &
               force_radiation_call,                                                               &
               id,                                                                                 &
               idown,                                                                              &
               ieast,                                                                              &
               inorth,                                                                             &
               isouth,                                                                             &
               iup,                                                                                &
               iwest,                                                                              &
               nd,                                                                                 &
               nz_urban_b,                                                                         &
               nz_urban_t,                                                                         &
               radiation_interaction,                                                              &
               radiation,                                                                          &
               rad_lw_in,                                                                          &
               rad_lw_out,                                                                         &
               rad_sw_in,                                                                          &
               rad_sw_out,                                                                         &
               unscheduled_radiation_calls

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_surface_filetypes,                                                        &
               rrd_mpi_io,                                                                         &
               rrd_mpi_io_surface,                                                                 &
               wrd_mpi_io,                                                                         &
               wrd_mpi_io_surface

    USE statistics,                                                                                &
        ONLY:  hom,                                                                                &
               statistic_regions

    USE surface_mod,                                                                               &
        ONLY:  ind_pav_green,                                                                      &
               ind_veg_wall,                                                                       &
               ind_wat_win,                                                                        &
               surf_type,                                                                          &
               surf_usm_h,                                                                         &
               surf_usm_v,                                                                         &
               surface_restore_elements


    IMPLICIT NONE

!
!-- USM model constants

    REAL(wp), PARAMETER ::  b_ch               = 6.04_wp    !< Clapp & Hornberger exponent
    REAL(wp), PARAMETER ::  lambda_h_green_dry = 0.19_wp    !< heat conductivity for dry soil
    REAL(wp), PARAMETER ::  lambda_h_green_sm  = 3.44_wp    !< heat conductivity of the soil matrix
    REAL(wp), PARAMETER ::  lambda_h_water     = 0.57_wp    !< heat conductivity of water
    REAL(wp), PARAMETER ::  psi_sat            = -0.388_wp  !< soil matrix potential at saturation
    REAL(wp), PARAMETER ::  rho_c_soil         = 2.19E6_wp  !< volumetric heat capacity of soil
    REAL(wp), PARAMETER ::  rho_c_water        = 4.20E6_wp  !< volumetric heat capacity of water
!    REAL(wp), PARAMETER ::  m_max_depth        = 0.0002_wp  !< Maximum capacity of the water reservoir (m)

!
!-- Soil parameters I           alpha_vg,      l_vg_green,    n_vg, gamma_w_green_sat
    REAL(wp), DIMENSION(0:3,1:7), PARAMETER ::  soil_pars = RESHAPE( (/     &
                                 3.83_wp,  1.250_wp, 1.38_wp,  6.94E-6_wp, &  !< soil 1
                                 3.14_wp, -2.342_wp, 1.28_wp,  1.16E-6_wp, &  !< soil 2
                                 0.83_wp, -0.588_wp, 1.25_wp,  0.26E-6_wp, &  !< soil 3
                                 3.67_wp, -1.977_wp, 1.10_wp,  2.87E-6_wp, &  !< soil 4
                                 2.65_wp,  2.500_wp, 1.10_wp,  1.74E-6_wp, &  !< soil 5
                                 1.30_wp,  0.400_wp, 1.20_wp,  0.93E-6_wp, &  !< soil 6
                                 0.00_wp,  0.00_wp,  0.00_wp,  0.57E-6_wp  &  !< soil 7
                                 /), (/ 4, 7 /) )

!
!-- Soil parameters II              swc_sat,     fc,   wilt,    swc_res
    REAL(wp), DIMENSION(0:3,1:7), PARAMETER ::  m_soil_pars = RESHAPE( (/ &
                                 0.403_wp, 0.244_wp, 0.059_wp, 0.025_wp, &  !< soil 1
                                 0.439_wp, 0.347_wp, 0.151_wp, 0.010_wp, &  !< soil 2
                                 0.430_wp, 0.383_wp, 0.133_wp, 0.010_wp, &  !< soil 3
                                 0.520_wp, 0.448_wp, 0.279_wp, 0.010_wp, &  !< soil 4
                                 0.614_wp, 0.541_wp, 0.335_wp, 0.010_wp, &  !< soil 5
                                 0.766_wp, 0.663_wp, 0.267_wp, 0.010_wp, &  !< soil 6
                                 0.472_wp, 0.323_wp, 0.171_wp, 0.000_wp  &  !< soil 7
                                 /), (/ 4, 7 /) )
!
!-- Value 9999999.9_wp -> Generic available or user-defined value must be set otherwise
!-- -> No generic variable and user setting is optional
    REAL(wp) ::  alpha_vangenuchten = 9999999.9_wp      !< NAMELIST alpha_vg
    REAL(wp) ::  field_capacity = 9999999.9_wp          !< NAMELIST fc
    REAL(wp) ::  hydraulic_conductivity = 9999999.9_wp  !< NAMELIST gamma_w_green_sat
    REAL(wp) ::  l_vangenuchten = 9999999.9_wp          !< NAMELIST l_vg
    REAL(wp) ::  n_vangenuchten = 9999999.9_wp          !< NAMELIST n_vg
    REAL(wp) ::  residual_moisture = 9999999.9_wp       !< NAMELIST m_res
    REAL(wp) ::  saturation_moisture = 9999999.9_wp     !< NAMELIST m_sat
    REAL(wp) ::  wilting_point = 9999999.9_wp           !< NAMELIST m_wilt

!
!-- Configuration parameters (they can be setup in PALM config)
    LOGICAL ::  force_radiation_call_l = .FALSE.   !< flag parameter for unscheduled radiation model calls
    LOGICAL ::  usm_wall_mod = .FALSE.             !< reduces conductivity of the first 2 wall layers by factor 0.1


    INTEGER(iwp) ::  building_type = 1               !< default building type (preleminary setting)
    INTEGER(iwp) ::  roof_category = 2               !< default category for root surface
    INTEGER(iwp) ::  wall_category = 2               !< default category for wall surface over pedestrian zone

    REAL(wp)     ::  d_roughness_concrete            !< inverse roughness length of average concrete surface
    REAL(wp)     ::  roughness_concrete = 0.001_wp   !< roughness length of average concrete surface

!
!-- Indices of input attributes in building_pars for (above) ground floor level
    INTEGER(iwp) ::  ind_alb_wall_agfl     = 38   !< index in input list for albedo_type of wall above ground floor level
    INTEGER(iwp) ::  ind_alb_wall_gfl      = 66   !< index in input list for albedo_type of wall ground floor level
    INTEGER(iwp) ::  ind_alb_wall_r        = 101  !< index in input list for albedo_type of wall roof
    INTEGER(iwp) ::  ind_alb_green_agfl    = 39   !< index in input list for albedo_type of green above ground floor level
    INTEGER(iwp) ::  ind_alb_green_gfl     = 78   !< index in input list for albedo_type of green ground floor level
    INTEGER(iwp) ::  ind_alb_green_r       = 117  !< index in input list for albedo_type of green roof
    INTEGER(iwp) ::  ind_alb_win_agfl      = 40   !< index in input list for albedo_type of window fraction above ground floor
                                                  !< level
    INTEGER(iwp) ::  ind_alb_win_gfl       = 77   !< index in input list for albedo_type of window fraction ground floor level
    INTEGER(iwp) ::  ind_alb_win_r         = 115  !< index in input list for albedo_type of window fraction roof
    INTEGER(iwp) ::  ind_emis_wall_agfl    = 14   !< index in input list for wall emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_wall_gfl     = 32   !< index in input list for wall emissivity, ground floor level
    INTEGER(iwp) ::  ind_emis_wall_r       = 100  !< index in input list for wall emissivity, roof
    INTEGER(iwp) ::  ind_emis_green_agfl   = 15   !< index in input list for green emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_green_gfl    = 34   !< index in input list for green emissivity, ground floor level
    INTEGER(iwp) ::  ind_emis_green_r      = 116  !< index in input list for green emissivity, roof
    INTEGER(iwp) ::  ind_emis_win_agfl     = 16   !< index in input list for window emissivity, above ground floor level
    INTEGER(iwp) ::  ind_emis_win_gfl      = 33   !< index in input list for window emissivity, ground floor level
    INTEGER(iwp) ::  ind_emis_win_r        = 113  !< index in input list for window emissivity, roof
    INTEGER(iwp) ::  ind_gflh              = 20   !< index in input list for ground floor level height
    INTEGER(iwp) ::  ind_green_frac_w_agfl = 2    !< index in input list for green fraction on wall, above ground floor level
    INTEGER(iwp) ::  ind_green_frac_w_gfl  = 23   !< index in input list for green fraction on wall, ground floor level
    INTEGER(iwp) ::  ind_green_frac_r_agfl = 3    !< index in input list for green fraction on roof, above ground floor level
    INTEGER(iwp) ::  ind_green_frac_r_gfl  = 24   !< index in input list for green fraction on roof, ground floor level
    INTEGER(iwp) ::  ind_green_type_roof   = 118  !< index in input list for type of green roof
    INTEGER(iwp) ::  ind_hc1_agfl          = 6    !< index in input list for heat capacity at first wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc1_gfl           = 26   !< index in input list for heat capacity at first wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc1_wall_r        = 94   !< index in input list for heat capacity at first wall layer, roof
    INTEGER(iwp) ::  ind_hc1_win_agfl      = 83   !< index in input list for heat capacity at first window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc1_win_gfl       = 71   !< index in input list for heat capacity at first window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_hc1_win_r         = 107  !< index in input list for heat capacity at first window layer, roof
    INTEGER(iwp) ::  ind_hc2_agfl          = 7    !< index in input list for heat capacity at second wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc2_gfl           = 27   !< index in input list for heat capacity at second wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc2_wall_r        = 95   !< index in input list for heat capacity at second wall layer, roof
    INTEGER(iwp) ::  ind_hc2_win_agfl      = 84   !< index in input list for heat capacity at second window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc2_win_gfl       = 72   !< index in input list for heat capacity at second window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_hc2_win_r         = 108  !< index in input list for heat capacity at second window layer, roof
    INTEGER(iwp) ::  ind_hc3_agfl          = 8    !< index in input list for heat capacity at third wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc3_gfl           = 28   !< index in input list for heat capacity at third wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc3_wall_r        = 96   !< index in input list for heat capacity at third wall layer, roof
    INTEGER(iwp) ::  ind_hc3_win_agfl      = 85   !< index in input list for heat capacity at third window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc3_win_gfl       = 73   !< index in input list for heat capacity at third window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_hc3_win_r         = 109  !< index in input list for heat capacity at third window layer, roof
    INTEGER(iwp) ::  ind_hc4_agfl          = 136  !< index in input list for heat capacity at fourth wall layer, above ground floor level
    INTEGER(iwp) ::  ind_hc4_gfl           = 138  !< index in input list for heat capacity at fourth wall layer, ground floor level
    INTEGER(iwp) ::  ind_hc4_wall_r        = 146  !< index in input list for heat capacity at fourth wall layer, roof
    INTEGER(iwp) ::  ind_hc4_win_agfl      = 144  !< index in input list for heat capacity at fourth window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_hc4_win_gfl       = 142  !< index in input list for heat capacity at fourth window layer, ground floor level
    INTEGER(iwp) ::  ind_hc4_win_r         = 148  !< index in input list for heat capacity at fourth window layer, roof
    INTEGER(iwp) ::  ind_indoor_target_temp_summer = 12  !<
    INTEGER(iwp) ::  ind_indoor_target_temp_winter = 13  !<
    INTEGER(iwp) ::  ind_lai_r_agfl        = 4    !< index in input list for LAI on roof, above ground floor level
    INTEGER(iwp) ::  ind_lai_r_gfl         = 4    !< index in input list for LAI on roof, ground floor level
    INTEGER(iwp) ::  ind_lai_w_agfl        = 5    !< index in input list for LAI on wall, above ground floor level
    INTEGER(iwp) ::  ind_lai_w_gfl         = 25   !< index in input list for LAI on wall, ground floor level
    INTEGER(iwp) ::  ind_tc1_agfl          = 9    !< index in input list for thermal conductivity at first wall layer, above ground floor level
    INTEGER(iwp) ::  ind_tc1_gfl           = 29   !< index in input list for thermal conductivity at first wall layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc1_wall_r        = 97   !< index in input list for thermal conductivity at first wall layer, roof
    INTEGER(iwp) ::  ind_tc1_win_agfl      = 86   !< index in input list for thermal conductivity at first window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc1_win_gfl       = 74   !< index in input list for thermal conductivity at first window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc1_win_r         = 110  !< index in input list for thermal conductivity at first window layer, roof
    INTEGER(iwp) ::  ind_tc2_agfl          = 10   !< index in input list for thermal conductivity at second wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc2_gfl           = 30   !< index in input list for thermal conductivity at second wall layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc2_wall_r        = 98   !< index in input list for thermal conductivity at second wall layer, roof
    INTEGER(iwp) ::  ind_tc2_win_agfl      = 87   !< index in input list for thermal conductivity at second window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc2_win_gfl       = 75   !< index in input list for thermal conductivity at second window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc2_win_r         = 111  !< index in input list for thermal conductivity at second window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc3_agfl          = 11   !< index in input list for thermal conductivity at third wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc3_gfl           = 31   !< index in input list for thermal conductivity at third wall layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc3_wall_r        = 99   !< index in input list for thermal conductivity at third wall layer, roof
    INTEGER(iwp) ::  ind_tc3_win_agfl      = 88   !< index in input list for thermal conductivity at third window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc3_win_gfl       = 76   !< index in input list for thermal conductivity at third window layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc3_win_r         = 112  !< index in input list for thermal conductivity at third window layer, roof
    INTEGER(iwp) ::  ind_tc4_agfl          = 137  !< index in input list for thermal conductivity at fourth wall layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc4_gfl           = 139  !< index in input list for thermal conductivity at fourth wall layer,
                                                  !< ground floor level
    INTEGER(iwp) ::  ind_tc4_wall_r        = 147  !< index in input list for thermal conductivity at fourth wall layer, roof
    INTEGER(iwp) ::  ind_tc4_win_agfl      = 145  !< index in input list for thermal conductivity at fourth window layer,
                                                  !< above ground floor level
    INTEGER(iwp) ::  ind_tc4_win_gfl       = 143  !< index in input list for thermal conductivity at first window layer, ground floor level
    INTEGER(iwp) ::  ind_tc4_win_r         = 149  !< index in input list for thermal conductivity at third window layer, roof
    INTEGER(iwp) ::  ind_thick_1_agfl      = 41   !< index for wall layer thickness - 1st layer above ground floor level
    INTEGER(iwp) ::  ind_thick_1_gfl       = 62   !< index for wall layer thickness - 1st layer ground floor level
    INTEGER(iwp) ::  ind_thick_1_wall_r    = 90   !< index for wall layer thickness - 1st layer roof
    INTEGER(iwp) ::  ind_thick_1_win_agfl  = 79   !< index for window layer thickness - 1st layer above ground floor level
    INTEGER(iwp) ::  ind_thick_1_win_gfl   = 67   !< index for window layer thickness - 1st layer ground floor level
    INTEGER(iwp) ::  ind_thick_1_win_r     = 103  !< index for window layer thickness - 1st layer roof
    INTEGER(iwp) ::  ind_thick_2_agfl      = 42   !< index for wall layer thickness - 2nd layer above ground floor level
    INTEGER(iwp) ::  ind_thick_2_gfl       = 63   !< index for wall layer thickness - 2nd layer ground floor level
    INTEGER(iwp) ::  ind_thick_2_wall_r    = 91   !< index for wall layer thickness - 2nd layer roof
    INTEGER(iwp) ::  ind_thick_2_win_agfl  = 80   !< index for window layer thickness - 2nd layer above ground floor level
    INTEGER(iwp) ::  ind_thick_2_win_gfl   = 68   !< index for window layer thickness - 2nd layer ground floor level
    INTEGER(iwp) ::  ind_thick_2_win_r     = 104  !< index for window layer thickness - 2nd layer roof
    INTEGER(iwp) ::  ind_thick_3_agfl      = 43   !< index for wall layer thickness - 3rd layer above ground floor level
    INTEGER(iwp) ::  ind_thick_3_gfl       = 64   !< index for wall layer thickness - 3rd layer ground floor level
    INTEGER(iwp) ::  ind_thick_3_wall_r    = 92   !< index for wall layer thickness - 3rd layer roof
    INTEGER(iwp) ::  ind_thick_3_win_agfl  = 81   !< index for window layer thickness - 3rd layer above ground floor level
    INTEGER(iwp) ::  ind_thick_3_win_gfl   = 69   !< index for window layer thickness - 3rd layer ground floor level
    INTEGER(iwp) ::  ind_thick_3_win_r     = 105  !< index for window layer thickness - 3rd layer roof
    INTEGER(iwp) ::  ind_thick_4_agfl      = 44   !< index for wall layer thickness - 4th layer above ground floor level
    INTEGER(iwp) ::  ind_thick_4_gfl       = 65   !< index for wall layer thickness - 4th layer ground floor level
    INTEGER(iwp) ::  ind_thick_4_wall_r    = 93   !< index for wall layer thickness - 4st layer roof
    INTEGER(iwp) ::  ind_thick_4_win_agfl  = 82   !< index for window layer thickness - 4th layer above ground floor level
    INTEGER(iwp) ::  ind_thick_4_win_gfl   = 70   !< index for window layer thickness - 4th layer ground floor level
    INTEGER(iwp) ::  ind_thick_4_win_r     = 106  !< index for window layer thickness - 4th layer roof
    INTEGER(iwp) ::  ind_trans_agfl        = 17   !< index in input list for window transmissivity, above ground floor level
    INTEGER(iwp) ::  ind_trans_gfl         = 35   !< index in input list for window transmissivity, ground floor level
    INTEGER(iwp) ::  ind_trans_r           = 114  !< index in input list for window transmissivity, roof
    INTEGER(iwp) ::  ind_wall_frac_agfl    = 0    !< index in input list for wall fraction, above ground floor level
    INTEGER(iwp) ::  ind_wall_frac_gfl     = 21   !< index in input list for wall fraction, ground floor level
    INTEGER(iwp) ::  ind_wall_frac_r       = 89   !< index in input list for wall fraction, roof
    INTEGER(iwp) ::  ind_win_frac_agfl     = 1    !< index in input list for window fraction, above ground floor level
    INTEGER(iwp) ::  ind_win_frac_gfl      = 22   !< index in input list for window fraction, ground floor level
    INTEGER(iwp) ::  ind_win_frac_r        = 102  !< index in input list for window fraction, roof
    INTEGER(iwp) ::  ind_z0_agfl           = 18   !< index in input list for z0, above ground floor level
    INTEGER(iwp) ::  ind_z0_gfl            = 36   !< index in input list for z0, ground floor level
    INTEGER(iwp) ::  ind_z0qh_agfl         = 19   !< index in input list for z0h / z0q, above ground floor level
    INTEGER(iwp) ::  ind_z0qh_gfl          = 37   !< index in input list for z0h / z0q, ground floor level
!
!-- Indices of input attributes in building_surface_pars (except for radiation-related, which are in
!-- radiation_model_mod)
    CHARACTER(37), DIMENSION(0:7), PARAMETER ::  building_type_name = (/     &
                                   'user-defined                         ', &  !< type 0
                                   'residential - 1950                   ', &  !< type  1
                                   'residential 1951 - 2000              ', &  !< type  2
                                   'residential 2001 -                   ', &  !< type  3
                                   'office - 1950                        ', &  !< type  4
                                   'office 1951 - 2000                   ', &  !< type  5
                                   'office 2001 -                        ', &  !< type  6
                                   'bridges                              '  &  !< type  7
                                                                     /)

    INTEGER(iwp) ::  ind_s_emis_green                = 14  !< index for emissivity of green fraction (0-1)
    INTEGER(iwp) ::  ind_s_emis_wall                 = 13  !< index for emissivity of wall fraction (0-1)
    INTEGER(iwp) ::  ind_s_emis_win                  = 15  !< index for emissivity o f window fraction (0-1)
    INTEGER(iwp) ::  ind_s_green_frac_r              = 3   !< index for green fraction on roof (0-1)
    INTEGER(iwp) ::  ind_s_green_frac_w              = 2   !< index for green fraction on wall (0-1)
    INTEGER(iwp) ::  ind_s_hc1                       = 5   !< index for heat capacity of wall layer 1
    INTEGER(iwp) ::  ind_s_hc2                       = 6   !< index for heat capacity of wall layer 2
    INTEGER(iwp) ::  ind_s_hc3                       = 7   !< index for heat capacity of wall layer 3
    INTEGER(iwp) ::  ind_s_indoor_target_temp_summer = 11  !< index for indoor target summer temperature
    INTEGER(iwp) ::  ind_s_indoor_target_temp_winter = 12  !< index for indoor target winter temperature
    INTEGER(iwp) ::  ind_s_lai_r                     = 4   !< index for leaf area index of green fraction
    INTEGER(iwp) ::  ind_s_tc1                       = 8   !< index for thermal conducivity of wall layer 1
    INTEGER(iwp) ::  ind_s_tc2                       = 9   !< index for thermal conducivity of wall layer 2
    INTEGER(iwp) ::  ind_s_tc3                       = 10  !< index for thermal conducivity of wall layer 3
    INTEGER(iwp) ::  ind_s_trans                     = 16  !< index for transmissivity of window fraction (0-1)
    INTEGER(iwp) ::  ind_s_wall_frac                 = 0   !< index for wall fraction (0-1)
    INTEGER(iwp) ::  ind_s_win_frac                  = 1   !< index for window fraction (0-1)
    INTEGER(iwp) ::  ind_s_z0                        = 17  !< index for roughness length for momentum (m)
    INTEGER(iwp) ::  ind_s_z0qh                      = 18  !< index for roughness length for heat (m)

    REAL(wp)  ::  ground_floor_level = 4.0_wp  !< default ground floor level

!
!-- Building facade/wall/green/window properties (partly according to PIDS).
!-- Initialization of building_pars is outsourced to usm_init_pars. This is needed because of the
!-- huge number of attributes given in building_pars (>700), while intel and gfortran compiler have
!-- hard limit of continuation lines of 511.
    REAL(wp), DIMENSION(0:149,1:7) ::  building_pars  !<
!
!-- Type for 1d surface variables as surface temperature and liquid water reservoir
    TYPE surf_type_1d_usm
       REAL(wp), DIMENSION(:), ALLOCATABLE         ::  val  !<
    END TYPE surf_type_1d_usm
!
!-- Type for 2d surface variables as wall temperature
    TYPE surf_type_2d_usm
       REAL(wp), DIMENSION(:,:), ALLOCATABLE       ::  val  !<
    END TYPE surf_type_2d_usm
!-- Wall surface model
!-- Wall surface model constants
    INTEGER(iwp), PARAMETER                        ::  nzb_wall = 0  !< inner side of the wall model (to be switched)
    INTEGER(iwp), PARAMETER                        ::  nzt_wall = 3  !< outer side of the wall model (to be switched)
    INTEGER(iwp), PARAMETER                        ::  nzw      = 4  !< number of wall layers (fixed for now)

    INTEGER(iwp)                                   ::  soil_type     !<

    REAL(wp)  ::  m_total                  = 0.0_wp    !< weighted total water content of the soil (m3/m3)
    REAL(wp)  ::  roof_inner_temperature   = 295.0_wp  !< temperature of the inner roof
                                                       !< surface (~22 degrees C) (K)
    REAL(wp)  ::  soil_inner_temperature   = 288.0_wp  !< temperature of the deep soil
                                                       !< (~15 degrees C) (K)
    REAL(wp)  ::  wall_inner_temperature   = 295.0_wp  !< temperature of the inner wall
                                                       !< surface (~22 degrees C) (K)
    REAL(wp)  ::  window_inner_temperature = 295.0_wp  !< temperature of the inner window
                                                       !< surface (~22 degrees C) (K)
!
!-- Surface and material model variables for walls, ground, roofs
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_green_h      !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_green_h_p    !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_wall_h       !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_wall_h_p     !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_window_h     !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_window_h_p   !<

    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  t_surf_green_h_1    !<
    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  t_surf_green_h_2    !<
    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  t_surf_wall_h_1     !<
    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  t_surf_wall_h_2     !<
    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  t_surf_window_h_1   !<
    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  t_surf_window_h_2   !<

    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_green_v      !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_green_v_p    !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_wall_v       !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_wall_v_p     !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_window_v     !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  t_surf_window_v_p   !<

    TYPE(surf_type_1d_usm), DIMENSION(0:3), TARGET  ::  t_surf_green_v_1    !<
    TYPE(surf_type_1d_usm), DIMENSION(0:3), TARGET  ::  t_surf_green_v_2    !<
    TYPE(surf_type_1d_usm), DIMENSION(0:3), TARGET  ::  t_surf_wall_v_1     !<
    TYPE(surf_type_1d_usm), DIMENSION(0:3), TARGET  ::  t_surf_wall_v_2     !<
    TYPE(surf_type_1d_usm), DIMENSION(0:3), TARGET  ::  t_surf_window_v_1   !<
    TYPE(surf_type_1d_usm), DIMENSION(0:3), TARGET  ::  t_surf_window_v_2   !<

!
!-- Energy balance variables
!-- Parameters of the land, roof and wall surfaces (only for horizontal upward surfaces)
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  fc_h          !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  rootfr_h      !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  swc_h         !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  swc_h_p       !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  swc_res_h     !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  swc_sat_h     !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_green_h     !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_green_h_p   !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_wall_h      !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_wall_h_p    !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  wilt_h        !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_window_h    !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_window_h_p  !<


    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  fc_h_1        !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  rootfr_h_1    !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  swc_h_1       !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  swc_h_2       !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  swc_res_h_1   !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  swc_sat_h_1   !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  t_green_h_1   !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  t_green_h_2   !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  t_wall_h_1    !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  t_wall_h_2    !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  wilt_h_1      !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  t_window_h_1  !<
    TYPE(surf_type_2d_usm), DIMENSION(0:1), TARGET  ::  t_window_h_2  !<

    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_green_v     !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_green_v_p   !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_wall_v      !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_wall_v_p    !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_window_v    !<
    TYPE(surf_type_2d_usm), DIMENSION(:), POINTER   ::  t_window_v_p  !<
    TYPE(surf_type_2d_usm), DIMENSION(0:3), TARGET  ::  t_green_v_1   !<
    TYPE(surf_type_2d_usm), DIMENSION(0:3), TARGET  ::  t_green_v_2   !<
    TYPE(surf_type_2d_usm), DIMENSION(0:3), TARGET  ::  t_wall_v_1    !<
    TYPE(surf_type_2d_usm), DIMENSION(0:3), TARGET  ::  t_wall_v_2    !<
    TYPE(surf_type_2d_usm), DIMENSION(0:3), TARGET  ::  t_window_v_1  !<
    TYPE(surf_type_2d_usm), DIMENSION(0:3), TARGET  ::  t_window_v_2  !<
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  m_liq_usm_h    !< liquid water reservoir (m), horizontal surface elements
    TYPE(surf_type_1d_usm), DIMENSION(:), POINTER   ::  m_liq_usm_h_p  !< progn. liquid water reservoir (m), horizontal surface elements
    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  m_liq_usm_h_1  !<
    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  m_liq_usm_h_2  !<
    TYPE(surf_type_1d_usm), DIMENSION(0:1), TARGET  ::  tm_liq_usm_h_m  !< liquid water reservoir tendency (m), horizontal surface elements
!-- Interfaces of subroutines accessed from outside of this module
    INTERFACE usm_3d_data_averaging
       MODULE PROCEDURE usm_3d_data_averaging
    END INTERFACE usm_3d_data_averaging

    INTERFACE usm_boundary_condition
       MODULE PROCEDURE usm_boundary_condition
    END INTERFACE usm_boundary_condition

    INTERFACE usm_check_data_output
       MODULE PROCEDURE usm_check_data_output
    END INTERFACE usm_check_data_output

    INTERFACE usm_check_parameters
       MODULE PROCEDURE usm_check_parameters
    END INTERFACE usm_check_parameters

    INTERFACE usm_data_output_3d
       MODULE PROCEDURE usm_data_output_3d
    END INTERFACE usm_data_output_3d

    INTERFACE usm_define_netcdf_grid
       MODULE PROCEDURE usm_define_netcdf_grid
    END INTERFACE usm_define_netcdf_grid

    INTERFACE usm_init
       MODULE PROCEDURE usm_init
    END INTERFACE usm_init

    INTERFACE usm_init_arrays
       MODULE PROCEDURE usm_init_arrays
    END INTERFACE usm_init_arrays

    INTERFACE usm_parin
       MODULE PROCEDURE usm_parin
    END INTERFACE usm_parin

    INTERFACE usm_rrd_local
       MODULE PROCEDURE usm_rrd_local_ftn
       MODULE PROCEDURE usm_rrd_local_mpi
    END INTERFACE usm_rrd_local

    INTERFACE usm_energy_balance
       MODULE PROCEDURE usm_energy_balance
    END INTERFACE usm_energy_balance

    INTERFACE usm_swap_timelevel
       MODULE PROCEDURE usm_swap_timelevel
    END INTERFACE usm_swap_timelevel

    INTERFACE usm_wrd_local
       MODULE PROCEDURE usm_wrd_local
    END INTERFACE usm_wrd_local


    SAVE

    PRIVATE

!
!-- Public functions
    PUBLIC usm_boundary_condition,                                                                 &
           usm_check_data_output,                                                                  &
           usm_check_parameters,                                                                   &
           usm_data_output_3d,                                                                     &
           usm_define_netcdf_grid,                                                                 &
           usm_init,                                                                               &
           usm_init_arrays,                                                                        &
           usm_parin,                                                                              &
           usm_rrd_local,                                                                          &
           usm_energy_balance,                                                                     &
           usm_swap_timelevel,                                                                     &
           usm_wrd_local,                                                                          &
           usm_3d_data_averaging

!
!-- Public parameters, constants and initial values
    PUBLIC building_type,                                                                          &
           building_pars,                                                                          &
           nzb_wall,                                                                               &
           nzt_wall,                                                                               &
           t_green_h,                                                                              &
           t_green_v,                                                                              &
           t_wall_h,                                                                               &
           t_wall_v,                                                                               &
           t_window_h,                                                                             &
           t_window_v,                                                                             &
           usm_wall_mod






 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This subroutine creates the necessary indices of the urban surfaces and plant canopy and it
!> allocates the needed arrays for USM
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_init_arrays

    IMPLICIT NONE

    INTEGER(iwp) ::  l  !<

    IF ( debug_output )  CALL debug_message( 'usm_init_arrays', 'start' )

!
!-- Allocate radiation arrays which are part of the new data type.
!-- For horizontal surfaces.
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%surfhf(1:surf_usm_h(l)%ns)    )
       ALLOCATE ( surf_usm_h(l)%rad_net_l(1:surf_usm_h(l)%ns) )
    ENDDO
!
!-- For vertical surfaces
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%surfhf(1:surf_usm_v(l)%ns)    )
       ALLOCATE ( surf_usm_v(l)%rad_net_l(1:surf_usm_v(l)%ns) )
    ENDDO

!
!-- Wall surface model
!-- Allocate arrays for wall surface model and define pointers
!-- Allocate array of wall types and wall parameters
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%surface_types(1:surf_usm_h(l)%ns)      )
       ALLOCATE ( surf_usm_h(l)%building_type(1:surf_usm_h(l)%ns)      )
       ALLOCATE ( surf_usm_h(l)%building_type_name(1:surf_usm_h(l)%ns) )
       surf_usm_h(l)%building_type      = 0
       surf_usm_h(l)%building_type_name = 'none'
    ENDDO
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%surface_types(1:surf_usm_v(l)%ns)      )
       ALLOCATE ( surf_usm_v(l)%building_type(1:surf_usm_v(l)%ns)      )
       ALLOCATE ( surf_usm_v(l)%building_type_name(1:surf_usm_v(l)%ns) )
       surf_usm_v(l)%building_type      = 0
       surf_usm_v(l)%building_type_name = 'none'
    ENDDO
!
!-- Allocate albedo_type and albedo. Each surface element has 3 values, 0: wall fraction,
!-- 1: green fraction, 2: window fraction.
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%albedo_type(1:surf_usm_h(l)%ns,0:2) )
       ALLOCATE ( surf_usm_h(l)%albedo(1:surf_usm_h(l)%ns,0:2)      )
       surf_usm_h(l)%albedo_type = albedo_type
    ENDDO
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%albedo_type(1:surf_usm_v(l)%ns,0:2) )
       ALLOCATE ( surf_usm_v(l)%albedo(1:surf_usm_v(l)%ns,0:2)      )
       surf_usm_v(l)%albedo_type = albedo_type
    ENDDO

!
!-- Allocate indoor target temperature for summer and winter
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%target_temp_summer(1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%target_temp_winter(1:surf_usm_h(l)%ns) )
    ENDDO
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%target_temp_summer(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%target_temp_winter(1:surf_usm_v(l)%ns) )
    ENDDO
!
!-- In case the indoor model is applied, allocate memory for waste heat and indoor temperature.
    IF ( indoor_model )  THEN
       DO  l = 0, 1
          ALLOCATE ( surf_usm_h(l)%waste_heat(1:surf_usm_h(l)%ns) )
          surf_usm_h(l)%waste_heat = 0.0_wp
       ENDDO
       DO  l = 0, 3
          ALLOCATE ( surf_usm_v(l)%waste_heat(1:surf_usm_v(l)%ns) )
          surf_usm_v(l)%waste_heat = 0.0_wp
       ENDDO
    ENDIF
!
!-- Allocate flag indicating ground floor level surface elements
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%ground_level(1:surf_usm_h(l)%ns) )
    ENDDO
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%ground_level(1:surf_usm_v(l)%ns) )
    ENDDO
!
!-- Allocate arrays for relative surface fraction.
!-- 0 - wall fraction, 1 - green fraction, 2 - window fraction
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%frac(1:surf_usm_h(l)%ns,0:2) )
       surf_usm_h(l)%frac = 0.0_wp
    ENDDO
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%frac(1:surf_usm_v(l)%ns,0:2) )
       surf_usm_v(l)%frac = 0.0_wp
    ENDDO

!
!-- Wall and roof surface parameters. First for horizontal surfaces
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%isroof_surf(1:surf_usm_h(l)%ns)        )
       ALLOCATE ( surf_usm_h(l)%lambda_surf(1:surf_usm_h(l)%ns)        )
       ALLOCATE ( surf_usm_h(l)%lambda_surf_window(1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%lambda_surf_green(1:surf_usm_h(l)%ns)  )
       ALLOCATE ( surf_usm_h(l)%c_surface(1:surf_usm_h(l)%ns)          )
       ALLOCATE ( surf_usm_h(l)%c_surface_window(1:surf_usm_h(l)%ns)   )
       ALLOCATE ( surf_usm_h(l)%c_surface_green(1:surf_usm_h(l)%ns)    )
       ALLOCATE ( surf_usm_h(l)%transmissivity(1:surf_usm_h(l)%ns)     )
       ALLOCATE ( surf_usm_h(l)%lai(1:surf_usm_h(l)%ns)                )
       ALLOCATE ( surf_usm_h(l)%emissivity(1:surf_usm_h(l)%ns,0:2)     )
       ALLOCATE ( surf_usm_h(l)%r_a(1:surf_usm_h(l)%ns)                )
       ALLOCATE ( surf_usm_h(l)%r_a_green(1:surf_usm_h(l)%ns)          )
       ALLOCATE ( surf_usm_h(l)%r_a_window(1:surf_usm_h(l)%ns)         )
       ALLOCATE ( surf_usm_h(l)%green_type_roof(1:surf_usm_h(l)%ns)    )
       ALLOCATE ( surf_usm_h(l)%r_s(1:surf_usm_h(l)%ns)                )
    ENDDO
!
!-- For vertical surfaces.
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%lambda_surf(1:surf_usm_v(l)%ns)        )
       ALLOCATE ( surf_usm_v(l)%c_surface(1:surf_usm_v(l)%ns)          )
       ALLOCATE ( surf_usm_v(l)%lambda_surf_window(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%c_surface_window(1:surf_usm_v(l)%ns)   )
       ALLOCATE ( surf_usm_v(l)%lambda_surf_green(1:surf_usm_v(l)%ns)  )
       ALLOCATE ( surf_usm_v(l)%c_surface_green(1:surf_usm_v(l)%ns)    )
       ALLOCATE ( surf_usm_v(l)%transmissivity(1:surf_usm_v(l)%ns)     )
       ALLOCATE ( surf_usm_v(l)%lai(1:surf_usm_v(l)%ns)                )
       ALLOCATE ( surf_usm_v(l)%emissivity(1:surf_usm_v(l)%ns,0:2)     )
       ALLOCATE ( surf_usm_v(l)%r_a(1:surf_usm_v(l)%ns)                )
       ALLOCATE ( surf_usm_v(l)%r_a_green(1:surf_usm_v(l)%ns)          )
       ALLOCATE ( surf_usm_v(l)%r_a_window(1:surf_usm_v(l)%ns)         )
       ALLOCATE ( surf_usm_v(l)%r_s(1:surf_usm_v(l)%ns)                )
    ENDDO

!
!-- Allocate wall and roof material parameters. First for horizontal surfaces
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%thickness_wall(1:surf_usm_h(l)%ns)                    )
       ALLOCATE ( surf_usm_h(l)%thickness_window(1:surf_usm_h(l)%ns)                  )
       ALLOCATE ( surf_usm_h(l)%thickness_green(1:surf_usm_h(l)%ns)                   )
       ALLOCATE ( surf_usm_h(l)%lambda_h(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)        )
       ALLOCATE ( surf_usm_h(l)%rho_c_wall(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)      )
       ALLOCATE ( surf_usm_h(l)%lambda_h_window(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%rho_c_window(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)    )
       ALLOCATE ( surf_usm_h(l)%lambda_h_green(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)  )
       ALLOCATE ( surf_usm_h(l)%rho_c_green(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)     )

       ALLOCATE ( surf_usm_h(l)%rho_c_total_green(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)    )
       ALLOCATE ( surf_usm_h(l)%n_vg_green(1:surf_usm_h(l)%ns)                             )
       ALLOCATE ( surf_usm_h(l)%alpha_vg_green(1:surf_usm_h(l)%ns)                         )
       ALLOCATE ( surf_usm_h(l)%l_vg_green(1:surf_usm_h(l)%ns)                             )
       ALLOCATE ( surf_usm_h(l)%gamma_w_green_sat(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)  )
       ALLOCATE ( surf_usm_h(l)%lambda_w_green(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)       )
       ALLOCATE ( surf_usm_h(l)%gamma_w_green(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)        )
       ALLOCATE ( surf_usm_h(l)%tswc_h_m(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)             )
    ENDDO
!
!-- For vertical surfaces.
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%thickness_wall(1:surf_usm_v(l)%ns)                    )
       ALLOCATE ( surf_usm_v(l)%thickness_window(1:surf_usm_v(l)%ns)                  )
       ALLOCATE ( surf_usm_v(l)%thickness_green(1:surf_usm_v(l)%ns)                   )
       ALLOCATE ( surf_usm_v(l)%lambda_h(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)        )
       ALLOCATE ( surf_usm_v(l)%rho_c_wall(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)      )
       ALLOCATE ( surf_usm_v(l)%lambda_h_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%rho_c_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)    )
       ALLOCATE ( surf_usm_v(l)%lambda_h_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
       ALLOCATE ( surf_usm_v(l)%rho_c_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)     )
    ENDDO

!
!-- Allocate green wall and roof vegetation and soil parameters. First horizontal surfaces
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%g_d(1:surf_usm_h(l)%ns)              )
       ALLOCATE ( surf_usm_h(l)%c_liq(1:surf_usm_h(l)%ns)            )
       ALLOCATE ( surf_usm_h(l)%qsws_liq(1:surf_usm_h(l)%ns)         )
       ALLOCATE ( surf_usm_h(l)%qsws_veg(1:surf_usm_h(l)%ns)         )
       ALLOCATE ( surf_usm_h(l)%r_canopy(1:surf_usm_h(l)%ns)         )
       ALLOCATE ( surf_usm_h(l)%r_canopy_min(1:surf_usm_h(l)%ns)     )
       ALLOCATE ( surf_usm_h(l)%pt_10cm(1:surf_usm_h(l)%ns)          )
    ENDDO
!
!-- For vertical surfaces.
    DO  l = 0, 3
      ALLOCATE ( surf_usm_v(l)%g_d(1:surf_usm_v(l)%ns)              )
      ALLOCATE ( surf_usm_v(l)%c_liq(1:surf_usm_v(l)%ns)            )
      ALLOCATE ( surf_usm_v(l)%qsws_liq(1:surf_usm_v(l)%ns)         )
      ALLOCATE ( surf_usm_v(l)%qsws_veg(1:surf_usm_v(l)%ns)         )
      ALLOCATE ( surf_usm_v(l)%r_canopy(1:surf_usm_v(l)%ns)         )
      ALLOCATE ( surf_usm_v(l)%r_canopy_min(1:surf_usm_v(l)%ns)     )
      ALLOCATE ( surf_usm_v(l)%pt_10cm(1:surf_usm_v(l)%ns)          )
    ENDDO

!
!-- Allocate wall and roof layers sizes. For horizontal surfaces.
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%dz_wall(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)       )
       ALLOCATE ( surf_usm_h(l)%dz_window(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)     )
       ALLOCATE ( surf_usm_h(l)%dz_green(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)      )
       ALLOCATE ( surf_usm_h(l)%ddz_wall(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)      )
       ALLOCATE ( surf_usm_h(l)%dz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)    )
       ALLOCATE ( surf_usm_h(l)%ddz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)   )
       ALLOCATE ( surf_usm_h(l)%zw(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)              )
       ALLOCATE ( surf_usm_h(l)%ddz_window(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)    )
       ALLOCATE ( surf_usm_h(l)%dz_window_stag(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)  )
       ALLOCATE ( surf_usm_h(l)%ddz_window_stag(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%zw_window(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)       )
       ALLOCATE ( surf_usm_h(l)%ddz_green(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)     )
       ALLOCATE ( surf_usm_h(l)%dz_green_stag(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)   )
       ALLOCATE ( surf_usm_h(l)%ddz_green_stag(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)  )
       ALLOCATE ( surf_usm_h(l)%zw_green(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns)        )
    ENDDO

!
!-- For vertical surfaces.
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%dz_wall(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)       )
       ALLOCATE ( surf_usm_v(l)%dz_window(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)     )
       ALLOCATE ( surf_usm_v(l)%dz_green(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)      )
       ALLOCATE ( surf_usm_v(l)%ddz_wall(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)      )
       ALLOCATE ( surf_usm_v(l)%dz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)    )
       ALLOCATE ( surf_usm_v(l)%ddz_wall_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)   )
       ALLOCATE ( surf_usm_v(l)%zw(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)              )
       ALLOCATE ( surf_usm_v(l)%ddz_window(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)    )
       ALLOCATE ( surf_usm_v(l)%dz_window_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
       ALLOCATE ( surf_usm_v(l)%ddz_window_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%zw_window(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)       )
       ALLOCATE ( surf_usm_v(l)%ddz_green(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns)     )
       ALLOCATE ( surf_usm_v(l)%dz_green_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)   )
       ALLOCATE ( surf_usm_v(l)%ddz_green_stag(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)  )
       ALLOCATE ( surf_usm_v(l)%zw_green(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns)        )
    ENDDO

!
!-- Allocate wall and roof temperature arrays, for horizontal walls.
!-- Allocate if required. Note, in case of restarts, some of these arrays might be already allocated.
    DO l = 0, 1
       IF ( .NOT. ALLOCATED( t_surf_wall_h_1(l)%val ) )                                                      &
          ALLOCATE ( t_surf_wall_h_1(l)%val(1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_wall_h_2(l)%val ) )                                                      &
          ALLOCATE ( t_surf_wall_h_2(l)%val(1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_wall_h_1(l)%val ) )                                                           &
          ALLOCATE ( t_wall_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_wall_h_2(l)%val ) )                                                           &
          ALLOCATE ( t_wall_h_2(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_window_h_1(l)%val ) )                                                    &
          ALLOCATE ( t_surf_window_h_1(l)%val(1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_window_h_2(l)%val ) )                                                    &
          ALLOCATE ( t_surf_window_h_2(l)%val(1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_window_h_1(l)%val ) )                                                         &
          ALLOCATE ( t_window_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_window_h_2(l)%val ) )                                                         &
          ALLOCATE ( t_window_h_2(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_green_h_1(l)%val ) )                                                     &
          ALLOCATE ( t_surf_green_h_1(l)%val(1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_green_h_2(l)%val ) )                                                     &
          ALLOCATE ( t_surf_green_h_2(l)%val(1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_green_h_1(l)%val ) )                                                          &
          ALLOCATE ( t_green_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( t_green_h_2(l)%val ) )                                                          &
          ALLOCATE ( t_green_h_2(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( swc_h_1(l)%val ) )                                                              &
          ALLOCATE ( swc_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( swc_sat_h_1(l)%val ) )                                                          &
          ALLOCATE ( swc_sat_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( swc_res_h_1(l)%val ) )                                                          &
          ALLOCATE ( swc_res_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( swc_h_2(l)%val ) )                                                              &
          ALLOCATE ( swc_h_2(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( rootfr_h_1(l)%val ) )                                                           &
          ALLOCATE ( rootfr_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( wilt_h_1(l)%val ) )                                                             &
          ALLOCATE ( wilt_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( fc_h_1(l)%val ) )                                                               &
          ALLOCATE ( fc_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )

       IF ( .NOT. ALLOCATED( m_liq_usm_h_1(l)%val ) )                                             &
          ALLOCATE ( m_liq_usm_h_1(l)%val(1:surf_usm_h(l)%ns) )
       IF ( .NOT. ALLOCATED( m_liq_usm_h_2(l)%val ) )                                             &
          ALLOCATE ( m_liq_usm_h_2(l)%val(1:surf_usm_h(l)%ns) )
    ENDDO
!
!-- Initial assignment of the pointers
    t_wall_h    => t_wall_h_1;   t_wall_h_p   => t_wall_h_2
    t_window_h  => t_window_h_1; t_window_h_p => t_window_h_2
    t_green_h   => t_green_h_1;  t_green_h_p  => t_green_h_2
    t_surf_wall_h   => t_surf_wall_h_1;   t_surf_wall_h_p   => t_surf_wall_h_2
    t_surf_window_h => t_surf_window_h_1; t_surf_window_h_p => t_surf_window_h_2
    t_surf_green_h  => t_surf_green_h_1;  t_surf_green_h_p  => t_surf_green_h_2
    m_liq_usm_h     => m_liq_usm_h_1;     m_liq_usm_h_p     => m_liq_usm_h_2
    swc_h     => swc_h_1; swc_h_p => swc_h_2
    swc_sat_h => swc_sat_h_1
    swc_res_h => swc_res_h_1
    rootfr_h  => rootfr_h_1
    wilt_h    => wilt_h_1
    fc_h      => fc_h_1

!
!-- Allocate wall and roof temperature arrays, for vertical walls if required.
!-- Allocate if required. Note, in case of restarts, some of these arrays might be already allocated.
    DO  l = 0, 3
       IF ( .NOT. ALLOCATED( t_surf_wall_v_1(l)%val ) )                                              &
          ALLOCATE ( t_surf_wall_v_1(l)%val(1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_wall_v_2(l)%val ) )                                              &
          ALLOCATE ( t_surf_wall_v_2(l)%val(1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_wall_v_1(l)%val ) )                                                   &
          ALLOCATE ( t_wall_v_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_wall_v_2(l)%val ) )                                                   &
          ALLOCATE ( t_wall_v_2(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_window_v_1(l)%val ) )                                            &
          ALLOCATE ( t_surf_window_v_1(l)%val(1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_window_v_2(l)%val ) )                                            &
          ALLOCATE ( t_surf_window_v_2(l)%val(1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_window_v_1(l)%val ) )                                                 &
          ALLOCATE ( t_window_v_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_window_v_2(l)%val ) )                                                 &
          ALLOCATE ( t_window_v_2(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_green_v_1(l)%val ) )                                             &
          ALLOCATE ( t_surf_green_v_1(l)%val(1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_surf_green_v_2(l)%val ) )                                             &
          ALLOCATE ( t_surf_green_v_2(l)%val(1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_green_v_1(l)%val ) )                                                  &
          ALLOCATE ( t_green_v_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       IF ( .NOT. ALLOCATED( t_green_v_2(l)%val ) )                                                  &
          ALLOCATE ( t_green_v_2(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
    ENDDO
!
!-- Initial assignment of the pointers
    t_wall_v        => t_wall_v_1;        t_wall_v_p        => t_wall_v_2
    t_surf_wall_v   => t_surf_wall_v_1;   t_surf_wall_v_p   => t_surf_wall_v_2
    t_window_v      => t_window_v_1;      t_window_v_p      => t_window_v_2
    t_green_v       => t_green_v_1;       t_green_v_p       => t_green_v_2
    t_surf_window_v => t_surf_window_v_1; t_surf_window_v_p => t_surf_window_v_2
    t_surf_green_v  => t_surf_green_v_1;  t_surf_green_v_p  => t_surf_green_v_2

!
!-- Allocate intermediate timestep arrays. For horizontal surfaces.
    DO  l = 0, 1
       ALLOCATE ( surf_usm_h(l)%tt_surface_wall_m(1:surf_usm_h(l)%ns)               )
       ALLOCATE ( surf_usm_h(l)%tt_wall_m(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)   )
       ALLOCATE ( surf_usm_h(l)%tt_surface_window_m(1:surf_usm_h(l)%ns)             )
       ALLOCATE ( surf_usm_h(l)%tt_window_m(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%tt_green_m(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns)  )
       ALLOCATE ( surf_usm_h(l)%tt_surface_green_m(1:surf_usm_h(l)%ns)              )
!
!--    Allocate intermediate timestep arrays
!--    Horizontal surfaces
       ALLOCATE ( tm_liq_usm_h_m(l)%val(1:surf_usm_h(l)%ns) )
       tm_liq_usm_h_m(l)%val = 0.0_wp
!
!--    Set inital values for prognostic quantities
       IF ( ALLOCATED( surf_usm_h(l)%tt_surface_wall_m )   )  surf_usm_h(l)%tt_surface_wall_m   = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%tt_wall_m )           )  surf_usm_h(l)%tt_wall_m           = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%tt_surface_window_m ) )  surf_usm_h(l)%tt_surface_window_m = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%tt_window_m    )      )  surf_usm_h(l)%tt_window_m         = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%tt_green_m    )       )  surf_usm_h(l)%tt_green_m          = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%tt_surface_green_m )  )  surf_usm_h(l)%tt_surface_green_m  = 0.0_wp
    END DO
!
!-- Now, for vertical surfaces
    DO  l = 0, 3
       ALLOCATE ( surf_usm_v(l)%tt_surface_wall_m(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%tt_wall_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       IF ( ALLOCATED( surf_usm_v(l)%tt_surface_wall_m ) )  surf_usm_v(l)%tt_surface_wall_m = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%tt_wall_m ) )  surf_usm_v(l)%tt_wall_m  = 0.0_wp
       ALLOCATE ( surf_usm_v(l)%tt_surface_window_m(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%tt_window_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       IF ( ALLOCATED( surf_usm_v(l)%tt_surface_window_m ) )  surf_usm_v(l)%tt_surface_window_m = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%tt_window_m ) )  surf_usm_v(l)%tt_window_m = 0.0_wp
       ALLOCATE ( surf_usm_v(l)%tt_surface_green_m(1:surf_usm_v(l)%ns) )
       IF ( ALLOCATED( surf_usm_v(l)%tt_surface_green_m ) )  surf_usm_v(l)%tt_surface_green_m = 0.0_wp
       ALLOCATE ( surf_usm_v(l)%tt_green_m(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       IF ( ALLOCATED( surf_usm_v(l)%tt_green_m ) )  surf_usm_v(l)%tt_green_m = 0.0_wp
    ENDDO
!
!-- Allocate wall heat flux output arrays and set initial values. For horizontal surfaces
    DO  l = 0, 1
!      ALLOCATE ( surf_usm_h(l)%wshf(1:surf_usm_h(l)%ns)    )  !can be removed
       ALLOCATE ( surf_usm_h(l)%ghf(1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%wshf_eb(1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%wghf_eb(1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%wghf_eb_window(1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%wghf_eb_green(1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%iwghf_eb(1:surf_usm_h(l)%ns) )
       ALLOCATE ( surf_usm_h(l)%iwghf_eb_window(1:surf_usm_h(l)%ns) )
       IF ( ALLOCATED( surf_usm_h(l)%ghf     ) )  surf_usm_h(l)%ghf     = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%wshf    ) )  surf_usm_h(l)%wshf    = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%wshf_eb ) )  surf_usm_h(l)%wshf_eb = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%wghf_eb ) )  surf_usm_h(l)%wghf_eb = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%wghf_eb_window ) )   surf_usm_h(l)%wghf_eb_window  = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%wghf_eb_green ) )    surf_usm_h(l)%wghf_eb_green   = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%iwghf_eb ) )         surf_usm_h(l)%iwghf_eb        = 0.0_wp
       IF ( ALLOCATED( surf_usm_h(l)%iwghf_eb_window ) )  surf_usm_h(l)%iwghf_eb_window = 0.0_wp
    ENDDO
!
!-- Now, for vertical surfaces
    DO  l = 0, 3
!       ALLOCATE ( surf_usm_v(l)%wshf(1:surf_usm_v(l)%ns)    )    ! can be removed
       ALLOCATE ( surf_usm_v(l)%ghf(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%wshf_eb(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%wghf_eb(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%wghf_eb_window(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%wghf_eb_green(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%iwghf_eb(1:surf_usm_v(l)%ns) )
       ALLOCATE ( surf_usm_v(l)%iwghf_eb_window(1:surf_usm_v(l)%ns) )
       IF ( ALLOCATED( surf_usm_v(l)%ghf     ) )  surf_usm_v(l)%ghf     = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%wshf    ) )  surf_usm_v(l)%wshf    = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%wshf_eb ) )  surf_usm_v(l)%wshf_eb = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%wghf_eb ) )  surf_usm_v(l)%wghf_eb = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%wghf_eb_window ) )   surf_usm_v(l)%wghf_eb_window  = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%wghf_eb_green ) )    surf_usm_v(l)%wghf_eb_green   = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%iwghf_eb ) )         surf_usm_v(l)%iwghf_eb        = 0.0_wp
       IF ( ALLOCATED( surf_usm_v(l)%iwghf_eb_window ) )  surf_usm_v(l)%iwghf_eb_window = 0.0_wp
    ENDDO
!
!-- Initialize building-surface properties, which are also required by other modules, e.g. the
!-- indoor model.
    CALL usm_define_pars

    IF ( debug_output )  CALL debug_message( 'usm_init_arrays', 'end' )

 END SUBROUTINE usm_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average urban surface output quantities as well as allocate the array necessary
!> for storing the average.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_3d_data_averaging( mode, variable )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) ::  variable  !<
    CHARACTER(LEN=*), INTENT(IN) ::  mode      !<

    INTEGER(iwp)                                       ::  i, j, k, l, m, ids, idsint, iwl, istat  !< runnin indices
    CHARACTER(LEN=varnamelength)                       ::  var                                     !< trimmed variable
    LOGICAL                                            ::  horizontal

    IF ( .NOT. variable(1:4) == 'usm_' )  RETURN  ! Is such a check really required?

!
!-- Find the real name of the variable
    ids = -1
    l = -1
    var = TRIM(variable)
    DO  i = 0, nd-1
       k = len( TRIM( var ) )
       j = len( TRIM( dirname(i) ) )
       IF ( TRIM( var(k-j+1:k) ) == TRIM( dirname(i) ) )  THEN
           ids = i
           idsint = dirint(ids)
           l = diridx(ids)  !> index of direction for _h and _v arrays
           var = var(:k-j)
           EXIT
       ENDIF
    ENDDO
    IF ( ids == -1 )  THEN
        var = TRIM( variable )
    ELSE
!--    Horizontal direction flag
       IF ( idsint == iup .OR. idsint == idown )  THEN
          horizontal = .TRUE.
       ELSE
          horizontal = .FALSE.
       ENDIF
    ENDIF
    IF ( var(1:11) == 'usm_t_wall_'  .AND.  len( TRIM( var ) ) >= 12 )  THEN
!
!--    Wall layers
       READ( var(12:12), '(I1)', iostat=istat ) iwl
       IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
          var = var(1:10)
       ELSE
!
!--       Wrong wall layer index
          RETURN
       ENDIF
    ENDIF
    IF ( var(1:13) == 'usm_t_window_'  .AND.  len( TRIM(var) ) >= 14 )  THEN
!
!--      Wall layers
        READ( var(14:14), '(I1)', iostat=istat ) iwl
        IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
            var = var(1:12)
        ELSE
!
!--         Wrong window layer index
            RETURN
        ENDIF
    ENDIF
    IF ( var(1:12) == 'usm_t_green_'  .AND.  len( TRIM( var ) ) >= 13 )  THEN
!
!--      Wall layers
        READ( var(13:13), '(I1)', iostat=istat ) iwl
        IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
            var = var(1:11)
        ELSE
!
!--         Wrong green layer index
            RETURN
        ENDIF
    ENDIF
    IF ( var(1:8) == 'usm_swc_'  .AND.  len( TRIM( var ) ) >= 9 )  THEN
!
!--      Swc layers
        READ( var(9:9), '(I1)', iostat=istat ) iwl
        IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
            var = var(1:7)
        ELSE
!
!--         Wrong swc layer index
            RETURN
        ENDIF
    ENDIF

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( var ) )

            CASE ( 'usm_wshf' )
!
!--             Array of sensible heat flux from surfaces
!--             Land surfaces
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%wshf_eb_av ) )  THEN
                      ALLOCATE ( surf_usm_h(l)%wshf_eb_av(1:surf_usm_h(l)%ns) )
                      surf_usm_h(l)%wshf_eb_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%wshf_eb_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%wshf_eb_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%wshf_eb_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_qsws' )
!
!--             Array of latent heat flux from surfaces
!--             Land surfaces
                IF ( horizontal .AND. .NOT.  ALLOCATED( surf_usm_h(l)%qsws_av ) )  THEN
                    ALLOCATE ( surf_usm_h(l)%qsws_av(1:surf_usm_h(l)%ns) )
                    surf_usm_h(l)%qsws_av = 0.0_wp
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%qsws_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%qsws_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%qsws_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_qsws_veg' )
!
!--             Array of latent heat flux from vegetation surfaces
!--             Land surfaces
                IF ( horizontal .AND. .NOT.  ALLOCATED( surf_usm_h(l)%qsws_veg_av ) )  THEN
                    ALLOCATE ( surf_usm_h(l)%qsws_veg_av(1:surf_usm_h(l)%ns) )
                    surf_usm_h(l)%qsws_veg_av = 0.0_wp
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%qsws_veg_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%qsws_veg_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%qsws_veg_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_qsws_liq' )
!
!--             Array of latent heat flux from surfaces with liquid
!--             Land surfaces
                IF ( horizontal .AND. .NOT.  ALLOCATED( surf_usm_h(l)%qsws_liq_av ) )  THEN
                    ALLOCATE ( surf_usm_h(l)%qsws_liq_av(1:surf_usm_h(l)%ns) )
                    surf_usm_h(l)%qsws_liq_av = 0.0_wp
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%qsws_liq_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%qsws_liq_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%qsws_liq_av = 0.0_wp
                   ENDIF
                ENDIF
!
!--         Please note, the following output quantities belongs to the individual tile fractions -
!--         ground heat flux at wall-, window-, and green fraction. Aggregated ground-heat flux is
!--         treated accordingly in average_3d_data, sum_up_3d_data, etc..
            CASE ( 'usm_wghf' )
!
!--             Array of heat flux from ground (wall, roof, land)
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%wghf_eb_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%wghf_eb_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%wghf_eb_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%wghf_eb_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%wghf_eb_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%wghf_eb_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_wghf_window' )
!
!--             Array of heat flux from window ground (wall, roof, land)
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%wghf_eb_window_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%wghf_eb_window_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%wghf_eb_window_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%wghf_eb_window_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%wghf_eb_window_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%wghf_eb_window_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_wghf_green' )
!
!--             Array of heat flux from green ground (wall, roof, land)
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%wghf_eb_green_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%wghf_eb_green_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%wghf_eb_green_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%wghf_eb_green_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%wghf_eb_green_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%wghf_eb_green_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_iwghf' )
!
!--             Array of heat flux from indoor ground (wall, roof, land)
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%iwghf_eb_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%iwghf_eb_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%iwghf_eb_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%iwghf_eb_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%iwghf_eb_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%iwghf_eb_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_iwghf_window' )
!
!--             Array of heat flux from indoor window ground (wall, roof, land)
                IF ( horizontal ) THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%iwghf_eb_window_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%iwghf_eb_window_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%iwghf_eb_window_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%iwghf_eb_window_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%iwghf_eb_window_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%iwghf_eb_window_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_t_surf_wall' )
!
!--             Surface temperature for surfaces
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%t_surf_wall_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%t_surf_wall_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%t_surf_wall_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%t_surf_wall_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%t_surf_wall_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%t_surf_wall_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_t_surf_window' )
!
!--             Surface temperature for window surfaces
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%t_surf_window_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%t_surf_window_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%t_surf_window_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%t_surf_window_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%t_surf_window_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%t_surf_window_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_t_surf_green' )
!
!--             Surface temperature for green surfaces
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%t_surf_green_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%t_surf_green_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%t_surf_green_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%t_surf_green_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%t_surf_green_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%t_surf_green_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_theta_10cm' )
!
!--             Near surface (10cm) temperature for whole surfaces
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%pt_10cm_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%pt_10cm_av(1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%pt_10cm_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%pt_10cm_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%pt_10cm_av(1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%pt_10cm_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_t_wall' )
!
!--             Wall temperature for iwl layer of walls and land
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%t_wall_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%t_wall_av(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%t_wall_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%t_wall_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%t_wall_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%t_wall_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_t_window' )
!
!--             Window temperature for iwl layer of walls and land
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%t_window_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%t_window_av(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%t_window_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%t_window_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%t_window_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%t_window_av = 0.0_wp
                   ENDIF
                ENDIF

            CASE ( 'usm_t_green' )
!
!--             Green temperature for iwl layer of walls and land
                IF ( horizontal )  THEN
                   IF ( .NOT.  ALLOCATED( surf_usm_h(l)%t_green_av ) )  THEN
                       ALLOCATE ( surf_usm_h(l)%t_green_av(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns) )
                       surf_usm_h(l)%t_green_av = 0.0_wp
                   ENDIF
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%t_green_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%t_green_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%t_green_av = 0.0_wp
                   ENDIF
                ENDIF
            CASE ( 'usm_swc' )
!
!--             Soil water content for iwl layer of walls and land
                IF ( horizontal .AND. .NOT.  ALLOCATED( surf_usm_h(l)%swc_av ) )  THEN
                    ALLOCATE ( surf_usm_h(l)%swc_av(nzb_wall:nzt_wall,1:surf_usm_h(l)%ns) )
                    surf_usm_h(l)%swc_av = 0.0_wp
                ELSE
                   IF ( .NOT.  ALLOCATED( surf_usm_v(l)%swc_av ) )  THEN
                       ALLOCATE ( surf_usm_v(l)%swc_av(nzb_wall:nzt_wall,1:surf_usm_v(l)%ns) )
                       surf_usm_v(l)%swc_av = 0.0_wp
                   ENDIF
                ENDIF

           CASE DEFAULT
               CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( var ) )

            CASE ( 'usm_wshf' )
!
!--             Array of sensible heat flux from surfaces (land, roof, wall)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%wshf_eb_av(m) = surf_usm_h(l)%wshf_eb_av(m) + surf_usm_h(l)%wshf_eb(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%wshf_eb_av(m) = surf_usm_v(l)%wshf_eb_av(m) +                  &
                                                    surf_usm_v(l)%wshf_eb(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_qsws' )
!
!--             Array of latent heat flux from surfaces (land, roof, wall)
                IF ( horizontal ) THEN
                DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%qsws_av(m) = surf_usm_h(l)%qsws_av(m) + surf_usm_h(l)%qsws(m) * l_v
                ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%qsws_av(m) =  surf_usm_v(l)%qsws_av(m) +                       &
                                                  surf_usm_v(l)%qsws(m) * l_v
                   ENDDO
                ENDIF

            CASE ( 'usm_qsws_veg' )
!
!--             Array of latent heat flux from vegetation surfaces (land, roof, wall)
                IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%qsws_veg_av(m) = surf_usm_h(l)%qsws_veg_av(m) + surf_usm_h(l)%qsws_veg(m)
                ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%qsws_veg_av(m) = surf_usm_v(l)%qsws_veg_av(m) +                &
                                                     surf_usm_v(l)%qsws_veg(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_qsws_liq' )
!
!--             Array of latent heat flux from surfaces with liquid (land, roof, wall)
                IF ( horizontal ) THEN
                DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%qsws_liq_av(m) = surf_usm_h(l)%qsws_liq_av(m) +                         &
                                               surf_usm_h(l)%qsws_liq(m)
                ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%qsws_liq_av(m) = surf_usm_v(l)%qsws_liq_av(m) +                &
                                                     surf_usm_v(l)%qsws_liq(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_wghf' )
!
!--             Array of heat flux from ground (wall, roof, land)
                IF ( horizontal ) THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%wghf_eb_av(m) = surf_usm_h(l)%wghf_eb_av(m) +                        &
                                                 surf_usm_h(l)%wghf_eb(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%wghf_eb_av(m) = surf_usm_v(l)%wghf_eb_av(m) +                  &
                                                    surf_usm_v(l)%wghf_eb(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_wghf_window' )
!
!--             Array of heat flux from window ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%wghf_eb_window_av(m) = surf_usm_h(l)%wghf_eb_window_av(m) +          &
                                                        surf_usm_h(l)%wghf_eb_window(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%wghf_eb_window_av(m) = surf_usm_v(l)%wghf_eb_window_av(m) +    &
                                                           surf_usm_v(l)%wghf_eb_window(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_wghf_green' )
!
!--             Array of heat flux from green ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%wghf_eb_green_av(m) = surf_usm_h(l)%wghf_eb_green_av(m) +            &
                                                       surf_usm_h(l)%wghf_eb_green(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%wghf_eb_green_av(m) = surf_usm_v(l)%wghf_eb_green_av(m) +      &
                                                          surf_usm_v(l)%wghf_eb_green(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_iwghf' )
!
!--             Array of heat flux from indoor ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%iwghf_eb_av(m) = surf_usm_h(l)%iwghf_eb_av(m) + surf_usm_h(l)%iwghf_eb(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%iwghf_eb_av(m) = surf_usm_v(l)%iwghf_eb_av(m) +                &
                                                     surf_usm_v(l)%iwghf_eb(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_iwghf_window' )
!
!--             Array of heat flux from indoor window ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%iwghf_eb_window_av(m) = surf_usm_h(l)%iwghf_eb_window_av(m) +        &
                                                         surf_usm_h(l)%iwghf_eb_window(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%iwghf_eb_window_av(m) = surf_usm_v(l)%iwghf_eb_window_av(m) +  &
                                                            surf_usm_v(l)%iwghf_eb_window(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_t_surf_wall' )
!
!--             Surface temperature for surfaces
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%t_surf_wall_av(m) = surf_usm_h(l)%t_surf_wall_av(m) + t_surf_wall_h(l)%val(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_surf_wall_av(m) = surf_usm_v(l)%t_surf_wall_av(m) +          &
                                                        t_surf_wall_v(l)%val(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_t_surf_window' )
!
!--             Surface temperature for window surfaces
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_surf_window_av(m) = surf_usm_h(l)%t_surf_window_av(m) +            &
                                                       t_surf_window_h(l)%val(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_surf_window_av(m) = surf_usm_v(l)%t_surf_window_av(m) +      &
                                                          t_surf_window_v(l)%val(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_t_surf_green' )
!
!--             Surface temperature for green surfaces
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_surf_green_av(m) = surf_usm_h(l)%t_surf_green_av(m) +              &
                                                      t_surf_green_h(l)%val(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_surf_green_av(m) = surf_usm_v(l)%t_surf_green_av(m) +        &
                                                         t_surf_green_v(l)%val(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_theta_10cm' )
!
!--             Near surface temperature for whole surfaces
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%pt_10cm_av(m) = surf_usm_h(l)%pt_10cm_av(m) +                        &
                                                 surf_usm_h(l)%pt_10cm(m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%pt_10cm_av(m) = surf_usm_v(l)%pt_10cm_av(m) +                  &
                                                    surf_usm_v(l)%pt_10cm(m)
                   ENDDO
                ENDIF

            CASE ( 'usm_t_wall' )
!
!--             Wall temperature for  iwl layer of walls and land
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_wall_av(iwl,m) = surf_usm_h(l)%t_wall_av(iwl,m) +                  &
                                                    t_wall_h(l)%val(iwl,m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_wall_av(iwl,m) = surf_usm_v(l)%t_wall_av(iwl,m) +            &
                                                       t_wall_v(l)%val(iwl,m)
                   ENDDO
                ENDIF

            CASE ( 'usm_t_window' )
!
!--             Window temperature for  iwl layer of walls and land
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_window_av(iwl,m) = surf_usm_h(l)%t_window_av(iwl,m) +              &
                                                         t_window_h(l)%val(iwl,m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_window_av(iwl,m) = surf_usm_v(l)%t_window_av(iwl,m) +        &
                                                         t_window_v(l)%val(iwl,m)
                   ENDDO
                ENDIF

            CASE ( 'usm_t_green' )
!
!--             Green temperature for  iwl layer of walls and land
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_green_av(iwl,m) = surf_usm_h(l)%t_green_av(iwl,m) + t_green_h(l)%val(iwl,m)
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_green_av(iwl,m) = surf_usm_v(l)%t_green_av(iwl,m) +          &
                                                        t_green_v(l)%val(iwl,m)
                   ENDDO
                ENDIF

            CASE ( 'usm_swc' )
!
!--             Soil water content for  iwl layer of walls and land
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%swc_av(iwl,m) = surf_usm_h(l)%swc_av(iwl,m) + swc_h(l)%val(iwl,m)
                   ENDDO
                ELSE
                ENDIF

            CASE DEFAULT
                CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( var ) )

            CASE ( 'usm_wshf' )
!
!--             Array of sensible heat flux from surfaces (land, roof, wall)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%wshf_eb_av(m) = surf_usm_h(l)%wshf_eb_av(m) /                        &
                                                 REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%wshf_eb_av(m) = surf_usm_v(l)%wshf_eb_av(m) /                  &
                                                    REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_qsws' )
!
!--             Array of latent heat flux from surfaces (land, roof, wall)
                IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%qsws_av(m) = surf_usm_h(l)%qsws_av(m) /                                 &
                                           REAL( average_count_3d, kind=wp )
                ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%qsws_av(m) = surf_usm_v(l)%qsws_av(m) /                        &
                                                 REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_qsws_veg' )
!
!--             Array of latent heat flux from vegetation surfaces (land, roof, wall)
                IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%qsws_veg_av(m) = surf_usm_h(l)%qsws_veg_av(m) /                         &
                                               REAL( average_count_3d, kind=wp )
                ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%qsws_veg_av(m) = surf_usm_v(l)%qsws_veg_av(m) /                &
                                                     REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_qsws_liq' )
!
!--             Array of latent heat flux from surfaces with liquid (land, roof, wall)
                IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%qsws_liq_av(m) = surf_usm_h(l)%qsws_liq_av(m) /                         &
                                               REAL( average_count_3d, kind=wp )
                ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%qsws_liq_av(m) = surf_usm_v(l)%qsws_liq_av(m) /                &
                                                     REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_wghf' )
!
!--             Array of heat flux from ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%wghf_eb_av(m) = surf_usm_h(l)%wghf_eb_av(m) /                        &
                                                 REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%wghf_eb_av(m) = surf_usm_v(l)%wghf_eb_av(m) /                  &
                                                    REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_wghf_window' )
!
!--             Array of heat flux from window ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%wghf_eb_window_av(m) = surf_usm_h(l)%wghf_eb_window_av(m) /          &
                                                        REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%wghf_eb_window_av(m) = surf_usm_v(l)%wghf_eb_window_av(m) /    &
                                                           REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_wghf_green' )
!
!--             Array of heat flux from green ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%wghf_eb_green_av(m) = surf_usm_h(l)%wghf_eb_green_av(m) /            &
                                                       REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%wghf_eb_green_av(m) = surf_usm_v(l)%wghf_eb_green_av(m) /      &
                                                          REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_iwghf' )
!
!--             Array of heat flux from indoor ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%iwghf_eb_av(m) = surf_usm_h(l)%iwghf_eb_av(m) /                      &
                                                  REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%iwghf_eb_av(m) = surf_usm_v(l)%iwghf_eb_av(m) /                &
                                                     REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_iwghf_window' )
!
!--             Array of heat flux from indoor window ground (wall, roof, land)
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%iwghf_eb_window_av(m) = surf_usm_h(l)%iwghf_eb_window_av(m) /        &
                                                         REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%iwghf_eb_window_av(m) = surf_usm_v(l)%iwghf_eb_window_av(m) /  &
                                                            REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_t_surf_wall' )
!
!--             Surface temperature for surfaces
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%t_surf_wall_av(m) = surf_usm_h(l)%t_surf_wall_av(m) /                   &
                                                  REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_surf_wall_av(m) = surf_usm_v(l)%t_surf_wall_av(m) /          &
                                                        REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_t_surf_window' )
!
!--             Surface temperature for window surfaces
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_surf_window_av(m) = surf_usm_h(l)%t_surf_window_av(m) /            &
                                                       REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_surf_window_av(m) = surf_usm_v(l)%t_surf_window_av(m) /      &
                                                          REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_t_surf_green' )
!
!--             Surface temperature for green surfaces
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_surf_green_av(m) = surf_usm_h(l)%t_surf_green_av(m) /              &
                                                      REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_surf_green_av(m) = surf_usm_v(l)%t_surf_green_av(m) /        &
                                                         REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_theta_10cm' )
!
!--             Near surface temperature for whole surfaces
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%pt_10cm_av(m) = surf_usm_h(l)%pt_10cm_av(m) /                        &
                                                 REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%pt_10cm_av(m) = surf_usm_v(l)%pt_10cm_av(m) /                  &
                                                    REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF


            CASE ( 'usm_t_wall' )
!
!--             Wall temperature for  iwl layer of walls and land
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_wall_av(iwl,m) = surf_usm_h(l)%t_wall_av(iwl,m) /                  &
                                                    REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_wall_av(iwl,m) = surf_usm_v(l)%t_wall_av(iwl,m) /            &
                                                       REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_t_window' )
!
!--             Window temperature for  iwl layer of walls and land
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_window_av(iwl,m) = surf_usm_h(l)%t_window_av(iwl,m) /              &
                                                      REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_window_av(iwl,m) = surf_usm_v(l)%t_window_av(iwl,m) /        &
                                                         REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_t_green' )
!
!--             Green temperature for  iwl layer of walls and land
                IF ( horizontal )  THEN
                   DO  m = 1, surf_usm_h(l)%ns
                      surf_usm_h(l)%t_green_av(iwl,m) = surf_usm_h(l)%t_green_av(iwl,m) /                &
                                                     REAL( average_count_3d, kind=wp )
                   ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%t_green_av(iwl,m) = surf_usm_v(l)%t_green_av(iwl,m) /          &
                                                        REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF

            CASE ( 'usm_swc' )
!
!--             Soil water content for  iwl layer of walls and land
                IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   surf_usm_h(l)%swc_av(iwl,m) = surf_usm_h(l)%swc_av(iwl,m) /                           &
                                              REAL( average_count_3d, kind=wp )
                ENDDO
                ELSE
                   DO  m = 1, surf_usm_v(l)%ns
                      surf_usm_v(l)%swc_av(iwl,m) = surf_usm_v(l)%swc_av(iwl,m) /                  &
                                                    REAL( average_count_3d, kind=wp )
                   ENDDO
                ENDIF


       END SELECT

    ENDIF

 END SUBROUTINE usm_3d_data_averaging



!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set internal Neumann boundary condition at outer soil grid points for temperature and humidity.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_boundary_condition

    IMPLICIT NONE

    INTEGER(iwp) ::  i      !< grid index x-direction
    INTEGER(iwp) ::  ioff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) ::  j      !< grid index y-direction
    INTEGER(iwp) ::  joff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) ::  k      !< grid index z-direction
    INTEGER(iwp) ::  koff   !< offset index x-direction indicating location of soil grid point
    INTEGER(iwp) ::  l      !< running index surface-orientation
    INTEGER(iwp) ::  m      !< running index surface elements

    DO  l = 0, 1
       koff = surf_usm_h(l)%koff
       DO  m = 1, surf_usm_h(l)%ns
          i = surf_usm_h(l)%i(m)
          j = surf_usm_h(l)%j(m)
          k = surf_usm_h(l)%k(m)
          pt(k+koff,j,i) = pt(k,j,i)
       ENDDO
    ENDDO

    DO  l = 0, 3
       ioff = surf_usm_v(l)%ioff
       joff = surf_usm_v(l)%joff
       DO  m = 1, surf_usm_v(l)%ns
          i = surf_usm_v(l)%i(m)
          j = surf_usm_v(l)%j(m)
          k = surf_usm_v(l)%k(m)
          pt(k,j+joff,i+ioff) = pt(k,j,i)
       ENDDO
    ENDDO

 END SUBROUTINE usm_boundary_condition


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Subroutine checks variables and assigns units.
!> It is called out from subroutine check_parameters.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_check_data_output( variable, unit )

    IMPLICIT NONE

    CHARACTER(LEN=*),INTENT(IN)    ::  variable   !<
    CHARACTER(LEN=*),INTENT(OUT)   ::  unit       !<

    CHARACTER(LEN=2)                              ::  ls            !<

    CHARACTER(LEN=varnamelength)                  ::  var           !< TRIM(variable)

    INTEGER(iwp)                                  ::  i,j,l         !< index

    INTEGER(iwp), PARAMETER                       ::  nl1 = 15      !< number of directional usm variables
    CHARACTER(LEN=varnamelength), DIMENSION(nl1)  ::  varlist1 = &  !< list of directional usm variables
              (/'usm_wshf                      ', &
                'usm_wghf                      ', &
                'usm_wghf_window               ', &
                'usm_wghf_green                ', &
                'usm_iwghf                     ', &
                'usm_iwghf_window              ', &
                'usm_surfz                     ', &
                'usm_surfwintrans              ', &
                'usm_surfcat                   ', &
                'usm_t_surf_wall               ', &
                'usm_t_surf_window             ', &
                'usm_t_surf_green              ', &
                'usm_t_green                   ', &
                'usm_qsws                      ', &
                'usm_theta_10cm                '/)

    INTEGER(iwp), PARAMETER                       ::  nl2 = 3       !< number of directional layer usm variables
    CHARACTER(LEN=varnamelength), DIMENSION(nl2)  ::  varlist2 = &  !< list of directional layer usm variables
              (/'usm_t_wall                    ', &
                'usm_t_window                  ', &
                'usm_t_green                   '/)

    LOGICAL                                       ::  lfound     !< flag if the variable is found


    lfound = .FALSE.

    var = TRIM( variable )

!
!-- Check if variable exists
!-- Directional variables
    DO  i = 1, nl1
       DO  j = 1, nd
          IF ( TRIM( var ) == TRIM( varlist1(i)) // TRIM( dirname(j) ) )  THEN
             lfound = .TRUE.
             EXIT
          ENDIF
          IF ( lfound )  EXIT
       ENDDO
    ENDDO
    IF ( lfound )  GOTO 10
!
!-- Directional layer variables
    DO  i = 1, nl2
       DO  j = 1, nd
          DO  l = nzb_wall, nzt_wall
             WRITE( ls,'(A1,I1)' ) '_', l
             IF ( TRIM( var ) == TRIM( varlist2(i) ) // TRIM( ls ) // TRIM( dirname(j) ) )  THEN
                lfound = .TRUE.
                EXIT
             ENDIF
          ENDDO
          IF ( lfound )  EXIT
       ENDDO
    ENDDO
    IF ( .NOT.  lfound )  THEN
       unit = 'illegal'
       RETURN
    ENDIF
10  CONTINUE

    IF ( var(1:9)  == 'usm_wshf_'  .OR.  var(1:9) == 'usm_wghf_' .OR.                              &
         var(1:16) == 'usm_wghf_window_' .OR. var(1:15) == 'usm_wghf_green_' .OR.                  &
         var(1:10) == 'usm_iwghf_' .OR. var(1:17) == 'usm_iwghf_window_'    .OR.                   &
         var(1:17) == 'usm_surfwintrans_' .OR.                                                     &
         var(1:9)  == 'usm_qsws_'  .OR.  var(1:13)  == 'usm_qsws_veg_'  .OR.                       &
         var(1:13) == 'usm_qsws_liq_' )  THEN
        unit = 'W/m2'
    ELSE IF ( var(1:15) == 'usm_t_surf_wall'   .OR.  var(1:10) == 'usm_t_wall' .OR.                &
              var(1:12) == 'usm_t_window' .OR. var(1:17) == 'usm_t_surf_window' .OR.               &
              var(1:16) == 'usm_t_surf_green'  .OR.                                                &
              var(1:11) == 'usm_t_green' .OR.  var(1:7) == 'usm_swc' .OR.                          &
              var(1:14) == 'usm_theta_10cm' )  THEN
        unit = 'K'
    ELSE IF ( var(1:9) == 'usm_surfz'  .OR.  var(1:11) == 'usm_surfcat' )  THEN
        unit = '1'
    ELSE
        unit = 'illegal'
    ENDIF

 END SUBROUTINE usm_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for urban surface model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_check_parameters

    USE control_parameters,                                                                        &
        ONLY:  bc_pt_b,                                                                            &
               bc_q_b,                                                                             &
               constant_flux_layer,                                                                &
               large_scale_forcing,                                                                &
               lsf_surf,                                                                           &
               topography

    USE netcdf_data_input_mod,                                                                     &
         ONLY:  building_type_f

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< running index, x-dimension
    INTEGER(iwp) ::  j  !< running index, y-dimension

!
!-- Dirichlet boundary conditions are required as the surface fluxes are calculated from the
!-- temperature/humidity gradients in the urban surface model
    IF ( bc_pt_b == 'neumann'   .OR.   bc_q_b == 'neumann' )  THEN
       message_string = 'urban surface model requires setting of bc_pt_b = "dirichlet" and '//     &
                        'bc_q_b  = "dirichlet"'
       CALL message( 'usm_check_parameters', 'PA0590', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( .NOT.  constant_flux_layer )  THEN
       message_string = 'urban surface model requires constant_flux_layer = .TRUE.'
       CALL message( 'usm_check_parameters', 'PA0084', 1, 2, 0, 6, 0 )
    ENDIF

    IF (  .NOT.  radiation )  THEN
       message_string = 'urban surface model requires the radiation model to be switched on'
       CALL message( 'usm_check_parameters', 'PA0084', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Surface forcing has to be disabled for LSF in case of enabled urban surface module
    IF ( large_scale_forcing )  THEN
       lsf_surf = .FALSE.
    ENDIF
!
!-- Topography
    IF ( topography == 'flat' )  THEN
       message_string = 'topography /= "flat" is required when using the urban surface model'
       CALL message( 'usm_check_parameters', 'PA0592', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check if building types are set within a valid range.
    IF ( building_type < LBOUND( building_pars, 2 )  .AND.                                         &
         building_type > UBOUND( building_pars, 2 ) )  THEN
       WRITE( message_string, * ) 'building_type = ', building_type, ' is out of the valid range'
       CALL message( 'usm_check_parameters', 'PA0529', 2, 2, 0, 6, 0 )
    ENDIF
    IF ( building_type_f%from_file )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
             IF ( building_type_f%var(j,i) /= building_type_f%fill  .AND.                          &
           ( building_type_f%var(j,i) < LBOUND( building_pars, 2 )  .OR.                           &
             building_type_f%var(j,i) > UBOUND( building_pars, 2 ) ) )  THEN
                WRITE( message_string, * ) 'building_type = is out of the valid range at (j,i) = ' &
                                           , j, i
                CALL message( 'usm_check_parameters', 'PA0529', 2, 2, myid, 6, 0 )
             ENDIF
          ENDDO
       ENDDO
    ENDIF
 END SUBROUTINE usm_check_parameters


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Output of the 3D-arrays in netCDF and/or AVS format for variables of urban_surface model.
!> It resorts the urban surface module output quantities from surf style indexing into temporary 3D
!> array with indices (i,j,k). It is called from subroutine data_output_3d.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)   ::  variable  !< variable name

    CHARACTER(LEN=varnamelength)   ::  var  !< trimmed variable name

    INTEGER(iwp), INTENT(IN)       ::  av        !< flag if averaged
    INTEGER(iwp), INTENT(IN)       ::  nzb_do    !< lower limit of the data output (usually 0)
    INTEGER(iwp), INTENT(IN)       ::  nzt_do    !< vertical upper limit of the data output (usually nz_do3d)

    INTEGER(iwp)  ::  ids, idsint, idsidx        !<
    INTEGER(iwp)  ::  i, j, k, iwl, istat, l, m  !< running indices
    LOGICAL       ::  horizontal                 !< horizontal upward or downeard facing surface

    LOGICAL, INTENT(OUT)           ::  found     !<

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !< sp - it has to correspond to module data_output_3d
    REAL(wp), DIMENSION(nzb:nzt+1,nys:nyn,nxl:nxr)     ::  temp_pf   !< temp array for urban surface output procedure

    found = .TRUE.
    temp_pf = -1._wp

    ids = -1
    var = TRIM( variable )
    DO i = 0, nd-1
        k = len( TRIM( var ) )
        j = len( TRIM( dirname(i) ) )
        IF ( TRIM( var(k-j+1:k) ) == TRIM( dirname(i) ) )  THEN
            ids = i
            idsint = dirint(ids)
            idsidx = diridx(ids)
            var = var(:k-j)
            EXIT
        ENDIF
    ENDDO
    horizontal = ( ( idsint == iup ) .OR. (idsint == idown ) )
    l = idsidx !< shorter direction index name

    IF ( ids == -1 )  THEN
        var = TRIM( variable )
    ENDIF
    IF ( var(1:11) == 'usm_t_wall_'  .AND.  len( TRIM( var ) ) >= 12 )  THEN
!
!--     Wall layers
        READ( var(12:12), '(I1)', iostat = istat ) iwl
        IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
            var = var(1:10)
        ENDIF
    ENDIF
    IF ( var(1:13) == 'usm_t_window_'  .AND.  len( TRIM( var ) ) >= 14 )  THEN
!
!--     Window layers
        READ( var(14:14), '(I1)', iostat = istat ) iwl
        IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
            var = var(1:12)
        ENDIF
    ENDIF
    IF ( var(1:12) == 'usm_t_green_'  .AND.  len( TRIM( var ) ) >= 13 )  THEN
!
!--     Green layers
        READ( var(13:13), '(I1)', iostat = istat ) iwl
        IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
            var = var(1:11)
        ENDIF
    ENDIF
    IF ( var(1:8) == 'usm_swc_'  .AND.  len( TRIM( var ) ) >= 9 )  THEN
!
!--     Green layers soil water content
        READ( var(9:9), '(I1)', iostat = istat ) iwl
        IF ( istat == 0  .AND.  iwl >= nzb_wall  .AND.  iwl <= nzt_wall )  THEN
            var = var(1:7)
        ENDIF
    ENDIF

    SELECT CASE ( TRIM( var ) )

      CASE ( 'usm_surfz' )
!
!--       Array of surface height (z)
          IF ( horizontal )  THEN
             DO  m = 1, surf_usm_h(l)%ns
                i = surf_usm_h(l)%i(m)
                j = surf_usm_h(l)%j(m)
                k = surf_usm_h(l)%k(m)
                temp_pf(0,j,i) = MAX( temp_pf(0,j,i), REAL( k, KIND = wp) )
             ENDDO
          ELSE
             DO  m = 1, surf_usm_v(l)%ns
                i = surf_usm_v(l)%i(m)
                j = surf_usm_v(l)%j(m)
                k = surf_usm_v(l)%k(m)
                temp_pf(0,j,i) = MAX( temp_pf(0,j,i), REAL( k, KIND = wp) + 1.0_sp )
             ENDDO
          ENDIF

      CASE ( 'usm_surfcat' )
!
!--       Surface category
          IF ( horizontal )  THEN
             DO  m = 1, surf_usm_h(l)%ns
                i = surf_usm_h(l)%i(m)
                j = surf_usm_h(l)%j(m)
                k = surf_usm_h(l)%k(m)
                temp_pf(k,j,i) = surf_usm_h(l)%surface_types(m)
             ENDDO
          ELSE
             DO  m = 1, surf_usm_v(l)%ns
                i = surf_usm_v(l)%i(m)
                j = surf_usm_v(l)%j(m)
                k = surf_usm_v(l)%k(m)
                temp_pf(k,j,i) = surf_usm_v(l)%surface_types(m)
             ENDDO
          ENDIF

      CASE ( 'usm_surfwintrans' )
!
!--       Transmissivity window tiles
          IF ( horizontal )  THEN
             DO  m = 1, surf_usm_h(l)%ns
                i = surf_usm_h(l)%i(m)
                j = surf_usm_h(l)%j(m)
                k = surf_usm_h(l)%k(m)
                temp_pf(k,j,i) = surf_usm_h(l)%transmissivity(m)
             ENDDO
          ELSE
             DO  m = 1, surf_usm_v(l)%ns
                i = surf_usm_v(l)%i(m)
                j = surf_usm_v(l)%j(m)
                k = surf_usm_v(l)%k(m)
                temp_pf(k,j,i) = surf_usm_v(l)%transmissivity(m)
             ENDDO
          ENDIF

      CASE ( 'usm_wshf' )
!
!--       Array of sensible heat flux from surfaces
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%wshf_eb(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%wshf_eb(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%wshf_eb_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%wshf_eb_av(m)
                ENDDO
             ENDIF
          ENDIF


      CASE ( 'usm_qsws' )
!
!--       Array of latent heat flux from surfaces
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%qsws(m) * l_v
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%qsws(m) * l_v
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%qsws_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%qsws_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_qsws_veg' )
!
!--       Array of latent heat flux from vegetation surfaces
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%qsws_veg(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%qsws_veg(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%qsws_veg_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%qsws_veg_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_qsws_liq' )
!
!--       Array of latent heat flux from surfaces with liquid
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%qsws_liq(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%qsws_liq(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%qsws_liq_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%qsws_liq_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_wghf' )
!
!--       Array of heat flux from ground (land, wall, roof)
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%wghf_eb(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%wghf_eb_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_wghf_window' )
!
!--       Array of heat flux from window ground (land, wall, roof)
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%wghf_eb_window(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_window(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%wghf_eb_window_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_window_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_wghf_green' )
!
!--       Array of heat flux from green ground (land, wall, roof)
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%wghf_eb_green(m)
                ENDDO
             ELSE
                l = idsidx
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_green(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%wghf_eb_green_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%wghf_eb_green_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_iwghf' )
!
!--       Array of heat flux from indoor ground (land, wall, roof)
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%iwghf_eb(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%iwghf_eb_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_iwghf_window' )
!
!--       Array of heat flux from indoor window ground (land, wall, roof)
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%iwghf_eb_window(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_window(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%iwghf_eb_window_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%iwghf_eb_window_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_t_surf_wall' )
!
!--       Surface temperature for surfaces
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = t_surf_wall_h(l)%val(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = t_surf_wall_v(l)%val(m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%t_surf_wall_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%t_surf_wall_av(m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_t_surf_window' )
!
!--       Surface temperature for window surfaces
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = t_surf_window_h(l)%val(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = t_surf_window_v(l)%val(m)
                ENDDO
             ENDIF

          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%t_surf_window_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%t_surf_window_av(m)
                ENDDO

             ENDIF

          ENDIF

      CASE ( 'usm_t_surf_green' )
!
!--       Surface temperature for green surfaces
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = t_surf_green_h(l)%val(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = t_surf_green_v(l)%val(m)
                ENDDO
             ENDIF

          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%t_surf_green_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%t_surf_green_av(m)
                ENDDO

             ENDIF

          ENDIF

      CASE ( 'usm_theta_10cm' )
!
!--       Near surface temperature for whole surfaces
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%pt_10cm(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%pt_10cm(m)
                ENDDO
             ENDIF


          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%pt_10cm_av(m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%pt_10cm_av(m)
                ENDDO

              ENDIF
          ENDIF

      CASE ( 'usm_t_wall' )
!
!--       Wall temperature for  iwl layer of walls and land
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = t_wall_h(l)%val(iwl,m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = t_wall_v(l)%val(iwl,m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%t_wall_av(iwl,m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%t_wall_av(iwl,m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_t_window' )
!
!--       Window temperature for iwl layer of walls and land
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = t_window_h(l)%val(iwl,m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = t_window_v(l)%val(iwl,m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%t_window_av(iwl,m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%t_window_av(iwl,m)
                ENDDO
             ENDIF
          ENDIF

      CASE ( 'usm_t_green' )
!
!--       Green temperature for  iwl layer of walls and land
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = t_green_h(l)%val(iwl,m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = t_green_v(l)%val(iwl,m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%t_green_av(iwl,m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%t_green_av(iwl,m)
                ENDDO
             ENDIF
          ENDIF

          CASE ( 'usm_swc' )
!
!--       Soil water content for  iwl layer of walls and land
          IF ( av == 0 )  THEN
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = swc_h(l)%val(iwl,m)
                ENDDO
             ENDIF
          ELSE
             IF ( horizontal )  THEN
                DO  m = 1, surf_usm_h(l)%ns
                   i = surf_usm_h(l)%i(m)
                   j = surf_usm_h(l)%j(m)
                   k = surf_usm_h(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_h(l)%swc_av(iwl,m)
                ENDDO
             ELSE
                DO  m = 1, surf_usm_v(l)%ns
                   i = surf_usm_v(l)%i(m)
                   j = surf_usm_v(l)%j(m)
                   k = surf_usm_v(l)%k(m)
                   temp_pf(k,j,i) = surf_usm_v(l)%swc_av(iwl,m)
                ENDDO
             ENDIF
          ENDIF


      CASE DEFAULT
          found = .FALSE.
          RETURN
    END SELECT

!
!-- Rearrange dimensions for NetCDF output
!-- FIXME: this may generate FPE overflow upon conversion from DP to SP
    DO  j = nys, nyn
        DO  i = nxl, nxr
            DO  k = nzb_do, nzt_do
                local_pf(i,j,k) = temp_pf(k,j,i)
            ENDDO
        ENDDO
    ENDDO

 END SUBROUTINE usm_data_output_3d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Soubroutine defines appropriate grid for netcdf variables.
!> It is called out from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  ::  variable  !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_x    !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_y    !<
    CHARACTER(LEN=*), INTENT(OUT) ::  grid_z    !<

    CHARACTER(LEN=varnamelength)  ::  var  !<

    LOGICAL, INTENT(OUT)  ::  found  !<

    var = TRIM( variable )
    IF ( var(1:9) == 'usm_wshf_'  .OR.  var(1:9) == 'usm_wghf_'  .OR.                              &
         var(1:16) == 'usm_wghf_window_'  .OR. var(1:15) == 'usm_wghf_green_' .OR.                 &
         var(1:10) == 'usm_iwghf_'  .OR. var(1:17) == 'usm_iwghf_window_' .OR.                     &
         var(1:9) == 'usm_qsws_'  .OR.  var(1:13) == 'usm_qsws_veg_'  .OR.                         &
         var(1:13) == 'usm_qsws_liq_' .OR.                                                         &
         var(1:15) == 'usm_t_surf_wall'  .OR.  var(1:10) == 'usm_t_wall'  .OR.                     &
         var(1:17) == 'usm_t_surf_window'  .OR.  var(1:12) == 'usm_t_window'  .OR.                 &
         var(1:16) == 'usm_t_surf_green'  .OR. var(1:11) == 'usm_t_green' .OR.                     &
         var(1:15) == 'usm_theta_10cm' .OR.                                                        &
         var(1:9) == 'usm_surfz'  .OR.  var(1:11) == 'usm_surfcat'  .OR.                           &
         var(1:16) == 'usm_surfwintrans'  .OR. var(1:7) == 'usm_swc' ) THEN

        found = .TRUE.
        grid_x = 'x'
        grid_y = 'y'
        grid_z = 'zu'
    ELSE
        found  = .FALSE.
        grid_x = 'none'
        grid_y = 'none'
        grid_z = 'none'
    ENDIF

 END SUBROUTINE usm_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the wall surface model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_init_wall_heat_model

    IMPLICIT NONE

    INTEGER(iwp) ::  k, l, m  !< running indices

    IF ( debug_output )  CALL debug_message( 'usm_init_wall_heat_model', 'start' )

!
!-- Calculate wall and window grid spacings. Wall temperature is defined at the center of the
!-- wall layers.
!-- First for horizontal surfaces:
    DO  l = 0, 1
       DO  m = 1, surf_usm_h(l)%ns

          surf_usm_h(l)%dz_wall(nzb_wall,m) = surf_usm_h(l)%zw(nzb_wall,m)
          DO k = nzb_wall+1, nzt_wall
             surf_usm_h(l)%dz_wall(k,m) = surf_usm_h(l)%zw(k,m) - surf_usm_h(l)%zw(k-1,m)
          ENDDO
          surf_usm_h(l)%dz_window(nzb_wall,m) = surf_usm_h(l)%zw_window(nzb_wall,m)
          DO  k = nzb_wall+1, nzt_wall
             surf_usm_h(l)%dz_window(k,m) = surf_usm_h(l)%zw_window(k,m) - surf_usm_h(l)%zw_window(k-1,m)
          ENDDO

          surf_usm_h(l)%dz_wall(nzt_wall+1,m) = surf_usm_h(l)%dz_wall(nzt_wall,m)

          DO k = nzb_wall, nzt_wall-1
            surf_usm_h(l)%dz_wall_stag(k,m) = 0.5 * ( surf_usm_h(l)%dz_wall(k+1,m) +                        &
                                                     surf_usm_h(l)%dz_wall(k,m) )
          ENDDO
          surf_usm_h(l)%dz_wall_stag(nzt_wall,m) = surf_usm_h(l)%dz_wall(nzt_wall,m)

          surf_usm_h(l)%dz_window(nzt_wall+1,m) = surf_usm_h(l)%dz_window(nzt_wall,m)

          DO  k = nzb_wall, nzt_wall-1
             surf_usm_h(l)%dz_window_stag(k,m) = 0.5 * ( surf_usm_h(l)%dz_window(k+1,m) +                   &
                                                      surf_usm_h(l)%dz_window(k,m) )
          ENDDO
          surf_usm_h(l)%dz_window_stag(nzt_wall,m) = surf_usm_h(l)%dz_window(nzt_wall,m)

          IF (surf_usm_h(l)%green_type_roof(m) == 2.0_wp )  THEN
!
!--          Extensive green roof
!--          Set ratio of substrate layer thickness, soil-type and LAI
             soil_type = 3
             surf_usm_h(l)%lai(m) = 2.0_wp

             surf_usm_h(l)%zw_green(nzb_wall,m)   = 0.05_wp
             surf_usm_h(l)%zw_green(nzb_wall+1,m) = 0.10_wp
             surf_usm_h(l)%zw_green(nzb_wall+2,m) = 0.15_wp
             surf_usm_h(l)%zw_green(nzb_wall+3,m) = 0.20_wp
          ELSE
!
!--          Intensiv green roof
!--          Set ratio of substrate layer thickness, soil-type and LAI
             soil_type = 6
             surf_usm_h(l)%lai(m) = 4.0_wp

             surf_usm_h(l)%zw_green(nzb_wall,m)   = 0.05_wp
             surf_usm_h(l)%zw_green(nzb_wall+1,m) = 0.10_wp
             surf_usm_h(l)%zw_green(nzb_wall+2,m) = 0.40_wp
             surf_usm_h(l)%zw_green(nzb_wall+3,m) = 0.80_wp
          ENDIF

          surf_usm_h(l)%dz_green(nzb_wall,m) = surf_usm_h(l)%zw_green(nzb_wall,m)
          DO k = nzb_wall+1, nzt_wall
              surf_usm_h(l)%dz_green(k,m) = surf_usm_h(l)%zw_green(k,m) - surf_usm_h(l)%zw_green(k-1,m)
          ENDDO
          surf_usm_h(l)%dz_green(nzt_wall+1,m) = surf_usm_h(l)%dz_green(nzt_wall,m)

          DO k = nzb_wall, nzt_wall-1
              surf_usm_h(l)%dz_green_stag(k,m) = 0.5 * ( surf_usm_h(l)%dz_green(k+1,m) +                    &
                                                      surf_usm_h(l)%dz_green(k,m) )
          ENDDO
          surf_usm_h(l)%dz_green_stag(nzt_wall,m) = surf_usm_h(l)%dz_green(nzt_wall,m)

         IF ( alpha_vangenuchten == 9999999.9_wp )  THEN
            alpha_vangenuchten = soil_pars(0,soil_type)
         ENDIF

         IF ( l_vangenuchten == 9999999.9_wp )  THEN
            l_vangenuchten = soil_pars(1,soil_type)
         ENDIF

         IF ( n_vangenuchten == 9999999.9_wp )  THEN
            n_vangenuchten = soil_pars(2,soil_type)
         ENDIF

         IF ( hydraulic_conductivity == 9999999.9_wp )  THEN
            hydraulic_conductivity = soil_pars(3,soil_type)
         ENDIF

         IF ( saturation_moisture == 9999999.9_wp )  THEN
            saturation_moisture = m_soil_pars(0,soil_type)
         ENDIF

         IF ( field_capacity == 9999999.9_wp )  THEN
            field_capacity = m_soil_pars(1,soil_type)
         ENDIF

         IF ( wilting_point == 9999999.9_wp )  THEN
            wilting_point = m_soil_pars(2,soil_type)
         ENDIF

         IF ( residual_moisture == 9999999.9_wp )  THEN
            residual_moisture = m_soil_pars(3,soil_type)
         ENDIF

         DO  k = nzb_wall, nzt_wall+1
            swc_h(l)%val(k,m) = field_capacity
            rootfr_h(l)%val(k,m) = 0.5_wp
            surf_usm_h(l)%alpha_vg_green(m)      = alpha_vangenuchten
            surf_usm_h(l)%l_vg_green(m)          = l_vangenuchten
            surf_usm_h(l)%n_vg_green(m)          = n_vangenuchten
            surf_usm_h(l)%gamma_w_green_sat(k,m) = hydraulic_conductivity
            swc_sat_h(l)%val(k,m)                = saturation_moisture
            fc_h(l)%val(k,m)                     = field_capacity
            wilt_h(l)%val(k,m)                   = wilting_point
            swc_res_h(l)%val(k,m)                = residual_moisture
         ENDDO

       ENDDO

       surf_usm_h(l)%ddz_wall        = 1.0_wp / surf_usm_h(l)%dz_wall
       surf_usm_h(l)%ddz_wall_stag   = 1.0_wp / surf_usm_h(l)%dz_wall_stag
       surf_usm_h(l)%ddz_window      = 1.0_wp / surf_usm_h(l)%dz_window
       surf_usm_h(l)%ddz_window_stag = 1.0_wp / surf_usm_h(l)%dz_window_stag
       surf_usm_h(l)%ddz_green       = 1.0_wp / surf_usm_h(l)%dz_green
       surf_usm_h(l)%ddz_green_stag  = 1.0_wp / surf_usm_h(l)%dz_green_stag
    ENDDO
!
!-- For vertical surfaces
    DO  l = 0, 3
       DO  m = 1, surf_usm_v(l)%ns
          surf_usm_v(l)%dz_wall(nzb_wall,m) = surf_usm_v(l)%zw(nzb_wall,m)
          DO  k = nzb_wall+1, nzt_wall
             surf_usm_v(l)%dz_wall(k,m) = surf_usm_v(l)%zw(k,m) - surf_usm_v(l)%zw(k-1,m)
          ENDDO
          surf_usm_v(l)%dz_window(nzb_wall,m) = surf_usm_v(l)%zw_window(nzb_wall,m)
          DO  k = nzb_wall+1, nzt_wall
             surf_usm_v(l)%dz_window(k,m) = surf_usm_v(l)%zw_window(k,m) -                         &
                                            surf_usm_v(l)%zw_window(k-1,m)
          ENDDO
          surf_usm_v(l)%dz_green(nzb_wall,m) = surf_usm_v(l)%zw_green(nzb_wall,m)
          DO  k = nzb_wall+1, nzt_wall
             surf_usm_v(l)%dz_green(k,m) = surf_usm_v(l)%zw_green(k,m) -                           &
                                           surf_usm_v(l)%zw_green(k-1,m)
          ENDDO

          surf_usm_v(l)%dz_wall(nzt_wall+1,m) = surf_usm_v(l)%dz_wall(nzt_wall,m)

          DO  k = nzb_wall, nzt_wall-1
             surf_usm_v(l)%dz_wall_stag(k,m) = 0.5 * ( surf_usm_v(l)%dz_wall(k+1,m) +              &
                                                       surf_usm_v(l)%dz_wall(k,m) )
          ENDDO
          surf_usm_v(l)%dz_wall_stag(nzt_wall,m) = surf_usm_v(l)%dz_wall(nzt_wall,m)
          surf_usm_v(l)%dz_window(nzt_wall+1,m)  = surf_usm_v(l)%dz_window(nzt_wall,m)

          DO  k = nzb_wall, nzt_wall-1
             surf_usm_v(l)%dz_window_stag(k,m) = 0.5 * ( surf_usm_v(l)%dz_window(k+1,m) +          &
                                                         surf_usm_v(l)%dz_window(k,m) )
          ENDDO
          surf_usm_v(l)%dz_window_stag(nzt_wall,m) = surf_usm_v(l)%dz_window(nzt_wall,m)
          surf_usm_v(l)%dz_green(nzt_wall+1,m)     = surf_usm_v(l)%dz_green(nzt_wall,m)

          DO  k = nzb_wall, nzt_wall-1
             surf_usm_v(l)%dz_green_stag(k,m) = 0.5 * ( surf_usm_v(l)%dz_green(k+1,m) +            &
                                                        surf_usm_v(l)%dz_green(k,m) )
          ENDDO
          surf_usm_v(l)%dz_green_stag(nzt_wall,m) = surf_usm_v(l)%dz_green(nzt_wall,m)
       ENDDO
       surf_usm_v(l)%ddz_wall        = 1.0_wp / surf_usm_v(l)%dz_wall
       surf_usm_v(l)%ddz_wall_stag   = 1.0_wp / surf_usm_v(l)%dz_wall_stag
       surf_usm_v(l)%ddz_window      = 1.0_wp / surf_usm_v(l)%dz_window
       surf_usm_v(l)%ddz_window_stag = 1.0_wp / surf_usm_v(l)%dz_window_stag
       surf_usm_v(l)%ddz_green       = 1.0_wp / surf_usm_v(l)%dz_green
       surf_usm_v(l)%ddz_green_stag  = 1.0_wp / surf_usm_v(l)%dz_green_stag
    ENDDO


    IF ( debug_output )  CALL debug_message( 'usm_init_wall_heat_model', 'end' )

 END SUBROUTINE usm_init_wall_heat_model


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the urban surface model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_init

    USE arrays_3d,                                                                                 &
        ONLY:  zw

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  albedo_type_f,                                                                      &
               building_pars_f,                                                                    &
               building_surface_pars_f,                                                            &
               building_type_f,                                                                    &
               terrain_height_f

    IMPLICIT NONE

    INTEGER(iwp) ::  i                 !< loop index x-dirction
    INTEGER(iwp) ::  ind_alb_green     !< index in input list for green albedo
    INTEGER(iwp) ::  ind_alb_wall      !< index in input list for wall albedo
    INTEGER(iwp) ::  ind_alb_win       !< index in input list for window albedo
    INTEGER(iwp) ::  ind_emis_wall     !< index in input list for wall emissivity
    INTEGER(iwp) ::  ind_emis_green    !< index in input list for green emissivity
    INTEGER(iwp) ::  ind_emis_win      !< index in input list for window emissivity
    INTEGER(iwp) ::  ind_green_frac_w  !< index in input list for green fraction on wall
    INTEGER(iwp) ::  ind_green_frac_r  !< index in input list for green fraction on roof
    INTEGER(iwp) ::  ind_hc1           !< index in input list for heat capacity at first wall layer
    INTEGER(iwp) ::  ind_hc1_win       !< index in input list for heat capacity at first window layer
    INTEGER(iwp) ::  ind_hc2           !< index in input list for heat capacity at second wall layer
    INTEGER(iwp) ::  ind_hc2_win       !< index in input list for heat capacity at second window layer
    INTEGER(iwp) ::  ind_hc3           !< index in input list for heat capacity at third wall layer
    INTEGER(iwp) ::  ind_hc3_win       !< index in input list for heat capacity at third window layer
    INTEGER(iwp) ::  ind_hc4           !< index in input list for heat capacity at fourth wall layer
    INTEGER(iwp) ::  ind_hc4_win       !< index in input list for heat capacity at fourth window layer
    INTEGER(iwp) ::  ind_lai_r         !< index in input list for LAI on roof
    INTEGER(iwp) ::  ind_lai_w         !< index in input list for LAI on wall
    INTEGER(iwp) ::  ind_tc1           !< index in input list for thermal conductivity at first wall layer
    INTEGER(iwp) ::  ind_tc1_win       !< index in input list for thermal conductivity at first window layer
    INTEGER(iwp) ::  ind_tc2           !< index in input list for thermal conductivity at second wall layer
    INTEGER(iwp) ::  ind_tc2_win       !< index in input list for thermal conductivity at second window layer
    INTEGER(iwp) ::  ind_tc3           !< index in input list for thermal conductivity at third wall layer
    INTEGER(iwp) ::  ind_tc3_win       !< index in input list for thermal conductivity at third window layer
    INTEGER(iwp) ::  ind_tc4           !< index in input list for thermal conductivity at fourth wall layer
    INTEGER(iwp) ::  ind_tc4_win       !< index in input list for thermal conductivity at fourth window layer
    INTEGER(iwp) ::  ind_thick_1       !< index in input list for thickness of first wall layer
    INTEGER(iwp) ::  ind_thick_1_win   !< index in input list for thickness of first window layer
    INTEGER(iwp) ::  ind_thick_2       !< index in input list for thickness of second wall layer
    INTEGER(iwp) ::  ind_thick_2_win   !< index in input list for thickness of second window layer
    INTEGER(iwp) ::  ind_thick_3       !< index in input list for thickness of third wall layer
    INTEGER(iwp) ::  ind_thick_3_win   !< index in input list for thickness of third window layer
    INTEGER(iwp) ::  ind_thick_4       !< index in input list for thickness of fourth wall layer
    INTEGER(iwp) ::  ind_thick_4_win   !< index in input list for thickness of fourth window layer
    INTEGER(iwp) ::  ind_trans         !< index in input list for window transmissivity
    INTEGER(iwp) ::  ind_wall_frac     !< index in input list for wall fraction
    INTEGER(iwp) ::  ind_win_frac      !< index in input list for window fraction
    INTEGER(iwp) ::  ind_z0            !< index in input list for z0
    INTEGER(iwp) ::  ind_z0qh          !< index in input list for z0h / z0q
    INTEGER(iwp) ::  is                !< loop index input surface element
    INTEGER(iwp) ::  j                 !< loop index y-dirction
    INTEGER(iwp) ::  k                 !< loop index z-dirction
    INTEGER(iwp) ::  l                 !< loop index surface orientation
    INTEGER(iwp) ::  m                 !< loop index surface element
    INTEGER(iwp) ::  st                !< dummy

    LOGICAL      ::  relative_fractions_corrected  !< flag indicating if relative surface fractions require normalization

    REAL(wp)     ::  c, tin, twin          !<
    REAL(wp)     ::  ground_floor_level_l  !< local height of ground floor level
    REAL(wp)     ::  sum_frac              !< sum of the relative material fractions at a surface element
    REAL(wp)     ::  z_agl                 !< height of the surface element above terrain

    IF ( debug_output )  CALL debug_message( 'usm_init', 'start' )

    CALL cpu_log( log_point_s(78), 'usm_init', 'start' )
!
!-- Surface forcing has to be disabled for LSF in case of enabled urban surface module
    IF ( large_scale_forcing )  THEN
        lsf_surf = .FALSE.
    ENDIF
!
!-- Calculate constant values
    d_roughness_concrete = 1.0_wp / roughness_concrete
!
!-- Flag surface elements belonging to the ground floor level. Therefore, use terrain height array
!-- from file, if available. This flag is later used to control initialization of surface attributes.
!-- Todo: for the moment disable initialization of building roofs with ground-floor-level properties.
    DO  l = 0, 1
      surf_usm_h(l)%ground_level = .FALSE.
    ENDDO

    DO  l = 0, 3
       surf_usm_v(l)%ground_level = .FALSE.
       DO  m = 1, surf_usm_v(l)%ns
          i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
          j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
          k = surf_usm_v(l)%k(m)
!
!--       Determine local ground level. Level 1 - default value, level 2 - initialization according
!--       to building type, level 3 - initialization from value read from file.
          ground_floor_level_l = ground_floor_level

          IF ( building_type_f%from_file )  THEN
              ground_floor_level_l = building_pars(ind_gflh,building_type_f%var(j,i))
          ENDIF

          IF ( building_pars_f%from_file )  THEN
             IF ( building_pars_f%pars_xy(ind_gflh,j,i) /=  building_pars_f%fill )                 &
                ground_floor_level_l = building_pars_f%pars_xy(ind_gflh,j,i)
          ENDIF
!
!--       Determine height of surface element above ground level. Please note, the height of a
!--       surface element is determined with respect to its height above ground of the reference
!--       grid point in the atmosphere. Therefore, substract the offset values when assessing the
!--       terrain height.
          IF ( terrain_height_f%from_file )  THEN
             z_agl = zw(k) - terrain_height_f%var(j-surf_usm_v(l)%joff, i-surf_usm_v(l)%ioff)
          ELSE
             z_agl = zw(k)
          ENDIF
!
!--       Set flag for ground level
          IF ( z_agl <= ground_floor_level_l ) surf_usm_v(l)%ground_level(m) = .TRUE.

       ENDDO
    ENDDO
!
!-- Initialization of resistances.
    DO  l = 0, 1
       DO  m = 1, surf_usm_h(l)%ns
          surf_usm_h(l)%r_a(m)        = 50.0_wp
          surf_usm_h(l)%r_a_green(m)  = 50.0_wp
          surf_usm_h(l)%r_a_window(m) = 50.0_wp
       ENDDO
    ENDDO
    DO  l = 0, 3
       DO  m = 1, surf_usm_v(l)%ns
          surf_usm_v(l)%r_a(m)        = 50.0_wp
          surf_usm_v(l)%r_a_green(m)  = 50.0_wp
          surf_usm_v(l)%r_a_window(m) = 50.0_wp
       ENDDO
    ENDDO

!
!-- Map values onto horizontal elemements
    DO  l = 0, 1
       DO  m = 1, surf_usm_h(l)%ns
          surf_usm_h(l)%r_canopy(m)     = 200.0_wp !< canopy_resistance
          surf_usm_h(l)%r_canopy_min(m) = 200.0_wp !< min_canopy_resistance
          surf_usm_h(l)%g_d(m)          = 0.0_wp   !< canopy_resistance_coefficient
       ENDDO
    ENDDO
!
!-- Map values onto vertical elements, even though this does not make much sense.
    DO  l = 0, 3
       DO  m = 1, surf_usm_v(l)%ns
          surf_usm_v(l)%r_canopy(m)     = 200.0_wp !< canopy_resistance
          surf_usm_v(l)%r_canopy_min(m) = 200.0_wp !< min_canopy_resistance
          surf_usm_v(l)%g_d(m)          = 0.0_wp   !< canopy_resistance_coefficient
       ENDDO
    ENDDO

!
!--  Initialize urban-type surface attribute. According to initialization in land-surface model,
!--  follow a 3-level approach.
!--  Level 1 - initialization via default attributes
     DO  l = 0, 1
        DO  m = 1, surf_usm_h(l)%ns
!
!--        Now, all horizontal surfaces are roof surfaces (?)
           surf_usm_h(l)%isroof_surf(m)   = .TRUE.
           surf_usm_h(l)%surface_types(m) = roof_category         !< default category for root surface
!
!--        In order to distinguish between ground floor level and above-ground-floor level surfaces,
!--        set input indices.

           ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl,                     &
                                     surf_usm_h(l)%ground_level(m) )
           ind_lai_r        = MERGE( ind_lai_r_gfl, ind_lai_r_agfl, surf_usm_h(l)%ground_level(m) )
           ind_z0           = MERGE( ind_z0_gfl, ind_z0_agfl, surf_usm_h(l)%ground_level(m) )
           ind_z0qh         = MERGE( ind_z0qh_gfl, ind_z0qh_agfl, surf_usm_h(l)%ground_level(m) )
!
!--        Store building type and its name on each surface element
           surf_usm_h(l)%building_type(m)      = building_type
           surf_usm_h(l)%building_type_name(m) = building_type_name(building_type)
!
!--        Initialize relatvie wall- (0), green- (1) and window (2) fractions
           surf_usm_h(l)%frac(m,ind_veg_wall)  = building_pars(ind_wall_frac_r,building_type)
           surf_usm_h(l)%frac(m,ind_pav_green) = building_pars(ind_green_frac_r,building_type)
           surf_usm_h(l)%frac(m,ind_wat_win)   = building_pars(ind_win_frac_r,building_type)
           surf_usm_h(l)%lai(m)                = building_pars(ind_lai_r,building_type)

           surf_usm_h(l)%rho_c_wall(nzb_wall,m)        = building_pars(ind_hc1_wall_r,building_type)
           surf_usm_h(l)%rho_c_wall(nzb_wall+1,m)      = building_pars(ind_hc2_wall_r,building_type)
           surf_usm_h(l)%rho_c_wall(nzb_wall+2,m)      = building_pars(ind_hc3_wall_r,building_type)
           surf_usm_h(l)%rho_c_wall(nzb_wall+3,m)      = building_pars(ind_hc4_wall_r,building_type)
           surf_usm_h(l)%lambda_h(nzb_wall,m)          = building_pars(ind_tc1_wall_r,building_type)
           surf_usm_h(l)%lambda_h(nzb_wall+1,m)        = building_pars(ind_tc2_wall_r,building_type)
           surf_usm_h(l)%lambda_h(nzb_wall+2,m)        = building_pars(ind_tc3_wall_r,building_type)
           surf_usm_h(l)%lambda_h(nzb_wall+3,m)        = building_pars(ind_tc4_wall_r,building_type)
           surf_usm_h(l)%rho_c_green(nzb_wall,m)       = rho_c_soil !building_pars(ind_hc1_wall_r,building_type)
           surf_usm_h(l)%rho_c_green(nzb_wall+1,m)     = rho_c_soil !building_pars(ind_hc1_wall_r,building_type)
           surf_usm_h(l)%rho_c_green(nzb_wall+2,m)     = rho_c_soil !building_pars(ind_hc2_wall_r,building_type)
           surf_usm_h(l)%rho_c_green(nzb_wall+3,m)     = rho_c_soil !building_pars(ind_hc3_wall_r,building_type)
           surf_usm_h(l)%lambda_h_green(nzb_wall,m)    = lambda_h_green_sm !building_pars(ind_tc1_wall_r,building_type)
           surf_usm_h(l)%lambda_h_green(nzb_wall+1,m)  = lambda_h_green_sm !building_pars(ind_tc1_wall_r,building_type)
           surf_usm_h(l)%lambda_h_green(nzb_wall+2,m)  = lambda_h_green_sm !building_pars(ind_tc2_wall_r,building_type)
           surf_usm_h(l)%lambda_h_green(nzb_wall+3,m)  = lambda_h_green_sm !building_pars(ind_tc3_wall_r,building_type)
           surf_usm_h(l)%rho_c_window(nzb_wall,m)      = building_pars(ind_hc1_win_r,building_type)
           surf_usm_h(l)%rho_c_window(nzb_wall+1,m)    = building_pars(ind_hc2_win_r,building_type)
           surf_usm_h(l)%rho_c_window(nzb_wall+2,m)    = building_pars(ind_hc3_win_r,building_type)
           surf_usm_h(l)%rho_c_window(nzb_wall+3,m)    = building_pars(ind_hc4_win_r,building_type)
           surf_usm_h(l)%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1_win_r,building_type)
           surf_usm_h(l)%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc2_win_r,building_type)
           surf_usm_h(l)%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc3_win_r,building_type)
           surf_usm_h(l)%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc4_win_r,building_type)

           surf_usm_h(l)%target_temp_summer(m)  = building_pars(ind_indoor_target_temp_summer,building_type)
           surf_usm_h(l)%target_temp_winter(m)  = building_pars(ind_indoor_target_temp_winter,building_type)
!
!--        Emissivity of wall-, green- and window fraction
           surf_usm_h(l)%emissivity(m,ind_veg_wall)  = building_pars(ind_emis_wall_r,building_type)
           surf_usm_h(l)%emissivity(m,ind_pav_green) = building_pars(ind_emis_green_r,building_type)
           surf_usm_h(l)%emissivity(m,ind_wat_win)   = building_pars(ind_emis_win_r,building_type)

           surf_usm_h(l)%transmissivity(m)           = building_pars(ind_trans_r,building_type)

           surf_usm_h(l)%z0(m)                       = building_pars(ind_z0,building_type)
           surf_usm_h(l)%z0h(m)                      = building_pars(ind_z0qh,building_type)
           surf_usm_h(l)%z0q(m)                      = building_pars(ind_z0qh,building_type)
!
!--        Albedo type for wall fraction, green fraction, window fraction
           surf_usm_h(l)%albedo_type(m,ind_veg_wall)  = INT( building_pars(ind_alb_wall_r,building_type) )
           surf_usm_h(l)%albedo_type(m,ind_pav_green) = INT( building_pars(ind_alb_green_r,building_type) )
           surf_usm_h(l)%albedo_type(m,ind_wat_win)   = INT( building_pars(ind_alb_win_r,building_type) )

           surf_usm_h(l)%zw(nzb_wall,m)         = building_pars(ind_thick_1_wall_r,building_type)
           surf_usm_h(l)%zw(nzb_wall+1,m)       = building_pars(ind_thick_2_wall_r,building_type)
           surf_usm_h(l)%zw(nzb_wall+2,m)       = building_pars(ind_thick_3_wall_r,building_type)
           surf_usm_h(l)%zw(nzb_wall+3,m)       = building_pars(ind_thick_4_wall_r,building_type)

           surf_usm_h(l)%zw_green(nzb_wall,m)         = building_pars(ind_thick_1_wall_r,building_type)
           surf_usm_h(l)%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2_wall_r,building_type)
           surf_usm_h(l)%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3_wall_r,building_type)
           surf_usm_h(l)%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4_wall_r,building_type)

           surf_usm_h(l)%zw_window(nzb_wall,m)         = building_pars(ind_thick_1_win_r,building_type)
           surf_usm_h(l)%zw_window(nzb_wall+1,m)       = building_pars(ind_thick_2_win_r,building_type)
           surf_usm_h(l)%zw_window(nzb_wall+2,m)       = building_pars(ind_thick_3_win_r,building_type)
           surf_usm_h(l)%zw_window(nzb_wall+3,m)       = building_pars(ind_thick_4_win_r,building_type)

           surf_usm_h(l)%green_type_roof(m)     = building_pars(ind_green_type_roof,building_type)
        ENDDO
     ENDDO

     DO  l = 0, 3
        DO  m = 1, surf_usm_v(l)%ns

           surf_usm_v(l)%surface_types(m) = wall_category     !< Default category for root surface
!
!--        In order to distinguish between ground floor level and above-ground-floor level surfaces,
!--        set input indices.
           ind_alb_green    = MERGE( ind_alb_green_gfl,    ind_alb_green_agfl,                     &
                                     surf_usm_v(l)%ground_level(m) )
           ind_alb_wall     = MERGE( ind_alb_wall_gfl,     ind_alb_wall_agfl,                      &
                                     surf_usm_v(l)%ground_level(m) )
           ind_alb_win      = MERGE( ind_alb_win_gfl,      ind_alb_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,                     &
                                     surf_usm_v(l)%ground_level(m) )
           ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,                      &
                                     surf_usm_v(l)%ground_level(m) )
           ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl,                  &
                                     surf_usm_v(l)%ground_level(m) )
           ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl,                  &
                                     surf_usm_v(l)%ground_level(m) )
           ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,                         &
                                     surf_usm_v(l)%ground_level(m) )
           ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,                         &
                                     surf_usm_v(l)%ground_level(m) )
           ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,                           &
                                     surf_usm_v(l)%ground_level(m) )
           ind_hc1_win      = MERGE( ind_hc1_win_gfl,      ind_hc1_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,                           &
                                     surf_usm_v(l)%ground_level(m) )
           ind_hc2_win      = MERGE( ind_hc2_win_gfl,      ind_hc2_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,                           &
                                     surf_usm_v(l)%ground_level(m) )
           ind_hc3_win      = MERGE( ind_hc3_win_gfl,      ind_hc3_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_hc4          = MERGE( ind_hc4_gfl,          ind_hc4_agfl,                           &
                                     surf_usm_v(l)%ground_level(m) )
           ind_hc4_win      = MERGE( ind_hc4_win_gfl,      ind_hc4_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,                           &
                                     surf_usm_v(l)%ground_level(m) )
           ind_tc1_win      = MERGE( ind_tc1_win_gfl,      ind_tc1_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,                           &
                                     surf_usm_v(l)%ground_level(m) )
           ind_tc2_win      = MERGE( ind_tc2_win_gfl,      ind_tc2_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,                           &
                                     surf_usm_v(l)%ground_level(m) )
           ind_tc3_win      = MERGE( ind_tc3_win_gfl,      ind_tc3_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_tc4          = MERGE( ind_tc4_gfl,          ind_tc4_agfl,                           &
                                     surf_usm_v(l)%ground_level(m) )
           ind_tc4_win      = MERGE( ind_tc4_win_gfl,      ind_tc4_win_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_thick_1      = MERGE( ind_thick_1_gfl,      ind_thick_1_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_thick_1_win  = MERGE( ind_thick_1_win_gfl,  ind_thick_1_win_agfl,                   &
                                     surf_usm_v(l)%ground_level(m) )
           ind_thick_2      = MERGE( ind_thick_2_gfl,      ind_thick_2_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_thick_2_win  = MERGE( ind_thick_2_win_gfl,  ind_thick_2_win_agfl,                   &
                                     surf_usm_v(l)%ground_level(m) )
           ind_thick_3      = MERGE( ind_thick_3_gfl,      ind_thick_3_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_thick_3_win  = MERGE( ind_thick_3_win_gfl,  ind_thick_3_win_agfl,                   &
                                     surf_usm_v(l)%ground_level(m) )
           ind_thick_4      = MERGE( ind_thick_4_gfl,      ind_thick_4_agfl,                       &
                                     surf_usm_v(l)%ground_level(m) )
           ind_thick_4_win  = MERGE( ind_thick_4_win_gfl,  ind_thick_4_win_agfl,                   &
                                     surf_usm_v(l)%ground_level(m) )
           ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,                     &
                                     surf_usm_v(l)%ground_level(m) )
           ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,                    &
                                     surf_usm_v(l)%ground_level(m) )
           ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,                      &
                                     surf_usm_v(l)%ground_level(m) )
           ind_trans        = MERGE( ind_trans_gfl,       ind_trans_agfl,                          &
                                     surf_usm_v(l)%ground_level(m) )
           ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,                            &
                                     surf_usm_v(l)%ground_level(m) )
           ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,                          &
                                     surf_usm_v(l)%ground_level(m) )
!
!--        Store building type and its name on each surface element
           surf_usm_v(l)%building_type(m)      = building_type
           surf_usm_v(l)%building_type_name(m) = building_type_name(building_type)
!
!--        Initialize relatvie wall- (0), green- (1) and window (2) fractions
           surf_usm_v(l)%frac(m,ind_veg_wall)   = building_pars(ind_wall_frac,building_type)
           surf_usm_v(l)%frac(m,ind_pav_green)  = building_pars(ind_green_frac_w,building_type)
           surf_usm_v(l)%frac(m,ind_wat_win)    = building_pars(ind_win_frac,building_type)
           surf_usm_v(l)%lai(m)                 = building_pars(ind_lai_w,building_type)

           surf_usm_v(l)%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1,building_type)
           surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc2,building_type)
           surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc3,building_type)
           surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc4,building_type)

           surf_usm_v(l)%rho_c_green(nzb_wall,m)   = rho_c_soil !building_pars(ind_hc1,building_type)
           surf_usm_v(l)%rho_c_green(nzb_wall+1,m) = rho_c_soil !building_pars(ind_hc1,building_type)
           surf_usm_v(l)%rho_c_green(nzb_wall+2,m) = rho_c_soil !building_pars(ind_hc2,building_type)
           surf_usm_v(l)%rho_c_green(nzb_wall+3,m) = rho_c_soil !building_pars(ind_hc3,building_type)

           surf_usm_v(l)%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1_win,building_type)
           surf_usm_v(l)%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc2_win,building_type)
           surf_usm_v(l)%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc3_win,building_type)
           surf_usm_v(l)%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc4_win,building_type)

           surf_usm_v(l)%lambda_h(nzb_wall,m)   = building_pars(ind_tc1,building_type)
           surf_usm_v(l)%lambda_h(nzb_wall+1,m) = building_pars(ind_tc2,building_type)
           surf_usm_v(l)%lambda_h(nzb_wall+2,m) = building_pars(ind_tc3,building_type)
           surf_usm_v(l)%lambda_h(nzb_wall+3,m) = building_pars(ind_tc4,building_type)

           surf_usm_v(l)%lambda_h_green(nzb_wall,m)   = lambda_h_green_sm !building_pars(ind_tc1,building_type)
           surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) = lambda_h_green_sm !building_pars(ind_tc1,building_type)
           surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) = lambda_h_green_sm !building_pars(ind_tc2,building_type)
           surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) = lambda_h_green_sm !building_pars(ind_tc3,building_type)

           surf_usm_v(l)%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1_win,building_type)
           surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc2_win,building_type)
           surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc3_win,building_type)
           surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc4_win,building_type)

           surf_usm_v(l)%target_temp_summer(m)  = building_pars(ind_indoor_target_temp_summer,building_type)
           surf_usm_v(l)%target_temp_winter(m)  = building_pars(ind_indoor_target_temp_winter,building_type)
!
!--        Emissivity of wall-, green- and window fraction
           surf_usm_v(l)%emissivity(m,ind_veg_wall)  = building_pars(ind_emis_wall,building_type)
           surf_usm_v(l)%emissivity(m,ind_pav_green) = building_pars(ind_emis_green,building_type)
           surf_usm_v(l)%emissivity(m,ind_wat_win)   = building_pars(ind_emis_win,building_type)

           surf_usm_v(l)%transmissivity(m)      = building_pars(ind_trans,building_type)

           surf_usm_v(l)%z0(m)                  = building_pars(ind_z0,building_type)
           surf_usm_v(l)%z0h(m)                 = building_pars(ind_z0qh,building_type)
           surf_usm_v(l)%z0q(m)                 = building_pars(ind_z0qh,building_type)

           surf_usm_v(l)%albedo_type(m,ind_veg_wall)  = INT( building_pars(ind_alb_wall,building_type) )
           surf_usm_v(l)%albedo_type(m,ind_pav_green) = INT( building_pars(ind_alb_green,building_type) )
           surf_usm_v(l)%albedo_type(m,ind_wat_win)   = INT( building_pars(ind_alb_win,building_type) )

           surf_usm_v(l)%zw(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
           surf_usm_v(l)%zw(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
           surf_usm_v(l)%zw(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
           surf_usm_v(l)%zw(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)

           surf_usm_v(l)%zw_green(nzb_wall,m)         = building_pars(ind_thick_1,building_type)
           surf_usm_v(l)%zw_green(nzb_wall+1,m)       = building_pars(ind_thick_2,building_type)
           surf_usm_v(l)%zw_green(nzb_wall+2,m)       = building_pars(ind_thick_3,building_type)
           surf_usm_v(l)%zw_green(nzb_wall+3,m)       = building_pars(ind_thick_4,building_type)

           surf_usm_v(l)%zw_window(nzb_wall,m)        = building_pars(ind_thick_1_win,building_type)
           surf_usm_v(l)%zw_window(nzb_wall+1,m)      = building_pars(ind_thick_2_win,building_type)
           surf_usm_v(l)%zw_window(nzb_wall+2,m)      = building_pars(ind_thick_3_win,building_type)
           surf_usm_v(l)%zw_window(nzb_wall+3,m)      = building_pars(ind_thick_4_win,building_type)

        ENDDO
     ENDDO
!
!--  Level 2 - initialization via building type read from file
     IF ( building_type_f%from_file )  THEN
        DO  l = 0, 1
           DO  m = 1, surf_usm_h(l)%ns
              i = surf_usm_h(l)%i(m)
              j = surf_usm_h(l)%j(m)
!
!--           For the moment, limit building type to 6 (to overcome errors in input file).
              st = building_type_f%var(j,i)
              IF ( st /= building_type_f%fill )  THEN

!
!--              In order to distinguish between ground floor level and above-ground-floor level
!--              surfaces, set input indices.

                 ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl,               &
                                           surf_usm_h(l)%ground_level(m) )
                 ind_lai_r        = MERGE( ind_lai_r_gfl, ind_lai_r_agfl, surf_usm_h(l)%ground_level(m) )
                 ind_z0           = MERGE( ind_z0_gfl, ind_z0_agfl, surf_usm_h(l)%ground_level(m) )
                 ind_z0qh         = MERGE( ind_z0qh_gfl, ind_z0qh_agfl, surf_usm_h(l)%ground_level(m) )
!
!--              Store building type and its name on each surface element
                 surf_usm_h(l)%building_type(m)      = st
                 surf_usm_h(l)%building_type_name(m) = building_type_name(st)
!
!--              Initialize relatvie wall- (0), green- (1) and window (2) fractions
                 surf_usm_h(l)%frac(m,ind_veg_wall)  = building_pars(ind_wall_frac_r,st)
                 surf_usm_h(l)%frac(m,ind_pav_green) = building_pars(ind_green_frac_r,st)
                 surf_usm_h(l)%frac(m,ind_wat_win)   = building_pars(ind_win_frac_r,st)
                 surf_usm_h(l)%lai(m)                = building_pars(ind_lai_r,st)

                 surf_usm_h(l)%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1_wall_r,st)
                 surf_usm_h(l)%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc2_wall_r,st)
                 surf_usm_h(l)%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc3_wall_r,st)
                 surf_usm_h(l)%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc4_wall_r,st)
                 surf_usm_h(l)%lambda_h(nzb_wall,m)     = building_pars(ind_tc1_wall_r,st)
                 surf_usm_h(l)%lambda_h(nzb_wall+1,m)   = building_pars(ind_tc2_wall_r,st)
                 surf_usm_h(l)%lambda_h(nzb_wall+2,m)   = building_pars(ind_tc3_wall_r,st)
                 surf_usm_h(l)%lambda_h(nzb_wall+3,m)   = building_pars(ind_tc4_wall_r,st)

                 surf_usm_h(l)%rho_c_green(nzb_wall,m)      = rho_c_soil !building_pars(ind_hc1_wall_r,st)
                 surf_usm_h(l)%rho_c_green(nzb_wall+1,m)    = rho_c_soil !building_pars(ind_hc1_wall_r,st)
                 surf_usm_h(l)%rho_c_green(nzb_wall+2,m)    = rho_c_soil !building_pars(ind_hc2_wall_r,st)
                 surf_usm_h(l)%rho_c_green(nzb_wall+3,m)    = rho_c_soil !building_pars(ind_hc3_wall_r,st)
                 surf_usm_h(l)%lambda_h_green(nzb_wall,m)   = lambda_h_green_sm !building_pars(ind_tc1_wall_r,st)
                 surf_usm_h(l)%lambda_h_green(nzb_wall+1,m) = lambda_h_green_sm !building_pars(ind_tc1_wall_r,st)
                 surf_usm_h(l)%lambda_h_green(nzb_wall+2,m) = lambda_h_green_sm !building_pars(ind_tc2_wall_r,st)
                 surf_usm_h(l)%lambda_h_green(nzb_wall+3,m) = lambda_h_green_sm !building_pars(ind_tc3_wall_r,st)

                 surf_usm_h(l)%rho_c_window(nzb_wall,m)      = building_pars(ind_hc1_win_r,st)
                 surf_usm_h(l)%rho_c_window(nzb_wall+1,m)    = building_pars(ind_hc2_win_r,st)
                 surf_usm_h(l)%rho_c_window(nzb_wall+2,m)    = building_pars(ind_hc3_win_r,st)
                 surf_usm_h(l)%rho_c_window(nzb_wall+3,m)    = building_pars(ind_hc4_win_r,st)
                 surf_usm_h(l)%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1_win_r,st)
                 surf_usm_h(l)%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc2_win_r,st)
                 surf_usm_h(l)%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc3_win_r,st)
                 surf_usm_h(l)%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc4_win_r,st)

                 surf_usm_h(l)%target_temp_summer(m)  = building_pars(ind_indoor_target_temp_summer,st)
                 surf_usm_h(l)%target_temp_winter(m)  = building_pars(ind_indoor_target_temp_winter,st)
!
!--              Emissivity of wall-, green- and window fraction
                 surf_usm_h(l)%emissivity(m,ind_veg_wall)  = building_pars(ind_emis_wall_r,st)
                 surf_usm_h(l)%emissivity(m,ind_pav_green) = building_pars(ind_emis_green_r,st)
                 surf_usm_h(l)%emissivity(m,ind_wat_win)   = building_pars(ind_emis_win_r,st)

                 surf_usm_h(l)%transmissivity(m)      = building_pars(ind_trans_r,st)

                 surf_usm_h(l)%z0(m)                  = building_pars(ind_z0,st)
                 surf_usm_h(l)%z0h(m)                 = building_pars(ind_z0qh,st)
                 surf_usm_h(l)%z0q(m)                 = building_pars(ind_z0qh,st)
!
!--              Albedo type for wall fraction, green fraction, window fraction
                 surf_usm_h(l)%albedo_type(m,ind_veg_wall)  = INT( building_pars(ind_alb_wall_r,st) )
                 surf_usm_h(l)%albedo_type(m,ind_pav_green) = INT( building_pars(ind_alb_green_r,st) )
                 surf_usm_h(l)%albedo_type(m,ind_wat_win)   = INT( building_pars(ind_alb_win_r,st) )

                 surf_usm_h(l)%zw(nzb_wall,m)   = building_pars(ind_thick_1_wall_r,st)
                 surf_usm_h(l)%zw(nzb_wall+1,m) = building_pars(ind_thick_2_wall_r,st)
                 surf_usm_h(l)%zw(nzb_wall+2,m) = building_pars(ind_thick_3_wall_r,st)
                 surf_usm_h(l)%zw(nzb_wall+3,m) = building_pars(ind_thick_4_wall_r,st)

                 surf_usm_h(l)%zw_green(nzb_wall,m)   = building_pars(ind_thick_1_wall_r,st)
                 surf_usm_h(l)%zw_green(nzb_wall+1,m) = building_pars(ind_thick_2_wall_r,st)
                 surf_usm_h(l)%zw_green(nzb_wall+2,m) = building_pars(ind_thick_3_wall_r,st)
                 surf_usm_h(l)%zw_green(nzb_wall+3,m) = building_pars(ind_thick_4_wall_r,st)

                 surf_usm_h(l)%zw_window(nzb_wall,m)   = building_pars(ind_thick_1_win_r,st)
                 surf_usm_h(l)%zw_window(nzb_wall+1,m) = building_pars(ind_thick_2_win_r,st)
                 surf_usm_h(l)%zw_window(nzb_wall+2,m) = building_pars(ind_thick_3_win_r,st)
                 surf_usm_h(l)%zw_window(nzb_wall+3,m) = building_pars(ind_thick_4_win_r,st)

                 surf_usm_h(l)%green_type_roof(m) = building_pars(ind_green_type_roof,st)

              ENDIF
           ENDDO
        ENDDO

        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
              j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
!
!--           For the moment, limit building type to 6 (to overcome errors in input file).

              st = building_type_f%var(j,i)
              IF ( st /= building_type_f%fill )  THEN

!
!--              In order to distinguish between ground floor level and above-ground-floor level
!--              surfaces, set input indices.
                 ind_alb_green    = MERGE( ind_alb_green_gfl,    ind_alb_green_agfl,               &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_alb_wall     = MERGE( ind_alb_wall_gfl,     ind_alb_wall_agfl,                &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_alb_win      = MERGE( ind_alb_win_gfl,      ind_alb_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_wall_frac    = MERGE( ind_wall_frac_gfl,    ind_wall_frac_agfl,               &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_win_frac     = MERGE( ind_win_frac_gfl,     ind_win_frac_agfl,                &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl,            &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl,            &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_lai_r        = MERGE( ind_lai_r_gfl,        ind_lai_r_agfl,                   &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_lai_w        = MERGE( ind_lai_w_gfl,        ind_lai_w_agfl,                   &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc1          = MERGE( ind_hc1_gfl,          ind_hc1_agfl,                     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc1_win      = MERGE( ind_hc1_win_gfl,      ind_hc1_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc2          = MERGE( ind_hc2_gfl,          ind_hc2_agfl,                     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc2_win      = MERGE( ind_hc2_win_gfl,      ind_hc2_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc3          = MERGE( ind_hc3_gfl,          ind_hc3_agfl,                     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc3_win      = MERGE( ind_hc3_win_gfl,      ind_hc3_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc4          = MERGE( ind_hc4_gfl,          ind_hc4_agfl,                     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_hc4_win      = MERGE( ind_hc4_win_gfl,      ind_hc4_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc1          = MERGE( ind_tc1_gfl,          ind_tc1_agfl,                     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc1_win      = MERGE( ind_tc1_win_gfl,      ind_tc1_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc2          = MERGE( ind_tc2_gfl,          ind_tc2_agfl,                     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc2_win      = MERGE( ind_tc2_win_gfl,      ind_tc2_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc3          = MERGE( ind_tc3_gfl,          ind_tc3_agfl,                     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc3_win      = MERGE( ind_tc3_win_gfl,      ind_tc3_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc4          = MERGE( ind_tc4_gfl,          ind_tc4_agfl,                     &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_tc4_win      = MERGE( ind_tc4_win_gfl,      ind_tc4_win_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_thick_1      = MERGE( ind_thick_1_gfl,      ind_thick_1_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_thick_1_win  = MERGE( ind_thick_1_win_gfl,  ind_thick_1_win_agfl,             &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_thick_2      = MERGE( ind_thick_2_gfl,      ind_thick_2_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_thick_2_win  = MERGE( ind_thick_2_win_gfl,  ind_thick_2_win_agfl,             &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_thick_3      = MERGE( ind_thick_3_gfl,      ind_thick_3_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_thick_3_win  = MERGE( ind_thick_3_win_gfl,  ind_thick_3_win_agfl,             &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_thick_4      = MERGE( ind_thick_4_gfl,      ind_thick_4_agfl,                 &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_thick_4_win  = MERGE( ind_thick_4_win_gfl,  ind_thick_4_win_agfl,             &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_emis_wall    = MERGE( ind_emis_wall_gfl,    ind_emis_wall_agfl,               &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_emis_green   = MERGE( ind_emis_green_gfl,   ind_emis_green_agfl,              &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_emis_win     = MERGE( ind_emis_win_gfl,     ind_emis_win_agfl,                &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_trans        = MERGE( ind_trans_gfl,       ind_trans_agfl,                    &
                                         surf_usm_v(l)%ground_level(m) )
                 ind_z0           = MERGE( ind_z0_gfl,           ind_z0_agfl,                      &
                                           surf_usm_v(l)%ground_level(m) )
                 ind_z0qh         = MERGE( ind_z0qh_gfl,         ind_z0qh_agfl,                    &
                                           surf_usm_v(l)%ground_level(m) )
!
!--              Store building type and its name on each surface element
                 surf_usm_v(l)%building_type(m)      = st
                 surf_usm_v(l)%building_type_name(m) = building_type_name(st)
!
!--              Initialize relatvie wall- (0), green- (1) and window (2) fractions
                 surf_usm_v(l)%frac(m,ind_veg_wall)  = building_pars(ind_wall_frac,st)
                 surf_usm_v(l)%frac(m,ind_pav_green) = building_pars(ind_green_frac_w,st)
                 surf_usm_v(l)%frac(m,ind_wat_win)   = building_pars(ind_win_frac,st)
                 surf_usm_v(l)%lai(m)                = building_pars(ind_lai_w,st)

                 surf_usm_v(l)%rho_c_wall(nzb_wall,m)   = building_pars(ind_hc1,st)
                 surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) = building_pars(ind_hc2,st)
                 surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) = building_pars(ind_hc3,st)
                 surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) = building_pars(ind_hc4,st)

                 surf_usm_v(l)%rho_c_green(nzb_wall,m)   = rho_c_soil !building_pars(ind_hc1,st)
                 surf_usm_v(l)%rho_c_green(nzb_wall+1,m) = rho_c_soil !building_pars(ind_hc1,st)
                 surf_usm_v(l)%rho_c_green(nzb_wall+2,m) = rho_c_soil !building_pars(ind_hc2,st)
                 surf_usm_v(l)%rho_c_green(nzb_wall+3,m) = rho_c_soil !building_pars(ind_hc3,st)

                 surf_usm_v(l)%rho_c_window(nzb_wall,m)   = building_pars(ind_hc1_win,st)
                 surf_usm_v(l)%rho_c_window(nzb_wall+1,m) = building_pars(ind_hc2_win,st)
                 surf_usm_v(l)%rho_c_window(nzb_wall+2,m) = building_pars(ind_hc3_win,st)
                 surf_usm_v(l)%rho_c_window(nzb_wall+3,m) = building_pars(ind_hc4_win,st)

                 surf_usm_v(l)%lambda_h(nzb_wall,m)   = building_pars(ind_tc1,st)
                 surf_usm_v(l)%lambda_h(nzb_wall+1,m) = building_pars(ind_tc2,st)
                 surf_usm_v(l)%lambda_h(nzb_wall+2,m) = building_pars(ind_tc3,st)
                 surf_usm_v(l)%lambda_h(nzb_wall+3,m) = building_pars(ind_tc4,st)

                 surf_usm_v(l)%lambda_h_green(nzb_wall,m)   = lambda_h_green_sm !building_pars(ind_tc1,st)
                 surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) = lambda_h_green_sm !building_pars(ind_tc1,st)
                 surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) = lambda_h_green_sm !building_pars(ind_tc2,st)
                 surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) = lambda_h_green_sm !building_pars(ind_tc3,st)

                 surf_usm_v(l)%lambda_h_window(nzb_wall,m)   = building_pars(ind_tc1_win,st)
                 surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) = building_pars(ind_tc2_win,st)
                 surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) = building_pars(ind_tc3_win,st)
                 surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) = building_pars(ind_tc4_win,st)

                 surf_usm_v(l)%target_temp_summer(m) = building_pars(ind_indoor_target_temp_summer,st)
                 surf_usm_v(l)%target_temp_winter(m) = building_pars(ind_indoor_target_temp_winter,st)
!
!--              Emissivity of wall-, green- and window fraction
                 surf_usm_v(l)%emissivity(m,ind_veg_wall)  = building_pars(ind_emis_wall,st)
                 surf_usm_v(l)%emissivity(m,ind_pav_green) = building_pars(ind_emis_green,st)
                 surf_usm_v(l)%emissivity(m,ind_wat_win)   = building_pars(ind_emis_win,st)

                 surf_usm_v(l)%transmissivity(m) = building_pars(ind_trans,st)

                 surf_usm_v(l)%z0(m)  = building_pars(ind_z0,st)
                 surf_usm_v(l)%z0h(m) = building_pars(ind_z0qh,st)
                 surf_usm_v(l)%z0q(m) = building_pars(ind_z0qh,st)

                 surf_usm_v(l)%albedo_type(m,ind_veg_wall)  = INT( building_pars(ind_alb_wall,st) )
                 surf_usm_v(l)%albedo_type(m,ind_pav_green) = INT( building_pars(ind_alb_green,st) )
                 surf_usm_v(l)%albedo_type(m,ind_wat_win)   = INT( building_pars(ind_alb_win,st) )

                 surf_usm_v(l)%zw(nzb_wall,m)   = building_pars(ind_thick_1,st)
                 surf_usm_v(l)%zw(nzb_wall+1,m) = building_pars(ind_thick_2,st)
                 surf_usm_v(l)%zw(nzb_wall+2,m) = building_pars(ind_thick_3,st)
                 surf_usm_v(l)%zw(nzb_wall+3,m) = building_pars(ind_thick_4,st)

                 surf_usm_v(l)%zw_green(nzb_wall,m)   = building_pars(ind_thick_1,st)
                 surf_usm_v(l)%zw_green(nzb_wall+1,m) = building_pars(ind_thick_2,st)
                 surf_usm_v(l)%zw_green(nzb_wall+2,m) = building_pars(ind_thick_3,st)
                 surf_usm_v(l)%zw_green(nzb_wall+3,m) = building_pars(ind_thick_4,st)

                 surf_usm_v(l)%zw_window(nzb_wall,m)   = building_pars(ind_thick_1_win,st)
                 surf_usm_v(l)%zw_window(nzb_wall+1,m) = building_pars(ind_thick_2_win,st)
                 surf_usm_v(l)%zw_window(nzb_wall+2,m) = building_pars(ind_thick_3_win,st)
                 surf_usm_v(l)%zw_window(nzb_wall+3,m) = building_pars(ind_thick_4_win,st)

              ENDIF
           ENDDO
        ENDDO
     ENDIF

!
!--  Level 3 - initialization via building_pars read from file. Note, only variables that are also
!--  defined in the input-standard can be initialized via file. Other variables will be initialized
!--  on level 1 or 2.
     IF ( building_pars_f%from_file )  THEN
        DO  l = 0, 1
           DO  m = 1, surf_usm_h(l)%ns
              i = surf_usm_h(l)%i(m)
              j = surf_usm_h(l)%j(m)

!
!--           In order to distinguish between ground floor level and above-ground-floor level surfaces,
!--           set input indices.
              ind_wall_frac    = MERGE( ind_wall_frac_gfl, ind_wall_frac_agfl,                        &
                                        surf_usm_h(l)%ground_level(m) )
              ind_green_frac_r = MERGE( ind_green_frac_r_gfl, ind_green_frac_r_agfl,                  &
                                        surf_usm_h(l)%ground_level(m) )
              ind_win_frac     = MERGE( ind_win_frac_gfl, ind_win_frac_agfl,                          &
                                        surf_usm_h(l)%ground_level(m) )
              ind_lai_r        = MERGE( ind_lai_r_gfl, ind_lai_r_agfl, surf_usm_h(l)%ground_level(m) )
              ind_z0           = MERGE( ind_z0_gfl, ind_z0_agfl, surf_usm_h(l)%ground_level(m) )
              ind_z0qh         = MERGE( ind_z0qh_gfl, ind_z0qh_agfl, surf_usm_h(l)%ground_level(m) )
              ind_hc1          = MERGE( ind_hc1_gfl, ind_hc1_agfl, surf_usm_h(l)%ground_level(m) )
              ind_hc2          = MERGE( ind_hc2_gfl, ind_hc2_agfl, surf_usm_h(l)%ground_level(m) )
              ind_hc3          = MERGE( ind_hc3_gfl, ind_hc3_agfl, surf_usm_h(l)%ground_level(m) )
              ind_hc4          = MERGE( ind_hc4_gfl, ind_hc4_agfl, surf_usm_h(l)%ground_level(m) )
              ind_tc1          = MERGE( ind_tc1_gfl, ind_tc1_agfl, surf_usm_h(l)%ground_level(m) )
              ind_tc2          = MERGE( ind_tc2_gfl, ind_tc2_agfl, surf_usm_h(l)%ground_level(m) )
              ind_tc3          = MERGE( ind_tc3_gfl, ind_tc3_agfl, surf_usm_h(l)%ground_level(m) )
              ind_tc4          = MERGE( ind_tc4_gfl, ind_tc4_agfl, surf_usm_h(l)%ground_level(m) )
              ind_emis_wall    = MERGE( ind_emis_wall_gfl, ind_emis_wall_agfl,                        &
                                        surf_usm_h(l)%ground_level(m) )
              ind_emis_green   = MERGE( ind_emis_green_gfl, ind_emis_green_agfl,                      &
                                        surf_usm_h(l)%ground_level(m) )
              ind_emis_win     = MERGE( ind_emis_win_gfl, ind_emis_win_agfl,                          &
                                        surf_usm_h(l)%ground_level(m) )
              ind_trans        = MERGE( ind_trans_gfl, ind_trans_agfl, surf_usm_h(l)%ground_level(m) )

!
!--           Initialize relatvie wall- (0), green- (1) and window (2) fractions
              IF ( building_pars_f%pars_xy(ind_wall_frac,j,i) /= building_pars_f%fill )               &
                 surf_usm_h(l)%frac(m,ind_veg_wall) = building_pars_f%pars_xy(ind_wall_frac,j,i)

              IF ( building_pars_f%pars_xy(ind_green_frac_r,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%frac(m,ind_pav_green) = building_pars_f%pars_xy(ind_green_frac_r,j,i)

              IF ( building_pars_f%pars_xy(ind_win_frac,j,i) /= building_pars_f%fill )                &
                 surf_usm_h(l)%frac(m,ind_wat_win) = building_pars_f%pars_xy(ind_win_frac,j,i)

              IF ( building_pars_f%pars_xy(ind_lai_r,j,i) /= building_pars_f%fill )                   &
                 surf_usm_h(l)%lai(m) = building_pars_f%pars_xy(ind_lai_r,j,i)

              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_wall(nzb_wall,m) = building_pars_f%pars_xy(ind_hc1,j,i)

              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_wall(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc2,j,i)

              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_wall(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc3,j,i)

              IF ( building_pars_f%pars_xy(ind_hc4,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_wall(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc4,j,i)

              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_green(nzb_wall,m) = building_pars_f%pars_xy(ind_hc1,j,i)

              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc2,j,i)

              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc3,j,i)

              IF ( building_pars_f%pars_xy(ind_hc4,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc4,j,i)

              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_window(nzb_wall,m) = building_pars_f%pars_xy(ind_hc1,j,i)

              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc2,j,i)

              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc3,j,i)

              IF ( building_pars_f%pars_xy(ind_hc4,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%rho_c_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc4,j,i)

              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h(nzb_wall,m)   = building_pars_f%pars_xy(ind_tc1,j,i)

              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc2,j,i)

              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc3,j,i)

              IF ( building_pars_f%pars_xy(ind_tc4,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc4,j,i)

              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h_green(nzb_wall,m) = building_pars_f%pars_xy(ind_tc1,j,i)

              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc2,j,i)

              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc3,j,i)

              IF ( building_pars_f%pars_xy(ind_tc4,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc4,j,i)

              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h_window(nzb_wall,m) = building_pars_f%pars_xy(ind_tc1,j,i)

              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc2,j,i)

              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc3,j,i)

              IF ( building_pars_f%pars_xy(ind_tc4,j,i) /= building_pars_f%fill )                     &
                 surf_usm_h(l)%lambda_h_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc4,j,i)

              IF ( building_pars_f%pars_xy(ind_indoor_target_temp_summer,j,i) /=                      &
                   building_pars_f%fill )                                                             &
                 surf_usm_h(l)%target_temp_summer(m) = building_pars_f%pars_xy(ind_indoor_target_temp_summer,j,i)

              IF ( building_pars_f%pars_xy(ind_indoor_target_temp_winter,j,i) /=                      &
                   building_pars_f%fill )                                                             &
                 surf_usm_h(l)%target_temp_winter(m) = building_pars_f%pars_xy(ind_indoor_target_temp_winter,j,i)

              IF ( building_pars_f%pars_xy(ind_emis_wall,j,i) /= building_pars_f%fill )               &
                 surf_usm_h(l)%emissivity(m,ind_veg_wall) = building_pars_f%pars_xy(ind_emis_wall,j,i)

              IF ( building_pars_f%pars_xy(ind_emis_green,j,i) /= building_pars_f%fill )              &
                 surf_usm_h(l)%emissivity(m,ind_pav_green) = building_pars_f%pars_xy(ind_emis_green,j,i)

              IF ( building_pars_f%pars_xy(ind_emis_win,j,i) /= building_pars_f%fill )                &
                 surf_usm_h(l)%emissivity(m,ind_wat_win) = building_pars_f%pars_xy(ind_emis_win,j,i)

              IF ( building_pars_f%pars_xy(ind_trans,j,i) /= building_pars_f%fill )                   &
                 surf_usm_h(l)%transmissivity(m) = building_pars_f%pars_xy(ind_trans,j,i)

              IF ( building_pars_f%pars_xy(ind_z0,j,i) /= building_pars_f%fill )                      &
                 surf_usm_h(l)%z0(m) = building_pars_f%pars_xy(ind_z0,j,i)

              IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /= building_pars_f%fill )                    &
                 surf_usm_h(l)%z0h(m) = building_pars_f%pars_xy(ind_z0qh,j,i)

              IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /= building_pars_f%fill )                    &
                 surf_usm_h(l)%z0q(m) = building_pars_f%pars_xy(ind_z0qh,j,i)

              IF ( building_pars_f%pars_xy(ind_alb_wall_agfl,j,i) /= building_pars_f%fill )           &
                 surf_usm_h(l)%albedo_type(m,ind_veg_wall)  = building_pars_f%pars_xy(ind_alb_wall_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_alb_green_agfl,j,i) /= building_pars_f%fill )          &
                 surf_usm_h(l)%albedo_type(m,ind_pav_green) = building_pars_f%pars_xy(ind_alb_green_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_alb_win_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%albedo_type(m,ind_wat_win) = building_pars_f%pars_xy(ind_alb_win_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_1_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%zw(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_2_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%zw(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_3_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%zw(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_4_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%zw(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_1_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%zw_green(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_2_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%zw_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_3_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%zw_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_4_agfl,j,i) /= building_pars_f%fill )            &
                 surf_usm_h(l)%zw_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4_agfl,j,i)
           ENDDO
        ENDDO

        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
              j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff

!
!--           In order to distinguish between ground floor level and above-ground-floor level
!--           surfaces, set input indices.
              ind_wall_frac    = MERGE( ind_wall_frac_gfl, ind_wall_frac_agfl, surf_usm_v(l)%ground_level(m) )
              ind_green_frac_w = MERGE( ind_green_frac_w_gfl, ind_green_frac_w_agfl, surf_usm_v(l)%ground_level(m) )
              ind_win_frac     = MERGE( ind_win_frac_gfl, ind_win_frac_agfl, surf_usm_v(l)%ground_level(m) )
              ind_lai_w        = MERGE( ind_lai_w_gfl, ind_lai_w_agfl, surf_usm_v(l)%ground_level(m) )
              ind_z0           = MERGE( ind_z0_gfl, ind_z0_agfl, surf_usm_v(l)%ground_level(m) )
              ind_z0qh         = MERGE( ind_z0qh_gfl, ind_z0qh_agfl, surf_usm_v(l)%ground_level(m) )
              ind_hc1          = MERGE( ind_hc1_gfl, ind_hc1_agfl, surf_usm_v(l)%ground_level(m) )
              ind_hc2          = MERGE( ind_hc2_gfl, ind_hc2_agfl, surf_usm_v(l)%ground_level(m) )
              ind_hc3          = MERGE( ind_hc3_gfl, ind_hc3_agfl, surf_usm_v(l)%ground_level(m) )
              ind_hc4          = MERGE( ind_hc4_gfl, ind_hc4_agfl, surf_usm_v(l)%ground_level(m) )
              ind_tc1          = MERGE( ind_tc1_gfl, ind_tc1_agfl, surf_usm_v(l)%ground_level(m) )
              ind_tc2          = MERGE( ind_tc2_gfl, ind_tc2_agfl, surf_usm_v(l)%ground_level(m) )
              ind_tc3          = MERGE( ind_tc3_gfl, ind_tc3_agfl, surf_usm_v(l)%ground_level(m) )
              ind_tc4          = MERGE( ind_tc4_gfl, ind_tc4_agfl, surf_usm_v(l)%ground_level(m) )
              ind_emis_wall    = MERGE( ind_emis_wall_gfl, ind_emis_wall_agfl, surf_usm_v(l)%ground_level(m) )
              ind_emis_green   = MERGE( ind_emis_green_gfl, ind_emis_green_agfl, surf_usm_v(l)%ground_level(m) )
              ind_emis_win     = MERGE( ind_emis_win_gfl, ind_emis_win_agfl, surf_usm_v(l)%ground_level(m) )
              ind_trans        = MERGE( ind_trans_gfl, ind_trans_agfl, surf_usm_v(l)%ground_level(m) )
!
!--           Initialize relatvie wall- (0), green- (1) and window (2) fractions
              IF ( building_pars_f%pars_xy(ind_wall_frac,j,i) /= building_pars_f%fill )            &
                 surf_usm_v(l)%frac(m,ind_veg_wall) = building_pars_f%pars_xy(ind_wall_frac,j,i)

              IF ( building_pars_f%pars_xy(ind_green_frac_w,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%frac(m,ind_pav_green) = building_pars_f%pars_xy(ind_green_frac_w,j,i)

              IF ( building_pars_f%pars_xy(ind_win_frac,j,i) /= building_pars_f%fill )             &
                 surf_usm_v(l)%frac(m,ind_wat_win) = building_pars_f%pars_xy(ind_win_frac,j,i)

              IF ( building_pars_f%pars_xy(ind_lai_w,j,i) /= building_pars_f%fill )                &
                 surf_usm_v(l)%lai(m) = building_pars_f%pars_xy(ind_lai_w,j,i)

              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_wall(nzb_wall,m) = building_pars_f%pars_xy(ind_hc1,j,i)

              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_wall(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc2,j,i)

              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc3,j,i)

              IF ( building_pars_f%pars_xy(ind_hc4,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc4,j,i)

              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_green(nzb_wall,m) = building_pars_f%pars_xy(ind_hc1,j,i)

              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc2,j,i)

              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc3,j,i)

              IF ( building_pars_f%pars_xy(ind_hc4,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc4,j,i)

              IF ( building_pars_f%pars_xy(ind_hc1,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_window(nzb_wall,m) = building_pars_f%pars_xy(ind_hc1,j,i)

              IF ( building_pars_f%pars_xy(ind_hc2,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_hc2,j,i)

              IF ( building_pars_f%pars_xy(ind_hc3,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_hc3,j,i)

              IF ( building_pars_f%pars_xy(ind_hc4,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%rho_c_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_hc4,j,i)

              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h(nzb_wall,m) = building_pars_f%pars_xy(ind_tc1,j,i)

              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc2,j,i)

              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc3,j,i)

              IF ( building_pars_f%pars_xy(ind_tc4,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc4,j,i)

              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h_green(nzb_wall,m) = building_pars_f%pars_xy(ind_tc1,j,i)

              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h_green(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc2,j,i)

              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc3,j,i)

              IF ( building_pars_f%pars_xy(ind_tc4,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc4,j,i)

              IF ( building_pars_f%pars_xy(ind_tc1,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h_window(nzb_wall,m) = building_pars_f%pars_xy(ind_tc1,j,i)

              IF ( building_pars_f%pars_xy(ind_tc2,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h_window(nzb_wall+1,m) = building_pars_f%pars_xy(ind_tc2,j,i)

              IF ( building_pars_f%pars_xy(ind_tc3,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) = building_pars_f%pars_xy(ind_tc3,j,i)

              IF ( building_pars_f%pars_xy(ind_tc4,j,i) /= building_pars_f%fill )                  &
                 surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) = building_pars_f%pars_xy(ind_tc4,j,i)

              IF ( building_pars_f%pars_xy(ind_indoor_target_temp_summer,j,i) /=                   &
                   building_pars_f%fill )                                                          &
                 surf_usm_v(l)%target_temp_summer(m) =                                             &
                 building_pars_f%pars_xy(ind_indoor_target_temp_summer,j,i)

              IF ( building_pars_f%pars_xy(ind_indoor_target_temp_winter,j,i) /=                   &
                   building_pars_f%fill )                                                          &
                 surf_usm_v(l)%target_temp_winter(m) =                                             &
                 building_pars_f%pars_xy(ind_indoor_target_temp_winter,j,i)

              IF ( building_pars_f%pars_xy(ind_emis_wall,j,i) /= building_pars_f%fill )            &
                 surf_usm_v(l)%emissivity(m,ind_veg_wall) =                                        &
                 building_pars_f%pars_xy(ind_emis_wall,j,i)

              IF ( building_pars_f%pars_xy(ind_emis_green,j,i) /= building_pars_f%fill )           &
                 surf_usm_v(l)%emissivity(m,ind_pav_green) =                                       &
                 building_pars_f%pars_xy(ind_emis_green,j,i)

              IF ( building_pars_f%pars_xy(ind_emis_win,j,i) /= building_pars_f%fill )             &
                 surf_usm_v(l)%emissivity(m,ind_wat_win)   =                                       &
                 building_pars_f%pars_xy(ind_emis_win,j,i)

              IF ( building_pars_f%pars_xy(ind_trans,j,i) /= building_pars_f%fill )                &
                 surf_usm_v(l)%transmissivity(m) =                                                 &
                 building_pars_f%pars_xy(ind_trans,j,i)

              IF ( building_pars_f%pars_xy(ind_z0,j,i) /= building_pars_f%fill )                   &
                 surf_usm_v(l)%z0(m) = building_pars_f%pars_xy(ind_z0,j,i)

              IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /= building_pars_f%fill )                 &
                 surf_usm_v(l)%z0h(m) = building_pars_f%pars_xy(ind_z0qh,j,i)

              IF ( building_pars_f%pars_xy(ind_z0qh,j,i) /= building_pars_f%fill )                 &
                 surf_usm_v(l)%z0q(m) = building_pars_f%pars_xy(ind_z0qh,j,i)

              IF ( building_pars_f%pars_xy(ind_alb_wall_agfl,j,i) /= building_pars_f%fill )        &
                 surf_usm_v(l)%albedo_type(m,ind_veg_wall)  =                                      &
                 building_pars_f%pars_xy(ind_alb_wall_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_alb_green_agfl,j,i) /= building_pars_f%fill )       &
                 surf_usm_v(l)%albedo_type(m,ind_pav_green) =                                      &
                 building_pars_f%pars_xy(ind_alb_green_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_alb_win_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%albedo_type(m,ind_wat_win)   =                                      &
                 building_pars_f%pars_xy(ind_alb_win_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_1_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%zw(nzb_wall,m) = building_pars_f%pars_xy(ind_thick_1_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_2_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%zw(nzb_wall+1,m) = building_pars_f%pars_xy(ind_thick_2_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_3_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%zw(nzb_wall+2,m) = building_pars_f%pars_xy(ind_thick_3_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_4_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%zw(nzb_wall+3,m) = building_pars_f%pars_xy(ind_thick_4_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_1_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%zw_green(nzb_wall,m) =                                              &
                 building_pars_f%pars_xy(ind_thick_1_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_2_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%zw_green(nzb_wall+1,m) =                                            &
                 building_pars_f%pars_xy(ind_thick_2_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_3_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%zw_green(nzb_wall+2,m) =                                            &
                 building_pars_f%pars_xy(ind_thick_3_agfl,j,i)

              IF ( building_pars_f%pars_xy(ind_thick_4_agfl,j,i) /= building_pars_f%fill )         &
                 surf_usm_v(l)%zw_green(nzb_wall+3,m) =                                            &
                 building_pars_f%pars_xy(ind_thick_4_agfl,j,i)

           ENDDO
        ENDDO
     ENDIF
!
!--  Read building surface pars. If present, they override LOD1-LOD3 building pars where applicable
     IF ( building_surface_pars_f%from_file )  THEN
        DO  l = 0, 1
           DO  m = 1, surf_usm_h(l)%ns
              i = surf_usm_h(l)%i(m)
              j = surf_usm_h(l)%j(m)
              k = surf_usm_h(l)%k(m)
!
!--           Iterate over surfaces in column, check height and orientation
              DO  is = building_surface_pars_f%index_ji(1,j,i), &
                       building_surface_pars_f%index_ji(2,j,i)
                 IF ( building_surface_pars_f%coords(4,is) == -surf_usm_h(l)%koff .AND.               &
                      building_surface_pars_f%coords(1,is) == k )  THEN

                    IF ( building_surface_pars_f%pars(ind_s_wall_frac,is) /=                          &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%frac(m,ind_veg_wall) =                                           &
                       building_surface_pars_f%pars(ind_s_wall_frac,is)

                    IF ( building_surface_pars_f%pars(ind_s_green_frac_w,is) /=                       &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%frac(m,ind_pav_green) =                                          &
                       building_surface_pars_f%pars(ind_s_green_frac_w,is)

                    IF ( building_surface_pars_f%pars(ind_s_green_frac_r,is) /=                       &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%frac(m,ind_pav_green) =                                          &
                       building_surface_pars_f%pars(ind_s_green_frac_r,is)
                       !TODO clarify: why should _w and _r be on the same surface?

                    IF ( building_surface_pars_f%pars(ind_s_win_frac,is) /=                           &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%frac(m,ind_wat_win) = building_surface_pars_f%pars(ind_s_win_frac,is)

                    IF ( building_surface_pars_f%pars(ind_s_lai_r,is) /=                              &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%lai(m) = building_surface_pars_f%pars(ind_s_lai_r,is)

                    IF ( building_surface_pars_f%pars(ind_s_hc1,is) /=                                &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h(l)%rho_c_wall(nzb_wall:nzb_wall+1,m) =                              &
                       building_surface_pars_f%pars(ind_s_hc1,is)
                       surf_usm_h(l)%rho_c_green(nzb_wall:nzb_wall+1,m) =                             &
                       building_surface_pars_f%pars(ind_s_hc1,is)
                       surf_usm_h(l)%rho_c_window(nzb_wall:nzb_wall+1,m) =                            &
                       building_surface_pars_f%pars(ind_s_hc1,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_hc2,is) /=                                &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h(l)%rho_c_wall(nzb_wall+2,m) =                                       &
                       building_surface_pars_f%pars(ind_s_hc2,is)
                       surf_usm_h(l)%rho_c_green(nzb_wall+2,m) =                                      &
                       building_surface_pars_f%pars(ind_s_hc2,is)
                       surf_usm_h(l)%rho_c_window(nzb_wall+2,m) =                                     &
                       building_surface_pars_f%pars(ind_s_hc2,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_hc3,is) /=                                &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h(l)%rho_c_wall(nzb_wall+3,m) =                                       &
                       building_surface_pars_f%pars(ind_s_hc3,is)
                       surf_usm_h(l)%rho_c_green(nzb_wall+3,m) =                                      &
                       building_surface_pars_f%pars(ind_s_hc3,is)
                       surf_usm_h(l)%rho_c_window(nzb_wall+3,m) =                                     &
                       building_surface_pars_f%pars(ind_s_hc3,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc1,is) /=                                &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h(l)%lambda_h(nzb_wall:nzb_wall+1,m) =                                &
                       building_surface_pars_f%pars(ind_s_tc1,is)
                       surf_usm_h(l)%lambda_h_green(nzb_wall:nzb_wall+1,m) =                          &
                       building_surface_pars_f%pars(ind_s_tc1,is)
                       surf_usm_h(l)%lambda_h_window(nzb_wall:nzb_wall+1,m) =                         &
                       building_surface_pars_f%pars(ind_s_tc1,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc2,is) /=                                &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h(l)%lambda_h(nzb_wall+2,m) =                                         &
                       building_surface_pars_f%pars(ind_s_tc2,is)
                       surf_usm_h(l)%lambda_h_green(nzb_wall+2,m) =                                   &
                       building_surface_pars_f%pars(ind_s_tc2,is)
                       surf_usm_h(l)%lambda_h_window(nzb_wall+2,m) =                                  &
                       building_surface_pars_f%pars(ind_s_tc2,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc3,is) /=                                &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h(l)%lambda_h(nzb_wall+3,m) =                                         &
                       building_surface_pars_f%pars(ind_s_tc3,is)
                       surf_usm_h(l)%lambda_h_green(nzb_wall+3,m) =                                   &
                       building_surface_pars_f%pars(ind_s_tc3,is)
                       surf_usm_h(l)%lambda_h_window(nzb_wall+3,m) =                                  &
                       building_surface_pars_f%pars(ind_s_tc3,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_indoor_target_temp_summer,is) /=          &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%target_temp_summer(m) =                                          &
                       building_surface_pars_f%pars(ind_s_indoor_target_temp_summer,is)

                    IF ( building_surface_pars_f%pars(ind_s_indoor_target_temp_winter,is) /=          &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%target_temp_winter(m) =                                          &
                       building_surface_pars_f%pars(ind_s_indoor_target_temp_winter,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_wall,is) /=                          &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%emissivity(m,ind_veg_wall) =                                     &
                       building_surface_pars_f%pars(ind_s_emis_wall,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_green,is) /=                         &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%emissivity(m,ind_pav_green) =                                    &
                       building_surface_pars_f%pars(ind_s_emis_green,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_win,is) /=                           &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%emissivity(m,ind_wat_win) =                                      &
                       building_surface_pars_f%pars(ind_s_emis_win,is)

                    IF ( building_surface_pars_f%pars(ind_s_trans,is) /=                              &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%transmissivity(m) = building_surface_pars_f%pars(ind_s_trans,is)

                    IF ( building_surface_pars_f%pars(ind_s_z0,is) /=                                 &
                         building_surface_pars_f%fill )                                               &
                       surf_usm_h(l)%z0(m) = building_surface_pars_f%pars(ind_s_z0,is)

                    IF ( building_surface_pars_f%pars(ind_s_z0qh,is) /=                               &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_h(l)%z0q(m) = building_surface_pars_f%pars(ind_s_z0qh,is)
                       surf_usm_h(l)%z0h(m) = building_surface_pars_f%pars(ind_s_z0qh,is)
                    ENDIF

                    EXIT ! Surface was found and processed
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              i = surf_usm_v(l)%i(m)
              j = surf_usm_v(l)%j(m)
              k = surf_usm_v(l)%k(m)
!
!--           Iterate over surfaces in column, check height and orientation
              DO  is = building_surface_pars_f%index_ji(1,j,i),                                    &
                       building_surface_pars_f%index_ji(2,j,i)
                 IF ( building_surface_pars_f%coords(5,is) == -surf_usm_v(l)%joff .AND.            &
                      building_surface_pars_f%coords(6,is) == -surf_usm_v(l)%ioff .AND.            &
                      building_surface_pars_f%coords(1,is) == k )  THEN

                    IF ( building_surface_pars_f%pars(ind_s_wall_frac,is) /=                       &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%frac(m,ind_veg_wall) =                                        &
                       building_surface_pars_f%pars(ind_s_wall_frac,is)

                    IF ( building_surface_pars_f%pars(ind_s_green_frac_w,is) /=                    &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%frac(m,ind_pav_green) =                                       &
                       building_surface_pars_f%pars(ind_s_green_frac_w,is)

                    IF ( building_surface_pars_f%pars(ind_s_green_frac_r,is) /=                    &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%frac(m,ind_pav_green) =                                       &
                       building_surface_pars_f%pars(ind_s_green_frac_r,is)
                       !TODO Clarify: why should _w and _r be on the same surface?

                    IF ( building_surface_pars_f%pars(ind_s_win_frac,is) /=                        &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%frac(m,ind_wat_win) =                                         &
                       building_surface_pars_f%pars(ind_s_win_frac,is)

                    IF ( building_surface_pars_f%pars(ind_s_lai_r,is) /=                           &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%lai(m) = building_surface_pars_f%pars(ind_s_lai_r,is)

                    IF ( building_surface_pars_f%pars(ind_s_hc1,is) /=                             &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_v(l)%rho_c_wall(nzb_wall:nzb_wall+1,m) =                           &
                       building_surface_pars_f%pars(ind_s_hc1,is)
                       surf_usm_v(l)%rho_c_green(nzb_wall:nzb_wall+1,m) =                          &
                       building_surface_pars_f%pars(ind_s_hc1,is)
                       surf_usm_v(l)%rho_c_window(nzb_wall:nzb_wall+1,m) =                         &
                       building_surface_pars_f%pars(ind_s_hc1,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_hc2,is) /=                             &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_v(l)%rho_c_wall(nzb_wall+2,m) =                                    &
                       building_surface_pars_f%pars(ind_s_hc2,is)
                       surf_usm_v(l)%rho_c_green(nzb_wall+2,m) =                                   &
                       building_surface_pars_f%pars(ind_s_hc2,is)
                       surf_usm_v(l)%rho_c_window(nzb_wall+2,m) =                                  &
                       building_surface_pars_f%pars(ind_s_hc2,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_hc3,is) /=                             &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_v(l)%rho_c_wall(nzb_wall+3,m) =                                    &
                       building_surface_pars_f%pars(ind_s_hc3,is)
                       surf_usm_v(l)%rho_c_green(nzb_wall+3,m) =                                   &
                       building_surface_pars_f%pars(ind_s_hc3,is)
                       surf_usm_v(l)%rho_c_window(nzb_wall+3,m) =                                  &
                       building_surface_pars_f%pars(ind_s_hc3,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc1,is) /=                             &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_v(l)%lambda_h(nzb_wall:nzb_wall+1,m) =                             &
                       building_surface_pars_f%pars(ind_s_tc1,is)
                       surf_usm_v(l)%lambda_h_green(nzb_wall:nzb_wall+1,m) =                       &
                       building_surface_pars_f%pars(ind_s_tc1,is)
                       surf_usm_v(l)%lambda_h_window(nzb_wall:nzb_wall+1,m) =                      &
                       building_surface_pars_f%pars(ind_s_tc1,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc2,is) /=                             &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_v(l)%lambda_h(nzb_wall+2,m) =                                      &
                       building_surface_pars_f%pars(ind_s_tc2,is)
                       surf_usm_v(l)%lambda_h_green(nzb_wall+2,m) =                                &
                       building_surface_pars_f%pars(ind_s_tc2,is)
                       surf_usm_v(l)%lambda_h_window(nzb_wall+2,m) =                               &
                       building_surface_pars_f%pars(ind_s_tc2,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_tc3,is) /=                             &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_v(l)%lambda_h(nzb_wall+3,m) =                                      &
                       building_surface_pars_f%pars(ind_s_tc3,is)
                       surf_usm_v(l)%lambda_h_green(nzb_wall+3,m) =                                &
                       building_surface_pars_f%pars(ind_s_tc3,is)
                       surf_usm_v(l)%lambda_h_window(nzb_wall+3,m) =                               &
                       building_surface_pars_f%pars(ind_s_tc3,is)
                    ENDIF

                    IF ( building_surface_pars_f%pars(ind_s_indoor_target_temp_summer,is) /=       &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%target_temp_summer(m) =                                       &
                       building_surface_pars_f%pars(ind_s_indoor_target_temp_summer,is)

                    IF ( building_surface_pars_f%pars(ind_s_indoor_target_temp_winter,is) /=       &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%target_temp_winter(m) =                                       &
                       building_surface_pars_f%pars(ind_s_indoor_target_temp_winter,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_wall,is) /=                       &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%emissivity(m,ind_veg_wall) =                                  &
                       building_surface_pars_f%pars(ind_s_emis_wall,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_green,is) /=                      &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%emissivity(m,ind_pav_green) =                                 &
                       building_surface_pars_f%pars(ind_s_emis_green,is)

                    IF ( building_surface_pars_f%pars(ind_s_emis_win,is) /=                        &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%emissivity(m,ind_wat_win) =                                   &
                       building_surface_pars_f%pars(ind_s_emis_win,is)

                    IF ( building_surface_pars_f%pars(ind_s_trans,is) /=                           &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%transmissivity(m) =                                           &
                       building_surface_pars_f%pars(ind_s_trans,is)

                    IF ( building_surface_pars_f%pars(ind_s_z0,is) /=                              &
                         building_surface_pars_f%fill )                                            &
                       surf_usm_v(l)%z0(m) = building_surface_pars_f%pars(ind_s_z0,is)

                    IF ( building_surface_pars_f%pars(ind_s_z0qh,is) /=                            &
                         building_surface_pars_f%fill )  THEN
                       surf_usm_v(l)%z0q(m) = building_surface_pars_f%pars(ind_s_z0qh,is)
                       surf_usm_v(l)%z0h(m) = building_surface_pars_f%pars(ind_s_z0qh,is)
                    ENDIF

                    EXIT ! Surface was found and processed
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
!
!--  Initialize albedo type via given type from static input file. Please note, even though
!--  the albedo type has been already given by the pars, albedo_type overwrites these values.
     IF ( albedo_type_f%from_file )  THEN
        DO  l = 0, 1
           DO  m = 1, surf_usm_h(l)%ns
              i = surf_usm_h(l)%i(m)
              j = surf_usm_h(l)%j(m)
              IF ( albedo_type_f%var(j,i) /= albedo_type_f%fill )                               &
                 surf_usm_h(l)%albedo_type(m,:) = albedo_type_f%var(j,i)
           ENDDO
        ENDDO
        DO  l = 0, 3
           DO  m = 1, surf_usm_v(l)%ns
              i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
              j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff
             IF ( albedo_type_f%var(j,i) /= albedo_type_f%fill )                                &
                 surf_usm_v(l)%albedo_type(m,:) = albedo_type_f%var(j,i)
           ENDDO
        ENDDO
     ENDIF
!
!--  Run further checks to ensure that the respecitve material fractions are prescribed properly.
!--  Start with horizontal surfaces (roofs).
     relative_fractions_corrected = .FALSE.
     DO  l = 0, 1
        DO  m = 1, surf_usm_h(l)%ns
           sum_frac = SUM( surf_usm_h(l)%frac(m,:) )
           IF ( sum_frac /= 1.0_wp )  THEN
              relative_fractions_corrected = .TRUE.
!
!--           Normalize relative fractions to 1. Deviations from 1 can arise, e.g. by rounding errors
!--           but also by inconsistent driver creation.
              IF ( sum_frac /= 0.0_wp )  THEN
                 surf_usm_h(l)%frac(m,:) = surf_usm_h(l)%frac(m,:) / sum_frac
!
!--           In case all relative fractions are erroneously set to zero, set wall fraction to 1.
              ELSE
                 surf_usm_h(l)%frac(m,ind_veg_wall)  = 1.0_wp
                 surf_usm_h(l)%frac(m,ind_wat_win)   = 0.0_wp
                 surf_usm_h(l)%frac(m,ind_pav_green) = 0.0_wp
              ENDIF
           ENDIF
        ENDDO
     ENDDO
!
!--  If fractions were normalized, give an informative message.
#if defined( __parallel )
     CALL MPI_ALLREDUCE( MPI_IN_PLACE, relative_fractions_corrected, 1,                            &
                         MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#endif
     IF ( relative_fractions_corrected )  THEN
        message_string = 'At some horizotal surfaces the relative material fractions do not ' //   &
                         'sum-up to one . Hence, the respective fractions were normalized.'
        CALL message( 'urban_surface_model_mod', 'PA0686', 0, 0, 0, 6, 0 )
     ENDIF
!
!--  Check relative fractions at vertical surfaces.
     relative_fractions_corrected = .FALSE.
     DO  l = 0, 3
        DO  m = 1, surf_usm_v(l)%ns
           sum_frac = SUM( surf_usm_v(l)%frac(m,:) )
           IF ( sum_frac /= 1.0_wp )  THEN
              relative_fractions_corrected = .TRUE.
!
!--           Normalize relative fractions to 1.
              IF ( sum_frac /= 0.0_wp )  THEN
                 surf_usm_v(l)%frac(m,:) = surf_usm_v(l)%frac(m,:) / sum_frac
!
!--           In case all relative fractions are erroneously set to zero, set wall fraction to 1.
              ELSE
                 surf_usm_v(l)%frac(m,ind_veg_wall)  = 1.0_wp
                 surf_usm_v(l)%frac(m,ind_wat_win)   = 0.0_wp
                 surf_usm_v(l)%frac(m,ind_pav_green) = 0.0_wp
              ENDIF
           ENDIF
        ENDDO
     ENDDO
!
!--  Also here, if fractions were normalized, give an informative message.
#if defined( __parallel )
     CALL MPI_ALLREDUCE( MPI_IN_PLACE, relative_fractions_corrected, 1,                            &
                         MPI_LOGICAL, MPI_LOR, comm2d, ierr )
#endif
     IF ( relative_fractions_corrected )  THEN
        message_string = 'At some vertical surfaces the relative material fractions do not ' //    &
                         'sum-up to one . Hence, the respective fractions were normalized.'
        CALL message( 'urban_surface_model_mod', 'PA0686', 0, 0, 0, 6, 0 )
     ENDIF
!
!--  Initialization of the wall/roof materials
     CALL usm_init_wall_heat_model()

!--  Init skin layer properties (can be done after initialization of wall layers)
     DO  l = 0, 1
        DO  m = 1, surf_usm_h(l)%ns
           i = surf_usm_h(l)%i(m)
           j = surf_usm_h(l)%j(m)

            surf_usm_h(l)%c_surface(m)           = surf_usm_h(l)%rho_c_wall(nzb_wall,m) *          &
                                                surf_usm_h(l)%dz_wall(nzb_wall,m) * 0.25_wp
            surf_usm_h(l)%lambda_surf(m)         = surf_usm_h(l)%lambda_h(nzb_wall,m) *            &
                                                surf_usm_h(l)%ddz_wall(nzb_wall,m) * 2.0_wp
            surf_usm_h(l)%c_surface_green(m)     = surf_usm_h(l)%rho_c_wall(nzb_wall,m) *          &
                                                surf_usm_h(l)%dz_wall(nzb_wall,m) * 0.25_wp
            surf_usm_h(l)%lambda_surf_green(m)   = surf_usm_h(l)%lambda_h_green(nzb_wall,m) *      &
                                                surf_usm_h(l)%ddz_green(nzb_wall,m) * 2.0_wp
            surf_usm_h(l)%c_surface_window(m)    = surf_usm_h(l)%rho_c_window(nzb_wall,m) *        &
                                                surf_usm_h(l)%dz_window(nzb_wall,m) * 0.25_wp
            surf_usm_h(l)%lambda_surf_window(m)  = surf_usm_h(l)%lambda_h_window(nzb_wall,m) *     &
                                                surf_usm_h(l)%ddz_window(nzb_wall,m) * 2.0_wp
        ENDDO
     ENDDO

     DO  l = 0, 3
         DO  m = 1, surf_usm_v(l)%ns
            i = surf_usm_v(l)%i(m) + surf_usm_v(l)%ioff
            j = surf_usm_v(l)%j(m) + surf_usm_v(l)%joff

             surf_usm_v(l)%c_surface(m) = surf_usm_v(l)%rho_c_wall(nzb_wall,m) *                   &
                                          surf_usm_v(l)%dz_wall(nzb_wall,m) * 0.25_wp
             surf_usm_v(l)%lambda_surf(m) = surf_usm_v(l)%lambda_h(nzb_wall,m) *                   &
                                            surf_usm_v(l)%ddz_wall(nzb_wall,m) * 2.0_wp
             surf_usm_v(l)%c_surface_green(m) = surf_usm_v(l)%rho_c_green(nzb_wall,m) *            &
                                                surf_usm_v(l)%dz_green(nzb_wall,m) * 0.25_wp
             surf_usm_v(l)%lambda_surf_green(m) = surf_usm_v(l)%lambda_h_green(nzb_wall,m) *       &
                                                  surf_usm_v(l)%ddz_green(nzb_wall,m) * 2.0_wp
             surf_usm_v(l)%c_surface_window(m) = surf_usm_v(l)%rho_c_window(nzb_wall,m) *          &
                                                    surf_usm_v(l)%dz_window(nzb_wall,m) * 0.25_wp
             surf_usm_v(l)%lambda_surf_window(m) = surf_usm_v(l)%lambda_h_window(nzb_wall,m) *     &
                                                    surf_usm_v(l)%ddz_window(nzb_wall,m) * 2.0_wp
         ENDDO
     ENDDO

!
!-- Check for consistent initialization.
!-- Check if roughness length for momentum, or heat, exceed surface-layer height and decrease local
!-- roughness length where necessary.
    DO  l = 0, 1
       DO  m = 1, surf_usm_h(l)%ns
          IF ( surf_usm_h(l)%z0(m) >= surf_usm_h(l)%z_mo(m) )  THEN

             surf_usm_h(l)%z0(m) = 0.9_wp * surf_usm_h(l)%z_mo(m)

             WRITE( message_string, * ) 'z0 exceeds surface-layer height at horizontal urban ' //     &
                                        'surface and is decreased appropriately at grid point ' //    &
                                        '(i,j) = ',  surf_usm_h(l)%i(m), surf_usm_h(l)%j(m)
             CALL message( 'urban_surface_model_mod', 'PA0503', 0, 0, myid, 6, 0 )
          ENDIF
          IF ( surf_usm_h(l)%z0h(m) >= surf_usm_h(l)%z_mo(m) )  THEN

             surf_usm_h(l)%z0h(m) = 0.9_wp * surf_usm_h(l)%z_mo(m)
             surf_usm_h(l)%z0q(m) = 0.9_wp * surf_usm_h(l)%z_mo(m)

             WRITE( message_string, * ) 'z0h exceeds surface-layer height at horizontal urban ' //    &
                                        'surface and is decreased appropriately at grid point ' //    &
                                        '(i,j) = ', surf_usm_h(l)%i(m), surf_usm_h(l)%j(m)
             CALL message( 'urban_surface_model_mod', 'PA0507', 0, 0, myid, 6, 0 )
          ENDIF
       ENDDO
    ENDDO

    DO  l = 0, 3
       DO  m = 1, surf_usm_v(l)%ns
          IF ( surf_usm_v(l)%z0(m) >= surf_usm_v(l)%z_mo(m) )  THEN

             surf_usm_v(l)%z0(m) = 0.9_wp * surf_usm_v(l)%z_mo(m)

             WRITE( message_string, * ) 'z0 exceeds surface-layer height at vertical urban ' //    &
                                        'surface and is decreased appropriately at grid point ' // &
                                        '(i,j) = ', surf_usm_v(l)%i(m)+surf_usm_v(l)%ioff,         &
                                         surf_usm_v(l)%j(m)+surf_usm_v(l)%joff
             CALL message( 'urban_surface_model_mod', 'PA0503', 0, 0, myid, 6, 0 )
          ENDIF
          IF ( surf_usm_v(l)%z0h(m) >= surf_usm_v(l)%z_mo(m) )  THEN

             surf_usm_v(l)%z0h(m) = 0.9_wp * surf_usm_v(l)%z_mo(m)
             surf_usm_v(l)%z0q(m) = 0.9_wp * surf_usm_v(l)%z_mo(m)

             WRITE( message_string, * ) 'z0h exceeds surface-layer height at vertical urban ' //   &
                                        'surface and is decreased appropriately at grid point ' // &
                                        '(i,j) = ', surf_usm_v(l)%i(m)+surf_usm_v(l)%ioff,         &
                                        surf_usm_v(l)%j(m)+surf_usm_v(l)%joff
             CALL message( 'urban_surface_model_mod', 'PA0507', 0, 0, myid, 6, 0 )
          ENDIF
       ENDDO
    ENDDO
!
!--  Intitialization of the surface and wall/ground/roof temperature
!
!--  Initialization for restart runs
     IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

!
!--     At horizontal surfaces. Please note, t_surf_wall_h is defined on a different data type,
!--     but with the same dimension.
         DO  l = 0, 1
            DO  m = 1, surf_usm_h(l)%ns
               i = surf_usm_h(l)%i(m)
               j = surf_usm_h(l)%j(m)
               k = surf_usm_h(l)%k(m)

               t_surf_wall_h(l)%val(m) = pt(k,j,i) * exner(k)
               t_surf_window_h(l)%val(m) = pt(k,j,i) * exner(k)
               t_surf_green_h(l)%val(m) = pt(k,j,i) * exner(k)
               surf_usm_h(l)%pt_surface(m) = pt(k,j,i) * exner(k)
            ENDDO
         ENDDO
!
!--      At vertical surfaces.
         DO  l = 0, 3
            DO  m = 1, surf_usm_v(l)%ns
               i = surf_usm_v(l)%i(m)
               j = surf_usm_v(l)%j(m)
               k = surf_usm_v(l)%k(m)

               t_surf_wall_v(l)%val(m) = pt(k,j,i) * exner(k)
               t_surf_window_v(l)%val(m) = pt(k,j,i) * exner(k)
               t_surf_green_v(l)%val(m) = pt(k,j,i) * exner(k)
               surf_usm_v(l)%pt_surface(m) = pt(k,j,i) * exner(k)
            ENDDO
         ENDDO

!
!--      For the sake of correct initialization, set also q_surface.
!--      Note, at urban surfaces q_surface is initialized with 0.
         IF ( humidity )  THEN
            DO  l = 0, 1
               DO  m = 1, surf_usm_h(l)%ns
                  surf_usm_h(l)%q_surface(m) = 0.0_wp
               ENDDO
            ENDDO
            DO  l = 0, 3
               DO  m = 1, surf_usm_v(l)%ns
                  surf_usm_v(l)%q_surface(m) = 0.0_wp
               ENDDO
            ENDDO
         ENDIF
!
!--      Initial values for t_wall
!--      Outer value is set to surface temperature, inner value is set to wall_inner_temperature
!--      and profile is logaritmic (linear in nz).
!--      Horizontal surfaces
         DO  l = 0, 1
            DO  m = 1, surf_usm_h(l)%ns
!
!--            Roof
               IF ( surf_usm_h(l)%isroof_surf(m) )  THEN
                   tin = roof_inner_temperature
                   twin = window_inner_temperature
!
!--            Normal land surface
               ELSE
                   tin = soil_inner_temperature
                   twin = window_inner_temperature
               ENDIF

               DO k = nzb_wall, nzt_wall+1
                   c = REAL( k - nzb_wall, wp ) / REAL( nzt_wall + 1 - nzb_wall , wp )

                   t_wall_h(l)%val(k,m) = ( 1.0_wp - c ) * t_surf_wall_h(l)%val(m) + c * tin
                   t_window_h(l)%val(k,m) = ( 1.0_wp - c ) * t_surf_window_h(l)%val(m) + c * twin
                   t_green_h(l)%val(k,m) = t_surf_wall_h(l)%val(m)
                   swc_h(l)%val(k,m) = 0.5_wp
                   swc_sat_h(l)%val(k,m) = 0.95_wp
                   swc_res_h(l)%val(k,m) = 0.05_wp
                   rootfr_h(l)%val(k,m) = 0.1_wp
                   wilt_h(l)%val(k,m) = 0.1_wp
                   fc_h(l)%val(k,m) = 0.9_wp
               ENDDO
            ENDDO
         ENDDO
!
!--      Vertical surfaces
         DO  l = 0, 3
            DO  m = 1, surf_usm_v(l)%ns
!
!--            Inner wall
               tin = wall_inner_temperature
               twin = window_inner_temperature

               DO k = nzb_wall, nzt_wall+1
                  c = REAL( k - nzb_wall, wp ) / REAL( nzt_wall + 1 - nzb_wall , wp )
                  t_wall_v(l)%val(k,m) = ( 1.0_wp - c ) * t_surf_wall_v(l)%val(m) + c * tin
                  t_window_v(l)%val(k,m) = ( 1.0_wp - c ) * t_surf_window_v(l)%val(m) + c * twin
                  t_green_v(l)%val(k,m) = t_surf_wall_v(l)%val(m)
               ENDDO
            ENDDO
         ENDDO
     ENDIF

!--
!--  Possibly DO user-defined actions (e.g. define heterogeneous wall surface)
     CALL user_init_urban_surface

!
!--  Initialize prognostic values for the first timestep
     t_surf_wall_h_p = t_surf_wall_h
     t_surf_wall_v_p = t_surf_wall_v
     t_surf_window_h_p = t_surf_window_h
     t_surf_window_v_p = t_surf_window_v
     t_surf_green_h_p = t_surf_green_h
     t_surf_green_v_p = t_surf_green_v

     t_wall_h_p = t_wall_h
     t_wall_v_p = t_wall_v
     t_window_h_p = t_window_h
     t_window_v_p = t_window_v
     t_green_h_p = t_green_h
     t_green_v_p = t_green_v

!
!-- Set initial values for prognostic soil quantities
    DO  l = 0, 1
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          m_liq_usm_h(l)%val  = 0.0_wp
       ENDIF
       m_liq_usm_h_p(l)%val = m_liq_usm_h(l)%val
!
!--    Set initial values for prognostic quantities
!--    Horizontal surfaces
       surf_usm_h(l)%c_liq     = 0.0_wp
       surf_usm_h(l)%qsws_liq  = 0.0_wp
       surf_usm_h(l)%qsws_veg  = 0.0_wp
    ENDDO

!
!-- Do the same for vertical surfaces
    DO  l = 0, 3
       surf_usm_v(l)%c_liq     = 0.0_wp
       surf_usm_v(l)%qsws_liq  = 0.0_wp
       surf_usm_v(l)%qsws_veg  = 0.0_wp
    ENDDO



    CALL cpu_log( log_point_s(78), 'usm_init', 'stop' )

    IF ( debug_output )  CALL debug_message( 'usm_init', 'end' )

 END SUBROUTINE usm_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Wall model as part of the urban surface model. The model predicts vertical and horizontal
!> wall / roof temperatures and window layer temperatures. No window layer temperature calculactions
!> during spinup to increase possible timestep.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_wall_heat_model( horizontal, l, during_spinup )


    IMPLICIT NONE

    LOGICAL       ::  during_spinup  !< if true, no calculation of window temperatures
    LOGICAL       ::  horizontal     !< Flag indicating horizontal or vertical surfaces

    INTEGER(iwp)  ::  kw             !< grid index - wall depth
    INTEGER(iwp)  ::  l              !< direction index
    INTEGER(iwp)  ::  m              !< running index for surface elements

    REAL(wp)  ::  win_absorp        !< absorption coefficient from transmissivity
    REAL(wp)  ::  win_nonrefl_1side !< non-reflected fraction after outer glass boundary

    REAL(wp), DIMENSION(nzb_wall:nzt_wall)  ::  wall_mod        !<
    REAL(wp), DIMENSION(nzb_wall:nzt_wall)  ::  wtend, wintend  !< tendency

    TYPE(surf_type),        POINTER         ::  surf            !< surface-date type variable
    TYPE(surf_type_2d_usm), POINTER         ::  t_wall
    TYPE(surf_type_2d_usm), POINTER         ::  t_wall_p
    TYPE(surf_type_2d_usm), POINTER         ::  t_window
    TYPE(surf_type_2d_usm), POINTER         ::  t_window_p
    TYPE(surf_type_2d_usm), POINTER         ::  t_green

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'usm_wall_heat_model: ', horizontal, l, during_spinup
       CALL debug_message( debug_string, 'start' )
    ENDIF

    wall_mod=1.0_wp
    IF ( usm_wall_mod  .AND.  during_spinup )  THEN
       DO  kw=nzb_wall, nzb_wall+1
          wall_mod(kw) = 0.1_wp
       ENDDO
    ENDIF

    IF ( horizontal )  THEN
       surf         => surf_usm_h(l)
       t_wall       => t_wall_h(l)
       t_wall_p     => t_wall_h_p(l)
       t_window     => t_window_h(l)
       t_window_p   => t_window_h_p(l)
       t_green      => t_green_h(l)
    ELSE
       surf         => surf_usm_v(l)
       t_wall       => t_wall_v(l)
       t_wall_p     => t_wall_v_p(l)
       t_window     => t_window_v(l)
       t_window_p   => t_window_v_p(l)
       t_green      => t_green_v(l)
    ENDIF
!
!-- Cycle for all surfaces in given direction
    !$OMP PARALLEL DO PRIVATE (m, kw, wtend, wintend, win_absorp) SCHEDULE (STATIC)
    DO  m = 1, surf%ns
!
!--    Prognostic equation for ground/roof temperature t_wall
       wtend(:) = 0.0_wp
       wtend(nzb_wall) = ( 1.0_wp / surf%rho_c_wall(nzb_wall,m) )                                 &
                          * ( surf%lambda_h(nzb_wall,m) * wall_mod(nzb_wall)                      &
                              * ( t_wall%val(nzb_wall+1,m) - t_wall%val(nzb_wall,m) )             &
                              * surf%ddz_wall(nzb_wall+1,m)                                       &
                              + surf%frac(m,ind_veg_wall)                                         &
                              / ( surf%frac(m,ind_veg_wall)                                       &
                                  + surf%frac(m,ind_pav_green) )                                  &
                              * surf%wghf_eb(m)                                                   &
                              - surf%frac(m,ind_pav_green)                                        &
                              / ( surf%frac(m,ind_veg_wall)                                       &
                                  + surf%frac(m,ind_pav_green) )                                  &
                              * ( surf%lambda_h_green(nzt_wall,m)                                 &
                              * wall_mod(nzt_wall)                                                &
                              * surf%ddz_green(nzt_wall,m)                                        &
                              + surf%lambda_h(nzb_wall,m)                                         &
                              * wall_mod(nzb_wall)                                                &
                              * surf%ddz_wall(nzb_wall,m) )                                       &
                              / ( surf%ddz_green(nzt_wall,m)                                      &
                              + surf%ddz_wall(nzb_wall,m) )                                       &
                              * ( t_wall%val(nzb_wall,m) - t_green%val(nzt_wall,m) )              &
                            ) * surf%ddz_wall_stag(nzb_wall,m)

!
!--    If indoor model is used inner wall layer is calculated by using iwghf (indoor
!--    wall ground heat flux)
       IF ( .NOT. indoor_model ) THEN
          surf%iwghf_eb(m) = surf%lambda_h(nzt_wall,m)  * wall_mod(nzt_wall)          &
                             * ( t_wall%val(nzt_wall+1,m) - t_wall%val(nzt_wall,m) )  &
                             * surf%ddz_wall(nzt_wall+1,m)
       ENDIF

       DO  kw = nzb_wall+1, nzt_wall-1
          wtend(kw) = ( 1.0_wp / surf%rho_c_wall(kw,m) )                                          &
                      * ( surf%lambda_h(kw,m) * wall_mod(kw)                                      &
                      * ( t_wall%val(kw+1,m) - t_wall%val(kw,m) )                                 &
                      * surf%ddz_wall(kw+1,m)                                                     &
                      - surf%lambda_h(kw-1,m) * wall_mod(kw-1)                                    &
                      * ( t_wall%val(kw,m) - t_wall%val(kw-1,m) )                                 &
                      * surf%ddz_wall(kw,m)                                                       &
                      ) * surf%ddz_wall_stag(kw,m)
       ENDDO
       wtend(nzt_wall) = ( 1.0_wp / surf%rho_c_wall(nzt_wall,m) )                                 &
                         * ( -surf%lambda_h(nzt_wall-1,m) * wall_mod(nzt_wall-1)                  &
                             * ( t_wall%val(nzt_wall,m) - t_wall%val(nzt_wall-1,m) )              &
                              * surf%ddz_wall(nzt_wall,m) + surf%iwghf_eb(m)                      &
                           ) * surf%ddz_wall_stag(nzt_wall,m)

       t_wall_p%val(nzb_wall:nzt_wall,m) = t_wall%val(nzb_wall:nzt_wall,m) + dt_3d                &
                                         * ( tsc(2) * wtend(nzb_wall:nzt_wall) + tsc(3)           &
                                             * surf%tt_wall_m(nzb_wall:nzt_wall,m) )

!
!--    During spinup the tempeature inside window layers is not calculated to make larger timesteps possible
       IF ( .NOT. during_spinup )  THEN
!
!--       Reflectivity in glass windows is considered as equal on frontal and rear side of the
!--       glass, which together make total reflectivity (albedo for win fraction).
          win_nonrefl_1side = 1.0 - ( surf%albedo(m,ind_wat_win) + surf%transmissivity(m)         &
                                      + 1.0_wp                                                    &
                                      - sqrt( ( surf%albedo(m,ind_wat_win)                        &
                                                + surf%transmissivity(m) + 1.0_wp ) ** 2          &
                                              - 4 * surf%albedo(m,ind_wat_win) ) ) / 2.0_wp
!
!--       Absorption coefficient is calculated using zw from internal tranmissivity, which only
!--       considers absorption without the effects of reflection.
          win_absorp = -log( ( surf%transmissivity(m) + surf%albedo(m,ind_wat_win) - 1.0_wp       &
                               + win_nonrefl_1side ) / win_nonrefl_1side )                        &
                       / surf%zw_window(nzt_wall,m)
!
!--       Prognostic equation for ground/roof window temperature t_window takes absorption
!--       of shortwave radiation into account
          wintend(:) = 0.0_wp
          wintend(nzb_wall) = ( 1.0_wp / surf%rho_c_window(nzb_wall,m) )                          &
                                * ( surf%lambda_h_window(nzb_wall,m)                              &
                                * ( t_window%val(nzb_wall+1,m) - t_window%val(nzb_wall,m) )       &
                                * surf%ddz_window(nzb_wall+1,m)                                   &
                                + surf%wghf_eb_window(m)                                          &
                                + surf%rad_sw_in(m) * win_nonrefl_1side                           &
                                * ( 1.0_wp - exp( -win_absorp                                     &
                                    * surf%zw_window(nzb_wall,m) ) )                              &
                              ) * surf%ddz_window_stag(nzb_wall,m)


          IF ( .NOT. indoor_model )  THEN
             surf%iwghf_eb_window(m) = surf%lambda_h_window(nzt_wall,m)                    &
                               * ( t_window%val(nzt_wall+1,m) - t_window%val(nzt_wall,m) ) &
                               * surf%ddz_window(nzt_wall+1,m)
          ENDIF


          DO  kw = nzb_wall+1, nzt_wall-1
             wintend(kw) = ( 1.0_wp / surf%rho_c_window(kw,m) )                                   &
                           * ( surf%lambda_h_window(kw,m)                                         &
                           * ( t_window%val(kw+1,m) - t_window%val(kw,m) )                        &
                           * surf%ddz_window(kw+1,m)                                              &
                           - surf%lambda_h_window(kw-1,m)                                         &
                           * ( t_window%val(kw,m) - t_window%val(kw-1,m) )                        &
                           * surf%ddz_window(kw,m)                                                &
                           + surf%rad_sw_in(m) * win_nonrefl_1side                                &
                           * ( exp( -win_absorp * surf%zw_window(kw-1,m) )                        &
                             - exp( -win_absorp * surf%zw_window(kw,  m) ) )                      &
                             ) * surf%ddz_window_stag(kw,m)

          ENDDO

          wintend(nzt_wall) = ( 1.0_wp / surf%rho_c_window(nzt_wall,m) )                          &
                                 * ( -surf%lambda_h_window(nzt_wall-1,m)                          &
                                 * ( t_window%val(nzt_wall,m) - t_window%val(nzt_wall-1,m) )      &
                                 * surf%ddz_window(nzt_wall,m)                                    &
                                 + surf%iwghf_eb_window(m)                                        &
                                 + surf%rad_sw_in(m) * win_nonrefl_1side                          &
                                 * ( exp( -win_absorp * surf%zw_window(nzt_wall-1,m) )            &
                                   - exp( -win_absorp * surf%zw_window(nzt_wall,  m) ) )          &
                                   ) * surf%ddz_window_stag(nzt_wall,m)

          t_window_p%val(nzb_wall:nzt_wall,m) = t_window%val(nzb_wall:nzt_wall,m) + dt_3d         &
                                              * ( tsc(2) * wintend(nzb_wall:nzt_wall) + tsc(3)    &
                                              * surf%tt_window_m(nzb_wall:nzt_wall,m) )

       ENDIF
!
!--    Calculate t_wall tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
           IF ( intermediate_timestep_count == 1 )  THEN
              DO  kw = nzb_wall, nzt_wall
                 surf%tt_wall_m(kw,m) = wtend(kw)
              ENDDO
           ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
               DO  kw = nzb_wall, nzt_wall
                  surf%tt_wall_m(kw,m) = -9.5625_wp * wtend(kw) +                                 &
                                               5.3125_wp * surf%tt_wall_m(kw,m)
               ENDDO
           ENDIF
       ENDIF

       IF ( .NOT. during_spinup )  THEN
!
!--       Calculate t_window tendencies for the next Runge-Kutta step
          IF ( timestep_scheme(1:5) == 'runge' )  THEN
              IF ( intermediate_timestep_count == 1 )  THEN
                 DO  kw = nzb_wall, nzt_wall
                    surf%tt_window_m(kw,m) = wintend(kw)
                 ENDDO
              ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
                  DO  kw = nzb_wall, nzt_wall
                     surf%tt_window_m(kw,m) = -9.5625_wp * wintend(kw) +                          &
                                                    5.3125_wp * surf%tt_window_m(kw,m)
                  ENDDO
              ENDIF
          ENDIF
       ENDIF
    ENDDO

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'usm_wall_heat_model: ', horizontal, l, during_spinup
       CALL debug_message( debug_string, 'end' )
    ENDIF

 END SUBROUTINE usm_wall_heat_model

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!
!> Green and substrate model as part of the urban surface model. The model predicts ground
!> temperatures.
!>
!> Important: green-heat model crashes due to unknown reason. Green fraction is thus set to zero
!> (in favor of wall fraction).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_green_heat_model( horizontal, l )


    IMPLICIT NONE

    LOGICAL       ::  horizontal          !< Flag indicating horizontal or vertical surfaces
    INTEGER(iwp)  ::  l                   !< direction index

    INTEGER(iwp)  ::  i, j, k, kw, m      !< running indices

    LOGICAL  ::  conserve_water_content = .TRUE.  !<

    REAL(wp)  ::  drho_l_lv               !< frequently used parameter
    REAL(wp)  ::  h_vg                    !< Van Genuchten coef. h
    REAL(wp)  ::  ke, lambda_h_green_sat  !< heat conductivity for saturated soil

    REAL(wp), DIMENSION(nzb_wall:nzt_wall)   ::  gtend,tend         !< tendency
    REAL(wp), DIMENSION(nzb_wall:nzt_wall)   ::  root_extr_green    !<

    REAL(wp), DIMENSION(nzb_wall:nzt_wall+1) ::  gamma_green_temp   !< temp. gamma
    REAL(wp), DIMENSION(nzb_wall:nzt_wall+1) ::  lambda_green_temp  !< temp. lambda

    TYPE(surf_type), POINTER                 ::  surf               !< surface-date type variable
    TYPE(surf_type_2d_usm), POINTER          ::  t_wall
    TYPE(surf_type_2d_usm), POINTER          ::  t_green

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'usm_green_heat_model: ', horizontal, l
       CALL debug_message( debug_string, 'start' )
    ENDIF

    drho_l_lv = 1.0_wp / (rho_l * l_v)

!
!-- For horizontal upward surfaces.
    IF ( horizontal .AND. l==0 )  THEN
       surf => surf_usm_h(l)
       t_wall       => t_wall_h(l)
       t_green      => t_green_h(l)

!--    Set tendency array for soil moisture to zero
       IF ( surf%ns > 0 )  THEN
          IF ( intermediate_timestep_count == 1 )  surf%tswc_h_m = 0.0_wp
       ENDIF

       !$OMP PARALLEL DO PRIVATE (m, i, j, k, kw, lambda_h_green_sat, ke, lambda_green_temp,      &
       !$OMP&  gtend, tend, h_vg, gamma_green_temp, m_total, root_extr_green) SCHEDULE (STATIC)
       DO  m = 1, surf%ns
          IF (surf%frac(m,ind_pav_green) > 0.0_wp)  THEN
!
!--         Obtain indices
             i = surf%i(m)
             j = surf%j(m)
             k = surf%k(m)

             DO  kw = nzb_wall, nzt_wall
!
!--             Calculate volumetric heat capacity of the soil, taking into account water content
                surf%rho_c_total_green(kw,m) = (surf%rho_c_green(kw,m)                            &
                                                     * (1.0_wp - swc_sat_h(l)%val(kw,m))          &
                                                     + rho_c_water * swc_h(l)%val(kw,m))

!
!--             Calculate soil heat conductivity at the center of the soil layers
                lambda_h_green_sat = lambda_h_green_sm ** ( 1.0_wp - swc_sat_h(l)%val(kw,m) )     &
                                     * lambda_h_water ** swc_h(l)%val(kw,m)

                ke = 1.0_wp + LOG10( MAX( 0.1_wp,swc_h(l)%val(kw,m) / swc_sat_h(l)%val(kw,m) ) )

                lambda_green_temp(kw) = ke * (lambda_h_green_sat - lambda_h_green_dry)            &
                                        + lambda_h_green_dry

             ENDDO
             lambda_green_temp(nzt_wall+1) = lambda_green_temp(nzt_wall)


!
!--          Calculate soil heat conductivity (lambda_h) at the _stag level using linear interpolation.
!--          For pavement surface, the true pavement depth is considered
             DO  kw = nzb_wall, nzt_wall
                surf%lambda_h_green(kw,m) = ( lambda_green_temp(kw+1) + lambda_green_temp(kw) )   &
                                                  * 0.5_wp
             ENDDO

             t_green_h(l)%val(nzt_wall+1,m) = t_wall_h(l)%val(nzb_wall,m)
!
!--          Prognostic equation for ground/roof temperature t_green_h
             gtend(:) = 0.0_wp
             gtend(nzb_wall) = ( 1.0_wp / surf%rho_c_total_green(nzb_wall,m) )                    &
                               * ( surf%lambda_h_green(nzb_wall,m)                                &
                                   * ( t_green_h(l)%val(nzb_wall+1,m)                             &
                                   - t_green_h(l)%val(nzb_wall,m) )                               &
                                   * surf%ddz_green(nzb_wall+1,m)                                 &
                                   + surf%wghf_eb_green(m)                                        &
                                 ) * surf%ddz_green_stag(nzb_wall,m)

             DO  kw = nzb_wall+1, nzt_wall
                gtend(kw) = ( 1.0_wp / surf%rho_c_total_green(kw,m) )                             &
                            * ( surf%lambda_h_green(kw,m)                                         &
                                * ( t_green_h(l)%val(kw+1,m) - t_green_h(l)%val(kw,m) )           &
                                * surf%ddz_green(kw+1,m)                                          &
                                - surf%lambda_h_green(kw-1,m)                                     &
                                * ( t_green_h(l)%val(kw,m) - t_green_h(l)%val(kw-1,m) )           &
                                * surf%ddz_green(kw,m)                                            &
                              ) * surf%ddz_green_stag(kw,m)
             ENDDO

             t_green_h_p(l)%val(nzb_wall:nzt_wall,m) = t_green_h(l)%val(nzb_wall:nzt_wall,m)      &
                            + dt_3d * ( tsc(2) * gtend(nzb_wall:nzt_wall) + tsc(3)                &
                            * surf%tt_green_m(nzb_wall:nzt_wall,m) )


!
!--          Calculate t_green tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                 IF ( intermediate_timestep_count == 1 )  THEN
                    DO  kw = nzb_wall, nzt_wall
                       surf%tt_green_m(kw,m) = gtend(kw)
                    ENDDO
                 ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
                     DO  kw = nzb_wall, nzt_wall
                        surf%tt_green_m(kw,m) = -9.5625_wp * gtend(kw) + 5.3125_wp                 &
                                                      * surf%tt_green_m(kw,m)
                     ENDDO
                 ENDIF
             ENDIF

             DO  kw = nzb_wall, nzt_wall

!
!--             Calculate soil diffusivity at the center of the soil layers
                lambda_green_temp(kw) = ( - b_ch * surf%gamma_w_green_sat(kw,m) * psi_sat          &
                                          / swc_sat_h(l)%val(kw,m) )                               &
                                        * ( MAX( swc_h(l)%val(kw,m), wilt_h(l)%val(kw,m) )         &
                                          / swc_sat_h(l)%val(kw,m) )**( b_ch + 2.0_wp )

!
!--             Parametrization of Van Genuchten
                IF ( soil_type /= 7 )  THEN
!
!--                Calculate the hydraulic conductivity after Van Genuchten (1980)
                   h_vg = ( ( (swc_res_h(l)%val(kw,m) - swc_sat_h(l)%val(kw,m))                    &
                              / ( swc_res_h(l)%val(kw,m) -                                         &
                              MAX( swc_h(l)%val(kw,m), wilt_h(l)%val(kw,m) ) ) )**                 &
                              ( surf%n_vg_green(m) / (surf%n_vg_green(m) - 1.0_wp ) )              &
                              - 1.0_wp                                                             &
                          )** ( 1.0_wp / surf%n_vg_green(m) ) / surf%alpha_vg_green(m)


                   gamma_green_temp(kw) = surf%gamma_w_green_sat(kw,m)                             &
                                          * ( ( ( 1.0_wp + ( surf%alpha_vg_green(m) * h_vg )**     &
                                                surf%n_vg_green(m) )**                             &
                                                ( 1.0_wp - 1.0_wp / surf%n_vg_green(m) )           &
                                                - ( surf%alpha_vg_green(m) * h_vg )**              &
                                                ( surf%n_vg_green(m) - 1.0_wp) )**2                &
                                            ) / ( ( 1.0_wp + ( surf%alpha_vg_green(m) * h_vg )**   &
                                                  surf%n_vg_green(m) )**                           &
                                                  ( ( 1.0_wp  - 1.0_wp / surf%n_vg_green(m) )      &
                                                    *( surf%l_vg_green(m) + 2.0_wp) )              &
                                                )

!
!--             Parametrization of Clapp & Hornberger
                ELSE
                   gamma_green_temp(kw) = surf%gamma_w_green_sat(kw,m) * ( swc_h(l)%val(kw,m)      &
                                          / swc_sat_h(l)%val(kw,m) )**( 2.0_wp * b_ch + 3.0_wp )
                ENDIF

             ENDDO

!
!--          Prognostic equation for soil moisture content. Only performed, when humidity is enabled in
!--          the atmosphere
             IF ( humidity )  THEN
!
!--             Calculate soil diffusivity (lambda_w) at the _stag level using linear interpolation.
!--             To do: replace this with ECMWF-IFS Eq. 8.81
                DO  kw = nzb_wall, nzt_wall-1

                   surf%lambda_w_green(kw,m) = ( lambda_green_temp(kw+1)                          &
                                                       + lambda_green_temp(kw) )                  &
                                                     * 0.5_wp
                   surf%gamma_w_green(kw,m)  = ( gamma_green_temp(kw+1)                           &
                                                       + gamma_green_temp(kw) )                   &
                                                     * 0.5_wp

                ENDDO

!
!--             In case of a closed bottom (= water content is conserved), set hydraulic conductivity
!--             to zero so that no water will be lost in the bottom layer.
                IF ( conserve_water_content )  THEN
                   surf%gamma_w_green(kw,m) = 0.0_wp
                ELSE
                   surf%gamma_w_green(kw,m) = gamma_green_temp(nzt_wall)
                ENDIF

!--             The root extraction (= root_extr * qsws_veg / (rho_l * l_v)) ensures the mass
!--             conservation for water. The transpiration of plants equals the cumulative withdrawals
!--             by the roots in the soil. The scheme takes into account the availability of water in
!--             the soil layers as well as the root fraction in the respective layer. Layer with
!--             moisture below wilting point will not contribute, which reflects the preference of
!--             plants to take water from moister layers.

!
!--             Calculate the root extraction (ECMWF 7.69, the sum of root_extr = 1). The energy
!--             balance solver guarantees a positive transpiration, so that there is no need for an
!--             additional check.
                m_total = 0.0_wp
                DO  kw = nzb_wall, nzt_wall
                    IF ( swc_h(l)%val(kw,m) > wilt_h(l)%val(kw,m) )  THEN
                       m_total = m_total + rootfr_h(l)%val(kw,m) * swc_h(l)%val(kw,m)
                    ENDIF
                ENDDO

                IF ( m_total > 0.0_wp )  THEN
                   DO  kw = nzb_wall, nzt_wall
                      IF ( swc_h(l)%val(kw,m) > wilt_h(l)%val(kw,m) )  THEN
                         root_extr_green(kw) = rootfr_h(l)%val(kw,m) * swc_h(l)%val(kw,m) / m_total
                      ELSE
                         root_extr_green(kw) = 0.0_wp
                      ENDIF
                   ENDDO
                ENDIF

!
!--             Prognostic equation for soil water content m_soil.
                tend(:) = 0.0_wp

                tend(nzb_wall) = ( surf_usm_h(l)%lambda_w_green(nzb_wall,m)                       &
                                 * ( swc_h(l)%val(nzb_wall+1,m) - swc_h(l)%val(nzb_wall,m) )      &
                                 * surf_usm_h(l)%ddz_green(nzb_wall+1,m)                          &
                                 - surf_usm_h(l)%gamma_w_green(nzb_wall,m)                        &
                                 - ( root_extr_green(nzb_wall) * surf_usm_h(l)%qsws_veg(m)        &
!                                    + surf_usm_h(l)%qsws_soil_green(m)                           &
                                   ) * drho_l_lv )                                                &
                                 * surf_usm_h(l)%ddz_green_stag(nzb_wall,m)

                DO  kw = nzb_wall+1, nzt_wall-1
                   tend(kw) = ( surf_usm_h(l)%lambda_w_green(kw,m)                                &
                                * ( swc_h(l)%val(kw+1,m) - swc_h(l)%val(kw,m) )                   &
                                * surf_usm_h(l)%ddz_green(kw+1,m)                                 &
                                - surf_usm_h(l)%gamma_w_green(kw,m)                               &
                                - surf_usm_h(l)%lambda_w_green(kw-1,m)                            &
                                * ( swc_h(l)%val(kw,m) - swc_h(l)%val(kw-1,m) )                   &
                                * surf_usm_h(l)%ddz_green(kw,m)                                   &
                                + surf_usm_h(l)%gamma_w_green(kw-1,m)                             &
                                - (root_extr_green(kw)                                            &
                                * surf_usm_h(l)%qsws_veg(m)                                       &
                                * drho_l_lv)                                                      &
                             ) * surf_usm_h(l)%ddz_green_stag(kw,m)

                ENDDO
                tend(nzt_wall) = ( - surf_usm_h(l)%gamma_w_green(nzt_wall,m)                      &
                                   - surf_usm_h(l)%lambda_w_green(nzt_wall-1,m)                   &
                                   * (swc_h(l)%val(nzt_wall,m)                                    &
                                   - swc_h(l)%val(nzt_wall-1,m))                                  &
                                   * surf_usm_h(l)%ddz_green(nzt_wall,m)                          &
                                   + surf_usm_h(l)%gamma_w_green(nzt_wall-1,m)                    &
                                   - ( root_extr_green(nzt_wall)                                  &
                                   * surf_usm_h(l)%qsws_veg(m)                                    &
                                   * drho_l_lv  )                                                 &
                                 ) * surf_usm_h(l)%ddz_green_stag(nzt_wall,m)

                swc_h_p(l)%val(nzb_wall:nzt_wall,m) = swc_h(l)%val(nzb_wall:nzt_wall,m) + dt_3d   &
                                               * ( tsc(2) * tend(:) + tsc(3)                      &
                                                   * surf_usm_h(l)%tswc_h_m(:,m)                  &
                                                  )

!
!--             Account for dry soils (find a better solution here!)
                DO  kw = nzb_wall, nzt_wall
                   IF ( swc_h_p(l)%val(kw,m) < 0.0_wp )  swc_h_p(l)%val(kw,m) = 0.0_wp
                ENDDO

!
!--             Calculate m_soil tendencies for the next Runge-Kutta step
                IF ( timestep_scheme(1:5) == 'runge' )  THEN
                   IF ( intermediate_timestep_count == 1 )  THEN
                      DO  kw = nzb_wall, nzt_wall
                         surf_usm_h(l)%tswc_h_m(kw,m) = tend(kw)
                      ENDDO
                   ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
                      DO  kw = nzb_wall, nzt_wall
                         surf_usm_h(l)%tswc_h_m(kw,m) = -9.5625_wp * tend(kw) + 5.3125_wp         &
                                                     * surf_usm_h(l)%tswc_h_m(kw,m)
                      ENDDO
                   ENDIF
                ENDIF
             ENDIF

          ENDIF
       ENDDO
    ELSE
       IF ( horizontal) THEN
!--       For horizontal downward surfaces
          surf         => surf_usm_h(l)
          t_wall       => t_wall_h(l)
          t_green      => t_green_h(l)
       ELSE
!--       For vertical surfaces
          surf         => surf_usm_v(l)
          t_wall       => t_wall_v(l)
          t_green      => t_green_v(l)
       ENDIF
       !$OMP PARALLEL DO PRIVATE (m, i, j, k, kw) SCHEDULE (STATIC)
       DO  m = 1, surf%ns
          IF (surf%frac(m,ind_pav_green) > 0.0_wp)  THEN
!
!-- No substrate layer for green walls / only groundbase green walls (ivy i.e.) -> Green layers get
!-- same temperature as first wall layer, therefore no temperature calculations for vertical green
!-- substrate layers now

!
! !
! !--          Obtain indices
!              i = surf%i(m)
!              j = surf%j(m)
!              k = surf%k(m)
!
!              t_green%val(nzt_wall+1,m) = t_wall%val(nzb_wall,m)
! !
! !--          Prognostic equation for green temperature t_green_v
!              gtend(:) = 0.0_wp
!              gtend(nzb_wall) = (1.0_wp / surf%rho_c_green(nzb_wall,m)) *                        &
!                                      ( surf%lambda_h_green(nzb_wall,m) *                        &
!                                        ( t_green%val(nzb_wall+1,m)                              &
!                                        - t_green%val(nzb_wall,m) ) *                            &
!                                        surf%ddz_green(nzb_wall+1,m)                             &
!                                      + surf%wghf_eb(m) ) *                                      &
!                                        surf%ddz_green_stag(nzb_wall,m)
!
!              DO  kw = nzb_wall+1, nzt_wall
!                 gtend(kw) = (1.0_wp / surf%rho_c_green(kw,m))                                   &
!                           * (   surf%lambda_h_green(kw,m)                                       &
!                             * ( t_green%val(kw+1,m) - t_green%val(kw,m) )                       &
!                             * surf%ddz_green(kw+1,m)                                            &
!                           - surf%lambda_h(kw-1,m)                                               &
!                             * ( t_green%val(kw,m) - t_green%val(kw-1,m) )                       &
!                             * surf%ddz_green(kw,m) )                                            &
!                           * surf%ddz_green_stag(kw,m)
!              ENDDO
!
!              t_green_v_p(l)%val(nzb_wall:nzt_wall,m) =                                          &
!                                   t_green%val(nzb_wall:nzt_wall,m)                              &
!                                 + dt_3d * ( tsc(2)                                              &
!                                 * gtend(nzb_wall:nzt_wall) + tsc(3)                             &
!                                 * surf%tt_green_m(nzb_wall:nzt_wall,m) )
!
! !
! !--          Calculate t_green tendencies for the next Runge-Kutta step
!              IF ( timestep_scheme(1:5) == 'runge' )  THEN
!                  IF ( intermediate_timestep_count == 1 )  THEN
!                     DO  kw = nzb_wall, nzt_wall
!                        surf%tt_green_m(kw,m) = gtend(kw)
!                     ENDDO
!                  ELSEIF ( intermediate_timestep_count <                                         &
!                           intermediate_timestep_count_max )  THEN
!                      DO  kw = nzb_wall, nzt_wall
!                         surf%tt_green_m(kw,m) =                                                 &
!                                     - 9.5625_wp * gtend(kw) +                                   &
!                                       5.3125_wp * surf%tt_green_m(kw,m)
!                      ENDDO
!                  ENDIF
!              ENDIF
             DO  kw = nzb_wall, nzt_wall+1
                 t_green%val(kw,m) = t_wall%val(nzb_wall,m)
             ENDDO
          ENDIF
       ENDDO
    ENDIF

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'usm_green_heat_model: ', horizontal, l
       CALL debug_message( debug_string, 'end' )
    ENDIF

 END SUBROUTINE usm_green_heat_model

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &usm_par for urban surface model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_parin

    IMPLICIT NONE

    CHARACTER(LEN=80)  ::  line  !< string containing current line of file PARIN

    NAMELIST /urban_surface_par/                                                                   &
                        building_type,                                                             &
                        roof_category,                                                             &
                        roof_inner_temperature,                                                    &
                        roughness_concrete,                                                        &
                        soil_inner_temperature,                                                    &
                        urban_surface,                                                             &
                        usm_wall_mod,                                                              &
                        wall_category,                                                             &
                        wall_inner_temperature,                                                    &
                        window_inner_temperature


    NAMELIST /urban_surface_parameters/                                                            &
                        building_type,                                                             &
                        roof_category,                                                             &
                        roof_inner_temperature,                                                    &
                        roughness_concrete,                                                        &
                        soil_inner_temperature,                                                    &
                        urban_surface,                                                             &
                        usm_wall_mod,                                                              &
                        wall_category,                                                             &
                        wall_inner_temperature,                                                    &
                        window_inner_temperature



!
!-- Try to find urban surface model package
    REWIND ( 11 )
    line = ' '
    DO WHILE ( INDEX( line, '&urban_surface_parameters' ) == 0 )
       READ ( 11, '(A)', END = 12 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read user-defined namelist
    READ ( 11, urban_surface_parameters, ERR = 10 )

!
!-- Set flag that indicates that the urban surface model is switched on
    urban_surface = .TRUE.

    GOTO 14

 10 BACKSPACE( 11 )
    READ( 11 , '(A)') line
    CALL parin_fail_message( 'urban_surface_parameters', line )
!
!-- Try to find old namelist
 12 REWIND ( 11 )
    line = ' '
    DO WHILE ( INDEX( line, '&urban_surface_par' ) == 0 )
       READ ( 11, '(A)', END = 14 )  line
    ENDDO
    BACKSPACE ( 11 )

!
!-- Read user-defined namelist
    READ ( 11, urban_surface_par, ERR = 13, END = 14 )

    message_string = 'namelist urban_surface_par is deprecated and will be removed in near ' //    &
                     'future. Please use namelist urban_surface_parameters instead'
    CALL message( 'usm_parin', 'PA0487', 0, 1, 0, 6, 0 )

!
!-- Set flag that indicates that the urban surface model is switched on
    urban_surface = .TRUE.

    GOTO 14

 13 BACKSPACE( 11 )
    READ( 11 , '(A)') line
    CALL parin_fail_message( 'urban_surface_par', line )


 14 CONTINUE


 END SUBROUTINE usm_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!> Soubroutine reads t_surf and t_wall.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxr_on_file, nynf, nyn_on_file,   &
                               nysf, nysc, nys_on_file, found )


    USE control_parameters,                                                                        &
        ONLY: length,                                                                              &
              restart_string

    IMPLICIT NONE

    INTEGER(iwp)  ::  k                 !< running index over previous input files covering current local domain
    INTEGER(iwp)  ::  l                 !< index variable for surface type
    INTEGER(iwp)  ::  nxlc              !< index of left boundary on current subdomain
    INTEGER(iwp)  ::  nxlf              !< index of left boundary on former subdomain
    INTEGER(iwp)  ::  nxl_on_file       !< index of left boundary on former local domain
    INTEGER(iwp)  ::  nxrf              !< index of right boundary on former subdomain
    INTEGER(iwp)  ::  nxr_on_file       !< index of right boundary on former local domain
    INTEGER(iwp)  ::  nynf              !< index of north boundary on former subdomain
    INTEGER(iwp)  ::  nyn_on_file       !< index of north boundary on former local domain
    INTEGER(iwp)  ::  nysc              !< index of south boundary on current subdomain
    INTEGER(iwp)  ::  nysf              !< index of south boundary on former subdomain
    INTEGER(iwp)  ::  nys_on_file       !< index of south boundary on former local domain
    INTEGER(iwp)  ::  ns_h_on_file_usm(0:1)  !< number of horizontal surface elements (urban type) on file
    INTEGER(iwp)  ::  ns_v_on_file_usm(0:3)  !< number of vertical surface elements (urban type) on file
!
!-- Note, the save attribute in the following array declaration is necessary, in order to keep the
!-- number of urban surface elements on file during rrd_local calls.
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  end_index_on_file    !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  start_index_on_file  !<

    LOGICAL, INTENT(OUT)  ::  found  !<

! MS: Why are there individual temporary arrays that all have the same size?
    TYPE( surf_type_1d_usm ), DIMENSION(0:1), SAVE ::  tmp_surf_green_h   !<
    TYPE( surf_type_1d_usm ), DIMENSION(0:1), SAVE ::  tmp_surf_mliq_h    !<
    TYPE( surf_type_1d_usm ), DIMENSION(0:1), SAVE ::  tmp_surf_wall_h    !<
    TYPE( surf_type_1d_usm ), DIMENSION(0:1), SAVE ::  tmp_surf_waste_h   !<
    TYPE( surf_type_1d_usm ), DIMENSION(0:1), SAVE ::  tmp_surf_window_h  !<

    TYPE( surf_type_2d_usm ), DIMENSION(0:1), SAVE ::  tmp_green_h   !<
    TYPE( surf_type_2d_usm ), DIMENSION(0:1), SAVE ::  tmp_wall_h    !<
    TYPE( surf_type_2d_usm ), DIMENSION(0:1), SAVE ::  tmp_window_h  !<

    TYPE( surf_type_1d_usm ), DIMENSION(0:3), SAVE ::  tmp_surf_green_v   !<
    TYPE( surf_type_1d_usm ), DIMENSION(0:3), SAVE ::  tmp_surf_wall_v    !<
    TYPE( surf_type_1d_usm ), DIMENSION(0:3), SAVE ::  tmp_surf_waste_v   !<
    TYPE( surf_type_1d_usm ), DIMENSION(0:3), SAVE ::  tmp_surf_window_v  !<

    TYPE( surf_type_2d_usm ), DIMENSION(0:3), SAVE ::  tmp_green_v   !<
    TYPE( surf_type_2d_usm ), DIMENSION(0:3), SAVE ::  tmp_wall_v    !<
    TYPE( surf_type_2d_usm ), DIMENSION(0:3), SAVE ::  tmp_window_v  !<


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'ns_h_on_file_usm')
          IF ( k == 1 )  THEN
             READ ( 13 ) ns_h_on_file_usm

             DO  l = 0, 1
                IF ( ALLOCATED( tmp_surf_wall_h(l)%val ) )    DEALLOCATE( tmp_surf_wall_h(l)%val )
                IF ( ALLOCATED( tmp_wall_h(l)%val ) )         DEALLOCATE( tmp_wall_h(l)%val )
                IF ( ALLOCATED( tmp_surf_window_h(l)%val ) )  DEALLOCATE( tmp_surf_window_h(l)%val )
                IF ( ALLOCATED( tmp_window_h(l)%val) )        DEALLOCATE( tmp_window_h(l)%val )
                IF ( ALLOCATED( tmp_surf_green_h(l)%val) )    DEALLOCATE( tmp_surf_green_h(l)%val )
                IF ( ALLOCATED( tmp_green_h(l)%val) )         DEALLOCATE( tmp_green_h(l)%val )
                IF ( ALLOCATED( tmp_surf_mliq_h(l)%val) )     DEALLOCATE( tmp_surf_mliq_h(l)%val )
                IF ( ALLOCATED( tmp_surf_waste_h(l)%val) )    DEALLOCATE( tmp_surf_waste_h(l)%val )
             ENDDO

!
!--          Allocate temporary arrays for reading data on file. Note, the size of allocated surface
!--          elements do not necessarily need  to match the size of present surface elements on
!--          current processor, as the number of processors between restarts can change.
             DO  l = 0, 1
                ALLOCATE( tmp_surf_wall_h(l)%val(1:ns_h_on_file_usm(l)) )
                ALLOCATE( tmp_wall_h(l)%val(nzb_wall:nzt_wall+1, 1:ns_h_on_file_usm(l) ) )
                ALLOCATE( tmp_surf_window_h(l)%val(1:ns_h_on_file_usm(l)) )
                ALLOCATE( tmp_window_h(l)%val(nzb_wall:nzt_wall+1, 1:ns_h_on_file_usm(l) ) )
                ALLOCATE( tmp_surf_green_h(l)%val(1:ns_h_on_file_usm(l)) )
                ALLOCATE( tmp_green_h(l)%val(nzb_wall:nzt_wall+1, 1:ns_h_on_file_usm(l) ) )
                ALLOCATE( tmp_surf_mliq_h(l)%val(1:ns_h_on_file_usm(l)) )
                ALLOCATE( tmp_surf_waste_h(l)%val(1:ns_h_on_file_usm(l)) )
             ENDDO

          ENDIF

       CASE ( 'ns_v_on_file_usm')
          IF ( k == 1 )  THEN
             READ ( 13 ) ns_v_on_file_usm

             DO  l = 0, 3
                IF ( ALLOCATED( tmp_surf_wall_v(l)%val ) )    DEALLOCATE( tmp_surf_wall_v(l)%val )
                IF ( ALLOCATED( tmp_wall_v(l)%val ) )         DEALLOCATE( tmp_wall_v(l)%val )
                IF ( ALLOCATED( tmp_surf_window_v(l)%val ) )  DEALLOCATE( tmp_surf_window_v(l)%val )
                IF ( ALLOCATED( tmp_window_v(l)%val ) )       DEALLOCATE( tmp_window_v(l)%val )
                IF ( ALLOCATED( tmp_surf_green_v(l)%val ) )   DEALLOCATE( tmp_surf_green_v(l)%val )
                IF ( ALLOCATED( tmp_green_v(l)%val ) )        DEALLOCATE( tmp_green_v(l)%val )
                IF ( ALLOCATED( tmp_surf_waste_v(l)%val ) )   DEALLOCATE( tmp_surf_waste_v(l)%val )
             ENDDO

!
!--          Allocate temporary arrays for reading data on file. Note, the size of allocated surface
!--          elements do not necessarily need to match the size of present surface elements on
!--          current processor, as the number of processors between restarts can change.
             DO  l = 0, 3
                ALLOCATE( tmp_surf_wall_v(l)%val(1:ns_v_on_file_usm(l)) )
                ALLOCATE( tmp_wall_v(l)%val(nzb_wall:nzt_wall+1, 1:ns_v_on_file_usm(l) ) )
                ALLOCATE( tmp_surf_window_v(l)%val(1:ns_v_on_file_usm(l)) )
                ALLOCATE( tmp_window_v(l)%val(nzb_wall:nzt_wall+1, 1:ns_v_on_file_usm(l) ) )
                ALLOCATE( tmp_surf_green_v(l)%val(1:ns_v_on_file_usm(l)) )
                ALLOCATE( tmp_green_v(l)%val(nzb_wall:nzt_wall+1, 1:ns_v_on_file_usm(l) ) )
                ALLOCATE( tmp_surf_waste_v(l)%val(1:ns_v_on_file_usm(l)) )
             ENDDO

          ENDIF

       CASE ( 'usm_start_index_h', 'usm_start_index_v'  )
          IF ( k == 1 )  THEN

             IF ( ALLOCATED( start_index_on_file ) )  DEALLOCATE( start_index_on_file )

             ALLOCATE ( start_index_on_file(nys_on_file:nyn_on_file, nxl_on_file:nxr_on_file) )

             READ ( 13 )  start_index_on_file

          ENDIF

       CASE ( 'usm_end_index_h', 'usm_end_index_v' )
          IF ( k == 1 )  THEN

             IF ( ALLOCATED( end_index_on_file ) )  DEALLOCATE( end_index_on_file )

             ALLOCATE ( end_index_on_file(nys_on_file:nyn_on_file, nxl_on_file:nxr_on_file) )

             READ ( 13 )  end_index_on_file

          ENDIF

       CASE ( 't_surf_wall_h(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_wall_h_1(0)%val ) )                                       &
                ALLOCATE( t_surf_wall_h_1(0)%val(1:surf_usm_h(0)%ns) )
             READ ( 13 )  tmp_surf_wall_h(0)%val
          ENDIF
          CALL surface_restore_elements( t_surf_wall_h_1(0)%val, tmp_surf_wall_h(0)%val,               &
                                         surf_usm_h(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_wall_h(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_wall_h_1(1)%val ) )                                       &
                ALLOCATE( t_surf_wall_h_1(1)%val(1:surf_usm_h(1)%ns) )
             READ ( 13 )  tmp_surf_wall_h(1)%val
          ENDIF
          CALL surface_restore_elements( t_surf_wall_h_1(1)%val, tmp_surf_wall_h(1)%val,               &
                                         surf_usm_h(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_wall_v(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(0)%val ) )                                       &
                ALLOCATE( t_surf_wall_v_1(0)%val(1:surf_usm_v(0)%ns) )
             READ ( 13 )  tmp_surf_wall_v(0)%val
          ENDIF
          CALL surface_restore_elements( t_surf_wall_v_1(0)%val, tmp_surf_wall_v(0)%val,               &
                                         surf_usm_v(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_wall_v(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(1)%val ) )                                       &
                ALLOCATE( t_surf_wall_v_1(1)%val(1:surf_usm_v(1)%ns) )
             READ ( 13 )  tmp_surf_wall_v(1)%val
          ENDIF
          CALL surface_restore_elements( t_surf_wall_v_1(1)%val, tmp_surf_wall_v(1)%val,               &
                                         surf_usm_v(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_wall_v(2)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(2)%val ) )                                       &
                ALLOCATE( t_surf_wall_v_1(2)%val(1:surf_usm_v(2)%ns) )
             READ ( 13 )  tmp_surf_wall_v(2)%val
          ENDIF
          CALL surface_restore_elements( t_surf_wall_v_1(2)%val, tmp_surf_wall_v(2)%val,               &
                                         surf_usm_v(2)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_wall_v(3)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(3)%val ) )                                       &
                ALLOCATE( t_surf_wall_v_1(3)%val(1:surf_usm_v(3)%ns) )
             READ ( 13 )  tmp_surf_wall_v(3)%val
          ENDIF
          CALL surface_restore_elements( t_surf_wall_v_1(3)%val, tmp_surf_wall_v(3)%val,               &
                                         surf_usm_v(3)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_window_h(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_window_h_1(0)%val ) )                                     &
                ALLOCATE( t_surf_window_h_1(0)%val(1:surf_usm_h(0)%ns) )
             READ ( 13 )  tmp_surf_window_h(0)%val
          ENDIF
          CALL surface_restore_elements( t_surf_window_h_1(0)%val, tmp_surf_window_h(0)%val,           &
                                         surf_usm_h(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_surf_window_h(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_window_h_1(1)%val ) )                                     &
                ALLOCATE( t_surf_window_h_1(1)%val(1:surf_usm_h(1)%ns) )
             READ ( 13 )  tmp_surf_window_h(1)%val
          ENDIF
          CALL surface_restore_elements( t_surf_window_h_1(1)%val, tmp_surf_window_h(1)%val,           &
                                         surf_usm_h(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_surf_window_v(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_window_v_1(0)%val ) )                                     &
                ALLOCATE( t_surf_window_v_1(0)%val(1:surf_usm_v(0)%ns) )
             READ ( 13 )  tmp_surf_window_v(0)%val
          ENDIF
          CALL surface_restore_elements( t_surf_window_v_1(0)%val, tmp_surf_window_v(0)%val,           &
                                         surf_usm_v(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_surf_window_v(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_window_v_1(1)%val ) )                                   &
                ALLOCATE( t_surf_window_v_1(1)%val(1:surf_usm_v(1)%ns) )
             READ ( 13 )  tmp_surf_window_v(1)%val
          ENDIF
          CALL surface_restore_elements( t_surf_window_v_1(1)%val, tmp_surf_window_v(1)%val,       &
                                         surf_usm_v(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_surf_window_v(2)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_window_v_1(2)%val ) )                                   &
                ALLOCATE( t_surf_window_v_1(2)%val(1:surf_usm_v(2)%ns) )
             READ ( 13 )  tmp_surf_window_v(2)%val
          ENDIF
          CALL surface_restore_elements( t_surf_window_v_1(2)%val, tmp_surf_window_v(2)%val,       &
                                         surf_usm_v(2)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_surf_window_v(3)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_window_v_1(3)%val ) )                                   &
                ALLOCATE( t_surf_window_v_1(3)%val(1:surf_usm_v(3)%ns) )
             READ ( 13 )  tmp_surf_window_v(3)%val
          ENDIF
          CALL surface_restore_elements( t_surf_window_v_1(3)%val, tmp_surf_window_v(3)%val,       &
                                         surf_usm_v(3)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_surf_green_h(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_green_h_1(0)%val ) )                                      &
                ALLOCATE( t_surf_green_h_1(0)%val(1:surf_usm_h(0)%ns) )
             READ ( 13 )  tmp_surf_green_h(0)%val
          ENDIF
          CALL surface_restore_elements( t_surf_green_h_1(0)%val, tmp_surf_green_h(0)%val,             &
                                         surf_usm_h(0)%start_index, start_index_on_file,         &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_green_h(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_green_h_1(1)%val ) )                                      &
                ALLOCATE( t_surf_green_h_1(1)%val(1:surf_usm_h(1)%ns) )
             READ ( 13 )  tmp_surf_green_h(1)%val
          ENDIF
          CALL surface_restore_elements( t_surf_green_h_1(1)%val, tmp_surf_green_h(1)%val,             &
                                         surf_usm_h(1)%start_index, start_index_on_file,         &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_green_v(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_green_v_1(0)%val ) )                                      &
                ALLOCATE( t_surf_green_v_1(0)%val(1:surf_usm_v(0)%ns) )
             READ ( 13 )  tmp_surf_green_v(0)%val
          ENDIF
          CALL surface_restore_elements( t_surf_green_v_1(0)%val, tmp_surf_green_v(0)%val,             &
                                         surf_usm_v(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_surf_green_v(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_green_v_1(1)%val ) )                                      &
                ALLOCATE( t_surf_green_v_1(1)%val(1:surf_usm_v(1)%ns) )
             READ ( 13 )  tmp_surf_green_v(1)%val
          ENDIF
          CALL surface_restore_elements( t_surf_green_v_1(1)%val, tmp_surf_green_v(1)%val,             &
                                         surf_usm_v(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_surf_green_v(2)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_green_v_1(2)%val ) )                                      &
                ALLOCATE( t_surf_green_v_1(2)%val(1:surf_usm_v(2)%ns) )
             READ ( 13 )  tmp_surf_green_v(2)%val
          ENDIF
          CALL surface_restore_elements( t_surf_green_v_1(2)%val, tmp_surf_green_v(2)%val,             &
                                         surf_usm_v(2)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_surf_green_v(3)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_surf_green_v_1(3)%val ) )                                      &
                ALLOCATE( t_surf_green_v_1(3)%val(1:surf_usm_v(3)%ns) )
             READ ( 13 )  tmp_surf_green_v(3)%val
          ENDIF
          CALL surface_restore_elements( t_surf_green_v_1(3)%val, tmp_surf_green_v(3)%val,             &
                                         surf_usm_v(3)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 'm_liq_usm_h(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( m_liq_usm_h(0)%val ) )                                         &
                ALLOCATE( m_liq_usm_h(0)%val(1:surf_usm_h(0)%ns) )
             READ ( 13 )  tmp_surf_mliq_h(0)%val
          ENDIF
          CALL surface_restore_elements( m_liq_usm_h(0)%val, tmp_surf_mliq_h(0)%val,               &
                                         surf_usm_h(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 'm_liq_usm_h(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( m_liq_usm_h(1)%val ) )                                         &
                ALLOCATE( m_liq_usm_h(1)%val(1:surf_usm_h(1)%ns) )
             READ ( 13 )  tmp_surf_mliq_h(1)%val
          ENDIF
          CALL surface_restore_elements( m_liq_usm_h(1)%val, tmp_surf_mliq_h(1)%val,               &
                                         surf_usm_h(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 'waste_heat_h(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm_h(0)%waste_heat ) )                                   &
                ALLOCATE( surf_usm_h(0)%waste_heat(1:surf_usm_h(0)%ns) )
             READ ( 13 )  tmp_surf_waste_h(0)%val
          ENDIF
          CALL surface_restore_elements( surf_usm_h(0)%waste_heat, tmp_surf_waste_h(0)%val,          &
                                         surf_usm_h(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 'waste_heat_h(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm_h(1)%waste_heat ) )                                   &
                ALLOCATE( surf_usm_h(1)%waste_heat(1:surf_usm_h(1)%ns) )
             READ ( 13 )  tmp_surf_waste_h(1)%val
          ENDIF
          CALL surface_restore_elements( surf_usm_h(1)%waste_heat, tmp_surf_waste_h(1)%val,          &
                                         surf_usm_h(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )


       CASE ( 'waste_heat_v(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm_v(0)%waste_heat ) )                                   &
                ALLOCATE( surf_usm_v(0)%waste_heat(1:surf_usm_v(0)%ns) )
             READ ( 13 )  tmp_surf_waste_v(0)%val
          ENDIF
          CALL surface_restore_elements( surf_usm_v(0)%waste_heat, tmp_surf_waste_v(0)%val,          &
                                         surf_usm_v(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 'waste_heat_v(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm_v(1)%waste_heat ) )                                   &
                ALLOCATE( surf_usm_v(1)%waste_heat(1:surf_usm_v(1)%ns) )
             READ ( 13 )  tmp_surf_waste_v(1)%val
          ENDIF
          CALL surface_restore_elements( surf_usm_v(1)%waste_heat, tmp_surf_waste_v(1)%val,          &
                                         surf_usm_v(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 'waste_heat_v(2)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm_v(2)%waste_heat ) )                                   &
                ALLOCATE( surf_usm_v(2)%waste_heat(1:surf_usm_v(2)%ns) )
             READ ( 13 )  tmp_surf_waste_v(2)%val
          ENDIF
          CALL surface_restore_elements( surf_usm_v(2)%waste_heat, tmp_surf_waste_v(2)%val,          &
                                         surf_usm_v(2)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 'waste_heat_v(3)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm_v(3)%waste_heat ) )                                   &
                ALLOCATE( surf_usm_v(3)%waste_heat(1:surf_usm_v(3)%ns) )
             READ ( 13 )  tmp_surf_waste_v(3)%val
          ENDIF
          CALL surface_restore_elements( surf_usm_v(3)%waste_heat, tmp_surf_waste_v(3)%val,          &
                                         surf_usm_v(3)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_wall_h(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_wall_h_1(0)%val ) )                                            &
                ALLOCATE( t_wall_h_1(0)%val(nzb_wall:nzt_wall+1, 1:surf_usm_h(0)%ns) )
             READ ( 13 )  tmp_wall_h(0)%val
          ENDIF
          CALL surface_restore_elements( t_wall_h_1(0)%val, tmp_wall_h(0)%val,                         &
                                         surf_usm_h(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_wall_h(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_wall_h_1(1)%val ) )                                            &
                ALLOCATE( t_wall_h_1(1)%val(nzb_wall:nzt_wall+1, 1:surf_usm_h(1)%ns) )
             READ ( 13 )  tmp_wall_h(1)%val
          ENDIF
          CALL surface_restore_elements( t_wall_h_1(1)%val, tmp_wall_h(1)%val,                         &
                                         surf_usm_h(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_wall_v(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_wall_v_1(0)%val ) )                                            &
                ALLOCATE( t_wall_v_1(0)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(0)%ns) )
             READ ( 13 )  tmp_wall_v(0)%val
          ENDIF
          CALL surface_restore_elements( t_wall_v_1(0)%val, tmp_wall_v(0)%val,                         &
                                         surf_usm_v(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_wall_v(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_wall_v_1(1)%val ) )                                            &
                ALLOCATE( t_wall_v_1(1)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(1)%ns) )
             READ ( 13 )  tmp_wall_v(1)%val
          ENDIF
          CALL surface_restore_elements( t_wall_v_1(1)%val, tmp_wall_v(1)%val,                         &
                                         surf_usm_v(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_wall_v(2)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_wall_v_1(2)%val ) )                                            &
                ALLOCATE( t_wall_v_1(2)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(2)%ns) )
             READ ( 13 )  tmp_wall_v(2)%val
          ENDIF
          CALL surface_restore_elements( t_wall_v_1(2)%val, tmp_wall_v(2)%val,                         &
                                         surf_usm_v(2)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_wall_v(3)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_wall_v_1(3)%val ) )                                            &
                ALLOCATE( t_wall_v_1(3)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(3)%ns) )
             READ ( 13 )  tmp_wall_v(3)%val
          ENDIF
          CALL surface_restore_elements( t_wall_v_1(3)%val, tmp_wall_v(3)%val,                         &
                                         surf_usm_v(3)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_window_h(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_window_h_1(0)%val ) )                                               &
                ALLOCATE( t_window_h_1(0)%val(nzb_wall:nzt_wall+1, 1:surf_usm_h(0)%ns) )
             READ ( 13 )  tmp_window_h(0)%val
          ENDIF
          CALL surface_restore_elements( t_window_h_1(0)%val, tmp_window_h(0)%val,                     &
                                         surf_usm_h(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_window_h(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_window_h_1(1)%val ) )                                               &
                ALLOCATE( t_window_h_1(1)%val(nzb_wall:nzt_wall+1, 1:surf_usm_h(1)%ns) )
             READ ( 13 )  tmp_window_h(1)%val
          ENDIF
          CALL surface_restore_elements( t_window_h_1(1)%val, tmp_window_h(1)%val,                     &
                                         surf_usm_h(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file, nxr_on_file )

       CASE ( 't_window_v(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_window_v_1(0)%val ) )                                          &
                ALLOCATE( t_window_v_1(0)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(0)%ns) )
             READ ( 13 )  tmp_window_v(0)%val
          ENDIF
          CALL surface_restore_elements( t_window_v_1(0)%val, tmp_window_v(0)%val,                     &
                                         surf_usm_v(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_window_v(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_window_v_1(1)%val ) )                                          &
                ALLOCATE( t_window_v_1(1)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(1)%ns) )
             READ ( 13 )  tmp_window_v(1)%val
          ENDIF
          CALL surface_restore_elements( t_window_v_1(1)%val, tmp_window_v(1)%val,                     &
                                         surf_usm_v(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_window_v(2)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_window_v_1(2)%val ) )                                          &
                ALLOCATE( t_window_v_1(2)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(2)%ns) )
             READ ( 13 )  tmp_window_v(2)%val
          ENDIF
          CALL surface_restore_elements( t_window_v_1(2)%val, tmp_window_v(2)%val,                     &
                                         surf_usm_v(2)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_window_v(3)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_window_v_1(3)%val ) )                                          &
                ALLOCATE( t_window_v_1(3)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(3)%ns) )
             READ ( 13 )  tmp_window_v(3)%val
          ENDIF
          CALL surface_restore_elements( t_window_v_1(3)%val, tmp_window_v(3)%val,                     &
                                         surf_usm_v(3)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_green_h(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_green_h_1(0)%val ) )                                           &
                ALLOCATE( t_green_h_1(0)%val(nzb_wall:nzt_wall+1, 1:surf_usm_h(0)%ns) )
             READ ( 13 )  tmp_green_h(0)%val
          ENDIF
          CALL surface_restore_elements( t_green_h_1(0)%val, tmp_green_h(0)%val,                       &
                                         surf_usm_h(0)%start_index, start_index_on_file,              &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_green_h(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_green_h_1(1)%val ) )                                           &
                ALLOCATE( t_green_h_1(1)%val(nzb_wall:nzt_wall+1, 1:surf_usm_h(1)%ns) )
             READ ( 13 )  tmp_green_h(1)%val
          ENDIF
          CALL surface_restore_elements( t_green_h_1(1)%val, tmp_green_h(1)%val,                       &
                                         surf_usm_h(1)%start_index, start_index_on_file,              &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_green_v(0)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_green_v_1(0)%val ) )                                           &
                ALLOCATE( t_green_v_1(0)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(0)%ns) )
             READ ( 13 )  tmp_green_v(0)%val
          ENDIF
          CALL surface_restore_elements( t_green_v_1(0)%val, tmp_green_v(0)%val,                       &
                                         surf_usm_v(0)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_green_v(1)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_green_v_1(1)%val ) )                                           &
                ALLOCATE( t_green_v_1(1)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(1)%ns) )
             READ ( 13 )  tmp_green_v(1)%val
          ENDIF
          CALL surface_restore_elements( t_green_v_1(1)%val, tmp_green_v(1)%val,                       &
                                         surf_usm_v(1)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_green_v(2)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_green_v_1(2)%val ) )                                           &
                ALLOCATE( t_green_v_1(2)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(2)%ns) )
             READ ( 13 )  tmp_green_v(2)%val
          ENDIF
          CALL surface_restore_elements( t_green_v_1(2)%val, tmp_green_v(2)%val,                       &
                                         surf_usm_v(2)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_green_v(3)' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( t_green_v_1(3)%val ) )                                           &
                ALLOCATE( t_green_v_1(3)%val(nzb_wall:nzt_wall+1, 1:surf_usm_v(3)%ns) )
             READ ( 13 )  tmp_green_v(3)%val
          ENDIF
          CALL surface_restore_elements( t_green_v_1(3)%val, tmp_green_v(3)%val,                       &
                                         surf_usm_v(3)%start_index, start_index_on_file,           &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE usm_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!> Soubroutine reads t_surf and t_wall.
!>
!> This read routine is a counterpart of usm_wrd_local.
!> In usm_wrd_local, all array are unconditionally written, therefore all arrays are read here.
!> This is a preliminary version of reading usm data. The final version has to be discussed with
!> the developers.
!>
!> If it is possible to call usm_allocate_surface before reading the restart file, this reading
!> routine would become much simpler, because no checking for allocation will be necessary any more.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_rrd_local_mpi


    CHARACTER(LEN=1) ::  dum  !< dummy string to create input-variable name

    INTEGER(iwp) ::  l  !< loop index for surface types

    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr) ::  global_start

    LOGICAL ::  ldum  !< dummy variable

    DO  l = 0, 1

       WRITE( dum, '(I1)' )  l

       CALL rrd_mpi_io( 'usm_start_index_h_' //dum,  surf_usm_h(l)%start_index )
       CALL rrd_mpi_io( 'usm_end_index_h_' //dum, surf_usm_h(l)%end_index )
       CALL rrd_mpi_io( 'usm_global_start_h_' //dum, global_start )

       CALL rd_mpi_io_surface_filetypes( surf_usm_h(l)%start_index, surf_usm_h(l)%end_index, ldum, &
                                         global_start )

       IF ( MAXVAL( surf_usm_h(l)%end_index ) <= 0 )  CYCLE

       IF ( .NOT.  ALLOCATED( t_surf_wall_h_1(l)%val ) )                                             &
          ALLOCATE( t_surf_wall_h_1(l)%val(1:surf_usm_h(l)%ns) )
       CALL rrd_mpi_io_surface( 't_surf_wall_h(' // dum // ')', t_surf_wall_h_1(l)%val )

       IF ( .NOT.  ALLOCATED( t_surf_window_h_1(l)%val ) )                                           &
          ALLOCATE( t_surf_window_h_1(l)%val(1:surf_usm_h(l)%ns) )
       CALL rrd_mpi_io_surface( 't_surf_window_h(' // dum // ')', t_surf_window_h_1(l)%val )

       IF ( .NOT.  ALLOCATED( t_surf_green_h_1(l)%val ) )                                            &
          ALLOCATE( t_surf_green_h_1(l)%val(1:surf_usm_h(l)%ns) )
       CALL rrd_mpi_io_surface( 't_surf_green_h(' // dum // ')', t_surf_green_h_1(l)%val )

       IF ( .NOT.  ALLOCATED( m_liq_usm_h_1(l)%val ) )                                            &
          ALLOCATE( m_liq_usm_h_1(l)%val(1:surf_usm_h(l)%ns) )
       CALL rrd_mpi_io_surface( 'm_liq_usm_h(' // dum // ')', m_liq_usm_h_1(l)%val )

       IF ( indoor_model )  THEN
          IF ( .NOT.  ALLOCATED( surf_usm_h(l)%waste_heat ) )                                     &
             ALLOCATE( surf_usm_h(l)%waste_heat(1:surf_usm_h(l)%ns) )
          CALL rrd_mpi_io_surface( 'waste_heat_h(' // dum // ')', surf_usm_h(l)%waste_heat )
       ENDIF

    ENDDO
    DO  l = 0, 3

       WRITE( dum, '(I1)' )  l

       CALL rrd_mpi_io( 'usm_start_index_v_' //dum, surf_usm_v(l)%start_index )
       CALL rrd_mpi_io( 'usm_end_index_v_' // dum, surf_usm_v(l)%end_index )
       CALL rrd_mpi_io( 'usm_global_start_v_' // dum, global_start )

       CALL rd_mpi_io_surface_filetypes( surf_usm_v(l)%start_index, surf_usm_v(l)%end_index, ldum, &
                                         global_start )

       IF ( MAXVAL( surf_usm_v(l)%end_index ) <= 0 )  CYCLE

       IF ( .NOT.  ALLOCATED( t_surf_wall_v_1(l)%val ) )                                             &
          ALLOCATE( t_surf_wall_v_1(l)%val(1:surf_usm_v(l)%ns) )
       CALL rrd_mpi_io_surface( 't_surf_wall_v(' // dum // ')', t_surf_wall_v_1(l)%val )

       IF ( .NOT.  ALLOCATED( t_surf_window_v_1(l)%val ) )                                           &
          ALLOCATE( t_surf_window_v_1(l)%val(1:surf_usm_v(l)%ns) )
       CALL rrd_mpi_io_surface( 't_surf_window_v(' // dum // ')', t_surf_window_v_1(l)%val )

       IF ( .NOT.  ALLOCATED( t_surf_green_v_1(l)%val ) )                                            &
          ALLOCATE( t_surf_green_v_1(l)%val(1:surf_usm_v(l)%ns) )
       CALL rrd_mpi_io_surface( 't_surf_green_v(' // dum // ')', t_surf_green_v_1(l)%val)

    ENDDO

    DO  l = 0, 1

       WRITE( dum, '(I1)' )  l

       CALL rrd_mpi_io( 'usm_start_index_h_2_' //dum,  surf_usm_h(l)%start_index )
       CALL rrd_mpi_io( 'usm_end_index_h_2_' //dum, surf_usm_h(l)%end_index )
       CALL rrd_mpi_io( 'usm_global_start_h_2_' //dum, global_start )

       CALL rd_mpi_io_surface_filetypes( surf_usm_h(l)%start_index, surf_usm_h(l)%end_index, ldum,          &
                                         global_start )

       IF ( MAXVAL( surf_usm_h(l)%end_index ) <= 0 )  CYCLE

       IF ( .NOT.  ALLOCATED( t_wall_h_1(l)%val ) )                                                          &
          ALLOCATE( t_wall_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       CALL rrd_mpi_io_surface( 't_wall_h(' // dum // ')', t_wall_h_1(l)%val )

       IF ( .NOT.  ALLOCATED( t_window_h_1(l)%val ) )                                                        &
          ALLOCATE( t_window_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       CALL rrd_mpi_io_surface( 't_window_h(' // dum // ')', t_window_h_1(l)%val )

       IF ( .NOT.  ALLOCATED( t_green_h_1(l)%val ) )                                                         &
          ALLOCATE( t_green_h_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_h(l)%ns) )
       CALL rrd_mpi_io_surface( 't_green_h(' // dum // ')', t_green_h_1(l)%val )

    ENDDO

    DO  l = 0, 3

       WRITE( dum, '(I1)' )  l

       CALL rrd_mpi_io( 'usm_start_index_v_2_' //dum, surf_usm_v(l)%start_index )
       CALL rrd_mpi_io( 'usm_end_index_v_2_' // dum, surf_usm_v(l)%end_index )
       CALL rrd_mpi_io( 'usm_global_start_v_2_' // dum, global_start )

       CALL rd_mpi_io_surface_filetypes( surf_usm_v(l)%start_index, surf_usm_v(l)%end_index, ldum, &
                                         global_start )

       IF ( MAXVAL( surf_usm_v(l)%end_index ) <= 0 )  CYCLE

       IF ( .NOT. ALLOCATED( t_wall_v_1(l)%val ) )                                                   &
          ALLOCATE ( t_wall_v_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       CALL rrd_mpi_io_surface( 't_wall_v(' // dum // ')', t_wall_v_1(l)%val )

       IF ( .NOT. ALLOCATED( t_window_v_1(l)%val ) )                                                 &
          ALLOCATE ( t_window_v_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       CALL rrd_mpi_io_surface( 't_window_v(' // dum // ')', t_window_v_1(l)%val )

       IF ( .NOT. ALLOCATED( t_green_v_1(l)%val ) )                                                  &
          ALLOCATE ( t_green_v_1(l)%val(nzb_wall:nzt_wall+1,1:surf_usm_v(l)%ns) )
       CALL rrd_mpi_io_surface( 't_green_v(' // dum // ')', t_green_v_1(l)%val )

    ENDDO

 END SUBROUTINE usm_rrd_local_mpi

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the ground/roof/wall surface. It follows the basic ideas and
!> structure of lsm_energy_balance with many simplifications and adjustments.
!> TODO better description
!> No calculation of window surface temperatures during spinup to increase maximum possible timstep
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_energy_balance( during_spinup )

    LOGICAL       ::  during_spinup  !< flag indicating soil/wall spinup phase
    INTEGER(iwp)  ::  l              !< loop index for surface types
!
!-- Call for horizontal surfaces
    DO l = 0, 1
       CALL usm_surface_energy_balance( .TRUE., l, during_spinup )
       CALL usm_green_heat_model( .TRUE., l )
       CALL usm_wall_heat_model( .TRUE., l, during_spinup )
    ENDDO
!
!--   Call for vertical surfaces
    DO l = 0, 3
       CALL usm_surface_energy_balance( .FALSE., l, during_spinup )
       CALL usm_green_heat_model( .FALSE., l )
       CALL usm_wall_heat_model( .FALSE., l, during_spinup )
    ENDDO

 END SUBROUTINE usm_energy_balance

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Solver for the energy balance at the ground/roof/wall surface. It follows the basic ideas and
!> structure of lsm_energy_balance with many simplifications and adjustments.
!> TODO better description
!> No calculation of window surface temperatures during spinup to increase maximum possible timstep
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_surface_energy_balance( horizontal, l, during_spinup )

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz


    IMPLICIT NONE

    LOGICAL       ::  horizontal     !< Flag indicating horizontal or vertical surfaces
    INTEGER(iwp)  ::  l              !< direction index
    LOGICAL       ::  during_spinup  !< flag indicating soil/wall spinup phase

    INTEGER(iwp)  ::  i, j, k, m     !< running indices
    INTEGER(iwp)  ::  i_off          !< offset to determine index of surface element, seen from atmospheric grid point, for x
    INTEGER(iwp)  ::  j_off          !< offset to determine index of surface element, seen from atmospheric grid point, for y
    INTEGER(iwp)  ::  k_off          !< offset to determine index of surface element, seen from atmospheric grid point, for z


    REAL(wp)  ::  coef_1             !< first coeficient for prognostic equation
    REAL(wp)  ::  coef_window_1      !< first coeficient for prognostic window equation
    REAL(wp)  ::  coef_green_1       !< first coeficient for prognostic green wall equation
    REAL(wp)  ::  coef_2             !< second  coeficient for prognostic equation
    REAL(wp)  ::  coef_window_2      !< second  coeficient for prognostic window equation
    REAL(wp)  ::  coef_green_2       !< second  coeficient for prognostic green wall equation
    REAL(wp)  ::  frac_win           !< window fraction, used to restore original values during spinup
    REAL(wp)  ::  frac_green         !< green fraction, used to restore original values during spinup
    REAL(wp)  ::  frac_wall          !< wall fraction, used to restore original values during spinup
    REAL(wp)  ::  f_shf              !< factor for shf_eb
    REAL(wp)  ::  f_shf_window       !< factor for shf_eb window
    REAL(wp)  ::  f_shf_green        !< factor for shf_eb green wall
    REAL(wp)  ::  lambda_surface     !< current value of lambda_surface (heat conductivity
                                     !<between air and wall)
    REAL(wp)  ::  lambda_surface_window  !< current value of lambda_surface (heat conductivity
                                         !< between air and window)
    REAL(wp)  ::  lambda_surface_green   !< current value of lambda_surface (heat conductivity
                                                                    !< between air and greeb wall)
    REAL(wp)  ::  rho_cp             !< rho_wall_surface * c_p
    REAL(wp)  ::  stend_wall         !< surface tendency
    REAL(wp)  ::  stend_window       !< surface tendency
    REAL(wp)  ::  stend_green        !< surface tendency


    REAL(wp)  ::  dq_s_dt,                 &  !< derivate of q_s with respect to T
                  drho_l_lv,               &  !< frequently used parameter for green layers
                  e,                       &  !< water vapour pressure
                  e_s,                     &  !< water vapour saturation pressure
                  e_s_dt,                  &  !< derivate of e_s with respect to T
                  f_qsws,                  &  !< factor for qsws
                  f_qsws_veg,              &  !< factor for qsws_veg
                  f_qsws_liq,              &  !< factor for qsws_liq
                  f1,                      &  !< resistance correction term 1
                  f2,                      &  !< resistance correction term 2
                  f3,                      &  !< resistance correction term 3
                  m_max_depth = 0.0002_wp, &  !< Maximum capacity of the water reservoir (m)
                  m_liq_max,               &  !< maxmimum value of the liq. water reservoir
                  qv1,                     &  !< specific humidity at first grid level
                  q_s,                     &  !< saturation specific humidity
                  rho_lv,                  &  !< frequently used parameter for green layers
                  tend,                    &  !< tendency
                  ueff                        !< limited near-surface wind speed - used for calculation of resistance

    TYPE(surf_type),        POINTER ::  surf              !< surface-date type variable
    TYPE(surf_type_1d_usm), POINTER ::  t_surf_green      !<
    TYPE(surf_type_1d_usm), POINTER ::  t_surf_green_p    !<
    TYPE(surf_type_1d_usm), POINTER ::  t_surf_wall       !<
    TYPE(surf_type_1d_usm), POINTER ::  t_surf_wall_p     !<
    TYPE(surf_type_1d_usm), POINTER ::  t_surf_window     !<
    TYPE(surf_type_1d_usm), POINTER ::  t_surf_window_p   !<
    TYPE(surf_type_2d_usm), POINTER ::  t_green           !<
    TYPE(surf_type_2d_usm), POINTER ::  t_wall            !<
    TYPE(surf_type_2d_usm), POINTER ::  t_window          !<

    LOGICAL   ::  upward                      !< Flag indicating upward horizontal surfaces

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'usm_surface_energy_balance:', horizontal, l, during_spinup
       CALL debug_message( debug_string, 'start' )
    ENDIF

    upward = .FALSE.
    IF ( horizontal )  THEN
       surf              => surf_usm_h(l)
       t_surf_wall       => t_surf_wall_h(l)
       t_surf_wall_p     => t_surf_wall_h_p(l)
       t_surf_window     => t_surf_window_h(l)
       t_surf_window_p   => t_surf_window_h_p(l)
       t_surf_green      => t_surf_green_h(l)
       t_surf_green_p    => t_surf_green_h_p(l)
       t_wall            => t_wall_h(l)
       t_window          => t_window_h(l)
       t_green           => t_green_h(l)
       IF ( l == 0 ) upward = .TRUE.
    ELSE
       surf              => surf_usm_v(l)
       t_surf_wall       => t_surf_wall_v(l)
       t_surf_wall_p     => t_surf_wall_v_p(l)
       t_surf_window     => t_surf_window_v(l)
       t_surf_window_p   => t_surf_window_v_p(l)
       t_surf_green      => t_surf_green_v(l)
       t_surf_green_p    => t_surf_green_v_p(l)
       t_wall            => t_wall_v(l)
       t_window          => t_window_v(l)
       t_green           => t_green_v(l)
    ENDIF
!
!-- Index offset of surface element point with respect to adjoining atmospheric grid point
    k_off = surf%koff
    j_off = surf%joff
    i_off = surf%ioff

!
!-- First, treat horizontal surface elements
    !$OMP PARALLEL DO PRIVATE (m, i, j, k, frac_win, frac_wall, frac_green, lambda_surface,       &
    !$OMP&             lambda_surface_window, lambda_surface_green, ueff, qv1, rho_cp, rho_lv,    &
    !$OMP&             drho_l_lv, f_shf, f_shf_window, f_shf_green, m_total, f1, f2, e_s, e,      &
    !$OMP&             f3, f_qsws_veg, q_s, f_qsws_liq, f_qsws, e_s_dt, dq_s_dt, coef_1,          &
    !$OMP&             coef_window_1, coef_green_1, coef_2, coef_window_2, coef_green_2,          &
    !$OMP&             stend_wall, stend_window, stend_green, tend, m_liq_max) SCHEDULE (STATIC)
    DO  m = 1, surf%ns
!
!--    During spinup set green and window fraction to zero and restore at the end of the loop.
!--    Note, this is a temporary fix and needs to be removed later.
       IF ( during_spinup )  THEN
          frac_win   = 0.0_wp
          frac_wall  = 1.0_wp
          frac_green = 0.0_wp
       ELSE
          frac_win   = surf%frac(m,ind_wat_win)
          frac_wall  = surf%frac(m,ind_veg_wall)
          frac_green = surf%frac(m,ind_pav_green)
       ENDIF
!
!--    Get indices of respective grid point
       i = surf%i(m)
       j = surf%j(m)
       k = surf%k(m)
!
!--    TODO - how to calculate lambda_surface for horizontal surfaces
!--    (lambda_surface shoud be set according to stratification in land surface model)
       lambda_surface        = surf%lambda_surf(m)
       lambda_surface_window = surf%lambda_surf_window(m)
       lambda_surface_green  = surf%lambda_surf_green(m)

       IF ( humidity )  THEN
          qv1 = q(k,j,i)
       ELSE
          qv1 = 0.0_wp
       ENDIF
!
!--    Calculate rho * c_p coefficient at surface layer
       rho_cp  = c_p * hyp(k) / ( r_d * surf%pt1(m) * exner(k) )

       IF ( frac_green > 0.0_wp )  THEN
!
!--       Calculate frequently used parameters
          rho_lv    = rho_cp / c_p * l_v
          drho_l_lv = 1.0_wp / ( rho_l * l_v )
       ENDIF

!
!--    Calculate aerodyamic resistance.
       IF ( upward ) THEN
!--       Calculation for horizontal upward facing surfaces follows LSM formulation
!--       pt, us, ts are not available for the prognostic time step, data from the
!--       last time step is used here.
!--       Workaround: use single r_a as stability is only treated for the average temperature
          surf%r_a(m)        = ( surf%pt1(m) - surf%pt_surface(m) ) /               &
                                     ( surf%ts(m) * surf%us(m) + 1.0E-20_wp )
       ELSE
!--           Calculation of r_a for vertical and downward facing horizontal surfaces
!--
!--           Heat transfer coefficient for forced convection along vertical walls follows formulation
!--           in TUF3d model (Krayenhoff & Voogt, 2006)
!--
!--           H = httc (Tsfc - Tair)
!--           httc = rw * (11.8 + 4.2 * Ueff) - 4.0
!--
!--                rw: Wall patch roughness relative to 1.0 for concrete
!--                Ueff: Effective wind speed
!--                - 4.0 is a reduction of Rowley et al (1930) formulation based on
!--                Cole and Sturrock (1977)
!--
!--                Ucan: Canyon wind speed
!--                wstar: Convective velocity
!--                Qs: Surface heat flux
!--                zH: Height of the convective layer
!--                wstar = (g/Tcan*Qs*zH)**(1./3.)
!--           Effective velocity components must always be defined at scalar grid point. The wall
!--           normal component is obtained by simple linear interpolation. (An alternative would be an
!--           logarithmic interpolation.) Parameter roughness_concrete (default value = 0.001) is used
!--           to calculation of roughness relative to concrete. Note, wind velocity is limited
!--           to avoid division by zero. The nominator can become <= 0.0 for values z0 < 3*10E-4.
              ueff        = MAX ( SQRT( ( ( u(k,j,i) + u(k,j,i+1) ) * 0.5_wp )**2 +               &
                                        ( ( v(k,j,i) + v(k,j+1,i) ) * 0.5_wp )**2 +               &
                                        ( ( w(k,j,i) + w(k-1,j,i) ) * 0.5_wp )**2 ),              &
                                  ( ( 4.0_wp + 0.1_wp )                                           &
                                    / ( surf%z0(m) * d_roughness_concrete )                       &
                                   - 11.8_wp ) / 4.2_wp                                           &
                                 )

              surf%r_a(m) = rho_cp / ( surf%z0(m) * d_roughness_concrete                          &
                                     * ( 11.8_wp + 4.2_wp * ueff )  - 4.0_wp  )
       ENDIF

!--    Make sure that the resistance does not drop to zero
!--    end does not exceed a maxmium value in case of zero velocities
       IF ( surf%r_a(m)        < 1.0_wp )  surf%r_a(m)        = 1.0_wp
       IF ( surf%r_a(m)        > 300.0_wp )  surf%r_a(m)        = 300.0_wp
!
!--    Aeorodynamical resistance for the window and green fractions are set to the same value
       surf%r_a_window(m) = surf%r_a(m)
       surf%r_a_green(m)  = surf%r_a(m)
!
!--    Factor for shf_eb
       f_shf         = rho_cp / surf%r_a(m)
       f_shf_window  = rho_cp / surf%r_a_window(m)
       f_shf_green   = rho_cp / surf%r_a_green(m)

       IF ( frac_green > 0.0_wp )  THEN
!
!--       Adapted from LSM:
!--       Second step: calculate canopy resistance r_canopy. f1-f3 here are defined as 1/f1-f3
!--       as in ECMWF documentation

!--       f1: Correction for incoming shortwave radiation (stomata close at night)
          f1 = MIN( 1.0_wp, ( 0.004_wp * surf%rad_sw_in(m) + 0.05_wp ) /                          &
                    (0.81_wp * ( 0.004_wp * surf%rad_sw_in(m) + 1.0_wp ) ) )
!
!--       f2: Correction for soil moisture availability to plants (the integrated soil moisture must
!--       thus be considered here) f2 = 0 for very dry soils
          IF ( upward ) THEN
             m_total = 0.0_wp
             DO  k = nzb_wall, nzt_wall+1
                 m_total = m_total + rootfr_h(l)%val(nzb_wall,m) &
                                     * MAX( swc_h(l)%val(nzb_wall,m),wilt_h(l)%val(nzb_wall,m) )
             ENDDO

             IF ( m_total > wilt_h(l)%val(nzb_wall,m)  .AND.  m_total < fc_h(l)%val(nzb_wall,m) )  THEN
                f2 = ( m_total - wilt_h(l)%val(nzb_wall,m) ) / (fc_h(l)%val(nzb_wall,m) - wilt_h(l)%val(nzb_wall,m) )
             ELSEIF ( m_total >= fc_h(l)%val(nzb_wall,m) )  THEN
                f2 = 1.0_wp
             ELSE
                f2 = 1.0E-20_wp
             ENDIF
          ELSE
             f2=1.0_wp
          ENDIF
!
!--       Calculate water vapour pressure at saturation
          e_s = 0.01_wp * magnus_tl( t_surf_green%val(m) )
!
!--       f3: Correction for vapour pressure deficit
          IF ( surf%g_d(m) /= 0.0_wp )  THEN
!--          Calculate vapour pressure
             e  = qv1 * surface_pressure / ( qv1 + 0.622_wp )
             f3 = EXP ( - surf%g_d(m) * (e_s - e) )
          ELSE
             f3 = 1.0_wp
          ENDIF
!
!--       Calculate canopy resistance. In case that c_veg is 0 (bare soils), this calculation is
!--       obsolete, as r_canopy is not used below.
!--       To do: check for very dry soil -> r_canopy goes to infinity
          surf%r_canopy(m) = surf%r_canopy_min(m) /                                               &
                                  ( surf%lai(m) * f1 * f2 * f3 + 1.0E-20_wp )
!
!--       Calculate saturation specific humidity
          q_s = 0.622_wp * e_s / ( surface_pressure - e_s )
!
!--       In case of dewfall, set evapotranspiration to zero
!--       All super-saturated water is then removed from the air
          IF ( humidity  .AND.  q_s <= qv1 )  THEN
             surf%r_canopy(m) = 0.0_wp
          ENDIF

          IF ( upward ) THEN
!--          Calculate the maximum possible liquid water amount on plants and bare surface. For
!--          vegetated surfaces, a maximum depth of 0.2 mm is assumed, while paved surfaces might hold
!--          up 1 mm of water. The liquid water fraction for paved surfaces is calculated after
!--          Noilhan & Planton (1989), while the ECMWF formulation is used for vegetated surfaces and
!--          bare soils.
             m_liq_max = m_max_depth * ( surf%lai(m) )
             surf%c_liq(m) = MIN( 1.0_wp, ( m_liq_usm_h(l)%val(m) / m_liq_max )**0.67 )

!
!--          Calculate coefficients for the total evapotranspiration
!--          In case of water surface, set vegetation and soil fluxes to zero.
!--          For pavements, only evaporation of liquid water is possible.
             f_qsws_veg  = rho_lv * ( 1.0_wp - surf%c_liq(m) ) /                                  &
                           ( surf%r_a_green(m) + surf%r_canopy(m) )
             f_qsws_liq  = rho_lv * surf%c_liq(m) / surf%r_a_green(m)
             f_qsws = f_qsws_veg + f_qsws_liq
          ELSE
             f_qsws_veg  = rho_lv * ( 1.0_wp - 0.0_wp ) / & !surf%c_liq(m) ) /                    &
                                ( surf%r_a_green(m) + surf%r_canopy(m) )
             f_qsws_liq  = 0._wp ! rho_lv * surf%c_liq(m) / surf%r_a_green(m)
             f_qsws = f_qsws_veg + f_qsws_liq
          ENDIF
!
!--       Calculate derivative of q_s for Taylor series expansion
          e_s_dt = e_s * ( 17.269_wp / ( t_surf_green%val(m) - 35.86_wp ) - 17.269_wp             &
                   * ( t_surf_green%val(m) - degc_to_k )                                          &
                   / ( t_surf_green%val(m) - 35.86_wp )**2 )
          dq_s_dt = 0.622_wp * e_s_dt / ( surface_pressure - e_s_dt )
       ENDIF
!
!--    Add LW up so that it can be removed in prognostic equation
       surf%rad_net_l(m) = surf%rad_sw_in(m)  - surf%rad_sw_out(m) +                              &
                           surf%rad_lw_in(m)  - surf%rad_lw_out(m)
!
!--    Numerator of the prognostic equation
!--    Todo: Adjust to tile approach. So far, emissivity for wall (element 0) is used
!--    Rem: Coef +1 corresponds to -lwout included in calculation of radnet_l
       coef_1 = surf%rad_net_l(m) + ( 3.0_wp + 1.0_wp )                                           &
               * surf%emissivity(m,ind_veg_wall) * sigma_sb * t_surf_wall%val(m)**4               &
               + f_shf * surf%pt1(m) + lambda_surface * t_wall%val(nzb_wall,m)

       IF ( ( .NOT. during_spinup ) .AND. (frac_win > 0.0_wp ) )  THEN
          coef_window_1 = surf%rad_net_l(m) +  ( 3.0_wp + 1.0_wp )                                &
                          * surf%emissivity(m,ind_wat_win) * sigma_sb                             &
                          * t_surf_window%val(m)**4 + f_shf_window * surf%pt1(m)                  &
                          + lambda_surface_window * t_window%val(nzb_wall,m)
       ENDIF
       IF ( ( humidity ) .AND. ( frac_green > 0.0_wp ) )  THEN
          coef_green_1 = surf%rad_net_l(m) + ( 3.0_wp + 1.0_wp )                                  &
                         * surf%emissivity(m,ind_pav_green) * sigma_sb                            &
                         * t_surf_green%val(m)**4 + f_shf_green * surf%pt1(m)                     &
                         + f_qsws * ( qv1 - q_s + dq_s_dt * t_surf_green%val(m) )                 &
                         + lambda_surface_green * t_green%val(nzb_wall,m)
       ELSE
          coef_green_1 = surf%rad_net_l(m) + ( 3.0_wp + 1.0_wp )                                  &
                         * surf%emissivity(m,ind_pav_green) * sigma_sb * t_surf_green%val(m)**4   &
                         + f_shf_green * surf%pt1(m) + lambda_surface_green                       &
                         * t_green%val(nzb_wall,m)
       ENDIF
!
!--    Denominator of the prognostic equation
       coef_2 = 4.0_wp * surf%emissivity(m,ind_veg_wall) * sigma_sb * t_surf_wall%val(m)**3       &
                + lambda_surface + f_shf / exner(k)
       IF ( ( .NOT. during_spinup ) .AND. ( frac_win > 0.0_wp ) )  THEN
          coef_window_2 = 4.0_wp * surf%emissivity(m,ind_wat_win) * sigma_sb *                    &
                          t_surf_window%val(m)**3 + lambda_surface_window + f_shf_window / exner(k)
       ENDIF
       IF ( ( humidity ) .AND. ( frac_green > 0.0_wp ) )  THEN
          coef_green_2 = 4.0_wp * surf%emissivity(m,ind_pav_green) * sigma_sb *                   &
                         t_surf_green%val(m)**3 + f_qsws * dq_s_dt + lambda_surface_green         &
                         + f_shf_green / exner(k)
       ELSE
          coef_green_2 = 4.0_wp * surf%emissivity(m,ind_pav_green) * sigma_sb                     &
                         * t_surf_green%val(m)**3 + lambda_surface_green + f_shf_green / exner(k)
       ENDIF
!
!--    Implicit solution when the surface layer has no heat capacity, otherwise use RK3 scheme.
       t_surf_wall_p%val(m) = ( coef_1 * dt_3d * tsc(2) + surf%c_surface(m)                       &
                            * t_surf_wall%val(m) )                                                &
                            / ( surf%c_surface(m) + coef_2 * dt_3d * tsc(2) )
       IF ( ( .NOT. during_spinup ) .AND. (frac_win > 0.0_wp) )  THEN
          t_surf_window_p%val(m) = ( coef_window_1 * dt_3d * tsc(2) +                             &
                                   surf%c_surface_window(m) * t_surf_window%val(m) ) /            &
                                ( surf%c_surface_window(m) + coef_window_2 * dt_3d * tsc(2) )
       ENDIF
       t_surf_green_p%val(m) = ( coef_green_1 * dt_3d * tsc(2) +                                  &
                               surf%c_surface_green(m) * t_surf_green%val(m) )                    &
                             / ( surf%c_surface_green(m) + coef_green_2 * dt_3d * tsc(2) )
!
!--    Add RK3 term
       t_surf_wall_p%val(m) = t_surf_wall_p%val(m) + dt_3d * tsc(3) *                             &
                            surf%tt_surface_wall_m(m)
       t_surf_window_p%val(m) = t_surf_window_p%val(m) + dt_3d * tsc(3) *                         &
                              surf%tt_surface_window_m(m)
       t_surf_green_p%val(m) = t_surf_green_p%val(m) + dt_3d * tsc(3) *                           &
                             surf%tt_surface_green_m(m)
!
!--    Store surface temperature on pt_surface. Further, in case humidity is used, store also
!--    vpt_surface, which is, due to the lack of moisture on roofs, simply assumed to be the surface
!--    temperature.
       surf%pt_surface(m) = ( frac_wall * t_surf_wall_p%val(m)                                    &
                                 + frac_win * t_surf_window_p%val(m)                              &
                                 + frac_green * t_surf_green_p%val(m)                             &
                             ) / exner(k)

       IF ( humidity )  surf%vpt_surface(m) = surf%pt_surface(m)
!
!--    Calculate true tendency
       stend_wall = ( t_surf_wall_p%val(m) - t_surf_wall%val(m) - dt_3d * tsc(3) *                &
                      surf%tt_surface_wall_m(m) ) / ( dt_3d  * tsc(2) )
       stend_window = ( t_surf_window_p%val(m) - t_surf_window%val(m) - dt_3d * tsc(3) *          &
                        surf%tt_surface_window_m(m) ) / ( dt_3d  * tsc(2) )
       stend_green = ( t_surf_green_p%val(m) - t_surf_green%val(m) - dt_3d * tsc(3) *             &
                       surf%tt_surface_green_m(m) ) / ( dt_3d  * tsc(2) )
!
!--    Calculate t_surf tendencies for the next Runge-Kutta step
       IF ( timestep_scheme(1:5) == 'runge' )  THEN
          IF ( intermediate_timestep_count == 1 )  THEN
             surf%tt_surface_wall_m(m)   = stend_wall
             surf%tt_surface_window_m(m) = stend_window
             surf%tt_surface_green_m(m)  = stend_green
          ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
             surf%tt_surface_wall_m(m) = -9.5625_wp * stend_wall +                                &
                                               5.3125_wp * surf%tt_surface_wall_m(m)
             surf%tt_surface_window_m(m) = -9.5625_wp * stend_window +                            &
                                                 5.3125_wp * surf%tt_surface_window_m(m)
             surf%tt_surface_green_m(m) = -9.5625_wp * stend_green +                              &
                                                5.3125_wp * surf%tt_surface_green_m(m)
          ENDIF
       ENDIF
!
!--    In case of fast changes in the skin temperature, it is required to update the radiative
!--    fluxes in order to keep the solution stable
       IF ( ( ( ABS( t_surf_wall_p%val(m)   - t_surf_wall%val(m) )   > 1.0_wp )   .OR.            &
            (   ABS( t_surf_green_p%val(m)  - t_surf_green%val(m) )  > 1.0_wp )   .OR.            &
            (   ABS( t_surf_window_p%val(m) - t_surf_window%val(m) ) > 1.0_wp ) )                 &
               .AND.  unscheduled_radiation_calls  )  THEN
          force_radiation_call_l = .TRUE.
       ENDIF
!
!--    Calculate new fluxes
!--    Rad_net_l is never used!
       surf%rad_net_l(m) = surf%rad_net_l(m) + frac_wall                                          &
                                 * sigma_sb * surf%emissivity(m,ind_veg_wall)                     &
                                 * ( t_surf_wall_p%val(m)**4 - t_surf_wall%val(m)**4 )            &
                                 + frac_win * sigma_sb                                            &
                                 * surf%emissivity(m,ind_wat_win)                                 &
                                 * ( t_surf_window_p%val(m)**4 - t_surf_window%val(m)**4 )        &
                                 + frac_green * sigma_sb                                          &
                                 * surf%emissivity(m,ind_pav_green)                               &
                                 * ( t_surf_green_p%val(m)**4 - t_surf_green%val(m)**4 )

       surf%wghf_eb(m)   = lambda_surface * ( t_surf_wall_p%val(m) - t_wall%val(nzb_wall,m) )
       surf%wghf_eb_green(m)  = lambda_surface_green                                              &
                                      * ( t_surf_green_p%val(m) - t_green%val(nzb_wall,m) )
       surf%wghf_eb_window(m) = lambda_surface_window                                             &
                                      * ( t_surf_window_p%val(m) - t_window%val(nzb_wall,m) )

!
!--    Ground/wall/roof surface heat flux
       surf%wshf_eb(m)   = - f_shf  * ( surf%pt1(m) - t_surf_wall_p%val(m) / exner(k) )           &
                                 * frac_wall - f_shf_window                                       &
                                 * ( surf%pt1(m) - t_surf_window_p%val(m) / exner(k) )            &
                                 * frac_win - f_shf_green                                         &
                                 * ( surf%pt1(m) - t_surf_green_p%val(m) / exner(k) )             &
                                 * frac_green
!
!--    Store kinematic surface heat fluxes for utilization in other processes diffusion_s,
!--    surface_layer_fluxes,...
       surf%shf(m) = surf%wshf_eb(m) / c_p
!
!--    If the indoor model is applied, further add waste heat from buildings to the kinematic flux.
       IF ( indoor_model )  THEN
          surf%shf(m) = surf%shf(m) + surf%waste_heat(m) / c_p
       ENDIF

       IF ( humidity .AND.  frac_green > 0.0_wp )  THEN!
!--       Calculate true surface resistance
          IF ( upward ) THEN
             surf%qsws(m)  = - f_qsws * ( qv1 - q_s + dq_s_dt * t_surf_green%val(m)               &
                                                - dq_s_dt * t_surf_green_p%val(m) )
             surf%qsws(m) = surf%qsws(m) / l_v
             surf%qsws_veg(m)  = - f_qsws_veg  * ( qv1 - q_s + dq_s_dt * t_surf_green%val(m)      &
                                                         - dq_s_dt * t_surf_green_p%val(m) )
             surf%qsws_liq(m)  = - f_qsws_liq  * ( qv1 - q_s + dq_s_dt * t_surf_green%val(m)      &
                                                         - dq_s_dt * t_surf_green_p%val(m) )
             surf%r_s(m) = - rho_lv * ( qv1 - q_s + dq_s_dt * t_surf_green%val(m)                 &
                                 - dq_s_dt * t_surf_green_p%val(m) ) /                            &
                                 (surf%qsws(m) + 1.0E-20) - surf%r_a_green(m)
             IF ( precipitation )  THEN
!--             Calculate change in liquid water reservoir due to dew fall or evaporation of liquid water
!--             If precipitation is activated, add rain water to qsws_liq and qsws_soil according the
!--             the vegetation coverage. Precipitation_rate is given in mm.
!--             Add precipitation to liquid water reservoir, if possible. Otherwise, add the water
!--             to soil. In case of pavements, the exceeding water amount is implicitely removed as
!--             runoff as qsws_soil is then not used in the soil model
                IF ( m_liq_usm_h(l)%val(m) /= m_liq_max )  THEN
                   surf%qsws_liq(m) = surf%qsws_liq(m)                                            &
                                            + frac_green                                          &
                                            * prr(k+k_off,j+j_off,i+i_off) * hyrho(k+k_off)       &
                                            * 0.001_wp * rho_l * l_v
                ENDIF
             ENDIF
!
!--          If the air is saturated, check the reservoir water level
             IF ( surf%qsws(m) < 0.0_wp )  THEN
!
!--                Check if reservoir is full (avoid values > m_liq_max) In that case, qsws_liq goes to
!--                qsws_soil. In this case qsws_veg is zero anyway (because c_liq = 1), so that tend is
!--                zero and no further check is needed
                IF ( m_liq_usm_h(l)%val(m) == m_liq_max )  THEN
                   surf%qsws_liq(m)  = 0.0_wp
                ENDIF
!--                In case qsws_veg becomes negative (unphysical behavior), let the water enter the
!--                liquid water reservoir as dew on the plant
                IF ( surf%qsws_veg(m) < 0.0_wp )  THEN
                   surf%qsws_liq(m) = surf%qsws_liq(m) + surf%qsws_veg(m)
                   surf%qsws_veg(m) = 0.0_wp
                ENDIF
             ENDIF

             tend = - surf%qsws_liq(m) * drho_l_lv
             m_liq_usm_h_p(l)%val(m) = m_liq_usm_h(l)%val(m) + dt_3d *                            &
                                       ( tsc(2) * tend + tsc(3) * tm_liq_usm_h_m(l)%val(m) )
!
!--             Check if reservoir is overfull -> reduce to maximum
!--             (conservation of water is violated here)
             m_liq_usm_h_p(l)%val(m) = MIN( m_liq_usm_h_p(l)%val(m), m_liq_max )
!
!--             Check if reservoir is empty (avoid values < 0.0) (conservation of water is violated here)
             m_liq_usm_h_p(l)%val(m) = MAX( m_liq_usm_h_p(l)%val(m), 0.0_wp )
!
!--             Calculate m_liq tendencies for the next Runge-Kutta step
             IF ( timestep_scheme(1:5) == 'runge' )  THEN
                IF ( intermediate_timestep_count == 1 )  THEN
                   tm_liq_usm_h_m(l)%val(m) = tend
                ELSEIF ( intermediate_timestep_count < intermediate_timestep_count_max )  THEN
                   tm_liq_usm_h_m(l)%val(m) = -9.5625_wp * tend +                                 &
                                                  5.3125_wp * tm_liq_usm_h_m(l)%val(m)
                ENDIF
             ENDIF
          ELSE
!--          Downward and vertical surfaces
             surf%qsws(m)  = - f_qsws * ( qv1 - q_s + dq_s_dt * t_surf_green%val(m)               &
                                                   - dq_s_dt * t_surf_green_p%val(m) )
             surf%qsws(m) = surf%qsws(m) / l_v
             surf%qsws_veg(m)  = - f_qsws_veg  * ( qv1 - q_s + dq_s_dt * t_surf_green%val(m)      &
                                                   - dq_s_dt * t_surf_green_p%val(m) )
             surf%r_s(m) = - rho_lv * ( qv1 - q_s + dq_s_dt * t_surf_green%val(m)                 &
                                                  - dq_s_dt * t_surf_green_p%val(m) ) /           &
                                                   (surf%qsws(m) + 1.0E-20)  - surf%r_a_green(m)
             surf%qsws_liq(m) = 0._wp  ! - f_qsws_liq  * ( qv1 - q_s + dq_s_dt * t_surf_green_h(m)&
!                                                  - dq_s_dt * t_surf_green_h_p(m) )
!--          If the air is saturated, check the reservoir water level
             IF ( surf%qsws(m) < 0.0_wp )  THEN
!--             In case qsws_veg becomes negative (unphysical behavior), let the water enter the
!--             liquid water reservoir as dew on the plant
                IF ( surf%qsws_veg(m) < 0.0_wp )  THEN
                   surf%qsws_veg(m) = 0.0_wp
                ENDIF
             ENDIF
          ENDIF
       ELSE
          surf%r_s(m) = 1.0E10_wp
       ENDIF

    ENDDO

!
!-- pt and shf are defined on nxlg:nxrg,nysg:nyng .Get the borders from neighbours.
    CALL exchange_horiz( pt, nbgp )
!
!-- Calculation of force_radiation_call:
!-- Make logical OR for all processes.
!-- Force radiation call if at least one processor forces it.
    IF ( intermediate_timestep_count == intermediate_timestep_count_max-1 )  THEN
#if defined( __parallel )
       IF ( .NOT. force_radiation_call ) THEN
          IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
          CALL MPI_ALLREDUCE( force_radiation_call_l, force_radiation_call,                       &
                              1, MPI_LOGICAL, MPI_LOR, comm2d, ierr )
       ENDIF
#else
       force_radiation_call = force_radiation_call .OR. force_radiation_call_l
#endif
       force_radiation_call_l = .FALSE.
    ENDIF

    IF ( debug_output_timestep )  THEN
       WRITE( debug_string, * ) 'usm_surface_energy_balance: ', horizontal, l, during_spinup
       CALL debug_message( debug_string, 'end' )
    ENDIF

 END SUBROUTINE usm_surface_energy_balance

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Swapping of time levels for t_surf and t_wall called out from subroutine swap_timelevel
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_swap_timelevel( mod_count )

    IMPLICIT NONE

    INTEGER(iwp), INTENT(IN)  ::  mod_count  !<


    SELECT CASE ( mod_count )

       CASE ( 0 )
!
!--      Horizontal surfaces
          t_surf_wall_h    => t_surf_wall_h_1;   t_surf_wall_h_p    => t_surf_wall_h_2
          t_wall_h         => t_wall_h_1;        t_wall_h_p         => t_wall_h_2
          t_surf_window_h  => t_surf_window_h_1; t_surf_window_h_p  => t_surf_window_h_2
          t_window_h       => t_window_h_1;      t_window_h_p       => t_window_h_2
          t_surf_green_h   => t_surf_green_h_1;  t_surf_green_h_p   => t_surf_green_h_2
          t_green_h        => t_green_h_1;       t_green_h_p        => t_green_h_2
!
!--      Vertical surfaces
          t_surf_wall_v    => t_surf_wall_v_1;   t_surf_wall_v_p    => t_surf_wall_v_2
          t_wall_v         => t_wall_v_1;        t_wall_v_p         => t_wall_v_2
          t_surf_window_v  => t_surf_window_v_1; t_surf_window_v_p  => t_surf_window_v_2
          t_window_v       => t_window_v_1;      t_window_v_p       => t_window_v_2
          t_surf_green_v   => t_surf_green_v_1;  t_surf_green_v_p   => t_surf_green_v_2
          t_green_v        => t_green_v_1;       t_green_v_p        => t_green_v_2
       CASE ( 1 )
!
!--      Horizontal surfaces
          t_surf_wall_h    => t_surf_wall_h_2;   t_surf_wall_h_p    => t_surf_wall_h_1
          t_wall_h         => t_wall_h_2;        t_wall_h_p         => t_wall_h_1
          t_surf_window_h  => t_surf_window_h_2; t_surf_window_h_p  => t_surf_window_h_1
          t_window_h       => t_window_h_2;      t_window_h_p       => t_window_h_1
          t_surf_green_h   => t_surf_green_h_2;  t_surf_green_h_p   => t_surf_green_h_1
          t_green_h        => t_green_h_2;       t_green_h_p        => t_green_h_1
!
!--      Vertical surfaces
          t_surf_wall_v    => t_surf_wall_v_2;   t_surf_wall_v_p    => t_surf_wall_v_1
          t_wall_v         => t_wall_v_2;        t_wall_v_p         => t_wall_v_1
          t_surf_window_v  => t_surf_window_v_2; t_surf_window_v_p  => t_surf_window_v_1
          t_window_v       => t_window_v_2;      t_window_v_p       => t_window_v_1
          t_surf_green_v   => t_surf_green_v_2;  t_surf_green_v_p   => t_surf_green_v_1
          t_green_v        => t_green_v_2;       t_green_v_p        => t_green_v_1
    END SELECT

 END SUBROUTINE usm_swap_timelevel

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine writes t_surf and t_wall data into restart files
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_wrd_local


    IMPLICIT NONE

    CHARACTER(LEN=1)  ::  dum  !< dummy string to create output-variable name

    INTEGER(iwp)  ::  l  !< index surface type orientation

    INTEGER(iwp), DIMENSION(nys:nyn,nxl:nxr)  ::  global_start_index  !< index for surface data (MPI-IO)

    LOGICAL  ::  surface_data_to_write  !< switch for MPI-I/O if PE has surface data to write


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'ns_h_on_file_usm' )
       WRITE ( 14 )  surf_usm_h(0:1)%ns

       CALL wrd_write_string( 'ns_v_on_file_usm' )
       WRITE ( 14 )  surf_usm_v(0:3)%ns

       DO  l = 0, 1

          CALL wrd_write_string( 'usm_start_index_h' )
          WRITE ( 14 )  surf_usm_h(l)%start_index

          CALL wrd_write_string( 'usm_end_index_h' )
          WRITE ( 14 )  surf_usm_h(l)%end_index

          WRITE( dum, '(I1)')  l

          CALL wrd_write_string( 't_surf_wall_h(' // dum // ')' )
          WRITE ( 14 )  t_surf_wall_h(l)%val

          CALL wrd_write_string( 't_surf_window_h(' // dum // ')' )
          WRITE ( 14 )  t_surf_window_h(l)%val

          CALL wrd_write_string( 't_surf_green_h(' // dum // ')' )
          WRITE ( 14 )  t_surf_green_h(l)%val

          CALL wrd_write_string( 'm_liq_usm_h(' // dum // ')' )
          WRITE ( 14 )  m_liq_usm_h(l)%val
!
!--       Write restart data which is especially needed for the urban-surface model. In order to do not
!--       fill up the restart routines in surface_mod. Output of waste heat from indoor model. Restart
!--       data is required in this special case, because the indoor model, where waste heat is
!--       computed, is called each hour (current default), so that waste heat would have zero value
!--       until next call of indoor model.
          IF ( indoor_model )  THEN
             CALL wrd_write_string( 'waste_heat_h(' // dum // ')' )
             WRITE ( 14 )  surf_usm_h(l)%waste_heat
          ENDIF
       ENDDO

       DO  l = 0, 3

          CALL wrd_write_string( 'usm_start_index_v' )
          WRITE ( 14 )  surf_usm_v(l)%start_index

          CALL wrd_write_string( 'usm_end_index_v' )
          WRITE ( 14 )  surf_usm_v(l)%end_index

          WRITE( dum, '(I1)')  l

          CALL wrd_write_string( 't_surf_wall_v(' // dum // ')' )
          WRITE ( 14 )  t_surf_wall_v(l)%val

          CALL wrd_write_string( 't_surf_window_v(' // dum // ')' )
          WRITE ( 14 ) t_surf_window_v(l)%val

          CALL wrd_write_string( 't_surf_green_v(' // dum // ')' )
          WRITE ( 14 ) t_surf_green_v(l)%val

          IF ( indoor_model )  THEN
             CALL wrd_write_string( 'waste_heat_v(' // dum // ')' )
             WRITE ( 14 )  surf_usm_v(l)%waste_heat
          ENDIF

       ENDDO

       DO  l = 0, 1

          CALL wrd_write_string( 'usm_start_index_h' )
          WRITE ( 14 )  surf_usm_h(l)%start_index

          CALL wrd_write_string( 'usm_end_index_h' )
          WRITE ( 14 )  surf_usm_h(l)%end_index

          WRITE( dum, '(I1)')  l

          CALL wrd_write_string( 't_wall_h(' // dum // ')' )
          WRITE ( 14 )  t_wall_h(l)%val

          CALL wrd_write_string( 't_window_h(' // dum // ')' )
          WRITE ( 14 )  t_window_h(l)%val

          CALL wrd_write_string( 't_green_h(' // dum // ')' )
          WRITE ( 14 )  t_green_h(l)%val

       ENDDO

       DO  l = 0, 3

          CALL wrd_write_string( 'usm_start_index_v' )
          WRITE ( 14 )  surf_usm_v(l)%start_index

          CALL wrd_write_string( 'usm_end_index_v' )
          WRITE ( 14 )  surf_usm_v(l)%end_index

          WRITE( dum, '(I1)')  l

          CALL wrd_write_string( 't_wall_v(' // dum // ')' )
          WRITE ( 14 )  t_wall_v(l)%val

          CALL wrd_write_string( 't_window_v(' // dum // ')' )
          WRITE ( 14 )  t_window_v(l)%val

          CALL wrd_write_string( 't_green_v(' // dum // ')' )
          WRITE ( 14 )  t_green_v(l)%val

       ENDDO

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
!
!--    There is no information about the PE-grid necessary because the restart files consists of the
!--    whole domain. Therefore, ns_h_on_file_usm and ns_v_on_file_usm are not used with MPI-IO.
       DO  l = 0, 1

          WRITE( dum, '(I1)')  l

          CALL rd_mpi_io_surface_filetypes( surf_usm_h(l)%start_index, surf_usm_h(l)%end_index,             &
                                            surface_data_to_write, global_start_index )

          CALL wrd_mpi_io( 'usm_start_index_h_' // dum,  surf_usm_h(l)%start_index )
          CALL wrd_mpi_io( 'usm_end_index_h_' // dum, surf_usm_h(l)%end_index )
          CALL wrd_mpi_io( 'usm_global_start_h_' // dum, global_start_index )

          IF ( .NOT. surface_data_to_write )  CYCLE

          CALL wrd_mpi_io_surface( 't_surf_wall_h(' // dum // ')',  t_surf_wall_h(l)%val )
          CALL wrd_mpi_io_surface( 't_surf_window_h(' // dum // ')', t_surf_window_h(l)%val )
          CALL wrd_mpi_io_surface( 't_surf_green_h(' // dum // ')', t_surf_green_h(l)%val )

          CALL wrd_mpi_io_surface( 'm_liq_usm_h(' // dum // ')', m_liq_usm_h(l)%val )
          IF ( indoor_model )  THEN
             CALL wrd_mpi_io_surface( 'waste_heat_h(' // dum // ')', surf_usm_h(l)%waste_heat ) ! NEED TO BE CHECKED!!!!!
          ENDIF

       ENDDO

       DO  l = 0, 3

          WRITE( dum, '(I1)')  l

          CALL rd_mpi_io_surface_filetypes( surf_usm_v(l)%start_index, surf_usm_v(l)%end_index,    &
                                            surface_data_to_write, global_start_index )

          CALL wrd_mpi_io( 'usm_start_index_v_' // dum, surf_usm_v(l)%start_index )
          CALL wrd_mpi_io( 'usm_end_index_v_' // dum, surf_usm_v(l)%end_index )
          CALL wrd_mpi_io( 'usm_global_start_v_' // dum, global_start_index )

          IF ( .NOT. surface_data_to_write )  CYCLE

          CALL wrd_mpi_io_surface( 't_surf_wall_v(' // dum // ')', t_surf_wall_v(l)%val )
          CALL wrd_mpi_io_surface( 't_surf_window_v(' // dum // ')', t_surf_window_v(l)%val )
          CALL wrd_mpi_io_surface( 't_surf_green_v(' // dum // ')', t_surf_green_v(l)%val )

       ENDDO

       DO  l = 0, 1

          WRITE( dum, '(I1)')  l

          CALL rd_mpi_io_surface_filetypes( surf_usm_h(l)%start_index, surf_usm_h(l)%end_index,             &
                                            surface_data_to_write, global_start_index )

          CALL wrd_mpi_io( 'usm_start_index_h_2_' // dum,  surf_usm_h(l)%start_index )
          CALL wrd_mpi_io( 'usm_end_index_h_2_' // dum, surf_usm_h(l)%end_index )
          CALL wrd_mpi_io( 'usm_global_start_h_2_' // dum, global_start_index )

          IF ( .NOT. surface_data_to_write )  CYCLE

          CALL wrd_mpi_io_surface( 't_wall_h(' // dum // ')', t_wall_h(l)%val )
          CALL wrd_mpi_io_surface( 't_window_h(' // dum // ')', t_window_h(l)%val )
          CALL wrd_mpi_io_surface( 't_green_h(' // dum // ')', t_green_h(l)%val )

       ENDDO

       DO  l = 0, 3

          WRITE( dum, '(I1)')  l

          CALL rd_mpi_io_surface_filetypes( surf_usm_v(l)%start_index, surf_usm_v(l)%end_index,    &
                                            surface_data_to_write, global_start_index )

          CALL wrd_mpi_io( 'usm_start_index_v_2_' //dum, surf_usm_v(l)%start_index )
          CALL wrd_mpi_io( 'usm_end_index_v_2_' // dum, surf_usm_v(l)%end_index )
          CALL wrd_mpi_io( 'usm_global_start_v_2_' // dum, global_start_index )

          IF ( .NOT. surface_data_to_write )  CYCLE

          CALL wrd_mpi_io_surface( 't_wall_v(' // dum // ')', t_wall_v(l)%val )
          CALL wrd_mpi_io_surface( 't_window_v(' // dum // ')', t_window_v(l)%val )
          CALL wrd_mpi_io_surface( 't_green_v(' // dum // ')', t_green_v(l)%val )

       ENDDO

    ENDIF

 END SUBROUTINE usm_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Define building properties
!> Parameters 12, 13, 119 - 135 exclusive used in indoor_model_mod.f90
!> Parameters 0-11, 14-118, 136 - 149 exclusive used in urban_surface_mod.f90
!> Parameters 31, 44 used in indoor_model_mod.f90 and urban_surface_mod.f90
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE usm_define_pars
!
!-- Define the building_pars
    building_pars(:,1) = (/                                                                        &
       0.82_wp,        &  !< parameter 0   - [-] wall fraction above ground floor level
       0.18_wp,        &  !< parameter 1   - [-] window fraction above ground floor level
       0.0_wp,         &  !< parameter 2   - [-] green fraction above ground floor level
       0.0_wp,         &  !< parameter 3   - [-] green fraction roof above ground floor level
       1.5_wp,         &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
       1.5_wp,         &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
       1520000.0_wp,   &  !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
       1512000.0_wp,   &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
       1512000.0_wp,   &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
       0.93_wp,        &  !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
       0.81_wp,        &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level       
       0.81_wp,        &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
       299.15_wp,      &  !< parameter 12  - [K] indoor target summer temperature
       293.15_wp,      &  !< parameter 13  - [K] indoor target winter temperature 
       0.93_wp,        &  !< parameter 14  - [-] wall emissivity above ground floor level
       0.86_wp,        &  !< parameter 15  - [-] green emissivity above ground floor level
       0.91_wp,        &  !< parameter 16  - [-] window emissivity above ground floor level
       0.7_wp,         &  !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
       0.001_wp,       &  !< parameter 18  - [m] z0 roughness above ground floor level
       0.0001_wp,      &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
       2.9_wp,         &  !< parameter 20  - [m] ground floor level height
       0.82_wp,        &  !< parameter 21  - [-] wall fraction ground floor level
       0.18_wp,        &  !< parameter 22  - [-] window fraction ground floor level
       0.0_wp,         &  !< parameter 23  - [-] green fraction ground floor level
       0.0_wp,         &  !< parameter 24  - [-] green fraction roof ground floor level
       1.5_wp,         &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
       1520000.0_wp,   &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
       1512000.0_wp,   &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
       1512000.0_wp,   &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
       0.81_wp,        &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
       0.81_wp,        &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 32  - [-] wall emissivity ground floor level
       0.91_wp,        &  !< parameter 33  - [-] window emissivity ground floor level
       0.86_wp,        &  !< parameter 34  - [-] green emissivity ground floor level
       0.7_wp,         &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
       0.001_wp,       &  !< parameter 36  - [m] z0 roughness ground floor level
       0.0001_wp,      &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
       36.0_wp,        &  !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
       37.0_wp,        &  !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 41  - [m] 1st wall layer thickness above ground floor level
       0.2_wp,         &  !< parameter 42  - [m] 2nd wall layer thickness above ground floor level
       0.38_wp,        &  !< parameter 43  - [m] 3rd wall layer thickness above ground floor level
       0.4_wp,         &  !< parameter 44  - [m] 4th wall layer thickness above ground floor level
       20000.0_wp,     &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
       23.0_wp,        &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
       20000.0_wp,     &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
       20000.0_wp,     &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
       23.0_wp,        &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
       10.0_wp,        &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
       1.0_wp,         &  !< parameter 51  - [-] wall fraction ground plate
       0.18_wp,        &  !< parameter 52  - [m] 1st wall layer thickness ground plate
       0.36_wp,        &  !< parameter 53  - [m] 2nd wall layer thickness ground plate
       0.42_wp,        &  !< parameter 54  - [m] 3rd wall layer thickness ground plate
       0.45_wp,        &  !< parameter 55  - [m] 4th wall layer thickness ground plate
       1512000.0_wp,   &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
       1512000.0_wp,   &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
       2112000.0_wp,   &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
       0.52_wp,        &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
       0.52_wp,        &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
       2.1_wp,         &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
       0.02_wp,        &  !< parameter 62  - [m] 1st wall layer thickness ground floor level
       0.2_wp,         &  !< parameter 63  - [m] 2nd wall layer thickness ground floor level
       0.38_wp,        &  !< parameter 64  - [m] 3rd wall layer thickness ground floor level
       0.4_wp,         &  !< parameter 65  - [m] 4th wall layer thickness ground floor level
       36.0_wp,        &  !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 67  - [m] 1st window layer thickness ground floor level
       0.04_wp,        &  !< parameter 68  - [m] 2nd window layer thickness ground floor level
       0.06_wp,        &  !< parameter 69  - [m] 3rd window layer thickness ground floor level
       0.08_wp,        &  !< parameter 70  - [m] 4th window layer thickness ground floor level
       1736000.0_wp,   &  !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
       1736000.0_wp,   &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
       1736000.0_wp,   &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
       0.45_wp,        &  !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
       0.45_wp,        &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
       0.45_wp,        &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
       37.0_wp,        &  !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 79  - [m] 1st window layer thickness above ground floor level
       0.04_wp,        &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
       0.06_wp,        &  !< parameter 81  - [m] 3rd window layer thickness above ground floor level
       0.08_wp,        &  !< parameter 82  - [m] 4th window layer thickness above ground floor level
       1736000.0_wp,   &  !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
       1736000.0_wp,   &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
       1736000.0_wp,   &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
       0.45_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
       0.45_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
       0.45_wp,        &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
       1.0_wp,         &  !< parameter 89  - [-] wall fraction roof
       0.02_wp,        &  !< parameter 90  - [m] 1st wall layer thickness roof
       0.06_wp,        &  !< parameter 91  - [m] 2nd wall layer thickness roof
       0.08_wp,        &  !< parameter 92  - [m] 3rd wall layer thickness roof
       0.1_wp,         &  !< parameter 93  - [m] 4th wall layer thickness roof
       1512000.0_wp,   &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
       709650.0_wp,    &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
       709650.0_wp,    &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
       0.52_wp,        &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
       0.12_wp,        &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
       0.12_wp,        &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
       0.90_wp,        &  !< parameter 100 - [-] wall emissivity roof
       42.0_wp,        &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 102 - [-] window fraction roof
       0.02_wp,        &  !< parameter 103 - [m] window 1st layer thickness roof
       0.04_wp,        &  !< parameter 104 - [m] window 2nd layer thickness roof
       0.06_wp,        &  !< parameter 105 - [m] window 3rd layer thickness roof
       0.08_wp,        &  !< parameter 106 - [m] window 4th layer thickness roof
       1736000.0_wp,   &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
       1736000.0_wp,   &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
       1736000.0_wp,   &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
       0.45_wp,        &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
       0.45_wp,        &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
       0.45_wp,        &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
       0.91_wp,        &  !< parameter 113 - [-] window emissivity roof
       0.7_wp,         &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
       37.0_wp,        &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
       0.86_wp,        &  !< parameter 116 - [-] green emissivity roof
       5.0_wp,         &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 118 - [-] green type roof
       0.75_wp,        &  !< parameter 119 - [-] shading factor
       0.8_wp,         &  !< parameter 120 - [-] g-value windows
       2.9_wp,         &  !< parameter 121 - [W/(m2*K)] u-value windows
       0.5_wp,         &  !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 0.5_wp, winter 0.5
       2.0_wp,         &  !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.5_wp, winter 0.0
       0.0_wp,         &  !< parameter 124 - [-] heat recovery efficiency
       3.0_wp,         &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
       260000.0_wp,    &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
       4.5_wp,         &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
       100.0_wp,       &  !< parameter 128 - [W] maximal heating capacity
       0.0_wp,         &  !< parameter 129 - [W] maximal cooling capacity
       0.0_wp,         &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
       4.2_wp,         &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
       2.9_wp,         &  !< parameter 132 - [m] storey height
       0.2_wp,         &  !< parameter 133 - [m] ceiling construction height
       0.1_wp,         &  !< parameter 134 - [-] anthropogenic heat output for heating
       1.333_wp,       &  !< parameter 135 - [-] anthropogenic heat output for cooling
       1526000.0_wp,   &  !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
       0.7_wp,         &  !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
       0.7_wp,         &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
       709650.0_wp,    &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
       0.12_wp,        &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
       1736000.0_wp,   &  !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
       0.45_wp,        &  !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
       1736000.0_wp,   &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
       0.45_wp,        &  !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
       0.7_wp,         &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
       1736000.0_wp,   &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
       0.45_wp         &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                        /)

    building_pars(:,2) = (/                                                                        &                        
       0.75_wp,        &  !< parameter 0   - [-] wall fraction above ground floor level
       0.25_wp,        &  !< parameter 1   - [-] window fraction above ground floor level
       0.0_wp,         &  !< parameter 2   - [-] green fraction above ground floor level
       0.0_wp,         &  !< parameter 3   - [-] green fraction roof above ground floor level
       1.5_wp,         &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
       1.5_wp,         &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
       1520000.0_wp,   &  !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
       79200.0_wp,     &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
       2112000.0_wp,   &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
       0.93_wp,        &  !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
       0.046_wp,       &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level       
       2.1_wp,         &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
       299.15_wp,      &  !< parameter 12  - [K] indoor target summer temperature
       293.15_wp,      &  !< parameter 13  - [K] indoor target winter temperature
       0.93_wp,        &  !< parameter 14  - [-] wall emissivity above ground floor level
       0.86_wp,        &  !< parameter 15  - [-] green emissivity above ground floor level
       0.87_wp,        &  !< parameter 16  - [-] window emissivity above ground floor level
       0.65_wp,        &  !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
       0.001_wp,       &  !< parameter 18  - [m] z0 roughness above ground floor level
       0.0001_wp,      &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
       2.5_wp,         &  !< parameter 20  - [m] ground floor level height
       0.75_wp,        &  !< parameter 21  - [-] wall fraction ground floor level
       0.25_wp,        &  !< parameter 22  - [-] window fraction ground floor level
       0.0_wp,         &  !< parameter 23  - [-] green fraction ground floor level
       0.0_wp,         &  !< parameter 24  - [-] green fraction roof ground floor level
       1.5_wp,         &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
       1520000.0_wp,   &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
       79200.0_wp,     &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
       2112000.0_wp,   &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
       0.046_wp,       &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
       2.1_wp,         &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 32  - [-] wall emissivity ground floor level
       0.87_wp,        &  !< parameter 33  - [-] window emissivity ground floor level
       0.86_wp,        &  !< parameter 34  - [-] green emissivity ground floor level
       0.65_wp,        &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
       0.001_wp,       &  !< parameter 36  - [m] z0 roughness ground floor level
       0.0001_wp,      &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
       36.0_wp,        &  !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
       37.0_wp,        &  !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 41  - [m] 1st wall layer thickness above ground floor level
       0.08_wp,        &  !< parameter 42  - [m] 2nd wall layer thickness above ground floor level
       0.32_wp,        &  !< parameter 43  - [m] 3rd wall layer thickness above ground floor level
       0.34_wp,        &  !< parameter 44  - [m] 4th wall layer thickness above ground floor level
       20000.0_wp,     &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
       23.0_wp,        &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
       20000.0_wp,     &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
       20000.0_wp,     &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
       23.0_wp,        &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
       10.0_wp,        &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
       1.0_wp,         &  !< parameter 51  - [-] wall fraction ground plate
       0.20_wp,        &  !< parameter 52  - [m] 1st wall layer thickness ground plate
       0.26_wp,        &  !< parameter 53  - [m] 2nd wall layer thickness ground plate
       0.32_wp,        &  !< parameter 54  - [m] 3rd wall layer thickness ground plate
       0.34_wp,        &  !< parameter 55  - [m] 4th wall layer thickness ground plate
       2112000.0_wp,   &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
       79200.0_wp,     &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
       2112000.0_wp,   &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
       2.1_wp,         &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
       0.05_wp,        &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
       2.1_wp,         &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
       0.02_wp,        &  !< parameter 62  - [m] 1st wall layer thickness ground floor level
       0.08_wp,        &  !< parameter 63  - [m] 2nd wall layer thickness ground floor level
       0.32_wp,        &  !< parameter 64  - [m] 3rd wall layer thickness ground floor level
       0.34_wp,        &  !< parameter 65  - [m] 4th wall layer thickness ground floor level
       36.0_wp,        &  !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 67  - [m] 1st window layer thickness ground floor level
       0.04_wp,        &  !< parameter 68  - [m] 2nd window layer thickness ground floor level
       0.06_wp,        &  !< parameter 69  - [m] 3rd window layer thickness ground floor level
       0.08_wp,        &  !< parameter 70  - [m] 4th window layer thickness ground floor level
       1736000.0_wp,   &  !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
       1736000.0_wp,   &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
       1736000.0_wp,   &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
       0.19_wp,        &  !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
       0.19_wp,        &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
       0.19_wp,        &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
       37.0_wp,        &  !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 79  - [m] 1st window layer thickness above ground floor level
       0.04_wp,        &  !< parameter 80  - [m] 2nd window layer thickness above ground floor level
       0.06_wp,        &  !< parameter 81  - [m] 3rd window layer thickness above ground floor level
       0.08_wp,        &  !< parameter 82  - [m] 4th window layer thickness above ground floor level
       1736000.0_wp,   &  !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
       1736000.0_wp,   &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
       1736000.0_wp,   &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
       0.19_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
       0.19_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
       0.19_wp,        &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
       1.0_wp,         &  !< parameter 89  - [-] wall fraction roof
       0.02_wp,        &  !< parameter 90  - [m] 1st wall layer thickness roof
       0.17_wp,        &  !< parameter 91  - [m] 2nd wall layer thickness roof
       0.37_wp,        &  !< parameter 92  - [m] 3rd wall layer thickness roof
       0.39_wp,        &  !< parameter 93  - [m] 4th wall layer thickness roof
       1700000.0_wp,   &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
       79200.0_wp,     &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
       2112000.0_wp,   &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
       0.16_wp,        &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
       0.046_wp,       &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
       2.1_wp,         &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
       0.93_wp,        &  !< parameter 100 - [-] wall emissivity roof
       42.0_wp,        &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 102 - [-] window fraction roof
       0.02_wp,        &  !< parameter 103 - [m] window 1st layer thickness roof
       0.04_wp,        &  !< parameter 104 - [m] window 2nd layer thickness roof
       0.06_wp,        &  !< parameter 105 - [m] window 3rd layer thickness roof
       0.08_wp,        &  !< parameter 106 - [m] window 4th layer thickness roof
       1736000.0_wp,   &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
       1736000.0_wp,   &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
       1736000.0_wp,   &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
       0.19_wp,        &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
       0.19_wp,        &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
       0.19_wp,        &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
       0.87_wp,        &  !< parameter 113 - [-] window emissivity roof
       0.65_wp,        &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
       37.0_wp,        &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
       0.86_wp,        &  !< parameter 116 - [-] green emissivity roof
       5.0_wp,         &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 118 - [-] green type roof
       0.75_wp,        &  !< parameter 119 - [-] shading factor
       0.7_wp,         &  !< parameter 120 - [-] g-value windows
       1.7_wp,         &  !< parameter 121 - [W/(m2*K)] u-value windows
       0.5_wp,         &  !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 0.5_wp, winter 0.5
       1.5_wp,         &  !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.5_wp, winter 0.0
       0.0_wp,         &  !< parameter 124 - [-] heat recovery efficiency
       3.5_wp,         &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
       370000.0_wp,    &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
       4.5_wp,         &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
       80.0_wp,        &  !< parameter 128 - [W] maximal heating capacity
       0.0_wp,         &  !< parameter 129 - [W] maximal cooling capacity
       0.0_wp,         &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
       4.2_wp,         &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
       2.5_wp,         &  !< parameter 132 - [m] storey height
       0.2_wp,         &  !< parameter 133 - [m] ceiling construction height
       0.0_wp,         &  !< parameter 134 - [-] anthropogenic heat output for heating
       2.54_wp,        &  !< parameter 135 - [-] anthropogenic heat output for cooling
       1526000.0_wp,   &  !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
       0.7_wp,         &  !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
       0.7_wp,         &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
       357200.0_wp,    &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
       0.04_wp,        &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
       1736000.0_wp,   &  !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
       0.19_wp,        &  !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
       1736000.0_wp,   &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
       0.19_wp,        &  !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
       0.7_wp,         &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
       1736000.0_wp,   &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
       0.19_wp         &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                        /)

    building_pars(:,3) = (/                                                                        &                        
       0.71_wp,        &  !< parameter 0   - [-] wall fraction above ground floor level
       0.29_wp,        &  !< parameter 1   - [-] window fraction above ground floor level
       0.0_wp,         &  !< parameter 2   - [-] green fraction above ground floor level
       0.0_wp,         &  !< parameter 3   - [-] green fraction roof above ground floor level
       1.5_wp,         &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
       1.5_wp,         &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
       1520000.0_wp,   &  !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
       79200.0_wp,     &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
       1344000.0_wp,   &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
       0.93_wp,        &  !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
       0.035_wp,       &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level       
       0.68_wp,        &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
       299.15_wp,      &  !< parameter 12  - [K] indoor target summer temperature
       293.15_wp,      &  !< parameter 13  - [K] indoor target winter temperature
       0.93_wp,        &  !< parameter 14  - [-] wall emissivity above ground floor level
       0.86_wp,        &  !< parameter 15  - [-] green emissivity above ground floor level
       0.8_wp,         &  !< parameter 16  - [-] window emissivity above ground floor level
       0.57_wp,        &  !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
       0.001_wp,       &  !< parameter 18  - [m] z0 roughness above ground floor level
       0.0001_wp,      &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
       2.7_wp,         &  !< parameter 20  - [m] ground floor level height
       0.71_wp,        &  !< parameter 21  - [-] wall fraction ground floor level
       0.29_wp,        &  !< parameter 22  - [-] window fraction ground floor level
       0.0_wp,         &  !< parameter 23  - [-] green fraction ground floor level
       0.0_wp,         &  !< parameter 24  - [-] green fraction roof ground floor level
       1.5_wp,         &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
       1520000.0_wp,   &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
       79200.0_wp,     &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
       1344000.0_wp,   &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
       0.035_wp,       &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
       0.68_wp,        &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 32  - [-] wall emissivity ground floor level
       0.8_wp,         &  !< parameter 33  - [-] window emissivity ground floor level
       0.86_wp,        &  !< parameter 34  - [-] green emissivity ground floor level
       0.57_wp,        &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
       0.001_wp,       &  !< parameter 36  - [m] z0 roughness ground floor level
       0.0001_wp,      &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
       36.0_wp,        &  !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
       38.0_wp,        &  !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 41  - [m] 1st wall layer thickness above ground floor level
       0.22_wp,        &  !< parameter 42  - [m] 2nd wall layer thickness above ground floor level
       0.58_wp,        &  !< parameter 43  - [m] 3rd wall layer thickness above ground floor level
       0.6_wp,         &  !< parameter 44  - [m] 4th wall layer thickness above ground floor level
       20000.0_wp,     &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
       23.0_wp,        &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
       20000.0_wp,     &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
       20000.0_wp,     &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
       23.0_wp,        &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
       10.0_wp,        &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
       1.0_wp,         &  !< parameter 51  - [-] wall fraction ground plate
       0.20_wp,        &  !< parameter 52  - [m] 1st wall layer thickness ground plate
       0.32_wp,        &  !< parameter 53  - [m] 2nd wall layer thickness ground plate
       0.38_wp,        &  !< parameter 54  - [m] 3rd wall layer thickness ground plate
       0.41_wp,        &  !< parameter 55  - [m] 4th wall layer thickness ground plate
       2112000.0_wp,   &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
       79200.0_wp,     &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
       2112000.0_wp,   &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
       2.1_wp,         &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
       0.05_wp,        &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
       2.1_wp,         &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
       0.02_wp,        &  !< parameter 62  - [m] 1st wall layer thickness ground floor level
       0.22_wp,        &  !< parameter 63  - [m] 2nd wall layer thickness ground floor level
       0.58_wp,        &  !< parameter 64  - [m] 3rd wall layer thickness ground floor level
       0.6_wp,         &  !< parameter 65  - [m] 4th wall layer thickness ground floor level
       36.0_wp,        &  !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
       0.03_wp,        &  !< parameter 67  - [m] 1st window layer thickness ground floor level
       0.06_wp,        &  !< parameter 68  - [m] 2nd window layer thickness ground floor level
       0.09_wp,        &  !< parameter 69  - [m] 3rd window layer thickness ground floor level
       0.012_wp,       &  !< parameter 70  - [m] 4th window layer thickness ground floor level
       1736000.0_wp,   &  !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
       1736000.0_wp,   &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
       1736000.0_wp,   &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
       0.11_wp,        &  !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
       0.11_wp,        &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
       0.11_wp,        &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
       38.0_wp,        &  !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
       0.03_wp,        &  !< parameter 79  - [m] 1st window layer thickness above ground floor level
       0.06_wp,        &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
       0.09_wp,        &  !< parameter 81  - [m] 3rd window layer thickness above ground floor level
       0.12_wp,        &  !< parameter 82  - [m] 4th window layer thickness above ground floor level
       1736000.0_wp,   &  !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
       1736000.0_wp,   &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
       1736000.0_wp,   &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
       0.11_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
       0.11_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
       0.11_wp,        &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
       1.0_wp,         &  !< parameter 89  - [-] wall fraction roof
       0.02_wp,        &  !< parameter 90  - [m] 1st wall layer thickness roof
       0.06_wp,        &  !< parameter 91  - [m] 2nd wall layer thickness roof
       0.36_wp,        &  !< parameter 92  - [m] 3rd wall layer thickness roof
       0.38_wp,        &  !< parameter 93  - [m] 4th wall layer thickness roof
       3753600.0_wp,   &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
       709650.0_wp,    &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
       79200.0_wp,     &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
       0.52_wp,        &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
       0.12_wp,        &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
       0.035_wp,       &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
       0.93_wp,        &  !< parameter 100 - [-] wall emissivity roof
       42.0_wp,        &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 102 - [-] window fraction roof
       0.03_wp,        &  !< parameter 103 - [m] window 1st layer thickness roof
       0.06_wp,        &  !< parameter 104 - [m] window 2nd layer thickness roof
       0.09_wp,        &  !< parameter 105 - [m] window 3rd layer thickness roof
       0.12_wp,        &  !< parameter 106 - [m] window 4th layer thickness roof
       1736000.0_wp,   &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
       1736000.0_wp,   &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
       1736000.0_wp,   &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
       0.11_wp,        &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
       0.11_wp,        &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
       0.11_wp,        &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
       0.8_wp,         &  !< parameter 113 - [-] window emissivity roof
       0.57_wp,        &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
       38.0_wp,        &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
       0.86_wp,        &  !< parameter 116 - [-] green emissivity roof
       5.0_wp,         &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 118 - [-] green type roof
       0.15_wp,        &  !< parameter 119 - [-] shading factor
       0.6_wp,         &  !< parameter 120 - [-] g-value windows
       0.8_wp,         &  !< parameter 121 - [W/(m2*K)] u-value windows
       0.5_wp,         &  !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 0.5_wp, winter 0.5_wp
       1.5_wp,         &  !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.5_wp, winter 0.0_wp
       0.8_wp,         &  !< parameter 124 - [-] heat recovery efficiency
       2.5_wp,         &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
       165000.0_wp,    &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
       4.5_wp,         &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
       40.0_wp,        &  !< parameter 128 - [W] maximal heating capacity
       0.0_wp,         &  !< parameter 129 - [W] maximal cooling capacity 
       0.0_wp,         &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
       4.2_wp,         &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
       2.7_wp,         &  !< parameter 132 - [m] storey height
       0.2_wp,         &  !< parameter 133 - [m] ceiling construction height
       -2.0_wp,        &  !< parameter 134 - [-] anthropogenic heat output for heating
       1.25_wp,        &  !< parameter 135 - [-] anthropogenic heat output for cooling
       1526000.0_wp,   &  !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
       0.7_wp,         &  !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
       0.7_wp,         &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
       709650.0_wp,    &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
       0.12_wp,        &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
       1736000.0_wp,   &  !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
       0.11_wp,        &  !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
       1736000.0_wp,   &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
       0.11_wp,        &  !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
       0.7_wp,         &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
       1736000.0_wp,   &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
       0.11_wp         &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                        /)

   building_pars(:,4) = (/                                                                        &
       0.82_wp,        &  !< parameter 0   - [-] wall fraction above ground floor level
       0.18_wp,        &  !< parameter 1   - [-] window fraction above ground floor level
       0.0_wp,         &  !< parameter 2   - [-] green fraction above ground floor level
       0.0_wp,         &  !< parameter 3   - [-] green fraction roof above ground floor level
       1.5_wp,         &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
       1.5_wp,         &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
       1520000.0_wp,   &  !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
       1512000.0_wp,   &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
       1512000.0_wp,   &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
       0.93_wp,        &  !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
       0.81_wp,        &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level       
       0.81_wp,        &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
       299.15_wp,      &  !< parameter 12  - [K] indoor target summer temperature
       293.15_wp,      &  !< parameter 13  - [K] indoor target winter temperature
       0.93_wp,        &  !< parameter 14  - [-] wall emissivity above ground floor level
       0.86_wp,        &  !< parameter 15  - [-] green emissivity above ground floor level
       0.91_wp,        &  !< parameter 16  - [-] window emissivity above ground floor level
       0.7_wp,         &  !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
       0.001_wp,       &  !< parameter 18  - [m] z0 roughness above ground floor level
       0.0001_wp,      &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
       2.9_wp,         &  !< parameter 20  - [m] ground floor level height
       0.82_wp,        &  !< parameter 21  - [-] wall fraction ground floor level
       0.18_wp,        &  !< parameter 22  - [-] window fraction ground floor level
       0.0_wp,         &  !< parameter 23  - [-] green fraction ground floor level
       0.0_wp,         &  !< parameter 24  - [-] green fraction roof ground floor level
       1.5_wp,         &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
       1520000.0_wp,   &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
       1512000.0_wp,   &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
       1512000.0_wp,   &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
       0.81_wp,        &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
       0.81_wp,        &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 32  - [-] wall emissivity ground floor level
       0.91_wp,        &  !< parameter 33  - [-] window emissivity ground floor level
       0.86_wp,        &  !< parameter 34  - [-] green emissivity ground floor level
       0.7_wp,         &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
       0.001_wp,       &  !< parameter 36  - [m] z0 roughness ground floor level
       0.0001_wp,      &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
       36.0_wp,        &  !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
       37.0_wp,        &  !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 41  - [m] 1st wall layer thickness above ground floor level
       0.2_wp,         &  !< parameter 42  - [m] 2nd wall layer thickness above ground floor level
       0.38_wp,        &  !< parameter 43  - [m] 3rd wall layer thickness above ground floor level
       0.4_wp,         &  !< parameter 44  - [m] 4th wall layer thickness above ground floor level
       20000.0_wp,     &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
       23.0_wp,        &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
       20000.0_wp,     &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
       20000.0_wp,     &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
       23.0_wp,        &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
       10.0_wp,        &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
       1.0_wp,         &  !< parameter 51  - [-] wall fraction ground plate
       0.18_wp,        &  !< parameter 52  - [m] 1st wall layer thickness ground plate
       0.36_wp,        &  !< parameter 53  - [m] 2nd wall layer thickness ground plate
       0.42_wp,        &  !< parameter 54  - [m] 3rd wall layer thickness ground plate
       0.45_wp,        &  !< parameter 55  - [m] 4th wall layer thickness ground plate
       1512000.0_wp,   &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
       1512000.0_wp,   &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
       2112000.0_wp,   &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
       0.52_wp,        &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
       0.52_wp,        &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
       2.1_wp,         &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
       0.02_wp,        &  !< parameter 62  - [m] 1st wall layer thickness ground floor level
       0.2_wp,         &  !< parameter 63  - [m] 2nd wall layer thickness ground floor level
       0.38_wp,        &  !< parameter 64  - [m] 3rd wall layer thickness ground floor level
       0.4_wp,         &  !< parameter 65  - [m] 4th wall layer thickness ground floor level
       36.0_wp,        &  !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 67  - [m] 1st window layer thickness ground floor level
       0.04_wp,        &  !< parameter 68  - [m] 2nd window layer thickness ground floor level
       0.06_wp,        &  !< parameter 69  - [m] 3rd window layer thickness ground floor level
       0.08_wp,        &  !< parameter 70  - [m] 4th window layer thickness ground floor level
       1736000.0_wp,   &  !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
       1736000.0_wp,   &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
       1736000.0_wp,   &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
       0.45_wp,        &  !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
       0.45_wp,        &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
       0.45_wp,        &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
       37.0_wp,        &  !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 79  - [m] 1st window layer thickness above ground floor level
       0.04_wp,        &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
       0.06_wp,        &  !< parameter 81  - [m] 3rd window layer thickness above ground floor level
       0.08_wp,        &  !< parameter 82  - [m] 4th window layer thickness above ground floor level
       1736000.0_wp,   &  !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
       1736000.0_wp,   &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
       1736000.0_wp,   &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
       0.45_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
       0.45_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
       0.45_wp,        &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
       1.0_wp,         &  !< parameter 89  - [-] wall fraction roof
       0.02_wp,        &  !< parameter 90  - [m] 1st wall layer thickness roof
       0.06_wp,        &  !< parameter 91  - [m] 2nd wall layer thickness roof
       0.08_wp,        &  !< parameter 92  - [m] 3rd wall layer thickness roof
       0.1_wp,         &  !< parameter 93  - [m] 4th wall layer thickness roof
       1512000.0_wp,   &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
       709650.0_wp,    &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
       709650.0_wp,    &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
       0.52_wp,        &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
       0.12_wp,        &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
       0.12_wp,        &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
       0.90_wp,        &  !< parameter 100 - [-] wall emissivity roof
       42.0_wp,        &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 102 - [-] window fraction roof
       0.02_wp,        &  !< parameter 103 - [m] window 1st layer thickness roof
       0.04_wp,        &  !< parameter 104 - [m] window 2nd layer thickness roof
       0.06_wp,        &  !< parameter 105 - [m] window 3rd layer thickness roof
       0.08_wp,        &  !< parameter 106 - [m] window 4th layer thickness roof
       1736000.0_wp,   &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
       1736000.0_wp,   &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
       1736000.0_wp,   &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
       0.45_wp,        &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
       0.45_wp,        &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
       0.45_wp,        &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
       0.91_wp,        &  !< parameter 113 - [-] window emissivity roof
       0.7_wp,         &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
       37.0_wp,        &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
       0.86_wp,        &  !< parameter 116 - [-] green emissivity roof
       5.0_wp,         &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 118 - [-] green type roof
       0.75_wp,        &  !< parameter 119 - [-] shading factor
       0.8_wp,         &  !< parameter 120 - [-] g-value windows
       2.9_wp,         &  !< parameter 121 - [W/(m2*K)] u-value windows
       1.0_wp,         &  !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 1.0_wp, winter 0.2
       1.0_wp,         &  !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.0_wp, winter 0.8  
       0.0_wp,         &  !< parameter 124 - [-] heat recovery efficiency
       3.0_wp,         &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
       260000.0_wp,    &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
       4.5_wp,         &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
       100.0_wp,       &  !< parameter 128 - [W] maximal heating capacity
       -100.0_wp,      &  !< parameter 129 - [W] maximal cooling capacity
       7.0_wp,         &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
       3.0_wp,         &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
       2.9_wp,         &  !< parameter 132 - [m] storey height
       0.2_wp,         &  !< parameter 133 - [m] ceiling construction height
       0.1_wp,         &  !< parameter 134 - [-] anthropogenic heat output for heating
       1.333_wp,       &  !< parameter 135 - [-] anthropogenic heat output for cooling
       1526000.0_wp,   &  !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
       0.7_wp,         &  !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
       0.7_wp,         &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
       709650.0_wp,    &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
       0.12_wp,        &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
       1736000.0_wp,   &  !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
       0.45_wp,        &  !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
       1736000.0_wp,   &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
       0.45_wp,        &  !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
       0.7_wp,         &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
       1736000.0_wp,   &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
       0.45_wp         &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                        /)

    building_pars(:,5) = (/                                                                        &
       0.75_wp,        &  !< parameter 0   - [-] wall fraction above ground floor level
       0.25_wp,        &  !< parameter 1   - [-] window fraction above ground floor level
       0.0_wp,         &  !< parameter 2   - [-] green fraction above ground floor level
       0.0_wp,         &  !< parameter 3   - [-] green fraction roof above ground floor level
       1.5_wp,         &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
       1.5_wp,         &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
       1520000.0_wp,   &  !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
       79200.0_wp,     &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
       2112000.0_wp,   &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
       0.93_wp,        &  !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
       2.1_wp,         &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level       
       0.046_wp,       &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
       299.15_wp,      &  !< parameter 12  - [K] indoor target summer temperature
       293.15_wp,      &  !< parameter 13  - [K] indoor target winter temperature
       0.93_wp,        &  !< parameter 14  - [-] wall emissivity above ground floor level
       0.86_wp,        &  !< parameter 15  - [-] green emissivity above ground floor level
       0.87_wp,        &  !< parameter 16  - [-] window emissivity above ground floor level
       0.65_wp,        &  !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
       0.001_wp,       &  !< parameter 18  - [m] z0 roughness above ground floor level
       0.0001_wp,      &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
       2.5_wp,         &  !< parameter 20  - [m] ground floor level height
       0.75_wp,        &  !< parameter 21  - [-] wall fraction ground floor level
       0.25_wp,        &  !< parameter 22  - [-] window fraction ground floor level
       0.0_wp,         &  !< parameter 23  - [-] green fraction ground floor level
       0.0_wp,         &  !< parameter 24  - [-] green fraction roof ground floor level
       1.5_wp,         &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
       1520000.0_wp,   &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
       79200.0_wp,     &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
       2112000.0_wp,   &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
       0.046_wp,       &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
       2.1_wp,         &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 32  - [-] wall emissivity ground floor level
       0.87_wp,        &  !< parameter 33  - [-] window emissivity ground floor level
       0.86_wp,        &  !< parameter 34  - [-] green emissivity ground floor level
       0.65_wp,        &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
       0.001_wp,       &  !< parameter 36  - [m] z0 roughness ground floor level
       0.0001_wp,      &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
       36.0_wp,        &  !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
       37.0_wp,        &  !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 41  - [m] 1st wall layer thickness above ground floor level
       0.08_wp,        &  !< parameter 42  - [m] 2nd wall layer thickness above ground floor level
       0.32_wp,        &  !< parameter 43  - [m] 3rd wall layer thickness above ground floor level
       0.34_wp,        &  !< parameter 44  - [m] 4th wall layer thickness above ground floor level
       20000.0_wp,     &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
       23.0_wp,        &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
       20000.0_wp,     &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
       20000.0_wp,     &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
       23.0_wp,        &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
       10.0_wp,        &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
       1.0_wp,         &  !< parameter 51  - [-] wall fraction ground plate
       0.20_wp,        &  !< parameter 52  - [m] 1st wall layer thickness ground plate
       0.26_wp,        &  !< parameter 53  - [m] 2nd wall layer thickness ground plate
       0.32_wp,        &  !< parameter 54  - [m] 3rd wall layer thickness ground plate
       0.34_wp,        &  !< parameter 55  - [m] 4th wall layer thickness ground plate
       2112000.0_wp,   &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
       79200.0_wp,     &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
       2112000.0_wp,   &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
       2.1_wp,         &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
       0.05_wp,        &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
       2.1_wp,         &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
       0.02_wp,        &  !< parameter 62  - [m] 1st wall layer thickness ground floor level
       0.08_wp,        &  !< parameter 63  - [m] 2nd wall layer thickness ground floor level
       0.32_wp,        &  !< parameter 64  - [m] 3rd wall layer thickness ground floor level
       0.34_wp,        &  !< parameter 65  - [m] 4th wall layer thickness ground floor level
       36.0_wp,        &  !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 67  - [m] 1st window layer thickness ground floor level
       0.04_wp,        &  !< parameter 68  - [m] 2nd window layer thickness ground floor level
       0.06_wp,        &  !< parameter 69  - [m] 3rd window layer thickness ground floor level
       0.08_wp,        &  !< parameter 70  - [m] 4th window layer thickness ground floor level
       1736000.0_wp,   &  !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
       1736000.0_wp,   &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
       1736000.0_wp,   &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
       0.19_wp,        &  !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
       0.19_wp,        &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
       0.19_wp,        &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
       37.0_wp,        &  !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 79  - [m] 1st window layer thickness above ground floor level
       0.04_wp,        &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
       0.06_wp,        &  !< parameter 81  - [m] 3rd window layer thickness above ground floor level
       0.08_wp,        &  !< parameter 82  - [m] 4th window layer thickness above ground floor level
       1736000.0_wp,   &  !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
       1736000.0_wp,   &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
       1736000.0_wp,   &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
       0.19_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
       0.19_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
       0.19_wp,        &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
       1.0_wp,         &  !< parameter 89  - [-] wall fraction roof
       0.02_wp,        &  !< parameter 90  - [m] 1st wall layer thickness roof
       0.17_wp,        &  !< parameter 91  - [m] 2nd wall layer thickness roof
       0.37_wp,        &  !< parameter 92  - [m] 3rd wall layer thickness roof
       0.39_wp,        &  !< parameter 93  - [m] 4th wall layer thickness roof
       1700000.0_wp,   &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
       79200.0_wp,     &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
       2112000.0_wp,   &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
       0.16_wp,        &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
       0.046_wp,       &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
       2.1_wp,         &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
       0.93_wp,        &  !< parameter 100 - [-] wall emissivity roof
       42.0_wp,        &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 102 - [-] window fraction roof
       0.02_wp,        &  !< parameter 103 - [m] window 1st layer thickness roof
       0.04_wp,        &  !< parameter 104 - [m] window 2nd layer thickness roof
       0.06_wp,        &  !< parameter 105 - [m] window 3rd layer thickness roof
       0.08_wp,        &  !< parameter 106 - [m] window 4th layer thickness roof
       1736000.0_wp,   &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
       1736000.0_wp,   &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
       1736000.0_wp,   &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
       0.19_wp,        &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
       0.19_wp,        &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
       0.19_wp,        &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
       0.87_wp,        &  !< parameter 113 - [-] window emissivity roof
       0.65_wp,        &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
       37.0_wp,        &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
       0.86_wp,        &  !< parameter 116 - [-] green emissivity roof
       5.0_wp,         &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 118 - [-] green type roof
       0.75_wp,        &  !< parameter 119 - [-] shading factor
       0.7_wp,         &  !< parameter 120 - [-] g-value windows
       1.7_wp,         &  !< parameter 121 - [W/(m2*K)] u-value windows
       1.0_wp,         &  !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 1.0_wp, winter 0.2 
       1.0_wp,         &  !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.0_wp, winter 0.8
       0.0_wp,         &  !< parameter 124 - [-] heat recovery efficiency
       3.5_wp,         &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
       370000.0_wp,    &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
       4.5_wp,         &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
       80.0_wp,        &  !< parameter 128 - [W] maximal heating capacity
       -120.0_wp,      &  !< parameter 129 - [W] maximal cooling capacity
       7.0_wp,         &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
       3.0_wp,         &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
       2.5_wp,         &  !< parameter 132 - [m] storey height
       0.2_wp,         &  !< parameter 133 - [m] ceiling construction height
       0.0_wp,         &  !< parameter 134 - [-] anthropogenic heat output for heating
       2.54_wp,        &  !< parameter 135 - [-] anthropogenic heat output for cooling
       1526000.0_wp,   &  !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
       0.7_wp,         &  !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
       0.7_wp,         &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
       357200.0_wp,    &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
       0.04_wp,        &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
       1736000.0_wp,   &  !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
       0.19_wp,        &  !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
       1736000.0_wp,   &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
       0.19_wp,        &  !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
       0.7_wp,         &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
       1736000.0_wp,   &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
       0.19_wp         &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                        /)

    building_pars(:,6) = (/                                                                        &
       0.71_wp,        &  !< parameter 0   - [-] wall fraction above ground floor level
       0.29_wp,        &  !< parameter 1   - [-] window fraction above ground floor level
       0.0_wp,         &  !< parameter 2   - [-] green fraction above ground floor level
       0.0_wp,         &  !< parameter 3   - [-] green fraction roof above ground floor level
       1.5_wp,         &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
       1.5_wp,         &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
       1520000.0_wp,   &  !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (outside) above ground floor level
       79200.0_wp,     &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
       1344000.0_wp,   &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
       0.93_wp,        &  !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (outside) above ground floor level
       0.035_wp,       &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level       
       0.68_wp,        &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
       299.15_wp,      &  !< parameter 12  - [K] indoor target summer temperature
       293.15_wp,      &  !< parameter 13  - [K] indoor target winter temperature
       0.93_wp,        &  !< parameter 14  - [-] wall emissivity above ground floor level
       0.86_wp,        &  !< parameter 15  - [-] green emissivity above ground floor level
       0.8_wp,         &  !< parameter 16  - [-] window emissivity above ground floor level
       0.57_wp,        &  !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
       0.001_wp,       &  !< parameter 18  - [m] z0 roughness above ground floor level
       0.0001_wp,      &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
       2.7_wp,         &  !< parameter 20  - [m] ground floor level height
       0.71_wp,        &  !< parameter 21  - [-] wall fraction ground floor level
       0.29_wp,        &  !< parameter 22  - [-] window fraction ground floor level
       0.0_wp,         &  !< parameter 23  - [-] green fraction ground floor level
       0.0_wp,         &  !< parameter 24  - [-] green fraction roof ground floor level
       1.5_wp,         &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
       1520000.0_wp,   &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground floor level
       79200.0_wp,     &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
       1344000.0_wp,   &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (outside) ground floor level
       0.035_wp,       &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
       0.68_wp,        &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
       0.93_wp,        &  !< parameter 32  - [-] wall emissivity ground floor level
       0.8_wp,         &  !< parameter 33  - [-] window emissivity ground floor level
       0.86_wp,        &  !< parameter 34  - [-] green emissivity ground floor level
       0.57_wp,        &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
       0.001_wp,       &  !< parameter 36  - [m] z0 roughness ground floor level
       0.0001_wp,      &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
       36.0_wp,        &  !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
       38.0_wp,        &  !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
       0.02_wp,        &  !< parameter 41  - [m] 1st wall layer thickness above ground floor level
       0.22_wp,        &  !< parameter 42  - [m] 2nd wall layer thickness above ground floor level
       0.58_wp,        &  !< parameter 43  - [m] 3rd wall layer thickness above ground floor level
       0.6_wp,         &  !< parameter 44  - [m] 4th wall layer thickness above ground floor level
       20000.0_wp,     &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
       23.0_wp,        &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
       20000.0_wp,     &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
       20000.0_wp,     &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
       23.0_wp,        &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
       10.0_wp,        &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
       1.0_wp,         &  !< parameter 51  - [-] wall fraction ground plate
       0.20_wp,        &  !< parameter 52  - [m] 1st wall layer thickness ground plate
       0.32_wp,        &  !< parameter 53  - [m] 2nd wall layer thickness ground plate
       0.38_wp,        &  !< parameter 54  - [m] 3rd wall layer thickness ground plate
       0.41_wp,        &  !< parameter 55  - [m] 4th wall layer thickness ground plate
       2112000.0_wp,   &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (outside) ground plate
       79200.0_wp,     &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
       2112000.0_wp,   &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
       2.1_wp,         &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (oustide) ground plate
       0.05_wp,        &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
       2.1_wp,         &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
       0.02_wp,        &  !< parameter 62  - [m] 1st wall layer thickness ground floor level
       0.22_wp,        &  !< parameter 63  - [m] 2nd wall layer thickness ground floor level
       0.58_wp,        &  !< parameter 64  - [m] 3rd wall layer thickness ground floor level
       0.6_wp,         &  !< parameter 65  - [m] 4th wall layer thickness ground floor level
       36.0_wp,        &  !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
       0.03_wp,        &  !< parameter 67  - [m] 1st window layer thickness ground floor level
       0.06_wp,        &  !< parameter 68  - [m] 2nd window layer thickness ground floor level
       0.09_wp,        &  !< parameter 69  - [m] 3rd window layer thickness ground floor level
       0.012_wp,       &  !< parameter 70  - [m] 4th window layer thickness ground floor level
       1736000.0_wp,   &  !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
       1736000.0_wp,   &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
       1736000.0_wp,   &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
       0.11_wp,        &  !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
       0.11_wp,        &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
       0.11_wp,        &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
       38.0_wp,        &  !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
       5.0_wp,         &  !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
       0.03_wp,        &  !< parameter 79  - [m] 1st window layer thickness above ground floor level
       0.06_wp,        &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
       0.09_wp,        &  !< parameter 81  - [m] 3rd window layer thickness above ground floor level
       0.12_wp,        &  !< parameter 82  - [m] 4th window layer thickness above ground floor level
       1736000.0_wp,   &  !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
       1736000.0_wp,   &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
       1736000.0_wp,   &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
       0.11_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
       0.11_wp,        &  !< parameter 86  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
       0.11_wp,        &  !< parameter 87  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
       1.0_wp,         &  !< parameter 89  - [-] wall fraction roof
       0.02_wp,        &  !< parameter 90  - [m] 1st wall layer thickness roof
       0.06_wp,        &  !< parameter 91  - [m] 2nd wall layer thickness roof
       0.36_wp,        &  !< parameter 92  - [m] 3rd wall layer thickness roof
       0.38_wp,        &  !< parameter 93  - [m] 4th wall layer thickness roof
       3753600.0_wp,   &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
       709650.0_wp,    &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
       79200.0_wp,     &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
       0.52_wp,        &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (outside) roof
       0.12_wp,        &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
       0.035_wp,       &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
       0.93_wp,        &  !< parameter 100 - [-] wall emissivity roof
       42.0_wp,        &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 102 - [-] window fraction roof
       0.03_wp,        &  !< parameter 103 - [m] window 1st layer thickness roof
       0.06_wp,        &  !< parameter 104 - [m] window 2nd layer thickness roof
       0.09_wp,        &  !< parameter 105 - [m] window 3rd layer thickness roof
       0.12_wp,        &  !< parameter 106 - [m] window 4th layer thickness roof
       1736000.0_wp,   &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
       1736000.0_wp,   &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
       1736000.0_wp,   &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
       0.11_wp,        &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
       0.11_wp,        &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
       0.11_wp,        &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
       0.8_wp,         &  !< parameter 113 - [-] window emissivity roof
       0.57_wp,        &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
       38.0_wp,        &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
       0.86_wp,        &  !< parameter 116 - [-] green emissivity roof
       5.0_wp,         &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
       0.0_wp,         &  !< parameter 118 - [-] green type roof
       0.15_wp,        &  !< parameter 119 - [-] shading factor
       0.6_wp,         &  !< parameter 120 - [-] g-value windows
       0.8_wp,         &  !< parameter 121 - [W/(m2*K)] u-value windows
       1.0_wp,         &  !< parameter 122 - [1/h] basic airflow without occupancy of the room for - summer 1.0_wp, winter 0.2
       1.0_wp,         &  !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room for - summer 1.0_wp, winter 0.8
       0.8_wp,         &  !< parameter 124 - [-] heat recovery efficiency
       2.5_wp,         &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
       165000.0_wp,    &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heat storage
       4.5_wp,         &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
       40.0_wp,        &  !< parameter 128 - [W] maximal heating capacity
       -80.0_wp,       &  !< parameter 129 - [W] maximal cooling capacity
       7.0_wp,         &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
       3.0_wp,         &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
       2.7_wp,         &  !< parameter 132 - [m] storey height
       0.2_wp,         &  !< parameter 133 - [m] ceiling construction height
       -2.0_wp,        &  !< parameter 134 - [-] anthropogenic heat output for heating
       1.25_wp,        &  !< parameter 135 - [-] anthropogenic heat output for cooling
       1526000.0_wp,   &  !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (inside) above ground floor level
       0.7_wp,         &  !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 138 - [J/(m3*K)] capacity 4th wall layer (inside) ground floor level
       0.7_wp,         &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground floor level
       709650.0_wp,    &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (inside) ground plate
       0.12_wp,        &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (inside) ground plate
       1736000.0_wp,   &  !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
       0.11_wp,        &  !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
       1736000.0_wp,   &  !< parameter 144 - [J/(m3*K)] heat capacity 4th layer (inside) above ground floor level
       0.11_wp,        &  !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
       1526000.0_wp,   &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
       0.7_wp,         &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (inside) roof
       1736000.0_wp,   &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
       0.11_wp         &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                        /)

    building_pars(:,7) = (/                                                                        &
      1.0_wp,          &  !< parameter 0   - [-] wall fraction above ground floor level
      0.0_wp,          &  !< parameter 1   - [-] window fraction above ground floor level
      0.0_wp,          &  !< parameter 2   - [-] green fraction above ground floor level
      0.0_wp,          &  !< parameter 3   - [-] green fraction roof above ground floor level
      1.5_wp,          &  !< parameter 4   - [m2/m2] LAI (Leaf Area Index) roof
      1.5_wp,          &  !< parameter 5   - [m2/m2] LAI (Leaf Area Index) on wall above ground floor level
      1950400.0_wp,    &  !< parameter 6   - [J/(m3*K)] heat capacity 1st wall layer (upside) above ground floor level
      1848000.0_wp,    &  !< parameter 7   - [J/(m3*K)] heat capacity 2nd wall layer above ground floor level
      1848000.0_wp,    &  !< parameter 8   - [J/(m3*K)] heat capacity 3rd wall layer above ground floor level
      0.7_wp,          &  !< parameter 9   - [W/(m*K)] thermal conductivity 1st wall layer (upside) above ground floor level
      1.0_wp,          &  !< parameter 10  - [W/(m*K)] thermal conductivity 2nd wall layer above ground floor level
      1.0_wp,          &  !< parameter 11  - [W/(m*K)] thermal conductivity 3rd wall layer above ground floor level
      372.15_wp,       &  !< parameter 12  - [K] indoor target summer temperature
      293.15_wp,       &  !< parameter 13  - [K] indoor target winter temperature
      0.93_wp,         &  !< parameter 14  - [-] wall emissivity above ground floor level
      0.86_wp,         &  !< parameter 15  - [-] green emissivity above ground floor level
      0.8_wp,          &  !< parameter 16  - [-] window emissivity above ground floor level
      0.7_wp,          &  !< parameter 17  - [-] window transmissivity (not visual transmissivity) above ground floor level
      0.001_wp,        &  !< parameter 18  - [m] z0 roughness above ground floor level
      0.0001_wp,       &  !< parameter 19  - [m] z0h/z0g roughness heat/humidity above ground floor level
      4.0_wp,          &  !< parameter 20  - [m] ground floor level height
      1.0_wp,          &  !< parameter 21  - [-] wall fraction ground floor level
      0.0_wp,          &  !< parameter 22  - [-] window fraction ground floor level
      0.0_wp,          &  !< parameter 23  - [-] green fraction ground floor level
      0.0_wp,          &  !< parameter 24  - [-] green fraction roof ground floor level
      1.5_wp,          &  !< parameter 25  - [m2/m2] LAI (Leaf Area Index) on wall ground floor level
      1950400.0_wp,    &  !< parameter 26  - [J/(m3*K)] heat capacity 1st wall layer (upside) ground floor level
      1848000.0_wp,    &  !< parameter 27  - [J/(m3*K)] heat capacity 2nd wall layer ground floor level
      1848000.0_wp,    &  !< parameter 28  - [J/(m3*K)] heat capacity 3rd wall layer ground floor level
      0.7_wp,          &  !< parameter 29  - [W/(m*K)] thermal conductivity 1st wall layer (upside) ground floor level
      1.0_wp,          &  !< parameter 30  - [W/(m*K)] thermal conductivity 2nd wall layer ground floor level
      1.0_wp,          &  !< parameter 31  - [W/(m*K)] thermal conductivity 3rd wall layer ground floor level
      0.93_wp,         &  !< parameter 32  - [-] wall emissivity ground floor level
      0.8_wp,          &  !< parameter 33  - [-] window emissivity ground floor level
      0.86_wp,         &  !< parameter 34  - [-] green emissivity ground floor level
      0.7_wp,          &  !< parameter 35  - [-] window transmissivity (not visual transmissivity) ground floor level
      0.001_wp,        &  !< parameter 36  - [m] z0 roughness ground floor level
      0.0001_wp,       &  !< parameter 37  - [m] z0h/z0q roughness heat/humidity
      20.0_wp,         &  !< parameter 38  - [-] wall albedo_type above ground floor level  (albedo_type specified in radiation model)
      5.0_wp,          &  !< parameter 39  - [-] green albedo_type above ground floor level  (albedo_type specified in radiation model)
      37.0_wp,         &  !< parameter 40  - [-] window albedo_type above ground floor level  (albedo_type specified in radiation model)
      0.29_wp,         &  !< parameter 41  - [m] 1st wall layer thickness above ground floor level
      0.4_wp,          &  !< parameter 42  - [m] 2nd wall layer thickness above ground floor level
      0.695_wp,        &  !< parameter 43  - [m] 3rd wall layer thickness above ground floor level
      0.985_wp,        &  !< parameter 44  - [m] 4th wall layer thickness above ground floor level
      20000.0_wp,      &  !< parameter 45  - [J/(m2*K)] heat capacity wall surface (1 cm air)
      23.0_wp,         &  !< parameter 46  - [W/(m2*K)] thermal conductivity of wall surface (1 cm air)
      20000.0_wp,      &  !< parameter 47  - [J/(m2*K)] heat capacity of window surface (1 cm air)
      20000.0_wp,      &  !< parameter 48  - [J/(m2*K)] heat capacity of green surface
      23.0_wp,         &  !< parameter 49  - [W/(m2*K)] thermal conductivity of window surface (1 cm air)
      10.0_wp,         &  !< parameter 50  - [W/(m2*K)] thermal conductivty of green surface
      1.0_wp,          &  !< parameter 51  - [-] wall fraction ground plate
      0.29_wp,         &  !< parameter 52  - [m] 1st wall layer thickness ground plate
      0.4_wp,          &  !< parameter 53  - [m] 2nd wall layer thickness ground plate
      0.695_wp,        &  !< parameter 54  - [m] 3rd wall layer thickness ground plate
      0.985_wp,        &  !< parameter 55  - [m] 4th wall layer thickness ground plate
      1950400.0_wp,    &  !< parameter 56  - [J/(m3*K)] heat capacity 1st wall layer (upside) ground plate
      1848000.0_wp,    &  !< parameter 57  - [J/(m3*K)] heat capacity 2nd wall layer ground plate
      1848000.0_wp,    &  !< parameter 58  - [J/(m3*K)] heat capacity 3rd wall layer ground plate
      0.7_wp,          &  !< parameter 59  - [W/(m*K)] thermal conductivity 1st wall layer (upside) ground plate
      1.0_wp,          &  !< parameter 60  - [W/(m*K)] thermal conductivity 2nd wall layer ground plate
      1.0_wp,          &  !< parameter 61  - [W/(m*K)] thermal conductivity 3rd wall layer ground plate
      0.29_wp,         &  !< parameter 62  - [m] 1st wall layer thickness ground floor level
      0.4_wp,          &  !< parameter 63  - [m] 2nd wall layer thickness ground floor level
      0.695_wp,        &  !< parameter 64  - [m] 3rd wall layer thickness ground floor level
      0.985_wp,        &  !< parameter 65  - [m] 4th wall layer thickness ground floor level
      20.0_wp,         &  !< parameter 66  - [-] wall albedo_type ground floor level (albedo_type specified in radiation model)
      0.003_wp,        &  !< parameter 67  - [m] 1st window layer thickness ground floor level
      0.006_wp,        &  !< parameter 68  - [m] 2nd window layer thickness ground floor level
      0.012_wp,        &  !< parameter 69  - [m] 3rd window layer thickness ground floor level
      0.018_wp,        &  !< parameter 70  - [m] 4th window layer thickness ground floor level
      1736000.0_wp,    &  !< parameter 71  - [J/(m3*K)] heat capacity 1st window layer (outside) ground floor level
      1736000.0_wp,    &  !< parameter 72  - [J/(m3*K)] heat capacity 2nd window layer ground floor level
      1736000.0_wp,    &  !< parameter 73  - [J/(m3*K)] heat capacity 3rd window layer ground floor level
      0.57_wp,         &  !< parameter 74  - [W/(m*K)] thermal conductivity 1st window layer (outside) ground floor level
      0.57_wp,         &  !< parameter 75  - [W/(m*K)] thermal conductivity 2nd window layer ground floor level
      0.57_wp,         &  !< parameter 76  - [W/(m*K)] thermal conductivity 3rd window layer ground floor level
      37.0_wp,         &  !< parameter 77  - [-] window albedo_type ground floor level (albedo_type specified in radiation model)
      5.0_wp,          &  !< parameter 78  - [-] green albedo_type ground floor level (albedo_type specified in radiation model)
      0.003_wp,        &  !< parameter 79  - [m] 1st window layer thickness above ground floor level
      0.006_wp,        &  !< parameter 80  - [m] 2nd thickness window layer above ground floor level
      0.012_wp,        &  !< parameter 81  - [m] 3rd window layer thickness above ground floor level
      0.018_wp,        &  !< parameter 82  - [m] 4th window layer thickness above ground floor level
      1736000.0_wp,    &  !< parameter 83  - [J/(m3*K)] heat capacity 1st window layer (outside) above ground floor level
      1736000.0_wp,    &  !< parameter 84  - [J/(m3*K)] heat capacity 2nd window layer above ground floor level
      1736000.0_wp,    &  !< parameter 85  - [J/(m3*K)] heat capacity 3rd window layer above ground floor level
      0.57_wp,         &  !< parameter 86  - [W/(m*K)] thermal conductivity 1st window layer (outside) above ground floor level
      0.57_wp,         &  !< parameter 87  - [W/(m*K)] thermal conductivity 2nd window layer above ground floor level
      0.57_wp,         &  !< parameter 88  - [W/(m*K)] thermal conductivity 3rd window layer above ground floor level
      1.0_wp,          &  !< parameter 89  - [-] wall fraction roof
      0.29_wp,         &  !< parameter 90  - [m] 1st wall layer thickness roof
      0.4_wp,          &  !< parameter 91  - [m] 2nd wall layer thickness roof
      0.695_wp,        &  !< parameter 92  - [m] 3rd wall layer thickness roof
      0.985_wp,        &  !< parameter 93  - [m] 4th wall layer thickness roof
      1950400.0_wp,    &  !< parameter 94  - [J/(m3*K)] heat capacity 1st wall layer (outside) roof
      1848000.0_wp,    &  !< parameter 95  - [J/(m3*K)] heat capacity 2nd wall layer roof
      1848000.0_wp,    &  !< parameter 96  - [J/(m3*K)] heat capacity 3rd wall layer roof
      0.7_wp,          &  !< parameter 97  - [W/(m*K)] thermal conductivity 1st wall layer (upside) roof
      1.0_wp,          &  !< parameter 98  - [W/(m*K)] thermal conductivity 2nd wall layer roof
      1.0_wp,          &  !< parameter 99  - [W/(m*K)] thermal conductivity 3rd wall layer roof
      0.93_wp,         &  !< parameter 100 - [-] wall emissivity roof
      19.0_wp,         &  !< parameter 101 - [-] wall albedo_type roof (albedo_type specified in radiation model)
      0.0_wp,          &  !< parameter 102 - [-] window fraction roof
      0.003_wp,        &  !< parameter 103 - [m] window 1st layer thickness roof
      0.006_wp,        &  !< parameter 104 - [m] window 2nd layer thickness roof
      0.012_wp,        &  !< parameter 105 - [m] window 3rd layer thickness roof
      0.018_wp,        &  !< parameter 106 - [m] window 4th layer thickness roof
      1736000.0_wp,    &  !< parameter 107 - [J/(m3*K)] heat capacity 1st window layer (outside) roof
      1736000.0_wp,    &  !< parameter 108 - [J/(m3*K)] heat capacity 2nd window layer roof
      1736000.0_wp,    &  !< parameter 109 - [J/(m3*K)] heat capacity 3rd window layer roof
      0.57_wp,         &  !< parameter 110 - [W/(m*K)] thermal conductivity 1st window layer (outside) roof
      0.57_wp,         &  !< parameter 111 - [W/(m*K)] thermal conductivity 2nd window layer roof
      0.57_wp,         &  !< parameter 112 - [W/(m*K)] thermal conductivity 3rd window layer roof
      0.8_wp,          &  !< parameter 113 - [-] window emissivity roof
      0.7_wp,          &  !< parameter 114 - [-] window transmissivity (not visual transmissivity) roof
      37.0_wp,         &  !< parameter 115 - [-] window albedo_type roof (albedo_type specified in radiation model)
      0.86_wp,         &  !< parameter 116 - [-] green emissivity roof
      5.0_wp,          &  !< parameter 117 - [-] green albedo_type roof (albedo_type specified in radiation model)
      0.0_wp,          &  !< parameter 118 - [-] green type roof
      0.8_wp,          &  !< parameter 119 - [-] shading factor
      100.0_wp,        &  !< parameter 120 - [-] g-value windows
      100.0_wp,        &  !< parameter 121 - [W/(m2*K)] u-value windows
      20.0_wp,         &  !< parameter 122 - [1/h] basic airflow without occupancy of the room
      20.0_wp,         &  !< parameter 123 - [1/h] additional airflow dependent on occupancy of the room
      0.0_wp,          &  !< parameter 124 - [-] heat recovery efficiency
      1.0_wp,          &  !< parameter 125 - [m2/m2] dynamic parameter specific effective surface
      1.0_wp,          &  !< parameter 126 - [J/(m2*K)] dynamic parameter innner heatstorage
      4.5_wp,          &  !< parameter 127 - [m2/m2] ratio internal surface/floor area
      100000.0_wp,     &  !< parameter 128 - [W] maximal heating capacity
      0.0_wp,          &  !< parameter 129 - [W] maximal cooling capacity
      0.0_wp,          &  !< parameter 130 - [W/m2] additional internal heat gains dependent on occupancy of the room
      0.0_wp,          &  !< parameter 131 - [W/m2] basic internal heat gains without occupancy of the room
      3.0_wp,          &  !< parameter 132 - [m] storey height
      0.2_wp,          &  !< parameter 133 - [m] ceiling construction height
      0.0_wp,          &  !< parameter 134 - [-] anthropogenic heat output for heating
      0.0_wp,          &  !< parameter 135 - [-] anthropogenic heat output for cooling
      1848000.0_wp,    &  !< parameter 136 - [J/(m3*K)] heat capacity 4th wall layer (downside) above ground floor level
      1.0_wp,          &  !< parameter 137 - [W/(m*K)] thermal conductivity 4th wall layer (downside) above ground floor level
      1848000.0_wp,    &  !< parameter 138 - [J/(m3*K)] heat capacity 4th wall layer (downside) ground floor level 
      1.0_wp,          &  !< parameter 139 - [W/(m*K)] thermal conductivity 4th wall layer (downside) ground floor level
      1848000.0_wp,    &  !< parameter 140 - [J/(m3*K)] heat capacity 4th wall layer (downside) ground plate
      1.0_wp,          &  !< parameter 141 - [W/(m*K)] thermal conductivity 4th wall layer (downside) ground plate
      1736000.0_wp,    &  !< parameter 142 - [J/(m3*K)] heat capacity 4th window layer (inside) ground floor level
      0.57_wp,         &  !< parameter 143 - [W/(m*K)] thermal conductivity 4th window layer (inside) ground floor level
      1736000.0_wp,    &  !< parameter 144 - [J/(m3*K)] heat capacity 4th window layer (inside) above ground floor level
      0.57_wp,         &  !< parameter 145 - [W/(m*K)] thermal conductivity 4th window layer (inside) above ground floor level
      1848000.0_wp,    &  !< parameter 146 - [J/(m3*K)] heat capacity 4th wall layer (inside) roof
      1.0_wp,          &  !< parameter 147 - [W/(m*K)] thermal conductivity 4th wall layer (downside) roof
      1736000.0_wp,    &  !< parameter 148 - [J/(m3*K)] heat capacity 4th window layer (inside) roof
      0.57_wp          &  !< parameter 149 - [W/(m*K)] thermal conductivity 4th window layer (inside) roof
                        /)

 END SUBROUTINE usm_define_pars


 END MODULE urban_surface_mod
