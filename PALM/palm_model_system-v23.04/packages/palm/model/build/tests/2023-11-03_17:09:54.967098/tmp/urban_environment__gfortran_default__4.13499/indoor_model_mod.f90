!> @file indoor_model_mod.f90
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
! Copyright 2018-2021 Leibniz Universitaet Hannover
! Copyright 2018-2021 Hochschule Offenburg
!--------------------------------------------------------------------------------------------------!
!
! Authors:
! --------
! @author Tobias Lang
! @author Jens Pfafferott
! @author Farah Kanani-Suehring
! @author Matthias Suehring
! @author Sascha Rissmann
! @author Bjoern Maronga
!
!
! Description:
! ------------
!> Module for Indoor Climate Model (ICM)
!> The module is based on the DIN EN ISO 13790 with simplified hour-based procedure.
!> This model is a equivalent circuit diagram of a three-point RC-model (5R1C).
!> This module differs between indoor-air temperature an average temperature of indoor surfaces which make it prossible to determine
!> thermal comfort
!> the heat transfer between indoor and outdoor is simplified

!> @todo Many statement comments that are given as doxygen comments need to be changed to normal comments
!> @todo Replace window_area_per_facade by %frac(1,m) for window
!> @todo emissivity change for window blinds if solar_protection_on=1

!> @note Do we allow use of integer flags, or only logical flags? (concerns e.g. cooling_on, heating_on)
!> @note How to write indoor temperature output to pt array?
!>
!> @bug  Calculation of iwghf_eb and iwghf_eb_window is faulty
!--------------------------------------------------------------------------------------------------!
 MODULE indoor_model_mod

#if defined( __parallel )
    USE MPI
#endif

    USE arrays_3d,                                                                                 &
        ONLY:  dzw,                                                                                &
               pt

    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  kappa

    USE control_parameters,                                                                        &
        ONLY:  initializing_actions,                                                               &
               restart_data_format_output

    USE grid_variables,                                                                            &
        ONLY:  dx,                                                                                 &
               dy

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys

    USE kinds

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  building_id_f,                                                                      &
               building_type_f

    USE palm_date_time_mod,                                                                        &
        ONLY:  get_date_time,                                                                      &
               northward_equinox,                                                                  &
               seconds_per_hour,                                                                   &
               southward_equinox

    USE pegrid

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array,                                                              &
               rd_mpi_io_surface_filetypes,                                                        &
               rrd_mpi_io,                                                                         &
               rrd_mpi_io_surface,                                                                 &
               wrd_mpi_io,                                                                         &
               wrd_mpi_io_surface

    USE surface_layer_fluxes_mod,                                                                  &
        ONLY:  psi_h

    USE surface_mod,                                                                               &
        ONLY:  ind_pav_green,                                                                      &
               ind_veg_wall,                                                                       &
               ind_wat_win,                                                                        &
               surf_usm,                                                                           &
               surface_restore_elements

    IMPLICIT NONE

!
!-- Define data structure for buidlings.
    TYPE build

       INTEGER(iwp) ::  id                                !< building ID
       INTEGER(iwp) ::  kb_max                            !< highest vertical index of a building
       INTEGER(iwp) ::  kb_min                            !< lowest vertical index of a building
       INTEGER(iwp) ::  num_facades_per_building = 0      !< total number of horizontal and vertical facade elements
       INTEGER(iwp) ::  num_facades_per_building_h = 0    !< total number of horizontal facades elements
       INTEGER(iwp) ::  num_facades_per_building_l = 0    !< number of horizontal and vertical facade elements on local subdomain
       INTEGER(iwp) ::  num_facades_per_building_v = 0    !< total number of vertical facades elements
       INTEGER(iwp) ::  ventilation_int_loads             !< [-] allocation of activity in the building

       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  m              !< index array linking surface-element index with building
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_facade_h   !< number of horizontal facade elements per buidling
                                                                  !< and height level
       INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_facade_v   !< number of vertical facades elements per buidling
                                                                  !< and height level


       LOGICAL ::  on_pe = .FALSE.   !< flag indicating whether a building with certain ID is on local subdomain

       REAL(wp) ::  air_change_high       !< [1/h] air changes per time_utc_hour
       REAL(wp) ::  air_change_low        !< [1/h] air changes per time_utc_hour
       REAL(wp) ::  area_facade           !< [m2] area of total facade
       REAL(wp) ::  building_height       !< building height
       REAL(wp) ::  eta_ve                !< [-] heat recovery efficiency
       REAL(wp) ::  factor_a              !< [-] Dynamic parameters specific effective surface according to Table 12; 2.5
                                          !< (very light, light and medium), 3.0 (heavy), 3.5 (very heavy)
       REAL(wp) ::  factor_c              !< [J/(m2 K)] Dynamic parameters inner heatstorage according to Table 12; 80000
                                          !< (very light), 110000 (light), 165000 (medium), 260000 (heavy), 370000 (very heavy)
       REAL(wp) ::  f_c_win               !< [-] shading factor
       REAL(wp) ::  fapf                  !< [m2/m2] floor area per facade
       REAL(wp) ::  g_value_win           !< [-] SHGC factor
       REAL(wp) ::  h_es                  !< [W/(m2 K)] surface-related heat transfer coefficient between extern and surface
       REAL(wp) ::  height_cei_con        !< [m] ceiling construction heigth
       REAL(wp) ::  height_storey         !< [m] storey heigth
       REAL(wp) ::  params_waste_heat_c   !< [-] anthropogenic heat outputs for cooling e.g. 1.33 for KKM with COP = 3
       REAL(wp) ::  params_waste_heat_h   !< [-] anthropogenic heat outputs for heating e.g. 1 - 0.9 = 0.1 for combustion with
                                          !< eta = 0.9 or -2 for WP with COP = 3
       REAL(wp) ::  phi_c_max             !< [W] Max. Cooling capacity (negative)
       REAL(wp) ::  phi_h_max             !< [W] Max. Heating capacity (positive)
       REAL(wp) ::  q_c_max               !< [W/m2] Max. Cooling heat flux per netto floor area (negative)
       REAL(wp) ::  q_h_max               !< [W/m2] Max. Heating heat flux per netto floor area (positive)
       REAL(wp) ::  qint_high             !< [W/m2] internal heat gains, option Database qint_0-23
       REAL(wp) ::  qint_low              !< [W/m2] internal heat gains, option Database qint_0-23
       REAL(wp) ::  lambda_at             !< [-] ratio internal surface/floor area chap. 7.2.2.2.
       REAL(wp) ::  lambda_layer3         !< [W/(m*K)] Thermal conductivity of the inner layer
       REAL(wp) ::  net_floor_area        !< [m2] netto ground area
       REAL(wp) ::  s_layer3              !< [m] half thickness of the inner layer (layer_3)
       REAL(wp) ::  theta_int_c_set       !< [degree_C] Max. Setpoint temperature (summer)
       REAL(wp) ::  theta_int_h_set       !< [degree_C] Max. Setpoint temperature (winter)
       REAL(wp) ::  u_value_win           !< [W/(m2*K)] transmittance
       REAL(wp) ::  vol_tot               !< [m3] total building volume

       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_in           !< mean building indoor temperature, height dependent
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_in_l         !< mean building indoor temperature on local subdomain, height dependent
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  theta_m_t_prev !< [degree_C] value of theta_m_t from previous time step
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  volume         !< total building volume, height dependent
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vol_frac       !< fraction of local on total building volume, height dependent
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  vpf            !< building volume volume per facade element, height dependent

    END TYPE build

    TYPE(build), DIMENSION(:), ALLOCATABLE ::  buildings   !< building array

    INTEGER(iwp) ::  cooling_on              !< Indoor cooling flag (0=off, 1=on)
    INTEGER(iwp) ::  heating_on              !< Indoor heating flag (0=off, 1=on)
    INTEGER(iwp) ::  num_build               !< total number of buildings in domain
    INTEGER(iwp) ::  solar_protection_off    !< Solar protection off
    INTEGER(iwp) ::  solar_protection_on     !< Solar protection on

    LOGICAL ::  indoor_during_spinup = .FALSE.      !< namelist parameter used to switch-off/on the indoor model during spinup
!
!-- Declare all global variables within the module

    REAL(wp), PARAMETER ::  dt_indoor = 3600.0_wp                  !< [s] time interval for indoor-model application, fixed to
                                                                   !< 3600.0 due to model requirements
    REAL(wp), PARAMETER ::  h_is                     = 3.45_wp     !< [W/(m2 K)] surface-related heat transfer coefficient between
                                                                   !< surface and air (chap. 7.2.2.2)
    REAL(wp), PARAMETER ::  h_ms                     = 9.1_wp      !< [W/(m2 K)] surface-related heat transfer coefficient between
                                                                   !< component and surface (chap. 12.2.2)
    REAL(wp), PARAMETER ::  params_f_f               = 0.3_wp      !< [-] frame ratio chap. 8.3.2.1.1 for buildings with mostly
                                                                   !< cooling 2.0_wp
    REAL(wp), PARAMETER ::  params_f_w               = 0.9_wp      !< [-] correction factor (fuer nicht senkrechten Stahlungseinfall
                                                                   !< DIN 4108-2 chap.8, (hier konstant, keine Winkelabhängigkeit)
    REAL(wp), PARAMETER ::  params_f_win             = 0.5_wp      !< [-] proportion of window area, Database A_win aus
                                                                   !< Datenbank 27 window_area_per_facade_percent
    REAL(wp), PARAMETER ::  params_solar_protection  = 300.0_wp    !< [W/m2] chap. G.5.3.1 sun protection closed, if the radiation
                                                                   !< on facade exceeds this value

    REAL(wp) ::  a_m                                 !< [m2] the effective mass-related area
    REAL(wp) ::  air_change                          !< [1/h] Airflow
    REAL(wp) ::  c_m                                 !< [J/K] internal heat storage capacity
    REAL(wp) ::  facade_element_area                 !< [m2_facade] building surface facade
    REAL(wp) ::  floor_area_per_facade               !< [m2/m2] floor area per facade area
    REAL(wp) ::  h_t_1                               !< [W/K] Heat transfer coefficient auxiliary variable 1
    REAL(wp) ::  h_t_2                               !< [W/K] Heat transfer coefficient auxiliary variable 2
    REAL(wp) ::  h_t_3                               !< [W/K] Heat transfer coefficient auxiliary variable 3
    REAL(wp) ::  h_t_es                              !< [W/K] heat transfer coefficient of doors, windows, curtain walls and
                                                     !< glazed walls (assumption: thermal mass=0)
    REAL(wp) ::  h_t_is                              !< [W/K] thermal coupling conductance (Thermischer Kopplungsleitwert)
    REAL(wp) ::  h_t_ms                              !< [W/K] Heat transfer conductance term (got with h_t_wm the thermal mass)
    REAL(wp) ::  h_t_wall                            !< [W/K] heat transfer coefficient of opaque components (assumption: got all
                                                     !< thermal mass) contains of h_t_wm and h_t_ms
    REAL(wp) ::  h_t_wm                              !< [W/K] Heat transfer coefficient of the emmision (got with h_t_ms the
                                                     !< thermal mass)
    REAL(wp) ::  h_v                                 !< [W/K] heat transfer of ventilation
    REAL(wp) ::  indoor_volume_per_facade            !< [m3] indoor air volume per facade element
    REAL(wp) ::  initial_indoor_temperature = 293.15 !< [K] initial indoor temperature (namelist parameter)
    REAL(wp) ::  net_sw_in                           !< [W/m2] net short-wave radiation
    REAL(wp) ::  phi_hc_nd                           !< [W] heating demand and/or cooling demand
    REAL(wp) ::  phi_hc_nd_10                        !< [W] heating demand and/or cooling demand for heating or cooling
    REAL(wp) ::  phi_hc_nd_ac                        !< [W] actual heating demand and/or cooling demand
    REAL(wp) ::  phi_hc_nd_un                        !< [W] unlimited heating demand and/or cooling demand which is necessary to
                                                     !< reach the demanded required temperature (heating is positive,
                                                     !< cooling is negative)
    REAL(wp) ::  phi_ia                              !< [W] internal air load = internal loads * 0.5, Eq. (C.1)
    REAL(wp) ::  phi_m                               !< [W] mass specific thermal load (internal and external)
    REAL(wp) ::  phi_mtot                            !< [W] total mass specific thermal load (internal and external)
    REAL(wp) ::  phi_sol                             !< [W] solar loads
    REAL(wp) ::  phi_st                              !< [W] mass specific thermal load implied non thermal mass
    REAL(wp) ::  q_wall                              !< [W/m2]heat flux from indoor into wall
    REAL(wp) ::  q_win                               !< [W/m2]heat flux from indoor into window
    REAL(wp) ::  q_waste_heat                        !< [W/m2]waste heat, sum of waste heat over the roof to Palm

    REAL(wp) ::  q_c_m                               !< [W] Energy of thermal storage mass specific thermal load for internal
                                                     !< and external heatsources (for energy bilanz)
    REAL(wp) ::  q_c_st                              !< [W] Energy of thermal storage mass specific thermal load implied non
                                                     !< thermal mass (for energy bilanz)
    REAL(wp) ::  q_int                               !< [W] Energy of internal air load (for energy bilanz)
    REAL(wp) ::  q_sol                               !< [W] Energy of solar (for energy bilanz)
    REAL(wp) ::  q_vent                              !< [W] Energy of ventilation (for energy bilanz)

    REAL(wp) ::  schedule_d                          !< [-] activation for internal loads (low or high + low)
    REAL(wp) ::  skip_time_do_indoor = 0.0_wp        !< [s] Indoor model is not called before this time
    REAL(wp) ::  theta_air                           !< [degree_C] air temperature of the RC-node
    REAL(wp) ::  theta_air_0                         !< [degree_C] air temperature of the RC-node in equilibrium
    REAL(wp) ::  theta_air_10                        !< [degree_C] air temperature of the RC-node from a heating capacity
                                                     !< of 10 W/m2
    REAL(wp) ::  theta_air_ac                        !< [degree_C] actual room temperature after heating/cooling
    REAL(wp) ::  theta_air_set                       !< [degree_C] Setpoint_temperature for the room
    REAL(wp) ::  theta_m                             !< [degree_C} inner temperature of the RC-node
    REAL(wp) ::  theta_m_t                           !< [degree_C] (Fictive) component temperature during timestep
    REAL(wp) ::  theta_op                            !< [degree_C] operative temperature
    REAL(wp) ::  theta_s                             !< [degree_C] surface temperature of the RC-node
    REAL(wp) ::  time_indoor = 0.0_wp                !< [s] time since last call of indoor model
    REAL(wp) ::  total_area                          !< [m2] area of all surfaces pointing to zone
    REAL(wp) ::  window_area_per_facade              !< [m2] window area per facade element

!
!-- Definition of seasonal parameters, summer and winter, for different building types
    REAL(wp), DIMENSION(0:1,1:7) ::  summer_pars = RESHAPE( (/                & ! building_type 1
                                          0.5_wp,                             & ! basical airflow without occupancy of the room
                                          1.5_wp,                             & ! additional airflow depend of occupancy of the room
                                          0.5_wp,                             & ! building_type 2: basical airflow without occupancy
                                                                                ! of the room
                                          1.5_wp,                             & ! additional airflow depend of occupancy of the room
                                          0.5_wp,                             & ! building_type 3: basical airflow without occupancy
                                                                                ! of the room
                                          1.5_wp,                             & ! additional airflow depend of occupancy of the room
                                          1.0_wp,                             & ! building_type 4: basical airflow without occupancy
                                                                                ! of the room
                                          1.0_wp,                             & ! additional airflow depend of occupancy of the room
                                          1.0_wp,                             & ! building_type 5: basical airflow without occupancy
                                                                                ! of the room
                                          1.0_wp,                             & ! additional airflow depend of occupancy of the room
                                          1.0_wp,                             & ! building_type 6: basical airflow without occupancy
                                                                                ! of the room
                                          1.0_wp,                             & ! additional airflow depend of occupancy of the room
                                          1.0_wp,                             & ! building_type 7: basical airflow without occupancy
                                                                                ! of the room
                                          1.0_wp                              & ! additional airflow depend of occupancy of the room
                                                           /), (/ 2, 7 /) )

    REAL(wp), DIMENSION(0:1,1:7) ::  winter_pars = RESHAPE( (/                & ! building_type 1
                                          0.5_wp,                             & ! basical airflow without occupancy of the room
                                          0.0_wp,                             & ! additional airflow depend of occupancy of the room
                                          0.5_wp,                             & ! building_type 2: basical airflow without occupancy
                                                                                ! of the room
                                          0.0_wp,                             & ! additional airflow depend of occupancy of the room
                                          0.5_wp,                             & ! building_type 3: basical airflow without occupancy
                                                                                ! of the room
                                          0.0_wp,                             & ! additional airflow depend of occupancy of the room
                                          0.2_wp,                             & ! building_type 4: basical airflow without occupancy
                                                                                ! of the room
                                          0.8_wp,                             & ! additional airflow depend of occupancy of the room
                                          0.2_wp,                             & ! building_type 5: basical airflow without occupancy
                                                                                ! of the room
                                          0.8_wp,                             & ! additional airflow depend of occupancy of the room
                                          0.2_wp,                             & ! building_type 6: basical airflow without occupancy
                                                                                ! of the room
                                          0.8_wp,                             & ! additional airflow depend of occupancy of the room
                                          0.2_wp,                             & ! building_type 7: basical airflow without occupancy
                                                                                ! of the room
                                          0.8_wp                              & ! additional airflow depend of occupancy of the room
                                                           /), (/ 2, 7 /) )

    TYPE surf_type_1d_im
       REAL(wp), DIMENSION(:), ALLOCATABLE ::  val  !<
    END TYPE surf_type_1d_im

    SAVE


    PRIVATE

!
!-- Add INTERFACES that must be available to other modules
    PUBLIC im_define_netcdf_grid,                                                                  &
           im_check_data_output,                                                                   &
           im_check_parameters,                                                                    &
           im_data_output_3d,                                                                      &
           im_init,                                                                                &
           im_init_arrays,                                                                         &
           im_main_heatcool,                                                                       &
           im_parin,                                                                               &
           im_rrd_local,                                                                           &
           im_wrd_local



!
!-- Add VARIABLES that must be available to other modules
    PUBLIC dt_indoor,                                                                              &
           indoor_during_spinup,                                                                   &
           skip_time_do_indoor,                                                                    &
           time_indoor

!
!-- PALM interfaces:
!-- Data output checks for 2D/3D data to be done in check_parameters
     INTERFACE im_check_data_output
        MODULE PROCEDURE im_check_data_output
     END INTERFACE im_check_data_output
!
!-- Input parameter checks to be done in check_parameters
     INTERFACE im_check_parameters
        MODULE PROCEDURE im_check_parameters
     END INTERFACE im_check_parameters
!
!-- Data output of 3D data
     INTERFACE im_data_output_3d
        MODULE PROCEDURE im_data_output_3d
     END INTERFACE im_data_output_3d

!
!-- Definition of data output quantities
     INTERFACE im_define_netcdf_grid
        MODULE PROCEDURE im_define_netcdf_grid
     END INTERFACE im_define_netcdf_grid
!
! !
! !-- Output of information to the header file
!     INTERFACE im_header
!        MODULE PROCEDURE im_header
!     END INTERFACE im_header
!
!-- Initialization actions
    INTERFACE im_init
       MODULE PROCEDURE im_init
    END INTERFACE im_init
!
!-- Array allocation
    INTERFACE im_init_arrays
       MODULE PROCEDURE im_init_arrays
    END INTERFACE im_init_arrays
!
!-- Main part of indoor model
    INTERFACE im_main_heatcool
       MODULE PROCEDURE im_main_heatcool
    END INTERFACE im_main_heatcool
!
!-- Reading of NAMELIST parameters
    INTERFACE im_parin
       MODULE PROCEDURE im_parin
    END INTERFACE im_parin

    INTERFACE im_rrd_local
       MODULE PROCEDURE im_rrd_local_ftn
       MODULE PROCEDURE im_rrd_local_mpi
    END INTERFACE im_rrd_local

    INTERFACE im_wrd_local
       MODULE PROCEDURE im_wrd_local
    END INTERFACE im_wrd_local

 CONTAINS

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!< Calculation of the air temperatures and mean radiation temperature.
!< This is basis for the operative temperature.
!< Based on a Crank-Nicholson scheme with a timestep of a hour.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_calc_temperatures ( i, j, k, indoor_wall_temperature,                               &
                                   near_facade_temperature, phi_hc_nd_dummy, theta_m_t_prev )

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j
    INTEGER(iwp) ::  k

    REAL(wp) ::  indoor_wall_temperature   !< temperature of innermost wall layer evtl in im_calc_temperatures einfügen
    REAL(wp) ::  near_facade_temperature
    REAL(wp) ::  phi_hc_nd_dummy
    REAL(wp), INTENT(IN) :: theta_m_t_prev
!
!-- Calculation of total mass specific thermal load (internal and external)   [degree_C] Eq. (C.5)
    phi_mtot = phi_m + h_t_wm * indoor_wall_temperature                                          &
                     + h_t_3  * ( phi_st + h_t_es * pt(k,j,i)                                    &
                                + h_t_1 * ( ( ( phi_ia + phi_hc_nd_dummy ) / h_v )               &
                                            + near_facade_temperature                            &
                                          )                                                      &
                                ) / h_t_2
!
!-- Calculation of component temperature at current timestep     [degree_C] Eq. (C.4)
    theta_m_t = ( theta_m_t_prev  * ( ( c_m / 3600.0_wp ) - 0.5_wp * ( h_t_3 + h_t_wm ) )          &
                  + phi_mtot                                                                       &
                ) / ( ( c_m / 3600.0_wp ) + 0.5_wp * ( h_t_3 + h_t_wm ) )

!
!-- Calculation of mean inner temperature for the RC-node in current timestep   [degree_C] Eq. (C.9)
    theta_m = ( theta_m_t + theta_m_t_prev ) * 0.5_wp

!
!-- Calculation of mean surface temperature of the RC-node in current timestep [degree_C] Eq. (C.10)
    theta_s = ( h_t_ms * theta_m + phi_st + h_t_es * pt(k,j,i)                                     &
                + h_t_1  * ( near_facade_temperature + ( phi_ia + phi_hc_nd_dummy ) / h_v )        &
              ) / ( h_t_ms + h_t_es + h_t_1 )


!
!-- Calculation of the air temperature of the RC-node     [degree_C] Eq. (C.11)
    theta_air = ( h_t_is * theta_s + h_v * near_facade_temperature + phi_ia + phi_hc_nd_dummy ) /  &
                ( h_t_is + h_v )

 END SUBROUTINE im_calc_temperatures


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the indoor model.
!> Static information are calculated here, e.g. building parameters and geometrical information,
!> anything that doesn't change in time.
!
!-- Input values
!-- Input datas from Palm, M4
!     i_global             -->  net_sw_in                         !< global radiation [W/m2]
!     theta_e              -->  pt(k,j,i)                         !< undisturbed outside temperature, 1. PALM volume, for windows
!     theta_sup = theta_f  -->  surf_usm%pt_10cm(m)
!                               surf_usm%pt_10cm(m)               !< Air temperature, facade near (10cm) air temperature from
                                                                  !< 1. Palm volume
!     theta_node           -->  t_wall(nzt_wall,m)
!                               t_wall%t(nzt_wall,m)              !< Temperature of innermost wall layer, for opaque wall
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_init

    USE control_parameters,                                                                        &
        ONLY:  message_string,                                                                     &
               time_since_reference_point

    USE indices,                                                                                   &
        ONLY:  nxl,                                                                                &
               nxr,                                                                                &
               nyn,                                                                                &
               nys,                                                                                &
               nzb,                                                                                &
               nzt,                                                                                &
               topo_flags

    USE urban_surface_mod,                                                                         &
        ONLY:  building_pars,                                                                      &
               building_type

    INTEGER(iwp) ::  bt          !< local building type
    INTEGER(iwp) ::  day_of_year !< day of the year
    INTEGER(iwp) ::  fa          !< running index for facade elements of each building
    INTEGER(iwp) ::  i           !< running index along x-direction
    INTEGER(iwp) ::  j           !< running index along y-direction
    INTEGER(iwp) ::  k           !< running index along z-direction
    INTEGER(iwp) ::  m           !< running index surface elements
    INTEGER(iwp) ::  n           !< building index
    INTEGER(iwp) ::  nb          !< building index

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids           !< building IDs on entire model domain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_final     !< building IDs on entire model domain,
                                                                    !< multiple occurences are sorted out
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_final_tmp !< temporary array used for resizing
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_l         !< building IDs on local subdomain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  build_ids_l_tmp     !< temporary array used to resize array of building IDs
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  displace_dum        !< displacements of start addresses, used for MPI_ALLGATHERV
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_max_l             !< highest vertical index of a building on subdomain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  k_min_l             !< lowest vertical index of a building on subdomain
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  n_fa                !< counting array
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_facades_h       !< dummy array used for summing-up total number of
                                                                    !< horizontal facade elements
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  num_facades_v       !< dummy array used for summing-up total number of
                                                                    !< vertical facade elements
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  receive_dum_h       !< dummy array used for MPI_ALLREDUCE
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  receive_dum_v       !< dummy array used for MPI_ALLREDUCE

    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_buildings         !< number of buildings with different ID on entire model domain
    INTEGER(iwp), DIMENSION(0:numprocs-1) ::  num_buildings_l       !< number of buildings with different ID on local subdomain

    REAL(wp) ::  u_tmp                                     !< dummy for temporary calculation of u-value without h_is
    REAL(wp) ::  du_tmp                                    !< 1/u_tmp
    REAL(wp) ::  du_win_tmp                                !< 1/building(nb)%u_value_win
    REAL(wp) ::  facade_area_v                             !< dummy to compute the total facade area from vertical walls

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  volume         !< total building volume at each discrete height level
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  volume_l       !< total building volume at each discrete height level,
                                                           !< on local subdomain

    CALL location_message( 'initializing indoor model', 'start' )
!
!-- Initializing of indoor model is only possible if buildings can be distinguished by their IDs.
    IF ( .NOT. building_id_f%from_file )  THEN
       message_string = 'indoor model requires information about building_id'
       CALL message( 'im_init', 'IDM0001', 1, 2, 0, 6, 0  )
    ENDIF
!
!-- Determine number of different building IDs on local subdomain.
    num_buildings_l = 0
    num_buildings   = 0
    ALLOCATE( build_ids_l(1) )
    DO  i = nxl, nxr
       DO  j = nys, nyn
          IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
             IF ( num_buildings_l(myid) > 0 )  THEN
                IF ( ANY( building_id_f%var(j,i) == build_ids_l ) )  THEN
                   CYCLE
                ELSE
                   num_buildings_l(myid) = num_buildings_l(myid) + 1
!
!--                Resize array with different local building ids
                   ALLOCATE( build_ids_l_tmp(1:SIZE(build_ids_l)) )
                   build_ids_l_tmp = build_ids_l
                   DEALLOCATE( build_ids_l )
                   ALLOCATE( build_ids_l(1:num_buildings_l(myid)) )
                   build_ids_l(1:num_buildings_l(myid)-1) =                                        &
                                                          build_ids_l_tmp(1:num_buildings_l(myid)-1)
                   build_ids_l(num_buildings_l(myid)) = building_id_f%var(j,i)
                   DEALLOCATE( build_ids_l_tmp )
                ENDIF
!
!--          First occuring building id on PE
             ELSE
                num_buildings_l(myid) = num_buildings_l(myid) + 1
                build_ids_l(1) = building_id_f%var(j,i)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
!
!-- Determine number of building IDs for the entire domain. (Note, building IDs can appear multiple
!-- times as buildings might be distributed over several PEs.)
#if defined( __parallel )
    CALL MPI_ALLREDUCE( num_buildings_l, num_buildings, numprocs, MPI_INTEGER, MPI_SUM, comm2d,    &
                        ierr )
#else
    num_buildings = num_buildings_l
#endif
    ALLOCATE( build_ids(1:SUM(num_buildings)) )
!
!-- Gather building IDs. Therefore, first, determine displacements used required for MPI_GATHERV
!-- call.
    ALLOCATE( displace_dum(0:numprocs-1) )
    displace_dum(0) = 0
    DO i = 1, numprocs-1
       displace_dum(i) = displace_dum(i-1) + num_buildings(i-1)
    ENDDO

#if defined( __parallel )
    CALL MPI_ALLGATHERV( build_ids_l(1:num_buildings_l(myid)),                                     &
                         num_buildings(myid),                                                      &
                         MPI_INTEGER,                                                              &
                         build_ids,                                                                &
                         num_buildings,                                                            &
                         displace_dum,                                                             &
                         MPI_INTEGER,                                                              &
                         comm2d, ierr )

    DEALLOCATE( displace_dum )

#else
    build_ids = build_ids_l
#endif
!
!-- Note: in parallel mode, building IDs can occur mutliple times, as each PE has send its own ids.
!-- Therefore, sort out building IDs which appear multiple times.
    num_build = 0
    DO  n = 1, SIZE(build_ids)

       IF ( ALLOCATED(build_ids_final) )  THEN
          IF ( ANY( build_ids(n) == build_ids_final ) )  THEN
             CYCLE
          ELSE
             num_build = num_build + 1
!
!--          Resize
             ALLOCATE( build_ids_final_tmp(1:num_build) )
             build_ids_final_tmp(1:num_build-1) = build_ids_final(1:num_build-1)
             DEALLOCATE( build_ids_final )
             ALLOCATE( build_ids_final(1:num_build) )
             build_ids_final(1:num_build-1) = build_ids_final_tmp(1:num_build-1)
             build_ids_final(num_build) = build_ids(n)
             DEALLOCATE( build_ids_final_tmp )
          ENDIF
       ELSE
          num_build = num_build + 1
          ALLOCATE( build_ids_final(1:num_build) )
          build_ids_final(num_build) = build_ids(n)
       ENDIF
    ENDDO

!
!-- Allocate building-data structure array. Note, this is a global array and all building IDs on
!-- domain are known by each PE. Further attributes, e.g. height-dependent arrays, however, are only
!-- allocated on PEs where  the respective building is present (in order to reduce memory demands).
    ALLOCATE( buildings(1:num_build) )

!
!-- Store building IDs and check if building with certain ID is present on subdomain.
    DO  nb = 1, num_build
       buildings(nb)%id = build_ids_final(nb)

       IF ( ANY( building_id_f%var(nys:nyn,nxl:nxr) == buildings(nb)%id ) )                        &
          buildings(nb)%on_pe = .TRUE.
    ENDDO
!
!-- Determine the maximum vertical dimension occupied by each building.
    ALLOCATE( k_min_l(1:num_build) )
    ALLOCATE( k_max_l(1:num_build) )
    k_min_l = nzt + 1
    k_max_l = 0

    DO  i = nxl, nxr
       DO  j = nys, nyn
          IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
             nb = MINLOC( ABS( buildings(:)%id - building_id_f%var(j,i) ), DIM=1 )
             DO  k = nzb, nzt+1
!
!--             Check if grid point belongs to a building.
                IF ( BTEST( topo_flags(k,j,i), 6 ) )  THEN
                   k_min_l(nb) = MIN( k_min_l(nb), k )
                   k_max_l(nb) = MAX( k_max_l(nb), k )
                ENDIF

             ENDDO
          ENDIF
       ENDDO
    ENDDO

#if defined( __parallel )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, k_min_l, num_build, MPI_INTEGER, MPI_MIN, comm2d, ierr )
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, k_max_l, num_build, MPI_INTEGER, MPI_MAX, comm2d, ierr )
#endif
    buildings(:)%kb_min = k_min_l(:)
    buildings(:)%kb_max = k_max_l(:)

    DEALLOCATE( k_min_l )
    DEALLOCATE( k_max_l )
!
!-- Calculate building height.
    DO  nb = 1, num_build
       buildings(nb)%building_height = 0.0_wp
       DO  k = buildings(nb)%kb_min, buildings(nb)%kb_max
          buildings(nb)%building_height = buildings(nb)%building_height + dzw(k+1)
       ENDDO
    ENDDO
!
!-- Calculate building volume
    DO  nb = 1, num_build
!
!--    Allocate temporary array for summing-up building volume
       ALLOCATE( volume(buildings(nb)%kb_min:buildings(nb)%kb_max)   )
       ALLOCATE( volume_l(buildings(nb)%kb_min:buildings(nb)%kb_max) )
       volume   = 0.0_wp
       volume_l = 0.0_wp
!
!--    Calculate building volume per height level on each PE where these building is present.
       IF ( buildings(nb)%on_pe )  THEN

          ALLOCATE( buildings(nb)%volume(buildings(nb)%kb_min:buildings(nb)%kb_max)   )
          ALLOCATE( buildings(nb)%vol_frac(buildings(nb)%kb_min:buildings(nb)%kb_max) )
          buildings(nb)%volume   = 0.0_wp
          buildings(nb)%vol_frac = 0.0_wp

          IF ( ANY( building_id_f%var(nys:nyn,nxl:nxr) == buildings(nb)%id ) )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = buildings(nb)%kb_min, buildings(nb)%kb_max
                      IF ( building_id_f%var(j,i) /= building_id_f%fill )                          &
                         volume_l(k) = volume_l(k) + dx * dy * dzw(k+1)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
!
!--    Sum-up building volume from all subdomains
#if defined( __parallel )
       CALL MPI_ALLREDUCE( volume_l, volume, SIZE(volume), MPI_REAL, MPI_SUM, comm2d, ierr )
#else
       volume = volume_l
#endif
!
!--    Save total building volume as well as local fraction on volume on building data structure.
       IF ( ALLOCATED( buildings(nb)%volume ) )  buildings(nb)%volume = volume
!
!--    Determine fraction of local on total building volume
       IF ( buildings(nb)%on_pe )  buildings(nb)%vol_frac = volume_l / volume
!
!--    Calculate total building volume
       IF ( ALLOCATED( buildings(nb)%volume ) )  buildings(nb)%vol_tot = SUM( buildings(nb)%volume )

       DEALLOCATE( volume   )
       DEALLOCATE( volume_l )

    ENDDO
!
!-- Allocate arrays for indoor temperature.
    DO  nb = 1, num_build
       IF ( buildings(nb)%on_pe )  THEN
          ALLOCATE( buildings(nb)%t_in(buildings(nb)%kb_min:buildings(nb)%kb_max)   )
          ALLOCATE( buildings(nb)%t_in_l(buildings(nb)%kb_min:buildings(nb)%kb_max) )
          buildings(nb)%t_in   = 0.0_wp
          buildings(nb)%t_in_l = 0.0_wp
       ENDIF
    ENDDO
!
!-- Allocate arrays for number of facades per height level. Distinguish between horizontal and
!-- vertical facades.
    DO  nb = 1, num_build
       IF ( buildings(nb)%on_pe )  THEN
          ALLOCATE( buildings(nb)%num_facade_h(buildings(nb)%kb_min:buildings(nb)%kb_max) )
          ALLOCATE( buildings(nb)%num_facade_v(buildings(nb)%kb_min:buildings(nb)%kb_max) )

          buildings(nb)%num_facade_h = 0
          buildings(nb)%num_facade_v = 0
       ENDIF
    ENDDO
!
!-- Determine number of facade elements per building on local subdomain.
    buildings(:)%num_facades_per_building_l = 0
    DO  m = 1, surf_usm%ns
!
!--    For the current facade element determine corresponding building index.
!--    First, obtain j,j,k indices of the building. Please note the offset between facade/surface
!--    element and building location (for horizontal surface elements the horizontal offsets are
!--    zero).
       i = surf_usm%i(m) + surf_usm%ioff(m)
       j = surf_usm%j(m) + surf_usm%joff(m)
       k = surf_usm%k(m) + surf_usm%koff(m)
!
!--    Determine building index and check whether building is on PE.
       nb = MINLOC( ABS( buildings(:)%id - building_id_f%var(j,i) ), DIM=1 )

       IF ( buildings(nb)%on_pe )  THEN
!
!--       Horizontal facade elements (roofs, etc.)
          IF ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) )  THEN
!
!--          Count number of facade elements at each height level.
             buildings(nb)%num_facade_h(k) = buildings(nb)%num_facade_h(k) + 1
!
!--          Moreover, sum up number of local facade elements per building.
             buildings(nb)%num_facades_per_building_l = buildings(nb)%num_facades_per_building_l + 1
!
!--       Vertical facades
          ELSE
             buildings(nb)%num_facade_v(k) = buildings(nb)%num_facade_v(k) + 1
             buildings(nb)%num_facades_per_building_l = buildings(nb)%num_facades_per_building_l + 1
          ENDIF
       ENDIF
    ENDDO
!
!-- Determine total number of facade elements per building and assign number to building data type.
    DO  nb = 1, num_build
!
!--    Allocate dummy array used for summing-up facade elements.
!--    Please note, dummy arguments are necessary as building-date type arrays are not necessarily
!--    allocated on all PEs.
       ALLOCATE( num_facades_h(buildings(nb)%kb_min:buildings(nb)%kb_max) )
       ALLOCATE( num_facades_v(buildings(nb)%kb_min:buildings(nb)%kb_max) )
       ALLOCATE( receive_dum_h(buildings(nb)%kb_min:buildings(nb)%kb_max) )
       ALLOCATE( receive_dum_v(buildings(nb)%kb_min:buildings(nb)%kb_max) )
       num_facades_h = 0
       num_facades_v = 0
       receive_dum_h = 0
       receive_dum_v = 0

       IF ( buildings(nb)%on_pe )  THEN
          num_facades_h = buildings(nb)%num_facade_h
          num_facades_v = buildings(nb)%num_facade_v
       ENDIF

#if defined( __parallel )
       CALL MPI_ALLREDUCE( num_facades_h,                                                          &
                           receive_dum_h,                                                          &
                           buildings(nb)%kb_max - buildings(nb)%kb_min + 1,                        &
                           MPI_INTEGER,                                                            &
                           MPI_SUM,                                                                &
                           comm2d,                                                                 &
                           ierr )

       CALL MPI_ALLREDUCE( num_facades_v,                                                          &
                           receive_dum_v,                                                          &
                           buildings(nb)%kb_max - buildings(nb)%kb_min + 1,                        &
                           MPI_INTEGER,                                                            &
                           MPI_SUM,                                                                &
                           comm2d,                                                                 &
                           ierr )
       IF ( ALLOCATED( buildings(nb)%num_facade_h ) )  buildings(nb)%num_facade_h = receive_dum_h
       IF ( ALLOCATED( buildings(nb)%num_facade_v ) )  buildings(nb)%num_facade_v = receive_dum_v
#else
       buildings(nb)%num_facade_h = num_facades_h
       buildings(nb)%num_facade_v = num_facades_v
#endif

!
!--    Deallocate dummy arrays
       DEALLOCATE( num_facades_h )
       DEALLOCATE( num_facades_v )
       DEALLOCATE( receive_dum_h )
       DEALLOCATE( receive_dum_v )
!
!--    Allocate index arrays which link facade elements with surface-data type.
!--    Please note, no height levels are considered here (information is stored in surface-data type
!--    itself).
       IF ( buildings(nb)%on_pe )  THEN
!
!--       Determine number of facade elements per building.
          buildings(nb)%num_facades_per_building_h = SUM( buildings(nb)%num_facade_h )
          buildings(nb)%num_facades_per_building_v = SUM( buildings(nb)%num_facade_v )
          buildings(nb)%num_facades_per_building   = buildings(nb)%num_facades_per_building_h +    &
                                                     buildings(nb)%num_facades_per_building_v
!
!--       Allocate array that links the building with the urban-type surfaces.
!--       Please note, the linking array is allocated over all facade elements, which is
!--       required in case a building is located at the subdomain boundaries, where the building and
!--       the corresponding surface elements are located on different subdomains.
          ALLOCATE( buildings(nb)%m(1:buildings(nb)%num_facades_per_building_l) )

          ALLOCATE( buildings(nb)%theta_m_t_prev(1:buildings(nb)%num_facades_per_building_l) )
       ENDIF

       IF ( buildings(nb)%on_pe )  THEN
          ALLOCATE( buildings(nb)%vpf(buildings(nb)%kb_min:buildings(nb)%kb_max) )
          buildings(nb)%vpf = 0.0_wp

          facade_area_v = 0.0_wp
          DO  k = buildings(nb)%kb_min, buildings(nb)%kb_max
             facade_area_v = facade_area_v + buildings(nb)%num_facade_v(k) * dzw(k+1) * dx
          ENDDO
!
!--       Determine volume per total facade area (vpf). For the horizontal facade area
!--       num_facades_per_building_h can be taken, multiplied with dx*dy.
!--       However, due to grid stretching, vertical facade elements must be summed-up vertically.
!--       Please note, if dx /= dy, an error is made!
          buildings(nb)%vpf = buildings(nb)%vol_tot /                                              &
                              ( buildings(nb)%num_facades_per_building_h * dx * dy + facade_area_v )
!
!--       Determine floor-area-per-facade.
          buildings(nb)%fapf = buildings(nb)%num_facades_per_building_h     * dx * dy              &
                               / ( buildings(nb)%num_facades_per_building_h * dx * dy              &
                                   + facade_area_v )
       ENDIF
    ENDDO
!
!-- Link facade elements with surface data type.
!-- Allocate array for counting.
    ALLOCATE( n_fa(1:num_build) )
    n_fa = 1

    DO  m = 1, surf_usm%ns
       i = surf_usm%i(m) + surf_usm%ioff(m)
       j = surf_usm%j(m) + surf_usm%joff(m)

       nb = MINLOC( ABS( buildings(:)%id - building_id_f%var(j,i) ), DIM=1 )

       IF ( buildings(nb)%on_pe )  THEN
          buildings(nb)%m(n_fa(nb)) = m
          n_fa(nb) = n_fa(nb) + 1
       ENDIF
    ENDDO

    DEALLOCATE( n_fa )
!
!-- Initialize building parameters, first by mean building type. Note, in this case all buildings
!-- have the same type.
!-- In a second step initialize with building tpyes from static input file, where building types can
!-- be individual for each building.
    buildings(:)%lambda_layer3       = building_pars(31,building_type)
    buildings(:)%s_layer3            = building_pars(44,building_type)
    buildings(:)%f_c_win             = building_pars(119,building_type)
    buildings(:)%g_value_win         = building_pars(120,building_type)
    buildings(:)%u_value_win         = building_pars(121,building_type)
    buildings(:)%eta_ve              = building_pars(124,building_type)
    buildings(:)%factor_a            = building_pars(125,building_type)
    buildings(:)%factor_c            = building_pars(126,building_type)
    buildings(:)%lambda_at           = building_pars(127,building_type)
    buildings(:)%theta_int_h_set     = building_pars(13,building_type)
    buildings(:)%theta_int_c_set     = building_pars(12,building_type)
    buildings(:)%q_h_max             = building_pars(128,building_type)
    buildings(:)%q_c_max             = building_pars(129,building_type)
    buildings(:)%qint_high           = building_pars(130,building_type)
    buildings(:)%qint_low            = building_pars(131,building_type)
    buildings(:)%height_storey       = building_pars(132,building_type)
    buildings(:)%height_cei_con      = building_pars(133,building_type)
    buildings(:)%params_waste_heat_h = building_pars(134,building_type)
    buildings(:)%params_waste_heat_c = building_pars(135,building_type)
!
!-- Initialize seasonal dependent parameters, depending on day of the year.
!-- First, calculated day of the year.
    CALL get_date_time( time_since_reference_point, day_of_year = day_of_year )
!
!-- Summer is defined in between northward- and southward equinox.
    IF ( day_of_year >= northward_equinox  .AND.  day_of_year <= southward_equinox )  THEN
       buildings(:)%air_change_low      = summer_pars(0,building_type)
       buildings(:)%air_change_high     = summer_pars(1,building_type)
    ELSE
       buildings(:)%air_change_low      = winter_pars(0,building_type)
       buildings(:)%air_change_high     = winter_pars(1,building_type)
    ENDIF
!
!-- Initialize ventilation load. Please note, building types > 7 are actually not allowed (check
!-- already in urban_surface_mod and netcdf_data_input_mod.
!-- However, the building data base may be later extended.
    IF ( building_type ==  1  .OR.  building_type ==  2  .OR.                                      &
         building_type ==  3  .OR.  building_type == 10  .OR.                                      &
         building_type == 11  .OR.  building_type == 12 )  THEN
       buildings(:)%ventilation_int_loads = 1
!
!-- Office, building with large windows
    ELSEIF ( building_type ==  4  .OR.  building_type ==  5  .OR.                                  &
             building_type ==  6  .OR.  building_type ==  7  .OR.                                  &
             building_type ==  8  .OR.  building_type ==  9)  THEN
       buildings(:)%ventilation_int_loads = 2
!
!-- Industry, hospitals
    ELSEIF ( building_type == 13  .OR.  building_type == 14  .OR.                                  &
             building_type == 15  .OR.  building_type == 16  .OR.                                  &
             building_type == 17  .OR.  building_type == 18 )  THEN
       buildings(:)%ventilation_int_loads = 3
    ENDIF
!
!-- Initialization of building parameters - level 2
    IF ( building_type_f%from_file )  THEN
       DO  i = nxl, nxr
          DO  j = nys, nyn
              IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
                 nb = MINLOC( ABS( buildings(:)%id - building_id_f%var(j,i) ), DIM=1 )
                 bt = building_type_f%var(j,i)

                 buildings(nb)%lambda_layer3       = building_pars(31,bt)
                 buildings(nb)%s_layer3            = building_pars(44,bt)
                 buildings(nb)%f_c_win             = building_pars(119,bt)
                 buildings(nb)%g_value_win         = building_pars(120,bt)
                 buildings(nb)%u_value_win         = building_pars(121,bt)
                 buildings(nb)%eta_ve              = building_pars(124,bt)
                 buildings(nb)%factor_a            = building_pars(125,bt)
                 buildings(nb)%factor_c            = building_pars(126,bt)
                 buildings(nb)%lambda_at           = building_pars(127,bt)
                 buildings(nb)%theta_int_h_set     = building_pars(13,bt)
                 buildings(nb)%theta_int_c_set     = building_pars(12,bt)
                 buildings(nb)%q_h_max             = building_pars(128,bt)
                 buildings(nb)%q_c_max             = building_pars(129,bt)
                 buildings(nb)%qint_high           = building_pars(130,bt)
                 buildings(nb)%qint_low            = building_pars(131,bt)
                 buildings(nb)%height_storey       = building_pars(132,bt)
                 buildings(nb)%height_cei_con      = building_pars(133,bt)
                 buildings(nb)%params_waste_heat_h = building_pars(134,bt)
                 buildings(nb)%params_waste_heat_c = building_pars(135,bt)

              IF ( day_of_year >= northward_equinox  .AND.  day_of_year <= southward_equinox )  THEN
                 buildings(nb)%air_change_low      = summer_pars(0,bt)
                 buildings(nb)%air_change_high     = summer_pars(1,bt)
              ELSE
                 buildings(nb)%air_change_low      = winter_pars(0,bt)
                 buildings(nb)%air_change_high     = winter_pars(1,bt)
              ENDIF

!
!--              Initialize ventilaation load. Please note, building types > 7
!--              are actually not allowed (check already in urban_surface_mod
!--              and netcdf_data_input_mod. However, the building data base may
!--              be later extended.
                 IF ( bt ==  1  .OR.  bt ==  2  .OR.                                               &
                      bt ==  3  .OR.  bt == 10  .OR.                                               &
                      bt == 11  .OR.  bt == 12 )  THEN
                    buildings(nb)%ventilation_int_loads = 1
!
!--              Office, building with large windows
                 ELSEIF ( bt ==  4  .OR.  bt ==  5  .OR.                                           &
                          bt ==  6  .OR.  bt ==  7  .OR.                                           &
                          bt ==  8  .OR.  bt ==  9)  THEN
                    buildings(nb)%ventilation_int_loads = 2
!
!--              Industry, hospitals
                 ELSEIF ( bt == 13  .OR.  bt == 14  .OR.                                           &
                          bt == 15  .OR.  bt == 16  .OR.                                           &
                          bt == 17  .OR.  bt == 18 )  THEN
                    buildings(nb)%ventilation_int_loads = 3
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
    ENDIF
!
!-- Calculation of surface-related heat transfer coeffiecient out of standard u-values from building
!-- database.
!-- Only amount of extern and surface is used.
!-- Otherwise amount between air and surface taken account twice.
    DO nb = 1, num_build
       IF ( buildings(nb)%on_pe ) THEN
          du_win_tmp = 1.0_wp / buildings(nb)%u_value_win
          u_tmp = buildings(nb)%u_value_win * ( du_win_tmp / ( du_win_tmp -                        &
                  0.125_wp + ( 1.0_wp / h_is ) ) )

          du_tmp = 1.0_wp / u_tmp

          buildings(nb)%h_es = 1.0_wp / ( du_tmp - ( 1.0_wp / h_is ) )

       ENDIF
    ENDDO
!
!-- Initialize indoor temperature. Actually only for output at initial state.
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
       DO  nb = 1, num_build
          IF ( buildings(nb)%on_pe )  THEN
             buildings(nb)%t_in(:) = initial_indoor_temperature

!
!--          (after first loop, use theta_m_t as theta_m_t_prev)
             buildings(nb)%theta_m_t_prev(:) = initial_indoor_temperature

          ENDIF
       ENDDO
!
!-- Initialize indoor temperature at previous timestep.
    ELSE
       DO  nb = 1, num_build
          IF ( buildings(nb)%on_pe )  THEN
!
!--          Mean indoor temperature can be initialized with initial value. This is just
!--          used for output.
             buildings(nb)%t_in(:) = initial_indoor_temperature
!
!--          Initialize theta_m_t_prev arrays. The respective data during the restart mechanism
!--          is stored on the surface-data array.
             DO  fa = 1, buildings(nb)%num_facades_per_building_l
!
!--             Determine indices where corresponding surface-type information is stored.
                m = buildings(nb)%m(fa)
                buildings(nb)%theta_m_t_prev(fa) = surf_usm%t_prev(m)
             ENDDO
          ENDIF
       ENDDO
    ENDIF
!
!-- Initially, compute 10-cm temperature.
    CALL calc_pt_10cm

    CALL location_message( 'initializing indoor model', 'finished' )

 END SUBROUTINE im_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of arrays related to the indoor model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_init_arrays

!
!-- Allocate memory for 10-cm temperature, waste heat and indoor temperature.
    ALLOCATE ( surf_usm%pt_10cm(1:surf_usm%ns) )
    ALLOCATE ( surf_usm%t_prev(1:surf_usm%ns) )
    ALLOCATE ( surf_usm%waste_heat(1:surf_usm%ns) )

    surf_usm%pt_10cm    = 0.0_wp
    surf_usm%t_prev     = 0.0_wp
    surf_usm%waste_heat = 0.0_wp

 END SUBROUTINE im_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Main part of the indoor model.
!> Computation of energy demand of the building volumes, waste heat, mean indoor temperature and
!> inner wall temperature.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_main_heatcool

    USE control_parameters,                                                                        &
        ONLY:  time_since_reference_point

    USE urban_surface_mod,                                                                         &
        ONLY:  nzt_wall,                                                                           &
               t_green,                                                                            &
               t_wall,                                                                             &
               t_window


    INTEGER(iwp) ::  fa   !< running index for facade elements of each building
    INTEGER(iwp) ::  i    !< index of facade-adjacent atmosphere grid point in x-direction
    INTEGER(iwp) ::  j    !< index of facade-adjacent atmosphere grid point in y-direction
    INTEGER(iwp) ::  k    !< index of facade-adjacent atmosphere grid point in z-direction
    INTEGER(iwp) ::  kk   !< vertical index of indoor grid point adjacent to facade
    INTEGER(iwp) ::  m    !< running index surface elements
    INTEGER(iwp) ::  nb   !< running index for buildings

    LOGICAL  ::  during_spinup                    !< flag indicating that the simulation is still in wall/soil spinup

    REAL(wp) ::  frac_green                       !< dummy for green fraction
    REAL(wp) ::  frac_wall                        !< dummy for wall fraction
    REAL(wp) ::  frac_win                         !< dummy for window fraction
!     REAL(wp) ::  indoor_wall_window_temperature   !< weighted temperature of innermost wall/window layer
    REAL(wp) ::  indoor_wall_temperature          !< temperature of innermost wall layer evtl in im_calc_temperatures einfügen
    REAL(wp) ::  near_facade_temperature          !< outside air temperature 10cm away from facade
    REAL(wp) ::  second_of_day                    !< second of the current day
    REAL(wp) ::  time_utc_hour                    !< time of day (hour UTC)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_in_l_send   !< dummy send buffer used for summing-up indoor temperature per kk-level
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  t_in_recv     !< dummy recv buffer used for summing-up indoor temperature per kk-level
!
!-- Determine time of day in hours.
    CALL get_date_time( time_since_reference_point, second_of_day=second_of_day )
    time_utc_hour = second_of_day / seconds_per_hour
!
!-- Check if the simulation is still in wall/soil spinup mode
    during_spinup = MERGE( .TRUE., .FALSE., time_since_reference_point < 0.0_wp )
!
!-- Update 10-cm temperature.
    CALL calc_pt_10cm
!
!-- Following calculations must be done for each facade element.
    DO  nb = 1, num_build
!
!--    First, check whether building is present on local subdomain.
       IF ( buildings(nb)%on_pe )  THEN
!
!--       Determine daily schedule. 08:00-18:00 = 1, other hours = 0.
!--       Residental Building, panel WBS 70
          IF ( buildings(nb)%ventilation_int_loads == 1 )  THEN
             IF ( time_utc_hour >= 8.0_wp  .AND.  time_utc_hour <= 18.0_wp )  THEN
                schedule_d = 0
             ELSE
                schedule_d = 1
             ENDIF
          ENDIF
!
!--       Office, building with large windows
          IF ( buildings(nb)%ventilation_int_loads == 2 )  THEN
             IF ( time_utc_hour >= 8.0_wp  .AND.  time_utc_hour <= 18.0_wp )  THEN
                schedule_d = 1
             ELSE
                schedule_d = 0
             ENDIF
          ENDIF
!
!--       Industry, hospitals
          IF ( buildings(nb)%ventilation_int_loads == 3 )  THEN
             IF ( time_utc_hour >= 6.0_wp  .AND.  time_utc_hour <= 22.0_wp )  THEN
                schedule_d = 1
             ELSE
                schedule_d = 0
             ENDIF
          ENDIF
!
!--       Initialize/reset indoor temperature - note, this is only for output
          buildings(nb)%t_in_l = 0.0_wp
!
!--       Horizontal surfaces
          DO  fa = 1, buildings(nb)%num_facades_per_building_l
!
!--          Determine indices where corresponding surface-type information is stored.
             m = buildings(nb)%m(fa)
!
!--          During spinup set window fraction to zero and add these to wall fraction.
             frac_win   = MERGE( surf_usm%frac(m,ind_wat_win), 0.0_wp, .NOT. during_spinup )
             frac_wall  = MERGE( surf_usm%frac(m,ind_veg_wall),                                    &
                                 surf_usm%frac(m,ind_veg_wall) +                                   &
                                 surf_usm%frac(m,ind_wat_win),                                     &
                                 .NOT. during_spinup )
             frac_green = surf_usm%frac(m,ind_pav_green)
!
!--          Determine building height level index.
             kk = surf_usm%k(m) + surf_usm%koff(m)
!
!--          Building geometries --> not time-dependent
             IF ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) )  THEN
                facade_element_area = dx * dy                                     !< [m2] surface area per facade element
             ELSEIF( surf_usm%northward(m)  .OR.  surf_usm%southward(m) )  THEN
                facade_element_area = dx * dzw(kk+1)                              !< [m2] surface area per facade element
             ELSEIF( surf_usm%eastward(m)  .OR.  surf_usm%westward(m) )  THEN
                facade_element_area = dy * dzw(kk+1)                              !< [m2] surface area per facade element
             ENDIF

             floor_area_per_facade        = buildings(nb)%fapf                    !< [m2/m2] floor area per facade area
             indoor_volume_per_facade     = buildings(nb)%vpf(kk)                 !< [m3/m2] indoor air volume per facade area
             buildings(nb)%area_facade    = facade_element_area *                                  &
                                            ( buildings(nb)%num_facades_per_building_h +           &
                                              buildings(nb)%num_facades_per_building_v ) !< [m2] area of total facade
             window_area_per_facade       = frac_win  * facade_element_area       !< [m2] window area per facade element

             buildings(nb)%net_floor_area = buildings(nb)%vol_tot / buildings(nb)%height_storey
             total_area                   = buildings(nb)%net_floor_area          !< [m2] area of all surfaces
                                                                                  !< pointing to zone  Eq. (9) according to section 7.2.2.2
             a_m                          = buildings(nb)%factor_a * total_area *                  &
                                            ( facade_element_area / buildings(nb)%area_facade ) *  &
                                            buildings(nb)%lambda_at                                 !< [m2] standard values
                                                                                                    !< according to Table 12 section 12.3.1.2  (calculate over Eq. (65) according to section 12.3.1.2)
             c_m                          = buildings(nb)%factor_c * total_area *                  &
                                            ( facade_element_area / buildings(nb)%area_facade )     !< [J/K] standard values
                                                                                                    !< according to table 12 section 12.3.1.2 (calculate over Eq. (66) according to section 12.3.1.2)
!
!--          Calculation of heat transfer coefficient for transmission --> not time-dependent
             h_t_es   = window_area_per_facade * buildings(nb)%h_es                                   !< [W/K] only for windows

             h_t_is  = buildings(nb)%area_facade * h_is                                               !< [W/K] with h_is = 3.45 W /
                                                                                                      !< (m2 K) between surface and air, Eq. (9)
             h_t_ms  = a_m * h_ms                                                                     !< [W/K] with h_ms = 9.10 W /
                                                                                                      !< (m2 K) between component and surface, Eq. (64)
             h_t_wall  = 1.0_wp / ( 1.0_wp / ( ( facade_element_area - window_area_per_facade )    &  !< [W/K]
                                    * buildings(nb)%lambda_layer3 / buildings(nb)%s_layer3 * 0.5_wp&
                                             ) + 1.0_wp / h_t_ms )                                    !< [W/K] opaque components
             h_t_wm  = 1.0_wp / ( 1.0_wp / h_t_wall - 1.0_wp / h_t_ms )                               !< [W/K] emmision Eq. (63),
                                                                                                      !< Section 12.2.2
!
!--          Internal air loads dependent on the occupacy of the room.
!--          Basical internal heat gains (qint_low) with additional internal heat gains by occupancy
!--          (qint_high) (0.5*phi_int).
             phi_ia = 0.5_wp * ( ( buildings(nb)%qint_high * schedule_d + buildings(nb)%qint_low ) &
                              * floor_area_per_facade )

             IF ( surf_usm%upward(m)  .OR.  surf_usm%downward(m) )  THEN
                q_int = phi_ia / total_area
             ELSE
                q_int = phi_ia
             ENDIF
!
!--          Airflow dependent on the occupacy of the room.
!--          Basical airflow (air_change_low) with additional airflow gains by occupancy
!--          (air_change_high)
             air_change = ( buildings(nb)%air_change_high * schedule_d +                           &
                            buildings(nb)%air_change_low )  !< [1/h]?
!
!--          Heat transfer of ventilation.
!--          Not less than 0.01 W/K to avoid division by 0 in further calculations with heat
!--          capacity of air 0.33 Wh/m2K.
             h_v   = MAX( 0.01_wp , ( air_change * indoor_volume_per_facade *                      &
                                      0.33_wp * (1.0_wp - buildings(nb)%eta_ve ) ) )    !< [W/K] from ISO 13789 Eq.(10)

!--          Heat transfer coefficient auxiliary variables
             h_t_1 = 1.0_wp / ( ( 1.0_wp / h_v )   + ( 1.0_wp / h_t_is ) )  !< [W/K] Eq. (C.6)
             h_t_2 = h_t_1 + h_t_es                                         !< [W/K] Eq. (C.7)
             h_t_3 = 1.0_wp / ( ( 1.0_wp / h_t_2 ) + ( 1.0_wp / h_t_ms ) )  !< [W/K] Eq. (C.8)
!
!--          Net short-wave radiation through window area (was i_global)
             net_sw_in = surf_usm%rad_sw_in(m) - surf_usm%rad_sw_out(m)
!
!--          Quantities needed for im_calc_temperatures
             i = surf_usm%i(m)
             j = surf_usm%j(m)
             k = surf_usm%k(m)
             near_facade_temperature = surf_usm%pt_10cm(m)

             indoor_wall_temperature = frac_wall  * t_wall%val(nzt_wall,m)                         &
                                     + frac_win   * t_window%val(nzt_wall,m)                       &
                                     + frac_green * t_green%val(nzt_wall,m)
!
!--          Solar thermal gains. If net_sw_in larger than sun-protection threshold parameter
!--          (params_solar_protection), sun protection will be activated.
             IF ( net_sw_in <= params_solar_protection )  THEN
                solar_protection_off = 1
                solar_protection_on  = 0
             ELSE
                solar_protection_off = 0
                solar_protection_on  = 1
             ENDIF
!
!--          Calculation of total heat gains from net_sw_in through windows [W] in respect on
!--          automatic sun protection.
!--          DIN 4108 - 2 chap.8
             phi_sol = (   window_area_per_facade * net_sw_in * solar_protection_off               &
                         + window_area_per_facade * net_sw_in * buildings(nb)%f_c_win *            &
                           solar_protection_on )                                                   &
                       * buildings(nb)%g_value_win * ( 1.0_wp - params_f_f ) * params_f_w
             q_sol = phi_sol
!
!--          Calculation of the mass specific thermal load for internal and external heatsources of
!--          the inner node.
             phi_m   = (a_m / total_area) * ( phi_ia + phi_sol )                                    !< [W] Eq. (C.2) with
                                                                                                    !< phi_ia=0,5*phi_int
             q_c_m = phi_m
!
!--          Calculation mass specific thermal load implied non thermal mass
             phi_st  =   ( 1.0_wp - ( a_m / total_area ) - ( h_t_es / ( 9.1_wp * total_area ) ) )  &
                       * ( phi_ia + phi_sol )                                                       !< [W] Eq. (C.3) with
                                                                                                    !< phi_ia=0,5*phi_int
             q_c_st = phi_st
!
!--          Calculations for deriving indoor temperature and heat flux into the wall
!--          Step 1: indoor temperature without heating and cooling
!--          section C.4.1 Picture C.2 zone 3)
             phi_hc_nd = 0.0_wp

             CALL  im_calc_temperatures ( i, j, k, indoor_wall_temperature,                        &
                                          near_facade_temperature, phi_hc_nd,                      &
                                          buildings(nb)%theta_m_t_prev(fa) )
!
!--          If air temperature between border temperatures of heating and cooling, assign output
!--          variable, then ready.
             IF ( buildings(nb)%theta_int_h_set <= theta_air  .AND.                                &
                  theta_air <= buildings(nb)%theta_int_c_set )  THEN
                phi_hc_nd_ac = 0.0_wp
                phi_hc_nd    = phi_hc_nd_ac
                theta_air_ac = theta_air
!
!--          Step 2: Else, apply 10 W/m2 heating/cooling power and calculate indoor temperature
!--          again.
             ELSE
!
!--             Temperature not correct, calculation method according to section C4.2
                theta_air_0 = theta_air       !< temperature without heating/cooling
!
!--             Heating or cooling?
                IF ( theta_air_0 > buildings(nb)%theta_int_c_set )  THEN
                   theta_air_set = buildings(nb)%theta_int_c_set
                ELSE
                   theta_air_set = buildings(nb)%theta_int_h_set
                ENDIF
!
!--             Calculate the temperature with phi_hc_nd_10
                phi_hc_nd_10 = 10.0_wp * floor_area_per_facade
                phi_hc_nd    = phi_hc_nd_10

                CALL  im_calc_temperatures ( i, j, k, indoor_wall_temperature,                     &
                                             near_facade_temperature, phi_hc_nd,                   &
                                             buildings(nb)%theta_m_t_prev(fa) )
                theta_air_10 = theta_air                                                !< temperature with 10 W/m2 of heating
!
!--             Avoid division by zero at first timestep where the denominator can become zero.
                IF ( ABS( theta_air_10  - theta_air_0 ) < 1E-10_wp )  THEN
                   phi_hc_nd_un = phi_hc_nd_10
                ELSE
                   phi_hc_nd_un = phi_hc_nd_10 * ( theta_air_set - theta_air_0 )                   &
                                               / ( theta_air_10  - theta_air_0 )        !< Eq. (C.13)
                ENDIF
!
!--             Step 3: with temperature ratio to determine the heating or cooling capacity.
!--             If necessary, limit the power to maximum power.
!--             section C.4.1 Picture C.2 zone 2) and 4)
                buildings(nb)%phi_c_max = buildings(nb)%q_c_max * floor_area_per_facade
                buildings(nb)%phi_h_max = buildings(nb)%q_h_max * floor_area_per_facade
                IF ( buildings(nb)%phi_c_max < phi_hc_nd_un  .AND.                                 &
                     phi_hc_nd_un < buildings(nb)%phi_h_max )  THEN
                   phi_hc_nd_ac = phi_hc_nd_un
                   phi_hc_nd = phi_hc_nd_un
                ELSE
!
!--             Step 4: inner temperature with maximum heating (phi_hc_nd_un positive) or cooling
!--                     (phi_hc_nd_un negative)
!--             section C.4.1 Picture C.2 zone 1) and 5)
                   IF ( phi_hc_nd_un > 0.0_wp )  THEN
                      phi_hc_nd_ac = buildings(nb)%phi_h_max                            !< Limit heating
                   ELSE
                      phi_hc_nd_ac = buildings(nb)%phi_c_max                            !< Limit cooling
                   ENDIF
                ENDIF
                phi_hc_nd = phi_hc_nd_ac
!
!--             Calculate the temperature with phi_hc_nd_ac (new)
                CALL  im_calc_temperatures ( i, j, k, indoor_wall_temperature,                     &
                                             near_facade_temperature, phi_hc_nd,                   &
                                             buildings(nb)%theta_m_t_prev(fa) )
                theta_air_ac = theta_air
             ENDIF
!
!--          Update theta_m_t_prev
             buildings(nb)%theta_m_t_prev(fa) = theta_m_t


             q_vent = h_v * ( theta_air - near_facade_temperature )
!
!--          Calculate the operating temperature with weighted mean temperature of air and mean
!--          solar temperature.
!--          Will be used for thermal comfort calculations.
             theta_op     = 0.3_wp * theta_air_ac + 0.7_wp * theta_s          !< [degree_C] operative Temperature Eq. (C.12)

!              surf_usm%t_indoor(m) = theta_op                               !< not integrated now
!
!--          Heat flux into the wall. Value needed in urban_surface_mod to
!--          calculate heat transfer through wall layers towards the facade
!--          (use c_p * rho_surface to convert [W/m2] into [K m/s])
             IF ( (facade_element_area - window_area_per_facade) > 0.0_wp )  THEN
                q_wall = h_t_wm * ( indoor_wall_temperature - theta_m )                            &
                                / ( facade_element_area - window_area_per_facade )
             ELSE
                q_wall = 0.0_wp
             ENDIF

             IF ( window_area_per_facade > 0.0_wp )  THEN
                q_win = h_t_es * ( pt(k,j,i) - theta_s ) / ( window_area_per_facade )
             ELSE
                q_win = 0.0_wp
             ENDIF
!
!--          Transfer q_wall & q_win back to USM (innermost wall/window layer)
             surf_usm%iwghf_eb(m)        = -q_wall
             surf_usm%iwghf_eb_window(m) = -q_win
!
!--          Sum up operational indoor temperature per kk-level. Further below, this temperature is
!--          reduced by MPI to one temperature per kk-level and building (processor overlapping).
             buildings(nb)%t_in_l(kk) = buildings(nb)%t_in_l(kk) + theta_op
!
!--          Calculation of waste heat.
!--          Anthropogenic heat output.
             IF ( phi_hc_nd_ac > 0.0_wp )  THEN
                heating_on = 1
                cooling_on = 0
             ELSE
                heating_on = 0
                cooling_on = -1
             ENDIF

             q_waste_heat = ( phi_hc_nd * (                                                        &
                              buildings(nb)%params_waste_heat_h * heating_on +                     &
                              buildings(nb)%params_waste_heat_c * cooling_on )                     &
                            ) / facade_element_area  !< [W/m2] , observe the directional convention in PALM!
!
!--          Store waste heat and previous previous indoor temperature on surface-data type.
!--          These will be used in the urban-surface model.
             surf_usm%t_prev(m) = buildings(nb)%theta_m_t_prev(fa)
             surf_usm%waste_heat(m) = q_waste_heat
          ENDDO
       ENDIF !< buildings(nb)%on_pe
    ENDDO !< buildings loop

!
!-- Determine the mean building temperature.
    DO  nb = 1, num_build
!
!--    Allocate dummy array used for summing-up facade elements.
!--    Please note, dummy arguments are necessary as building-date type arrays are not necessarily
!--    allocated on all PEs.
       ALLOCATE( t_in_l_send(buildings(nb)%kb_min:buildings(nb)%kb_max) )
       ALLOCATE( t_in_recv(buildings(nb)%kb_min:buildings(nb)%kb_max) )
       t_in_l_send = 0.0_wp
       t_in_recv   = 0.0_wp

       IF ( buildings(nb)%on_pe )  THEN
          t_in_l_send = buildings(nb)%t_in_l
       ENDIF


#if defined( __parallel )
       CALL MPI_ALLREDUCE( t_in_l_send,                                                            &
                           t_in_recv,                                                              &
                           buildings(nb)%kb_max - buildings(nb)%kb_min + 1,                        &
                           MPI_REAL,                                                               &
                           MPI_SUM,                                                                &
                           comm2d,                                                                 &
                           ierr )

       IF ( ALLOCATED( buildings(nb)%t_in ) )  buildings(nb)%t_in = t_in_recv
#else
       IF ( ALLOCATED( buildings(nb)%t_in ) )  buildings(nb)%t_in = buildings(nb)%t_in_l
#endif

       IF ( ALLOCATED( buildings(nb)%t_in ) )  THEN
!
!--       Average indoor temperature. Note, in case a building is completely surrounded by higher
!--       buildings, it may have no facade elements at some height levels, which will lead to a
!--       division by zero.
          DO  k = buildings(nb)%kb_min, buildings(nb)%kb_max
             IF ( buildings(nb)%num_facade_h(k) + buildings(nb)%num_facade_v(k) > 0 )  THEN
                buildings(nb)%t_in(k) = buildings(nb)%t_in(k) /                                    &
                                        REAL( buildings(nb)%num_facade_h(k) +                      &
                                              buildings(nb)%num_facade_v(k), KIND = wp )
             ENDIF
          ENDDO
!
!--       If indoor temperature is not defined because of missing facade elements, the values from
!--       the above-lying level will be taken.
!--       At least at the top of the buildings facades are defined, so that at least there an indoor
!--       temperature is defined. This information will propagate downwards the building.
          DO  k = buildings(nb)%kb_max-1, buildings(nb)%kb_min, -1
             IF ( buildings(nb)%num_facade_h(k) + buildings(nb)%num_facade_v(k) <= 0 )  THEN
                buildings(nb)%t_in(k) = buildings(nb)%t_in(k+1)
             ENDIF
          ENDDO
       ENDIF


!
!--    Deallocate dummy arrays
       DEALLOCATE( t_in_l_send )
       DEALLOCATE( t_in_recv )

    ENDDO

 END SUBROUTINE im_main_heatcool


!--------------------------------------------------------------------------------------------------!
! Description:
!-------------
!> Check data output for plant canopy model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_check_data_output( var, unit )

    CHARACTER (LEN=*) ::  unit   !<
    CHARACTER (LEN=*) ::  var    !<

    SELECT CASE ( TRIM( var ) )


        CASE ( 'im_hf_roof')
           unit = 'W m-2'

        CASE ( 'im_hf_wall_win' )
           unit = 'W m-2'

        CASE ( 'im_hf_wall_win_waste' )
           unit = 'W m-2'

        CASE ( 'im_hf_roof_waste' )
           unit = 'W m-2'

        CASE ( 'im_theta_10cm_roof', 'im_theta_10cm_wall' )
           unit = 'K'

        CASE ( 'im_t_indoor_mean' )
           unit = 'K'

        CASE ( 'im_t_indoor_roof' )
           unit = 'K'

        CASE ( 'im_t_indoor_wall_win' )
           unit = 'K'

        CASE ( 'im_t_indoor_wall' )
           unit = 'K'

        CASE DEFAULT
           unit = 'illegal'

    END SELECT

 END SUBROUTINE


!--------------------------------------------------------------------------------------------------!
! Description:
!-------------
!> Check parameters routine for plant canopy model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_check_parameters

!   USE control_parameters,
!       ONLY: message_string

 END SUBROUTINE im_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
!-------------
!> Subroutine defining appropriate grid for netcdf variables.
!> It is called from subroutine netcdf.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )

    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z
    CHARACTER (LEN=*), INTENT(IN)  ::  var

    LOGICAL, INTENT(OUT)           ::  found


    found   = .TRUE.
!
!-- Check for the grid
    SELECT CASE ( TRIM( var ) )

       CASE ( 'im_hf_roof', 'im_hf_roof_waste' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'
!
!--    Heat fluxes at vertical walls are actually defined on stagged grid, i.e. xu, yv.
       CASE ( 'im_hf_wall_win', 'im_hf_wall_win_waste' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zu'

       CASE ( 'im_t_indoor_mean', 'im_t_indoor_roof', 'im_t_indoor_wall_win', 'im_t_indoor_wall', &
              'im_theta_10cm_roof', 'im_theta_10cm_wall' )
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'
    END SELECT

 END SUBROUTINE im_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine defining 3D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )

    USE indices

    USE kinds

    CHARACTER (LEN=*) ::  variable !<

    INTEGER(iwp) ::  av    !<
    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  k     !<
    INTEGER(iwp) ::  m     !<
    INTEGER(iwp) ::  nb    !< index of the building in the building data structure
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found !<

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf !<

    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )
!
!--     Output of indoor temperature. All grid points within the building are filled with values,
!--     while atmospheric grid points are set to _FillValues.
        CASE ( 'im_t_indoor_mean' )
           IF ( av == 0 ) THEN
              DO  i = nxl, nxr
                 DO  j = nys, nyn
                    IF ( building_id_f%var(j,i) /= building_id_f%fill )  THEN
!
!--                    Determine index of the building within the building data structure.
                       nb = MINLOC( ABS( buildings(:)%id - building_id_f%var(j,i) ), DIM=1 )
                       IF ( buildings(nb)%on_pe )  THEN
!
!--                       Write mean building temperature onto output array. Please note, in
!--                       contrast to many other loops in the output, the vertical bounds are
!--                       determined by the lowest and hightest vertical index occupied by the
!--                       building.
                          DO  k = buildings(nb)%kb_min, buildings(nb)%kb_max
                             local_pf(i,j,k) = buildings(nb)%t_in(k)
                          ENDDO
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
           ENDIF

        CASE ( 'im_hf_roof' )
           IF ( av == 0 )  THEN
              DO  m = 1, surf_usm%ns
                 i = surf_usm%i(m)
                 j = surf_usm%j(m)
                 k = surf_usm%k(m)
!
!--              Only values at upward-facing surfaces are taken.
                 local_pf(i,j,k) = MERGE( surf_usm%iwghf_eb(m), local_pf(i,j,k),                   &
                                          surf_usm%upward(m) )
              ENDDO
           ENDIF

        CASE ( 'im_hf_roof_waste' )
           IF ( av == 0 )  THEN
              DO m = 1, surf_usm%ns
                 i = surf_usm%i(m)
                 j = surf_usm%j(m)
                 k = surf_usm%k(m)
!
!--              Only values at upward-facing surfaces are taken.
                 local_pf(i,j,k) = MERGE( surf_usm%waste_heat(m), local_pf(i,j,k),                 &
                                          surf_usm%upward(m) )
              ENDDO
           ENDIF

       CASE ( 'im_hf_wall_win' )
           IF ( av == 0 )  THEN
              DO m = 1, surf_usm%ns
                 i = surf_usm%i(m)
                 j = surf_usm%j(m)
                 k = surf_usm%k(m)
!
!--              Only values at vertical surfaces are taken.
                 local_pf(i,j,k) = MERGE( surf_usm%iwghf_eb(m), local_pf(i,j,k),                   &
                                          surf_usm%northward(m)  .OR.  surf_usm%southward(m)  .OR. &
                                          surf_usm%eastward(m)   .OR.  surf_usm%westward(m) )
              ENDDO
           ENDIF

        CASE ( 'im_hf_wall_win_waste' )
           IF ( av == 0 )  THEN
              DO m = 1, surf_usm%ns
                 i = surf_usm%i(m)
                 j = surf_usm%j(m)
                 k = surf_usm%k(m)
!
!--              Only values at vertical surfaces are taken.
                 local_pf(i,j,k) = MERGE( surf_usm%waste_heat(m), local_pf(i,j,k),                 &
                                          surf_usm%northward(m)  .OR.  surf_usm%southward(m)  .OR. &
                                          surf_usm%eastward(m)   .OR.  surf_usm%westward(m) )
              ENDDO
           ENDIF

        CASE ( 'im_theta_10cm_wall' )
           IF ( av == 0 )  THEN
              DO m = 1, surf_usm%ns
                 i = surf_usm%i(m)
                 j = surf_usm%j(m)
                 k = surf_usm%k(m)
!
!--              Only values at vertical surfaces are taken.
                 local_pf(i,j,k) = MERGE( surf_usm%pt_10cm(m), local_pf(i,j,k),                    &
                                          surf_usm%northward(m)  .OR.  surf_usm%southward(m)  .OR. &
                                          surf_usm%eastward(m)   .OR.  surf_usm%westward(m) )
              ENDDO
           ENDIF

        CASE ( 'im_theta_10cm_roof' )
           IF ( av == 0 )  THEN
              DO m = 1, surf_usm%ns
                 i = surf_usm%i(m)
                 j = surf_usm%j(m)
                 k = surf_usm%k(m)
!
!--              Only values at vertical surfaces are taken.
                 local_pf(i,j,k) = MERGE( surf_usm%pt_10cm(m), local_pf(i,j,k), surf_usm%upward(m) )
              ENDDO
           ENDIF

!
!< NOTE im_t_indoor_roof and im_t_indoor_wall_win not work yet

!         CASE ( 'im_t_indoor_roof' )
!            IF ( av == 0 )  THEN
!               DO  m = 1, surf_usm_h%ns
!                   i = surf_usm_h%i(m) !+ surf_usm_h%ioff(m)
!                   j = surf_usm_h%j(m) !+ surf_usm_h%joff(m)
!                   k = surf_usm_h%k(m) !+ surf_usm_h%koff(m)
!                   local_pf(i,j,k) = surf_usm_h%t_indoor(m)
!               ENDDO
!            ENDIF
!
!         CASE ( 'im_t_indoor_wall_win' )
!            IF ( av == 0 )  THEN
!               DO l = 0, 3
!                  DO m = 1, surf_usm_v(l)%ns
!                     i = surf_usm_v(l)%i(m) !+ surf_usm_v(l)%ioff(m)
!                     j = surf_usm_v(l)%j(m) !+ surf_usm_v(l)%joff(m)
!                     k = surf_usm_v(l)%k(m) !+ surf_usm_v(l)%koff(m)
!                     local_pf(i,j,k) = surf_usm_v(l)%t_indoor(m)
!                  ENDDO
!               ENDDO
!            ENDIF

        CASE DEFAULT
           found = .FALSE.

    END SELECT

 END SUBROUTINE im_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &indoor_parameters for indoor model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_parin

    USE control_parameters,                                                                        &
        ONLY:  indoor_model


    CHARACTER(LEN=100) ::  line  !< string containing current line of file PARIN

    INTEGER(iwp)  ::  io_status  !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /indoor_parameters/  indoor_during_spinup,                                            &
                                  initial_indoor_temperature,                                      &
                                  switch_off_module

!
!-- Move to the beginning of the namelist file and try to find and read the namelist.
    REWIND( 11 )
    READ( 11, indoor_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    indoor_parameters namelist was found and read correctly. Set flag that indicates that the
!--    indoor model is switched on.
       IF ( .NOT. switch_off_module )  indoor_model = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    indoor_parameters namelist was found, but contained errors. Print an error message including
!--    the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'indoor_parameters', line )

    ENDIF

!
!--    Activate spinup (maybe later
!        IF ( spinup_time > 0.0_wp )  THEN
!           coupling_start_time = spinup_time
!           end_time = end_time + spinup_time
!           IF ( spinup_pt_mean == 9999999.9_wp )  THEN
!              spinup_pt_mean = pt_surface
!           ENDIF
!           spinup = .TRUE.
!        ENDIF

 END SUBROUTINE im_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!> Soubroutine reads t_surf and t_wall.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxr_on_file, nynf, nyn_on_file,    &
                              nysf, nysc, nys_on_file, found )


    USE control_parameters,                                                                        &
        ONLY: length,                                                                              &
              restart_string

    IMPLICIT NONE

    INTEGER(iwp) ::  k                 !< running index over previous input files covering current local domain
    INTEGER(iwp) ::  nxlc              !< index of left boundary on current subdomain
    INTEGER(iwp) ::  nxlf              !< index of left boundary on former subdomain
    INTEGER(iwp) ::  nxl_on_file       !< index of left boundary on former local domain
    INTEGER(iwp) ::  nxrf              !< index of right boundary on former subdomain
    INTEGER(iwp) ::  nxr_on_file       !< index of right boundary on former local domain
    INTEGER(iwp) ::  nynf              !< index of north boundary on former subdomain
    INTEGER(iwp) ::  nyn_on_file       !< index of north boundary on former local domain
    INTEGER(iwp) ::  nysc              !< index of south boundary on current subdomain
    INTEGER(iwp) ::  nysf              !< index of south boundary on former subdomain
    INTEGER(iwp) ::  nys_on_file       !< index of south boundary on former local domain
    INTEGER(iwp) ::  ns_on_file_usm_im !< number of surface elements (urban type) on file
!
!-- Note, the save attribute in the following array declaration is necessary, in order to keep the
!-- number of urban surface elements on file during rrd_local calls.
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  end_index_on_file    !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE, SAVE ::  start_index_on_file  !<

    LOGICAL, INTENT(OUT)  ::  found  !<

    TYPE( surf_type_1d_im ), SAVE ::  tmp_surf   !< temporary variable to read surface data

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

       CASE ( 'ns_on_file_usm_im')
          IF ( k == 1 )  THEN
             READ ( 13 ) ns_on_file_usm_im
!
!--          In case of changing mpi topology, this routine could be called more than once.
!--          Hence, arrays need to be deallocated before allocated again.
             IF ( ALLOCATED( tmp_surf%val ) )  DEALLOCATE( tmp_surf%val )
!
!--          Allocate temporary arrays for reading data on file. Note, the size of allocated surface
!--          elements do not necessarily need to match the size of present surface elements on
!--          current processor, as the number of processors between restarts can change.
             ALLOCATE( tmp_surf%val(1:ns_on_file_usm_im) )
          ENDIF

       CASE ( 'im_start_index' )
          IF ( k == 1 )  THEN

             IF ( ALLOCATED( start_index_on_file ) )  DEALLOCATE( start_index_on_file )

             ALLOCATE ( start_index_on_file(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) )

             READ ( 13 )  start_index_on_file

          ENDIF

       CASE ( 'im_end_index' )
          IF ( k == 1 )  THEN

             IF ( ALLOCATED( end_index_on_file ) )  DEALLOCATE( end_index_on_file )

             ALLOCATE ( end_index_on_file(nys_on_file:nyn_on_file,nxl_on_file:nxr_on_file) )

             READ ( 13 )  end_index_on_file

          ENDIF

       CASE ( 'waste_heat' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm%waste_heat ) )                                        &
                ALLOCATE( surf_usm%waste_heat(1:surf_usm%ns) )
             READ ( 13 )  tmp_surf%val
          ENDIF
          CALL surface_restore_elements( surf_usm%waste_heat, tmp_surf%val,                        &
                                         surf_usm%start_index, start_index_on_file,                &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE ( 't_prev' )
          IF ( k == 1 )  THEN
             IF ( .NOT.  ALLOCATED( surf_usm%t_prev ) )  ALLOCATE( surf_usm%t_prev(1:surf_usm%ns) )
             READ ( 13 )  tmp_surf%val
          ENDIF
          CALL surface_restore_elements( surf_usm%t_prev, tmp_surf%val,                            &
                                         surf_usm%start_index, start_index_on_file,                &
                                         end_index_on_file, nxlc, nysc, nxlf, nxrf, nysf, nynf,    &
                                         nys_on_file, nyn_on_file, nxl_on_file,nxr_on_file )

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE im_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!> Soubroutine reads t_surf and t_wall.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_rrd_local_mpi

    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_end_index
    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_start_index

    LOGICAL ::  array_found  !<
    LOGICAL ::  data_to_read !< dummy variable


    CALL rd_mpi_io_check_array( 'usm_global_start', found = array_found )
    IF ( array_found )  CALL rrd_mpi_io( 'usm_global_start', global_start_index )

    CALL rd_mpi_io_check_array( 'usm_global_end', found = array_found )
    IF ( array_found )  CALL rrd_mpi_io( 'usm_global_end', global_end_index )
!
!-- Check if data input for surface-type variables is required. Note, only invoke routine if USM
!-- surface restart data is on file. In case of cyclic fill initialization this is not necessarily
!-- guaranteed. To check this use the array_found control flag.
    IF ( array_found )  THEN
       CALL rd_mpi_io_surface_filetypes( surf_usm%start_index, surf_usm%end_index,                 &
                                         data_to_read, global_start_index, global_end_index )
    ELSE
       data_to_read = .FALSE.
    ENDIF

    IF ( data_to_read )  THEN
       CALL rd_mpi_io_check_array( 'waste_heat', found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_usm%waste_heat ) )                                            &
             ALLOCATE( surf_usm%waste_heat(1:surf_usm%ns) )
          CALL rrd_mpi_io_surface( 'waste_heat', surf_usm%waste_heat )
       ENDIF

       CALL rd_mpi_io_check_array( 't_prev', found = array_found )
       IF ( array_found )  THEN
          IF ( .NOT. ALLOCATED( surf_usm%t_prev ) )  ALLOCATE( surf_usm%t_prev(1:surf_usm%ns) )
          CALL rrd_mpi_io_surface( 't_prev', surf_usm%t_prev )
       ENDIF
    ENDIF

 END SUBROUTINE im_rrd_local_mpi

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Subroutine to writes restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE im_wrd_local

    IMPLICIT NONE

    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_end_index    !< end index for surface data (MPI-IO)
    INTEGER(idp), DIMENSION(nys:nyn,nxl:nxr) ::  global_start_index  !< start index for surface data (MPI-IO)

    LOGICAL  ::  surface_data_to_write  !< switch for MPI-I/O if PE has surface data to write


    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

       CALL wrd_write_string( 'ns_on_file_usm_im' )
       WRITE ( 14 )  surf_usm%ns

       CALL wrd_write_string( 'waste_heat' )
       WRITE ( 14 )  surf_usm%waste_heat

       CALL wrd_write_string( 't_prev' )
       WRITE ( 14 )  surf_usm%t_prev

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
!
!--    There is no information about the PE-grid necessary because the restart files consists of the
!--    whole domain. Therefore, ns_on_file_usm are not used with MPI-IO.
       CALL rd_mpi_io_surface_filetypes( surf_usm%start_index, surf_usm%end_index,                 &
                                         surface_data_to_write, global_start_index,                &
                                         global_end_index )

       CALL wrd_mpi_io( 'im_global_start', global_start_index )
       CALL wrd_mpi_io( 'im_global_end', global_end_index )

       IF ( surface_data_to_write )  THEN
          CALL wrd_mpi_io_surface( 'waste_heat', surf_usm%waste_heat )
          CALL wrd_mpi_io_surface( 't_prev',     surf_usm%t_prev     )
       ENDIF

    ENDIF

 END SUBROUTINE im_wrd_local

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculates temperature near surface (10 cm) for indoor model.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE calc_pt_10cm

    INTEGER(iwp) ::  m  !< running index for surface elements

    DO  m = 1, surf_usm%ns
!
!--    For horizontal upward-facing surfaces 10-cm temperature can be calculated
!--    using MOST. At vertical surfaces 10-cm temperature is not calculated via MOST
!--    but set to the grid-cell temperature.
       surf_usm%pt_10cm(m) = MERGE(                                                                &
                             surf_usm%pt_surface(m) + surf_usm%ts(m) / kappa                       &
                          * ( LOG( 0.1_wp /  surf_usm%z0h(m) ) - psi_h( 0.1_wp / surf_usm%ol(m) )  &
                              + psi_h( surf_usm%z0h(m) / surf_usm%ol(m) ) ),                       &
                              pt(surf_usm%k(m),surf_usm%j(m),surf_usm%i(m)),                       &
                              surf_usm%upward(m) )
    ENDDO

 END SUBROUTINE calc_pt_10cm


END MODULE indoor_model_mod
