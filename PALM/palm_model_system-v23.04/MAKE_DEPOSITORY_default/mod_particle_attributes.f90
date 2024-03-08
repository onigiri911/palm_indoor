!> @file mod_particle_attributes.f90
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
! Copyright 1997-2020 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
!
!
! Former revisions:
! -----------------
! $Id: mod_particle_attributes.f90 4778 2020-11-09 13:40:05Z raasch $
! variables for particle output renamed, id counter removed from particle type
!
! 4677 2020-09-14 07:55:28Z raasch
! file re-formatted to follow the PALM coding standard
!
! 4628 2020-07-29 07:23:03Z raasch
! extensions required for MPI-I/O of particle data to restart files
!
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
!
! 4043 2019-06-18 16:59:00Z schwenkel
! Remove min_nr_particle
!
! 4017 2019-06-06 12:16:46Z schwenkel
! interoperable C datatypes introduced in particle type to avoid compiler warnings
!
! 3720 2019-02-06 13:19:55Z knoop
! time_prel replaced by last_particle_release_time
!
! 1359 2014-04-11 17:15:14Z hoffmann
! new module containing all particle related variables
! -dt_sort_particles
!
! Description:
! ------------
!> Definition of variables used to compute particle transport
!--------------------------------------------------------------------------------------------------!
 MODULE particle_attributes

    USE, INTRINSIC ::  ISO_C_BINDING

    USE control_parameters,                                                                        &
        ONLY: varnamelength

    USE kinds

    INTEGER(iwp), PARAMETER ::  max_number_of_particle_groups = 10 !< maximum allowed number of particle groups

    CHARACTER(LEN=64)                           ::  pts_id_file = ''  !< namelist parameter
    CHARACTER(LEN=varnamelength), DIMENSION(50) ::  data_output_pts = ''    !< namelist parameter

    INTEGER(iwp) ::  dissipation_classes = 10                     !< namelist parameter (see documentation)
    INTEGER(iwp) ::  ibc_par_b                                    !< particle bottom boundary condition dummy
    INTEGER(iwp) ::  ibc_par_lr                                   !< particle left/right boundary condition dummy
    INTEGER(iwp) ::  ibc_par_ns                                   !< particle north/south boundary condition dummy
    INTEGER(iwp) ::  ibc_par_t                                    !< particle top boundary condition dummy
    INTEGER(iwp) ::  number_of_output_particles = 0               !< number of output particles
    INTEGER(iwp) ::  number_of_particles = 0                      !< number of particles for each grid box (3d array is saved on
                                                                  !< prt_count)
    INTEGER(iwp) ::  number_of_particle_groups = 1                !< namelist parameter (see documentation)
    INTEGER(iwp) ::  pts_increment = 1                            !< increment of particles in output file

    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  prt_count  !< 3d array of number of particles of every grid box

    LOGICAL ::  particle_advection = .FALSE.              !< parameter to steer the advection of particles
    LOGICAL ::  unlimited_dimension = .TRUE.              !< umlimited dimension for particle output
    LOGICAL ::  use_sgs_for_particles = .FALSE.           !< namelist parameter (see documentation)
    LOGICAL ::  wang_kernel = .FALSE.                     !< flag for collision kernel

    REAL(wp) ::  alloc_factor = 20.0_wp                   !< namelist parameter (see documentation)
    REAL(wp) ::  oversize = 100.0_wp                      !< reserve spare particles in output file (in % relative to initial
                                                          !< number)
    REAL(wp) ::  particle_advection_start = 0.0_wp        !< namelist parameter (see documentation)
    REAL(wp) ::  pts_percentage = 100.0_wp                !< percentage of particles in output file

    TYPE, PUBLIC ::  particle_type
        REAL(wp)     ::  aux1             !< auxiliary multi-purpose feature
        REAL(wp)     ::  aux2             !< auxiliary multi-purpose feature
        REAL(wp)     ::  radius           !< radius of particle
        REAL(wp)     ::  age              !< age of particle
        REAL(wp)     ::  age_m            !<
        REAL(wp)     ::  dt_sum           !<
        REAL(wp)     ::  e_m              !< interpolated sgs tke
        REAL(wp)     ::  origin_x         !< origin x-position of particle (changed cyclic bc)
        REAL(wp)     ::  origin_y         !< origin y-position of particle (changed cyclic bc)
        REAL(wp)     ::  origin_z         !< origin z-position of particle (changed cyclic bc)
        REAL(wp)     ::  rvar1            !<
        REAL(wp)     ::  rvar2            !<
        REAL(wp)     ::  rvar3            !<
        REAL(wp)     ::  speed_x          !< speed of particle in x
        REAL(wp)     ::  speed_y          !< speed of particle in y
        REAL(wp)     ::  speed_z          !< speed of particle in z
        REAL(wp)     ::  weight_factor    !< weighting factor
        REAL(wp)     ::  x                !< x-position
        REAL(wp)     ::  y                !< y-position
        REAL(wp)     ::  z                !< z-position
        INTEGER(iwp) ::  class            !< radius class needed for collision
        INTEGER(iwp) ::  group            !< number of particle group
        INTEGER(idp) ::  id               !< particle ID (64 bit integer)
        LOGICAL      ::  particle_mask    !< if this parameter is set to false the particle will be deleted
        INTEGER(iwp) ::  block_nr         !< number for sorting (removable?)
        INTEGER(iwp) ::  particle_nr=-1   !< particle number for particle IO (increment one
    END TYPE particle_type

    TYPE(particle_type), DIMENSION(:), POINTER ::  particles      !< Particle array for this grid cell
    TYPE(particle_type)                        ::  zero_particle  !< zero particle to avoid weird things

    TYPE particle_groups_type
        SEQUENCE
        REAL(wp) ::  density_ratio  !< density ratio of the fluid and the particles
        REAL(wp) ::  radius         !< radius of particle
        REAL(wp) ::  exp_arg        !< exponential term of particle inertia
        REAL(wp) ::  exp_term       !< exponential term of particle inertia
    END TYPE particle_groups_type

    TYPE(particle_groups_type), DIMENSION(max_number_of_particle_groups) ::  particle_groups

    TYPE  grid_particle_def
        INTEGER(iwp), DIMENSION(0:7)               ::  start_index     !< start particle index for current block
        INTEGER(iwp), DIMENSION(0:7)               ::  end_index       !< end particle index for current block
        LOGICAL                                    ::  time_loop_done  !< timestep loop for particle advection
        TYPE(particle_type), POINTER, DIMENSION(:) ::  particles       !< Particle array for this grid cell
    END TYPE grid_particle_def

    TYPE(grid_particle_def), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  grid_particles

    TYPE block_offset_def          !<
        INTEGER(iwp) ::  i_off     !<
        INTEGER(iwp) ::  j_off     !<
        INTEGER(iwp) ::  k_off     !<
    END TYPE block_offset_def

    TYPE(block_offset_def), DIMENSION(0:7)         ::  block_offset

    SAVE


 END MODULE particle_attributes
