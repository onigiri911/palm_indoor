! !> @file palm.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2020 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
!
! Former revisions:
! -----------------
! $Id: palm.f90 4539 2020-05-18 14:05:17Z raasch $
! log point name changed
! 
! 4535 2020-05-15 12:07:23Z raasch
! bugfix for restart data format query
! 
! 4496 2020-04-15 08:37:26Z raasch
! bugfix: coupling character added to restart output filename
! 
! 4495 2020-04-13 20:11:20Z raasch
! restart data handling with MPI-IO added
! 
! 4457 2020-03-11 14:20:43Z raasch
! use statement for exchange horiz added
! 
! 4444 2020-03-05 15:59:50Z raasch
! bugfix: cpp-directives for serial mode added
! 
! 4414 2020-02-19 20:16:04Z suehring
! Call to module_interface_init_numerics
! 
! 4400 2020-02-10 20:32:41Z suehring
! Add interface to initialize data output with dom
! 
! 4360 2020-01-07 11:25:50Z suehring
! implement new palm_date_time_mod
! 
! 4094 2019-07-12 09:24:21Z gronemeier
! Corrected "Former revisions" section
! 
! 4039 2019-06-18 10:32:41Z suehring
! Rename subroutines in module for diagnostic quantities
! 
! 4017 2019-06-06 12:16:46Z schwenkel
! new module for calculation and output of diagnostic quantities added
! 
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction 
! of additional debug messages
! 
! 3761 2019-02-25 15:31:42Z raasch
! unused variable removed
! 
! 3719 2019-02-06 13:10:18Z kanani
! Included cpu measurement for wall/soil spinup
! 
! 3703 2019-01-29 16:43:53Z knoop
! Some interface calls moved to module_interface + cleanup
! 
! 3648 2019-01-02 16:35:46Z suehring
! Rename subroutines for surface-data output
!
! Revision 1.1  1997/07/24 11:23:35  raasch
! Initial revision
!
!
! Description:
! ------------
!> Large-Eddy Simulation (LES) model for atmospheric and oceanic boundary-layer
!> flows
!> see the PALM homepage https://palm-model.org for further information
!------------------------------------------------------------------------------!
 PROGRAM palm
 

    USE arrays_3d

#if defined( __parallel )
    USE bulk_cloud_model_mod,                                                  &
        ONLY: bulk_cloud_model, microphysics_morrison, microphysics_seifert
#endif

    USE control_parameters,                                                    &
        ONLY:  coupling_char, do2d_at_begin, do3d_at_begin, io_blocks,         &
               io_group, message_string, restart_data_format_output, runnr, simulated_time_chr, spinup,   &
               time_since_reference_point, user_interface_current_revision,    &
               user_interface_required_revision, version, write_binary, pt_surface

#if defined( __parallel )
    USE control_parameters,                                                    &
        ONLY:  child_domain, constant_diffusion, humidity,                     &
               initializing_actions, neutral, passive_scalar
#endif

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, cpu_statistics

#if defined( __parallel )
    USE cpulog,                                                                &
        ONLY:  log_point_s
#endif

    USE diagnostic_output_quantities_mod,                                      &
        ONLY:  doq_calculate

#if defined( __parallel )
    USE exchange_horiz_mod,                                                    &
        ONLY:  exchange_horiz

    USE indices,                                                               &
        ONLY:  nbgp
#endif
    
!-- COVID-19 specific code 
!       
    USE indices,                                                               &
        ONLY: nxl, nxr, nys, nyn, nzb, nzt, wall_flags_total_0
!
!-- COVID-19 specific code ends   

    USE kinds

    USE module_interface,                                                      &
        ONLY:  module_interface_init_numerics,                                 &
               module_interface_init_output,                                   &
               module_interface_last_actions


    USE multi_agent_system_mod,                                                &
        ONLY:  agents_active, mas_last_actions

    USE netcdf_data_input_mod,                                                 &
        ONLY:  netcdf_data_input_inquire_file, netcdf_data_input_init,         &
               netcdf_data_input_surface_data, netcdf_data_input_topo

    USE pegrid

#if defined( __parallel )
    USE pmc_particle_interface,                                                &
        ONLY: pmcp_g_alloc_win

    USE pmc_interface,                                                         &
        ONLY:  nested_run, pmci_child_initialize, pmci_init,                   &
               pmci_modelconfiguration, pmci_parent_initialize
#endif

    USE restart_data_mpi_io_mod,                                               &
        ONLY:  rd_mpi_io_close, rd_mpi_io_open

    USE surface_data_output_mod,                                               &
        ONLY:  surface_data_output_last_action
!
!-- COVID-19 specific code 
    USE user,                                                                  &
        ONLY:  reset_scalar_at_restart
!
!-- COVID-19 specific code ends

    USE write_restart_data_mod,                                                &
        ONLY:  wrd_global, wrd_local

#if defined( __parallel )  &&  defined( _OPENACC )
    USE openacc
#endif


    IMPLICIT NONE

!
!-- Local variables
    CHARACTER(LEN=9) ::  time_to_string  !<
    INTEGER(iwp)     ::  i               !< loop counter for blocked I/O
    
!-- COVID-19 specific code 
!
    INTEGER(iwp)     ::  j               !< loop counter for pt override
    INTEGER(iwp)     ::  k               !< loop counter for pt override 
!
!-- COVID-19 specific code ends
       
#if defined( __parallel) && defined( _OPENACC )
    INTEGER(iwp)     :: local_comm       !< local communicator (shared memory)
    INTEGER(iwp)     :: local_num_procs  !< local number of processes
    INTEGER(iwp)     :: local_id         !< local id
    INTEGER(acc_device_kind) :: device_type !< device type for OpenACC
    INTEGER(iwp)     ::  num_devices     !< number of devices visible to OpenACC
    INTEGER(iwp)     ::  my_device       !< device used by this process
#endif

    version = 'PALM 6.0'
    user_interface_required_revision = 'r4495'

#if defined( __parallel )
!
!-- MPI initialisation. comm2d is preliminary set, because
!-- it will be defined in init_pegrid but is used before in cpu_log.
    CALL MPI_INIT( ierr )

!
!-- Initialize the coupling for nested-domain runs
!-- comm_palm is the communicator which includes all PEs (MPI processes)
!-- available for this (nested) model. If it is not a nested run, comm_palm
!-- is returned as MPI_COMM_WORLD
    CALL cpu_log( log_point_s(70), 'pmci_init', 'start' )
    CALL pmci_init( comm_palm )
    CALL cpu_log( log_point_s(70), 'pmci_init', 'stop' )
    comm2d = comm_palm
!
!-- Get the (preliminary) number of MPI processes and the local PE-id (in case
!-- of a further communicator splitting in init_coupling, these numbers will
!-- be changed in init_pegrid).
    IF ( nested_run )  THEN

       CALL MPI_COMM_SIZE( comm_palm, numprocs, ierr )
       CALL MPI_COMM_RANK( comm_palm, myid, ierr )

    ELSE

       CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
       CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
!
!--    Initialize PE topology in case of coupled atmosphere-ocean runs (comm_palm
!--    will be splitted in init_coupling)
       CALL init_coupling
    ENDIF

#ifdef _OPENACC
!
!-- Select OpenACC device to use in this process. For this find out how many
!-- neighbors there are running on the same node and which id this process is.
    IF ( nested_run )  THEN
       CALL MPI_COMM_SPLIT_TYPE( comm_palm, MPI_COMM_TYPE_SHARED, 0,           &
                                 MPI_INFO_NULL, local_comm, ierr )
    ELSE
       CALL MPI_COMM_SPLIT_TYPE( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,      &
                                 MPI_INFO_NULL, local_comm, ierr )
    ENDIF
    CALL MPI_COMM_SIZE( local_comm, local_num_procs, ierr )
    CALL MPI_COMM_RANK( local_comm, local_id, ierr )

!
!-- This loop including the barrier is a workaround for PGI compiler versions
!-- up to and including 18.4. Later releases are able to select their GPUs in
!-- parallel, without running into spurious errors.
    DO i = 0, local_num_procs-1
       CALL MPI_BARRIER( local_comm, ierr )

       IF ( i == local_id )  THEN
          device_type = acc_get_device_type()
          num_devices = acc_get_num_devices( device_type )
          my_device = MOD( local_id, num_devices )
          CALL acc_set_device_num( my_device, device_type )
       ENDIF
    ENDDO

    CALL MPI_COMM_FREE( local_comm, ierr )
#endif
#endif

!
!-- Initialize measuring of the CPU-time remaining to the run
    CALL local_tremain_ini

!
!-- Start of total CPU time measuring.
    CALL cpu_log( log_point(1), 'total', 'start' )
    CALL cpu_log( log_point(2), 'initialisation', 'start' )

!
!-- Open a file for debug output
    WRITE (myid_char,'(''_'',I6.6)')  myid
    OPEN( 9, FILE='DEBUG'//TRIM( coupling_char )//myid_char, FORM='FORMATTED' )

!
!-- Initialize dvrp logging. Also, one PE maybe split from the global
!-- communicator for doing the dvrp output. In that case, the number of
!-- PEs available for PALM is reduced by one and communicator comm_palm
!-- is changed respectively.
#if defined( __parallel )
    CALL MPI_COMM_RANK( comm_palm, myid, ierr )
#endif

!
!-- Read control parameters from NAMELIST files and read environment-variables
    CALL parin

!
!-- Check for the user's interface version
    IF ( user_interface_current_revision /= user_interface_required_revision )  &
    THEN
       message_string = 'current user-interface revision "' //                  &
                        TRIM( user_interface_current_revision ) // '" does ' // &
                        'not match the required revision ' //                   &
                        TRIM( user_interface_required_revision )
        CALL message( 'palm', 'PA0169', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine processor topology and local array indices
    CALL init_pegrid
!
!-- Check if input file according to input-data standard exists
    CALL netcdf_data_input_inquire_file
!
!-- Read topography input data if required. This is required before the 
!-- numerical grid is finally created in init_grid
    CALL netcdf_data_input_topo  
!
!-- Generate grid parameters, initialize generic topography and further process
!-- topography information if required
    CALL init_grid
!
!-- Initialize boundary conditions and numerics such as the multigrid solver or 
!-- the advection routine
    CALL module_interface_init_numerics
!
!-- Read global attributes if available.  
    CALL netcdf_data_input_init 
!
!-- Read surface classification data, e.g. vegetation and soil types, water 
!-- surfaces, etc., if available. Some of these data is required before 
!-- check parameters is invoked.     
    CALL netcdf_data_input_surface_data
!
!-- Check control parameters and deduce further quantities
    CALL check_parameters

    CALL init_3d_model
!                                                
!-- COVID-19 specific code
!
!-- Resetting of passive scalar after restart
    IF ( TRIM( initializing_actions ) == 'read_restart_data' .OR. TRIM( initializing_actions ) == 'cyclic_fill' ) THEN
       IF ( reset_scalar_at_restart )  THEN
          s(:,:,:) = 0.0_wp
       ENDIF
    ENDIF
    
!    IF ( TRIM( initializing_actions ) == 'cyclic_fill' ) THEN
!       IF ( reset_scalar_at_restart )  THEN
!          
!          DO  i = nxl, nxr
!             DO  j = nys, nyn
!                DO  k = nzb+1, nzt
!		
!		   IF( BTEST(wall_flags_total_0(k,j,i), 0) ) THEN
!                      pt(k,j,i) = pt_surface
!		   ENDIF
!		   
!                ENDDO
!             ENDDO
!          ENDDO
!	  
!       ENDIF
!    ENDIF  
!
!-- COVID-19 specific code ends

    CALL module_interface_init_output

#if defined( __parallel )
!
!-- Coupling protocol setup for nested-domain runs
    IF ( nested_run )  THEN
       CALL pmci_modelconfiguration
!
!--    Receive and interpolate initial data on children.
!--    Child initialization must be made first if the model is both child and
!--    parent if necessary
       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN
          CALL pmci_child_initialize
!
!--       Send initial condition data from parent to children
          CALL pmci_parent_initialize
!
!--       Exchange_horiz is needed after the nest initialization
          IF ( child_domain )  THEN
             CALL exchange_horiz( u, nbgp )
             CALL exchange_horiz( v, nbgp )
             CALL exchange_horiz( w, nbgp )
             IF ( .NOT. neutral )  THEN
                CALL exchange_horiz( pt, nbgp )
             ENDIF
             IF ( .NOT. constant_diffusion )  CALL exchange_horiz( e, nbgp )
             IF ( humidity )  THEN
                CALL exchange_horiz( q, nbgp )
                IF ( bulk_cloud_model  .AND.  microphysics_morrison )  THEN
                  CALL exchange_horiz( qc, nbgp )
                  CALL exchange_horiz( nc, nbgp )
                ENDIF
                IF ( bulk_cloud_model  .AND.  microphysics_seifert )  THEN
                   CALL exchange_horiz( qr, nbgp ) 
                   CALL exchange_horiz( nr, nbgp )
                ENDIF
             ENDIF
             IF ( passive_scalar )  CALL exchange_horiz( s, nbgp )
          ENDIF
       ENDIF

       CALL pmcp_g_alloc_win                    ! Must be called after pmci_child_initialize and pmci_parent_initialize
    ENDIF
#endif

!
!-- Output of program header
    IF ( myid == 0 )  CALL header

    CALL cpu_log( log_point(2), 'initialisation', 'stop' )

!
!-- Integration of the non-atmospheric equations (land surface model, urban
!-- surface model)
    IF ( spinup )  THEN
       CALL cpu_log( log_point(41), 'wall/soil spinup', 'start' )
       CALL time_integration_spinup
       CALL cpu_log( log_point(41), 'wall/soil spinup', 'stop' )
    ENDIF

!
!-- Set start time in format hh:mm:ss
    simulated_time_chr = time_to_string( time_since_reference_point )

!
!-- If required, output of initial arrays
    IF ( do2d_at_begin )  THEN
       CALL doq_calculate    !TODO, will be called twice

       CALL data_output_2d( 'xy', 0 )
       CALL data_output_2d( 'xz', 0 )
       CALL data_output_2d( 'yz', 0 )
    ENDIF

    IF ( do3d_at_begin )  THEN
       CALL doq_calculate    !TODO, will be called twice

       CALL data_output_3d( 0 )
    ENDIF

!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

!
!-- If required, write binary data for restart runs
    IF ( write_binary )  THEN

       CALL cpu_log( log_point(22), 'write-restart-data', 'start' )

       CALL location_message( 'writing restart data', 'start' )

       IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

          DO  i = 0, io_blocks-1
             IF ( i == io_group )  THEN

!
!--             Open binary file
                CALL check_open( 14 )
!
!--             Write control parameters and other global variables for restart.
                IF ( myid == 0 )  CALL wrd_global
!
!--             Write processor specific flow field data for restart runs
                CALL wrd_local
!
!--             Close binary file
                CALL close_file( 14 )

             ENDIF
#if defined( __parallel )
             CALL MPI_BARRIER( comm2d, ierr )
#endif
          ENDDO

       ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN
!
!--       Open MPI-IO restart file
          CALL rd_mpi_io_open( 'write', 'BINOUT' // TRIM( coupling_char ) )
!
!--       Write control parameters and other global variables for restart.
          CALL wrd_global
!
!--       Write processor specific flow field data for restart runs
          CALL wrd_local
!
!--       Close restart File
          CALL rd_mpi_io_close

       ENDIF

       CALL location_message( 'writing restart data', 'finished' )

       CALL cpu_log( log_point(22), 'write-restart-data', 'stop' )
       
    ENDIF
!
!-- Last actions for surface output, for instantaneous and time-averaged data
    CALL surface_data_output_last_action( 0 )
    CALL surface_data_output_last_action( 1 )

!
!-- If required, repeat output of header including the required CPU-time
    IF ( myid == 0 )  CALL header
!
!-- Perform module specific last actions
    CALL cpu_log( log_point(4), 'last actions', 'start' )

    IF ( myid == 0 .AND. agents_active ) CALL mas_last_actions ! ToDo: move to module_interface

    CALL module_interface_last_actions

    CALL cpu_log( log_point(4), 'last actions', 'stop' )

!
!-- Close files
    CALL close_file( 0 )

!
!-- Write run number to file (used by palmrun to create unified cycle numbers
!-- for output files
    IF ( myid == 0  .AND.  runnr > 0 )  THEN
       OPEN( 90, FILE='RUN_NUMBER', FORM='FORMATTED' )
       WRITE( 90, '(I4)' )  runnr
       CLOSE( 90 )
    ENDIF

!
!-- Take final CPU-time for CPU-time analysis
    CALL cpu_log( log_point(1), 'total', 'stop' )
    CALL cpu_statistics

#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif

 END PROGRAM palm
