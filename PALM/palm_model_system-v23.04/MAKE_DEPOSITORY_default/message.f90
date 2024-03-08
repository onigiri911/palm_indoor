!> @file message.f90
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
! -----------------
!
!
! Former revisions:
! -----------------
! $Id: message.f90 4677 2020-09-14 07:55:28Z raasch $
! file re-formatted to follow the PALM coding standard
!
! 4580 2020-06-29 07:54:21Z raasch
! bugfix for aborts in case of nested runs
!
! 4578 2020-06-25 15:43:32Z gronemeier
! bugfix : do not save input values from last call of routines debug_message and location_message
! changes: layout changes according to PALM coding standards
!
! 4536 2020-05-17 17:24:13Z raasch
! location message format changed
!
! 4360 2020-01-07 11:25:50Z suehring
! Corrected "Former revisions" section
!
! 4097 2019-07-15 11:59:11Z suehring
! Avoid overlong lines - limit is 132 characters per line
!
! 3987 2019-05-22 09:52:13Z kanani
! Improved formatting of job logfile output,
! changed output of DEBUG file
!
! 3885 2019-04-11 11:29:34Z kanani
! Changes related to global restructuring of location messages and introduction of additional debug
! messages
!
! 3655 2019-01-07 16:51:22Z knoop
! Minor formating changes
!
! 213 2008-11-13 10:26:18Z raasch
! Initial revision
!
! Description:
! ------------
!> Handling of the different kinds of messages.
!> Meaning of formal parameters:
!> requested_action: 0 - continue, 1 - abort by stop, 2 - abort by mpi_abort
!>                   3 - abort by mpi_abort using MPI_COMM_WORLD
!> message_level: 0 - informative, 1 - warning, 2 - error
!> output_on_pe: -1 - all, else - output on specified PE
!> file_id: 6 - stdout (*)
!> flush_file: 0 - no action, 1 - flush the respective output buffer
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE message( routine_name, message_identifier, requested_action, message_level,            &
                     output_on_pe, file_id, flush_file )

    USE control_parameters,                                                                        &
        ONLY:  abort_mode, message_string

    USE kinds

    USE pegrid

    USE pmc_interface,                                                                             &
        ONLY:  cpl_id, nested_run

    IMPLICIT NONE

    CHARACTER(LEN=6)   ::  message_identifier            !<
    CHARACTER(LEN=20)  ::  nest_string                   !< nest id information
    CHARACTER(LEN=*)   ::  routine_name                  !<
    CHARACTER(LEN=200) ::  header_string                 !<
    CHARACTER(LEN=200) ::  header_string_2               !< for message ID and routine name
    CHARACTER(LEN=200) ::  information_string_1          !<
    CHARACTER(LEN=200) ::  information_string_2          !<

    INTEGER(iwp) ::  file_id                             !<
    INTEGER(iwp) ::  flush_file                          !<
    INTEGER(iwp) ::  i                                   !<
    INTEGER(iwp) ::  message_level                       !<
    INTEGER(iwp) ::  output_on_pe                        !<
    INTEGER(iwp) ::  requested_action                    !<

    LOGICAL ::  do_output                                !<
    LOGICAL ::  pe_out_of_range                          !<


    do_output       = .FALSE.
    pe_out_of_range = .FALSE.

!
!-- In case of nested runs create the nest id informations
    IF ( nested_run )  THEN
       WRITE( nest_string, '(1X,A,I2.2)' )  'from nest-id ', cpl_id
    ELSE
       nest_string = ''
    ENDIF
!
!-- Create the complete output string, starting with the message level
    IF ( message_level == 0 )  THEN
       header_string = '--- informative message' // TRIM(nest_string) // ' ---'
    ELSEIF ( message_level == 1 )  THEN
       header_string = '+++ warning message' // TRIM(nest_string) // ' ---'
    ELSEIF ( message_level == 2 )  THEN
       header_string = '+++ error message' // TRIM(nest_string) // ' ---'
    ELSE
       WRITE( header_string,'(A,I2)' )  '+++ unknown message level' //                             &
                                        TRIM(nest_string) // ': ', message_level
    ENDIF

!
!-- Add the message identifier and the generating routine
    header_string_2 = 'ID: ' // message_identifier //                                              &
                      '  generated by routine: ' // TRIM( routine_name )

    information_string_1 = 'Further information can be found at'
    IF ( message_identifier(1:2) == 'NC' )  THEN
       information_string_2 = 'http://palm.muk.uni-hannover.de/trac/wiki/doc/app/errmsg#NC'
    ELSE
       information_string_2 = 'http://palm.muk.uni-hannover.de/trac/wiki/doc/app/errmsg#' //       &
                              message_identifier
    ENDIF


!
!-- Output the output string and the corresponding message string which had been already assigned in
!-- the calling subroutine.
!
!-- First find out if output shall be done on this PE.
    IF ( output_on_pe == -1 )  THEN
       do_output = .TRUE.
    ELSEIF ( myid == output_on_pe )  THEN
       do_output = .TRUE.
    ENDIF
#if defined( __parallel )
!
!-- In case of illegal pe number output on pe0
    IF ( output_on_pe > numprocs-1 )  THEN
       pe_out_of_range = .TRUE.
       IF ( myid == 0 )  do_output = .TRUE.
    ENDIF
#endif

!
!-- Now do the output
    IF ( do_output )  THEN

       IF ( file_id == 6 )  THEN
!
!--       Output on stdout
          WRITE( *, '(16X,A)' )  TRIM( header_string )
          WRITE( *, '(20X,A)' )  TRIM( header_string_2 )
!
!--       Cut message string into pieces and output one piece per line.
!--       Remove leading blanks.
          message_string = ADJUSTL( message_string )
          i = INDEX( message_string, '&' )
          DO WHILE ( i /= 0 )
             WRITE( *, '(20X,A)' )  ADJUSTL( message_string(1:i-1) )
             message_string = ADJUSTL( message_string(i+1:) )
             i = INDEX( message_string, '&' )
          ENDDO
          WRITE( *, '(20X,A)' )  ''
          WRITE( *, '(20X,A)' )  TRIM( message_string )
          WRITE( *, '(20X,A)' )  ''
          WRITE( *, '(20X,A)' )  TRIM( information_string_1 )
          WRITE( *, '(20X,A)' )  TRIM( information_string_2 )
          WRITE( *, '(20X,A)' )  ''

       ELSE
!
!--       Output on requested file id (file must have been opened elsewhere!)
          WRITE( file_id, '(A/)' )  TRIM( header_string )
!
!--       Cut message string into pieces and output one piece per line.
!--       Remove leading blanks.
          message_string = ADJUSTL( message_string )
          i = INDEX( message_string, '&' )
          DO WHILE ( i /= 0 )
             WRITE( file_id, '(4X,A)' )  ADJUSTL( message_string(1:i-1) )
             message_string = ADJUSTL( message_string(i+1:) )
             i = INDEX( message_string, '&' )
          ENDDO
          WRITE( file_id, '(4X,A)' )  TRIM( message_string )
          WRITE( file_id, '(4X,A)' )  ''
          WRITE( file_id, '(4X,A)' )  TRIM( information_string_1 )
          WRITE( file_id, '(4X,A)' )  TRIM( information_string_2 )
          WRITE( file_id, '(4X,A)' )  ''
!
!--       Flush buffer, if requested
          IF ( flush_file == 1 )  FLUSH( file_id )
       ENDIF

       IF ( pe_out_of_range )  THEN
          WRITE ( *, '(A)' )  '+++ WARNING from routine message:'
          WRITE ( *, '(A,I6,A)' )  '    PE ', output_on_pe, ' choosed for output is larger '
          WRITE ( *, '(A,I6)' )  '    than the maximum number of used PEs', numprocs-1
          WRITE ( *, '(A)' )  '    Output is done on PE0 instead'
       ENDIF

    ENDIF

!
!-- Abort execution, if requested
    IF ( requested_action > 0 )  THEN
       abort_mode = requested_action
!
!--    Since nested runs always use MPI_ABORT, let only the PE which output a message initiate the
!--    abort. Others just wait.
       IF ( nested_run  .AND.  requested_action == 1  .AND.  .NOT. do_output )  THEN
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ELSE
          CALL local_stop
       ENDIF
    ENDIF

 END SUBROUTINE message


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints out the given location on stdout
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE location_message( location, message_type )

    USE, INTRINSIC ::  ISO_FORTRAN_ENV,                                                            &
        ONLY:  OUTPUT_UNIT

    USE pegrid,                                                                                    &
        ONLY:  myid

    USE pmc_interface,                                                                             &
        ONLY:  cpl_id

    IMPLICIT NONE

    CHARACTER(LEN=*)  ::  location             !< text to be output on stdout
    CHARACTER(LEN=60) ::  location_trimmed     !< trimmed text to be output on stdout
    CHARACTER(LEN=*)  ::  message_type         !< type of message; supported values: 'start', 'finished'
    CHARACTER(LEN=10) ::  message_type_string  !< formatted message-type string for output
    CHARACTER(LEN=8)  ::  system_time          !< formatted system clock time
    CHARACTER(LEN=10) ::  time                 !< current time of system

!
!-- Output for nested runs only on the root domain
    IF ( cpl_id /= 1 )  RETURN

    IF ( myid == 0 )  THEN
!
!--    Get system time for debug info output (helpful to estimate the required computing time for
!--    specific parts of code)
       CALL date_and_time( TIME=time )
       system_time = time(1:2) // ':' // time(3:4) // ':' // time(5:6)
!
!--    Write message-type string depending on message_type
       message_type_string = REPEAT( '-', 10 )
       IF ( TRIM( message_type ) == 'start' )                                                      &
          message_type_string(2:) = TRIM( message_type ) // '----'
       IF ( TRIM( message_type ) == 'finished' )                                                   &
          message_type_string(2:) = TRIM( message_type ) // '-'
!
!--    Trim location text to a maximum of 60 chars
!--    Note: if the length is set within the write format, the string is right-aligned; to trim and
!--    left-align the output, we need to use this detour
       WRITE( location_trimmed, '(A)' )  ADJUSTL( TRIM( location ) )
!
!--    Write and flush debug location or info message to file
       WRITE( OUTPUT_UNIT, 200 )  system_time, message_type_string, TRIM( location_trimmed )
       FLUSH( OUTPUT_UNIT )
!
!--    Message formats
200    FORMAT ( 3X, A, 3x, A, 3X, A )

    ENDIF

 END SUBROUTINE location_message


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Prints out the given debug information to unit 9 (DEBUG files in temporary directory)
!> for each PE on each domain.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE debug_message( debug_string, message_type )

    USE control_parameters,                                                                        &
        ONLY:  time_since_reference_point

    IMPLICIT NONE

    CHARACTER(LEN=*)  ::  debug_string         !< debug message to be output to debug_output_unit
    CHARACTER(LEN=*)  ::  message_type         !< type of message; supported values: 'start', 'end', 'info'
    CHARACTER(LEN=7)  ::  message_type_string  !< formatted message-type string for output
    CHARACTER(LEN=8)  ::  system_time          !< formatted system clock time
    CHARACTER(LEN=10) ::  time                 !< current time of system

    INTEGER, PARAMETER ::  debug_output_unit = 9

!
!-- Get system time for debug info output (helpful to estimate the required computing time for
!-- specific parts of code)
    CALL date_and_time( TIME=time )
    system_time = time(1:2) // ':' // time(3:4) // ':' // time(5:6)
!
!-- Write message-type string depending on message_type
    message_type_string = REPEAT( '-', 7 )
    IF ( TRIM( message_type ) == 'start' )  message_type_string(2:) = TRIM( message_type ) // '-'
    IF ( TRIM( message_type ) == 'end' )    message_type_string(2:) = TRIM( message_type ) // '---'
    IF ( TRIM( message_type ) == 'info' )   message_type_string(2:) = TRIM( message_type ) // '--'
!
!-- Write and flush debug location or info message to file
    WRITE( debug_output_unit, 201 )    system_time, time_since_reference_point, &
                                       message_type_string, TRIM( debug_string )
    FLUSH( debug_output_unit )
!
!-- Message formats
201 FORMAT ( 'System time: ', A, ' | simulated time (s): ', F12.3, ' | ', A, ' ', A )

 END SUBROUTINE debug_message


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Abort routine for failures durin reading of namelists
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE parin_fail_message( location, line )

    USE control_parameters,                                                                        &
        ONLY:  message_string

    USE kinds

    IMPLICIT NONE

    CHARACTER(LEN=*) ::  location !< text to be output on stdout
    CHARACTER(LEN=*) ::  line

    CHARACTER(LEN=80) ::  line_dum

    INTEGER(iwp) ::  line_counter

    line_dum = ' '
    line_counter = 0

    REWIND( 11 )
    DO WHILE ( INDEX( line_dum, TRIM(line) ) == 0 )
       READ ( 11, '(A)', END=20 )  line_dum
       line_counter = line_counter + 1
    ENDDO

 20 WRITE( message_string, '(A,I3,A)' ) 'Error(s) in NAMELIST '// TRIM(location) //                &
                                        '&Reading fails on line ', line_counter, ' at&' // line
    CALL message( 'parin', 'PA0271', 1, 2, 0, 6, 0 )

 END SUBROUTINE parin_fail_message
