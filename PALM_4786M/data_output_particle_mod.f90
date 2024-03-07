!> @file data_output_particle_mod.f90
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
! $Id: data_output_particle_mod.f90 4779 2020-11-09 17:45:22Z suehring $
! Avoid overlong lines and unused variables
!
! 4778 2020-11-09 13:40:05Z raasch
! Initial implementation (K. Ketelsen)
!
!
! Description:
! ------------
!> Output of particle time series
!--------------------------------------------------------------------------------------------------!
 MODULE data_output_particle_mod

#if defined( __parallel )
   USE MPI
#endif

#if defined( __netcdf4 )
   USE NETCDF
#endif

   USE, INTRINSIC ::  ISO_C_BINDING

   USE kinds,                                                                                     &
      ONLY: wp, iwp, sp, dp, idp, isp

   USE indices,                                                                                   &
       ONLY:  nbgp, nnx, nny, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz, nzb, nzt

   USE control_parameters,                                                                        &
       ONLY:   end_time, simulated_time, dt_dopts

   USE pegrid,                                                                                    &
       ONLY:  comm1dx, comm1dy, comm2d, myid, myidx, myidy, npex, npey, numprocs, pdims

   USE particle_attributes,                                                                       &
       ONLY: particle_type, particles, grid_particle_def, grid_particles, prt_count,              &
             unlimited_dimension, pts_id_file, data_output_pts, pts_increment, pts_percentage, &
             oversize, number_of_output_particles

   USE shared_memory_io_mod,                                                                      &
       ONLY: sm_class

   USE data_output_netcdf4_module,                                                                   &
       ONLY: netcdf4_init_dimension, netcdf4_get_error_message, netcdf4_stop_file_header_definition, &
             netcdf4_init_module, netcdf4_init_variable, netcdf4_finalize,                           &
             netcdf4_open_file, netcdf4_write_attribute, netcdf4_write_variable

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

   IMPLICIT NONE

   PRIVATE
   SAVE


   LOGICAL, PUBLIC      :: dop_active = .FALSE.


!  Variables for restart

   INTEGER(iwp), PUBLIC    :: dop_prt_axis_dimension
   INTEGER(iwp), PUBLIC    :: dop_last_active_particle




!kk Private module variables can be declared inside the submodule

   INTEGER,PARAMETER    :: MAX_NR_VARIABLES    = 128
   INTEGER,PARAMETER    :: NR_FIXED_VARIABLES  = 16

   CHARACTER(LEN=32)                         :: file_name            !< Name of NetCDF file
   INTEGER(iwp)                              :: nr_time_values       !< Number of values on time axis
   REAL(sp),ALLOCATABLE,DIMENSION(:)         :: time_axis_values     !< time axis Values
   INTEGER(iwp),ALLOCATABLE,DIMENSION(:,:)   :: rma_particles        !< Start address and number of remote particles
   INTEGER(iwp)                              :: nr_particles_PE      !< Number of particles assigned for output on this thread
   INTEGER(iwp)                              :: nr_particles_out     !< total number od particles assigned for output

   INTEGER(iwp),ALLOCATABLE,DIMENSION(:,:)   :: sh_indices                   !< Indices in shared memory group
   INTEGER(iwp),ALLOCATABLE,DIMENSION(:,:)   :: io_indices                   !< Indices on IO processes
   INTEGER(iwp),ALLOCATABLE,DIMENSION(:,:)   :: mo_indices                   !< Indices for model communicator
   INTEGER(iwp),ALLOCATABLE,DIMENSION(:,:)   :: remote_nr_particles
   INTEGER(iwp)                              :: start_local_numbering        !< start increment 1 numbering on this thread
   INTEGER(iwp)                              :: end_local_numbering
   INTEGER(iwp)                              :: pe_start_index               !< start index of the output area on this thread
   INTEGER(iwp)                              :: pe_end_index
   INTEGER(iwp)                              :: io_start_index               !< start index of the output area on IO thread
   INTEGER(iwp)                              :: io_end_index
   INTEGER(iwp)                              :: nr_particles_rest            !< Numbere of rest Particle (ireg. distribution)
   LOGICAL                                   :: irregular_distribubtion      !< irregular distribution of output particlesexit

   TYPE var_def
      INTEGER(iwp)         :: var_id
      CHARACTER(len=32)    :: name
      CHARACTER(len=32)    :: units
      LOGICAL              :: is_integer = .FALSE.
   END TYPE var_def

   TYPE(var_def),DIMENSION(MAX_NR_VARIABLES)    :: variables
   TYPE(var_def),DIMENSION(NR_FIXED_VARIABLES)  :: fix_variables

   TYPE(sm_class)  :: part_io                 !< manage communicator for particle IO

!
!  NetCDF

   INTEGER(iwp)                              :: file_id = -1         !< id of Netcdf file
   INTEGER(iwp)                              :: nr_fix_variables     !< Number of fixed variables scheduled for output
   INTEGER(iwp)                              :: nr_variables         !< Number of variables  scheduled for output

   TYPE dimension_id
      INTEGER(iwp)         :: prt
      INTEGER(iwp)         :: time
   END TYPE dimension_id

   TYPE variable_id
      INTEGER(iwp)         :: prt
      INTEGER(iwp)         :: time
   END TYPE variable_id

   TYPE(dimension_id)                        :: did
   TYPE(variable_id)                         :: var_id

!  shared memory buffer
   INTEGER(iwp)    :: win_prt_i = -1                                 ! integer MPI shared Memory window
   INTEGER(iwp)    :: win_prt_r = -1                                 !< real  MPI shared Memory window
   INTEGER(iwp), POINTER, CONTIGUOUS, DIMENSION(:) :: out_buf_i      !< integer output buffer
   REAL(sp), POINTER, CONTIGUOUS, DIMENSION(:)     :: out_buf_r      !< real output buffer

!
!  Particle list in file

   LOGICAL         :: part_list_in_file
   INTEGER(idp),ALLOCATABLE, DIMENSION(:)   :: part_id_list_file
!
!  RMA window
#if defined( __parallel )
   INTEGER(iwp)    ::  win_rma_buf_i
   INTEGER(iwp)    ::  win_rma_buf_r
#endif

   INTEGER(iwp),POINTER,DIMENSION(:)         :: transfer_buffer_i    !< rma window to provide data, which can be fetch via MPI_Get
   REAL(sp),POINTER,DIMENSION(:)             :: transfer_buffer_r    !< Same for REAL
   INTEGER(iwp),ALLOCATABLE,DIMENSION(:)     :: remote_indices       !< particle nubmer of the remodte partices, used as indices in the output array
   INTEGER(iwp)                              :: initial_number_of_active_particles

!  Public subroutine Interface

   INTERFACE dop_init
      MODULE PROCEDURE dop_init
   END INTERFACE dop_init

   INTERFACE dop_output_tseries
      MODULE PROCEDURE dop_output_tseries
   END INTERFACE dop_output_tseries

   INTERFACE dop_finalize
      MODULE PROCEDURE dop_finalize
   END INTERFACE  dop_finalize

#if defined( __parallel )
   INTERFACE dop_alloc_rma_mem
      MODULE PROCEDURE dop_alloc_rma_mem_i1
      MODULE PROCEDURE dop_alloc_rma_mem_r1
   END INTERFACE dop_alloc_rma_mem
#endif

   PUBLIC dop_init, dop_output_tseries, dop_finalize
#if defined( __parallel )
   PUBLIC dop_alloc_rma_mem       ! Must be PUBLIC on NEC, although it is only used in Submodule
#endif



 CONTAINS

   SUBROUTINE dop_init (read_restart)
      IMPLICIT NONE
      LOGICAL,INTENT(IN)    :: read_restart

      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: nr_particles_local          !< total number of particles scheduled for output on this thread
#if defined( __parallel )
      INTEGER(iwp) :: ierr                        !< MPI error code
#endif
      INTEGER(idp) :: nr_particles_8              !< Total number of particles in 64 bit
      REAL(dp)     :: xnr_part                    !< Must be 64 Bit REAL
      INTEGER(iwp) :: nr_local_last_pe               !< Number of output particles on myid == numprocs-2

      INTEGER(iwp),DIMENSION(0:numprocs-1)     :: nr_particles_all_s
      INTEGER(iwp),DIMENSION(0:numprocs-1)     :: nr_particles_all_r
      INTEGER(idp),DIMENSION(0:numprocs-1)     :: nr_particles_all_8

      INTEGER(iwp),ALLOCATABLE,DIMENSION(:,:)  :: sh_indices_s
      INTEGER(iwp),ALLOCATABLE,DIMENSION(:,:)  :: io_indices_s
      INTEGER(iwp),ALLOCATABLE,DIMENSION(:,:)  :: mo_indices_s

      part_list_in_file = (pts_id_file /= ' ')

      IF(.NOT. part_list_in_file)   THEN
         IF(.NOT. read_restart)   THEN

            CALL set_indef_particle_nr

            CALL count_output_particles (nr_particles_local)

            nr_particles_all_s = 0
            nr_particles_all_s(myid) = nr_particles_local

#if defined( __parallel )
            CALL MPI_ALLREDUCE( nr_particles_all_s, nr_particles_all_r, SIZE( nr_particles_all_s ), MPI_INTEGER,    &
               MPI_SUM, comm2d, ierr )
#else
            nr_particles_all_r = nr_particles_all_s
#endif

            initial_number_of_active_particles = SUM(nr_particles_all_r)

            start_local_numbering = 1
            end_local_numbering   = nr_particles_all_r(0)

            IF ( myid > 0)   THEN
               DO i=1,numprocs-1
                  start_local_numbering = start_local_numbering+nr_particles_all_r(i-1)
                  end_local_numbering   = start_local_numbering+nr_particles_all_r(i)-1
                  IF(myid == i)  EXIT
               ENDDO
            END IF

            nr_particles_all_8 = nr_particles_all_r              ! Use 64 INTEGER, maybe more than 2 GB particles
            nr_particles_8 = SUM(nr_particles_all_8)
            IF (nr_particles_8 > 2000000000_8)   THEN            ! This module works only with 32 Bit INTEGER
               write(9,*) 'Number of particles too large ',nr_particles_8; flush(9)
#if defined( __parallel )
               CALL MPI_ABORT (MPI_COMM_WORLD, 1, ierr)
#endif
            ENDIF
            nr_particles_out = nr_particles_8                    ! Total number of particles scheduled for output

!
!--         reserve space for additional particle
            IF(number_of_output_particles > 0)   THEN
               nr_particles_out = MAX(number_of_output_particles, nr_particles_out)
            ELSE IF(oversize > 100.)   THEN
               xnr_part = nr_particles_out*oversize/100.
               nr_particles_out = xnr_part
            ENDIF
         ELSE
            nr_particles_out = dop_prt_axis_dimension            ! Get Number from restart file
            initial_number_of_active_particles = dop_last_active_particle
         ENDIF

      ELSE
         CALL set_indef_particle_nr
         CALL dop_read_output_particle_list (nr_particles_out)
      ENDIF
      dop_prt_axis_dimension = nr_particles_out
!
!--   The number of particles must be at least the number of MPI processes

      nr_particles_out = MAX(nr_particles_out, numprocs)


      nr_particles_PE = (nr_particles_out+numprocs-1)/numprocs   ! Number of paricles scheduled for ouput on this thread

      pe_start_index = myid*nr_particles_PE+1                              !Output numberimng on this thread
      pe_end_index   = MIN((myid+1)*nr_particles_PE,nr_particles_out)

      irregular_distribubtion = .FALSE.
#if defined( __parallel )
!
!--   In case of few particles, it can happen that not only the last thread gets fewer output particles.
!--   In this case, the local number of particles on thread numprocs-1 will be < 1
!--   If this happens, irregular distribution of output particles will be used

      IF(myid == numprocs-1)   Then
         nr_local_last_pe = pe_end_index-pe_start_index+1
      ELSE
         nr_local_last_pe = 0
      ENDIF
      CALL MPI_BCAST (nr_local_last_pe, 1, MPI_INTEGER, numprocs-1, comm2d, ierr)
#else
      nr_local_last_pe = nr_particles_PE
#endif
      IF(nr_local_last_pe < 1)   THEN
         irregular_distribubtion = .TRUE.
         CALL dop_setup_ireg_distribution
      ENDIF


      IF(.NOT. read_restart .AND. .NOT. part_list_in_file)   THEN
         CALL set_particle_number
      ENDIF

      IF(part_list_in_file)   THEN
         CALL dop_find_particle_in_outlist
      ENDIF

      CALL part_io%sm_init_data_output_particles ()

      ALLOCATE (sh_indices_s(2,0:part_io%sh_npes-1))
      ALLOCATE (sh_indices(2,0:part_io%sh_npes-1))


#if defined( __parallel )
      sh_indices_s = 0
      sh_indices_s(1,part_io%sh_rank) = pe_start_index
      sh_indices_s(2,part_io%sh_rank) = pe_end_index

      CALL MPI_ALLREDUCE( sh_indices_s, sh_indices, 2*part_io%sh_npes, MPI_INTEGER,    &
                                              MPI_SUM, part_io%comm_shared, ierr )

      io_start_index = sh_indices(1,0)                         ! output numbering on actual IO thread
      io_end_index   = sh_indices(2,part_io%sh_npes-1)
#else
      io_start_index = pe_start_index                          ! output numbering
      io_end_index   = pe_end_index
#endif


#if defined( __parallel )
      CALL MPI_BCAST (part_io%io_npes, 1, MPI_INTEGER, 0,  part_io%comm_shared, ierr)
#endif
      ALLOCATE (io_indices(2,0:part_io%io_npes-1))
      IF (part_io%iam_io_pe)   THEN
         ALLOCATE (io_indices_s(2,0:part_io%io_npes-1))


         io_indices_s = 0
         io_indices_s(1,part_io%io_rank) = io_start_index
         io_indices_s(2,part_io%io_rank) = io_end_index

#if defined( __parallel )
         CALL MPI_ALLREDUCE( io_indices_s, io_indices, 2*part_io%io_npes, MPI_INTEGER,    &
            MPI_SUM, part_io%comm_io, ierr )
#else
         io_indices = io_indices_s
#endif

      ENDIF

#if defined( __parallel )
      CALL MPI_BCAST (io_indices, size(io_indices), MPI_INTEGER, 0,  part_io%comm_shared, ierr)
#endif

      ALLOCATE (remote_nr_particles(2,0:numprocs-1))
      ALLOCATE (rma_particles(2,0:numprocs-1))

      ALLOCATE (mo_indices(2,0:numprocs-1))
      ALLOCATE (mo_indices_s(2,0:numprocs-1))


      mo_indices_s = 0
      mo_indices_s(1,myid) = pe_start_index
      mo_indices_s(2,myid) = pe_end_index

#if defined( __parallel )
      CALL MPI_ALLREDUCE( mo_indices_s, mo_indices, 2*numprocs, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
      mo_indices = mo_indices_s
#endif

!--   Allocate output buffer

#if defined( __parallel )
      CALL part_io%sm_allocate_shared (out_buf_r, io_start_index, io_end_index, win_prt_r)
      CALL part_io%sm_allocate_shared (out_buf_i, io_start_index, io_end_index, win_prt_i)
#else
      ALLOCATE(out_buf_r(io_start_index:io_end_index))
      ALLOCATE(out_buf_i(io_start_index:io_end_index))
#endif

!--   NetCDF

      CALL dop_netcdf_setup ()

#if defined( __parallel )
      CALL MPI_BCAST (nr_fix_variables, 1, MPI_INTEGER, 0,  part_io%comm_shared, ierr)
      CALL MPI_BCAST (nr_variables, 1, MPI_INTEGER, 0,  part_io%comm_shared, ierr)
#endif

      CALL dop_count_remote_particles

      CALL dop_write_fixed_variables

      CALL deallocate_and_free

      CALL dop_output_tseries

      RETURN
   CONTAINS
      SUBROUTINE dop_setup_ireg_distribution
         IMPLICIT NONE


         nr_particles_PE   = (nr_particles_out)/numprocs   ! Number of paricles scheduled for ouput on this thread

         nr_particles_rest = nr_particles_out - numprocs * nr_particles_PE

         nr_particles_PE = nr_particles_PE+1
         IF(myid < nr_particles_rest)   THEN
            pe_start_index = myid*nr_particles_PE+1                              !Output numberimng on this thread
            pe_end_index   = MIN((myid+1)*nr_particles_PE,nr_particles_out)
         ELSE
            pe_start_index = nr_particles_rest*(nr_particles_PE)+(myid-nr_particles_rest)*(nr_particles_PE-1)+1                              !Output numberimng on this thread
            pe_end_index   = MIN(pe_start_index+nr_particles_PE-2,nr_particles_out)
         ENDIF


      END SUBROUTINE dop_setup_ireg_distribution

   END SUBROUTINE dop_init
!
!  particle output on selected time steps

   SUBROUTINE dop_output_tseries
      IMPLICIT NONE

      INTEGER(iwp)             :: i                           !<
      INTEGER(iwp)             :: return_value                !< Return value data_output_netcdf4 .. routines
#if defined( __parallel )
      INTEGER(iwp)             :: ierr                        !< MPI error code
#endif
      INTEGER(iwp),SAVE        :: icount=0                    !< count output steps

      INTEGER(iwp),DIMENSION(2)      :: bounds_origin, bounds_start, value_counts
      REAl(wp),POINTER, CONTIGUOUS, DIMENSION(:)  :: my_time

      icount = icount+1

      CALL dop_delete_particle_number

      IF(part_list_in_file)   THEN

         CALL dop_find_particle_in_outlist
      ELSE
         CALL dop_newly_generated_particles
      ENDIF

      CALL dop_count_remote_particles

      bounds_origin    = 1

      bounds_start(1)  = io_start_index
      bounds_start(2)  = icount

      value_counts(1)  = io_end_index-io_start_index+1
      value_counts(2)  = 1

      DO i=1,nr_variables
#if defined( __netcdf4 )
         IF(variables(i)%is_integer)   THEN
            out_buf_i = NF90_FILL_INT
         ELSE
            out_buf_r = NF90_FILL_REAL
         ENDIF
#endif

         CALL cpu_log( log_point_s(99), 'dop_fill_out_buf', 'start' )
         CALL dop_fill_out_buf (variables(i))
         CALL cpu_log( log_point_s(99), 'dop_fill_out_buf', 'stop' )

         CALL cpu_log( log_point_s(88), 'dop_get_remote_particle', 'start' )
#if defined( __parallel )
         CALL  dop_get_remote_particle (variables(i)%is_integer)
#endif
         CALL cpu_log( log_point_s(88), 'dop_get_remote_particle', 'stop' )

         CALL cpu_log( log_point_s(89), 'particle NetCDF output', 'start' )
         CALL  part_io%sm_node_barrier()
         IF (part_io%iam_io_pe)   THEN
            IF(variables(i)%is_integer)   THEN
               CALL  netcdf4_write_variable('parallel', file_id, variables(i)%var_id, bounds_start,&
                                             value_counts, bounds_origin,                          &
                                             .FALSE., values_int32_1d=out_buf_i, return_value=return_value)
            ELSE
               CALL  netcdf4_write_variable('parallel', file_id, variables(i)%var_id, bounds_start,&
                                             value_counts, bounds_origin,                          &
                                            .FALSE., values_real32_1d=out_buf_r, return_value=return_value)
            ENDIF
         ENDIF
         CALL  part_io%sm_node_barrier()
         CALL cpu_log( log_point_s(89), 'particle NetCDF output', 'stop' )

#if defined( __parallel )
         CALL MPI_BARRIER(comm2d, ierr)               !kk This Barrier is necessary, not sure why
#endif
      END DO

      CALL deallocate_and_free

!     write Time value
      IF(myid == 0)   THEN
         ALLOCATE(my_time(1))
         bounds_start(1)  = icount
         value_counts(1)  = 1
         my_time(1)       = simulated_time
         IF(unlimited_dimension)   THEN               ! For workaround described in dop finalize
            time_axis_values(icount) = my_time(1)
         ELSE
            CALL  netcdf4_write_variable('parallel', file_id, var_id%time, bounds_start(1:1),      &
                                          value_counts(1:1), bounds_origin(1:1),                   &
                                         .TRUE., values_realwp_1d=my_time, return_value=return_value)
         ENDIF
         DEALLOCATE(my_time)
      ENDIF

      RETURN
   END SUBROUTINE dop_output_tseries

   SUBROUTINE dop_finalize
      IMPLICIT NONE
#if defined( __netcdf4 )
      INTEGER(iwp) :: var_len
#endif
      INTEGER(iwp) :: return_value                !< Return value data_output_netcdf4 .. routines
#if defined( __netcdf4 )
      INTEGER(iwp) :: ierr                        !< MPI error code
#endif

      IF( win_prt_i /= -1)   THEN
         CALL part_io%sm_free_shared (win_prt_i)
      ENDIF
      IF( win_prt_r /= -1)   THEN
         CALL part_io%sm_free_shared (win_prt_r)
      ENDIF

      IF( file_id /= -1 .AND. part_io%iam_io_pe)   THEN
         CALL netcdf4_finalize( 'parallel', file_id, return_value )
         file_id = -1

#if defined( __netcdf4 )
!
!--      For yet unknown reasons it is not possible to write in parallel mode the time values to NetCDF variable time
!--      This workaround closes the parallel file and writes the time values in sequential mode on PE0
!kk      This is a real quick and dirty workaround and the problem should be solved in the final version !!!

         If(myid == 0 .AND. unlimited_dimension)   THEN
            ierr = nf90_open (TRIM(file_name),NF90_WRITE, file_id)
            ierr = nf90_inquire_dimension(file_id, did%time, len = var_len)

            ierr = nf90_put_var (file_id, var_id%time, time_axis_values(1:var_len))

            ierr = nf90_close (file_id)
         ENDIF
#endif
      ENDIF

      RETURN
   END SUBROUTINE dop_finalize

!  Private subroutines


!
!  kk    Not sure if necessary, but set ALL particle numnber to -1

   SUBROUTINE set_indef_particle_nr

      IMPLICIT NONE
      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: k                           !<
      INTEGER(iwp) :: n                           !<

      DO i=nxl,nxr
         DO j=nys,nyn
            DO k=nzb+1,nzt
               DO n=1,SIZE(grid_particles(k,j,i)%particles)
                  grid_particles(k,j,i)%particles(n)%particle_nr = -1
               END DO
            END DO
         ENDDO
      ENDDO

      RETURN
   END SUBROUTINE set_indef_particle_nr
!
!-- Count particles scheduled for output
!-- here are pts_increment and pts_percentage are used to select output particles

   SUBROUTINE count_output_particles (pcount)
      IMPLICIT NONE

      INTEGER(iwp), INTENT(OUT)  :: pcount        !<

      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: k                           !<
      INTEGER(iwp) :: n                           !<
      INTEGER(iwp) :: n_all                       !< count all particles for MOD function
      REAL(dp)     :: fcount
      REAL(dp)     :: finc

      pcount = 0
      IF(pts_increment == 1)   THEN
         pcount = SUM(prt_count)
      ELSE
         n_all = 0
         DO i=nxl,nxr
            DO j=nys,nyn
               DO k=nzb+1,nzt
                  DO n=1,prt_count(k,j,i)
                     IF(MOD(n_all,pts_increment) == 0)   THEN
                        pcount = pcount+1
                     ENDIF
                     n_all = n_all+1
                  END DO
               END DO
            ENDDO
         ENDDO
      ENDIF

      IF(pts_percentage < 100. )   THEN

         finc   = pts_percentage/100

         pcount = 0
         fcount = 0.0

         DO i=nxl,nxr
            DO j=nys,nyn
               DO k=nzb+1,nzt
                  DO n=1,prt_count(k,j,i)
                     fcount = fcount + finc
                     IF(pcount < int(fcount) )  THEN
                        pcount = pcount+1
                     ENDIF
                  END DO
               END DO
            ENDDO
         ENDDO

      ENDIF

      RETURN
   END SUBROUTINE count_output_particles

   SUBROUTINE dop_read_output_particle_list (nr_particles_out)
      IMPLICIT NONE
      INTEGER(iwp),INTENT(OUT)           :: nr_particles_out

      INTEGER(iwp)             :: i      !<
      INTEGER(iwp)             :: iu     !<
      INTEGER(iwp)             :: istat  !<
      INTEGER(idp)             :: dummy  !<

      iu    = 345
      istat = 0
      nr_particles_out = 0

      OPEN(unit=iu, file=TRIM(pts_id_file))      !kk should be changed to check_open

!     First stridem cout output particle

      DO WHILE (istat == 0)
         READ(iu,*,iostat=istat)  dummy
         nr_particles_out = nr_particles_out+1
      END DO

      nr_particles_out = nr_particles_out-1             ! subtract 1 for end of file read

      ALLOCATE(part_id_list_file(nr_particles_out))

      REWIND(iu)

!--   second stride, read particle ids for scheduled output particle

      DO i=1,nr_particles_out
         READ(iu,*) part_id_list_file(i)
      END DO


      CLOSE (iu)

   END SUBROUTINE dop_read_output_particle_list
!
!-- Setb output particle number for selected active particles

   SUBROUTINE set_particle_number
      IMPLICIT NONE

      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: k                           !<
      INTEGER(iwp) :: n                           !<
      INTEGER(iwp) :: n_all                       !< count all particles for MOD function
      INTEGER(iwp) :: particle_nr                 !< output particle number
      INTEGER(iwp) :: pcount                      !< local particle count in case of pts_percentage
      REAL(dp)     :: fcount                      !< partical progress in %/100
      REAL(dp)     :: finc                        !< increment of particle

      pcount = 0
      fcount = 0.0

      particle_nr = start_local_numbering
      n_all       = 0
      DO i=nxl,nxr
         DO j=nys,nyn
            DO k=nzb+1,nzt
               IF(pts_increment > 1)   THEN
                  DO n=1,prt_count(k,j,i)
                     IF(MOD(n_all,pts_increment) == 0)   THEN
                        grid_particles(k,j,i)%particles(n)%particle_nr = particle_nr
                        particle_nr = particle_nr+1
                     ELSE
                        grid_particles(k,j,i)%particles(n)%particle_nr = -2
                     ENDIF
                     n_all = n_all+1
                  END DO
               ELSE IF(pts_percentage < 100. )   THEN
                  finc   = pts_percentage/100

                  DO n=1,prt_count(k,j,i)
!
!--                  Every particle move fraction on particle axis (i.e part percent == 80; move 0.8)
!--                  if increases next whole number, the particle is taken for output
                     fcount = fcount + finc
                     IF(pcount < int(fcount) )  THEN
                        pcount = pcount+1
                        grid_particles(k,j,i)%particles(n)%particle_nr = particle_nr
                        particle_nr = particle_nr+1
                     ELSE
                        grid_particles(k,j,i)%particles(n)%particle_nr = -2
                     ENDIF
                  END DO
               ELSE
                  DO n=1,prt_count(k,j,i)
                     grid_particles(k,j,i)%particles(n)%particle_nr = particle_nr
                     particle_nr = particle_nr+1
                  END DO
               ENDIF
            END DO
         ENDDO
      ENDDO

      RETURN
   END SUBROUTINE set_particle_number

   SUBROUTINE dop_find_particle_in_outlist
      IMPLICIT NONE

      INTEGER(iwp) :: i                           !<
#if defined( __parallel )
      INTEGER(iwp) :: ierr                        !< MPI error code
#endif
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: k                           !<
      INTEGER(iwp) :: l                           !<
      INTEGER(iwp) :: n                           !<
      INTEGER(iwp) :: nr_part                     !<
!      INTEGER, save :: icount=0

      nr_part = 0

!
!--   If there is a long particle output list, for performance reason it may become necessary to optimize the
!--   following loop.
!
!--       Decode the particle id in i,j,k,n
!--       Split serach, i.e search first in i, then in j ...

      DO i=nxl,nxr
         DO j=nys,nyn
            DO k=nzb+1,nzt
               DO n=1,prt_count(k,j,i)
                  DO l=1,SIZE(part_id_list_file)
                     IF(grid_particles(k,j,i)%particles(n)%id == part_id_list_file(l))   THEN
                        grid_particles(k,j,i)%particles(n)%particle_nr =  l
                        nr_part = nr_part+1
                     ENDIF
                  END DO
               END DO
            END DO
         ENDDO
      ENDDO

#if defined( __parallel )
      CALL MPI_ALLREDUCE( nr_part, initial_number_of_active_particles, 1, MPI_INTEGER, MPI_SUM, comm2d, ierr )
#else
      initial_number_of_active_particles = nr_part
#endif

   END SUBROUTINE dop_find_particle_in_outlist

!-- Netcdf Setup
!
!--   Open NetCDF File DATA_1D_PTS_NETCDF
!--   Define Dimensions and variables
!--   Write constant variables

   SUBROUTINE dop_netcdf_setup
      IMPLICIT NONE

      INTEGER,PARAMETER        :: global_id_in_file = -1
      INTEGER(iwp)             :: i
      INTEGER(iwp)             :: fix_ind
      INTEGER(iwp)             :: var_ind

      INTEGER(iwp)             :: return_value
      LOGICAL                  :: const_flag

      INTEGER, DIMENSION(2)    :: dimension_ids

      nr_time_values = end_time/dt_dopts+1               !kk has to be adapted to formular of 3d output

      IF (part_io%iam_io_pe)   THEN
         CALL netcdf4_init_module( "", part_io%comm_io, 0, 9, .TRUE., -1 )

         file_name = 'DATA_1D_PTS_NETCDF'
#if defined( __parallel )
         CALL netcdf4_open_file( 'parallel', trim(file_name), file_id, return_value )
#else
         CALL netcdf4_open_file( 'serial', trim(file_name), file_id, return_value )
#endif
!
!--      global attributes

         CALL netcdf4_write_attribute( 'parallel', file_id, global_id_in_file, 'comment', &
                     'Particle ouput created by PALM module data_output_particle', return_value=return_value )

         CALL netcdf4_write_attribute( 'parallel', file_id, global_id_in_file, 'initial_nr_particles', &
                     value_int32=initial_number_of_active_particles, return_value=return_value )
!
!--      define dimensions
         CALL netcdf4_init_dimension( 'parallel', file_id, did%prt, var_id%prt, &
               'prt','int32' , nr_particles_out, .TRUE., return_value )

         IF(unlimited_dimension)   THEN
            CALL netcdf4_init_dimension( 'parallel', file_id, did%time, var_id%time, &
                                 'time','real32' , -1, .TRUE., return_value )
            ALLOCATE (time_axis_values(nr_time_values))
         ELSE
            CALL netcdf4_init_dimension( 'parallel', file_id, did%time, var_id%time, &
                                 'time','real32' , nr_time_values, .TRUE., return_value )
         ENDIF
      END IF
!
!--   Variables without time axis
!--   These variables will always be written only once at the beginning of the file

      dimension_ids(1) = did%prt

      fix_ind = 1
      fix_variables(fix_ind)%name  = 'origin_x'
      fix_variables(fix_ind)%units = 'meter'

      IF (part_io%iam_io_pe)   THEN
         CALL netcdf4_init_variable( 'parallel', file_id, fix_variables(fix_ind)%var_id, fix_variables(fix_ind)%name, 'real32',  &
                                   dimension_ids(1:1), .FALSE., return_value )

         CALL netcdf4_write_attribute( 'parallel', file_id, fix_variables(fix_ind)%var_id, 'units', &
                                   value_char=trim(fix_variables(fix_ind)%units), return_value=return_value )
      ENDIF

      fix_ind = fix_ind+1
      fix_variables(fix_ind)%name  = 'origin_y'
      fix_variables(fix_ind)%units = 'meter'

      IF (part_io%iam_io_pe)   THEN
         CALL netcdf4_init_variable( 'parallel', file_id, fix_variables(fix_ind)%var_id, fix_variables(fix_ind)%name, 'real32',  &
                                   dimension_ids(1:1), .FALSE., return_value )

         CALL netcdf4_write_attribute( 'parallel', file_id, fix_variables(fix_ind)%var_id, 'units', &
                                   value_char=trim(fix_variables(fix_ind)%units), return_value=return_value )

      ENDIF

      fix_ind = fix_ind+1
      fix_variables(fix_ind)%name  = 'origin_z'
      fix_variables(fix_ind)%units = 'meter'

      IF (part_io%iam_io_pe)   THEN
         CALL netcdf4_init_variable( 'parallel', file_id, fix_variables(fix_ind)%var_id, fix_variables(fix_ind)%name, 'real32',  &
                                   dimension_ids(1:1), .FALSE., return_value )

         CALL netcdf4_write_attribute( 'parallel', file_id, fix_variables(fix_ind)%var_id, 'units', &
                                   value_char=trim(fix_variables(fix_ind)%units), return_value=return_value )
      ENDIF

!
!--   These variables are written if name end with _const'
      DO i=1,size(data_output_pts)
         const_flag = (INDEX(TRIM(data_output_pts(i)),'_const') > 0)
         IF(LEN(TRIM(data_output_pts(i))) > 0 .AND. const_flag)    THEN
            fix_ind = fix_ind+1
            fix_variables(fix_ind)%name  = TRIM(data_output_pts(i))

            SELECT CASE (TRIM(fix_variables(fix_ind)%name))
               CASE('radius_const')
                  fix_variables(fix_ind)%units = 'meter'
                  fix_variables(fix_ind)%is_integer = .FALSE.
               CASE('aux1_const')
                  fix_variables(fix_ind)%units = 'depend_on_setup'
                  fix_variables(fix_ind)%is_integer = .FALSE.
               CASE('aux2_const')
                  fix_variables(fix_ind)%units = 'depend_on_setup'
                  fix_variables(fix_ind)%is_integer = .FALSE.
               CASE('rvar1_const')
                  fix_variables(fix_ind)%units = 'depend_on_setup'
                  fix_variables(fix_ind)%is_integer = .FALSE.
               CASE('rvar2_const')
                  fix_variables(fix_ind)%units = 'depend_on_setup'
                  fix_variables(fix_ind)%is_integer = .FALSE.
               CASE('rvar3_const')
                  fix_variables(fix_ind)%units = 'depend_on_setup'
                  fix_variables(fix_ind)%is_integer = .FALSE.
            END SELECT

            IF (part_io%iam_io_pe)   THEN
               IF(fix_variables(fix_ind)%is_integer)   THEN
                  CALL netcdf4_init_variable( 'parallel', file_id, fix_variables(fix_ind)%var_id,  &
                                              fix_variables(fix_ind)%name, 'int32',                &
                                              dimension_ids(1:1), .FALSE., return_value )
               ELSE
                  CALL netcdf4_init_variable( 'parallel', file_id, fix_variables(fix_ind)%var_id,  &
                                               fix_variables(fix_ind)%name, 'real32',              &
                                               dimension_ids(1:1), .FALSE., return_value )
               ENDIF

               CALL netcdf4_write_attribute( 'parallel', file_id, fix_variables(fix_ind)%var_id, 'units', &
                  value_char=trim(fix_variables(fix_ind)%units), return_value=return_value )
            ENDIF

         ENDIF
      ENDDO

      nr_fix_variables = fix_ind
!
!--   Variables time axis
!--   These variables will always be written in the time loop

      dimension_ids(1) = did%prt
      dimension_ids(2) = did%time

      var_ind = 0

      DO i=1,size(data_output_pts)
         const_flag = (INDEX(TRIM(data_output_pts(i)),'_const') > 0)
         IF(LEN(TRIM(data_output_pts(i))) > 0 .AND. .NOT. const_flag)    THEN
            var_ind = var_ind+1
            variables(var_ind)%name  = TRIM(data_output_pts(i))

            SELECT CASE (TRIM(variables(var_ind)%name))
               CASE('id')
                  variables(var_ind)%name  = TRIM(data_output_pts(i)) // '_low'
                  variables(var_ind)%units = 'Number'
                  variables(var_ind)%is_integer = .TRUE.
                  IF (part_io%iam_io_pe)   THEN
                     CALL netcdf4_init_variable( 'parallel', file_id, variables(var_ind)%var_id, variables(var_ind)%name,  &
                                                'int32', dimension_ids(1:2), .FALSE., return_value )
                     CALL netcdf4_write_attribute( 'parallel', file_id, variables(var_ind)%var_id, 'units', &
                                                value_char=trim(variables(var_ind)%units), return_value=return_value )
                  ENDIF

                  var_ind = var_ind+1
                  variables(var_ind)%name  = TRIM(data_output_pts(i)) // '_high'
                  variables(var_ind)%units = 'Number'
                  variables(var_ind)%is_integer = .TRUE.
               CASE('particle_nr')
                  variables(var_ind)%units = 'Number'
                  variables(var_ind)%is_integer = .TRUE.
               CASE('class')
                  variables(var_ind)%units = 'Number'
                  variables(var_ind)%is_integer = .TRUE.
               CASE('group')
                  variables(var_ind)%units = 'Number'
                  variables(var_ind)%is_integer = .TRUE.
               CASE('x')
                  variables(var_ind)%units = 'meter'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('y')
                  variables(var_ind)%units = 'meter'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('z')
                  variables(var_ind)%units = 'meter'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('speed_x')
                  variables(var_ind)%units = 'm/s'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('speed_y')
                  variables(var_ind)%units = 'm/s'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('speed_z')
                  variables(var_ind)%units = 'm/s'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('radius')
                  variables(var_ind)%units = 'meter'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('age')
                  variables(var_ind)%units = 'sec'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('age_m')
                  variables(var_ind)%units = 'sec'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('dt_sum')
                  variables(var_ind)%units = 'sec'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('e_m')
                  variables(var_ind)%units = 'Ws'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('weight_factor')
                  variables(var_ind)%units = 'factor'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('aux1')
                  variables(var_ind)%units = 'depend_on_setup'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('aux2')
                  variables(var_ind)%units = 'depend_on_setup'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('rvar1')
                  variables(var_ind)%units = 'depend_on_setup'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('rvar2')
                  variables(var_ind)%units = 'depend_on_setup'
                  variables(var_ind)%is_integer = .FALSE.
               CASE('rvar3')
                  variables(var_ind)%units = 'depend_on_setup'
                  variables(var_ind)%is_integer = .FALSE.
            END SELECT

            IF (part_io%iam_io_pe)   THEN
               IF(variables(var_ind)%is_integer)   THEN
                  CALL netcdf4_init_variable( 'parallel', file_id, variables(var_ind)%var_id, variables(var_ind)%name, 'int32',  &
                                                  dimension_ids(1:2), .FALSE., return_value )
               ELSE
                  CALL netcdf4_init_variable( 'parallel', file_id, variables(var_ind)%var_id, variables(var_ind)%name, 'real32',  &
                                                  dimension_ids(1:2), .FALSE., return_value )
               ENDIF

               CALL netcdf4_write_attribute( 'parallel', file_id, variables(var_ind)%var_id, 'units', &
                  value_char=trim(variables(var_ind)%units), return_value=return_value )
            ENDIF

         ENDIF
      ENDDO

      nr_variables = var_ind

      IF (part_io%iam_io_pe)   THEN
         CALL netcdf4_stop_file_header_definition( 'parallel', file_id, return_value )
      ENDIF

      CALL dop_write_axis

      RETURN

    CONTAINS
       SUBROUTINE dop_write_axis
          IMPLICIT NONE
          INTEGER(iwp)         :: i
          INTEGER,DIMENSION(1) :: bounds_start, value_counts, bounds_origin
          INTEGER(iwp)         :: return_value                !< Return value data_output_netcdf4 .. routines

          INTEGER, POINTER, CONTIGUOUS, DIMENSION(:)   :: prt_val

          bounds_origin = 1
          bounds_start(1) = 1

          IF(myid == 0)   THEN

             ALLOCATE(prt_val(nr_particles_out))
             DO i=1,nr_particles_out
                prt_val(i) = i
             ENDDO
             value_counts(1) = nr_particles_out

             CALL  netcdf4_write_variable('parallel', file_id, var_id%prt, bounds_start, value_counts, bounds_origin, &
                .TRUE., values_int32_1d=prt_val, return_value=return_value)

             DEALLOCATE(prt_val)
          END IF

          RETURN
       END SUBROUTINE dop_write_axis

   END SUBROUTINE dop_netcdf_setup
!
!- write constant variables

   SUBROUTINE dop_write_fixed_variables
      IMPLICIT NONE

      INTEGER(iwp)             :: i                              !
      INTEGER(iwp)             :: return_value                   !

      INTEGER(iwp),DIMENSION(1)   :: bounds_origin, bounds_start, value_counts

      bounds_origin(1) = 1
      bounds_start(1)  = io_start_index
      value_counts(1)  = io_end_index-io_start_index+1


      DO i=1,nr_fix_variables

#if defined( __netcdf4 )
         IF(fix_variables(i)%is_integer)   THEN
            out_buf_i = NF90_FILL_INT
         ELSE
            out_buf_r = NF90_FILL_REAL
         ENDIF
#endif

         CALL dop_fill_out_buf (fix_variables(i))

         CALL  dop_get_remote_particle (fix_variables(i)%is_integer)

         CALL  part_io%sm_node_barrier()
         IF (part_io%iam_io_pe)   THEN
            IF(fix_variables(i)%is_integer)   THEN
               CALL  netcdf4_write_variable('parallel', file_id, fix_variables(i)%var_id, bounds_start, value_counts, &
                                            bounds_origin, .FALSE., values_int32_1d=out_buf_i, return_value=return_value)
            ELSE
               CALL  netcdf4_write_variable('parallel', file_id, fix_variables(i)%var_id, bounds_start, value_counts, &
                                            bounds_origin, .FALSE., values_real32_1d=out_buf_r, return_value=return_value)
            ENDIF
         ENDIF
         CALL  part_io%sm_node_barrier()
      ENDDO

      RETURN
   END SUBROUTINE dop_write_fixed_variables

   SUBROUTINE dop_delete_particle_number
      IMPLICIT NONE

      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: k                           !<
      INTEGER(iwp) :: n                           !<

!
!--   delete inactive particles, i.e. all particles with particle_mask == .FALSE.
!--   get an output particle nr of -1
!kk   Not sure if it is required here or already done in lagrangian_particle_model_mod before this call

      remote_nr_particles = 0
      DO i=nxl,nxr
         DO j=nys,nyn
            DO k=nzb+1,nzt
               DO n=1,prt_count(k,j,i)
                  IF( .NOT. grid_particles(k,j,i)%particles(n)%particle_mask)   THEN
                     grid_particles(k,j,i)%particles(n)%particle_nr = -1
                  END IF
               END DO
            END DO
         ENDDO
      ENDDO

      RETURN
   END SUBROUTINE dop_delete_particle_number

   SUBROUTINE dop_newly_generated_particles
      IMPLICIT NONE

      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: k                           !<
      INTEGER(iwp) :: n                           !<
      INTEGER(iwp) :: n_all                       !< count all particles for MOD function
      INTEGER(iwp) :: nr_new_particle             !<
      INTEGER(iwp) :: particle_nr                 !<
#if defined( __parallel )
      INTEGER(iwp) :: ierr                        !< MPI error code
#endif
      REAL(dp)     :: fcount
      REAL(dp)     :: finc
      INTEGER(iwp),DIMENSION(0:numprocs-1)     :: nr_particles_new_s
      INTEGER(iwp),DIMENSION(0:numprocs-1)     :: nr_particles_new_r
      INTEGER(iwp)                             :: start_new_numbering


!
!--   Count Number of Newly Generated particles
!--   Condition for a newparticle: particle_mask = .TRUE. and particle nr of -1
!
!--   For performance reasons, this subroutine may be combined later with dop_delete_particle_number

      nr_new_particle     = 0

      IF(pts_increment > 1)   THEN
         n_all = 0
         DO i=nxl,nxr
            DO j=nys,nyn
               DO k=nzb+1,nzt
                  DO n=1,prt_count(k,j,i)
                     IF( grid_particles(k,j,i)%particles(n)%particle_mask)   THEN
                        IF(grid_particles(k,j,i)%particles(n)%particle_nr == -1)   THEN
                           IF(MOD(n_all,pts_increment) == 0)   THEN
                              nr_new_particle = nr_new_particle+1
                           ENDIF
                           n_all = n_all+1
                        ENDIF
                     END IF
                  END DO
               END DO
            ENDDO
         ENDDO
      ELSEIF(pts_percentage < 100. )   THEN
         finc   = pts_percentage/100
         fcount = 0.0

         DO i=nxl,nxr
            DO j=nys,nyn
               DO k=nzb+1,nzt
                  DO n=1,prt_count(k,j,i)
                     IF( grid_particles(k,j,i)%particles(n)%particle_mask)   THEN
                        IF(grid_particles(k,j,i)%particles(n)%particle_nr == -1)   THEN
                           fcount = fcount + finc
                           IF(nr_new_particle < int(fcount) )  THEN
                              nr_new_particle = nr_new_particle+1
                           ENDIF
                        ENDIF
                     END IF
                  END DO
               END DO
            ENDDO
         ENDDO

      ELSE
         DO i=nxl,nxr
            DO j=nys,nyn
               DO k=nzb+1,nzt
                  DO n=1,prt_count(k,j,i)
                     IF( grid_particles(k,j,i)%particles(n)%particle_mask)   THEN
                        IF(grid_particles(k,j,i)%particles(n)%particle_nr == -1)   THEN
                           nr_new_particle = nr_new_particle+1
                        ENDIF
                     END IF
                  END DO
               END DO
            ENDDO
         ENDDO
      ENDIF
!
!--   Determine start number of new particles on every thread

      nr_particles_new_s = 0
      nr_particles_new_s(myid) = nr_new_particle

#if defined( __parallel )
      CALL MPI_ALLREDUCE( nr_particles_new_s, nr_particles_new_r, SIZE( nr_particles_new_s ), MPI_INTEGER,    &
                                              MPI_SUM, comm2d, ierr )
#else
      nr_particles_new_r = nr_particles_new_s
#endif
!
!--   Abortion if selected particles from new particle set would exceed particle axis of output
!--   changed by JS
      IF ( ( SUM(nr_particles_new_r) + initial_number_of_active_particles)  >  nr_particles_out ) THEN
         RETURN
      ENDIF

      start_new_numbering = initial_number_of_active_particles+1

      IF ( myid > 0)   THEN
         DO i=1,numprocs-1
            start_new_numbering = start_new_numbering+nr_particles_new_r(i-1)
            IF(myid == i)  EXIT
         ENDDO
      END IF

      initial_number_of_active_particles = initial_number_of_active_particles+SUM(nr_particles_new_r)

      dop_last_active_particle = initial_number_of_active_particles
!
!--   Set number of new particles

      particle_nr         = start_new_numbering
      nr_new_particle     = 0

      IF(pts_increment > 1)   THEN
         n_all = 0
         DO i=nxl,nxr
            DO j=nys,nyn
               DO k=nzb+1,nzt
                  DO n=1,prt_count(k,j,i)
                     IF( grid_particles(k,j,i)%particles(n)%particle_mask)   THEN
                        IF(grid_particles(k,j,i)%particles(n)%particle_nr == -1)   THEN
                           IF(MOD(n_all,pts_increment) == 0)   THEN
                              grid_particles(k,j,i)%particles(n)%particle_nr = particle_nr
                              particle_nr = particle_nr+1
                           ELSE
                              grid_particles(k,j,i)%particles(n)%particle_nr = -2
                           ENDIF
                           n_all = n_all+1
                        ENDIF
                     END IF
                  END DO
               END DO
            ENDDO
         ENDDO
      ELSEIF(pts_percentage < 100. )   THEN
         finc   = pts_percentage/100
         fcount = 0.0

         DO i=nxl,nxr
            DO j=nys,nyn
               DO k=nzb+1,nzt
                  DO n=1,prt_count(k,j,i)
                     IF( grid_particles(k,j,i)%particles(n)%particle_mask)   THEN
                        IF(grid_particles(k,j,i)%particles(n)%particle_nr == -1)   THEN
                           fcount = fcount + finc
                           IF(nr_new_particle < int(fcount) )  THEN
                              grid_particles(k,j,i)%particles(n)%particle_nr = particle_nr
                              particle_nr = particle_nr+1
                              nr_new_particle = nr_new_particle+1
                           ELSE
                              grid_particles(k,j,i)%particles(n)%particle_nr = -2
                           ENDIF
                        ENDIF
                     END IF
                  END DO
               END DO
            ENDDO
         ENDDO
      ELSE
         DO i=nxl,nxr
            DO j=nys,nyn
               DO k=nzb+1,nzt
                  DO n=1,prt_count(k,j,i)
                     IF( grid_particles(k,j,i)%particles(n)%particle_mask)   THEN
                        IF(grid_particles(k,j,i)%particles(n)%particle_nr == -1)   THEN
                            grid_particles(k,j,i)%particles(n)%particle_nr = particle_nr
                            particle_nr = particle_nr+1
                        ENDIF
                     END IF
                  END DO
               END DO
            ENDDO
         ENDDO


      ENDIF



      RETURN
   END SUBROUTINE dop_newly_generated_particles

   SUBROUTINE dop_count_remote_particles
      IMPLICIT NONE

#if defined( __parallel )
      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: k                           !<
      INTEGER(iwp) :: n                           !<
      INTEGER(iwp) :: iop                         !<
      INTEGER(iwp) :: particle_nr
      INTEGER(iwp) :: pe_nr
      INTEGER(iwp) :: win_size
      INTEGER(iwp) :: ierr                        !< MPI error code
      INTEGER(iwp), DIMENSION(0:numprocs-1) :: part_ind

!
!--   Count remote particles

      remote_nr_particles = 0
      DO i=nxl,nxr
         DO j=nys,nyn
            DO k=nzb+1,nzt
               DO n=1,prt_count(k,j,i)
                  particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                  IF ( particle_nr > 0)   THEN
                     DO iop=0,numprocs-1                  !kk this loop has to be optimized
!
!--                     Although the counting is local PE based, the following if is io processor based
!--                     because particles in MPI shared memory do not have to be transfered
                        IF(particle_nr < io_start_index .OR. particle_nr > io_end_index)   THEN
                           IF(particle_nr >= mo_indices(1,iop) .AND. particle_nr <= mo_indices(2,iop))   THEN
                              remote_nr_particles(2,iop) = remote_nr_particles(2,iop)+1
                           ENDIF
                        ENDIF
                     ENDDO
                  END IF
               END DO
            END DO
         ENDDO
      ENDDO

      remote_nr_particles(1,0) = 0
      DO i=1,numprocs-1
         remote_nr_particles(1,i) = remote_nr_particles(1,i-1) + remote_nr_particles(2,i-1)
      END DO

      win_size = sum(remote_nr_particles(2,:))
      CALL dop_alloc_rma_mem (transfer_buffer_i, win_size, win_rma_buf_i)
      CALL dop_alloc_rma_mem (transfer_buffer_r, win_size, win_rma_buf_r)

      CALL MPI_ALLTOALL (remote_nr_particles, 2, MPI_INTEGER, rma_particles, 2, MPI_INTEGER, comm2d, ierr)

!
!--   The particles indices are the same for all output variables during one time step
!--   therefore, the indices are transfered only once here

      part_ind = remote_nr_particles(1,:)
      transfer_buffer_i = -9999

      DO i=nxl,nxr
         DO j=nys,nyn
            DO k=nzb+1,nzt

               DO n=1,prt_count(k,j,i)
                  particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                  IF(particle_nr < io_start_index .OR. particle_nr > io_end_index)   THEN
                     IF (particle_nr > 0)   THEN
                        pe_nr = find_pe_from_particle_nr (particle_nr)
                        transfer_buffer_i (part_ind(pe_nr)) = particle_nr
                        part_ind(pe_nr) = part_ind(pe_nr)+1
                     ENDIF
                  ENDIF
               END DO
            END DO
         ENDDO
      ENDDO

      CALL MPI_Barrier (MPI_COMM_WORLD, ierr)

      CALL dop_get_remote_indices

#endif
      RETURN
   END SUBROUTINE dop_count_remote_particles

!
!- Fill ouput buffer
!
!    local variables values are copied into output buffer
!               local here means local to shared memory group
!    remote variable values are copied into transfer buffer
!               this is done by all threads

   SUBROUTINE dop_fill_out_buf (var)
      IMPLICIT NONE

      TYPE(var_def),INTENT(IN)    :: var

      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: k                           !<
      INTEGER(iwp) :: n                           !<
      INTEGER(idp) :: pval                        !<
      INTEGER(iwp) :: particle_nr
      INTEGER(iwp) :: pe_nr
      INTEGER(iwp) :: local_len
      CHARACTER(len=32) :: local_name
      INTEGER(iwp), DIMENSION(0:numprocs-1) :: part_ind

      part_ind = remote_nr_particles(1,:)
      transfer_buffer_i = -9998
!
!--   filling output buffer is the same for variable name and variable name_const
!--   therefore set local_name without

      local_len = INDEX(TRIM(var%name),'_const')
      IF(local_len == 0)   THEN
         local_name = var%name
      ELSE
         local_name = var%name(1:local_len-1)
      END IF
!
!--   In this subroutine the particles are seperated:
!
!--   All particles which are located in the share memory area of the respective IO thread are copied into
!--   the output buffer. The other output particle are copied into the transfer buffer.

      SELECT CASE (TRIM(local_name))
         CASE('origin_x')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%origin_x
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%origin_x
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('origin_y')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%origin_y
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%origin_y
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('origin_z')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%origin_z
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%origin_z
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('id_low')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        pval = IBITS(grid_particles(k,j,i)%particles(n)%id,0,32)
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_i(particle_nr) = INT(pval,4)
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_i (part_ind(pe_nr)) = INT(pval,4)
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('id_high')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        pval = IBITS(grid_particles(k,j,i)%particles(n)%id,32,32)
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_i(particle_nr) = INT(pval,4)
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_i (part_ind(pe_nr)) = INT(pval,4)
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('particle_nr')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_i(particle_nr) = grid_particles(k,j,i)%particles(n)%particle_nr
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_i (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%particle_nr
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('class')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_i(particle_nr) = grid_particles(k,j,i)%particles(n)%class
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_i (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%class
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('group')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_i(particle_nr) = grid_particles(k,j,i)%particles(n)%group
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_i (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%group
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('x')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%x
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%x
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('y')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%y
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%y
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('z')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%z
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%z
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('speed_x')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%speed_x
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%speed_x
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('speed_y')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%speed_y
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%speed_y
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('speed_z')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%speed_z
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%speed_z
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('radius')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%radius
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%radius
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('age')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%age
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%age
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('age_m')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%age_m
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%age_m
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('dt_sum')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%dt_sum
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%dt_sum
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('e_m')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%e_m
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%e_m
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('weight_factor')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%weight_factor
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%weight_factor
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('aux1')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%aux1
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%aux1
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('aux2')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%aux2
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%aux2
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('rvar1')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%rvar1
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%rvar1
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('rvar2')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%rvar2
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%rvar2
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
         CASE('rvar3')
            DO i=nxl,nxr
               DO j=nys,nyn
                  DO k=nzb+1,nzt
                     DO n=1,prt_count(k,j,i)
                        particle_nr = grid_particles(k,j,i)%particles(n)%particle_nr
                        IF(particle_nr >= io_start_index .AND. particle_nr <= io_end_index)   THEN
                           out_buf_r(particle_nr) = grid_particles(k,j,i)%particles(n)%rvar3
                        ELSE IF (particle_nr > 0)   THEN
                           pe_nr = find_pe_from_particle_nr (particle_nr)
                           transfer_buffer_r (part_ind(pe_nr)) = grid_particles(k,j,i)%particles(n)%rvar3
                           part_ind(pe_nr) = part_ind(pe_nr)+1
                        ENDIF
                     END DO
                  END DO
               ENDDO
            ENDDO
      END SELECT

      RETURN
   END SUBROUTINE dop_fill_out_buf

#if defined( __parallel )
!
!- Get indices (displacement) of remot particles
   SUBROUTINE dop_get_remote_indices

      IMPLICIT NONE

      INTEGER(iwp) :: i                           !<

      INTEGER(iwp) :: bufsize                     !< size of remote indices array
      INTEGER(iwp) :: ind_local                   !< index in remore indices array
      INTEGER(iwp) :: ierr                        !< MPI error code
      INTEGER(KIND=MPI_ADDRESS_KIND)   :: disp    !< displacement in RMA window

      bufsize = SUM(rma_particles(2,:))
      ALLOCATE (remote_indices(0:bufsize))
      remote_indices = -1

      ind_local = 0
      CALL MPI_WIN_FENCE( 0, win_rma_buf_i, ierr )
      DO i=0,numprocs-1
         IF (rma_particles(2,i) > 0)   THEN
            disp = rma_particles(1,i)
            IF(rma_particles(2,i) > 0)   THEN
               CALL MPI_GET (remote_indices(ind_local), rma_particles(2,i), MPI_INTEGER, i, disp,     &
                                           rma_particles(2,i), MPI_INTEGER,  win_rma_buf_i, ierr)
               ind_local = ind_local + rma_particles(2,i)
            END IF
         END IF
      ENDDO
      CALL MPI_WIN_FENCE( 0, win_rma_buf_i, ierr )

      RETURN
   END SUBROUTINE dop_get_remote_indices
#endif


   SUBROUTINE dop_get_remote_particle (is_integer)

      IMPLICIT NONE

      LOGICAL,INTENT(IN)   :: is_integer

#if defined( __parallel )
      INTEGER(iwp) :: i                           !<
      INTEGER(iwp) :: j                           !<
      INTEGER(iwp) :: bufsize                     !< size of remote data array
      INTEGER(iwp) :: ind_local                   !< index in remore indices array
      INTEGER(iwp) :: particle_nr                 !< particle number
      INTEGER(iwp) :: ierr                        !< MPI error code
      INTEGER(KIND=MPI_ADDRESS_KIND)           :: disp        !< displacement in RMA window
      REAL(sp),ALLOCATABLE, DIMENSION(:)       :: rma_buf_r   !< buffer to receive remote data (REAL)
      INTEGER(iwp),ALLOCATABLE, DIMENSION(:)   :: rma_buf_i   !< buffer to receive remote data (INTEGER)

      bufsize   = sum(rma_particles(2,:))
      ind_local = 0
      ALLOCATE (rma_buf_r(0:bufsize-1))
      ALLOCATE (rma_buf_i(0:bufsize-1))

      IF(is_integer)   THEN
         CALL MPI_WIN_FENCE( 0, win_rma_buf_i, ierr )
      ELSE
         CALL MPI_WIN_FENCE( 0, win_rma_buf_r, ierr )
      ENDIF
      DO i=0,numprocs-1
         IF (rma_particles(2,i) > 0)   THEN
            IF(is_integer)   THEN
               disp = rma_particles(1,i)
               CALL MPI_GET (rma_buf_i(ind_local), rma_particles(2,i), MPI_INTEGER, i, disp,     &
                                        rma_particles(2,i), MPI_INTEGER,  win_rma_buf_i, ierr)
               ind_local = ind_local + rma_particles(2,i)
            ELSE
               disp = rma_particles(1,i)
               CALL MPI_GET (rma_buf_r(ind_local), rma_particles(2,i), MPI_real, i, disp,     &
                                         rma_particles(2,i), MPI_real,  win_rma_buf_r, ierr)
               ind_local = ind_local + rma_particles(2,i)
            ENDIF
         END IF
      ENDDO
      IF(is_integer)   THEN
         CALL MPI_WIN_FENCE( 0, win_rma_buf_i, ierr )
      ELSE
         CALL MPI_WIN_FENCE( 0, win_rma_buf_r, ierr )
      ENDIF

      ind_local = 0

      DO i=0,numprocs-1
         IF (rma_particles(2,i) > 0)   THEN
            IF(is_integer)   THEN
!
!--            Copy data from remote PEs into output array
               DO j=0,rma_particles(2,i)-1
                  particle_nr = remote_indices(ind_local)
                  out_buf_i(particle_nr) = rma_buf_i(ind_local)
                  ind_local = ind_local+1
               END DO
            ELSE
!
!--            Copy data from remote PEs into output array

               DO j=0,rma_particles(2,i)-1
                  particle_nr = remote_indices(ind_local)
                  out_buf_r(particle_nr) = rma_buf_r(ind_local)
                  ind_local = ind_local+1
               END DO
            ENDIF
         END IF
      ENDDO

      IF(ALLOCATED(rma_buf_r)) DEALLOCATE(rma_buf_r)
      IF(ALLOCATED(rma_buf_i)) DEALLOCATE(rma_buf_i)
#else
   IF (is_integer)  THEN
   ENDIF
#endif

      RETURN
   END SUBROUTINE dop_get_remote_particle

#if defined( __parallel )
!
!- Allocate memory and cread window for one-sided communication (INTEGER 1-D array)
   SUBROUTINE dop_alloc_rma_mem_i1( array, idim1, win )
      IMPLICIT NONE

      INTEGER(isp), DIMENSION(:), POINTER, INTENT(INOUT)   ::  array  !<
      INTEGER(iwp), INTENT(IN)                             ::  idim1   !<
      INTEGER(iwp), INTENT(OUT)                            ::  win     !<

      INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !< size of RMA window
      INTEGER                        ::  ierr     !< MPI error code

      winsize = max(idim1, 2)

      ALLOCATE(array(0:winsize-1))

      winsize = winsize * isp

      CALL MPI_WIN_CREATE( array, winsize, isp, MPI_INFO_NULL, comm2d, win, ierr )

      array = -1

      CALL MPI_WIN_FENCE( 0, win, ierr )

   END SUBROUTINE dop_alloc_rma_mem_i1
#endif

#if defined( __parallel )
!
!- Allocate memory and cread window for one-sided communication (REAL 1-D array)
   SUBROUTINE dop_alloc_rma_mem_r1( array, idim1, win )
      IMPLICIT NONE

      REAL(sp), DIMENSION(:), POINTER, INTENT(INOUT)     ::  array   !<
      INTEGER(iwp), INTENT(IN)                           ::  idim1   !<
      INTEGER(iwp), INTENT(OUT)                          ::  win     !<

      INTEGER(KIND=MPI_ADDRESS_KIND) ::  winsize  !< size of RMA window
      INTEGER                        ::  ierr     !< MPI error code


      winsize = max(idim1, 2)


      ALLOCATE(array(0:winsize-1))

      winsize = winsize * sp

      CALL MPI_WIN_CREATE( array, winsize, sp, MPI_INFO_NULL, comm2d, win, ierr )

      array = -1.0

      CALL MPI_WIN_FENCE( 0, win, ierr )

   END SUBROUTINE dop_alloc_rma_mem_r1
#endif


   SUBROUTINE deallocate_and_free
      IMPLICIT NONE

#if defined( __parallel )
      INTEGER                        ::  ierr     !< MPI error code
#endif

#if defined( __parallel )
      CALL MPI_Win_free (win_rma_buf_i, ierr)
      CALL MPI_Win_free (win_rma_buf_r, ierr)
#endif
      IF (ALLOCATED(remote_indices)) DEALLOCATE(remote_indices)

      DEALLOCATE(transfer_buffer_i)
      DEALLOCATE(transfer_buffer_r)

      RETURN

   END SUBROUTINE deallocate_and_free

   FUNCTION find_pe_from_particle_nr (particle_nr)  RESULT (pe_nr)
      IMPLICIT NONE

      INTEGER(iwp), INTENT(IN)       :: particle_nr
      INTEGER(iwp)                   :: pe_nr             !<
      INTEGER(iwp)                   :: base              !<
      INTEGER(iwp)                   :: pnr               !<

      IF(irregular_distribubtion)   THEN
         IF(particle_nr <= nr_particles_rest*nr_particles_PE)   THEN
            pe_nr = (particle_nr-1)/nr_particles_PE
         ELSE
            base  = nr_particles_rest*nr_particles_PE
            pnr   = particle_nr - base
            pe_nr = (pnr-1)/(nr_particles_PE-1)
            pe_nr = pe_nr+nr_particles_rest
         ENDIF
      ELSE
         pe_nr = (particle_nr-1)/nr_particles_PE
      ENDIF


!kk   This error test is to detect programming errors. For performance reasons it can be removed in
!kk   the final, stable version

   END FUNCTION find_pe_from_particle_nr

END MODULE data_output_particle_mod
