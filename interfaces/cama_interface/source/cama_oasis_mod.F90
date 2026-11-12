!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief OASIS3-MCT coupling layer for the CaMa-Flood river routing model.
!>
!> @details CaMa-Flood takes the place of the JULES river-routing executable
!>          (river.exe) in the LFRic coupled configuration. The exchanges and
!>          field names are deliberately identical to those of river.exe so
!>          that neither LFRic nor NEMO need to change:
!>
!>            received (2D, on the CaMa-Flood map grid) [kg m-2 s-1]
!>              surf_roff      surface runoff
!>              sub_surf_roff  sub-surface runoff
!>            sent (1D, one value per NEMO river number) [kg s-1]
!>              outflow_per_river
!>
!>          The CaMa-Flood river mouths are summed into the NEMO river
!>          numbers using a map of river numbers on the CaMa-Flood grid.
module cama_oasis_mod

  use mod_oasis,  only: oasis_init_comp, oasis_get_localcomm,                  &
                        oasis_def_partition, oasis_def_var, oasis_enddef,      &
                        oasis_get, oasis_put, oasis_terminate, oasis_abort,    &
                        oasis_in, oasis_out, oasis_real, oasis_ok,             &
                        oasis_sent, oasis_loctrans, oasis_torest,              &
                        oasis_output, oasis_sentout, oasis_torestout,          &
                        oasis_waitgroup, oasis_recvd, oasis_fromrest,          &
                        oasis_input, oasis_recvout, oasis_fromrestout
  use mpi,        only: mpi_comm_size
  use, intrinsic :: iso_fortran_env, only: int64

  implicit none

  private

  public :: cama_cpl_type, map_mouths_to_rivers

  integer, parameter, public :: wp = selected_real_kind(15, 307)

  integer, parameter :: str_len = 256

  ! Component and field names. The field names must match the namcouple.
  character(len=*), parameter :: comp_name          = 'cama'
  character(len=*), parameter :: surf_roff_name     = 'surf_roff'
  character(len=*), parameter :: sub_surf_roff_name = 'sub_surf_roff'
  character(len=*), parameter :: outflow_name       = 'outflow_per_river'

  !> @brief Holds the state of the OASIS coupling of CaMa-Flood
  type :: cama_cpl_type
    private
    integer :: comp_id    = -1
    integer :: id_surf    = -1
    integer :: id_sub     = -1
    integer :: id_outflow = -1
    integer :: nx         = 0
    integer :: ny         = 0
    !> Number of rivers in the NEMO river number space
    integer :: n_rivers   = 0
    !> Coupling period [s]
    integer :: cpl_freq   = 3600
    !> Density used to convert volume flux to mass flux [kg m-3]
    real(kind=wp) :: water_density = 1000.0_wp
    character(len=str_len) :: river_number_file = 'cama_river_number.bin'
    !> NEMO river number of each CaMa-Flood river mouth (0 = not mapped)
    integer, allocatable :: mouth_river(:)
  contains
    procedure :: initialise_oasis
    procedure :: read_config
    procedure :: get_cpl_freq
    procedure :: get_n_rivers
    procedure :: define_coupling
    procedure :: send_outflow
    procedure :: receive_runoff
    procedure :: finalise
    procedure :: abort_run
  end type cama_cpl_type

contains

  !> @brief Initialise OASIS for this component and check it runs serially
  subroutine initialise_oasis(self)

    implicit none

    class(cama_cpl_type), intent(inout) :: self

    integer :: ierror, local_comm, n_ranks

    call oasis_init_comp(self%comp_id, comp_name, ierror)
    if (ierror /= oasis_ok) then
      call oasis_abort(self%comp_id, 'initialise_oasis',                       &
                       'oasis_init_comp failed')
    end if

    call oasis_get_localcomm(local_comm, ierror)
    if (ierror /= oasis_ok) then
      call oasis_abort(self%comp_id, 'initialise_oasis',                       &
                       'oasis_get_localcomm failed')
    end if

    ! The OASIS partitions and the CaMa-Flood build used here are serial
    call mpi_comm_size(local_comm, n_ranks, ierror)
    if (n_ranks /= 1) then
      call oasis_abort(self%comp_id, 'initialise_oasis',                       &
                       'The CaMa-Flood component must run on a single rank')
    end if

  end subroutine initialise_oasis

  !> @brief Read the coupling configuration namelist
  !> @param[in] nml_file  Name of the namelist file
  subroutine read_config(self, nml_file)

    implicit none

    class(cama_cpl_type), intent(inout) :: self
    character(len=*),     intent(in)    :: nml_file

    integer :: unit, ierror

    character(len=str_len) :: river_number_file
    integer                :: cpl_freq
    integer                :: river_count
    real(kind=wp)          :: water_density

    namelist /cama_oasis_nml/ river_number_file, cpl_freq, river_count,        &
                              water_density

    river_number_file = self%river_number_file
    cpl_freq          = self%cpl_freq
    river_count       = self%n_rivers
    water_density     = self%water_density

    open(newunit=unit, file=trim(nml_file), status='old', action='read',       &
         iostat=ierror)
    if (ierror /= 0) then
      call oasis_abort(self%comp_id, 'read_config',                            &
                       'Cannot open namelist file '//trim(nml_file))
    end if
    read(unit, nml=cama_oasis_nml, iostat=ierror)
    close(unit)
    if (ierror /= 0) then
      call oasis_abort(self%comp_id, 'read_config',                            &
                       'Cannot read cama_oasis_nml from '//trim(nml_file))
    end if

    if (river_count < 1 .or. cpl_freq < 1 .or. water_density <= 0.0_wp) then
      call oasis_abort(self%comp_id, 'read_config',                            &
                       'Invalid values in cama_oasis_nml')
    end if

    self%river_number_file = river_number_file
    self%cpl_freq          = cpl_freq
    self%n_rivers          = river_count
    self%water_density     = water_density

  end subroutine read_config

  !> @brief Coupling period [s]
  function get_cpl_freq(self) result(cpl_freq)

    implicit none

    class(cama_cpl_type), intent(in) :: self
    integer                          :: cpl_freq

    cpl_freq = self%cpl_freq

  end function get_cpl_freq

  !> @brief Number of rivers in the NEMO river number space
  function get_n_rivers(self) result(n_rivers)

    implicit none

    class(cama_cpl_type), intent(in) :: self
    integer                          :: n_rivers

    n_rivers = self%n_rivers

  end function get_n_rivers

  !> @brief Map every CaMa-Flood river mouth onto a NEMO river number
  !> @param[in]  river_number  NEMO river number for each CaMa-Flood grid cell
  !> @param[in]  nseqriv       Index of the last river-link cell
  !> @param[in]  nseqall       Index of the last river-mouth cell
  !> @param[in]  seq_x         Grid x index of each cell in the 1D sequence
  !> @param[in]  seq_y         Grid y index of each cell in the 1D sequence
  !> @param[in]  n_rivers      Number of rivers in the NEMO river number space
  !> @param[out] mouth_river   River number of each mouth cell; the array is
  !>                           indexed by position in the 1D sequence
  !> @param[out] n_unmapped    Number of mouth cells without a valid river
  subroutine map_mouths_to_rivers(river_number, nseqriv, nseqall,              &
                                  seq_x, seq_y, n_rivers,                      &
                                  mouth_river, n_unmapped)

    implicit none

    integer, intent(in)  :: river_number(:,:)
    integer, intent(in)  :: nseqriv
    integer, intent(in)  :: nseqall
    integer, intent(in)  :: seq_x(:)
    integer, intent(in)  :: seq_y(:)
    integer, intent(in)  :: n_rivers
    integer, allocatable, intent(out) :: mouth_river(:)
    integer, intent(out) :: n_unmapped

    integer :: iseq, ix, iy, river

    allocate(mouth_river(nseqriv+1:nseqall))
    mouth_river(:) = 0
    n_unmapped = 0

    do iseq = nseqriv + 1, nseqall
      ix = seq_x(iseq)
      iy = seq_y(iseq)
      river = 0
      if (ix >= 1 .and. ix <= size(river_number, 1) .and.                      &
          iy >= 1 .and. iy <= size(river_number, 2)) then
        river = river_number(ix, iy)
      end if
      if (river >= 1 .and. river <= n_rivers) then
        mouth_river(iseq) = river
      else
        n_unmapped = n_unmapped + 1
      end if
    end do

  end subroutine map_mouths_to_rivers

  !> @brief Define the OASIS partitions and variables, and end the
  !>        definition phase
  !> @param[in] nx       Number of CaMa-Flood runoff grid points in x
  !> @param[in] ny       Number of CaMa-Flood runoff grid points in y
  !> @param[in] nseqriv  Index of the last river-link cell
  !> @param[in] nseqall  Index of the last river-mouth cell
  !> @param[in] seq_x    Grid x index of each cell in the 1D sequence
  !> @param[in] seq_y    Grid y index of each cell in the 1D sequence
  subroutine define_coupling(self, nx, ny, nseqriv, nseqall, seq_x, seq_y)

    implicit none

    class(cama_cpl_type), intent(inout) :: self
    integer,              intent(in)    :: nx
    integer,              intent(in)    :: ny
    integer,              intent(in)    :: nseqriv
    integer,              intent(in)    :: nseqall
    integer,              intent(in)    :: seq_x(:)
    integer,              intent(in)    :: seq_y(:)

    integer, allocatable :: river_number(:,:)
    integer :: ig_paral(3)
    integer :: part_2d, part_1d
    integer :: shape_2d(4), shape_1d(2)
    integer :: n_unmapped, ierror

    self%nx = nx
    self%ny = ny

    call read_river_numbers(self, river_number)
    call map_mouths_to_rivers(river_number, nseqriv, nseqall, seq_x, seq_y,    &
                              self%n_rivers, self%mouth_river, n_unmapped)
    write(*,'(A,I0,A,I0,A)') 'cama_oasis: ', nseqall - nseqriv,                &
        ' river mouths, ', n_unmapped, ' not mapped to a NEMO river'
    deallocate(river_number)

    ! Serial partitions: the component runs on a single rank
    ig_paral(1) = 0
    ig_paral(2) = 0
    ig_paral(3) = nx * ny
    call oasis_def_partition(part_2d, ig_paral, ierror)
    call check_ok(self, ierror, 'oasis_def_partition (2d)')

    ig_paral(3) = self%n_rivers
    call oasis_def_partition(part_1d, ig_paral, ierror)
    call check_ok(self, ierror, 'oasis_def_partition (1d)')

    shape_2d = [1, nx, 1, ny]
    shape_1d = [1, self%n_rivers]

    call oasis_def_var(self%id_surf, surf_roff_name, part_2d, [2, 1],          &
                       oasis_in, shape_2d, oasis_real, ierror)
    call check_ok(self, ierror, 'oasis_def_var '//surf_roff_name)
    call oasis_def_var(self%id_sub, sub_surf_roff_name, part_2d, [2, 1],       &
                       oasis_in, shape_2d, oasis_real, ierror)
    call check_ok(self, ierror, 'oasis_def_var '//sub_surf_roff_name)
    call oasis_def_var(self%id_outflow, outflow_name, part_1d, [1, 1],         &
                       oasis_out, shape_1d, oasis_real, ierror)
    call check_ok(self, ierror, 'oasis_def_var '//outflow_name)

    call oasis_enddef(ierror)
    call check_ok(self, ierror, 'oasis_enddef')

    ! A variable that is absent from the namcouple is given an id of -1
    if (self%id_surf < 0 .or. self%id_sub < 0 .or. self%id_outflow < 0) then
      call oasis_abort(self%comp_id, 'define_coupling',                        &
                       'A coupling field was not found in the namcouple')
    end if

  end subroutine define_coupling

  !> @brief Read the NEMO river number on the CaMa-Flood grid. The file holds
  !>        one 4 byte integer per grid point, x varying fastest, with the
  !>        first row in the north (the same layout as CaMa-Flood binary maps).
  subroutine read_river_numbers(self, river_number)

    implicit none

    class(cama_cpl_type),  intent(in)  :: self
    integer, allocatable,  intent(out) :: river_number(:,:)

    integer :: unit, ierror
    integer(kind=int64) :: file_size

    inquire(file=trim(self%river_number_file), size=file_size)
    if (file_size /= 4_int64 * self%nx * self%ny) then
      call oasis_abort(self%comp_id, 'read_river_numbers',                     &
                       'River number file '//trim(self%river_number_file)//    &
                       ' does not have nx*ny 4 byte integers')
    end if

    allocate(river_number(self%nx, self%ny))
    open(newunit=unit, file=trim(self%river_number_file), access='stream',     &
         form='unformatted', status='old', action='read', iostat=ierror)
    if (ierror == 0) read(unit, iostat=ierror) river_number
    if (ierror /= 0) then
      call oasis_abort(self%comp_id, 'read_river_numbers',                     &
                       'Cannot read '//trim(self%river_number_file))
    end if
    close(unit)

  end subroutine read_river_numbers

  !> @brief Send the discharge of every river to the ocean
  !> @param[in] kstep    OASIS time [s]
  !> @param[in] rivout   River outflow of each cell [m3 s-1]
  !> @param[in] fldout   Floodplain outflow of each cell [m3 s-1]
  !> @param[in] nseqriv  Index of the last river-link cell
  !> @param[in] nseqall  Index of the last river-mouth cell
  subroutine send_outflow(self, kstep, rivout, fldout, nseqriv, nseqall)

    implicit none

    class(cama_cpl_type), intent(in) :: self
    integer,              intent(in) :: kstep
    real(kind=wp),        intent(in) :: rivout(:)
    real(kind=wp),        intent(in) :: fldout(:)
    integer,              intent(in) :: nseqriv
    integer,              intent(in) :: nseqall

    real(kind=wp) :: outflow(self%n_rivers)
    integer :: iseq, river, ierror

    outflow(:) = 0.0_wp
    do iseq = nseqriv + 1, nseqall
      river = self%mouth_river(iseq)
      if (river > 0) then
        outflow(river) = outflow(river) +                                      &
                         (rivout(iseq) + fldout(iseq)) * self%water_density
      end if
    end do

    call oasis_put(self%id_outflow, kstep, outflow, ierror)
    if (ierror /= oasis_sent .and. ierror /= oasis_loctrans .and.              &
        ierror /= oasis_torest .and. ierror /= oasis_output .and.              &
        ierror /= oasis_sentout .and. ierror /= oasis_torestout .and.          &
        ierror /= oasis_waitgroup .and. ierror /= oasis_ok) then
      call oasis_abort(self%comp_id, 'send_outflow',                           &
                       'oasis_put failed for '//outflow_name)
    end if

  end subroutine send_outflow

  !> @brief Receive the surface and sub-surface runoff from LFRic
  !> @param[in]  kstep  OASIS time [s]
  !> @param[out] surf   Surface runoff [kg m-2 s-1]
  !> @param[out] sub    Sub-surface runoff [kg m-2 s-1]
  subroutine receive_runoff(self, kstep, surf, sub)

    implicit none

    class(cama_cpl_type), intent(in)  :: self
    integer,              intent(in)  :: kstep
    real(kind=wp),        intent(out) :: surf(:,:)
    real(kind=wp),        intent(out) :: sub(:,:)

    integer :: ierror

    call oasis_get(self%id_surf, kstep, surf, ierror)
    call check_received(self, ierror, surf_roff_name)
    call oasis_get(self%id_sub, kstep, sub, ierror)
    call check_received(self, ierror, sub_surf_roff_name)

  end subroutine receive_runoff

  !> @brief Terminate OASIS
  subroutine finalise(self)

    implicit none

    class(cama_cpl_type), intent(in) :: self

    integer :: ierror

    call oasis_terminate(ierror)
    if (ierror /= oasis_ok) then
      call oasis_abort(self%comp_id, 'finalise', 'oasis_terminate failed')
    end if

  end subroutine finalise

  !> @brief Abort every component of the coupled model
  !> @param[in] message  Reason for the abort
  subroutine abort_run(self, message)

    implicit none

    class(cama_cpl_type), intent(in) :: self
    character(len=*),     intent(in) :: message

    call oasis_abort(self%comp_id, 'cama_oasis_driver', message)

  end subroutine abort_run

  subroutine check_ok(self, ierror, message)

    implicit none

    class(cama_cpl_type), intent(in) :: self
    integer,              intent(in) :: ierror
    character(len=*),     intent(in) :: message

    if (ierror /= oasis_ok) then
      call oasis_abort(self%comp_id, 'cama_oasis_mod', message//' failed')
    end if

  end subroutine check_ok

  subroutine check_received(self, ierror, name)

    implicit none

    class(cama_cpl_type), intent(in) :: self
    integer,              intent(in) :: ierror
    character(len=*),     intent(in) :: name

    if (ierror /= oasis_recvd .and. ierror /= oasis_fromrest .and.             &
        ierror /= oasis_input .and. ierror /= oasis_recvout .and.              &
        ierror /= oasis_fromrestout .and. ierror /= oasis_ok) then
      call oasis_abort(self%comp_id, 'receive_runoff',                         &
                       'oasis_get failed for '//name)
    end if

  end subroutine check_received

end module cama_oasis_mod
