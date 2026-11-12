!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Minimal stand-in for the OASIS3-MCT mod_oasis module.
!>
!> @details Only for testing the CaMa-Flood OASIS driver on a machine without
!>          OASIS3-MCT (see tests/run_stub_test.sh). Received fields are
!>          filled with a uniform value taken from the environment variable
!>          STUB_RUNOFF [kg m-2 s-1] and sent fields are summed and appended
!>          to stub_put.log.
module mod_oasis

  implicit none

  private

  public :: oasis_init_comp, oasis_get_localcomm, oasis_def_partition,          &
            oasis_def_var, oasis_enddef, oasis_get, oasis_put,                 &
            oasis_terminate, oasis_abort

  integer, parameter, public :: oasis_ok = 0, oasis_in = 1, oasis_out = 2,     &
                                oasis_real = 4, oasis_sent = 5,                &
                                oasis_loctrans = 6, oasis_torest = 7,          &
                                oasis_output = 8, oasis_sentout = 9,           &
                                oasis_torestout = 10, oasis_waitgroup = 11,    &
                                oasis_recvd = 12, oasis_fromrest = 13,         &
                                oasis_input = 14, oasis_recvout = 15,          &
                                oasis_fromrestout = 16

  integer, parameter :: wp = selected_real_kind(15, 307)

  interface oasis_get
    module procedure oasis_get_2d
  end interface

  interface oasis_put
    module procedure oasis_put_1d
  end interface

contains

  subroutine oasis_init_comp(comp_id, comp_name, ierror)
    integer,          intent(out) :: comp_id
    character(len=*), intent(in)  :: comp_name
    integer,          intent(out) :: ierror
    comp_id = 1
    ierror  = oasis_ok
    write(*,'(A)') 'stub oasis: init component '//trim(comp_name)
  end subroutine oasis_init_comp

  subroutine oasis_get_localcomm(comm, ierror)
    integer, intent(out) :: comm
    integer, intent(out) :: ierror
    comm   = 0
    ierror = oasis_ok
  end subroutine oasis_get_localcomm

  subroutine oasis_def_partition(part_id, ig_paral, ierror)
    integer, intent(out) :: part_id
    integer, intent(in)  :: ig_paral(3)
    integer, intent(out) :: ierror
    part_id = ig_paral(3)
    ierror  = oasis_ok
  end subroutine oasis_def_partition

  subroutine oasis_def_var(var_id, name, part_id, var_nodims, kinout,          &
                           var_shape, var_type, ierror)
    integer,          intent(out) :: var_id
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: part_id
    integer,          intent(in)  :: var_nodims(2)
    integer,          intent(in)  :: kinout
    integer,          intent(in)  :: var_shape(:)
    integer,          intent(in)  :: var_type
    integer,          intent(out) :: ierror
    var_id = 1
    if (kinout == oasis_out) var_id = 2
    ierror = oasis_ok
    write(*,'(A,I0)') 'stub oasis: define variable '//trim(name)//' size ',    &
                      part_id
  end subroutine oasis_def_var

  subroutine oasis_enddef(ierror)
    integer, intent(out) :: ierror
    ierror = oasis_ok
  end subroutine oasis_enddef

  subroutine oasis_get_2d(var_id, kstep, field, ierror)
    integer,       intent(in)  :: var_id
    integer,       intent(in)  :: kstep
    real(kind=wp), intent(out) :: field(:,:)
    integer,       intent(out) :: ierror
    character(len=64) :: text
    real(kind=wp)     :: value
    integer           :: ios
    value = 1.0e-5_wp
    call get_environment_variable('STUB_RUNOFF', text)
    if (len_trim(text) > 0) read(text, *, iostat=ios) value
    field(:,:) = value
    if (var_id < 0) field(:,:) = 0.0_wp
    ierror = oasis_recvd
  end subroutine oasis_get_2d

  subroutine oasis_put_1d(var_id, kstep, field, ierror)
    integer,       intent(in)  :: var_id
    integer,       intent(in)  :: kstep
    real(kind=wp), intent(in)  :: field(:)
    integer,       intent(out) :: ierror
    integer :: unit
    open(newunit=unit, file='stub_put.log', position='append', action='write')
    write(unit,'(A,I0,A,ES14.6,A,I0,A)') 'time ', kstep, ' sum ',              &
        sum(field), ' nonzero ', count(field /= 0.0_wp), ''
    close(unit)
    ierror = oasis_sent
  end subroutine oasis_put_1d

  subroutine oasis_terminate(ierror)
    integer, intent(out) :: ierror
    ierror = oasis_ok
    write(*,'(A)') 'stub oasis: terminate'
  end subroutine oasis_terminate

  subroutine oasis_abort(comp_id, routine, message)
    integer,          intent(in) :: comp_id
    character(len=*), intent(in) :: routine
    character(len=*), intent(in) :: message
    write(*,'(A)') 'stub oasis ABORT in '//trim(routine)//': '//trim(message)
    error stop 1
  end subroutine oasis_abort

end module mod_oasis
