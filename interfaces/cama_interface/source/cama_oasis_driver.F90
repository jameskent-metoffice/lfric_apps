!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Runs the CaMa-Flood river routing model as a component of the
!>        OASIS3-MCT coupled LFRic model.
!>
!> @details The exchange sequence follows that of the JULES river.exe: at the
!>          start of every coupling period the river discharge is sent
!>          and the runoff is received, and then CaMa-Flood is advanced over
!>          the coupling period.
!>
!>          CaMa-Flood is configured through its usual namelist file
!>          (input_cmf.nam). It must be configured to read runoff on its own
!>          map grid in units of mm s-1 (kg m-2 s-1):
!>            &NFORCE   LINTERP = .FALSE., DROFUNIT = 1000
!>            &NDIMTIME DT = IFRQ_INP * 3600 = the coupling period
!>
!>          The coupling configuration is read from the namelist file named
!>          by the first command line argument (default: cama_oasis.nml).
program cama_oasis_driver

  use cama_oasis_mod,         only: cama_cpl_type, wp
  use cmf_ctrl_forcing_mod,   only: cmf_forcing_put, drofunit
  use cmf_drv_advance_mod,    only: cmf_drv_advance
  use cmf_drv_control_mod,    only: cmf_drv_input, cmf_drv_init, cmf_drv_end
  use parkind1,               only: jprb, jpim
  use yos_cmf_input,          only: nxin, nyin, dt, lrosplit
  use yos_cmf_map,            only: nseqriv, nseqall, i1seqx, i1seqy
  use yos_cmf_prog,           only: d2rivout, d2fldout
  use yos_cmf_time,           only: nsteps

  implicit none

  type(cama_cpl_type) :: cpl

  character(len=256) :: nml_file
  character(len=256) :: message
  real(kind=wp), allocatable    :: surf(:,:), sub(:,:)
  real(kind=jprb), allocatable  :: runoff(:,:,:)
  integer(kind=jpim) :: istep, istepadv
  integer :: oasis_time, cpl_freq, model_step

  nml_file = 'cama_oasis.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, nml_file)

  call cpl%initialise_oasis()
  call cpl%read_config(trim(nml_file))

  call cmf_drv_input()
  call cmf_drv_init()

  ! Runoff is received in kg m-2 s-1, which is mm s-1 as CaMa-Flood sees it
  if (abs(drofunit - 1000.0_jprb) > 1.0_jprb) then
    call cpl%abort_run('DROFUNIT must be 1000 (runoff in mm s-1)')
  end if

  cpl_freq   = cpl%get_cpl_freq()
  model_step = nint(dt)
  istepadv   = int(cpl_freq / model_step, jpim)
  if (istepadv < 1 .or. istepadv * model_step /= cpl_freq) then
    write(message, '(A,I0,A,I0,A)') 'The coupling period ', cpl_freq,          &
        ' s is not a multiple of the CaMa-Flood time step ', model_step, ' s'
    call cpl%abort_run(trim(message))
  end if

  call cpl%define_coupling(nxin, nyin, nseqriv, nseqall, i1seqx, i1seqy)

  allocate(surf(nxin, nyin), sub(nxin, nyin))
  allocate(runoff(nxin, nyin, 2))

  do istep = 1, nsteps, istepadv
    oasis_time = int(istep - 1) * model_step

    ! Send first and then receive to avoid deadlocks, as the JULES river does
    call cpl%send_outflow(oasis_time, real(d2rivout(:,1), wp),                 &
                          real(d2fldout(:,1), wp), nseqriv, nseqall)
    call cpl%receive_runoff(oasis_time, surf, sub)

    if (lrosplit) then
      runoff(:,:,1) = real(surf, jprb)
      runoff(:,:,2) = real(sub, jprb)
    else
      runoff(:,:,1) = real(surf + sub, jprb)
      runoff(:,:,2) = 0.0_jprb
    end if

    call cmf_forcing_put(runoff)
    call cmf_drv_advance(istepadv)
  end do

  deallocate(surf, sub, runoff)

  call cmf_drv_end()
  call cpl%finalise()

end program cama_oasis_driver
