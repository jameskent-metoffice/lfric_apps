!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Minimal stand-in for the MPI module, for testing the CaMa-Flood
!>        OASIS driver without MPI (see tests/run_stub_test.sh).
module mpi

  implicit none

  private

  public :: mpi_comm_size

contains

  subroutine mpi_comm_size(comm, size, ierror)
    integer, intent(in)  :: comm
    integer, intent(out) :: size
    integer, intent(out) :: ierror
    size   = 1
    ierror = 0
  end subroutine mpi_comm_size

end module mpi
