!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   AdHImEx kernel in the z direction.
!> @details Adapted from the Piecewise Constant Method applied to the z direction.
!!          Original file: ffsl_flux_z_constant_kernel_mod.F90

module ffsl_flux_z_adhimex_kernel_mod

use argument_mod,                   only : arg_type,              &
                                           GH_FIELD, GH_REAL,     &
                                           GH_READ, GH_WRITE,     &
                                           GH_LOGICAL, GH_SCALAR, &
                                           CELL_COLUMN
use fs_continuity_mod,              only : W3, W2v
use constants_mod,                  only : r_tran, i_def, l_def, EPS_R_TRAN
use kernel_mod,                     only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: ffsl_flux_z_adhimex_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                  &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! flux
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! detj
       arg_type(GH_SCALAR, GH_REAL,    GH_READ),       & ! dt
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)        & ! monotonicity
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: ffsl_flux_z_adhimex_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: ffsl_flux_z_adhimex_code
public :: fifth_order_interp
public :: gcrk
public :: solve_fifth_order_matrix
public :: fluxdiv
public :: advdiff
public :: fct
public :: adimex_upwind
public :: first_order_interp
public :: solve_first_order_matrix
public :: set_extrema

contains

!> @brief Computes the mass flux for FFSL using AdHImEx in the z direction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] flux      The flux to be computed
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
!> @param[in]     field     The field to construct the flux
!> @param[in]     detj      Volume of cells
!> @param[in]     dt        Time step
!> @param[in]     monotonicity Whether to apply the FCT limiter
!> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine ffsl_flux_z_adhimex_code( nlayers,       &
                                      flux,         &
                                      dep_dist,     &
                                      field,        &
                                      detj,         &
                                      dt,           &
                                      monotonicity, &
                                      ndf_w2v,      &
                                      undf_w2v,     &
                                      map_w2v,      &
                                      ndf_w3,       &
                                      undf_w3,      &
                                      map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: undf_w2v
  integer(kind=i_def), intent(in)    :: ndf_w2v
  integer(kind=i_def), intent(in)    :: undf_w3
  integer(kind=i_def), intent(in)    :: ndf_w3
  real(kind=r_tran),   intent(inout) :: flux(undf_w2v)
  real(kind=r_tran),   intent(in)    :: field(undf_w3)
  real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2v)
  real(kind=r_tran),   intent(in)    :: detj(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  real(kind=r_tran),   intent(in)    :: dt
  logical(kind=l_def), intent(in)    :: monotonicity

  ! Internal variables
  integer(kind=i_def) :: k, s, i_s, w2v_idx, w3_idx
  logical(kind=l_def) :: gcrk_fct

  ! Parameters
  integer(kind=i_def), parameter :: nstages = 5
  real(kind=r_tran),   parameter :: zero = 0.0_r_tran

  ! Internal Fields
  real(kind=r_tran)   :: courant(nlayers)             ! Abs Courant number at cell centres
  real(kind=r_tran)   :: cfl_w3(nlayers)              ! Courant number at cell centres
  real(kind=r_tran)   :: implness_w2v(nlayers + 1)    ! implicitness at faces (div 1D)
  real(kind=r_tran)   :: onemimplness_w2v(nlayers + 1)! 1-implicitness at faces (div 1D)
  real(kind=r_tran)   :: onemimplness_w3(nlayers)     ! 1-implicitness at centres (div 1D)
  real(kind=r_tran)   :: ones(nlayers + 1)            ! 1D array of ones (used in fluxdiv)
  real(kind=r_tran)   :: implness_w3(nlayers)         ! implicitness at centres
  real(kind=r_tran)   :: a_ex(nstages, nstages)       ! explicit Butcher tableau
  real(kind=r_tran)   :: a_im(nstages, nstages)       ! implicit Butcher tableau
  real(kind=r_tran)   :: field_s(nlayers)             ! stage field
  real(kind=r_tran)   :: rhs(nlayers)                 ! rhs (b) in Ax=b
  real(kind=r_tran)   :: field_s_w2v(nlayers + 1)     ! interpolated field_s at faces
  real(kind=r_tran)   :: c_field_s_w2v(nlayers + 1)   ! interpolated field_s*C at faces
  real(kind=r_tran)   :: f_ex_adv(nstages, nlayers)   ! ex advective divergence
  real(kind=r_tran)   :: f_im_adv(nstages, nlayers)   ! im advective divergence
  real(kind=r_tran)   :: detj_upwind

  ! Map indices and constants
  w2v_idx = map_w2v(1)
  w3_idx = map_w3(1)
  ones = 1.0_r_tran
  gcrk_fct = .false.

  ! Calculate absolute Courant number and implicitness - assumes uniform vertical grid
  courant(1) = 0.5_r_tran*(ABS(dep_dist(w2v_idx)) + ABS(dep_dist(w2v_idx + 1)))
  implness_w3(1) = 1.0_r_tran - 1.0_r_tran / &
                   (1.0_r_tran + 0.7_r_tran*(MAX( 1.4_r_tran, courant(1) ) - 1.4_r_tran))
  implness_w2v(1) = zero
  do k = 1, nlayers - 1
    courant(k + 1) = 0.5_r_tran*(ABS(dep_dist(w2v_idx + k)) + ABS(dep_dist(w2v_idx + k + 1)))
    implness_w3(k + 1) = 1.0_r_tran - 1.0_r_tran/(1.0_r_tran + 0.7_r_tran*(MAX( &
                          1.4_r_tran, courant(k + 1)) - 1.4_r_tran))
    implness_w2v(k + 1) = MAX(implness_w3(k), implness_w3(k+1))
  end do
  implness_w2v(nlayers + 1) = zero
  ! One minus implness at W2v and W3
  onemimplness_w2v = ones - implness_w2v
  onemimplness_w3  = (ones(1:nlayers) - implness_w3)
  ! Signed Courant number at W3
  do k = 1, nlayers
     cfl_w3(k) = (dep_dist(w2v_idx + k - 1)+dep_dist(w2v_idx + k))/2.0_r_tran
  end do

  ! Set up Butcher tableau (remember column-major order of reshape)
  a_ex = reshape((/ zero, zero, zero, zero, zero,                               &
                    zero, zero, 1.0_r_tran, 0.25_r_tran, 1.0_r_tran/6.0_r_tran, &
                    zero, zero, zero, 0.25_r_tran, 1.0_r_tran/6.0_r_tran,       &
                    zero, zero, zero, zero, 2.0_r_tran/3.0_r_tran,              &
                    zero, zero, zero, zero, zero /), shape(a_ex))
  a_im = reshape((/ zero, 0.5_r_tran, 0.5_r_tran, 0.5_r_tran, 0.5_r_tran,       &
                    zero, zero, zero, zero, zero,                               &
                    zero, zero, zero, zero, zero,                               &
                    zero, zero, zero, zero, zero,                               &
                    zero, zero, zero, zero, 0.5_r_tran /), shape(a_im))

  ! Setting elements to zero
  f_ex_adv = zero
  f_im_adv = zero
  rhs = zero
  flux(w2v_idx : w2v_idx + nlayers) = zero

  ! Loop over number of RK stages
  do s = 1, nstages
    ! Calculate rhs using advective form
    rhs = field(w3_idx : w3_idx + nlayers - 1)
    do i_s = 1, s
      rhs = rhs + f_ex_adv(i_s,:)*a_ex(s, i_s) + f_im_adv(i_s,:)*a_im(s, i_s)
    end do

    ! Call matrix solver for last stage if required
    if ( (s == 5) .and. (any(implness_w2v /= zero)) ) then
      call gcrk( nlayers, rhs, field_s, field(w3_idx : w3_idx + nlayers - 1), &
                 a_im(s,s), dep_dist(w2v_idx : w2v_idx + nlayers),            &
                 implness_w2v, gcrk_fct)
    else
      field_s = rhs
    end if

    ! Interpolate field_s to faces
    call fifth_order_interp( nlayers, &
                    field_s_w2v, field_s, dep_dist(w2v_idx : w2v_idx + nlayers) )
    c_field_s_w2v = dep_dist(w2v_idx : w2v_idx + nlayers)*field_s_w2v

    ! Calculate various differences (f) based on the new stage field
    call advdiff( nlayers, field_s_w2v(2:nlayers+1), field_s_w2v(1:nlayers), &
                  cfl_w3, f_ex_adv(s,:) )
    f_ex_adv(s,:) = onemimplness_w3 * f_ex_adv(s,:)
    call advdiff( nlayers, field_s_w2v(2:nlayers+1), field_s_w2v(1:nlayers), &
                  cfl_w3, f_im_adv(s,:) )
    f_im_adv(s,:) = implness_w3 * f_im_adv(s,:)

    ! Update total flux
    flux(w2v_idx) = 0.0_r_tran
    do k = 2, nlayers
      detj_upwind = MAX(0.0_r_tran, SIGN(1.0_r_tran, c_field_s_w2v(k)))   &
                    * detj(w3_idx + k - 2)                                &
                    - MIN(0.0_r_tran, SIGN(1.0_r_tran, c_field_s_w2v(k))) &
                    * detj(w3_idx + k - 1)
      flux(w2v_idx + k - 1) = flux(w2v_idx + k - 1) +                     &
                              ( a_ex(nstages,s) * onemimplness_w2v(k) +   &
                                a_im(nstages,s) * implness_w2v(k) )       &
                              * c_field_s_w2v(k) * detj_upwind / dt
    end do
    flux(w2v_idx + nlayers) = 0.0_r_tran
  end do

  if (monotonicity) then
    ! Apply FCT limiter
    call fct( nlayers, flux(w2v_idx : w2v_idx + nlayers), & 
              field(w3_idx : w3_idx + nlayers - 1),       &
              dep_dist(w2v_idx : w2v_idx + nlayers),      &
              detj(w3_idx : w3_idx + nlayers - 1), dt)
  end if

end subroutine ffsl_flux_z_adhimex_code


!> @brief Calculates fifth-order interpolated field at faces for AdHImEx
!> @param[in]     nl        Number of layers (needs to be at least 6)
!> @param[in,out] fieldh    The flux to be computed
!> @param[in]     field     The field to construct the flux
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
subroutine fifth_order_interp( nl,        &
                               fieldh,    &
                               field,     &
                               dep_dist )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                ! nlayers
  real(kind=r_tran),   intent(inout) :: fieldh(nl + 1)    ! interpolated field
  real(kind=r_tran),   intent(in)    :: field(nl)         ! field to interpolate
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)  ! Courant (used for sign)

  ! Internal variables
  integer(kind=i_def) :: k

  ! Lower boundary
  fieldh(1) = field(1)                                                   ! u>=0
  fieldh(2) = MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(2)))           & ! u>=0
        * field(1)                                                     &
        - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(2)))               & ! u<0
        * field(2)
  fieldh(3) = MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(3)))           & ! u>=0
        * field(2)                                                     &
        - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(3)))               & ! u<0
        * field(3)

  ! Loop over most of the spatial domain (excluding boundaries)
  do k = 3, nl - 3
    fieldh(k + 1) = (MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(k + 1))) & ! u>=0
      * ( 2.0_r_tran*field(k - 2) - 13.0_r_tran*field(k - 1)            &
      + 47.0_r_tran*field(k)                                            &
      + 27.0_r_tran*field(k + 1) - 3.0_r_tran*field(k + 2) )            &
      - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(k + 1)))              & ! u<0
      * ( -3.0_r_tran*field(k - 1) + 27.0_r_tran*field(k)               &
      + 47.0_r_tran*field(k + 1) - 13.0_r_tran*field(k + 2)             &
      + 2.0_r_tran*field(k + 3) ))/60.0_r_tran
  end do

  ! Upper boundary
  fieldh(nl + 1) = field(nl)                                              ! u<=0
  fieldh(nl) = MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl)))          & ! u>=0
                  * field(nl - 1)                                       &
                  - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl)))     & ! u<0
                  * field(nl) 
  fieldh(nl - 1) = MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl - 1)))  & ! u>=0
                  * field(nl - 2)                                       &
                  - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl - 1))) & ! u<0
                  * field(nl - 1)

end subroutine fifth_order_interp


!> @brief     Solves the matrix with restarted GCR
!> @details   Generalised conjugate gradient iterative matrix solver
!> @param[in]     nl           Number of layers
!> @param[in]     rhs          b in Ax = b
!> @param[in,out] field        The field to be computed
!> @param[in]     initialguess The initial guess
!> @param[in]     a_im         The diagonal Butcher tableau coefficient for the stage
!> @param[in]     dep_dist     The vertical departure points (signed C at faces)
!> @param[in]     implness_w2v Implicitness
!> @param[in]     gcrk_fct     Logical to determine whether this is for FCT (true) or not
subroutine gcrk( nl,           &
                 rhs,          &
                 field,        &
                 initialguess, &
                 a_im,         &
                 dep_dist,     &
                 implness_w2v, &
                 gcrk_fct)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                   ! nlayers
  real(kind=r_tran),   intent(in)    :: rhs(nl)              ! b of Ax=b equation
  real(kind=r_tran),   intent(inout) :: field(nl)            ! x (field to solve for)
  real(kind=r_tran),   intent(in)    :: initialguess(nl)     ! initial guess
  real(kind=r_tran),   intent(in)    :: a_im                 ! diagonal Butcher coefficient
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)     ! Courant number
  real(kind=r_tran),   intent(in)    :: implness_w2v(nl + 1) ! implicitness
  logical(kind=l_def), intent(in)    :: gcrk_fct        ! solving for FCT or not

  ! Internal variables
  integer(kind=i_def) :: k, m, mrestart, j, i
  real(kind=r_tran)   :: tol, reltol, Avj2_sum, Avi2_sum, alpha, r2max, rmx, zero

  integer(kind=i_def), parameter :: jiters = 5

  real(kind=r_tran) :: v(jiters + 1, nl)
  real(kind=r_tran) :: r(nl)
  real(kind=r_tran) :: r2(nl)
  real(kind=r_tran) :: Avj(nl)
  real(kind=r_tran) :: Avi(nl)
  real(kind=r_tran) :: Ar(nl)
  real(kind=r_tran) :: beta(jiters + 1)
  real(kind=r_tran) :: lhs(nl)
  real(kind=r_tran) :: bv(nl)
  real(kind=r_tran) :: guess(nl)

  zero = 0.0_r_tran
  guess = initialguess

  if (gcrk_fct) then
    ! Solve for Upwind scheme
    mrestart = 100_i_def
    tol = 1.0E-15_r_tran
    call solve_first_order_matrix( nl, guess, lhs, dep_dist, implness_w2v)
  else
    ! Solve for fifth-order scheme
    mrestart = 20_i_def
    tol = 1.0E-6_r_tran
    call solve_fifth_order_matrix( nl, guess, lhs, a_im, dep_dist, implness_w2v)
  end if

  r = rhs - lhs
  reltol = tol*SQRT(SUM(rhs*rhs))

  outer: do m = 1, mrestart
    v = zero
    v(1,:) = r

    do j = 1, jiters
  
      if (gcrk_fct) then
        ! Solve for Upwind scheme
        call solve_first_order_matrix(nl, v(j,:), Avj, dep_dist, implness_w2v)
      else
        ! Solve for fifth-order scheme
        call solve_fifth_order_matrix(nl, v(j,:), Avj, a_im, dep_dist, implness_w2v)
      end if

      Avj2_sum = zero
      do k = 1, nl
        Avj2_sum = Avj2_sum + Avj(k)*Avj(k)
      end do
      Avj2_sum = MAX(Avj2_sum, 1.0E-15_r_tran)
      alpha = zero
      do k = 1, nl
        alpha = alpha + r(k)*Avj(k)
      end do
      alpha = alpha/Avj2_sum
      guess = guess + alpha*v(j,:)
      r = r - alpha*Avj
  
      if (gcrk_fct) then
        ! Solve for Upwind scheme
        call solve_first_order_matrix(nl, r, Ar, dep_dist, implness_w2v)
      else
        ! Solve for fifth-order scheme
        call solve_fifth_order_matrix(nl, r, Ar, a_im, dep_dist, implness_w2v)
      end if

      beta = zero
      do i = 1, j
        if (gcrk_fct) then
          ! Solve for Upwind scheme
          call solve_first_order_matrix(nl, v(i,:), Avi, dep_dist, implness_w2v)
        else
          ! Solve for fifth-order scheme
          call solve_fifth_order_matrix(nl, v(i,:), Avi, a_im, dep_dist, implness_w2v)
        end if
        Avi2_sum = zero
        do k = 1, nl
          Avi2_sum = Avi2_sum + Avi(k)*Avi(k)
        end do
        Avi2_sum = MAX(Avi2_sum, 1.0E-15_r_tran)
        do k = 1, nl
          beta(i) = beta(i) - Ar(k)*Avi(k)
        end do
        beta(i) = beta(i)/Avi2_sum
      end do

      bv = zero
      do i = 1, j ! + 1
        bv = bv + beta(i)*v(i,:)
      end do

      v(j+1,:) = r + bv ! no need to explicitly ignore beyond j+1 elements as beta is zero there

      r2 = r*r
      r2max = MAXVAL(r2)
      rmx = r2max**0.5_r_tran ! root max square error

      if ( rmx < reltol ) exit outer
    end do
  end do outer

  field = guess ! todo: perhaps add warning if not converged

end subroutine gcrk


!> @brief Applies the fifth-order matrix to a field (guess)
!> @param[in]     nl        Number of layers
!> @param[in]     guess     The initial guess
!> @param[in,out] lhs       Left-hand side result of applying matrix
!> @param[in]     a_im      The diagonal Butcher tableau coefficient for the stage
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
!> @param[in]     implness_w2v  Implicitness at faces
subroutine solve_fifth_order_matrix( nl,           &
                                     guess,        &
                                     lhs,          &
                                     a_im,         &
                                     dep_dist,     &
                                     implness_w2v )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                   ! nlayers
  real(kind=r_tran),   intent(in)    :: guess(nl)            ! matrix input field
  real(kind=r_tran),   intent(inout) :: lhs(nl)              ! lhs field after applying matrix
  real(kind=r_tran),   intent(in)    :: a_im                 ! diagonal Butcher coefficient
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)     ! Courant number
  real(kind=r_tran),   intent(in)    :: implness_w2v(nl + 1) ! implicitness

  ! Internal variables
  real(kind=r_tran)   :: fieldh(nl + 1)
  real(kind=r_tran)   :: c_fieldh(nl + 1)
  real(kind=r_tran)   :: div(nl)

  ! Interpolate guess to faces
  call fifth_order_interp(nl, fieldh, guess, dep_dist) ! out=fieldh (interp guess at faces)

  c_fieldh = dep_dist*fieldh

  call fluxdiv(nl, c_fieldh, div, implness_w2v) ! outputs div

  lhs = guess - a_im*div

end subroutine solve_fifth_order_matrix


!> @brief Outputs the divergence of the input fields
!> @param[in]     nl        Number of layers
!> @param[in]     c_fieldh  Courant number times field at faces
!> @param[in,out] div       The divergence to be computed
!> @param[in]     implfac   Implicitness factor (implness_w2v (im) or onemimplness_w2v (ex))
subroutine fluxdiv( nl,           &
                    c_fieldh,     &
                    div,          &
                    implfac )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                ! nlayers
  real(kind=r_tran),   intent(in)    :: c_fieldh(nl + 1)  ! C*field at faces
  real(kind=r_tran),   intent(inout) :: div(nl)           ! divergence
  real(kind=r_tran),   intent(in)    :: implfac(nl + 1)   ! implicitness factor

  ! Internal variables
  integer(kind=i_def) :: k

  ! (assume uniform grid)
  do k = 1, nl
    div(k) = - implfac(k + 1)*c_fieldh(k + 1) + implfac(k)*c_fieldh(k)
  end do

end subroutine fluxdiv


!> @brief Outputs the advective difference of the input fields
!> @param[in]     nl          Number of layers
!> @param[in]     fieldh_up   Field at face above
!> @param[in]     fieldh_down Field at face below
!> @param[in]     cfl         Courant number at cell centre
!> @param[in,out] diff        The difference to be computed
subroutine advdiff( nl,          &
                    fieldh_up,   &
                    fieldh_down, &
                    cfl,         &
                    diff )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl              ! nlayers
  real(kind=r_tran),   intent(in)    :: fieldh_up(nl)   ! field at above faces
  real(kind=r_tran),   intent(in)    :: fieldh_down(nl) ! field at below faces
  real(kind=r_tran),   intent(in)    :: cfl(nl)         ! courant at cell centres
  real(kind=r_tran),   intent(inout) :: diff(nl)        ! difference

  ! Internal variables
  integer(kind=i_def) :: k

  ! (assume uniform grid)
  do k = 1, nl
    diff(k) = - fieldh_up(k) + fieldh_down(k)
    diff(k) = cfl(k) * diff(k)
  end do

end subroutine advdiff

!> @brief Limits the high-order flux with flux-corrected transport (Zalesak 1979)
!> @param[in]     nl        Number of layers
!> @param[in,out] flux      High-order flux to be limited 
!> @param[in]     field     Field at previous time step, needed for low-order solution
!> @param[in]     dep_dist  Courant number at faces
!> @param[in]     detj      Cell volume
!> @param[in]     dt        Time step
subroutine fct( nl,       &
                flux,     &
                field,    &
                dep_dist, &
                detj,     &
                dt )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                ! nlayers
  real(kind=r_tran),   intent(inout) :: flux(nl + 1)      ! high-order flux
  real(kind=r_tran),   intent(in)    :: field(nl)         ! previous field
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)  ! Courant
  real(kind=r_tran),   intent(in)    :: detj(nl)          ! volume
  real(kind=r_tran),   intent(in)    :: dt                ! time step

  ! Internal variables
  integer(kind=i_def) :: k
  real(kind=r_tran)   :: field_lo(nl)      ! low-order solution
  real(kind=r_tran)   :: flux_lo(nl + 1)   ! low-order flux
  real(kind=r_tran)   :: min_allowed(nl)   ! minimum allowable values
  real(kind=r_tran)   :: max_allowed(nl)   ! maximum allowable values
  real(kind=r_tran)   :: corr(nl + 1)      ! flux correction (high-order - low-order)
  real(kind=r_tran)   :: qp(nl)
  real(kind=r_tran)   :: qm(nl)
  real(kind=r_tran)   :: pp(nl)
  real(kind=r_tran)   :: pm(nl)
  real(kind=r_tran)   :: rp(nl)
  real(kind=r_tran)   :: rm(nl)
  real(kind=r_tran)   :: lim(nl + 1)

  ! Calculate low-order solution (AdImEx upwind with 1-1/(2C))
  call adimex_upwind(nl, field_lo, flux_lo, field, dep_dist, detj, dt)

  ! Calculate allowable extrema
  call set_extrema(nl, min_allowed, max_allowed, field_lo, field, dep_dist)

  !============ FCT algorithm (corrects the flux) ============

  corr = flux - flux_lo ! flux has units field*w*dx*dy

  ! Calculate allowable mass in/out for max rise and fall
  qp = detj*(max_allowed - field_lo)
  qm = detj*(field_lo - min_allowed)

  ! Calculate in/out fluxes at cell centers
  do k = 1, nl
    pp(k) = dt*MAX(0.0_r_tran, corr(k)) - MIN(0.0_r_tran, corr(k + 1))
    pm(k) = dt*MAX(0.0_r_tran, corr(k + 1)) - MIN(0.0_r_tran, corr(k))
  end do

  ! Calculate ratios of allowable (Q) to existing high-order (P) fluxes
  do k = 1, nl
    if (pp(k) .gt. 1.0E-15_r_tran) then
      rp(k) = MIN(1.0_r_tran, qp(k)/MAX(pp(k),1.0E-15_r_tran))
    else
      rp(k) = 0.0_r_tran
    end if
    if (pm(k) .gt. 1.0E-15_r_tran) then
      rm(k) = MIN(1.0_r_tran, qm(k)/MAX(pm(k),1.0E-15_r_tran))
    else
      rm(k) = 0.0_r_tran
    end if
  end do

  lim(1) = 0.0_r_tran
  do k = 2, nl
    if (corr(k) .ge. 0.0_r_tran) then
      lim(k) = MIN(rp(k), rm(k - 1))
    else
      lim(k) = MIN(rm(k), rp(k - 1))
    end if
  end do
  lim(nl + 1) = 0.0_r_tran

  flux = flux_lo + lim*corr

end subroutine fct


!> @brief Calculates first-order AdImEx upwind solution and fluxes
!> @details Used for FCT and assuming divergent flow (implness=1-1/(2C))
!> @param[in]     nl        Number of layers
!> @param[in,out] field_lo  AdImEx upwind solution
!> @param[in,out] flux_lo   Flux that gives the AdImEx upwind solution
!> @param[in]     field     Field at previous time step
!> @param[in]     dep_dist  Courant number at faces
!> @param[in]     detj      Cell-centred volume
!> @param[in]     dt        Time step
subroutine adimex_upwind( nl,        &
                           field_lo, &
                           flux_lo,  &
                           field,    &
                           dep_dist, &
                           detj,     &
                           dt )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                ! nlayers
  real(kind=r_tran),   intent(inout) :: field_lo(nl)      ! low-order field
  real(kind=r_tran),   intent(inout) :: flux_lo(nl + 1)   ! low-order flux
  real(kind=r_tran),   intent(in)    :: field(nl)         ! previous field
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)  ! Courant
  real(kind=r_tran),   intent(in)    :: detj(nl)          ! cell volume
  real(kind=r_tran),   intent(in)    :: dt                ! time step

  ! Internal variables
  integer(kind=i_def)  :: k
  real(kind=r_tran)    :: zero
  real(kind=r_tran)    :: ones(nl + 1)
  real(kind=r_tran)    :: courant(nl)
  real(kind=r_tran)    :: implness_1st_w3(nl)
  real(kind=r_tran)    :: implness_1st_w2v(nl + 1)
  real(kind=r_tran)    :: rhs(nl)
  real(kind=r_tran)    :: fieldh_ex(nl + 1)
  real(kind=r_tran)    :: c_fieldh_ex(nl + 1)
  real(kind=r_tran)    :: div_ex(nl)
  logical(kind=l_def)  :: bool_gcrk_fct

  zero = 0.0_r_tran
  ones = 1.0_r_tran

  ! Calculate Courant number and implicitness - assumes uniform vertical grid
  courant(1) = 0.5_r_tran*(ABS(dep_dist(1)) + ABS(dep_dist(2)))
  implness_1st_w3(1) = 1.0_r_tran - 1.0_r_tran/MAX(1.0_r_tran, 2.0_r_tran*courant(1))
  implness_1st_w2v(1) = zero
  do k = 1, nl - 1
    courant(k + 1) = 0.5_r_tran*(ABS(dep_dist(k + 1)) + ABS(dep_dist(k + 2)))
    implness_1st_w3(k + 1) = 1.0_r_tran - 1.0_r_tran/MAX(1.0_r_tran, 2.0_r_tran*courant(k + 1))
    implness_1st_w2v(k + 1) = MAX(implness_1st_w3(k), implness_1st_w3(k+1))
  end do
  implness_1st_w2v(nl + 1) = zero

  ! Calculate rhs (explicit part)
  call first_order_interp(nl, fieldh_ex, field, dep_dist) ! out=fieldh_ex (interp field at faces)

  c_fieldh_ex = dep_dist*fieldh_ex

  ! Calculate explicit flux divergence
  call fluxdiv(nl, c_fieldh_ex, div_ex, 1.0_r_tran - implness_1st_w2v)

  ! Calculate first-order explicit RHS
  rhs = field + div_ex

  ! Solve matrix and find field_lo
  if (any(implness_1st_w2v /= zero)) then
    bool_gcrk_fct = .true.
    call gcrk( nl, rhs, field_lo, field, zero, dep_dist, implness_1st_w2v, bool_gcrk_fct)
  else
    field_lo = rhs
  end if

  ! Find flux_lo based on field_lo
  flux_lo(1) = zero
  do k = 2, nl ! detj index points to upwind cell
    flux_lo(k) = (MAX(zero, dep_dist(k))*(ones(k) - implness_1st_w2v(k))*field(k - 1)         &
                 + MAX(zero, dep_dist(k))*implness_1st_w2v(k)*field_lo(k - 1))*detj(k - 1)/dt &
                 + (MIN(zero, dep_dist(k))*(ones(k) - implness_1st_w2v(k))*field(k)           &
                 + MIN(zero, dep_dist(k))*implness_1st_w2v(k)*field_lo(k))*detj(k)/dt
  end do
  flux_lo(nl + 1) = zero

end subroutine adimex_upwind


!> @brief Calculates first-order interpolated field at faces for AdImEx upwind (for FCT)
!> @param[in]     nl        Number of layers (needs to be at least 6)
!> @param[in,out] fieldh    The flux to be computed
!> @param[in]     field     The field to construct the flux
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
subroutine first_order_interp( nl,     &
                               fieldh, &
                               field,  &
                               dep_dist )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                ! nlayers
  real(kind=r_tran),   intent(inout) :: fieldh(nl + 1)    ! interpolated field
  real(kind=r_tran),   intent(in)    :: field(nl)         ! field to interpolate
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)  ! Courant (used for sign)

  ! Internal variables
  integer(kind=i_def) :: k

  fieldh(1) = field(1) ! lower boundary

  ! Loop over most of the spatial domain (excluding boundaries)
  do k = 1, nl - 1
    fieldh(k + 1) = MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(k + 1)))*field(k) &   ! u>=0
                  - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(k + 1)))*field(k + 1) ! u<0
  end do

  fieldh(nl + 1) = field(nl) ! upper boundary

end subroutine first_order_interp


!> @brief Applies the first-order matrix to a field (guess)
!> @param[in]     nl        Number of layers
!> @param[in]     guess     The initial guess
!> @param[in,out] lhs       Left-hand side result of applying matrix
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
!> @param[in]     implness_1st_w2v  Implicitness at faces
subroutine solve_first_order_matrix( nl,           &
                                     guess,        &
                                     lhs,          &
                                     dep_dist,     &
                                     implness_1st_w2v )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                       ! nlayers
  real(kind=r_tran),   intent(in)    :: guess(nl)                ! matrix input field
  real(kind=r_tran),   intent(inout) :: lhs(nl)                  ! lhs field after applying matrix
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)         ! Courant number
  real(kind=r_tran),   intent(in)    :: implness_1st_w2v(nl + 1) ! implicitness


  ! Internal variables
  real(kind=r_tran)   :: fieldh(nl + 1)
  real(kind=r_tran)   :: c_fieldh(nl + 1)
  real(kind=r_tran)   :: div(nl)

  ! Interpolate guess to faces
  call first_order_interp(nl, fieldh, guess, dep_dist) ! out=fieldh (interp guess at faces)

  c_fieldh = dep_dist*fieldh

  call fluxdiv(nl, c_fieldh, div, implness_1st_w2v) ! outputs div

  lhs = guess - div

end subroutine solve_first_order_matrix


!> @brief Find min and max allowable values in each cell for FCT
!> @param[in]     nl          Number of layers
!> @param[in,out] min_allowed Minimum value allowed in cell
!> @param[in,out] max_allowed Maximum value allowed in cell
!> @param[in]     field_lo    Low-order solution
!> @param[in]     field       Field at previous time step, needed for low-order solution
!> @param[in]     dep_dist    Courant number at faces
subroutine set_extrema( nl,          &
                        min_allowed, &
                        max_allowed, &
                        field_lo,    &
                        field,       &
                        dep_dist )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                ! nlayers
  real(kind=r_tran),   intent(inout) :: min_allowed(nl)   ! min values
  real(kind=r_tran),   intent(inout) :: max_allowed(nl)   ! max values
  real(kind=r_tran),   intent(in)    :: field_lo(nl)      ! low-order field
  real(kind=r_tran),   intent(in)    :: field(nl)         ! previous field
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)  ! Courant

  ! Internal variables
  integer(kind=i_def)  :: k

  if ( MAXVAL(dep_dist) .gt. 1.0_r_tran ) then
    min_allowed(1) = MIN(field_lo(1), field_lo(2))
    max_allowed(1) = MAX(field_lo(1), field_lo(2))
    do k = 2, nl - 1
      min_allowed(k) = MIN(field_lo(k - 1), field_lo(k), field_lo(k + 1))
      max_allowed(k) = MAX(field_lo(k - 1), field_lo(k), field_lo(k + 1))
    end do
    min_allowed(nl) = MIN(field_lo(nl - 1), field_lo(nl))
    max_allowed(nl) = MAX(field_lo(nl - 1), field_lo(nl))
  else
    min_allowed(1) = MIN(field_lo(1), field_lo(2), field(1), field(2))
    max_allowed(1) = MAX(field_lo(1), field_lo(2), field(1), field(2))
    do k = 2, nl - 1
      min_allowed(k) = MIN(field_lo(k - 1), field_lo(k), field_lo(k + 1), &
                           field(k - 1), field(k), field(k + 1))
      max_allowed(k) = MAX(field_lo(k - 1), field_lo(k), field_lo(k + 1), &
                           field(k - 1), field(k), field(k + 1))
    end do
    min_allowed(nl) = MIN(field_lo(nl - 1), field_lo(nl), field(nl - 1), field(nl))
    max_allowed(nl) = MAX(field_lo(nl - 1), field_lo(nl), field(nl - 1), field(nl))
  end if

end subroutine set_extrema

end module ffsl_flux_z_adhimex_kernel_mod
