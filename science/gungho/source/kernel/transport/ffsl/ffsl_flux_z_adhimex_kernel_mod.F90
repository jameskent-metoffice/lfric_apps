!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
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
                                           GH_SCALAR, CELL_COLUMN
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
  type(arg_type) :: meta_args(5) = (/                  &
       arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2v), & ! flux
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2v), & ! dep pts
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! field
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3),  & ! detj
       arg_type(GH_SCALAR, GH_REAL,    GH_READ)        & ! dt
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: ffsl_flux_z_adhimex_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: ffsl_flux_z_adhimex_code
public :: fifth_order_adhimex
public :: gcrk
public :: solve_fifth_order_matrix
public :: fluxdiv
public :: advdiv

contains

!> @brief Computes the mass flux for FFSL using AdHImEx in the z direction.
!> @param[in]     nlayers   Number of layers
!> @param[in,out] flux      The flux to be computed
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
!> @param[in]     field     The field to construct the flux
!> @param[in]     detj      Volume of cells
!> @param[in]     dt        Time step
!> @param[in]     ndf_w2v   Number of degrees of freedom for W2v per cell
!> @param[in]     undf_w2v  Number of unique degrees of freedom for W2v
!> @param[in]     map_w2v   The dofmap for the W2v cell at the base of the column
!> @param[in]     ndf_w3    Number of degrees of freedom for W3 per cell
!> @param[in]     undf_w3   Number of unique degrees of freedom for W3
!> @param[in]     map_w3    The dofmap for the cell at the base of the column
subroutine ffsl_flux_z_adhimex_code( nlayers,    &
                                      flux,      &
                                      dep_dist,  &
                                      field,     &
                                      detj,      &
                                      dt,        &
                                      ndf_w2v,   &
                                      undf_w2v,  &
                                      map_w2v,   &
                                      ndf_w3,    &
                                      undf_w3,   &
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

  ! Internal variables
  integer(kind=i_def) :: k, s, i_s, w2v_idx, w3_idx

  integer(kind=i_def), parameter :: nstages = 5

  real(kind=r_tran)   :: courant(nlayers)             ! Abs Courant number at cell centres
  real(kind=r_tran)   :: cfl_w2v(nlayers+1)           ! Courant number at cell faces
  real(kind=r_tran)   :: cfl_w3(nlayers)           ! Courant number at cell centres
  real(kind=r_tran)   :: implness_w2v(nlayers + 1)    ! implicitness at faces (div 1D)
  real(kind=r_tran)   :: onemimplness_w2v(nlayers + 1)! 1-implicitness at faces (div 1D)
  real(kind=r_tran)   :: ones(nlayers + 1)            ! 1D array of ones (used in fluxdiv)
  real(kind=r_tran)   :: implness_w3(nlayers)         ! implicitness at centres
  real(kind=r_tran)   :: a_ex(nstages, nstages)       ! explicit Butcher tableau
  real(kind=r_tran)   :: a_im(nstages, nstages)       ! implicit Butcher tableau
  real(kind=r_tran)   :: field_s(nlayers)             ! stage field
  real(kind=r_tran)   :: rhs(nlayers)                 ! rhs (b) in Ax=b
  real(kind=r_tran)   :: field_s_w2v(nlayers + 1)     ! interpolated field_s at faces
  real(kind=r_tran)   :: field_s_up(nlayers)          ! interpolated field_s at face above cell
  real(kind=r_tran)   :: field_s_down(nlayers)        ! interpolated field_s at face below cell
  real(kind=r_tran)   :: c_field_s_w2v(nlayers + 1)   ! interpolated field_s*C at faces
  real(kind=r_tran)   :: f_ex_adv(nstages, nlayers)   ! ex advective divergence
  real(kind=r_tran)   :: f_im_adv(nstages, nlayers)   ! im advective divergence
  real(kind=r_tran)   :: f_ex_con(nstages, nlayers)   ! ex conservative divergence
  real(kind=r_tran)   :: f_im_con(nstages, nlayers)   ! im conservative divergence
  real(kind=r_tran)   :: g_ex_adv(nstages, nlayers)   ! ex advective difference
  real(kind=r_tran)   :: g_im_adv(nstages, nlayers)   ! im advective difference
  real(kind=r_tran)   :: rhs_adv(nlayers)                 ! rhs (b) in Ax=b
  real(kind=r_tran)   :: zero
  real(kind=r_tran)   :: detj_upwind


  w2v_idx = map_w2v(1)
  w3_idx = map_w3(1)
  zero = 0.0_r_tran
  ones = 1.0_r_tran

  ! Calculate Courant number and implicitness - assumes uniform vertical grid
  courant(1) = 0.5_r_tran*(ABS(dep_dist(w2v_idx)) &
                               + ABS(dep_dist(w2v_idx + 1)))
  implness_w3(1) = 1.0_r_tran - 1.0_r_tran/(1.0_r_tran + 0.7_r_tran*(MAX( &
                          1.4_r_tran, courant(1)) - 1.4_r_tran))
  implness_w2v(1) = zero
  do k = 1, nlayers - 1
    courant(k + 1) = 0.5_r_tran*(ABS(dep_dist(w2v_idx + k)) + ABS(dep_dist(w2v_idx + k + 1)))
    implness_w3(k + 1) = 1.0_r_tran - 1.0_r_tran/(1.0_r_tran + 0.7_r_tran*(MAX( &
                          1.4_r_tran, courant(k + 1)) - 1.4_r_tran))
    implness_w2v(k + 1) = MAX(implness_w3(k), implness_w3(k+1))
  end do
  implness_w2v(nlayers + 1) = zero

  do k = 1, nlayers+1
     cfl_w2v(k) = dep_dist(w2v_idx + k - 1)
  end do
  do k = 1, nlayers
     cfl_w3(k) = (cfl_w2v(k)+cfl_w2v(k+1))/2.0_r_tran
  end do

  ! ! Calculate Courant number and implicitness - assumes uniform vertical grid
  ! courant(1) = 0.5_r_tran*(ABS(dep_dist(w2v_idx)) &
  !                              + ABS(dep_dist(w2v_idx + 1)))
  ! implness_w3(1) = 1.0_r_tran - 1.0_r_tran/(MAX(1.0_r_tran, courant(1)))
  ! implness_w2v(1) = zero
  ! do k = 1, nlayers - 1
  !   courant(k + 1) = 0.5_r_tran*(ABS(dep_dist(w2v_idx + k)) + ABS(dep_dist(w2v_idx + k + 1)))
  !   implness_w3(k + 1) = 1.0_r_tran - 1.0_r_tran/(MAX(1.0_r_tran, courant(k + 1)))
  !   implness_w2v(k + 1) = MAX(implness_w3(k), implness_w3(k+1))
  ! end do
  ! implness_w2v(nlayers + 1) = zero

  ! implness_w2v = ones
  ! implness_w3 = ones(1 : nlayers)

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

  ! Setting all f_... elements to zero
  f_ex_adv = zero
  f_im_adv = zero
  f_ex_con = zero
  f_im_con = zero
  g_ex_adv = zero
  g_im_adv = zero
  rhs = zero
  rhs_adv = zero
  flux(w2v_idx : w2v_idx + nlayers) = zero

  do s = 1, nstages
    ! Calculate rhs (depends on stage for constancy)
    rhs = field(w3_idx : w3_idx + nlayers - 1)
    rhs_adv = field(w3_idx : w3_idx + nlayers - 1)
    if ( (s == 2) .or. (s == 3) ) then ! advective
      do i_s = 1, s
        rhs = rhs + f_ex_adv(i_s,:)*a_ex(s, i_s) &
                   + f_im_adv(i_s,:)*a_im(s, i_s)
        rhs_adv = rhs_adv + g_ex_adv(i_s,:)*a_ex(s, i_s) &
                   + g_im_adv(i_s,:)*a_im(s, i_s)
      end do
    else ! conservative
      do i_s = 1, s
        rhs = rhs + f_ex_con(i_s,:)*a_ex(s, i_s) &
                   + f_im_con(i_s,:)*a_im(s, i_s)
        rhs_adv = rhs_adv + g_ex_adv(i_s,:)*a_ex(s, i_s) &
                   + g_im_adv(i_s,:)*a_im(s, i_s)
      end do
    end if

    ! Call matrix solver for last stage if required
    if ( (s == 5) .and. (any(implness_w2v /= zero)) ) then
      call gcrk( nlayers, rhs_adv, field_s, field(w3_idx : w3_idx + nlayers - 1), a_im(s,s), &
                 cfl_w2v, implness_w2v)
    else
      if ( (s == 2) .or. (s == 3) ) then
        field_s = rhs_adv
      else
        field_s = rhs_adv
      end if
    end if

    ! Interpolate field_s to faces
    call fifth_order_adhimex( nlayers, &
                    field_s_w2v, field_s, cfl_w2v )
    !c_field_s_w2v = dep_dist(w2v_idx : w2v_idx + nlayers)*field_s_w2v
    c_field_s_w2v = cfl_w2v*field_s_w2v
    field_s_up(:) = field_s_w2v(2:nlayers+1)
    field_s_down(:) = field_s_w2v(1:nlayers)

    ! Calculate various divergences (f_...) based on the new stage field
    ! (new description with fluxdiv and I needed to change adv and con around)
    onemimplness_w2v = ones - implness_w2v
    call fluxdiv(nlayers, c_field_s_w2v, f_ex_con(s,:), onemimplness_w2v)
    ! f_ex_con(s,:) = f_ex_con(s,:)
    call fluxdiv(nlayers, c_field_s_w2v, f_im_con(s,:), implness_w2v)
    ! f_im_con(s,:) = f_im_con(s,:)
    call fluxdiv(nlayers, c_field_s_w2v, f_ex_adv(s,:), ones)
    f_ex_adv(s,:) = (ones(1:nlayers) - implness_w3)*f_ex_adv(s,:)
    call advdiv( nlayers, field_s_up, field_s_down, cfl_w3, g_ex_adv(s,:) )
    g_ex_adv(s,:) = (ones(1:nlayers) - implness_w3)*g_ex_adv(s,:)
    call fluxdiv(nlayers, c_field_s_w2v, f_im_adv(s,:), ones)
    f_im_adv(s,:) = implness_w3*f_im_adv(s,:)
    call advdiv( nlayers, field_s_up, field_s_down, cfl_w3, g_im_adv(s,:) )
    g_im_adv(s,:) = implness_w3*g_im_adv(s,:)

    ! Update total flux
    flux(w2v_idx) = 0.0_r_tran
    do k = 2, nlayers
      detj_upwind = MAX(0.0_r_tran, SIGN(1.0_r_tran, c_field_s_w2v(k)))   &
                    * detj(w3_idx + k - 2)                                &
                    - MIN(0.0_r_tran, SIGN(1.0_r_tran, c_field_s_w2v(k))) &
                    * detj(w3_idx + k - 1)
      flux(w2v_idx + k - 1) = flux(w2v_idx + k - 1) + (a_ex(nstages,s) &
                                           *(ones(k) - implness_w2v(k)) + &
          a_im(nstages,s)*implness_w2v(k))*c_field_s_w2v(k)*detj_upwind/dt
    end do
    flux(w2v_idx + nlayers) = 0.0_r_tran
  end do

end subroutine ffsl_flux_z_adhimex_code


!> @brief Calculates fifth-order interpolated field at faces for AdHImEx
!> @param[in]     nl        Number of layers (needs to be at least 6)
!> @param[in,out] fieldh    The flux to be computed
!> @param[in]     field     The field to construct the flux
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
subroutine fifth_order_adhimex( nl,              &
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

  ! Loop over most of the spatial domain (excluding boundaries)
  do k = 3, nl - 3
    fieldh(k + 1) = (MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(k + 1))) & ! u>=0
      *( 2.0_r_tran*field(k - 2) - 13.0_r_tran*field(k - 1)             &
      + 47.0_r_tran*field(k)                                            &
      + 27.0_r_tran*field(k + 1) - 3.0_r_tran*field(k + 2) )            &
      - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(k + 1)))              & ! u<0
      *( -3.0_r_tran*field(k - 1) + 27.0_r_tran*field(k)                &
      + 47.0_r_tran*field(k + 1) - 13.0_r_tran*field(k + 2)             &
      + 2.0_r_tran*field(k + 3) ))/60.0_r_tran
  end do

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

  ! Upper boundary
  fieldh(nl + 1) = field(nl)                                             ! u<=0
  fieldh(nl) = MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl)))         & ! u>=0
                  * field(nl - 1)                                      &
                  - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl)))    & ! u<0
                  * field(nl) 
  fieldh(nl - 1) = MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl - 1)))  & ! u>=0
                  * field(nl - 2)                                       &
                  - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl - 1))) & ! u<0
                  * field(nl - 1)

  !! Lower boundary
  !fieldh(1) = (21.0_r_tran*field(1) - field(2))/20.0_r_tran              ! u>=0
  !fieldh(2) = (MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(2)))          & ! u>=0
  !      *( 36.0_r_tran*field(1) + 27.0_r_tran*field(2)                 &
  !      - 3.0_r_tran*field(3) )                                        &
  !      - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(2)))               & ! u<0
  !      *( 24.0_r_tran*field(1) + 47.0_r_tran*field(2)                 &
  !      - 13.0_r_tran*field(3) + 2.0_r_tran*field(4) ))/60.0_r_tran
  !fieldh(3) = (MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(3)))          & ! u>=0
  !      *( -11.0_r_tran*field(1) + 47.0_r_tran*field(2)                &
  !      + 27.0_r_tran*field(3) - 3.0_r_tran*field(4) )                 &
  !      - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(3)))               & ! u<0
  !      *( -3.0_r_tran*field(1) + 27.0_r_tran*field(2)                 &
  !      + 47.0_r_tran*field(3) - 13.0_r_tran*field(4)                  &
  !      + 2.0_r_tran*field(5) ))/60.0_r_tran

  !! Upper boundary
  !fieldh(nl + 1) = (-field(nl - 1) + 21.0_r_tran*field(nl))/20.0_r_tran  ! u<=0
  !fieldh(nl) = (MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl)))        & ! u>=0
  !                *( 2.0_r_tran*field(nl - 3)                          &
  !                - 13.0_r_tran*field(nl - 2)                          &
  !                + 47.0_r_tran*field(nl - 1)                          &
  !                + 24.0_r_tran*field(nl) )                            &
  !                - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl)))    & ! u<0
  !                *( -3.0_r_tran*field(nl - 2)                         &
  !                + 27.0_r_tran*field(nl - 1)                          &
  !                + 36.0_r_tran*field(nl) ))/60.0_r_tran
  !fieldh(nl - 1) = (MAX(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl - 1))) & ! u>=0
  !                *( 2.0_r_tran*field(nl - 4)                           &
  !                - 13.0_r_tran*field(nl - 3)                           &
  !                + 47.0_r_tran*field(nl - 2)                           &
  !                + 27.0_r_tran*field(nl - 1)                           &
  !                - 3.0_r_tran*field(nl) )                              &
  !                - MIN(0.0_r_tran, SIGN(1.0_r_tran, dep_dist(nl - 1))) & ! u<0
  !                *( -3.0_r_tran*field(nl - 3)                          &
  !                + 27.0_r_tran*field(nl - 2)                           &
  !                + 47.0_r_tran*field(nl - 1)                           &
  !                - 11.0_r_tran*field(nl)))/60.0_r_tran

end subroutine fifth_order_adhimex


!> @brief     Solves the matrix with restarted GCR
!> @details   Generalised conjugate gradient iterative matrix solver
!> @param[in]     nl        Number of layers
!> @param[in]     rhs       The flux to be computed
!> @param[in,out] field     The field to be computed
!> @param[in]     initialguess     The initial guess
!> @param[in]     a_im      The diagonal Butcher tableau coefficient for the stage
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
!> @param[in]     implness_w2v  Implicitness
subroutine gcrk( nl,                  &
                        rhs,          &
                        field,        &
                        initialguess, &
                        a_im,         &
                        dep_dist,     &
                        implness_w2v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                   ! nlayers
  real(kind=r_tran),   intent(in)    :: rhs(nl)              ! b of Ax=b equation
  real(kind=r_tran),   intent(inout) :: field(nl)            ! x (field to solve for)
  real(kind=r_tran),   intent(in)    :: initialguess(nl)     ! initial guess
  real(kind=r_tran),   intent(in)    :: a_im                 ! diagonal Butcher coefficient
  real(kind=r_tran),   intent(in)    :: dep_dist(nl + 1)     ! Courant number
  real(kind=r_tran),   intent(in)    :: implness_w2v(nl + 1) ! implicitness

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

  tol = 1.0E-6
  mrestart = 20
  zero = 0.0_r_tran
  guess = initialguess

  call solve_fifth_order_matrix( nl, guess, lhs, a_im, dep_dist, implness_w2v)
  r = rhs - lhs
  reltol = tol*SQRT(SUM(rhs*rhs))

  outer: do m = 1, mrestart
    v = zero
    v(1,:) = r

    do j = 1, jiters
      call solve_fifth_order_matrix(nl, v(j,:), Avj, a_im, dep_dist, implness_w2v)
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

      call solve_fifth_order_matrix(nl, r, Ar, a_im, dep_dist, implness_w2v)
      beta = zero
      do i = 1, j ! + 1
        call solve_fifth_order_matrix(nl, v(i,:), Avi, a_im, dep_dist, implness_w2v)
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
      rmx = r2max**0.5 ! root max square error

      if ( rmx < reltol ) exit outer
    end do
  end do outer

  field = guess ! todo: perhaps add warning if not converged

end subroutine gcrk


!> @brief         Applies the fifth-order matrix to a field (guess)
!> @param[in]     nl        Number of layers
!> @param[in]     guess     The initial guess
!> @param[in,out] lhs       Left-hand side result of applying matrix
!> @param[in]     a_im      The diagonal Butcher tableau coefficient for the stage
!> @param[in]     dep_dist  The vertical departure points (signed C at faces)
!> @param[in]     implness_w2v  Implicitness at faces
subroutine solve_fifth_order_matrix( nl,            &
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
  call fifth_order_adhimex(nl, fieldh, guess, dep_dist) ! out=fieldh (interp guess at faces)

  c_fieldh = dep_dist*fieldh

  call fluxdiv(nl, c_fieldh, div, implness_w2v) ! outputs div

  lhs = guess - a_im*div

end subroutine solve_fifth_order_matrix


!> @brief     Outputs the divergence of the input fields
!> @param[in]     nl        Number of layers
!> @param[in]     c_fieldh  Courant number times field at faces
!> @param[in,out] div       The divergence to be computed
!> @param[in]     implfac   Implicitness factor (implness_w2v (im) or onemimplness_w2v (ex))
subroutine fluxdiv( nl,                  &
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


!> @brief     Outputs the advective difference of the input fields
!> @param[in]     nl        Number of layers
!> @param[in]     c_fieldh  Courant number times field at faces
!> @param[in,out] div       The divergence to be computed
!> @param[in]     implfac   Implicitness factor (implness_w2v (im) or onemimplness_w2v (ex))
subroutine advdiv( nl,                  &
                   fieldh_up,     &
                   fieldh_down,     &
                   cfl,       &
                   diff )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nl                ! nlayers
  real(kind=r_tran),   intent(in)    :: fieldh_up(nl)  ! field at above faces
  real(kind=r_tran),   intent(in)    :: fieldh_down(nl)  ! field at below faces
  real(kind=r_tran),   intent(in)    :: cfl(nl)  ! courant at cell centres
  real(kind=r_tran),   intent(inout) :: diff(nl)           ! difference

  ! Internal variables
  integer(kind=i_def) :: k

  ! (assume uniform grid)
  do k = 1, nl
    diff(k) = - fieldh_up(k) + fieldh_down(k)
    diff(k) = cfl(k) * diff(k)
  end do

end subroutine advdiv


end module ffsl_flux_z_adhimex_kernel_mod
