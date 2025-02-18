!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the horizontal mass flux for FFSL, with special treatment
!!       at edges of cubed-sphere panels. This is intended for only a single
!!       MPI rank, as if-checks are embedded deep in the code to identify when
!!       a calculation has crossed particular panel boundaries.
!> @details This kernel computes the flux in the x and y directions. A choice of
!!          constant, Nirvana, or PPM is used to compute the edge reconstruction.
!!          This is multiplied by the fractional velocity to give
!!          the flux. For CFL > 1 the field values multiplied by the mass
!!          are summed between the flux point and the departure cell.
!!          Special treatment is used at the edges of cubed-sphere panels, to
!!          shift the reconstruction to improve accuracy.
!!          This kernel is used for the FFSL horizontal transport scheme.
!!
!> @note This kernel only works when field is a W3 field at lowest order since
!!       it is assumed that ndf_w3 = 1 with stencil_map(1,:) containing the
!!       relevant dofmaps.

module ffsl_flux_xy_1rank_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_REAL,         &
                                    CELL_COLUMN, GH_WRITE,     &
                                    GH_READ, GH_SCALAR,        &
                                    STENCIL, X1D, GH_INTEGER,  &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    Y1D, GH_LOGICAL
  use constants_mod,         only : i_def, r_tran, r_def, l_def
  use fs_continuity_mod,     only : W3, W2h
  use kernel_mod,            only : kernel_type
  use reference_element_mod, only : W, E, N, S

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: ffsl_flux_xy_1rank_kernel_type
    private
    type(arg_type) :: meta_args(18) = (/                                       &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),          & ! field_for_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(X1D)),          & ! dry_mass_for_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(Y1D)),          & ! field_for_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W3, STENCIL(Y1D)),          & ! dry_mass_for_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! dep_dist
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! frac_dry_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                                STENCIL(X1D)), & ! panel_id_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                                STENCIL(Y1D)), & ! panel_id_y
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ew
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! face selector ns
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ     ),                        & ! high_order_edges
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                        & ! order
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                        & ! monotone
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                        & ! edge_monotone
         arg_type(GH_SCALAR, GH_REAL,    GH_READ     ),                        & ! min_val
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                        & ! extent_size
         arg_type(GH_SCALAR, GH_REAL,    GH_READ     )                         & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ffsl_flux_xy_1rank_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: ffsl_flux_xy_1rank_code

contains

  !> @brief Compute the horizontal fluxes for FFSL, on a single MPI rank.
  !> @param[in]     nlayers             Number of layers
  !> @param[in,out] flux                The output flux
  !> @param[in]     field_for_x         Field to use in evaluating x-flux
  !> @param[in]     stencil_size_x      Local length of field_for_x W3 stencil
  !> @param[in]     stencil_map_x       Dofmap for the field_for_x stencil
  !> @param[in]     dry_mass_for_x      Volume or dry mass field at W3 points
  !!                                    for use in evaluating x-flux
  !> @param[in]     stencil_size_mass_x Local length of mass x-stencil
  !> @param[in]     stencil_map_mass_x  Dofmap for the mass x-stencil
  !> @param[in]     field_for_y         Field to use in evaluating x-flux
  !> @param[in]     stencil_size_y      Local length of field_for_y W3 stencil
  !> @param[in]     stencil_map_y       Dofmap for the field_for_y stencil
  !> @param[in]     dry_mass_for_y      Volume or dry mass field at W3 points
  !!                                    for use in evaluating x-flux
  !> @param[in]     stencil_size_mass_y Local length of mass y-stencil
  !> @param[in]     stencil_map_mass_y  Dofmap for the mass y-stencil
  !> @param[in]     dep_dist            Horizontal departure distances
  !> @param[in]     frac_dry_flux       Fractional part of the dry flux or wind
  !> @param[in]     panel_id_x          Panel ID used for the x-flux
  !> @param[in]     stencil_size_px     Local length of panel_id_x stencil
  !> @param[in]     stencil_map_px      Dofmap for the panel_id_x stencil
  !> @param[in]     panel_id_y          Panel ID used for the y-flux
  !> @param[in]     stencil_size_py     Local length of panel_id_y stencil
  !> @param[in]     stencil_map_py      Dofmap for the panel_id_y stencil
  !> @param[in]     face_selector_ew    2D field indicating which W/E faces to
  !!                                    loop over for this column
  !> @param[in]     face_selector_ns    2D field indicating which N/S faces to
  !!                                    loop over for this column
  !> @param[in]     high_order_edges    Flag to use high order with special edges
  !> @param[in]     order               Order of reconstruction
  !> @param[in]     monotone            Horizontal monotone option for FFSL
  !> @param[in]     edge_monotone       Monotone option for FFSL at panel edges
  !> @param[in]     min_val             Minimum value to enforce when using
  !!                                    quasi-monotone limiter
  !> @param[in]     extent_size         Stencil extent needed for the LAM edge
  !> @param[in]     dt                  Time step
  !> @param[in]     ndf_w2h             Num of DoFs for W2h per cell
  !> @param[in]     undf_w2h            Num of DoFs for W2h in this partition
  !> @param[in]     map_w2h             Map for W2h
  !> @param[in]     ndf_w3              Num of DoFs for W3 per cell
  !> @param[in]     undf_w3             Num of DoFs for W3 in this partition
  !> @param[in]     map_w3              Map for W3
  !> @param[in]     ndf_pid             Number of DoFs for panel ID per cell
  !> @param[in]     undf_pid            Num DoFs for this partition for panel ID
  !> @param[in]     map_pid             Map for panel ID
  !> @param[in]     ndf_w3_2d           Num of DoFs for 2D W3 per cell
  !> @param[in]     undf_w3_2d          Num of DoFs for this partition for 2D W3
  !> @param[in]     map_w3_2d           Map for 2D W3
  subroutine ffsl_flux_xy_1rank_code( nlayers,             &
                                      flux,                &
                                      field_for_x,         &
                                      stencil_size_x,      &
                                      stencil_map_x,       &
                                      dry_mass_for_x,      &
                                      stencil_size_mass_x, &
                                      stencil_map_mass_x,  &
                                      field_for_y,         &
                                      stencil_size_y,      &
                                      stencil_map_y,       &
                                      dry_mass_for_y,      &
                                      stencil_size_mass_y, &
                                      stencil_map_mass_y,  &
                                      dep_dist,            &
                                      frac_dry_flux,       &
                                      panel_id_x,          &
                                      stencil_size_px,     &
                                      stencil_map_px,      &
                                      panel_id_y,          &
                                      stencil_size_py,     &
                                      stencil_map_py,      &
                                      face_selector_ew,    &
                                      face_selector_ns,    &
                                      high_order_edges,    &
                                      order,               &
                                      monotone,            &
                                      edge_monotone,       &
                                      min_val,             &
                                      extent_size,         &
                                      dt,                  &
                                      ndf_w2h,             &
                                      undf_w2h,            &
                                      map_w2h,             &
                                      ndf_w3,              &
                                      undf_w3,             &
                                      map_w3,              &
                                      ndf_pid,             &
                                      undf_pid,            &
                                      map_pid,             &
                                      ndf_w3_2d,           &
                                      undf_w3_2d,          &
                                      map_w3_2d )

    use subgrid_horizontal_support_mod, only: horizontal_nirvana_case,     &
                                              horizontal_ppm_case,         &
                                              nirvana_special_edge,        &
                                              linear_special_edge,         &
                                              fourth_order_special_edge,   &
                                              fourth_nirvana_special_edge, &
                                              monotonic_horizontal_edge,   &
                                              subgrid_quadratic_recon

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: undf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_w3_2d
    integer(kind=i_def), intent(in) :: ndf_pid
    integer(kind=i_def), intent(in) :: undf_pid
    integer(kind=i_def), intent(in) :: stencil_size_x
    integer(kind=i_def), intent(in) :: stencil_size_mass_x
    integer(kind=i_def), intent(in) :: stencil_size_y
    integer(kind=i_def), intent(in) :: stencil_size_mass_y
    integer(kind=i_def), intent(in) :: stencil_size_px
    integer(kind=i_def), intent(in) :: stencil_size_py

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_pid(ndf_pid)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: map_w3_2d(ndf_w3_2d)
    integer(kind=i_def), intent(in) :: stencil_map_x(ndf_w3,stencil_size_x)
    integer(kind=i_def), intent(in) :: stencil_map_mass_x(ndf_w3,stencil_size_mass_x)
    integer(kind=i_def), intent(in) :: stencil_map_y(ndf_w3,stencil_size_y)
    integer(kind=i_def), intent(in) :: stencil_map_mass_y(ndf_w3,stencil_size_mass_y)
    integer(kind=i_def), intent(in) :: stencil_map_px(ndf_pid,stencil_size_px)
    integer(kind=i_def), intent(in) :: stencil_map_py(ndf_pid,stencil_size_py)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: flux(undf_w2h)
    real(kind=r_tran),   intent(in)    :: field_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: field_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_for_x(undf_w3)
    real(kind=r_tran),   intent(in)    :: dry_mass_for_y(undf_w3)
    real(kind=r_tran),   intent(in)    :: dep_dist(undf_w2h)
    real(kind=r_tran),   intent(in)    :: frac_dry_flux(undf_w2h)
    logical(kind=l_def), intent(in)    :: high_order_edges
    integer(kind=i_def), intent(in)    :: order
    integer(kind=i_def), intent(in)    :: monotone
    integer(kind=i_def), intent(in)    :: edge_monotone
    real(kind=r_tran),   intent(in)    :: min_val
    integer(kind=i_def), intent(in)    :: extent_size
    real(kind=r_tran),   intent(in)    :: dt
    real(kind=r_def),    intent(in)    :: panel_id_x(undf_pid)
    real(kind=r_def),    intent(in)    :: panel_id_y(undf_pid)
    integer(kind=i_def), intent(in)    :: face_selector_ew(undf_w3_2d)
    integer(kind=i_def), intent(in)    :: face_selector_ns(undf_w3_2d)

    ! Integers
    integer(kind=i_def) :: local_dofs_x(2), local_dofs_y(2)
    integer(kind=i_def) :: dof_iterator, order_case, special_edge_case
    integer(kind=i_def) :: stencil_idx, col_idx, pid_idx
    integer(kind=i_def) :: sign_displacement, int_displacement
    integer(kind=i_def) :: rel_idx, rel_dep_cell_idx
    integer(kind=i_def) :: dof_offset, sign_offset, edge_offset
    integer(kind=i_def) :: k, j, idx_3d
    integer(kind=i_def) :: stencil_half, lam_edge_size, recon_size
    integer(kind=i_def) :: mono_option
    integer(kind=i_def) :: owned_panel, panel_edge

    ! Reals
    real(kind=r_tran)   :: displacement, frac_dist
    real(kind=r_tran)   :: recon_field, int_mass
    real(kind=r_tran)   :: field_edge_left, field_edge_right

    real(kind=r_tran),   allocatable :: field_local(:)
    integer(kind=i_def), allocatable :: ipanel_local(:)

    ! x-direction
    local_dofs_x = (/ W, E /)
    local_dofs_y = (/ S, N /)

    ! Set size of array for reconstruction
    if (high_order_edges) then
      ! Increase stencil size by one for special edges
      recon_size  = 3 + 2*order
      ! The order_case determines which 1D scheme to use
      ! and is set equal to the order argument
      order_case  = order
      edge_offset = 2
    else
      ! Use same stencil size as without special edges
      recon_size  = 3 + 2*(order-1)
      ! The order_case determines which 1D scheme to use
      ! and the +10 is to indicate we reduce the order
      ! at special edges
      order_case  = order + 10
      edge_offset = 1
    end if
    allocate(field_local(recon_size))
    allocate(ipanel_local(recon_size))

    ! panel_id_x and panel_id_y are identical fields, so we can use either to
    ! determine the panel for this column
    owned_panel = INT(panel_id_x(map_pid(1)), i_def)

    ! Get size the stencil should be to check if we are at the edge of a LAM domain
    lam_edge_size = 2_i_def*extent_size+1_i_def

    ! = X Calculation ======================================================== !
    stencil_half = (stencil_size_x + 1_i_def) / 2_i_def

    if ( lam_edge_size > stencil_size_x ) then

      ! At edge of LAM, so set output to zero ----------------------------------
      do k = 0, nlayers - 1
        do dof_iterator = 1, face_selector_ew(map_w3_2d(1))
          flux( map_w2h(local_dofs_x(dof_iterator)) + k ) = 0.0_r_tran
        end do
      end do

    else

      ! Not at edge of LAM so compute fluxes -----------------------------------

      ! Loop over the x direction dofs to compute flux at each dof
      do dof_iterator = 1, face_selector_ew(map_w3_2d(1))

        ! Pull out index to avoid multiple indirections
        idx_3d = map_w2h(local_dofs_x(dof_iterator))

        ! Set a local offset, dependent on the face we are looping over
        select case (local_dofs_x(dof_iterator))
        case ( W, S )
          dof_offset = 0
        case ( E, N )
          dof_offset = 1
        end select

        ! Loop over vertical levels
        do k = 0, nlayers - 1

          ! Flux calculation for a single face =================================

          ! Pull out departure point, and separate into integer / frac parts
          displacement = dep_dist(idx_3d + k)
          int_displacement = INT(displacement, i_def)
          frac_dist = displacement - REAL(int_displacement, r_tran)
          sign_displacement = INT(SIGN(1.0_r_tran, displacement))

          ! Set an offset for the stencil index, based on dep point sign
          sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

          ! The local index of the departure cell
          rel_dep_cell_idx = dof_offset - int_displacement + sign_offset - 1

          ! Integer part =======================================================
          int_mass = 0.0_r_tran
          do j = 1, ABS(int_displacement)
            ! If this column has idx 0, find relative index along column of
            ! the departure cell, between - stencil_half and stencil_half
            rel_idx = dof_offset - (j - 1) * sign_displacement + sign_offset - 1

            ! Determine the index in the stencil from rel_idx
            ! e.g. for extent 4:
            ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
            ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
            stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2

            col_idx = stencil_map_x(1,stencil_idx)
            pid_idx = stencil_map_px(1,stencil_idx)

            ! Compute mass over integer number of cells
            panel_edge = 10*owned_panel + INT(panel_id_x(pid_idx), i_def)
            select case (panel_edge)
            case (14, 41, 25, 52, 23, 32, 16, 61, 46, 64, 35, 53)
              ! Calculation has crossed a panel rotation boundary, so use the
              ! opposite fields
              int_mass = int_mass + field_for_y(col_idx + k) * dry_mass_for_y(col_idx + k)
            case default
              int_mass = int_mass + field_for_x(col_idx + k) * dry_mass_for_x(col_idx + k)
            end select
          end do

          ! Fractional part ====================================================
          ! Extract reconstruction data
          do j = 1, recon_size
            ! If this column has idx 0, find relative index along column of
            ! the departure cell, between - stencil_half and stencil_half
            rel_idx = rel_dep_cell_idx + j - order - edge_offset

            ! Determine the index in the stencil from rel_idx
            ! e.g. for extent 4:
            ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
            ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
            stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2

            col_idx = stencil_map_x(1,stencil_idx)
            pid_idx = stencil_map_px(1,stencil_idx)

            ! Populate small local array of field values
            ipanel_local(j) = int(panel_id_x(pid_idx), i_def)
            panel_edge = 10*owned_panel + INT(panel_id_x(pid_idx), i_def)
            select case (panel_edge)
            case (14, 41, 25, 52, 23, 32, 16, 61, 46, 64, 35, 53)
              ! Calculation has crossed a panel rotation boundary, so use the
              ! opposite field
              field_local(j) = field_for_y(col_idx + k)
            case default
              field_local(j) = field_for_x(col_idx + k)
            end select
          end do

          select case ( order_case )
          case ( 0 )
            ! Constant reconstruction
            recon_field = field_local(2)
          case ( 10 )
            ! Constant reconstruction
            recon_field = field_local(1)
          case ( 1 )
            ! Nirvana reconstruction
            ! Get special edge case
            call horizontal_nirvana_case(ipanel_local(2:4), special_edge_case)

            ! Compute edge values using shifted stencil Nirvana at panel edges
            call nirvana_special_edge(field_local, special_edge_case,          &
                                      field_edge_left, field_edge_right)

            ! Apply monotonicity to edges if required
            if (special_edge_case==1_i_def) then
              mono_option = monotone
            else
              mono_option = edge_monotone
            end if
            call monotonic_horizontal_edge(field_local(2:4), mono_option,      &
                                           min_val, field_edge_left, field_edge_right)

            ! Compute reconstruction using field edge values
            ! and quadratic subgrid reconstruction
            call subgrid_quadratic_recon(recon_field, frac_dist,               &
                                         field_local(3), field_edge_left,      &
                                         field_edge_right, mono_option)
          case ( 11 )
            ! Nirvana reconstruction reverting to linear edges at panel boundaries
            ! Get special edge case
            call horizontal_nirvana_case(ipanel_local(1:3), special_edge_case)

            ! Compute edge values reverting to linear interpolation at panel edges
            call linear_special_edge(field_local, special_edge_case,           &
                                     field_edge_left, field_edge_right)

            ! Apply monotonicity to edges if required
            if (special_edge_case==1_i_def) then
              mono_option = monotone
            else
              mono_option = edge_monotone
            end if
            call monotonic_horizontal_edge(field_local(1:3), mono_option,      &
                                          min_val, field_edge_left, field_edge_right)

            ! Compute reconstruction using field edge values
            ! and quadratic subgrid reconstruction
            call subgrid_quadratic_recon(recon_field, frac_dist,               &
                                         field_local(2), field_edge_left,      &
                                         field_edge_right, mono_option)
          case ( 2 )
            ! PPM reconstruction
            ! Get special edge case
            call horizontal_ppm_case(ipanel_local(2:6), special_edge_case)

            ! Compute edge values using shifted stencil PPM at panel edges
            call fourth_order_special_edge(field_local, special_edge_case,     &
                                           field_edge_left, field_edge_right)

            ! Apply monotonicity to edges if required
            if (special_edge_case==1_i_def) then
              mono_option = monotone
            else
              mono_option = edge_monotone
            end if
            call monotonic_horizontal_edge(field_local(3:5), mono_option,      &
                                           min_val, field_edge_left, field_edge_right)

            ! Compute reconstruction using field edge values
            ! and quadratic subgrid reconstruction
            call subgrid_quadratic_recon(recon_field, frac_dist,               &
                                         field_local(4), field_edge_left,      &
                                         field_edge_right, mono_option)

          case ( 12 )
            ! PPM reconstruction reverting to Nirvana edges at panel boundaries
            ! Get special edge case
            call horizontal_ppm_case(ipanel_local(1:5), special_edge_case)

            ! Compute edge values reverting to Nirvana at panel edges
            call fourth_nirvana_special_edge(field_local, special_edge_case,   &
                                             field_edge_left, field_edge_right)

            ! Apply monotonicity to edges if required
            if (special_edge_case==1_i_def) then
              mono_option = monotone
            else
              mono_option = edge_monotone
            end if
            call monotonic_horizontal_edge(field_local(2:4), mono_option,      &
                                           min_val, field_edge_left, field_edge_right)

            ! Compute reconstruction using field edge values
            ! and quadratic subgrid reconstruction
            call subgrid_quadratic_recon(recon_field, frac_dist,               &
                                         field_local(3), field_edge_left,      &
                                         field_edge_right, mono_option)
          end select

          ! Assign flux ========================================================
          flux(idx_3d + k) = (recon_field * frac_dry_flux(idx_3d + k)          &
                              + int_mass * SIGN(1.0_r_tran, displacement)) / dt
        end do ! vertical levels k

      end do ! dof_iterator
    end if

    ! = Y Calculation ======================================================== !
    stencil_half = (stencil_size_y + 1_i_def) / 2_i_def

    if ( lam_edge_size > stencil_size_y ) then

      ! At edge of LAM, so set output to zero ----------------------------------
      do k = 0, nlayers - 1
        do dof_iterator = 1, face_selector_ns(map_w3_2d(1))
          flux( map_w2h(local_dofs_y(dof_iterator)) + k ) = 0.0_r_tran
        end do
      end do

    else

      ! Not at edge of LAM so compute fluxes -----------------------------------

      ! Loop over the y direction dofs to compute flux at each dof
      do dof_iterator = 1, face_selector_ns(map_w3_2d(1))

        ! Pull out index to avoid multiple indirections
        idx_3d = map_w2h(local_dofs_y(dof_iterator))

        ! Set a local offset, dependent on the face we are looping over
        select case (local_dofs_y(dof_iterator))
        case ( W, S )
          dof_offset = 0
        case ( E, N )
          dof_offset = 1
        end select

        ! Loop over vertical levels
        do k = 0, nlayers - 1

          ! Flux calculation for a single face =================================

          ! Pull out departure point, and separate into integer / frac parts
          ! NB: minus signs are needed because the Y1D stencil runs from S to N
          ! but winds blowing S to N have negative sign
          displacement = -dep_dist(idx_3d + k)     ! NB: minus sign
          int_displacement = INT(displacement, i_def)
          frac_dist = displacement - REAL(int_displacement, r_tran)
          sign_displacement = INT(SIGN(1.0_r_tran, displacement))

          ! Set an offset for the stencil index, based on dep point sign
          sign_offset = (1 - sign_displacement) / 2   ! 0 if sign == 1, 1 if sign == -1

          ! The local index of the departure cell
          rel_dep_cell_idx = dof_offset - int_displacement + sign_offset - 1

          ! Integer part =======================================================
          int_mass = 0.0_r_tran
          do j = 1, ABS(int_displacement)
            ! If this column has idx 0, find relative index along column of
            ! the departure cell, between - stencil_half and stencil_half
            rel_idx = dof_offset - (j - 1) * sign_displacement + sign_offset - 1

            ! Determine the index in the stencil from rel_idx
            ! e.g. for extent 4:
            ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
            ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
            stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2

            col_idx = stencil_map_y(1,stencil_idx)
            pid_idx = stencil_map_py(1,stencil_idx)

            ! Compute mass over integer number of cells
            panel_edge = 10*owned_panel + INT(panel_id_y(pid_idx), i_def)
            select case (panel_edge)
            case (14, 41, 25, 52, 23, 32, 16, 61, 46, 64, 35, 53)
              ! Calculation has crossed a panel rotation boundary, so use the
              ! opposite fields
              int_mass = int_mass + field_for_x(col_idx + k) * dry_mass_for_x(col_idx + k)
            case default
              int_mass = int_mass + field_for_y(col_idx + k) * dry_mass_for_y(col_idx + k)
            end select
          end do

          ! Fractional part ====================================================
          ! Extract reconstruction data
          do j = 1, recon_size
            ! If this column has idx 0, find relative index along column of
            ! the departure cell, between - stencil_half and stencil_half
            rel_idx = rel_dep_cell_idx + j - order - edge_offset

            ! Determine the index in the stencil from rel_idx
            ! e.g. for extent 4:
            ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
            ! Stencil has order |  5 |  4 |  3 |  2 |  1 |  6 |  7 |  8 |  9 |
            stencil_idx = 1 + ABS(rel_idx) + (stencil_half - 1)*(1 - SIGN(1, -rel_idx))/2

            col_idx = stencil_map_y(1,stencil_idx)
            pid_idx = stencil_map_py(1,stencil_idx)

            ! Populate small local array of field values
            ipanel_local(j) = int(panel_id_y(pid_idx), i_def)
            panel_edge = 10*owned_panel + INT(panel_id_y(pid_idx), i_def)
            select case (panel_edge)
            case (14, 41, 25, 52, 23, 32, 16, 61, 46, 64, 35, 53)
              ! Calculation has crossed a panel rotation boundary, so use the
              ! opposite field
              field_local(j) = field_for_x(col_idx + k)
            case default
              field_local(j) = field_for_y(col_idx + k)
            end select
          end do

          select case ( order_case )
          case ( 0 )
            ! Constant reconstruction
            recon_field = field_local(2)
          case ( 10 )
            ! Constant reconstruction
            recon_field = field_local(1)
          case ( 1 )
            ! Nirvana reconstruction
            ! Get special edge case
            call horizontal_nirvana_case(ipanel_local(2:4), special_edge_case)

            ! Compute edge values using shifted stencil Nirvana at panel edges
            call nirvana_special_edge(field_local, special_edge_case,          &
                                      field_edge_left, field_edge_right)

            ! Apply monotonicity to edges if required
            if (special_edge_case==1_i_def) then
              mono_option = monotone
            else
              mono_option = edge_monotone
            end if
            call monotonic_horizontal_edge(field_local(2:4), mono_option,      &
                                           min_val, field_edge_left, field_edge_right)

            ! Compute reconstruction using field edge values
            ! and quadratic subgrid reconstruction
            call subgrid_quadratic_recon(recon_field, frac_dist,               &
                                         field_local(3), field_edge_left,      &
                                         field_edge_right, mono_option)
          case ( 11 )
            ! Nirvana reconstruction reverting to linear edges at panel boundaries
            ! Get special edge case
            call horizontal_nirvana_case(ipanel_local(1:3), special_edge_case)

            ! Compute edge values reverting to linear interpolation at panel edges
            call linear_special_edge(field_local, special_edge_case,           &
                                     field_edge_left, field_edge_right)

            ! Apply monotonicity to edges if required
            if (special_edge_case==1_i_def) then
              mono_option = monotone
            else
              mono_option = edge_monotone
            end if
            call monotonic_horizontal_edge(field_local(1:3), mono_option,      &
                                           min_val, field_edge_left, field_edge_right)

            ! Compute reconstruction using field edge values
            ! and quadratic subgrid reconstruction
            call subgrid_quadratic_recon(recon_field, frac_dist,               &
                                         field_local(2), field_edge_left,      &
                                         field_edge_right, mono_option)
          case ( 2 )
            ! PPM reconstruction
            ! Get special edge case
            call horizontal_ppm_case(ipanel_local(2:6), special_edge_case)

            ! Compute edge values using shifted stencil PPM at panel edges
            call fourth_order_special_edge(field_local, special_edge_case,     &
                                           field_edge_left, field_edge_right)

            ! Apply monotonicity to edges if required
            if (special_edge_case==1_i_def) then
              mono_option = monotone
            else
              mono_option = edge_monotone
            end if
            call monotonic_horizontal_edge(field_local(3:5), mono_option,      &
                                           min_val, field_edge_left, field_edge_right)

            ! Compute reconstruction using field edge values
            ! and quadratic subgrid reconstruction
            call subgrid_quadratic_recon(recon_field, frac_dist,               &
                                         field_local(4), field_edge_left,      &
                                         field_edge_right, mono_option)

          case ( 12 )
            ! PPM reconstruction reverting to Nirvana edges at panel boundaries
            ! Get special edge case
            call horizontal_ppm_case(ipanel_local(1:5), special_edge_case)

            ! Compute edge values reverting to Nirvana at panel edges
            call fourth_nirvana_special_edge(field_local, special_edge_case,   &
                                             field_edge_left, field_edge_right)

            ! Apply monotonicity to edges if required
            if (special_edge_case==1_i_def) then
              mono_option = monotone
            else
              mono_option = edge_monotone
            end if
            call monotonic_horizontal_edge(field_local(2:4), mono_option,      &
                                           min_val, field_edge_left, field_edge_right)

            ! Compute reconstruction using field edge values
            ! and quadratic subgrid reconstruction
            call subgrid_quadratic_recon(recon_field, frac_dist,               &
                                         field_local(3), field_edge_left,      &
                                         field_edge_right, mono_option)
          end select

          ! Assign flux ========================================================
          ! NB: minus sign before integer component of flux returns us to the
          ! usual y-direction for the rest of the model
          ! This is not needed for the fractional part, as frac_dry_flux has
          ! the correct sign already
          flux(idx_3d + k) = (recon_field * frac_dry_flux(idx_3d + k)        &
                              - int_mass * SIGN(1.0_r_tran, displacement)) / dt
        end do ! vertical levels k

      end do ! dof_iterator

    end if

    deallocate(field_local)
    deallocate(ipanel_local)

  end subroutine ffsl_flux_xy_1rank_code

end module ffsl_flux_xy_1rank_kernel_mod
