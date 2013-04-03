Module Run_Info
  Implicit None
  Real*8, Public :: erhx, erhy, erhxc, erhyc, forward_weight
  Real*8, Allocatable, Public :: rhs_x(:), rhs_y(:)
  Real*8, Allocatable, Public :: atax(:,:), atay(:,:), regm(:,:), ataxt(:,:), atayt(:,:), ataxm(:,:), ataym(:,:)
  Real*8, Allocatable, Public :: ataxmt(:,:), ataymt(:,:), rhsxy(:,:), rhsxyt(:,:),rhsxyzero(:,:,:), rhsxyiter(:,:,:)
  Real*8, Allocatable, Public :: L(:,:), Lx(:,:), Ly(:,:), Lz(:,:)

  Integer*8, Public :: n_lons, n_lats
  Integer*8, Public :: i_min, i_max, j_min, j_max
  Integer*8, Public :: inv_dir,myrow,mycol
  Integer, Public :: nprow, npcol, context, iam, nprocs, max_iter
  Integer, Public :: x_block_size, y_block_size, x_block_size_a
  Logical, Public :: Stitch_Boundaries_Only = .False.
  Logical, Public :: Stitch_Regions = .False.
  Logical, Public :: Propagate_Errors = .True.
  Logical, Public :: compute_coarse_kernels = .False.
  Logical, Public :: compute_avg_kernels = .True.
	Logical, Public :: onesided_boundary_stitching = .True.

End Module Run_Info
