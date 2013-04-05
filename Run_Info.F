MODULE Run_Info
	IMPLICIT NONE
	INCLUDE "mpif.h"
	REAL*8, PUBLIC :: erhx, erhy, erhxc, erhyc, forward_weight
	REAL*8, ALLOCATABLE, PUBLIC :: rhs_x(:), rhs_y(:)
	REAL*8, ALLOCATABLE, PUBLIC :: atax(:,:), atay(:,:), regm(:,:), ataxt(:,:), atayt(:,:), ataxm(:,:), ataym(:,:)
	REAL*8, ALLOCATABLE, PUBLIC :: ataxmt(:,:), ataymt(:,:), rhsxy(:,:), rhsxyt(:,:),rhsxyzero(:,:,:), rhsxyiter(:,:,:)
	REAL*8, ALLOCATABLE, PUBLIC :: L(:,:), Lx(:,:), Ly(:,:), Lz(:,:)

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

	! MPI stuff
	INTEGER, PUBLIC :: myid, numproc
	! myid is the MPI process ID, iam is the BLACS process ID. they are probably
	! the same
	DOUBLE PRECISION :: btime0, btime1 = 0

END MODULE Run_Info
