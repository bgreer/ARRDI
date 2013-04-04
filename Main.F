PROGRAM Main
	USE Data_Management
	USE Run_Info
	USE Fitslib_NF
	USE Matrix_Magic
	IMPLICIT NONE
	INTEGER :: i,cct, j, ii, jj, ierr
	INTEGER*8 :: nnw, regions_per_pe, n_rounds, last_id, r_id, iter_reg
	INTEGER :: nrow, ncol, cxt, serial_inversion, write_time, rs, cs, last_region, left_over
	INTEGER :: iter
	REAL*8 :: del_lat, del_lon
	CHARACTER*120 :: savefile
	CHARACTER*50 :: istring
	! For timing with MPI_Wtime:
	DOUBLE PRECISION :: time0, time1, time2, time3, time4, time_build, time_invert

  ! Stitching version
  ! Basic Plan
  ! will have N processors and N' identical regions
  ! Each processor gets some number of regions to build the matrix for
  ! All processors invert each region in parallell after it is distributed
  ! 

  ! For the duration of the program each processor holds :
  !	ata_small(my_nrow,my_ncol,N')  - its piece of each subregion
  !     aug_small(my_nrow,my_ncol,N')  - its piece of the augmenting matrix to ata that describes the stitching regularization
  !     regm                           - the full regularization matrix, which is the SAME for each IDENTICAL subregion

	! I want some MPI calls later, so init MPI
	CALL MPI_Init(ierr)
	CALL MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
	CALL MPI_Comm_size(MPI_COMM_WORLD, numproc, ierr)
	time0 = MPI_Wtime(ierr) ! start time

	! We also will use BLACS, so init that
	CALL blacs_pinfo(iam, nprocs)
	CALL blacs_setup ( iam, nprocs )


	IF (myid == 0) WRITE(*,*) 'Starting...'

	! Data_management.F
	CALL Main_Initialization()      ! Load in the solution grid parameters and
				  ! initialize the primary matrices

	serial_inversion = 0

	! Decide nprow, npcol base on numproc
	! These are needed for BLACS
	npcol = INT(FLOOR(SQRT(numproc*1D0)))
	nprow = numproc/npcol
	! Find a better way to do this later
	IF (nprow*npcol /= numproc) THEN
		WRITE(*,*) "ERROR: nprow * npcol /= numproc, pick a different numproc"
		STOP
	ENDIF
  
	CALL blacs_get ( -1, 0, context )
	CALL blacs_gridinit ( context, 'Row-major', nprow, npcol )
	CALL blacs_gridinfo ( context, nprow, npcol, myrow, mycol )
!  write(6,*)'Grid initialization complete for processor :', iam, n_regions
  


	CALL Get_Processor_Info()
	CALL Main_Allocation()
	IF (myid == 0) WRITE(*,*) 'Computing Regularization Matrices....'

	CALL Compute_Regularization_Matrices()
	time1 = MPI_Wtime(ierr) ! end of init
	IF (myid == 0) WRITE(*,*) 'Initialization Complete... Entering Main Loop'
	CALL Findmns()
	DO i=1, n_regions
		r_id = i-1
    	
		time2 = MPI_Wtime(ierr)
		Call Build_Region(r_id)
		time_build = time_build + MPI_Wtime(ierr) - time2

		rhsxy(:,1) = rhs_x(:)
		rhsxy(:,2) = rhs_y(:)

		DO ii = 1, my_nrow
			rhsxyzero(my_rows(ii),1,r_id+1) = rhsxy(ii,1)
			rhsxyzero(my_rows(ii),2,r_id+1) = rhsxy(ii,2)
		ENDDO

		IF (n_regions > 1) THEN
			CALL Add_Stitch(r_id)
		ENDIF
       
		! Everyone has built their piece of these matrices.  Wait on everyone to come together.
		CALL Blacs_Barrier(context,'A') ! Hold Everyone Up
		time2 = MPI_Wtime(ierr)
	    ! next, we have 4 matrices per region to invert
	    ! inverse matrices will be saved
		inv_dir = 1
		CALL pinv2_parallel(0,r_id) ! x-direction, zero-order inversion
		IF (n_regions > 1) THEN
			CALL pinv2_parallel(1,r_id) ! x-direction, stitching matrix inversion
		ENDIF
	    inv_dir = 2
        CALL pinv2_parallel(0,r_id) ! y-direction, zero-order inversion
		IF (n_regions > 1) THEN
			CALL pinv2_parallel(1,r_id) ! y-direction, stitching matrix inversion
		ENDIF
	    time_invert = time_invert+(MPI_Wtime(ierr)-time2)
	ENDDO
	CALL blacs_barrier(context,'A')

  If (n_regions .gt. 1) Then
	  DeAllocate(ataxm)
	  DeAllocate(ataym)
	  If (have_diagonal) Then
		  DeAllocate(my_diagonal)
	  Endif
  Endif

  DeAllocate(nkindex)
  

	time3 = MPI_Wtime(ierr)
  ! Now iterate for max number of iterations
  if (n_regions .eq. 1) then
	max_iter = 0
  endif 
  Do iter = 0, max_iter
	esqx_last(:,:) = esqx(:,:)
	esqy_last(:,:) = esqy(:,:)
	lastvx(:,:) = vx(:,:)
	lastvy(:,:) = vy(:,:)
	if (iam .eq. 0) then
		write(6,*)'iter = ', iter
	endif

	Do r_id = 1, n_regions
	  	Call Update_RHS(iter,r_id)  ! this will do nothing when iter is zero
		inv_dir = 1
	  	Call RHS_Solve_Parallel(iter,r_id)
		inv_dir = 2
	  	Call RHS_Solve_Parallel(iter,r_id)
	Enddo
	If (iter .gt. 0) then
		If (propagate_errors == .True.) then
			call advance_errors()
		Endif
	endif

	if (iam .eq. 0) then  !output the results for this iteration
		!open the file
		! write nsolution
		! write n_regions
		WRITE(istring,'(i10)')iter
		savefile = 'results_'//TRIM(ADJUSTL(istring))//'.dat'
		open(50,FILE=savefile,status='replace',form='formatted',&
                  & access='sequential')
		write(50,180)nsolution
		write(50,180)n_regions
		do r_id = 1, n_regions
			del_lat = lat_lims(1,r_id)-lat_lims(1,1)
			del_lon = lon_lims(1,r_id)-lon_lims(1,1)
			do i = 1, nsolution
				write(50,179)(lons_pack(i)+del_lon),(lats_pack(i)+del_lat),zs_pack(i),vx(i,r_id),sqrt(esqx(i,r_id)),vy(i,r_id),sqrt(esqy(i,r_id))
			enddo
		enddo
		close(50)
180     format(i10)
179     format(7E19.12)
	endif
  Enddo

	time4 = MPI_Wtime(ierr) ! end of everything
	IF (myid == 0) THEN
		WRITE(*,'(A)') '------------------ Run Timing -------------------------'
		WRITE(*,'(A,F8.2,A)') ' Total time elapsed = ', (time4-time0)/60D0, ' min'
		WRITE(*,'(A,F8.2,A)') '          Init Time = ', (time1-time0)/60D0, ' min'
		WRITE(*,'(A,F8.2,A)') '         Build Time = ', time_build/60D0, ' min'
		WRITE(*,'(A,F8.2,A)') '     Inversion Time = ', time_invert/60D0, ' min'
		WRITE(*,'(A,F8.2,A)') '     Stitching Time = ', (time4-time3)/60D0, ' min'
	ENDIF

  if (serial_inversion == 0) then
	  call blacs_exit (0)
  endif
!	CALL MPI_Finalize(ierr) ! this segfaults. what.
END PROGRAM Main
