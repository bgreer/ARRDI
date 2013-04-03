Program Main
  Use Data_Management
  Use Run_Info
  Use Fitslib_NF
  Use Matrix_Magic
  Implicit None
  Integer :: i,cct, j, ii, jj
  Integer*8 :: n1,n2,n3,n4
  Integer*8 :: nnw, regions_per_pe, n_rounds, last_id, r_id, iter_reg
  integer :: nrow, ncol, cxt, serial_inversion, write_time, rs, cs, last_region, left_over
  integer :: time0, time1, time2, count_rate, time3, iter, time_invert, time_build, time_init, time4, time_solve
  Integer*8 :: fits_test4(2,2,2,2)
  Real*8  :: elapsed_time, init_time, build_time, inversion_time, del_lat, del_lon, solve_time
  Real*8, Allocatable :: reg_check(:), test_arr(:,:)
  Real*8 :: e_save(1:2,1:2)
  Character*120 :: check_file, savefile, filename
  Character*50 :: istring

!  call test_sort()
!  stop
  Call System_Clock(time0,count_rate)
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



  Write(6,*)'Starting....'

  Call Main_Initialization()      ! Load in the solution grid parameters and
				  ! initialize the primary matrices

  
  

  write(6,*)'calling gridinit'
  serial_inversion = 0
  
  call blacs_pinfo(iam, nprocs)
  call blacs_setup ( iam, nprocs )
  call blacs_get ( -1, 0, context )
  call blacs_gridinit ( context, 'Row-major', nprow, npcol )
  call blacs_gridinfo ( context, nprow, npcol, myrow, mycol )
!  write(6,*)'Grid initialization complete for processor :', iam, n_regions
  


  Call Get_Processor_Info()
  Call Main_Allocation()
  Write(6,*)'Computing Regularization Matrices....'

  Call Compute_Regularization_Matrices()
  Call System_Clock(time1,count_rate)
  time_init = time1-time0
  Write(6,*)'Initialization Complete... Entering Main Loop'
  Call Findmns()
  Do i = 1, n_regions
    r_id = i-1
    
       Call Build_Region(r_id)
	rhsxy(:,1) = rhs_x(:)
	rhsxy(:,2) = rhs_y(:)

	Do ii = 1, my_nrow
	  rhsxyzero(my_rows(ii),1,r_id+1) = rhsxy(ii,1)
	  rhsxyzero(my_rows(ii),2,r_id+1) = rhsxy(ii,2)
	Enddo

       If (n_regions .gt. 1) Then
          Call Add_Stitch(r_id)
       Endif
       
    ! Everyone has built their piece of these matrices.  Wait on everyone to come together.
    Call Blacs_Barrier(context,'A') ! Hold Everyone Up
  
	    Call System_Clock(time2,count_rate)
	    ! next, we have 4 matrices per region to invert
	    ! inverse matrices will be saved
	     inv_dir = 1
    	    call pinv2_parallel(0,r_id) ! x-direction, zero-order inversion
	    If (n_regions .gt. 1) Then
	            call pinv2_parallel(1,r_id) ! x-direction, stitching matrix inversion
	    Endif
	     inv_dir = 2
            call pinv2_parallel(0,r_id) ! y-direction, zero-order inversion
	    If (n_regions .gt. 1) Then
	            call pinv2_parallel(1,r_id) ! y-direction, stitching matrix inversion
	    Endif
	    Call System_Clock(time3,count_rate)
	    time_invert = time_invert+(time3-time2)
  Enddo
  Call System_Clock(time4,count_rate)
  time_build = time4-time1-time_invert
  call blacs_barrier(context,'A')

  If (n_regions .gt. 1) Then
	  DeAllocate(ataxm)
	  DeAllocate(ataym)
	  If (have_diagonal) Then
		  DeAllocate(my_diagonal)
	  Endif
  Endif

  DeAllocate(nkindex)
  


  Call System_Clock(time2,count_rate)

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
 Call System_Clock(time3,count_rate)
 time_solve = time3-time2
  write_time = 1
  if (serial_inversion == 0) then
	if (iam .ne. 0) then
		write_time = 0
	endif
  endif

  if (write_time .eq. 1) then
	  Call System_Clock(time3,count_rate)

	  elapsed_time = Real((time3-time0))/Real(count_rate)
	  inversion_time = Real((time_invert))/Real(count_rate)
	  build_time = Real((time_build))/Real(count_rate)
	  init_time = Real((time_init))/Real(count_rate)
	  solve_time = Real((time_solve))/Real(count_rate)

    Write(6,*) '------------------ Full Run Timings -------------------------'
    Write(6,*) 'Total time elapsed (min) = ',elapsed_time /60.d0
    Write(6,*) '         Init time (min) = ',init_time / 60.d0
    Write(6,*) '         Build Time (min) = ',build_time / 60.d0
    Write(6,*) '         Inversion Time (min) = ',inversion_time / 60.d0
    Write(6,*) '         Solve Time (min) = ',solve_time / 60.d0

  endif
  if (serial_inversion == 0) then
	  call blacs_exit (0)
  endif
End Program Main
