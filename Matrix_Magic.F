Module Matrix_Magic
  Use Run_Info
  Use Fitslib_NF
  Use Solution_Grid
  Use Data_Management
  Use Region_Info
  Use Kernel_Management
  Implicit None
    integer :: n_blocks_x, n_blocks_y,nb,mb,nba, n_blocks_xa
    integer*8 :: my_num_sol_pts, nmeasurements, next_big_column, next_local_column
    integer :: x_block_per_pe, y_block_per_pe, x_block_per_pea
    integer :: my_nblocks_x, my_nblocks_y, my_ncol, my_nrow, y_block_mod, x_block_mod
    integer :: my_nblocks_xa, my_ncola, x_block_moda
    integer, Allocatable :: block_size_x(:), block_size_y(:), block_size_xa(:), have_column(:), isboundary(:)
    integer, Allocatable :: my_blocks_x(:), my_blocks_xa(:), my_blocks_y(:), my_rows(:), my_cols(:), my_colsa(:)
	integer, Allocatable :: rowmap(:), colmap(:)
    Integer, Allocatable :: true_overlap(:,:), found_partner(:,:)
    integer :: have_diagonal
    integer*8 :: my_n_diagonal
    integer*8, Allocatable :: my_diagonal(:,:), xncount(:), yncount(:), x_neighbors(:,:), y_neighbors(:,:)
    Real*8 :: lh_amp
    Real*8, Allocatable :: stitch_weight(:,:)
    Real*8, Allocatable :: ainv0(:,:,:,:), ainvn(:,:,:,:), covar0(:,:,:,:), rhsxysave(:,:,:), vx(:,:), vy(:,:)
    Real*8, Allocatable :: lastvx(:,:), lastvy(:,:)
    Real*8, Allocatable :: atransposex(:,:), atransposey(:,:), reskercoefs(:,:), all_diagonals(:,:)
    Real*8, Allocatable :: esqx(:,:), esqy(:,:), xerh(:), yerh(:), esqx_last(:,:), esqy_last(:,:)
    Real*8, Allocatable :: rows_and_cols_unsorted(:), rows_and_cols(:), esqx_iter(:,:), esqy_iter(:,:)
	Real*8, Allocatable :: boundary_factor(:)
Contains
  Subroutine Main_Allocation()
     Implicit None
     ! Allocate the primary arrays needed for the code to run
     erhx = 0.0     ! hrms error
     erhy = 0.0     ! hrms error

      Allocate(atax(1:my_nrow,1:my_ncol))
      Allocate(atay(1:my_nrow,1:my_ncol))
!     Allocate(atax(1:nsolution,1:nsolution))  ! Matrix to invert - x
!     Allocate(atay(1:nsolution,1:nsolution))  ! ' ' -y
     If (n_regions .gt. 1) Then 
!	     Allocate(ataxm(1:nsolution,1:nsolution))  ! Stitching Matrix to invert - x
!	     Allocate(ataym(1:nsolution,1:nsolution))  ! ' ' -y
	     Allocate(ataxm(1:my_nrow,1:my_ncol))
	     Allocate(ataym(1:my_nrow,1:my_ncol))	
     Endif
	
     ! temporary matrices to hold copies of the above
!     Allocate(ataxt(1:nsolution,1:nsolution))  ! Matrix to invert - x
!     Allocate(atayt(1:nsolution,1:nsolution))  ! ' ' -y
!     If (n_regions .gt. 1) Then
!	     Allocate(ataxmt(1:nsolution,1:nsolution))  ! Stitching Matrix to invert - x
!	     Allocate(ataymt(1:nsolution,1:nsolution))  ! ' ' -y
!     Endif
     Allocate(rhsxy(1:my_nrow,2))
!     Allocate(rhsxyt(1:nsolution,2))
     
     Allocate(regm(1:my_nrow,1:my_ncol))  ! regularization matrix

     Allocate(rhs_x(1:my_nrow))	      ! right hand side of matrix equation - x
     Allocate(rhs_y(1:my_nrow))	      ! ' ' - y



     atax(:,:) = 0.0D0
     atay(:,:) = 0.0D0
     rhs_x(:) = 0.0D0
     rhs_y(:) = 0.0D0


  End Subroutine Main_Allocation
  Subroutine Get_Processor_Info
	Implicit None
	Integer :: x_mod, y_mod, thiscount
	Integer :: ii, jj,i,j,nsol,m,n
	!Here we figure out which piece of the global matrices our processor owns

	nb = x_block_size	
        mb = y_block_size  
	nsol = nsolution
        n = nsol
        m = nsol
  
        if ( nprocs < 1 ) then
          if ( iam == 0 ) then
            nprocs = nprow * npcol
          end if
        end if


	call blacs_gridinfo ( context, nprow, npcol, myrow, mycol )


	if ( myrow == -1 ) then
	  call blacs_exit ( 0 )
          stop
	end if


	! first, we identify the number of blocks that our grid
	! is partitioned into base off of the block size
 	! nb is the x block size and mb is the y block size (i.e. nb is column block size and mb is row block size)

	n_blocks_x = nsol/nb
	x_mod = MOD(nsol,nb)

	if (x_mod .gt. 0) then
	  n_blocks_x = n_blocks_x+1
	endif

	n_blocks_y = nsol/nb
	y_mod = MOD(nsol,mb)
	if (y_mod .gt. 0) Then
	  n_blocks_y = n_blocks_y+1
	endif

	Allocate(block_size_x(1:n_blocks_x))
	Allocate(block_size_y(1:n_blocks_y))


	block_size_x(:) = nb
	block_size_y(:) = mb


	If (x_mod .gt. 0) Then    ! take care of the remainder
   	  block_size_x(n_blocks_x) = x_mod
	Endif

	If (y_mod .gt. 0) Then
 	  block_size_y(n_blocks_y) = y_mod
	Endif



	! Now we need to know who has which blocks
 
	! First, figure out how many blocks each processor has
	x_block_per_pe = n_blocks_x/npcol
	y_block_per_pe = n_blocks_y/nprow

	x_block_mod = MOD(n_blocks_x,npcol)
	y_block_mod = MOD(n_blocks_y,nprow)

 
	my_nblocks_x = x_block_per_pe  ! everyone gets at least this many blocks
	my_nblocks_y = y_block_per_pe  ! note that it could be zero in theory.
 
	If (x_block_mod .gt. 0) Then
	  If (mycol .lt. x_block_mod) Then ! Not .le. because mycol starts at 0
	    my_nblocks_x = my_nblocks_x+1
	  Endif
	Endif

	If (y_block_mod .gt. 0) Then
	  If (myrow .lt. y_block_mod) Then ! Not .le. because myrow starts at 0
	    my_nblocks_y = my_nblocks_y+1
	  Endif
	Endif



	! Now, the blocks are delt out like a deck of cards as noted in
	! the scalapack tutorials.   In doing so, we'll figure out how many
	! rows and columns are in our piece of the distributed matrix BigA

	Allocate(my_blocks_x(1:my_nblocks_x))
	Allocate(my_blocks_y(1:my_nblocks_y))

	my_nrow = 0
	my_ncol = 0
	do i = 1, my_nblocks_x
	  my_blocks_x(i) = mycol+1 +(i-1)*npcol ! plus 1 since mycol starts at 0
  					  ! but block inds start at 1
	  my_ncol = my_ncol+block_size_x(my_blocks_x(i))
	enddo

	do i = 1, my_nblocks_y
	  my_blocks_y(i) = myrow +1 +(i-1)*nprow ! plus 1 since myrow starts at 0
 					  ! but block inds start at 1
	  my_nrow = my_nrow+block_size_y(my_blocks_y(i))
	enddo


!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!	DDDDDDDDDDDDDDDD
	Allocate(my_cols(1:my_ncol))
	Allocate(my_rows(1:my_nrow))
	Allocate(colmap(1:nsolution))
	Allocate(rowmap(1:nsolution))
	rowmap(:) = -1
	colmap(:) = -1
!	my_num_sol_pts = my_ncol+my_nrow
!	Allocate(rows_and_cols_unsorted(1:my_num_sol_pts))
!	Allocate(rows_and_cols(1:my_num_sol_pts))
	thiscount =1
	do j = 1, my_nblocks_y
	  do i = 1, block_size_y(my_blocks_y(j))
		my_rows(thiscount) = (my_blocks_y(j)-1)*mb+i
		! rowmap will map a global index to a local index
		rowmap((my_blocks_y(j)-1)*mb+i) = thiscount
		thiscount = thiscount+1
	  enddo
	enddo
!	Write(6,*)'check rows : ', my_nrow, (thiscount-1)
	thiscount =1

	do j = 1, my_nblocks_x
	  do i = 1, block_size_x(my_blocks_x(j))
		my_cols(thiscount) = (my_blocks_x(j)-1)*nb+i
		colmap((my_blocks_x(j)-1)*nb+i) = thiscount
		thiscount = thiscount+1
	  enddo
	enddo
	If (n_regions .gt. 1) Then
		If (Allocated(stitch_weight)) DeAllocate(stitch_weight)
		Allocate(stitch_weight(2,1:nsolution))
		stitch_weight(:,:) = 0.0D0
		call find_my_diagonal()
	Endif
!	Write(6,*)'check cols : ', my_ncol, (thiscount-1)


	!  Now that we know the dimensions of our piece of the big matrix,
	!  We can allocate the smaller matrices that we actually do the work on.


 
	Allocate(ainv0(1:my_nrow,1:my_ncol,n_regions,1:2))
	Allocate(covar0(1:my_nrow,1:my_ncol, n_regions,1:2))
	If (n_regions .gt. 1) Then
		Allocate(ainvn(1:my_nrow,1:my_ncol,n_regions,1:2))
	Endif
	Allocate(rhsxyiter(1:nsolution,2,1:n_regions))
	Allocate(rhsxyzero(1:nsolution,2,1:n_regions))
	Allocate(vx(1:nsolution,n_regions))
	Allocate(vy(1:nsolution,n_regions))
	Allocate(lastvx(1:nsolution,n_regions))
	Allocate(lastvy(1:nsolution,n_regions))
	Allocate(esqy(1:nsolution,n_regions))
	Allocate(esqx(1:nsolution,n_regions))
	Allocate(esqx_iter(1:nsolution,n_regions))
	Allocate(esqy_iter(1:nsolution,n_regions))
	esqx(:,:) = 0.0D0
	esqy(:,:) = 0.0D0
	Allocate(esqx_last(1:nsolution,n_regions))
	Allocate(esqy_last(1:nsolution,n_regions))
	Allocate(xerh(1:n_regions))
	Allocate(yerh(1:n_regions))

  End Subroutine Get_Processor_Info
  Subroutine Find_My_Diagonal()
	Integer*8 :: i, j, dd
	Integer*8, Allocatable :: temp(:,:)
	! Using my_rows and my_cols, identify which
	! elements of the diagonal of the global matrix 
	! this processor possesses.  Store those elements in
	! the matrix my_diagonal([dd,i,j]) such that
	! my_cols(i) = dd and my_rows(j) = dd
	! I think i should equal j due to the blocking, but just to be safe
	!
	! Variables set in this subroutine are my_diagonal, my_n_diagonal, have_diagonal
	Allocate(temp(1:nsolution,1:3))
	my_n_diagonal = 0
	have_diagonal = 0
	do j = 1, my_ncol
	  do i = 1, my_nrow
		if (my_rows(i) .eq. my_cols(j)) then
			dd = my_rows(i)
			my_n_diagonal = my_n_diagonal+1
			temp(my_n_diagonal,1) = dd
			temp(my_n_diagonal,2) = i
			temp(my_n_diagonal,3) = j
		endif
	  enddo
	enddo
	if (my_n_diagonal .gt. 0) then
		have_diagonal = 1
		Allocate(my_diagonal(1:my_n_diagonal,1:3))
		my_diagonal(1:my_n_diagonal,1:3) = temp(1:my_n_diagonal,1:3)
	endif
	DeAllocate(temp)
  End Subroutine Find_My_Diagonal


  Subroutine Get_Atranspose_Info()
	Implicit None
	Integer :: x_mod, y_mod, thiscount
	Integer :: ii, jj,i,j,nsol,m,n


	nba = x_block_size_a
	nsol = nsolution
        n = nsol
        m = nsol
  
        if ( nprocs < 1 ) then
          if ( iam == 0 ) then
            nprocs = nprow * npcol
          end if
        end if


	call blacs_gridinfo ( context, nprow, npcol, myrow, mycol )

!	Write(6,*)'made it past gridinfo'
	if ( myrow == -1 ) then
	  call blacs_exit ( 0 )
          stop
	end if


	! first, we identify the number of blocks that our grid
	! is partitioned into base off of the block size
 	! nb is the x block size and mb is the y block size (i.e. nb is column block size and mb is row block size)

	n_blocks_xa = nmeasurements/nba
	x_mod = MOD(nmeasurements,nba)

	if (x_mod .gt. 0) then
	  n_blocks_xa = n_blocks_xa+1
	endif


	If (Allocated(block_size_xa)) DeAllocate(block_size_xa)
	Allocate(block_size_xa(1:n_blocks_xa))
!	Write(6,*)'past block size xa alloaction'
	block_size_xa(:) = nba
	If (x_mod .gt. 0) Then    ! take care of the remainder
   	  block_size_xa(n_blocks_xa) = x_mod
	Endif


	! Now we need to know who has which blocks
 
	! First, figure out how many blocks each processor has
	x_block_per_pea = n_blocks_xa/npcol

	x_block_moda = MOD(n_blocks_xa,npcol)

 
	my_nblocks_xa = x_block_per_pea  ! everyone gets at least this many blocks
	                                 ! note that it could be zero in theory.
 
	If (x_block_moda .gt. 0) Then
	  If (mycol .lt. x_block_moda) Then ! Not .le. because mycol starts at 0
	    my_nblocks_xa = my_nblocks_xa+1
	  Endif
	Endif


	! Now, the blocks are delt out like a deck of cards as noted in
	! the scalapack tutorials.   In doing so, we'll figure out how many
	! rows and columns are in our piece of the distributed matrix BigA

	If (Allocated(my_blocks_xa)) DeAllocate(my_blocks_xa)
	Allocate(my_blocks_xa(1:my_nblocks_xa))
!	Write(6,*)'past my blocks xa'


	my_ncola = 0
	do i = 1, my_nblocks_xa
	  my_blocks_xa(i) = mycol+1 +(i-1)*npcol ! plus 1 since mycol starts at 0
  					  ! but block inds start at 1
	  my_ncola = my_ncola+block_size_xa(my_blocks_xa(i))
	enddo


	If (Allocated(my_colsa)) DeAllocate(my_colsa)
	Allocate(my_colsa(1:my_ncola))
!	write(6,*)'past my colsa'

	thiscount =1

	do j = 1, my_nblocks_xa
	  do i = 1, block_size_xa(my_blocks_xa(j))
		my_colsa(thiscount) = (my_blocks_xa(j)-1)*nba+i
		thiscount = thiscount+1
	  enddo
	enddo

!	Write(6,*)'check cols : ', my_ncola, (thiscount-1)

!	If (Allocated(have_column)) DeAllocate(have_column)
!	Allocate(have_column(1:nmeasurements))
!	have_column(:) = 0
!	do j = 1, my_ncola
!		have_column(my_colsa(j)) = 1
!	enddo


 	If (Allocated(atransposex)) DeAllocate(atransposex)
 	If (Allocated(atransposey)) DeAllocate(atransposey)
	If (Allocated(reskercoefs)) DeAllocate(reskercoefs)
!	write(6,*)'atranspose in'
	Allocate(atransposex(1:my_nrow,1:my_ncola))
	Allocate(atransposey(1:my_nrow,1:my_ncola))
	Allocate(reskercoefs(1:my_nrow,1:my_ncola))
	! write(6,*)'atranspose out'
	atransposex(:,:) = 0.0D0
	atransposey(:,:) = 0.0D0
!	write(6,*)'atranspose out', myrow,mycol
    Call Blacs_Barrier(context,'A') ! Hold Everyone Up
  End Subroutine Get_Atranspose_Info


	SUBROUTINE Build_Region(region_id)
		IMPLICIT NONE
		INTEGER*8 :: region_id, i, cct, tsize_pix, delta_pix, ndtemp
		REAL*8 :: tsize_deg, apod_factor
		CHARACTER*50 :: rstring, tstring
		CHARACTER*120 :: parfile, dfile
		INTEGER :: ierr
		! for now, everyone reads in the same files
		! these filenames should depend on the region id, however
		WRITE(rstring,'(i10)')region_id

		cct = 0
		atax(:,:) = 0.0
		rhs_x(:) = 0.0
		atay(:,:)= 0.0
		rhs_y(:) = 0.0
		erhx = 0.0
		erhxc = 0.0
		erhy = 0.0
		erhyc = 0.0
		! Here we go ahead and get the total number of measurements
		!	If(Allocated(nmeasurements)) Then DeAllocate(nmeasurements)
		!	Allocate(nmeasurements
		nmeasurements = 0
		DO i = 1, ntsizes
			nmeasurements = nmeasurements + ndata(i)
		ENDDO
		!	write(6,*)'nmeasurements', nmeasurements
		CALL Get_Atranspose_Info()
		next_big_column = 1
		next_local_column = 1
		DO i = 1, ntsizes
			tsize_pix = tsizes(i)
			tsize_deg = tsize_pix/8.0
			apod_factor = 7.5/8.0
			delta_pix = (nxk-int(tsize_deg*apod_factor/dlonlat)-1)/2  ! split the difference
			i_min = 1+delta_pix
			i_max = nxk-delta_pix
			j_min = i_min
			j_max = i_max
			! Loading kernels and data does not take much time
			! Maybe a few seconds, but this could be significant for tiny
			! inversions
			CALL Load_Kernels(tsize_pix) ! TODO: make this parallel
			smeared_kers = smeared_kers*dlonlat*dlonlat
			CALL Load_Dataset(region_id,i)
			CALL Create_Nonzero_Table(i)  ! more bookeeping for parallel build of matrix
			btime0 = MPI_Wtime(ierr)
			CALL Augment_A_Parallel(i)   ! load the A_transpose#A and A_tranpose#data for each tilesize
			btime1 = btime1 + MPI_Wtime(ierr) - btime0
			! ^ takes the majority of the time

!			DEALLOCATE(uxy)
!			DEALLOCATE(n)
!			DEALLOCATE(ki)
			DEALLOCATE(kxy)
			DEALLOCATE(indic)
			DEALLOCATE(indicn)
			DEALLOCATE(nktindex)
			DEALLOCATE(ata_indices)
			DEALLOCATE(ata_num)
		ENDDO
		erhx = erhx/erhxc
		erhx = sqrt(1.0/erhx)
		erhy = erhy/erhyc
		erhy = sqrt(1.0/erhy)

		IF (myid == 0) WRITE(*,'(A,F8.3,A,F8.3)')' Error Harmonic Mean is :', erhx, ', ', erhy
		atax = atax*(erhx*erhx) ! +regm    - regularization is now done in pinv()
		atay = atay*(erhy*erhy) ! +regm
		atransposex = atransposex*(erhx*erhx)
		atransposey = atransposey*(erhy*erhy)
		rhs_x = rhs_x*(erhx*erhx)
		rhs_y = rhs_y*(erhy*erhy)	
	END SUBROUTINE Build_Region

	SUBROUTINE Create_Nonzero_Table(ts_index)
		! no implicit none?
		INTEGER*8 :: nloc, ilon, ilat, ii, jj, ind, ts_index
		INTEGER*8, ALLOCATABLE :: nz_locs(:)
		ALLOCATE(nz_locs(1:nsolution))		

		DO ilat = 1, nlat(ts_index)
			DO ilon = 1, nlon(ts_index)
				nloc = ata_num(ilon,ilat)
				nz_locs(:) = 0
				DO ii = 1, nloc
					nz_locs(ata_indices(ii,ilon,ilat)+1) = 1
				ENDDO
				ind = 0
				DO jj = 1, my_ncol
					IF (nz_locs(my_cols(jj)) == 1) THEN
						ind = ind+1
						jlocal_indices(ind,ilon,ilat) = jj
						nzc_indices(ind,ilon,ilat) = my_cols(jj)
					ENDIF
				ENDDO
				nznc_num(ilon,ilat) = ind
				ind = 0
				DO ii = 1, my_nrow
					IF (nz_locs(my_rows(ii)) == 1) THEN
					ind = ind+1
					ilocal_indices(ind,ilon,ilat) = ii
					nzr_indices(ind,ilon,ilat) = my_rows(ii)
					ENDIF
				ENDDO
				nznr_num(ilon,ilat) = ind
			ENDDO
		ENDDO

		DEALLOCATE(nz_locs)
	END SUBROUTINE Create_Nonzero_Table

	SUBROUTINE Augment_A_Parallel(ts_index)
    ! Different from the original Augment A in that each processors builds a piece of
    ! each matrix - rather than a whole matrix which is then communicates.  Should allow
    ! for more memory efficiency and should help with the computation of the averaging
    ! kernel coefficients.
		IMPLICIT NONE
		REAL*8, ALLOCATABLE :: augres(:,:,:)
		REAL*8, ALLOCATABLE :: result(:,:,:)
		REAL*8, ALLOCATABLE :: spread_ker(:), temp(:)
		REAL*8, ALLOCATABLE :: slice_new(:,:)
		REAL*8 :: dx,dy,ex,ey,esqx,esqy, my_total, prod

		INTEGER*8 :: ii, ik, el, icn, iz, im, jm, ilat, ilon, this_j, this_i, pack_ind, ts_index
		INTEGER*8 :: el_num, s_num, el_ind, nnk, this_index, aindex, iin, num, i, j, offset
		INTEGER*8 :: some_i, some_j
		INTEGER*8 :: nshftx, nshfty,loc_num, nstop, test, my_j, my_i, my_z, nx0, ny0, s_ind, this_ind
		INTEGER*8 :: nznc, nznr
		INTEGER*8, ALLOCATABLE :: spread_locs(:), nzrows(:), nzcols(:), ilocal(:), jlocal(:)
		CHARACTER*120 :: tempfile
		n_lats = nlat(ts_index)   ! this is a revision to Augment_A to 
		n_lons = nlon(ts_index)   ! make this routine mesh with the data_management routines.
		ALLOCATE(augres(1:nxk,1:nyk,1:nzs))  ! augres gets dimensions of kernels horizontally -
                                    ! not of full solution grid
		ALLOCATE(result(1:nxk,1:nyk,1:nzk))
		ALLOCATE(temp(1:nzk))

		nx0 = nxsz(1)
		ny0 = nysz(1)

		DO ii = 1, n_n(ts_index)
			IF (indicn(ii) > 0) THEN  ! do we process this n?
				DO ik = 1, n_k(ts_index)
					IF (indic(ik,ii) > 0) THEN !do we process this k?

						! first, find the kernel for this mode
						this_index = 0
						nnk = 1
						DO WHILE (this_index == 0)
							IF (nkindex(1,nnk) == n(ii,ts_index) ) THEN
								IF (nkindex(2,nnk) == ki(ik,ts_index)) THEN
									this_index = nnk
								ENDIF
							ENDIF
							nnk = nnk+1
							IF (nnk > nkers) THEN
								WRITE(*,*)'Fatal Error :  kernel not found for this n k combo : ', n(ii,ts_index), ki(ik,ts_index)
								STOP
							ENDIF
						ENDDO

						! now we have identified where in the kernel array we are
						! loop over solution grid depths and compact the kernel onto the solution
						! grid in depth
						result(:,:,:) = smeared_kers(:,:,:,this_index)

						DO  iz = 1, nzs
							temp(:) = z_owner(iz,:)*wz(:)
							DO im = i_min, i_max
								DO jm = j_min, j_max
									augres(im,jm,iz) = SUM(result(im,jm,:)*temp(:))
								ENDDO    ! jm
							ENDDO      ! im
						ENDDO   ! iz
						! The basic kernel for one n+k combination is stored in augres
						! at this point, augres is on the proper z solution grid,
						! but is still uniform in x and y.

						! We now read this kernel into the inversion matrix and the rhs
						! at every location that a measurement of this n+k was made.
						DO ilat = 1, n_lats
							DO ilon = 1, n_lons
								loc_num = ata_num(ilon,ilat)
								ALLOCATE(spread_locs(1:loc_num))
								ALLOCATE(spread_ker(1:nsolution))
								spread_locs(:) = ata_indices(1:loc_num,ilon,ilat)  ! non zero locations
                                                                   ! of the kernel
								! We are very concerned about where the nonzero locations are because
								! the kernel matrix A is sparse.  The bulk of this loop is essentially
								! bookeeping for a sparse matrix multiply.

								!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
								!  ---Parallell Modification ----
								nznr = nznr_num(ilon,ilat)  ! number of rows in this processor that correspond to nonzero kernel locations
								nznc = nznc_num(ilon,ilat)  ! number of columns
								ALLOCATE(nzrows(1:nznr))    ! The nonzero rows  (indexed with respect to the big matrix - i.e. 
								ALLOCATE(nzcols(1:nznc))    ! The nonzero columns (these rows and columns can run from 1 to nsolution)

								ALLOCATE(ilocal(1:nznr))    ! Indices in the local, small array that correspond to the nonzero
								ALLOCATE(jlocal(1:nznc))    ! rows and columns above

								nzrows(:) = nzr_indices(1:nznr,ilon,ilat)
								nzcols(:) = nzc_indices(1:nznc,ilon,ilat)
								ilocal(:) = ilocal_indices(1:nznr,ilon,ilat)
								jlocal(:) = jlocal_indices(1:nznc,ilon,ilat)
								!PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

								aindex = nktindex(ii,ik,ilat,ilon)+1 !idl stuff starts at 0
								IF (aindex >= 1) THEN ! we measured this mode at this locale
									nshftx = kxy(1,ilon,ilat)  ! Shift kernel into proper place
									nshfty = kxy(2,ilon,ilat)  ! on the uniform grid

									! just do it here - not modular, but it's pretty unique code anyway
									! this could be optimized further, but we'll stick to the idl format

									! we need an approach suitable for fortran
									! that ensure all of a kernel is distributed among its owners
									! Let us build a reverse lookup table instead and iterate over the kernel points
									! 16 deg : i_min = 1, i_max = 17   - all points
									! 4 deg  : i_min = 7, i_max = 11   - 5 square
									! 2 deg  : i_min = 8, i_max = 10   - 3 square
									! indarr says : given a location on the regular grid, what indices on the nonuniform grid "own" me.
									! weight_arr says : these indices each own me by this amount
									spread_ker(:) = 0.0
									offset = 0
									DO iz = 1, nzs
										DO my_j = j_min, j_max
											DO my_i = i_min, i_max
												num = num_arr(my_i+nshftx,my_j+nshfty,iz)
												my_total = 0.0
												DO iin = 1, num
													this_ind = ind_arr(iin,my_i+nshftx,my_j+nshfty,iz)+1  ! index on 2d slice (with nxsz0 and nysz0 points)
													this_j = (this_ind/nx0)   ! +1 left out because I would subtract it below anyway
													this_i = MOD(this_ind,nx0)
													s_ind = offset+this_i+(this_j)*nxsz(iz)
													spread_ker(s_ind) = spread_ker(s_ind)+weight_arr(iin,my_i+nshftx,my_j+nshfty,iz)*augres(my_i,my_j,iz)
												ENDDO
											ENDDO
										ENDDO
										offset = offset+nxsz(iz)*nysz(iz)
									ENDDO
									!########   We have just computed spread_ker.  This is the 3-D kernel, represented in 1-D
									!#######    on the inversion grid (it has been converted from the native
									!#######    kernel grid using weight_arr
									dx = uxy(1,aindex,ts_index)
									ex = uxy(2,aindex,ts_index)
									dy = uxy(3,aindex,ts_index)
									ey = uxy(4,aindex,ts_index)
									esqx = ex*ex
									esqy = ey*ey

									erhx = erhx+1.0/esqx
									erhy = erhy+1.0/esqy
									erhxc = erhxc+1.0
									erhyc = erhyc+1.0

									! Parallel Modification.  Because we are no longer doing the entire
									! matrix, the manner in which we load our matrices changes
									DO my_i = 1, nznr
										some_i = ilocal(my_i)
										rhs_x(some_i) = rhs_x(some_i)+dx/(esqx)*spread_ker(nzrows(my_i))
										rhs_y(some_i) = rhs_y(some_i)+dy/(esqy)*spread_ker(nzrows(my_i))
									ENDDO

									DO my_j = 1, nznc
										DO my_i = 1, nznr
											some_i = ilocal(my_i)
											some_j = jlocal(my_j)
											prod = spread_ker(nzrows(my_i))*spread_ker(nzcols(my_j))
											atax(some_i,some_j) = atax(some_i,some_j)+prod/esqx
											atay(some_i,some_j) = atay(some_i,some_j)+prod/esqy
										ENDDO
									ENDDO
									IF (next_local_column <= my_ncola) THEN
										IF (my_colsa(next_local_column) == next_big_column) THEN
											DO my_i = 1, nznr		! into that matrix as well
												some_i = ilocal(my_i)
												atransposex(some_i,next_local_column) = spread_ker(nzrows(my_i))/esqx
												atransposey(some_i,next_local_column) = spread_ker(nzrows(my_i))/esqy	
											ENDDO
											next_local_column = next_local_column+1
										ENDIF
									ENDIF
									next_big_column = next_big_column+1
									!PPPPPPPPPPPPPPPPPPPPPPPPPPP
								ENDIF
								DEALLOCATE(nzrows)
								DEALLOCATE(nzcols)
								DEALLOCATE(ilocal)
								DEALLOCATE(jlocal)
								DEALLOCATE(spread_ker)
								DEALLOCATE(spread_locs)
							ENDDO  ! ilon
						ENDDO ! ilat
					ENDIF
				ENDDO  ! ik
			ENDIF
		ENDDO    ! ii
		DEALLOCATE(result)
		DEALLOCATE(augres)
		DEALLOCATE(temp)
	END SUBROUTINE Augment_A_Parallel


	Subroutine Add_Stitch(r_id)
		Implicit None
		Integer*8 :: r_id, dd, indi, indj
		Integer :: i, j, k,kk, ncount, cc
		real*8 :: max_x, max_y, test_x, test_y, stitch_min, sweight_factor
		real*8 :: sigma2
		stitch_min = 0.0d0
		ataxm(:,:) = atax
		ataym(:,:) = atay
		If (.not. Allocated(boundary_factor)) Then
			Allocate(boundary_factor(1:nsolution))
			boundary_factor(:) = 1.0d0
			If (onesided_boundary_stitching) Then
				indi = 1
				Do k = 1, nzs
					Do j = 1, nysz(k)
						Do i = 1, nxsz(k)
							! wrap sigma into boundary_factor
							boundary_factor(indi) = sigma / (1D0 + zs(k))
							If (i .eq. 1) 			boundary_factor(indi) = 0.0d0
							If (i .eq. nxsz(k)) 	boundary_factor(indi) = 0.0d0
							If (j .eq. 1) 			boundary_factor(indi) = 0.0d0
							If (j .eq. nysz(k)) 	boundary_factor(indi) = 0.0d0
							indi = indi+1
						Enddo
					Enddo
				Enddo
			Endif
		Endif
		If (stitch_regions == .True.) Then
			If (have_diagonal .eq. 1) Then

				max_y = 0.0D0
				max_x = 0.0D0
				do dd = 1, my_n_diagonal
					indi = my_diagonal(dd,2)
					indj = my_diagonal(dd,3)
					test_x = abs(atax(indi,indj))
					test_y = abs(atay(indi,indj))
					if (test_x .gt. max_x) max_x = test_x
					if (test_y .gt. max_y) max_y = test_y
				enddo
				do dd = 1, my_n_diagonal
					k = my_diagonal(dd,1)
					stitch_weight(1,k) = 1.0D0-dabs(atax(indi,indj)/max_x)
					stitch_weight(2,k) = 1.0D0-dabs(atay(indi,indj)/max_y)
					!	stitch_weight(1:2,k) = stitch_weight(1:2,k)*(dvolume(nsolution)/dvolume(k))**(2.0D0)
					!	stitch_weight(1:2,k) = 1.0D0/dvolume(k)*sum(dvolume)
					!	stitch_weight(1,K) = (stitch_weight(1,k)*(1.0-stitch_min)+stitch_min)/sum(dvolume)
					!	stitch_weight(2,k) = (stitch_weight(2,k)*(1.0-stitch_min)+stitch_min)/sum(dvolume)
					stitch_weight(1:2,k) = 1.0D0*(dvolume(k)/sum(dvolume))*(zs_pack(k)+1.0D0)**2.0D0
					!	stitch_weight(1:2,k) = stitch_weight(1:2,k)*sqrt(zs_pack(k)+0.1)
				enddo



			Endif
			lh_amp = 0.0D0


			If (have_diagonal .eq. 1) Then  ! No need to do anything if we don't have any elements of the diagonal
				do i = 1, n_regions	! iterate over all the regions
					if (neighbors(r_id+1,i) .eq. 1) then ! this region (i) is my (r_id) neighbor
						ncount = p_count(r_id+1,i)  ! how many points overlap with this particular neighbor
						do j = 1, ncount
							k = partner_points(1,j,r_id+1,i) +1 ! this point of mine overlaps with one in region i
							        ! so augment my diagonal there.  Add 1 since IDL to Fortran.
							kk = partner_points(2,j,r_id+1,i) +1	! This is the index of the point in region i that I overlap with
							do dd = 1, my_n_diagonal
								If (k .eq. my_diagonal(dd,1)) Then
									indi = my_diagonal(dd,2)
									indj = my_diagonal(dd,3)
									! We multiply by boundary_factor(kk), which is zero if point
									! kk is a boundary point in region r_id
									ataxm(indi,indj) = ataxm(indi,indj)+2.0*stitch_weight(1,k)*boundary_factor(kk)
									ataym(indi,indj) = ataym(indi,indj)+2.0*stitch_weight(2,k)*boundary_factor(kk)
								Endif
							enddo
						enddo
					endif
				enddo
			Endif	! have_diagonal
		Endif	! stitch interior
	End Subroutine Add_Stitch

  Subroutine Update_RHS(iter,r_id)
	Implicit None
	Integer :: i, j, k, ncount, kk, iter, cc, ccc
	Integer*8 :: r_id
	Real*8 :: vtx, vty, vtx2, vty2, fweight, bweight
	Real*8, Allocatable :: stitch_count(:)
	rhsxyiter(:,:,:) = rhsxyzero(:,:,:)
	found_partner(:,r_id) = 0
	fweight = forward_weight
	bweight = 1.0D0-fweight
	if (iter .gt. 0) then
	! Assumes first regions is labeled 1 !  Sorry - some parts of code are zero, some are 1.  Keep track
	! I hate me.
	! here we augument the rhs using the partner points array

	if (stitch_regions == .True.) then
		!write(6,*)'Stitching Interior'
	do i = 1, n_regions	! iterate over all the regions
	   	if (neighbors(r_id,i) .eq. 1) then ! this region (i) is my (r_id) neighbor
			ncount = p_count(r_id,i)  ! how many points overlap with this particular neighbor
			do j = 1, ncount
			    k = partner_points(1,j,r_id,i) +1  ! this point of mine overlaps with one in region i
							        ! so augment my diagonal there.  Add 1 for IDL to Fortran
			    kk = partner_points(2,j,r_id,i) +1  ! this point in the neighboring region corresponds to my point k			
		       rhsxyiter(k,1,r_id) = rhsxyiter(k,1,r_id)+2.0*lastvx(kk,i)*stitch_weight(1,k)*boundary_factor(kk)	    
		       rhsxyiter(k,2,r_id) = rhsxyiter(k,2,r_id)+2.0*lastvy(kk,i)*stitch_weight(2,k)*boundary_factor(kk)
			enddo
		endif
	enddo


	endif ! stitch_regions eq true

	if (stitch_boundaries_only == .True.) then
		!write(6,*)' Stitching Boundaries'
	found_partner(:,r_id) = 0
	do i = 1, nsolution	! zero out the rhs wherever we have overlap in this region with interior point of another
		if (true_overlap(i,r_id) .eq. 1) then
		   rhsxyiter(i,1:2,r_id) = 0.0D0
		endif
	enddo
	allocate(stitch_count(1:nsolution))
	stitch_count(:) = 0.0D0
	do i = 1, n_regions	! iterate over all the regions
	   	if (neighbors(r_id,i) .eq. 1) then ! this region (i) is my (r_id) neighbor
			ncount = p_count(r_id,i)  ! how many points overlap with this particular neighbor
			do j = 1, ncount
			    k = partner_points(1,j,r_id,i) +1  ! this point of mine overlaps with one in region i
							        ! so augment my diagonal there.  Add 1 for IDL to Fortran
			    kk = partner_points(2,j,r_id,i) +1  ! this point in the neighboring region corresponds to my point k			
			    vtx = (lastvx(kk,i)*fweight+lastvx(k,r_id)*bweight)  ! "t" is for "target"
			    vty = (lastvy(kk,i)*fweight+lastvy(k,r_id)*bweight)  ! aim for a midpoint 			
			    if (isboundary(k)) then
				if (isboundary(kk) .eq. 0) then
					rhsxyiter(k,1,r_id) = rhsxyiter(k,1,r_id)+vtx
					rhsxyiter(k,2,r_id) = rhsxyiter(k,2,r_id)+vty
					stitch_count(k) = stitch_count(k)+1.0D0
					found_partner(k,r_id) = 1
				endif
			    endif
			enddo
		endif
	enddo
	do i = 1, nsolution	! Turn the sum over overlapping target velocities into an average
		if (true_overlap(i,r_id) .eq. 1) then
		   if(found_partner(i,r_id) .eq. 0) then
			Write(6,*)'There is a bug....'
		   endif
		   rhsxyiter(i,1,r_id) = rhsxyiter(i,1,r_id)/stitch_count(i)
		   rhsxyiter(i,2,r_id) = rhsxyiter(i,2,r_id)/stitch_count(i)
		endif
	enddo
	deallocate(stitch_count)
	endif ! Stitch_Regions eq true



	endif  ! iter gt 0


  End Subroutine Update_RHS

Subroutine Check_Bounds()
	Integer :: i, j, k, kk, ncount, r_id
	allocate(true_overlap(1:nsolution,1:n_regions))
	allocate(found_partner(1:nsolution,1:n_regions))
	true_overlap(:,:) = 0
	do r_id = 1, n_regions
	do i = 1, n_regions
		if (neighbors(r_id,i) .eq. 1) then ! region i is a neighbor of r_id
			ncount = p_count(r_id,i) ! this many points overlap between neighbors
			do j = 1, ncount
				 k = partner_points(1,j,r_id,i)+1
				kk = partner_points(2,j,r_id,i)+1
				if (isboundary(k) .eq. 1) then
				if (isboundary(kk) .eq. 0) then
					true_overlap(k,r_id) = 1
					!point k in r_id overlaps with a point in
					! another region which is NOT also a boundary point
				endif
				endif
			enddo
		endif
	enddo
	enddo

END Subroutine Check_Bounds

Subroutine Findmns
	Integer :: ind, i, j, k
	! Identify the points which are ahead of and behind each mesh point in the solution grid - 
	! both in x and y.  This will be the same for each sub-region since the grid's are identical
	! in structure.
	! first find the x neighbors
	ind = 1
	Allocate(x_neighbors(2,nsolution))
	Allocate(y_neighbors(2,nsolution))
	Allocate(xncount(nsolution))
	Allocate(yncount(nsolution))
	Allocate(deltax(nsolution))
	Allocate(deltay(nsolution))
	Allocate(isboundary(nsolution))
	xncount(:) = 0
	yncount(:) = 0
	isboundary(:) = 0
	do k = 1, nzs
		do j = 1, nysz(k)
			do i = 1, nxsz(k)
				if (i .eq. 1) isboundary(ind) = 1
				if (i .eq. nxsz(k)) isboundary(ind) = 1
				if (i .lt. nxsz(k)) then ! search for forward neighbors
				  x_neighbors(xncount(ind)+1,ind) = ind+1
				  xncount(ind) = xncount(ind)+1
				endif
				if (i .gt. 1) then
				  x_neighbors(xncount(ind)+1,ind) = ind-1
				  xncount(ind) = xncount(ind)+1
				endif
				deltax(ind) = deg_factor(k)*(lons_pack(2)-lons_pack(1))
				deltay(ind) = deg_factor(k)*(lats_pack(nxsz(1)+1)-lats_pack(1))
				ind = ind+1
			enddo
		enddo
	enddo

	! now find the y neighbors
	ind = 1
	do k = 1, nzs
		do j = 1, nysz(k)
			do i = 1, nxsz(k)
				if (j .eq. 1) isboundary(ind) = 1
				if (j .eq. nysz(k)) isboundary(ind) = 1
				if (j .lt. nysz(k)) then ! search for forward neighbors
				  !this point is an x neighbor
				  y_neighbors(yncount(ind)+1,ind) = ind+nxsz(k)
				  yncount(ind) = yncount(ind)+1
				endif
				if (j .gt. 1) then ! search for the reverse neighbors
				  y_neighbors(yncount(ind)+1,ind) = ind-nxsz(k)
				  yncount(ind) = yncount(ind)+1
				endif
				ind = ind+1
			enddo
		enddo
	enddo	
	if (n_regions .gt. 1) then
	Call Check_Bounds() ! establish whether points that overlap are boundary or interior points
	endif
End Subroutine Findmns

!!! New here
Subroutine Compute_Regularization_Matrices
	Implicit None
	Integer*8 :: xtotal, ytotal, ztotal
	Integer*8 :: ix,iy, iz, i, ile, i0, ile0,j, k 
	Integer*8 :: i2, ij_offset0, ij_offset1, ij_offset2
	Integer*8 :: deg_now, deg_next, deg_diff
	Integer*8 :: offset, offset0, ii, jj
        Integer*8 :: nsolnsol, ind1, ind2, inds(1:9), nind
	Real*8, Allocatable :: dxz(:), dyz(:), zss(:), reg_temp(:,:), temp_arr(:), bregm(:,:)
	Real*8 :: delz,delz4,delz6,delz9, zweight, ts(1:9)
	Real*8 :: f1, f2, t1,t2,t3, lv2, lh2
	xtotal = 0
	ytotal = 0
	ztotal = 0
	do iz = 1, nzs
		xtotal = xtotal+(nxsz(iz)-1)*nysz(iz)
	enddo
	do iz = 1, nzs
		ytotal = ytotal+(nysz(iz)-1)*nxsz(iz)
	enddo
	do iz = 1, nzs
		ztotal = ztotal+nysz(iz)*nxsz(iz)
	enddo
	
!	write(*,*) "WARNING: Allocating large matrix (bregm)"
!	Allocate(bregm(1:nsolution,1:nsolution))
!	Allocate(bregm(my_rows(1):my_rows(my_nrow),my_cols(1):my_cols(my_ncol)))
	Allocate(dyz(1:nzs))
	Allocate(dxz(1:nzs))
	Allocate(zss(1:nzs))
	do iz = 1, nzs
		zss(iz) = iz
	enddo
	dxz(:) = 6.10719*deg_factor(:)   ! hard coded in - might need to change that number
	dxz(:) = 1.0/dxz(:)*sqrt(dxz(:)/dxz(1))
	dyz(:) = dxz			! not sure what I'm really doing here.
!	dxz(:) = dxz(1)			!  This is a bit convoluted - check later
!	dyz(:) = dyz(1)			! for now, just match to idl for testing purposes

!	Write(6,*)'check 1 :', xtotal, dxz(1)
	i = 1
	ile = 1
!	Lx(:,:) = 0.0D0
!	bregm(:,:) = 0.0D0
	regm(:,:) = 0D0
	! x matrix
	Do iz = 1, nzs
		zweight = Sum(z_owner(iz,:)*wz(:))
!		write(6,*) 'dz ', zweight
		!lh2 = lh + deltalh*zs(iz)
		!lh2 = MAX(lh2,1D-5)
!		lh2 = 13D0 * lh / (zs(iz) + 4D0)**2D0
!		lh2 = lh*exp(-zs(iz)*0.3)
		lh2 = lh
		Do iy =1 , nysz(iz)
			Do ix = 1, nxsz(iz)-1
	!			Lx(i,ile) = -dxz(iz)*zweight
	!			Lx(i,ile+1) = dxz(iz)*zweight
				ind1 = ile
				t1 = -dxz(iz)*zweight
				ind2 = ile+1
				t2= dxz(iz)*zweight
				IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
					!bregm(ind1,ind1) = bregm(ind1,ind1)+lh*t1*t1
					regm(rowmap(ind1),colmap(ind1)) = regm(rowmap(ind1),colmap(ind1))+lh2*t1*t1
				ENDIF
				IF (rowmap(ind2).GE.1 .AND. rowmap(ind2).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
					regm(rowmap(ind2),colmap(ind2)) = regm(rowmap(ind2),colmap(ind2))+lh2*t2*t2
				ENDIF
				IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
					regm(rowmap(ind1),colmap(ind2)) = regm(rowmap(ind1),colmap(ind2))+lh2*t1*t2
				ENDIF
				IF (rowmap(ind2).GE.1 .AND. rowmap(ind2).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
					regm(rowmap(ind2),colmap(ind1)) = regm(rowmap(ind2),colmap(ind1))+lh2*t1*t2
				ENDIF
				i = i+1
				ile = ile+1
			Enddo
			ile = ile+1
		Enddo
	Enddo
	
	i = 1
	ile = 1

	! y matrix
	Do iz = 1, nzs
		zweight = Sum(z_owner(iz,:)*wz(:))
		!lh2 = lh + deltalh*zs(iz)
		!lh2 = MAX(lh2,1D-5)
!		lh2 = 13D0 * lh / (zs(iz) + 4D0)**2D0
!		lh2 = lh*exp(-zs(iz)*0.3)
		lh2 = lh
		Do iy =1 , nysz(iz)-1
			Do ix = 1, nxsz(iz)
			!	Ly(i,ile) = -dyz(iz)*zweight
			!	Ly(i,ile+nxsz(iz)) = dyz(iz)*zweight
			!	regm(ile,ile+nxsz(iz)) = regm(ile,ile+nxsz(iz)) + Ly(i,ile)*Ly(i,ile+nxsz(iz))
			!	regm(ile+nxsz(iz),ile) = regm(ile+nxsz(iz),ile)+Ly(i,ile)*Ly(i,ile+nxsz(iz))
				ind1 = ile
				t1 = -dyz(iz)*zweight
				ind2 = ile+nxsz(iz)
				t2= dyz(iz)*zweight
				IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
					regm(rowmap(ind1),colmap(ind1)) = regm(rowmap(ind1),colmap(ind1))+lh2*t1*t1
				ENDIF
				IF (rowmap(ind2).GE.1 .AND. rowmap(ind2).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
					regm(rowmap(ind2),colmap(ind2)) = regm(rowmap(ind2),colmap(ind2))+lh2*t2*t2
				ENDIF
				IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
					regm(rowmap(ind1),colmap(ind2)) = regm(rowmap(ind1),colmap(ind2))+lh2*t1*t2
				ENDIF
				IF (rowmap(ind2).GE.1 .AND. rowmap(ind2).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
					regm(rowmap(ind2),colmap(ind1)) = regm(rowmap(ind2),colmap(ind1))+lh2*t1*t2
				ENDIF

				i = i+1
				ile = ile+1
			Enddo
		Enddo
		ile = ile+nxsz(iz)
	Enddo

	! z matrix - more complicated
	If (nzs .gt. 1) Then
		i0 = ztotal
		ile0 = nsolution
		offset = 0
		do iz = nzs,2,-1
			delz = zss(iz)-zss(iz-1)
			delz6 = delz*6.0
			delz9 = delz*9.0
			delz4 = delz*4.0
			offset = offset+nxsz(iz)*nysz(iz)
			i = i0-offset+1           ! the row of the matrix for the this point
			ile = ile0-offset+1       ! the first point (x1,y1) at the deepest depth  : nsolution-nx(z=zmax)*ny(z=zmax) +1 (plus 1 because array starts at 1)
			deg_now = deg_factor(iz)
			deg_next = deg_factor(iz-1)
			deg_diff = deg_next-deg_now
			! changing vertical regularization
			!lv2 = lv + deltalv*zs(iz)
			!lv2 = MAX(lv2,1D-5)
!			lv2 = 13D0 * lv / (zs(iz) + 4D0)**2D0
!			lv2 = lv*exp(-zs(iz)*0.10)
			lv2 = lv
			do iy = 1, nysz(iz)
				do ix = 1, nxsz(iz)
!					Lz(i,ile) = 1.0/delz
					ind2 = ile
					t2 = 1.0/delz
					If (deg_diff .eq. 0) Then 
		!				Lz(i,ile-nxsz(iz)*nysz(iz))=-1.0/delz   ! if there's no difference in the grid, the points are separated by nx*ny
						ind1 = ile-nxsz(iz)*nysz(iz)
						t1 = -1.0/delz
						IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
							regm(rowmap(ind1),colmap(ind1)) = regm(rowmap(ind1),colmap(ind1))+lv2*t1*t1
						ENDIF
						IF (rowmap(ind2).GE.1 .AND. rowmap(ind2).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
							regm(rowmap(ind2),colmap(ind2)) = regm(rowmap(ind2),colmap(ind2))+lv2*t2*t2
						ENDIF
						IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
							regm(rowmap(ind1),colmap(ind2)) = regm(rowmap(ind1),colmap(ind2))+lv2*t1*t2
						ENDIF
						IF (rowmap(ind2).GE.1 .AND. rowmap(ind2).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
							regm(rowmap(ind2),colmap(ind1)) = regm(rowmap(ind2),colmap(ind1))+lv2*t1*t2
						ENDIF
					Else
						!scale back from previous layer
						i2 = ile0-offset-nxsz(iz-1)*nysz(iz-1)+1  ! i2 is the first point (x1,y1) at the next depth up
						! first assume we are at an interior point and
						! don't need to worry about the edges
						! now get the corresponding x,y one point up - current scheme says grid spacing changes by a factor of 2 when
						!   it changes - hence ix*2+1 (always need to move one extra point in - similar reasoning for the (2*(iy-1)+1) type terms.  
						! there will be 9 points 
						ij_offset0 = i2+(2*(iy-1)+1-1)*nxsz(iz-1)+(ix-1)*2 ! +1 !  row of the corresponding point
						ij_offset1 = i2+(2*(iy-1)+2-1)*nxsz(iz-1)+(ix-1)*2 ! +1 !nextrow
						ij_offset2 = i2+(2*(iy-1)-1)*nxsz(iz-1)+(ix-1)*2 ! +1 !previous row
						If ((ix .ne. 1) .and. (ix .ne. nxsz(iz)) )Then
							! interior x
							If ((iy .ne. 1) .and. (iy .ne. nysz(iz))) Then						
								!point interior in x and in y
!								Lz(i, ij_offset0-1) = -1.0/delz9
!								Lz(i, ij_offset0) = -1.0/delz9
!								Lz(i, ij_offset0+1) = -1.0/delz9
!								Lz(i, ij_offset1-1) = -1.0/delz9
!								Lz(i, ij_offset1) = -1.0/delz9
!								Lz(i, ij_offset1+1) = -1.0/delz9
!								Lz(i, ij_offset2-1) = -1.0/delz9
!								Lz(i, ij_offset2) = -1.0/delz9
!								Lz(i, ij_offset2+1) = -1.0/delz9


								nind = 9
								ts(1:9) = -1.0/delz9
								inds(1) = ij_offset0-1
								inds(2) = ij_offset0
								inds(3) = ij_offset0+1
								inds(4) = ij_offset1-1
								inds(5) = ij_offset1
								inds(6) = ij_offset1+1
								inds(7) = ij_offset2-1
								inds(8) = ij_offset2
								inds(9) = ij_offset2+1

							Endif
							If (iy .eq. 1) Then
								!interior x, lower row y
!								Lz(i, ij_offset0-1) = -1.0/delz6
!								Lz(i, ij_offset0) = -1.0/delz6
!								Lz(i, ij_offset0+1) = -1.0/delz6
!								Lz(i, ij_offset1-1) = -1.0/delz6
!								Lz(i, ij_offset1) = -1.0/delz6
!								Lz(i, ij_offset1+1) = -1.0/delz6

								nind = 6
								ts(1:6) = -1.0/delz6
								inds(1) = ij_offset0-1
								inds(2) = ij_offset0
								inds(3) = ij_offset0+1
								inds(4) = ij_offset1-1
								inds(5) = ij_offset1
								inds(6) = ij_offset1+1			
							Endif
							If (iy .eq. nysz(iz)) Then
								!interior x, upper row y
!								Lz(i, ij_offset0-1) = -1.0/delz6
!								Lz(i, ij_offset0) = -1.0/delz6
!								Lz(i, ij_offset0+1) = -1.0/delz6
!								Lz(i, ij_offset2-1) = -1.0/delz6
!								Lz(i, ij_offset2) = -1.0/delz6
!								Lz(i, ij_offset2+1) = -1.0/delz6

								nind = 6
								ts(1:6) = -1.0/delz6
								inds(1) = ij_offset0-1
								inds(2) = ij_offset0
								inds(3) = ij_offset0+1
								inds(4) = ij_offset2-1
								inds(5) = ij_offset2
								inds(6) = ij_offset2+1

							Endif
						Endif
						If (ix .eq. 1) Then
							!x end point - left edge
							If ((iy .ne. 1) .and. (iy .ne. nysz(iz))) Then
								!interior y point
!								Lz(i,ij_offset0) = -1.0/delz6
!								Lz(i,ij_offset0+1) = -1.0/delz6
!								Lz(i,ij_offset1) = -1.0/delz6
!								Lz(i,ij_offset1+1) = -1.0/delz6
!								Lz(i,ij_offset2) = -1.0/delz6
!								Lz(i,ij_offset2+1) = -1.0/delz6

								nind = 6
								ts(1:6) = -1.0/delz6
								inds(1) = ij_offset0
								inds(2) = ij_offset0+1
								inds(3) = ij_offset1
								inds(4) = ij_offset1+1
								inds(5) = ij_offset2
								inds(6) = ij_offset2+1
			
							Endif
							If (iy .eq. 1) Then
								! y end point - lower left corner
!								Lz(i,ij_offset0) = -1.0/delz4
!								Lz(i,ij_offset0+1) = -1.0/delz4
!								Lz(i,ij_offset1) = -1.0/delz4
!								Lz(i,ij_offset1+1) = -1.0/delz4								

								nind = 4
								ts(1:4) = -1.0/delz4
								inds(1) = ij_offset0
								inds(2) = ij_offset0+1
								inds(3) = ij_offset1
								inds(4) = ij_offset1+1								

							Endif
							If (iy .eq. nysz(iz)) Then
								! y end point - upper left corner
!								Lz(i,ij_offset0) = -1.0/delz4
!								Lz(i,ij_offset0+1) = -1.0/delz4
!								Lz(i,ij_offset2) = -1.0/delz4
!								Lz(i,ij_offset2+1) = -1.0/delz4								

								nind = 4
								ts(1:4) = -1.0/delz4
								inds(1) = ij_offset0
								inds(2) = ij_offset0+1
								inds(3) = ij_offset2
								inds(4) = ij_offset2+1

							Endif

						Endif
						If (ix .eq. nxsz(iz)) Then
							! x end point - right edge
							If ((iy .ne. 1) .and. (iy .ne. nysz(iz))) Then
								!interior y point
!								Lz(i,ij_offset0-1) = -1.0/delz6
!								Lz(i,ij_offset0) = -1.0/delz6
!								Lz(i,ij_offset1-1) = -1.0/delz6
!								Lz(i,ij_offset1) = -1.0/delz6	
!								Lz(i,ij_offset2-1) = -1.0/delz6
!								Lz(i,ij_offset2) = -1.0/delz6								

								nind = 6
								ts(1:6) = -1.0/delz6
								inds(1) = ij_offset0-1
								inds(2) = ij_offset0
								inds(3) = ij_offset1-1
								inds(4) = ij_offset1
								inds(5) = ij_offset2-1
								inds(6) = ij_offset2

							Endif
							If (iy .eq. 1) Then
								! y end point - lower right corner
!								Lz(i,ij_offset0-1) = -1.0/delz4
!								Lz(i,ij_offset0) = -1.0/delz4
!								Lz(i,ij_offset1-1) = -1.0/delz4
!								Lz(i,ij_offset1) = -1.0/delz4	

								nind = 4
								ts(1:4) = -1.0/delz4
								inds(1) = ij_offset0-1
								inds(2) = ij_offset0
								inds(3) = ij_offset1-1
								inds(4) = ij_offset1

							Endif
							If (iy .eq. nysz(iz)) Then
								!y end point - upper right corner
!								Lz(i,ij_offset0-1) = -1.0/delz4
!								Lz(i,ij_offset0) = -1.0/delz4
!								Lz(i,ij_offset2-1) = -1.0/delz4
!								Lz(i,ij_offset2) = -1.0/delz4	

								nind = 4
								ts(1:4) = -1.0/delz4
								inds(1) = ij_offset0-1
								inds(2) = ij_offset0
								inds(3) = ij_offset2-1
								inds(4) = ij_offset2

							Endif

						Endif
						IF (rowmap(ind2).GE.1 .AND. rowmap(ind2).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
							regm(rowmap(ind2),colmap(ind2)) = regm(rowmap(ind2),colmap(ind2))+lv2*t2*t2
						ENDIF
						do ii = 1, nind-1
						     ind1 = inds(ii)
							 IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
							     regm(rowmap(ind1),colmap(ind1)) = regm(rowmap(ind1),colmap(ind1))+lv2*ts(ii)*ts(ii)
							 ENDIF
							 IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
							     regm(rowmap(ind1),colmap(ind2)) = regm(rowmap(ind1),colmap(ind2))+lv2*ts(ii)*t2
							 ENDIF
							 IF (rowmap(ind2).GE.1 .AND. rowmap(ind2).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
							     regm(rowmap(ind2),colmap(ind1)) = regm(rowmap(ind2),colmap(ind1))+lv2*ts(ii)*t2
							 ENDIF
						     do jj = ii+1, nind
								 IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(inds(jj)).GE.1 .AND. colmap(inds(jj)).LE.my_ncol) THEN
						     		 regm(rowmap(ind1),colmap(inds(jj))) = regm(rowmap(ind1),colmap(inds(jj)))+lv2*ts(ii)*ts(jj)
								 ENDIF
								 IF (rowmap(inds(jj)).GE.1 .AND. rowmap(inds(jj)).LE.my_nrow .AND. colmap(ind1).GE.1 .AND. colmap(ind1).LE.my_ncol) THEN
									 regm(rowmap(inds(jj)),colmap(ind1)) = regm(rowmap(inds(jj)),colmap(ind1))+lv2*ts(ii)*ts(jj)
								 ENDIF
						     enddo
					
						enddo
						ii = nind
						ind1 = inds(nind)
						IF (rowmap(ind1).GE.1.AND.rowmap(ind1).LE.my_nrow.AND.colmap(ind1).GE.1.AND.colmap(ind1).LE.my_ncol) THEN
							regm(rowmap(ind1),colmap(ind1)) = regm(rowmap(ind1),colmap(ind1))+lv2*ts(ii)*ts(ii)
						ENDIF
						IF(rowmap(ind1).GE.1.AND.rowmap(ind1).LE.my_nrow.AND.colmap(ind2).GE.1.AND.colmap(ind2).LE.my_ncol)THEN
							regm(rowmap(ind1),colmap(ind2)) = regm(rowmap(ind1),colmap(ind2))+lv2*ts(ii)*t2
						ENDIF
						IF(rowmap(ind2).GE.1.AND.rowmap(ind2).LE.my_nrow.AND.colmap(ind1).GE.1.AND.colmap(ind1).LE.my_ncol)THEN
							regm(rowmap(ind2),colmap(ind1)) = regm(rowmap(ind2),colmap(ind1))+lv2*ts(ii)*t2
						ENDIF


					Endif

					i = i+1			! move to next point in at this depth
					ile = ile+1
				enddo
			enddo
		enddo

	Endif
!	regm(:,:) = 0.0D0
!	Do j = 1, my_ncol
!		Do i = 1, my_nrow
!			if (colmap(my_cols(j)) .ne. j) then
!				write(*,*), j, colmap(j), my_cols(j), "STOP"
!				stop
!			endif
			!regm(i,j) = bregm(my_rows(i), my_cols(j))
!		Enddo
!	Enddo

!	DeAllocate(bregm)
!	write(*,*) "bregm deallocated"
	DeAllocate(dyz)
	DeAllocate(dxz)
End Subroutine Compute_Regularization_Matrices

SUBROUTINE IncrementRegm(ind1, ind2, val)
	INTEGER*8 :: ind1, ind2
	REAL*8 :: val
	IF (rowmap(ind1).GE.1 .AND. rowmap(ind1).LE.my_nrow .AND. colmap(ind2).GE.1 .AND. colmap(ind2).LE.my_ncol) THEN
		regm(rowmap(ind1),colmap(ind2)) = regm(rowmap(ind1),colmap(ind2))+val
	ENDIF
END SUBROUTINE IncrementRegm



	Subroutine Write_Matrix_To_Disk(prefix,mat,mlocal,nlocal,rows_local,columns_local, mtotal,ntotal,rid)
	  Real*8, Intent(In) :: mat(:,:)
	  Integer, Intent(In) :: rows_local(:), columns_local(:), mlocal,nlocal
	  Integer*8, Intent(In) :: mtotal, ntotal
	  Integer*8 :: i, j, myid, rid
	  Character*120 :: prefix, mfile
	  Character*4 :: idstring,ridstring
	  myid = myrow+nprow*mycol
	  write(idstring,'(i4.4)')myid
	  write(ridstring,'(i4.4)')rid
	  mfile= TRIM(ADJUSTL(prefix))//'_cpu'//TRIM(ADJUSTL(idstring))//'_r'//TRIM(ADJUSTL(ridstring))

          Open(unit=13, file=mfile,form='unformatted')
	  If (myid .eq. 0) Then
		write(13)nprocs
		write(13)mtotal, ntotal
	  Endif
          Write(13) mlocal, nlocal
          Write(13) ((mat(i,j),i=1,mlocal),j=1,nlocal)
          Write(13) (rows_local(i),i=1,mlocal)
	  Write(13) (columns_local(j),j=1,nlocal)
	  close(13)
	End Subroutine Write_Matrix_To_Disk
  Subroutine Pinv2_Parallel(inv_order,r_id)
    implicit none
	real*8, Allocatable :: a(:,:)
        real*8, Allocatable :: b(:), y(:)
        real*8, Allocatable :: bigb(:), biga(:,:),btemp(:,:), little_ata(:,:), dtemp(:,:), etemp(:,:)
	real*8, Allocatable :: work(:), covar(:,:), covar_piece(:,:), res_ker(:,:)
	real*8, allocatable :: covardiag(:)
	real*8 :: erh2
    !	  integer ipiv(8000)  !10 + 2
    integer, Allocatable :: ipiv(:)
	character*120 :: ofits, res_ker_file
    integer :: lwork, inv_order, nipiv
    integer*8 :: r_id
    integer, Allocatable :: iwork(:)
    integer :: locc, locr
	Integer*8 :: ix,iy, iz, ile

    integer desca(9)
    integer descb(9)
    integer descy(9)
    integer desce(9)
    integer descd(9)
    integer descla(9)
    integer descrk(9)
    integer descat(9)  
    integer i

    integer ihi
    integer ilo
    integer info
    integer :: nsol

    integer istore
    integer j
    integer jhi
    integer jlo
    integer jstore
    integer m

    integer n, na
    integer nrhs
    integer nbrhs

    integer prow
    integer pcol
    integer liwork

    integer :: yr_nrow, yr_nblocks
    integer, Allocatable :: yr_blocks_y(:), xr_blocks_x(:)
    integer :: a_ind_i, a_ind_j, biga_ind_i, biga_ind_j,xr_ncol, xr_nblocks
    integer :: ii,jj, nrow, ncol, bind_j, bind_i, nr, nc
    Character*120 :: prefix

    nipiv = my_nrow+mb
    allocate(ipiv(nipiv))
    allocate(biga(1:my_nrow,1:my_ncol))	
		
	If (inv_dir .eq. 1) Then
		! x direction
		if (inv_order == 0) Then
		biga(:,:) = atax(:,:)
		else
		biga(:,:) = ataxm(:,:)
		endif
!		erh2 = xerh(r_id+1)*xerh(r_id+1)
		        erh2 = erhx**2     
        	
        Endif
	If (inv_dir .eq. 2) Then
		! y direction
		if (inv_order == 0) Then
		biga(:,:) = atay(:,:)
		else
		biga(:,:) = ataym(:,:)
		endif
!		erh2 = yerh(r_id+1)*yerh(r_id+1)
		erh2 = erhy**2
        Endif
!	erh2 = 1.0
	nsol = nsolution

        n = nsol
	na = nmeasurements
        m = nsol
  
        if ( nprocs < 1 ) then
          if ( iam == 0 ) then
            nprocs = nprow * npcol
          end if
        end if


	call blacs_gridinfo ( context, nprow, npcol, myrow, mycol )


	if ( myrow == -1 ) then
	  call blacs_exit ( 0 )
          stop
	end if

 
	Allocate(a(1:my_nrow,1:my_ncol))
	Allocate(dtemp(1:my_nrow,1:my_ncol))
	Allocate(etemp(1:my_nrow,1:my_ncol))
	Allocate(little_ata(1:my_nrow,1:my_ncol))
	
	little_ata(:,:) = 0.0
	
	! we need some work arrays for the matrix inversion
	locr = my_nrow*nsolution
	locc = my_ncol*nsolution

	!write(6,*)myrow, locr
	!write(6,*)mycol, locc	
	liwork = nsolution !only true if nprow = npcol
	lwork = nb*my_nrow
	allocate(work(1:lwork))
	allocate(iwork(1:liwork))


	!Everyone copies their piece of the big matrices into a smaller array

	ile = 1
	Do iz = 1, nzs
		!lh2 = lh / (zs(iz)+1.0)
		Do iy =1 , nysz(iz)
			Do ix = 1, nxsz(iz)-1
!				WRITE(*,'(I5,3E15.4)') ile, biga(ile,ile), regm(ile,ile), zs(iz)
				ile = ile + 1
			ENDDO
		ENDDO
	ENDDO

	do j = 1, my_ncol
		do i = 1, my_nrow
!			if (biga(i,j) .ne. 0.D0 .or. regm(i,j) .ne. 0D0) then
!			if (i .eq. j) then
!				WRITE(*,'(I5,3E15.4)') i, biga(i,j), regm(i,j), 
!			endif
			a(i,j) = biga(i,j)+regm(i,j)

			little_ata(i,j) = biga(i,j) 	
		enddo
	enddo

	! Replace rows of the matrix with 1.0 to enforce a boundary
	! condition of sorts
	If (inv_order ==1) Then
		do i = 1, my_nrow
		!------------------
		  if (isboundary(my_rows(i))==1) then
		 ! if (have_match(my_rows(i),r_id+1)==1) then
		  if (true_overlap(my_rows(i),r_id+1)==1) then
				a(i,:) = 0.0D0
				do j =1, my_ncol
					if (my_cols(j) .eq. my_rows(i)) then
						a(i,j) = 1.0D0
					endif
				enddo
		  endif
		 ! endif
		  endif
		enddo
	endif

	
	nrhs = 1

	!  Set the array descriptors
	!
	call descinit ( desca, m, n, mb, nb, 0, 0, context, my_nrow, info )
	call descinit ( descd, m, n, mb, nb, 0, 0, context, my_nrow, info )
	call descinit ( desce, m, n, mb, nb, 0, 0, context, my_nrow, info )
	call descinit ( descla, m, n, mb, nb, 0, 0, context, my_nrow, info )
	call descinit ( descrk, m, na, mb, nba, 0, 0, context, my_nrow, info)
	call descinit ( descat, m, na, mb, nba, 0, 0, context, my_nrow, info)
!	call descinit ( descb, n, nrhs, nb, 1, 0, 0, context, my_nrow, info )
!	call descinit ( descy, n, nrhs, nb, 1, 0, 0, context, my_nrow, info )


	 if ((myrow .eq. 0) .and. (mycol .eq. 0)) then
		write(6,*)'Performing parallel inversion.'
	 endif
	! Get the inverse matrix
	call pdgetrf( nsolution, nsolution, a, 1, 1, desca, ipiv, info )
	! at this point, i believe a() is correct. problem lies elsewhere

	!call pdpotrf( 'U', nsolution, a, 1, 1, desca, info )
	! write(6,*)'factor info is ', info
	lwork = -1
	liwork = -1
	! figure out what lwork and liwork should be
	call pdgetri(nsolution, a, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)  
	!call pdpotri( 'U', nsolution, a, 1, 1, desca, info )
	lwork = work(1)
	liwork = iwork(1)
	! write(6,*), 'lwork is ', lwork, nsolution
	! write(6,*), 'liwork is ', liwork, nsolution
	deallocate(iwork)
	allocate(iwork(liwork))
	deallocate(work)
	allocate(work(lwork))
	call pdgetri(nsolution, a, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)
	!call pdpotri( 'U', nsolution, a, 1, 1, desca, info )
	if ((myrow .eq. 0) .and. (mycol .eq. 0)) then
		write(6,*)'Inversion Complete. Inv info is ', info, r_id
		write(6,*)'Computing Covariance Matrix...'
	endif


	!  Compute the covariance Matrix
	! first get ata#transpose(a_inv)
	dtemp(:,:) = 0.0

	call PDGEMM( 'N','T', nsol, nsol,nsol, 1.0, little_ata, 1, 1, DESCLA, A, 1, 1, DESCA, 0.0, DTEMP, 1, 1, DESCD )
	! now get a_inv#ata#transpose(a_inv)
	etemp(:,:) = 0.0
	call PDGEMM( 'N','N', nsol, nsol,nsol, 1.0, a, 1, 1, DESCA, Dtemp, 1, 1, DESCD, 0.0, ETEMP, 1, 1, DESCE)
 	
	! covariance matrix is stored in etemp
	if ((myrow .eq. 0) .and. (mycol .eq. 0)) write(6,*)'...complete.'
	if (inv_order == 0) then
	    If (compute_coarse_kernels == .True.) Then
		if ((myrow .eq. 0) .and. (mycol .eq. 0)) write(6,*)'Computing Coarse Averaging Kernels'
		!  Compute the Resolution Kernels
		! these should just be a_inv#ata
	    	dtemp(:,:) = 0.0
	    	call PDGEMM( 'N','N', nsol, nsol,nsol, 1.0, a, 1, 1, DESCA, little_ata, 1, 1, DESCLA, 0.0, DTEMP, 1, 1, DESCD)
		If (inv_dir ==1 ) prefix = 'aux/avg_ker_unres_x'
	  	If (inv_dir ==2 ) prefix = 'aux/avg_ker_unres_y'
	        Call Write_Matrix_To_Disk(prefix,dtemp,my_nrow,my_ncol,my_rows,my_cols, nsolution,nsolution,r_id)
		if ((myrow .eq. 0) .and. (mycol .eq. 0)) write(6,*)'...complete'
	    Endif
	endif
	! save the results for this region  
	if (inv_order == 0) then
	  covar0(:,:,r_id+1,inv_dir) = Etemp(:,:)
	  ainv0(:,:,r_id+1,inv_dir) = a(:,:)	
	endif
	if (inv_order == 1) then
	  ainvn(:,:,r_id+1,inv_dir) = a(:,:)	
	endif

	!Compute the averaging kernel coefficients
	If (compute_avg_kernels == .True.) Then
	IF ((myrow .eq. 0) .and. (mycol .eq. 0) ) write(6,*)'Computing averaging kernel coefficients...'
	if (inv_order ==0) Then
	  reskercoefs(:,:) = 0.0D0
	  If (inv_dir == 1) Then
	   prefix = 'aux/avg_coefs_x'
	   call PDGEMM( 'N','N', nsol, na,nsol, 1.0, a, 1, 1, DESCA, atransposex, 1, 1, DESCAT, 0.0, reskercoefs, 1, 1, DESCRK)
	  Endif
	  If (inv_dir == 2) Then
	   prefix = 'aux/avg_coefs_y'
	   call PDGEMM( 'N','N', nsol, na,nsol, 1.0, a, 1, 1, DESCA, atransposey, 1, 1, DESCAT, 0.0, reskercoefs, 1, 1, DESCRK)
	  Endif
	  call Write_Matrix_To_Disk(prefix,reskercoefs,my_nrow,my_ncola,my_rows,my_colsa, nsolution,nmeasurements,r_id)

	Endif
        IF ((myrow .eq. 0) .and. (mycol .eq. 0) ) write(6,*)'...complete'
	! Everyone gets the covariance matrix and computes the errors for the zero order inversion
        ! we'll figure out how to get the errors for the non zero inversions later....
	Endif
!	write(*,*) "WARNING: Allocating large matrix (covar)"
!	allocate(covar(1:nsolution,1:nsolution))
	allocate(covardiag(1:nsolution))
	!if (inv_order == 0) then   ! covariance block --------------------------------------

	  do j = 1, my_ncol
		do i = 1, my_nrow
			if (my_rows(i) .eq. my_cols(j)) then
				covardiag(my_rows(i)) = etemp(i,j)
			endif
			!covar(my_rows(i),my_cols(j)) = etemp(i,j)
		enddo
	  enddo
	!----------------------------------------
	!		Distribute the Covariance Matrix
	!----------------------------
	if ( (myrow .eq. 0) .and. (mycol .eq. 0) ) write(6,*)'Distributing Covariance Matrix'
	do nr = 0, nprow-1
		do nc = 0, npcol-1
   			If ((myrow .eq. nr) .and. (mycol .eq. nc) ) Then
	    !broadcast my results to everyone else
	   !  write(6,*)'sending ', my_nrow, my_ncol, myrow, mycol
				call DGEBS2D( context, 'ALL', ' ',my_nrow, my_ncol, etemp(1,1), my_nrow )
			Else
	    		yr_nblocks = y_block_per_pe !number of rhs blocks this row has to contribute to the end result
     	    	xr_nblocks = x_block_per_pe !number of rhs blocks this row has to contribute to the end result
   	        	If (y_block_mod .gt. 0) Then
       	    		If (nr .lt. y_block_mod) Then ! Not .le. because myrow starts at 0
						yr_nblocks = yr_nblocks+1
           			Endif
           		Endif
           		If (x_block_mod .gt. 0) Then
           			If (nc .lt. x_block_mod) Then ! Not .le. because myrow starts at 0
        				xr_nblocks = xr_nblocks+1
           			Endif
           		Endif

           		Allocate(yr_blocks_y(1:yr_nblocks))
	            yr_nrow = 0 
	            Allocate(xr_blocks_x(1:xr_nblocks))
	            xr_ncol = 0 
    	        do j = 1, yr_nblocks
					yr_blocks_y(j) = nr +1 +(j-1)*nprow ! plus 1 since myrow starts at 0
					  ! but block inds start at 1
    	       		yr_nrow = yr_nrow+block_size_y(yr_blocks_y(j))
       	    	enddo
        	   	do j = 1, xr_nblocks
           			xr_blocks_x(j) = nc +1 +(j-1)*npcol ! plus 1 since myrow starts at 0
 					  ! but block inds start at 1
           			xr_ncol = xr_ncol+block_size_x(xr_blocks_x(j))
           		enddo

           		Allocate(covar_piece(1:yr_nrow,1:xr_ncol))
          !  write(6,*)'receiving ', yr_nrow, xr_ncol, nr, nc
           		call DGEBR2D( context, 'ALL', ' ',yr_nrow, xr_ncol, covar_piece(1,1), yr_nrow,nr, nc)

				do jj = 1, xr_nblocks
					do ii =1, yr_nblocks
	        			bind_j = xr_blocks_x(jj)
	        			bind_i = yr_blocks_y(ii)
	        			do j = 1, block_size_x(bind_j)
	           				do i = 1, block_size_y(bind_i)
								a_ind_j = (jj-1)*nb+j
								a_ind_i = (ii-1)*mb+i
								biga_ind_j = (xr_blocks_x(jj)-1)*nb+j
								biga_ind_i = (yr_blocks_y(ii)-1)*mb+i
								if (biga_ind_i .eq. biga_ind_j) then
									! load ement into diagonal of covar
									covardiag(biga_ind_i) = covar_piece(a_ind_i,a_ind_j)
								endif
								!covar(biga_ind_i,biga_ind_j) = covar_piece(a_ind_i,a_ind_j)
	          				enddo
	         			enddo
	       			enddo
	     		enddo
	     		DeAllocate(covar_piece)
	     		DeAllocate(xr_blocks_x)
	     		DeAllocate(yr_blocks_y)
	   		Endif
		enddo
	enddo

	!Now that we have the covariance matrix, copy the diagonal into esq
	If (inv_order == 0) then
		! These are the zero-order iteration errors
	if (inv_dir == 1) then
	  do j = 1, nsol
	      !esqx(j,r_id+1) = covar(j,j)*erh2
		  esqx(j,r_id+1) = covardiag(j)*erh2
	  enddo
	else
 	  do j = 1, nsol
   	  	  !esqy(j,r_id+1) = covar(j,j)*erh2
	      esqy(j,r_id+1) = covardiag(j)*erh2
	  enddo
	endif
	Endif
	If (inv_order == 1) then
		! These errors form the basis for the iterative errors
	if (inv_dir == 1) then
	  do j = 1, nsol
	      esqx_iter(j,r_id+1) = covardiag(j)*erh2
	  enddo
	else
 	  do j = 1, nsol
	      esqy_iter(j,r_id+1) = covardiag(j)*erh2
	  enddo
	endif
	Endif

	
	
!	DeAllocate(covar)
!	write(*,*) "covar deallocated"
	deallocate(covardiag)
	Deallocate(biga)
	Deallocate(a)
	Deallocate(little_ata)
	Deallocate(dtemp)
	Deallocate(etemp)
	Deallocate(work)
	Deallocate(iwork)
	DeAllocate(ipiv)
  End Subroutine Pinv2_Parallel
!  Subroutine Condense_Partner_Points()
!	! There are times when it is useful 
!  End Subroutine Condense_Partner_Points

	Subroutine Advance_Errors()
		Integer*8 :: r,rr,i,j,jj,k,kk,ncount
		Integer :: ii, bind_i, a_ind_i, biga_ind_i
		Integer, Allocatable :: yr_blocks_y(:)
		Integer :: yr_nrow, yr_nblocks
		Real*8, Allocatable :: my_modx(:,:), my_mody(:,:), emod_global(:,:), my_emod(:,:), btemp(:,:)
		Real*8 :: aweightx, aweighty
		Real*8 :: eweightx, eweighty
		write(6,*)'propagating errors'

		! take ainvn which we have saved and use the rows to propagate
		! the errors
		! here we propagate the errors

		! for this to function, we'll need to compute eweight(j,rr)
		!  with eweight(1:my_nrow,n_regions)

		! Step 1:  Because matrices are distributed across columns and rows
		!	   things are complicated a little bit.
		!	   compute the modification to esq_iter that each processor 
		!	   knows.
		!          Some organization of the partner points relative to which
		!	     processor has which columns would be useful here
		!	     But just brute force it for now to get it working.
		!	Allocate(my_modx(my_nrow,n_regions))	! my modification to the square of the x-error
		!	Allocate(my_mody(my_nrow,n_regions))	! my modification to the square of the y-error

		Allocate(my_emod(my_nrow*2,n_regions))	! actually just put the x and y errors in 1 array so we only have to
		! communicate one array.


		my_emod(:,:) = 0.0D0
		do r = 1, n_regions
			do rr = 1, n_regions
				if (neighbors(r,rr)) then
					ncount = p_count(r,rr)       
					do i = 1, my_nrow
						do j = 1, my_ncol
							aweightx = ainvn(i,j,r,1)**2
							aweighty = ainvn(i,j,r,2)**2
							do jj = 1, ncount

								kk = partner_points(2,jj,r,rr)+1	! point in subregion rr that overlaps with a points in region r
								k = partner_points(1,jj,r,rr)+1  ! This point in region r corresponds to point kk  in region rr
								eweightx = (2.0*stitch_weight(1,k))**2
								eweighty = (2.0*stitch_weight(2,k))**2
								if (my_cols(j) .eq. k) then			      
									! my_modx(i,r) = my_modx(i,r)+esqx_last(k,rr)*aweightx*eweight(k,rr)
									! my_mody(i,r) = my_mody(i,r)+esqy_last(k,rr)*aweighty*eweight(k,rr)
									my_emod(i,r)         =  my_emod(i,r)         + esqx_last(kk,rr)*(boundary_factor(kk)**2) * aweightx  * eweightx
									my_emod(i+my_nrow,r) =  my_emod(i+my_nrow,r) + esqy_last(kk,rr)*(boundary_factor(kk)**2) * aweighty  * eweighty
								endif
							enddo
						enddo
					enddo
				endif
			enddo
		enddo
		! Now, we need to communicate
		Allocate(emod_global(1:nsolution*2,n_regions))
		emod_global(:,:) = 0.0D0
		do i = 0, nprow-1
			do r = 0, npcol-1
				If ((myrow .eq. i) .and. (mycol .eq. r) ) Then
					!broadcast my results to everyone else
					allocate(btemp(1:my_nrow*2,1:n_regions))
					btemp(:,:) =  my_emod(:,:)  
					call DGEBS2D( context, 'ALL', ' ',2*my_nrow, n_regions, my_emod(1,1), 2*my_nrow )
					do rr = 1, my_nrow
						emod_global(my_rows(rr),:) = emod_global(my_rows(rr),:)+my_emod(rr,:)
						emod_global(my_rows(rr)+nsolution,:) = emod_global(my_rows(rr)+nsolution,:)+my_emod(rr+my_nrow,:)
					enddo
				Else
					! receive
					yr_nblocks = y_block_per_pe !number of rhs blocks this row has to contribute to the end result
					If (y_block_mod .gt. 0) Then
						If (i .lt. y_block_mod) Then ! Not .le. because myrow starts at 0
							yr_nblocks = yr_nblocks+1
						Endif
					Endif
					allocate(yr_blocks_y(1:yr_nblocks))
					yr_nrow = 0 
					do j = 1, yr_nblocks
						yr_blocks_y(j) = i +1 +(j-1)*nprow ! plus 1 since myrow starts at 0
						! but block inds start at 1
						yr_nrow = yr_nrow+block_size_y(yr_blocks_y(j))
					enddo
					Allocate(btemp(1:yr_nrow*2,n_regions))
					call DGEBR2D( context, 'ALL', ' ',yr_nrow*2, n_regions, btemp(1,1), yr_nrow*2,i, r)
					do ii =1, yr_nblocks
						bind_i = yr_blocks_y(ii)
						do j = 1, block_size_y(bind_i)     
							a_ind_i = (ii-1)*mb+j
							biga_ind_i = (yr_blocks_y(ii)-1)*mb+j
							emod_global(biga_ind_i,:) = emod_global(biga_ind_i,:)+btemp(a_ind_i,:)
							emod_global(biga_ind_i+nsolution,:) = emod_global(biga_ind_i+nsolution,:)+btemp(a_ind_i+yr_nrow,:)
						enddo
					enddo
					DeAllocate(yr_blocks_y)  
				Endif
				DeAllocate(btemp)
				call blacs_barrier(context,'A')
			enddo
		enddo

		esqx(1:nsolution,:) = esqx_iter(1:nsolution,:)+emod_global(1:nsolution,:)
		esqy(1:nsolution,:) = esqy_iter(1:nsolution,:)+emod_global(nsolution+1:2*nsolution,:)
		DeAllocate(emod_global)
		DeAllocate(my_emod)	

	End Subroutine Advance_Errors

	Subroutine RHS_Solve_Parallel(iter, r_ind)
	Implicit None
	Real*8, Allocatable :: bigb(:)
	Real*8, Allocatable :: a(:,:), y(:), btemp(:,:), b(:)
	Integer :: i,ii, bind_i, a_ind_i, biga_ind_i, nrhs, m, n, iter
	Integer :: j, info
	Integer*8 :: r_ind
	Integer, Allocatable :: yr_blocks_y(:)
	Integer :: yr_nrow, yr_nblocks
	Integer :: desca(9)
        Integer :: descb(9)
        Integer :: descy(9)
!	mb = 64   ! block size
!	nb = 64
	Allocate(bigb(1:nsolution))
	bigb(:) = rhsxyiter(:,inv_dir,r_ind)
	! write(6,*)'test ', rhsxysave(1,1,:)
	Allocate(a(1:my_nrow,1:my_ncol))
	Allocate(b(1:my_nrow))

	Allocate(y(1:my_nrow))
	if (iter == 0) then
		a(:,:) = ainv0(:,:,r_ind,inv_dir)		
	endif
	if (iter .gt. 0) then
		a(:,:) = ainvn(:,:,r_ind,inv_dir)

	endif
	
	
	do i = 1, my_nrow
		b(i) = bigb(my_rows(i))
	enddo


	nrhs = 1
	m = nsolution
        n = nsolution

	call descinit ( desca, m, n, mb, nb, 0, 0, context, my_nrow, info )
	call descinit ( descb, n, nrhs, nb, 1, 0, 0, context, my_nrow, info )
	call descinit ( descy, n, nrhs, nb, 1, 0, 0, context, my_nrow, info )

	  
        call PDGEMV( 'N', nsolution, nsolution, 1.0, A, 1, 1, DESCA, b, 1, 1, DESCB, 1, 0.0, Y, 1, 1, DESCY, 1 ) 


	do i = 1, my_nrow
		bigb(my_rows(i)) = y(i)
	enddo

	do i = 0, nprow-1
          If ((myrow .eq. i) .and. (mycol .eq. 0) ) Then
	    !broadcast my results to everyone else
	     allocate(btemp(1:my_nrow,2))
	     btemp(:,1) =  y(:)   ! made a 2d array because not sure if DGEBS2D works with 1d arrays
             call DGEBS2D( context, 'ALL', ' ',my_nrow, 2, btemp(1,1), my_nrow )
	  Else
	    ! receive
     	    yr_nblocks = y_block_per_pe !number of rhs blocks this row has to contribute to the end result
            If (y_block_mod .gt. 0) Then
              If (i .lt. y_block_mod) Then ! Not .le. because myrow starts at 0
	        yr_nblocks = yr_nblocks+1
              Endif
            Endif
            allocate(yr_blocks_y(1:yr_nblocks))
            yr_nrow = 0 
            do j = 1, yr_nblocks
              yr_blocks_y(j) = i +1 +(j-1)*nprow ! plus 1 since myrow starts at 0
 					  ! but block inds start at 1
              yr_nrow = yr_nrow+block_size_y(yr_blocks_y(j))
	    enddo
            Allocate(btemp(1:yr_nrow,2))
            call DGEBR2D( context, 'ALL', ' ',yr_nrow, 2, btemp(1,1), yr_nrow,i, 0)
	    do ii =1, yr_nblocks
              bind_i = yr_blocks_y(ii)
                do j = 1, block_size_y(bind_i)     
	  	  a_ind_i = (ii-1)*mb+j
		  biga_ind_i = (yr_blocks_y(ii)-1)*mb+j
		  bigb(biga_ind_i) = btemp(a_ind_i,1)
                enddo
            enddo
	    DeAllocate(yr_blocks_y)  
    	  Endif
	  DeAllocate(btemp)
          call blacs_barrier(context,'A')
	enddo
	if (inv_dir == 1) then
	  vx(:,r_ind) = bigb(:)
	else
	  vy(:,r_ind) = bigb(:)
	endif
	DeAllocate(a)
	DeAllocate(b)
	DeAllocate(y)
	DeAllocate(bigb)
	End Subroutine RHS_Solve_Parallel
End Module Matrix_Magic
