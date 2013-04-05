Module Data_Management
	USE Solution_Grid
	USE Region_Info
	USE Kernel_Management
	USE Run_Info
	! Contains global variables and subroutines related to the managment and
	!   organization of the ring analysis measurements (i.e. the data). 
	! Much of the organization is done prior to the fortran code via the HMI pipeline
	!   or (in the case of MDI data) IDL. 

	Implicit None
	Integer*8, Public :: ntsizes
	INTEGER*8, PUBLIC, ALLOCATABLE :: ndata(:), n_n(:), n_k(:)
	Logical, Public :: compute_res_ker_unres = .True.
	Integer*8, Allocatable, Public :: nktindex(:,:,:,:), indic(:,:), indicn(:), kxy(:,:,:), tsizes(:)
	Integer*8, Allocatable, Public :: ki(:,:), n(:,:), tile_sizes(:)
	CHARACTER*300, ALLOCATABLE, PUBLIC :: data_files(:)
	INTEGER*8, ALLOCATABLE, PUBLIC :: lon_all(:,:), lat_all(:,:)
        Integer*8, Allocatable :: nzr_indices(:,:,:), nzc_indices(:,:,:), ilocal_indices(:,:,:), jlocal_indices(:,:,:)
	Integer*8, Allocatable :: nznr_num(:,:), nznc_num(:,:)
	Integer*8, Public, ALLOCATABLE :: nlon(:), nlat(:)  !  , grid_scheme
	Integer*8, Allocatable, Private :: nvals(:,:), kvals(:,:)
	Real*8, Allocatable, Public :: latvals(:), lonvals(:)
	Real*8, Allocatable, Public :: uxy(:,:,:)
	real*8, Allocatable :: lon(:), lat(:)

	Real*8, Public :: lh, lv, deltalv, deltalh, sigma, lat_lims0(2), lon_lims0(2), max_tile_size, tsize_deg
    	Namelist /Input_Namelist/ x_block_size, y_block_size, x_block_size_a,n_regions, &
		kernel_set, ntsizes, lh, lv, deltalv, deltalh, sigma, dlonlat, grid_scheme, tile_sizes, max_iter, &
		stitch_boundaries_only, stitch_regions, propagate_errors, compute_avg_kernels, &
		compute_coarse_kernels, forward_weight, onesided_boundary_stitching, &
		grid_file, data_files

	CONTAINS
	
	SUBROUTINE Main_Initialization()
		IMPLICIT NONE
		! 1.  Gather the inversion parameters
		CALL Read_Run_Info()
		! 2. Initialize the kernels  - need kernel parameters to set and validate the solution grid
		CALL Initialize_Kernels()	
 		! 3.  Get the basic solution grid
		CALL Get_Solution_Grid(lon_lims0, lat_lims0, max_tile_size)

		CALL Pack_Grid()

		IF (n_regions > 1) THEN
			! 4. Of the different regions, figure out which ones overlap
			CALL Identify_Neighbors(lon_lims0,lat_lims0)
			CALL find_overlap(lon_lims0,lat_lims0)
		ENDIF
	END SUBROUTINE Main_Initialization
 
	! Reads run parameters from namelist
	! path is hard-coded
	SUBROUTINE Read_Run_Info()
		IMPLICIT NONE
		CHARACTER*120 :: run_file, region_file
		INTEGER*8 :: i, mxts = 20
		! Initialize lh, lv, lsigma, n_regions
		! then read the latitude and longitude limits
		ALLOCATE(tile_sizes(1:mxts))
		ALLOCATE(data_files(1:mxts))
		run_file = 'input/main_input'
		region_file = 'input/region_limits'
		OPEN(10, FILE=run_file)
		READ(10, NML=Input_Namelist)

		ALLOCATE(tsizes(ntsizes))	
		tsizes(1:ntsizes) = tile_sizes(1:ntsizes)
		DEALLOCATE(tile_sizes)
		CLOSE(10)
		max_tile_size = maxval(tsizes)/8.0

		ALLOCATE(lat_lims(2,n_regions))
		ALLOCATE(lon_lims(2,n_regions))

		OPEN(10, FILE=region_file)
		DO i = 1, n_regions
			READ(10,'(E12.4)')lon_lims(1,i)
			READ(10,'(E12.4)')lon_lims(2,i)
			READ(10,'(E12.4)')lat_lims(1,i)
			READ(10,'(E12.4)')lat_lims(2,i)
		ENDDO
		CLOSE(10)
		lat_lims0(:) = lat_lims(:,1)
		lon_lims0(:) = lon_lims(:,1)
	END SUBROUTINE Read_Run_Info
	
	SUBROUTINE Load_Dataset(region_id,ts_index)
		IMPLICIT NONE
		INTEGER*8 :: region_id, ts_index

		! The indice building routines need to know the 
		! solution grid for this region

		lats_pack(:) = lats_pack(:)+lat_lims(1,region_id+1)-lat_lims0(1)
		lats(:,:) = lats(:,:)+lat_lims(1,region_id+1)-lat_lims0(1)
		lons_pack(:) = lons_pack(:)+lon_lims(1,region_id+1)-lon_lims0(1)
		lons(:,:) = lons(:,:)+lon_lims(1,region_id+1)-lon_lims0(1)
		CALL Read_Raw_Data(region_id,ts_index)

		tsize_deg = tsizes(ts_index)/8.0
		CALL Build_Indices(ts_index)
		CALL Build_Shift(ts_index)
		CALL Generate_ATA_Indices(ts_index)

		! Revert in case we need to load another region
		lats_pack(:) = lats_pack(:)-lat_lims(1,region_id+1)+lat_lims0(1)
		lats(:,:) = lats(:,:)-lat_lims(1,region_id+1)+lat_lims0(1)
		lons_pack(:) = lons_pack(:)-lon_lims(1,region_id+1)+lon_lims0(1)
		lons(:,:) = lons(:,:)-lon_lims(1,region_id+1)+lon_lims0(1)

		DEALLOCATE(nvals)
		DEALLOCATE(kvals)
		DEALLOCATE(latvals)
		DEALLOCATE(lonvals)
		DEALLOCATE(lon)
		DEALLOCATE(lat)	
		IF (myid == 0) WRITE(*,'(A,I3,A,I4,A,I8,A)')" Data from region ",region_id," for tilesize ",ts_index,&
			" is loaded (ndata=", ndata(ts_index), ")."
	END SUBROUTINE Load_Dataset

	! load all of the data into the large array uxy
	! should do something similar to Read_Raw_Data
	! Added 04/2013 bgreer
	SUBROUTINE Read_All_Data()
		IMPLICIT NONE
		CHARACTER*120 :: dfile, indir, tstring
		INTEGER*8 :: maxndata
		INTEGER :: ii, ij, stat, ierr
		REAL*8 :: dummy
		REAL*8, ALLOCATABLE :: lonarr(:), latarr(:)
		REAL*8, ALLOCATABLE :: temparr(:)

		ALLOCATE(ndata(ntsizes))

		IF (myid == 0) THEN
			! find max ndata
			DO ii=1,ntsizes
				WRITE(tstring,'(I10)')tsizes(ii)
				dfile = 'input/data_'//TRIM(ADJUSTL(tstring))
				ndata(ii) = 0
				stat = 0
				OPEN(30,FILE=dfile)
				DO WHILE (stat == 0)
					READ(30,*,IOSTAT=stat), dummy
					ndata(ii) = ndata(ii) + 1
				ENDDO
				CLOSE(30)
				ndata(ii) = ndata(ii) - 1
			ENDDO
		ENDIF

		CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
		CALL MPI_Bcast(ndata, ntsizes, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)

		! allocate space for all data
		maxndata = MAXVAL(ndata)
		ALLOCATE(uxy(4,maxndata,ntsizes))
		ALLOCATE(nvals(maxndata,ntsizes))
		ALLOCATE(kvals(maxndata,ntsizes))
		ALLOCATE(lon_all(maxndata,ntsizes))
		ALLOCATE(lat_all(maxndata,ntsizes))
		ALLOCATE(nlon(ntsizes))
		ALLOCATE(nlat(ntsizes))
		ALLOCATE(n_n(ntsizes))
		ALLOCATE(n_k(ntsizes))

		! go back and read the data
		DO ii=1,ntsizes
			IF (myid == 0) THEN
				OPEN(30,FILE=data_files(ii))
				DO ij=1,ndata(ii)
					READ(30,*) dummy
					uxy(1,ij,ii) = dummy
					uxy(2,ij,ii) = dummy
					uxy(3,ij,ii) = dummy
					uxy(4,ij,ii) = dummy
					nvals(ij,ii) = dummy
					kvals(ij,ii) = dummy
					lat_all(ij,ii) = dummy
					lon_all(ij,ii) = dummy
				ENDDO
				CLOSE(30)
			ENDIF

			CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
			CALL MPI_Bcast(uxy, 4*maxndata*ntsizes, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
			CALL MPI_Bcast(nvals, maxndata*ntsizes, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
			CALL MPI_Bcast(kvals, maxndata*ntsizes, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
			CALL MPI_Bcast(lat_all, maxndata*ntsizes, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
			CALL MPI_Bcast(lon_all, maxndata*ntsizes, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

			! sort lon/lat arrays to get unique elements
			ALLOCATE(lonarr(ndata(ii)))
			ALLOCATE(latarr(ndata(ii)))
			ALLOCATE(temparr(ndata(ii)))
			temparr(:) = lon_all(1:ndata(ii),ii)
			lonarr = sort_array(temparr,ndata(ii))
			temparr(:) = lat_all(1:ndata(ii),ii)
			latarr = sort_array(temparr,ndata(ii))
			CALL find_unique(lonarr)
			CALL find_unique(latarr)
			nlon(ii) = SIZE(lonarr)
			nlat(ii) = SIZE(latarr)
			DEALLOCATE(temparr)
			DEALLOCATE(lonarr)
			DEALLOCATE(latarr)

			n_n(ii) = MAXVAL(nvals(:,ii))-MINVAL(nvals(:,ii))+1
			n_k(ii) = MAXVAL(kvals(:,ii))-MINVAL(kvals(:,ii))+1
			ALLOCATE(n(n_n(ii),ntsizes))
			ALLOCATE(ki(n_k(ii),ntsizes))
			ki(1,ii) = MINVAL(kvals(:,ii))
			DO ij = 2, n_k(ii)
				ki(ij,ii) = ki(ij-1,ii)+1
			ENDDO
			n(1,ii) = MINVAL(nvals(:,ii))
			DO ij = 2, n_n(ii)
				n(ij,ii) = n(ij-1,ii)+1
			ENDDO
		ENDDO

	END SUBROUTINE Read_All_Data
  
	SUBROUTINE Read_Raw_Data(r_id,t_id)
		IMPLICIT NONE
		REAL*8, ALLOCATABLE :: ex(:), ux(:), ey(:), uy(:)
		INTEGER*8 :: i, r_id, t_id
		CHARACTER*120 :: dfile, rstring,tstring, indir, file1
		! need to read ndata, nlon, and nlat
		! Otherwise will need to sort lonvals and latvals to get nlon, nlat, and lat and lon
		! lat and lon are the unique latitudes and longitudes for a given tile size

		! open 1 file first to get the dimensions

		WRITE(rstring,'(i10)')r_id
		WRITE(tstring,'(i10)')t_id
		indir = 'input/'

		IF (ALLOCATED(ex)) DEALLOCATE(ex)
		IF (ALLOCATED(ey)) DEALLOCATE(ey)
		IF (ALLOCATED(ux)) DEALLOCATE(ux)
		IF (ALLOCATED(uy)) DEALLOCATE(uy)
		IF (ALLOCATED(uxy)) DEALLOCATE(uxy)

		dfile = TRIM(ADJUSTL(indir))//'ex_'//TRIM(ADJUSTL(rstring))//'_'//TRIM(ADJUSTL(tstring))
		CALL Read_Ascii(ex,dfile)

		dfile = TRIM(ADJUSTL(indir))//'ey_'//TRIM(ADJUSTL(rstring))//'_'//TRIM(ADJUSTL(tstring))
		CALL Read_Ascii(ey,dfile)

		dfile = TRIM(ADJUSTL(indir))//'ux_'//TRIM(ADJUSTL(rstring))//'_'//TRIM(ADJUSTL(tstring))
		CALL Read_Ascii(ux,dfile)

		dfile = TRIM(ADJUSTL(indir))//'uy_'//TRIM(ADJUSTL(rstring))//'_'//TRIM(ADJUSTL(tstring))
		CALL Read_Ascii(uy,dfile)

!		ndata = SIZE(ux)
!		ALLOCATE(uxy(4,ndata))

!		uxy(1,:) = ux(:)
!		uxy(2,:) = ex(:)
!		uxy(3,:) = uy(:)
!		uxy(4,:) = ey(:)	

		DEALLOCATE(ex)
		DEALLOCATE(ey)
		DEALLOCATE(ux)
		DEALLOCATE(uy)

		IF (ALLOCATED(nvals)) DEALLOCATE(nvals)
		IF (ALLOCATED(kvals)) DEALLOCATE(kvals)
		IF (ALLOCATED(latvals)) DEALLOCATE(latvals)
		IF (ALLOCATED(lonvals)) DEALLOCATE(lonvals)
		IF (ALLOCATED(lat)) DEALLOCATE(lat)
		IF (ALLOCATED(lon)) DEALLOCATE(lon)

		dfile = TRIM(ADJUSTL(indir))//'n_'//TRIM(ADJUSTL(rstring))//'_'//TRIM(ADJUSTL(tstring))
!		CALL Read_Ascii_int(nvals,dfile)

		dfile = TRIM(ADJUSTL(indir))//'k_'//TRIM(ADJUSTL(rstring))//'_'//TRIM(ADJUSTL(tstring))
!		CALL Read_Ascii_int(kvals,dfile)

		dfile = TRIM(ADJUSTL(indir))//'lat_'//TRIM(ADJUSTL(rstring))//'_'//TRIM(ADJUSTL(tstring))
!		CALL Read_Ascii(latvals,dfile)

		dfile = TRIM(ADJUSTL(indir))//'lon_'//TRIM(ADJUSTL(rstring))//'_'//TRIM(ADJUSTL(tstring))
!		CALL Read_Ascii(lonvals,dfile)	

!		nlat = SIZE(latvals)
!		nlon = SIZE(lonvals)
	
!		ALLOCATE(lat(1:nlat))
!		ALLOCATE(lon(1:nlon))
!		lat = sort_array(latvals,nlat)
!		lon = sort_array(lonvals, nlon)
	
!		CALL find_unique(lat)
!		CALL find_unique(lon)

!		nlon = SIZE(lon)
!		nlat = SIZE(lat)

!		n_n = MAXVAL(nvals)-MINVAL(nvals)+1
!		n_k = MAXVAL(kvals)-MINVAL(kvals)+1
!		IF (ALLOCATED(n)) DEALLOCATE(n)
!		IF (ALLOCATED(ki)) DEALLOCATE(ki)
!		ALLOCATE(n(n_n))
!		ALLOCATE(ki(n_k))
!		ki(1) = MINVAL(kvals)
!		DO i = 2, n_k
!			ki(i) = ki(i-1)+1
!		ENDDO
!		n(1) = MINVAL(nvals)
!		DO i = 2, n_n
!			n(i) = n(i-1)+1
!		ENDDO
	END SUBROUTINE Read_Raw_Data



  Subroutine Build_Indices(ts_index)
	Implicit None
	Integer*8 :: i, ii, jj, this_n, this_k, found,nmin,kmin, ts_index
	Real*8 :: this_lat, this_lon
	Character*120 :: file1
	!  Builds nktindex, indicn, indic - to keep track of measurement modes and locations
	!  nktindex(n,k,ii,jj) = i if measurement i was made at n-1, k-1, lon(ii), lat(jj)
	!  If no measurements were made at that mode/location combination, then nktindex = 0 otherwise.
	!  indicn(n) = 1 if a measurement was made at n-1, and = 0 otherwise.
	!  indic(k,n) = 1 if a measurement was made at n-1, k-1, and = 0 otherwise.
	If (Allocated(nktindex)) DeAllocate(nktindex)
	If (Allocated(indicn)) DeAllocate(indicn)
	If (Allocated(indic)) DeAllocate(indic)

	
	! need to make this work for an arbitrary n_min (right now assumes 0 - which is fine) - same with k

	Allocate(nktindex(n_n(ts_index),n_k(ts_index),nlat(ts_index),nlon(ts_index)))   ! note that the indices are a little out of the sensible order here - this is correct
	Allocate(indicn(n_n(ts_index)))
	Allocate(indic(n_k(ts_index),n_n(ts_index)))
	nktindex(:,:,:,:) = int(-1)     
	indicn(:) = int(0)
	indic(:,:) = int(0)
!	Write(6,*), 'n_n is ',n_n, n(1)
!	Write(6,*), 'n_k is ', n_k, ki(1), maxval(ki)
	nmin = minval(n(:,ts_index))
	kmin = minval(ki(:,ts_index))
	Do i = 1, ndata(1)
		this_n = nvals(i,ts_index)-nmin
		this_k = kvals(i,ts_index)-kmin
		this_lat = latvals(i)
		this_lon = lonvals(i)
		found = 0
!		Do ii = 1, nxsz(1)
		Do ii = 1, nlon(ts_index)
		    if (found .eq. 0) then
!		    Do jj = 1, nysz(1)
		    Do jj = 1, nlat(ts_index)
			If ( (lon(ii) .eq. this_lon) .and. (lat(jj) .eq. this_lat) ) then
				found = 1
				nktindex(this_n+1,this_k+1,jj,ii) = i -1 ! see note in generate ATA indices for -1 explanation
				indicn(this_n+1) = 1
				indic(this_k+1,this_n+1) = 1
			endif
		    Enddo
		    endif
		Enddo
	Enddo
!	Write(6,*)'indic '
!  	do i = 1, 
!		write(6,'(6I9)')i, indic(i,:)
!	enddo
	file1 = 'nktindex.fits'
!	call writefits_4dl(nktindex,file1,n_n,n_k,nlat,nlon)
		
  End Subroutine Build_Indices
  
  Subroutine Build_Shift(ts_index)
	! builds the kxy array that describes how each measurement's kernel should be 
	! shifted into position onto the global grid.  kxy has been tested against
	! IDL and the output verified.  There is a precision issue in the area of the code
	! that uses int().   Hence the real*16 to be safe.
	Implicit None
	Integer ::i,j, xstart, xcen, ystart, ycen
	Integer*8 :: xoff, yoff, ts_index
	Real*16 :: lon_min, lat_min, dlon, dlat, this_lat, this_lon
	Character*120 :: file1 
	dlon = dlonlat
	dlat = dlonlat
	file1 = 'kxy.fits'
	If (Allocated(kxy)) DeAllocate(kxy)
	Allocate(kxy(2,nlon(ts_index),nlat(ts_index)))
	lon_min = lons(1,1)
	lat_min = lats(1,1)

	DO j = 1, nlat(ts_index)
	  DO i = 1, nlon(ts_index)
		xcen = lon(i)
		ycen = lat(j)
!		xstart = int(xcen/dlon)-(nxk-1)/2
!		ystart = int(ycen/dlat)-(nyk-1)/2
!		xstart = (nxk-1)/2
!		ystart = (nyk-1)/2
		this_lat = lat(j)
		this_lon = lon(i)
		xoff = int(((this_lon-lon_min)/dlon))  ! How far is the center of the measurement from the edge of the grid in pixels?
		yoff = int((this_lat-lat_min)/dlat)  ! Shift by that much, but
		kxy(1,i,j) = xoff-(nxk-1)/2        ! compensate for the size of the kernel tile to get the appropriate shift
		kxy(2,i,j) = yoff-(nyk-1)/2
	!	write(6,*)'kxy',kxy(, lon(i),lon_min
	  ENDDO
	ENDDO
!	call writefits_3dl(kxy,file1,2,nlon,nlat)
!	stop
!	Write(6,*)'------- Shifs -------'
!	DO j = 1, nlat
!	DO i = 1, nlon
!		Write(6,'(5E13.7)')real(kxy(1,i,1)), lon(i), lon_min,real(nxk), (lon(i)-lon_min)
!	Enddo
!	Enddo
	
!	stop
  ! This Array appears to be correct - NF 9/17/2010
  End Subroutine Build_Shift

  Subroutine Generate_ATA_Indices(ts_index)
		! Figures out the nonzero locations for kernels of different horizontal (tile) sizes at
		! each location on the solution grid.  Otherwise a lot of time can be wasted cross-correlating
		! zeros when the inversion matrix is assembled.

		! need 1D arrays for longitude and latitude from the input data - call them lon and lat
		Implicit None
		Integer*8 :: i,j,k, ioffset, nsolz, ii, iistart, iiend, ts_index
		Real*8 :: fudge, rad, radz, mxarr(2), lon_min, lon_max, lat_min, lat_max
		Character*120 :: file1, file2
		  If (Allocated(ata_indices)) DeAllocate(ata_indices)
		  If (Allocated(ata_num)) DeAllocate(ata_num)
		  If (Allocated(nzc_indices)) DeAllocate(nzc_indices)
	 	  If (Allocated(nzr_indices)) DeAllocate(nzr_indices)
		  If (Allocated(ilocal_indices)) DeAllocate(ilocal_indices)
		  If (Allocated(jlocal_indices)) DeAllocate(jlocal_indices)
		  If (Allocated(nznr_num)) DeAllocate(nznr_num)
		  If (Allocated(nznc_num)) DeAllocate(nznc_num)
!		  Allocate(ata_indices(1:nsolution,1:nxsz(1),1:nysz(1)))
		  Allocate(ata_indices(1:nsolution,nlon(ts_index),nlat(ts_index)))
		  Allocate(nzc_indices(1:nsolution,nlon(ts_index),nlat(ts_index)))
		  Allocate(nzr_indices(1:nsolution,nlon(ts_index),nlat(ts_index)))
		  Allocate(ilocal_indices(1:nsolution,nlon(ts_index),nlat(ts_index)))
		  Allocate(jlocal_indices(1:nsolution,nlon(ts_index),nlat(ts_index)))
		  Allocate(ata_num(nlon(ts_index),nlat(ts_index)))
		  Allocate(nznr_num(nlon(ts_index),nlat(ts_index)))
		  Allocate(nznc_num(nlon(ts_index),nlat(ts_index)))

		  ata_indices(:,:,:) = 0
		  ata_num(:,:) = 0
		fudge = 0.0  ! fudge factor since degraded resolution at depth can effectively smear
			     ! out small kernels.
		if (int(tsize_deg) .eq. 4) then 
			
			fudge = 1.0	
		endif
		if (int(tsize_deg) .eq. 2) then 
			fudge = 1.0
		endif
		rad = tsize_deg/2.0+fudge
		
		ioffset = 1
		mxarr(1) = rad
		Do k = 1, nzs
			nsolz = nxsz(k)*nysz(k)
			mxarr(2) = deg_factor(k)
			radz = maxval(mxarr)
			iistart = ioffset
			iiend = iistart+nsolz-1
			Do i = 1, nlon(ts_index)
				lon_min = lon(i)-radz
				lon_max = lon(i)+radz
				Do j = 1, nlat(ts_index)
					lat_min = lat(j)-radz
					lat_max = lat(j)+radz
					do ii = iistart, iiend
						if ( (lats_pack(ii) .ge. lat_min) .and. (lats_pack(ii) .le. lat_max) ) then
						if ( (lons_pack(ii) .ge. lon_min) .and. (lons_pack(ii) .le. lon_max) ) then
							ata_num(i,j) = ata_num(i,j)+1
							ata_indices(ata_num(i,j),i,j) = ii -1 ! NOTE: Adding -1 so that these indices are consistent with what would have			
						endif					      ! been generated by IDL (indices start at zero).  Augment_A and Add_stitch/Add_stitch2 use IDL-like indices
						endif					      ! and I would like to avoid modifying those routines.
					enddo
				Enddo
			Enddo
			ioffset = ioffset+nsolz
		Enddo		  
	!	WRITE(6,*)'ata num', tsize_deg
	!	do i = 1, nlat
	!		write(6,'(16I7)')i, ata_num(:,i)
	!	enddo
	!	file1 = 'ata_indices.fits'
	!	file2 = 'ata_num.fits'
	!	call writefits_2dl(ata_num,file2,nlon,nlat)
	!	call writefits_3dl(ata_indices,file1,nsolution,nlon,nlat)	
	!	stop
	!       ata_indices and ata_num have been verified against IDL's output
   End Subroutine Generate_ATA_Indices
   Subroutine test_sort()
	Real*8, allocatable :: sortme(:), sorted(:)
	integer*8 :: i, npts
	npts = 100
	allocate(sortme(npts))
	allocate(sorted(npts))
	do i = 1, npts 
		sortme(i) = i*((-1.0)**(i))	
!		write(6,*)sortme(i)
	enddo
	sortme(10) = sortme(1)
	sortme(20) = sortme(12)
!	call sort_array(sortme)
	sorted = sort_array(sortme,npts)
	do i = 1, npts
		write(6,*)sortme(i)
	enddo
	write(6,*)'******************************'
	do i = 1, npts
		write(6,*)sorted(i)
	enddo
!	deallocate(sortme)
!	deallocate(sorted)
	call find_unique(sorted)
	npts = size(sorted)
	write(6,*)'******************************'
	write(6,*)'unique number of points is ', npts
	do i = 1, npts+2
		write(6,*)sortme(i)
	enddo
	deallocate(sortme)
	deallocate(sorted) 
  end subroutine test_sort

   subroutine find_unique(allvals)
	Real*8, Allocatable :: allvals(:), valcopy(:)
	Integer*8 :: npts, ucount, i
	npts = size(allvals)
	allocate(valcopy(1:npts))
	ucount = 1
	valcopy(1) = allvals(1)
	do i = 2, npts
		if ( allvals(i) .gt. allvals(i-1) ) then
			ucount = ucount+1
			valcopy(ucount) = allvals(i)
		endif
	enddo
	deallocate(allvals)
	allocate(allvals(1:ucount))
	allvals(1:ucount) = valcopy(1:ucount)
	deallocate(valcopy)
   end subroutine find_unique

   Recursive Function Sort_Array(arr_unsorted,npts) Result(arr)
	Integer*8 :: npts
	Real*8 :: arr_unsorted(1:npts), tmp1, tmp2, tmp3, arr(1:npts)
	Real*8, Allocatable :: arrlt(:), arrgt(:), arreq(:), temp_arr(:), arrgtsort(:), arrltsort(:)
	Integer*8 :: ngt, nlt,neq, i
!	npts = size(arr_unsorted)
	allocate(arrlt(1:npts))
	allocate(arrgt(1:npts))
!	allocate(arr(1:npts))

	arr(1:npts) = arr_unsorted(1:npts)
!	write(6,*)arr(:)
	IF ( npts .eq. 1) Then
		write(6,*)arr(1:npts)
	Endif
	IF ( npts .eq. 2) Then
		If (arr(1) .gt. arr(2)) then
!			write(6,*)'swapping ', arr(1:2)
			tmp1 = arr(1)
			arr(1) = arr(2)
			arr(2) = tmp1
!			write(6,*)'swapped ', arr(1:2)
		Endif		
		
!		Return  ! arr
	Endif
	If (npts .gt. 2) Then
		tmp1 = minval(arr)
		tmp2 = maxval(arr)
		tmp3 = (tmp1+tmp2)/2.0
		ngt = 0
		nlt = 0
		neq = 0
		DO i = 1, npts
			if (arr(i) .gt.  tmp3) then
				ngt = ngt+1
				arrgt(ngt) = arr(i)
			endif
			if (arr(i) .lt. tmp3 ) then
				nlt = nlt+1
				arrlt(nlt) = arr(i)
			endif
!			if (arr(i) .eq. tmp3) then
!				neq = neq+1
!			endif
		ENDDO 
		neq = npts - nlt-ngt
!		write(6,*)neq, tmp3
		if (neq .eq. npts) then  ! array values are all the same - return
			return ! arr
		endif

		! next resize the arrays
		if (nlt .ge. 1) then
			allocate(temp_arr(nlt))
			temp_arr(1:nlt) = arrlt(1:nlt)
			deallocate(arrlt)
			allocate(arrlt(nlt))
			allocate(arrltsort(nlt))
			arrlt(:) = temp_arr(:)
			deallocate(temp_arr)
			if (nlt .gt. 1) then
			arrltsort=sort_array(arrlt,nlt)
			else
			arrltsort(1) = arrlt(1)
			endif
			arr(1:nlt) = arrltsort(1:nlt)
			deallocate(arrltsort)
		endif
		deallocate(arrlt)
		if (neq .ge. 1) then
			allocate(arreq(1:neq))
			arreq(1:neq) = tmp3
!			write(6,*)'neq is ', neq
!			write(6,*)'req ', arreq(1:neq)
			arr(nlt+1:nlt+1+neq-1) = arreq(1:neq)
			deallocate(arreq)
		endif	
		if (ngt .ge. 1) then
			allocate(temp_arr(ngt))
			temp_arr(1:ngt) = arrgt(1:ngt)
			deallocate(arrgt)
			allocate(arrgt(1:ngt))
			allocate(arrgtsort(1:ngt))
			arrgt(1:ngt) = temp_arr(1:ngt)
			deallocate(temp_arr)
			if (ngt .gt. 1) then
				arrgtsort = sort_array(arrgt,ngt)
			else
				arrgtsort(1) = arrgt(1)
			endif
			arr(nlt+neq+1:nlt+neq+1+ngt-1) = arrgtsort(1:ngt)
			deallocate(arrgtsort)
		endif
		deallocate(arrgt)
	!	return arr
	Endif
   End Function Sort_Array
End Module Data_Management
