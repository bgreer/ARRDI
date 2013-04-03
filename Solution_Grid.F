Module Solution_Grid
	Use Kernel_Management
	! This module contains global variables and subroutines related
	! to the setup of the solution (i.e. inversion) grid.
	! Contains Subroutines
	! I.   Get_Solution_Grid : general setup of solution grid
	!
	! II.  Generate_Weight_Arrays : Creates book-keeping arrays needed to manage the 
	!                               variable resolution with depth.
	!
	! III. Pick_Grid_Scheme :  Sets sdepths and the horizontals resolution at each depth
	!
	! IV.  Pack_Grid :  Takes the 3D arrays of lat,lon, and depth and stores them as 1D arrays.
	!
	Implicit None
	Integer :: nzs
	Integer, Allocatable :: deg_factor(:), nxsz(:), nysz(:)
	Integer*8, Allocatable, Public :: ata_indices(:,:,:), ata_num(:,:)
	Integer*8 :: nsolution, dims(3)
	Real*8, Allocatable :: weight_arr(:,:,:,:), weight_arro(:,:,:,:)
	Real*8, Allocatable :: zs(:), lats(:,:), lons(:,:), dzs(:), base_lat(:), base_lon(:)
	Real*8, Allocatable :: zs_pack(:), deltaz(:), lons_pack(:), lats_pack(:), dvolume(:)
	Real*8, Allocatable :: deltax(:), deltay(:)
	Real*8, Public :: dlonlat
	Integer*8, Public :: grid_scheme
	Integer*8, Allocatable :: ind_arro(:,:,:,:), num_arro(:,:,:), ind_arr(:,:,:,:), num_arr(:,:,:)
contains
	subroutine test_routine()
		Write(6,*)'I compiled'
	end subroutine test_routine
	Subroutine Get_Solution_Grid(lon_lims,lat_lims,max_tile_size)
	        Implicit None
	        Real*8 :: dlon, dlat, lat_lims(2), lon_lims(2)
		Real*8 :: lat_min, lat_max, lon_min, lon_max, icheck1
		Real*8 :: delta1, delta2, apod_factor,mtrad, max_tile_size, kerrad
		Integer*8 :: icheck2, max_tile_size_pix
		Integer*8 :: nxs0, nys0
		Integer*8 :: i,j
		dlon = dlonlat
		dlat = dlonlat
		Call Pick_Grid_Scheme()
		
		apod_factor = 7.5/8.0
		max_tile_size_pix = int(max_tile_size/dlonlat*apod_factor)
		mtrad = max_tile_size*apod_factor/2.0   ! radius of the apodized tile
		kerrad = (nxk-1)*dlonlat/2.0
		If (kerrad .lt. mtrad ) then
			Write(6,*)'Fatal Error : Kernel dimensions too small for largest tilesize'
			Write(6,*)'Check kernel set.'
			Write(6,*)'Largest tilesize (after apodization) is : ',max_tile_size_pix*dlonlat
			Write(6,*)'Kernel size is : ',(nxk-1)*dlonlat
			STOP
		Endif
			
		! max_tile_size_pix is the size of the *apodized* tile in (solution grid) pixels
		! For a solution grid with resolution of dlonlat = 7.5/8.0 degrees and the standard
		! apodization function (15 degrees in size), max_tile_size_pix will be 16 for a 16 degree tile.
!		Write(6,*)'max_tile_size ', max_tile_size
!		Write(6,*)'max_tile_size_pix ', max_tile_size_pix
!		Write(6,*)'kerrad is ',kerrad
		lat_min = lat_lims(1)-kerrad
		lat_max = lat_lims(2)+kerrad
		lon_min = lon_lims(1)-kerrad
		lon_max = lon_lims(2)+kerrad

		icheck1 = (lon_max-lon_min)/dlon
		icheck2 = int( (lon_max-lon_min)/dlon )
	
		If (icheck1 .ne. icheck2) then
			Write(6,*)'Range of longitudes must be integer multiple of ',dlon,' degrees.'
			Stop
		Endif
		icheck1 = (lat_max-lat_min)/dlat
		icheck2 = int( (lat_max-lat_min)/dlat )

		If (icheck1 .ne. icheck2) then
			Write(6,*)'Range of latitudes must be integer multiple of ',dlat,' degrees.'
			Stop
		Endif

	       ! note that in this notation, max_tile_size_pix is really the maximum tile size
	       ! divided by dlon or dlat - so a tile apodized to 15 deg at a resolution
	       ! of .9375 deg (7/8) would give max_tile_size_pix of 16.
	       nxs0=(lon_max-lon_min)/dlon + 1 !max_tile_size_pix + 1   ! Number of solution grid points at
               nys0=(lat_max-lat_min)/dlat + 1 ! max_tile_size_pix + 1   ! the uppermost depth
	       Allocate(nxsz(nzs))  ! # of solution grid points at each depth
	       Allocate(nysz(nzs))
!	write(6,*)'-----min max ------------'
!	write(6,*)lat_min, lat_max
!	write(6,*)lon_min, lon_max
        DO i = 1, nzs 
                nxsz(i) = (nxs0-1)/deg_factor(i)+1
                nysz(i) = (nys0-1)/deg_factor(i)+1
!		write(6,*)nxsz(i), nysz(i)
        ENDDO
!	Write(6,*)'------------'

	Allocate(lats(nxs0,nzs))
	Allocate(lons(nxs0,nzs))
	
	Do j = 1, nzs
		Do i = 1, nxsz(j)
			lons(i,j) = lon_min+(i-1)*dlon*deg_factor(j)
		EndDo
		Do i = 1, nysz(j)
			lats(i,j) = lat_min+(i-1)*dlat*deg_factor(j)
		EndDo
!		Write(6,*)'************** lons *************'
!		Write(6,*)lons(:,j)
!		Write(6,*)'************** lats *************'
!		Write(6,*)lats(:,j)
		
	EndDo

	Allocate(dzs(nzs))
	Do i = 1, nzs
		If ( (i .ne. 1) .and. (i .ne. nzs)) Then
			delta1 = zs(i)-zs(i-1)
			delta2 = zs(i+1)-zs(i)
			dzs(i) = 0.5*(delta1+delta2)
		Endif
		If (i .eq. 1) Then
			dzs(i) = 0.5*(zs(i+1)-zs(i))
		Endif
		If (i .eq. nzs) Then
			dzs(i) = 0.5*(zs(i)-zs(i-1))
		Endif
	EndDo

	dims(1) = nxs0
	dims(2) = nys0
	dims(3) = nzs

	nsolution = 0

	Do i = 1, nzs
		nsolution = nsolution+nxsz(i)*nysz(i)
	EndDo
	write(6,*)'Matrix is NxN with N = ', nsolution
	Allocate(base_lat(nys0))
	Allocate(base_lon(nxs0))
	base_lat(:) = lats(:,1)
	base_lon(:) = lons(:,1)
	Call Generate_Weight_Arrays()
	Allocate(zs_pack(nsolution))
	Allocate(lats_pack(nsolution))
	Allocate(lons_pack(nsolution))
	Allocate(deltaz(nsolution))
	Allocate(dvolume(nsolution))
	Call Build_Zweights()
	End Subroutine Get_Solution_Grid
	
	Subroutine Generate_Weight_Arrays
		Implicit None
		Integer*8 :: nxs0, nys0, mwtotal, icount
		Real*8, Allocatable :: wtotal(:,:,:), lats0(:), lons0(:), lon_temp(:), lat_temp(:)
		Real*8, Allocatable :: i_weight_arr(:,:,:,:,:)
		Real*8 :: dlon, dlat, dlonz, dlatz, lat_maxz, lat_minz
		Real*8 :: divtotal, wtmax, lon_minz, lon_maxz
		Integer :: i, j, k, ny_z, nx_z, ii, jj
		Integer*8, Allocatable :: num_test(:,:,:)
		Character*120 :: test_file1, test_file2, test_file3
		nxs0 = dims(1)
		nys0 = dims(2)

		Allocate(i_weight_arr(nxs0,nys0,nxs0,nys0,nzs))
		
		Allocate(wtotal(nxs0,nys0,nzs))	

		Allocate(num_test(nxs0,nys0,nzs))
		i_weight_arr = 0.0D0
		wtotal = 0.0D0
	
		dlon = lons(2,1)-lons(1,1)	
		dlat = lats(2,1)-lats(1,1)	
		dlon = dlon/deg_factor(1)
		dlat = dlat/deg_factor(1)

		Allocate(lats0(nys0))
		Allocate(lons0(nxs0))

		lats0(:) = base_lat(:)
		lons0(:) = base_lon(:)

		Allocate(lon_temp(nxs0))
		Allocate(lat_temp(nys0))
		do k = 1, nzs
			dlonz = dlon*deg_factor(k)
			dlatz = dlat*deg_factor(k)
			nx_z = nxsz(k)
			ny_z = nysz(k)
			lon_temp(:) = lons(:,k)
			lat_temp(:) = lats(:,k)
	!		write(6,*)'dlatz ', dlatz
	!		write(6,*)'dlonz ', dlonz
			Do j = 1, ny_z
				lat_minz = lat_temp(j)-dlatz/2.0
				lat_maxz = lat_temp(j)+dlatz/2.0
				IF ( ny_z .eq. 1) then
					lat_minz = lats0(1)
					lat_maxz = lats0(nys0)
				endif
				Do i = 1, nx_z
					lon_minz = lon_temp(i)-dlonz/2.0
					lon_maxz = lon_temp(i)+dlonz/2.0
					IF ( nx_z .eq. 1) then
						lon_minz = lons0(1)
						lon_maxz = lons0(nxs0)
					endif
					icount = 0
					do jj = 1, nys0
						if ( (lats0(jj) .ge. lat_minz) .and. (lats0(jj) .le. lat_maxz)) then
						do ii = 1, nxs0
						if ( (lons0(ii) .ge. lon_minz) .and. (lons0(ii) .le. lon_maxz)) then
								i_weight_arr(ii,jj,i,j,k) = 1.0	
								icount = icount+1
						endif
						enddo
						endif

					enddo
					num_test(i,j,k) = icount
					wtotal(i,j,k) = SUM(i_weight_arr(:,:,i,j,k))
				Enddo
			Enddo
			do j = 1, nys0
				do i = 1, nxs0
					divtotal = sum(i_weight_arr(i,j,:,:,k))
					i_weight_arr(i,j,:,:,k) = i_weight_arr(i,j,:,:,k)/divtotal
				enddo
			enddo			
		enddo
		test_file1 = 'i_weight_arr.fits'
!		call writefits_5dd(i_weight_arr,test_file1,nxs0,nys0,nxs0,nys0,nzs)
		wtmax = maxval(wtotal)
		mwtotal = int(wtmax)

		Allocate(ind_arr(mwtotal,nxs0,nys0,nzs))
		Allocate(num_arr(nxs0,nys0,nzs))
		Allocate(weight_arr(mwtotal,nxs0,nys0,nzs))
		do k = 1, nzs
			do j = 1, nys0
				do i = 1, nxs0
				  icount = 0
				  do jj = 1, nys0
				  do ii = 1, nxs0
					if (i_weight_arr(ii,jj,i,j,k) .gt. 0.0 ) then
						
						icount = icount+1
						ind_arr(icount,i,j,k) = ii+(jj-1)*nxs0-1  ! had ii+jj*(ii-1) here, but I think this is wrong.  we need an idl-like index on a 2d slice - hence the -1
					 	weight_arr(icount,i,j,k) = i_weight_arr(ii,jj,i,j,k)
					endif		
				  enddo
				  enddo
				  num_arr(i,j,k) = icount
				
				enddo
			enddo
		enddo



		Allocate(ind_arro(mwtotal,nxs0,nys0,nzs))
		Allocate(num_arro(nxs0,nys0,nzs))
		Allocate(weight_arro(mwtotal,nxs0,nys0,nzs))
		do k = 1, nzs
			do j = 1, nys0
				do i = 1, nxs0
				  icount = 0
				  do jj = 1, nys0
				  do ii = 1, nxs0
					if (i_weight_arr(i,j,ii,jj,k) .gt. 0.0 ) then
						
						icount = icount+1
						ind_arro(icount,i,j,k) = ii+(jj-1)*nxs0-1
					 	weight_arro(icount,i,j,k) = i_weight_arr(i,j,ii,jj,k)
					endif		
				  enddo
				  enddo
				  num_arro(i,j,k) = icount
				enddo
			enddo
		enddo
		DeAllocate(wtotal)
		! looks like in IDL, I had actually used ind_arro - not ind_arr
		ind_arr(:,:,:,:) = 0
		weight_arr(:,:,:,:) = 0.0D0
		num_arr(:,:,:) = 0
		num_arr = num_arro
		ind_arr = ind_arro
		weight_arr = weight_arro
		test_file1 = 'ind_arr.fits'
		test_file2 = 'num_arr.fits'
		test_file3 = 'weight_arr.fits'
!		call writefits_4dl(ind_arr,test_file1,mwtotal,nxs0,nys0,nzs)
!		call writefits_3dl(num_arr,test_file2,nxs0,nys0,nzs)
!		call writefits_4dd(weight_arr,test_file3,mwtotal,nxs0,nys0,nzs)
!		stop
!	        write_fits_3dl(ind_arr,test_file1,mwtotal,nxs0,nys0,nzs)
		! weight_arr, ind_arr, and num_arr have been checked against IDL's output and are correct.
	End Subroutine Generate_Weight_Arrays


	Subroutine Pick_Grid_Scheme()
		Implicit None
!		Integer*8 :: grid_scheme
		select case (grid_scheme)            
     			case (43)
			   write(6,*)'using grid scheme 43'
			   nzs = 19
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.25 
	                   zs(3) = 0.5 
	                   zs(4) = 1.0 
	                   zs(5) = 1.5 
	                   zs(6) = 2.0 
	                   zs(7) = 3.0 
	                   zs(8) = 4.0 
	                   zs(9) = 5.0 
	                   zs(10) = 6.0 
	                   zs(11) = 7.0 
	                   zs(12) = 9.0 
	                   zs(13) = 11.0 
	                   zs(14) = 13.0 
	                   zs(15) = 15.0 
	                   zs(16) = 17.0 
	                   zs(17) = 19.0 
	                   zs(18) = 22.0 
	                   zs(19) = 25.0 
		
			   deg_factor(1) = 1
			   deg_factor(2:11) = 2
			   deg_factor(12:15) = 4
			   deg_factor(16:19) = 4
     			case (42)
			   write(6,*)'using grid scheme 42  - this is most suited for a dlonlat of 1.25 deg'
			   nzs = 19
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.25 
	                   zs(3) = 0.5 
	                   zs(4) = 1.0 
	                   zs(5) = 1.5 
	                   zs(6) = 2.0 
	                   zs(7) = 3.0 
	                   zs(8) = 4.0 
	                   zs(9) = 5.0 
	                   zs(10) = 6.0 
	                   zs(11) = 7.0 
	                   zs(12) = 9.0 
	                   zs(13) = 11.0 
	                   zs(14) = 13.0 
	                   zs(15) = 15.0 
	                   zs(16) = 17.0 
	                   zs(17) = 19.0 
	                   zs(18) = 22.0 
	                   zs(19) = 25.0 
		
			   deg_factor(1:8) = 1
			   deg_factor(9:14) = 2
			   !deg_factor(9:11) = 1
			   deg_factor(15:18) = 3
			   deg_factor(19:19) = 3
     			case (421)
			   write(6,*)'using grid scheme 421  - this is most suited for a dlonlat of 0.9375 deg'
			   nzs = 19
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.25 
	                   zs(3) = 0.5 
	                   zs(4) = 1.0 
	                   zs(5) = 1.5 
	                   zs(6) = 2.0 
	                   zs(7) = 3.0 
	                   zs(8) = 4.0 
	                   zs(9) = 5.0 
	                   zs(10) = 6.0 
	                   zs(11) = 7.0 
	                   zs(12) = 9.0 
	                   zs(13) = 11.0 
	                   zs(14) = 13.0 
	                   zs(15) = 15.0 
	                   zs(16) = 17.0 
	                   zs(17) = 19.0 
	                   zs(18) = 22.0 
	                   zs(19) = 25.0 
		
			   deg_factor(1:8) = 1
			   deg_factor(9:14) = 2
			   !deg_factor(9:11) = 1
			   deg_factor(15:18) = 4
			   deg_factor(19:19) = 4
				case (422)
			   write(6,*)'using grid scheme 422  - same as 421 but goes deeper'
			   nzs = 24
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.25 
	                   zs(3) = 0.5 
	                   zs(4) = 1.0 
	                   zs(5) = 1.5 
	                   zs(6) = 2.0 
	                   zs(7) = 3.0 
	                   zs(8) = 4.0 
	                   zs(9) = 5.0 
	                   zs(10) = 6.0 
	                   zs(11) = 7.0 
	                   zs(12) = 9.0 
	                   zs(13) = 11.0 
	                   zs(14) = 13.0 
	                   zs(15) = 15.0 
	                   zs(16) = 17.0 
	                   zs(17) = 19.0 
	                   zs(18) = 22.0 
	                   zs(19) = 25.0 
					   zs(20) = 30.0
					   zs(21) = 35.0
					   zs(22) = 40.0
					   zs(23) = 45.0
					   zs(24) = 50.0
			   deg_factor(1:8) = 1
			   deg_factor(9:14) = 2
			   deg_factor(15:24) = 4
				case (4000)
			   write(6,*)'using grid scheme 4000 - 421 cut off at 10Mm'
			   nzs = 17
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.15 
	                   zs(3) = 0.25 
	                   zs(4) = 0.5 
	                   zs(5) = 0.75 
	                   zs(6) = 1.0 
	                   zs(7) = 1.25
	                   zs(8) = 1.5
	                   zs(9) = 2.0 
	                   zs(10) = 3.0 
	                   zs(11) = 4.0 
	                   zs(12) = 5.0 
	                   zs(13) = 6.0 
	                   zs(14) = 7.0 
	                   zs(15) = 8.0
					   zs(16) = 9.0
	                   zs(17) = 10.0 
			   deg_factor(1:11) = 1
			   deg_factor(12:17) = 2
				case (9999)
			   write(6,*)'using grid scheme 9999 i give up describing these'
			   nzs = 29
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
					   zs(2) = 0.15
	                   zs(3) = 0.25
					   zs(4) = 0.50
	                   zs(5) = 0.75 
	                   zs(6) = 1.0 
					   zs(7) = 1.25
	                   zs(8) = 1.5 
	                   zs(9) = 2.0 
	                   zs(10) = 3.0 
					   zs(11) = 3.5
	                   zs(12) = 4.0 
					   zs(13) = 4.5
	                   zs(14) = 5.0 
	                   zs(15) = 6.0 
	                   zs(16) = 7.0 
	                   zs(17) = 9.0 
	                   zs(18) = 11.0 
	                   zs(19) = 13.0 
	                   zs(20) = 15.0 
	                   zs(21) = 17.0 
	                   zs(22) = 19.0 
	                   zs(23) = 22.0 
	                   zs(24) = 25.0 
					   zs(25) = 30.0
					   zs(26) = 35.0
					   zs(27) = 40.0
					   zs(28) = 45.0
					   zs(29) = 50.0
			   deg_factor(1:1) = 1
			   deg_factor(2:19) = 2
			   deg_factor(20:29) = 4
				case (9998)
			   write(6,*)'using grid scheme 9998 i give up describing these'
			   nzs = 48
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
					   zs(2) = 0.15
	                   zs(3) = 0.25
					   zs(4) = 0.50
	                   zs(5) = 0.75 
	                   zs(6) = 1.0 
					   zs(7) = 1.25
	                   zs(8) = 1.5 
					   zs(9) = 1.75
	                   zs(10) = 2.0
					   zs(11) = 2.25
					   zs(12) = 2.5
					   zs(13) = 2.75
	                   zs(14) = 3.0
					   zs(15) = 3.25
					   zs(16) = 3.5
					   zs(17) = 3.75
	                   zs(18) = 4.0
					   zs(19) = 4.25
					   zs(20) = 4.5
					   zs(21) = 4.75
	                   zs(22) = 5.0
					   zs(23) = 5.5
	                   zs(24) = 6.0 
	                   zs(25) = 7.0 
	                   zs(26) = 9.0 
	                   zs(27) = 11.0 
	                   zs(28) = 13.0 
	                   zs(29) = 15.0 
	                   zs(30) = 17.0 
	                   zs(31) = 19.0 
	                   zs(32) = 21.0 
					   zs(33) = 22.5
	                   zs(34) = 25.0 
					   zs(35) = 27.5
					   zs(36) = 30.0
					   zs(37) = 32.5
					   zs(38) = 35.0
					   zs(39) = 37.5
					   zs(40) = 40.0
					   zs(41) = 42.5
					   zs(42) = 45.0
					   zs(43) = 50.0
					   zs(44) = 55.0
					   zs(45) = 60.0
					   zs(46) = 65.0
					   zs(47) = 70.0
					   zs(48) = 75.0
			   deg_factor(1:1) = 1
			   deg_factor(2:2) = 2
			   deg_factor(3:48) = 4
				case (9997)
			   write(6,*)'using grid scheme 9997 i give up describing these'
			   nzs = 43
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
					   zs(2) = 0.15
	                   zs(3) = 0.25
					   zs(4) = 0.50
	                   zs(5) = 0.75 
	                   zs(6) = 1.0 
					   zs(7) = 1.25
	                   zs(8) = 1.5 
					   zs(9) = 1.75
	                   zs(10) = 2.0
					   zs(11) = 2.25
					   zs(12) = 2.5
					   zs(13) = 2.75
	                   zs(14) = 3.0
					   zs(15) = 3.25
					   zs(16) = 3.5
					   zs(17) = 3.75
	                   zs(18) = 4.0
					   zs(19) = 4.25
					   zs(20) = 4.5
					   zs(21) = 4.75
	                   zs(22) = 5.0
					   zs(23) = 5.5
	                   zs(24) = 6.0 
	                   zs(25) = 7.0 
	                   zs(26) = 9.0 
	                   zs(27) = 11.0 
	                   zs(28) = 13.0 
	                   zs(29) = 15.0 
	                   zs(30) = 17.0 
	                   zs(31) = 19.0 
	                   zs(32) = 22.0 
	                   zs(33) = 25.0 
					   zs(34) = 30.0
					   zs(35) = 35.0
					   zs(36) = 40.0
					   zs(37) = 45.0
					   zs(38) = 50.0
					   zs(39) = 55.0
					   zs(40) = 60.0
					   zs(41) = 65.0
					   zs(42) = 70.0
					   zs(43) = 75.0
			   deg_factor(1:1) = 1
			   deg_factor(2:2) = 2
			   deg_factor(3:43) = 2
     			case (9001)
			   write(6,*)'using grid scheme 9001  - it is pretty great'
			   nzs = 24
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.25 
	                   zs(3) = 0.5 
	                   zs(4) = 1.0 
	                   zs(5) = 1.5 
	                   zs(6) = 2.0 
	                   zs(7) = 3.0 
	                   zs(8) = 4.0 
	                   zs(9) = 5.0 
	                   zs(10) = 6.0 
	                   zs(11) = 7.0 
	                   zs(12) = 9.0 
	                   zs(13) = 11.0 
	                   zs(14) = 13.0 
	                   zs(15) = 15.0 
	                   zs(16) = 17.0 
	                   zs(17) = 19.0 
	                   zs(18) = 22.0 
	                   zs(19) = 25.0 
					   zs(20) = 30.0
					   zs(21) = 35.0
					   zs(22) = 40.0
					   zs(23) = 45.0
					   zs(24) = 50.0
			   deg_factor(1:1) = 1
			   deg_factor(2:14) = 2
			   deg_factor(15:24) = 4
				case (9002)
			   write(6,*)'using grid scheme 9002  - it is pretty great'
			   nzs = 8
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 1.0 
	                   zs(3) = 4.0 
	                   zs(4) = 9.0 
	                   zs(5) = 15.0 
	                   zs(6) = 20.0 
	                   zs(7) = 30.0 
	                   zs(8) = 40.0 
			   deg_factor(1:1) = 1
			   deg_factor(2:4) = 2
			   deg_factor(5:8) = 4
				case (9003)
			   write(6,*)'using grid scheme 9003  - it is pretty great'
			   nzs = 8
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 1.0 
	                   zs(3) = 4.0 
	                   zs(4) = 9.0 
	                   zs(5) = 15.0 
	                   zs(6) = 20.0 
	                   zs(7) = 30.0 
	                   zs(8) = 40.0 
			   deg_factor(1:1) = 1
			   deg_factor(2:8) = 2
			case (4359)
			   write(6,*)'using grid scheme 4359'
			nzs = 15
			Allocate(zs(nzs))
			Allocate(deg_factor(nzs))
			zs(1) = 0.0
			zs(2) = 0.5
			zs(3) = 1.0
			zs(4) = 2.0
			zs(5) = 3.0
			zs(6) = 4.0
			zs(7) = 5.0
			zs(8) = 6.0
			zs(9) = 9.0
			zs(10) = 11.0
			zs(11) = 13.0
			zs(12) = 17.0
			zs(13) = 19.0
			zs(14) = 22.0
			zs(15) = 25.0
			deg_factor(1) = 1
			deg_factor(2:15) = 2

     			case (41)
			   write(6,*)'using grid scheme 41'
			   nzs = 19
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.25 
	                   zs(3) = 0.5 
	                   zs(4) = 1.0 
	                   zs(5) = 1.5 
	                   zs(6) = 2.0 
	                   zs(7) = 3.0 
	                   zs(8) = 4.0 
	                   zs(9) = 5.0 
	                   zs(10) = 6.0 
	                   zs(11) = 7.0 
	                   zs(12) = 9.0 
	                   zs(13) = 11.0 
	                   zs(14) = 13.0 
	                   zs(15) = 15.0 
	                   zs(16) = 17.0 
	                   zs(17) = 19.0 
	                   zs(18) = 22.0 
	                   zs(19) = 25.0 
		
			   deg_factor(1) = 1
			   deg_factor(2:11) = 2
			   deg_factor(12:19) = 3
     			case (39)
			   write(6,*)'using grid scheme 39'
			   nzs = 18
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.5 
	                   zs(3) = 1.0 
	                   zs(4) = 1.5 
	                   zs(5) = 2.0 
	                   zs(6) = 3.0 
	                   zs(7) = 4.0 
	                   zs(8) = 5.0 
	                   zs(9) = 6.0 
	                   zs(10) = 7.0 
	                   zs(11) = 9.0 
	                   zs(12) = 11.0 
	                   zs(13) = 13.0 
	                   zs(14) = 15.0 
	                   zs(15) = 17.0 
	                   zs(16) = 19.0 
	                   zs(17) = 22.0 
	                   zs(18) = 25.0 
		
			   deg_factor(1:3) = 1
			   deg_factor(4:11) = 2
			   deg_factor(12:18) = 3
     			case (38)
			   write(6,*)'using grid scheme 38'
			   nzs = 18
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.5 
	                   zs(3) = 1.0 
	                   zs(4) = 1.5 
	                   zs(5) = 2.0 
	                   zs(6) = 3.0 
	                   zs(7) = 4.0 
	                   zs(8) = 5.0 
	                   zs(9) = 6.0 
	                   zs(10) = 7.0 
	                   zs(11) = 9.0 
	                   zs(12) = 11.0 
	                   zs(13) = 13.0 
	                   zs(14) = 15.0 
	                   zs(15) = 17.0 
	                   zs(16) = 19.0 
	                   zs(17) = 22.0 
	                   zs(18) = 25.0 
		
			   deg_factor(1:8) = 1
			   deg_factor(9:11) = 2
			   deg_factor(12:18) = 3
     			case (37)
			   write(6,*)'using grid scheme 37'
			   nzs = 18
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.5 
	                   zs(3) = 1.0 
	                   zs(4) = 1.5 
	                   zs(5) = 2.0 
	                   zs(6) = 3.0 
	                   zs(7) = 4.0 
	                   zs(8) = 5.0 
	                   zs(9) = 6.0 
	                   zs(10) = 7.0 
	                   zs(11) = 9.0 
	                   zs(12) = 11.0 
	                   zs(13) = 13.0 
	                   zs(14) = 15.0 
	                   zs(15) = 17.0 
	                   zs(16) = 19.0 
	                   zs(17) = 22.0 
	                   zs(18) = 25.0 
		
			   deg_factor(1:6) = 1
			   deg_factor(7:13) = 2
			   deg_factor(14:18) = 3

     			case (371)
			   write(6,*)'using grid scheme 371'
			   nzs = 26
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 0.5 
	                   zs(3) = 1.0 
	                   zs(4) = 1.5 
	                   zs(5) = 2.0
			   zs(6) = 2.5 
	                   zs(7) = 3.0
			   zs(8) = 3.5 
	                   zs(9) = 4.0
			   zs(10) = 4.5 
	                   zs(11) = 5.0 
			   zs(12) = 5.5
	                   zs(13) = 6.0 
			   zs(14) = 6.5
	                   zs(15) = 7.0 
			   zs(16) = 8.0
	                   zs(17) = 9.0 
			   zs(18) = 10.0
	                   zs(19) = 11.0
			   zs(20) = 12.0 
	                   zs(21) = 13.0 
	                   zs(22) = 15.0 
	                   zs(23) = 17.0 
	                   zs(24) = 19.0 
	                   zs(25) = 22.0 
	                   zs(26) = 25.0 
		
			   deg_factor(1:7) = 1
			   deg_factor(8:18) = 2
			   deg_factor(19:26) = 3

     			case (89)
			   write(6,*)'using grid scheme 89'
			   nzs = 6
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
	                   zs(1) = 0.0 
	                   zs(2) = 1.25 
	                   zs(3) = 2.5 
	                   zs(4) = 10.0 
	                   zs(5) = 15.0 
	                   zs(6) = 17.0 
		
			   deg_factor(1) = 1
			   deg_factor(2:6) = 3
     			case default
			   write(6,*)'using default scheme'
        		   nzs = 6
			   Allocate(zs(nzs))
			   Allocate(deg_factor(nzs))
			   zs(1) = 0.0
			   zs(2) = 3.0
			   zs(3) = 5.0
			   zs(4) = 6.0
			   zs(5) = 8.0
			   zs(6) = 14.0
			   deg_factor(1) = 1
			   deg_factor(2) = 2
			   deg_factor(3) = 2
			   deg_factor(4) = 4
			   deg_factor(5) = 4
			   deg_factor(6) = 4
  			end select  
		!	Write(6,*)zs
		!	Write(6,*)deg_factor

	End Subroutine Pick_Grid_Scheme

	Subroutine Pack_Grid()
		Implicit None
		Integer :: iz, iy,ix,ii
		Integer*8 :: ind, nxy, indz
		Real*8 :: this_dxy
		ind = 1
		indz = 1
		do iz = 1, nzs
		  nxy = nxsz(iz)*nysz(iz)
		  do iy = 1, nysz(iz)
		
			lats_pack(ind:ind+nxsz(iz)-1) = lats(iy,iz)
			do ii = 1, nxsz(iz)
				lons_pack(ind+ii-1) = lons(ii,iz)
			enddo
			ind = ind+nxsz(iz)
		  enddo
		  zs_pack(indz:indz+nxy-1) = zs(iz)
		  deltaz(indz:indz+nxy-1) = dzs(iz)
		  indz = indz+nxy
		enddo
		ind = 1
		do iz = 1, nzs
		this_dxy = (lons_pack(2)-lons_pack(1))*deg_factor(iz)
		do iy = 1, nysz(iz)
		do ii = 1, nxsz(iz)
			dvolume(ind) = this_dxy*dzs(iz)
			ind = ind+1
		enddo
		enddo
		enddo

	End Subroutine Pack_grid

	Subroutine Build_Zweights()
		Implicit None
		! verified against IDL output - NF 9/20/2010
		Integer :: i, j
		Character*120 :: file1, file2
		Allocate(z_owner(nzs,nzk))
		Allocate(wz(nzk))
		z_owner(:,:) = 0
		wz(:) = 0.0D0
		! compute wz  - really dz
		! it's the dz of the kernel, not the target grid
		wz(1) = (z_kernel(2)-z_kernel(1))/2.0
		wz(nzk) = (z_kernel(nzk)-z_kernel(nzk-1))/2.0
		Do i = 2, nzk-1
			wz(i) = (z_kernel(i+1)-z_kernel(i-1))/2.0 
		Enddo
		
	
		! compute z_owner
		! z_owner is either 0 or 1, used as a mask later
		! 0 by default
		! 1 if kernel depth j is contained within target depths i and i+1
		Do i = 1, nzs-1 
			Do j = 1, nzk
				If  ( (z_kernel(j) .ge. zs(i)) .and. (z_kernel(j) .lt. zs(i+1)) ) Then
					z_owner(i,j) = 1
				Endif
			Enddo
		Enddo
		Do j = 1, nzk
			If (z_kernel(j) .ge. zs(nzs)) then
				z_owner(nzs,j) = 1
			Endif
		Enddo	
		DeAllocate(z_kernel)
!		Write(6,*)'z owner)'
!		do j = 1, nzk
!			write(6,*)wz(j)
!		enddo
	!	file1 = 'z_owner.fits'
	!	file2 = 'wz.fits'
	!	call writefits_2dl(z_owner,file1,nzs,nzk)
	!	call writefits_1dd(wz,file2,nzk)
	End Subroutine Build_Zweights


End Module Solution_Grid
