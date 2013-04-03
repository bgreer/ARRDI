Module Region_Info
	! This module contains global variables and subroutines used to manage
	! the book-keeping associated with the overlapping regions.
	USE Solution_Grid
	Integer*8, Allocatable, Public :: neighbors(:,:), p_count(:,:), partner_points(:,:,:,:)
	Integer, Allocatable, Public :: have_match(:,:)
	Integer*8, Public :: n_regions, mxcp
	Real*8, Allocatable, Public :: lat_lims(:,:), lon_lims(:,:)
contains
	Subroutine Write_LonLat_Lims()
		Integer :: i
		do i = 1, n_regions
			Write(6,*)'-----------------'
			Write(6,*)'   Region ', i
			Write(6,*)lon_lims(1,i), lon_lims(2,i)
			Write(6,*)lat_lims(1,i), lat_lims(2,i)
		enddo	

	End Subroutine Write_lonlat_lims

	Subroutine Identify_Neighbors(lon_lims0,lat_lims0)
	  Real*8, Allocatable :: temp_lats(:), temp_lons(:)
	  Real*8, Allocatable :: mylons(:), mylats(:), lonsn(:), latsn(:)
	  Real*8 :: mymin, mymax, lmin, lmax
	  Real*8 :: lon_lims0(2), lat_lims0(2)
	  Integer*8 :: lat_true, lon_true, i,j

	  Allocate(neighbors(n_regions,n_regions))
	  Allocate(lonsn(nxsz(1)))
	  Allocate(latsn(nysz(1)))
	  Allocate(mylons(nxsz(1)))
          Allocate(mylats(nysz(1)))
          Allocate(temp_lons(nxsz(1)))
          Allocate(temp_lats(nysz(1)))


          latsn(:) = lats(:,1)
          lonsn(:) = lons(:,1)

          do i = 1, n_regions
                mylats(:) = latsn(:)-lat_lims0(1)+lat_lims(1,i)
                mylons(:) = lonsn(:)-lon_lims0(1)+lon_lims(1,i)
                do j = 1, n_regions
                        neighbors(i,j) = 0
                        if (i .ne. j) then
				!write(6,*)'test 1 ', size(temp_lats), size(latsn)
                                temp_lats(:) = latsn(:)-lat_lims0(1)+lat_lims(1,j)
                                temp_lons(:) = lonsn(:)-lon_lims0(1)+lon_lims(1,j)

                                lat_true = 0
                                lmin = minval(temp_lats)
                                lmax = maxval(temp_lats)
                                mymax = maxval(mylats)
                                mymin = minval(mylats)
                                if ( (lmin .le. mymax) .and. (lmin .ge. mymin) ) then
                                        lat_true = 1
                                endif
                                if ( (lmax .le. mymax) .and. (lmax .ge. mymin) ) then
                                        lat_true = 1
                                endif


                                lon_true = 0
                                lmin = minval(temp_lons)
                                lmax = maxval(temp_lons)
                                mymax = maxval(mylons)
                                mymin = minval(mylons)
                                if ( (lmin .le. mymax) .and. (lmin .ge. mymin) ) then
                                        lon_true = 1
                                endif
                                if ( (lmax .le. mymax) .and. (lmax .ge. mymin) ) then
                                        lon_true = 1
                                endif

                                if ( (lon_true .eq. 1) .and. (lat_true .eq. 1) ) then
                                        neighbors(i,j) = 1
                                endif
                        endif
                enddo
          enddo
	  DeAllocate(temp_lons)
	  DeAllocate(temp_lats)
	  DeAllocate(mylons)
	  DeAllocate(mylats)
	  DeAllocate(latsn)
	  DeAllocate(lonsn)
	  	


	 End Subroutine Identify_Neighbors

	 Subroutine find_overlap(lon_lims0,lat_lims0)
	  Implicit None
          Real*8, Allocatable :: temp_lats(:), temp_lons(:)
          Real*8, Allocatable :: mylon(:), mylat(:), myz(:), azs(:), alats(:), alons(:)
          Real*8 :: lon_min,lon_max,lat_min,lat_max, this_lon, this_lat, this_z
          Real*8 :: lon_lims0(2), lat_lims0(2), az, alat, alon
          Integer*8 :: lat_true, lon_true, i,j, ii, jj, ind
          Integer*8, Allocatable :: temp_arr(:,:,:,:)
	  
	  Allocate(mylat(nsolution))
	  Allocate(mylon(nsolution))	  
	  Allocate(myz(nsolution))	  
	  Allocate(alats(nsolution))
	  Allocate(alons(nsolution))	  
	  Allocate(azs(nsolution))	  

	
	  Allocate(partner_points(2,nsolution,n_regions,n_regions))
	  Allocate(p_count(nsolution,nsolution))
	  Allocate(have_match(nsolution,n_regions))
	  do i = 1, n_regions
			mylat(:) = 0.0D0
			mylon(:) = 0.0D0
			myz(:) = 0.0D0
			mylat(:) = lats_pack(:)-lat_lims0(1)+lat_lims(1,i)
			mylon(:) = lons_pack(:)-lon_lims0(1)+lon_lims(1,i)
			myz(:) = zs_pack(:)
			lat_max = maxval(mylat)
			lat_min = minval(mylat)
			lon_max = maxval(mylon)
			lon_min = minval(mylon)
			Write(6,*)'------ Region limits ------'
			Write(6,*)'   Region ', i
			Write(6,*)'lon : ', lon_min, lon_max
			Write(6,*)'lat : ', lat_min, lat_max
		!	if (i .eq. 1) then
		!		write(6,*)'testing'
		!		write(6,*)mylon(1:127)
		!		write(6,*)mylat(1:127)
		!	endif						
			do j = 1, n_regions
				ind = 1
				if (neighbors(i,j) .eq. 1) then
				  alats(:) = 0.0D0
				  alons(:) = 0.0D0
				  azs(:) = 0.0D0
				  alats(:) = lats_pack(:)-lat_lims0(1)+lat_lims(1,j)
				  alons(:) = lons_pack(:)-lon_lims0(1)+lon_lims(1,j)
				  azs(:) = zs_pack(:)
			if (j .eq. 40) then
				write(6,*)'testing'
				write(6,*)alats(3993:4158)
				write(6,*)alons(3993:4158)
				write(6,*)azs(3993:4158)
				write(6,*)zs_pack(3993:4158)
			endif		
  				  do ii = 1, nsolution

				    this_lat = mylat(ii)
				    this_lon = mylon(ii)
				    this_z = myz(ii)
				    do jj = 1, nsolution
					    alat = alats(jj)
					    alon = alons(jj)
					    az = azs(jj)
					     if (az .eq. this_z) then
					     if (alat .eq. this_lat) then
					     if (alon .eq. this_lon) then
						 ! ii from data set i corresponds to
						 ! jj from data set j
					         partner_points(1,ind,i,j) = ii-1   ! Start indices at zero to be consistent with arrays 
					   	 partner_points(2,ind,i,j) = jj-1   ! that were generated by idl with indices starting at zero
										    ! and which are then consistent with what matrix_magic routines expect. 
						have_match(ii,i) = 1
						have_match(jj,j) = 1
						 ! correlation the other way
					  	 partner_points(2,ind,j,i) = ii-1
					   	 partner_points(1,ind,j,i) = jj-1
						 ind = ind+1
						! if ( (i .eq. 1) .and. (j .eq. 4) ) then
						!	 if (az .eq. zs(1)) then
						!		write(6,*) ii, jj
						!	endif
						! endif
						! if ( (i .eq. 2) .and. (j .eq. 3) ) then
						!	if (az .eq. zs(1)) then
						!		write(6,*)this_lon,this_lat
						!	endif
						! endif

					     endif
					     endif
					     endif
					     
				    enddo
				  enddo
				  p_count(i,j) = ind-1
				  p_count(j,i) = ind-1
			!	  write(6,*)p_count(i,j)
				endif
				
			enddo
	  enddo
	
	  DeAllocate(mylat)
	  DeAllocate(mylon)	
	  DeAllocate(myz)	
	  DeAllocate(alats)
	  DeAllocate(alons)	
	  DeAllocate(azs)	
	
	  ! Resize the partner points array to save some memory
	  mxcp = maxval(p_count)
	  Allocate(temp_arr(2,mxcp,n_regions,n_regions))
	  temp_arr(:,:,:,:) = partner_points(1:2,1:mxcp,1:n_regions,1:n_regions)
	  DeAllocate(partner_points)
	  Allocate(partner_points(2,mxcp,n_regions,n_regions))
	  partner_points(:,:,:,:) = temp_arr(:,:,:,:)
	  DeAllocate(temp_arr)

	 End Subroutine find_overlap

End Module Region_Info
