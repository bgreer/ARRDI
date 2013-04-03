Module Kernel_Management
	Use IO_Routines
	Use FitsLib_NF
	Implicit None
        Integer*8, Public :: nxk, nyk, nzk, nkers, kernel_set
	Integer*8, Allocatable, Public :: z_owner(:,:), nkindex(:,:)
	Real*8, Allocatable, Public :: wz(:), smeared_kers(:,:,:,:), z_kernel(:)
	Character*120 :: ker_dir

  Contains
	Subroutine Initialize_Kernels()
		Implicit None
		Character*120 :: ker_dim_file, ker_z_file, ker_set_string
		Write(ker_set_string,'(i10)')kernel_set
		ker_dir = '/polaris/kernel_sets/'//TRIM(ADJUSTL(ker_set_string))
		ker_z_file = TRIM(ADJUSTL(ker_dir))//'/z_grid'
		ker_dim_file = TRIM(ADJUSTL(ker_dir))//'/dimensions'
!		Write(6,*)'kernel_set is', kernel_set
!		Write(6,*)'ker_set_string is ', ker_set_string
!		Write(6,*)'ker_dir is ', ker_dir
!		Write(6,*)'ker_z_file is ', ker_z_file
	
		Call Read_Ascii(z_kernel,ker_z_file)

		Open (unit = 10, file = ker_dim_file)
		READ(10,'(I7)')nxk
		READ(10,'(I7)')nyk
		READ(10,'(I7)')nzk	
		READ(10,'(I7)')nkers	
	        Close (unit = 10)	
	!	write(6,*)nxk, nyk, nzk	
		! find nxk, nyk, nzk
		! read kernel zgrid


	End Subroutine Initialize_Kernels

	Subroutine Load_Kernels(tsize)
		Implicit None
		Integer*8 :: tsize
		Character*120 :: nk_file, tsize_string, ker_file
		Write(tsize_string,'(i10)')tsize
		Write(6,*)'Loading kernel for tilesize : ', tsize_string
		   nk_file = TRIM(ADJUSTL(ker_dir))//'/nk_'//TRIM(ADJUSTL(tsize_string))//'.fits'
		  ker_file = TRIM(ADJUSTL(ker_dir))//'/kers_'//TRIM(ADJUSTL(tsize_string))//'.fits'
		If (Allocated(smeared_kers)) DeAllocate(smeared_kers)
		If (Allocated(nkindex)) DeAllocate(nkindex)
		Write(6,*)'Number of kernels = ', nkers
		Allocate(nkindex(2,nkers))
		Allocate(smeared_kers(1:nxk,1:nyk,1:nzk,1:nkers)) ! The kernel array for all tile sizes
		smeared_kers(:,:,:,:) = 0.0D0
                call readfits_4dd(smeared_kers,ker_file,nxk,nyk,nzk,nkers)
	!	smeared_kers = smeared_kers   ! go ahead and multiply by dxdy		
		nkindex(:,:) = 0
		call readfits_2dl(nkindex,nk_file,2,nkers)
		

	End Subroutine Load_Kernels

End Module Kernel_Management
