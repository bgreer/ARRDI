MODULE Kernel_Management
	USE IO_Routines
	USE FitsLib_NF
	USE Run_Info
	IMPLICIT NONE
	INTEGER*8, PUBLIC :: nxk, nyk, nzk, nkers, kernel_set
	INTEGER*8, ALLOCATABLE, PUBLIC :: z_owner(:,:), nkindex(:,:)
	REAL*8, ALLOCATABLE, PUBLIC :: wz(:), smeared_kers(:,:,:,:), z_kernel(:)
	CHARACTER*120 :: ker_dir

	CONTAINS

	! Read kernel set dimensions
	! Path is hard-coded
	SUBROUTINE Initialize_Kernels()
		IMPLICIT NONE
		CHARACTER*120 :: ker_dim_file, ker_z_file, ker_set_string
		WRITE(ker_set_string,'(i10)')kernel_set
		ker_dir = '/polaris/kernel_sets/'//TRIM(ADJUSTL(ker_set_string))
		ker_z_file = TRIM(ADJUSTL(ker_dir))//'/z_grid'
		ker_dim_file = TRIM(ADJUSTL(ker_dir))//'/dimensions'
	
		CALL Read_Ascii(z_kernel,ker_z_file)

		OPEN(10, FILE=ker_dim_file)
		READ(10,'(I7)')nxk
		READ(10,'(I7)')nyk
		READ(10,'(I7)')nzk
		READ(10,'(I7)')nkers
		CLOSE(10)
	END SUBROUTINE Initialize_Kernels

	! read a kernel set (single tile size) into a 4d array
	SUBROUTINE Load_Kernels(tsize)
		IMPLICIT NONE
		INTEGER*8 :: tsize
		CHARACTER*120 :: nk_file, tsize_string, ker_file

		WRITE(tsize_string,'(i10)')tsize
		IF (myid == 0) WRITE(*,'(A,I4,A,A)')' Loading ',nkers,' kernels for tilesize : ', TRIM(tsize_string)

		nk_file = TRIM(ADJUSTL(ker_dir))//'/nk_'//TRIM(ADJUSTL(tsize_string))//'.fits'
		ker_file = TRIM(ADJUSTL(ker_dir))//'/kers_'//TRIM(ADJUSTL(tsize_string))//'.fits'
		IF (ALLOCATED(smeared_kers)) DEALLOCATE(smeared_kers)
		IF (ALLOCATED(nkindex)) DEALLOCATE(nkindex)
		ALLOCATE(nkindex(2,nkers))
		ALLOCATE(smeared_kers(1:nxk,1:nyk,1:nzk,1:nkers)) ! The kernel array for all tile sizes

		smeared_kers(:,:,:,:) = 0.0D0
		CALL readfits_4dd(smeared_kers,ker_file,nxk,nyk,nzk,nkers)
		nkindex(:,:) = 0
		CALL readfits_2dl(nkindex,nk_file,2,nkers)
	END SUBROUTINE Load_Kernels

END MODULE Kernel_Management
