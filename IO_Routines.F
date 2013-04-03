Module IO_Routines
  Contains
  Subroutine Read_Ascii_defunct(arr,dfile)
	Implicit None
	Character*120 :: dfile
	Real*8, Intent(InOut) :: arr(:)
	Integer*8 :: nd, i
	open (unit = 10, file = dfile)
	READ(10,'(I7)')nd	
	DO i = 1, nd
		READ(10,'(E13.5)')arr(i)
	ENDDO
	close (unit = 10)	
  End Subroutine Read_Ascii_defunct

  Subroutine Read_Ascii(arr,dfile)
	Implicit None
	Character*120 :: dfile
	Real*8, Allocatable, Intent(InOut) :: arr(:)
	Integer*8 :: nd, i
	open (unit = 10, file = dfile)
	READ(10,'(I7)')nd
	IF (.not. Allocated(arr)) Then	
	     Allocate(arr(nd))
	Endif
	DO i = 1, nd
		READ(10,'(E19.12)')arr(i)
	ENDDO
	close (unit = 10)	
  End Subroutine Read_Ascii

  Subroutine Read_Ascii_Int(arr,dfile)
	Implicit None
	Character*120 :: dfile
	Real*8 :: temp
	Integer*8, Allocatable, Intent(InOut) :: arr(:)
	Integer*8 :: nd, i
	open (unit = 10, file = dfile)
	READ(10,'(I7)')nd	
	IF (.not. Allocated(arr)) Then	
	     Allocate(arr(nd))
	Endif
	DO i = 1, nd
		READ(10,'(E19.12)')temp
		arr(i) = int(temp)
	ENDDO
	close (unit = 10)	
  End Subroutine Read_Ascii_Int


End Module IO_Routines
