
c===============================================================c
c                                                               c
c                       Prterr	                        	c
c                                                               c
c	This routine prints error messages when a FITSIO error	c
c   occurs.  The appropriate error message is printed as well	c
c   as a marker number to determine where in the calling code	c
c   the error occured.						c
c                                                               c
c===============================================================c


	SUBROUTINE Prterr(status,marker)

      	INTEGER status,marker
      	CHARACTER errtext*30,errmessage*80


	IF (status .EQ. 0) RETURN

      	END


