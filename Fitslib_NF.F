Module Fitslib_NF

Contains
	SUBROUTINE Prterr(ier,marker)

      	INTEGER ier,marker
      	CHARACTER errtext*30,errmessage*80


	IF (status .EQ. 0) RETURN

      	END SUBROUTINE Prterr

      SUBROUTINE readfits_1dl(data,filename,nx)

      implicit none
      integer*8 nx,ny,nz
      integer ier, naxes(1), j, i, hdu, blksize,nfnd
      integer nbuf, ifp, jfp,it
      character*80 comment, record
      character*120 filename
      integer*8 data(0:nx-1),g
      integer*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf


      ier = 0
      hdu = 1
      nv  = 0.0
      naxes(1) = nx
      call ftopen(hdu,filename,0,blksize,ier)
      call prterr(ier,10)
      call ftgknj(hdu,'NAXIS',1,1,naxes,nfnd,ier)
      call prterr(ier,20)
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
         call ftgpvk(hdu,g,ifp,nbuf,nv,buf,anyf,ier)
         call prterr(ier,30)
         ifp=ifp+nbuf
         do i=1,nbuf
           data(i-1) = buf(i)
         enddo
      call prterr(ier,100)
      call ftclos(hdu,ier)
      call prterr(ier,113)

      end subroutine readfits_1dl

      SUBROUTINE readfits_1dd(data,filename,nx)

      implicit none
      integer*8 nx,ny,nz
      integer ier, naxes(1), j, i, hdu, blksize,nfnd
      integer nbuf, ifp, jfp,it
      character*80 comment, record
      character*120 filename
      Real*8 data(0:nx-1),g
      Real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf


      ier = 0
      hdu = 1
      nv  = 0.0
      naxes(1) = nx
      call ftopen(hdu,filename,0,blksize,ier)
      call prterr(ier,10)
      call ftgknj(hdu,'NAXIS',1,1,naxes,nfnd,ier)
      call prterr(ier,20)
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
         call ftgpvd(hdu,g,ifp,nbuf,nv,buf,anyf,ier)
         call prterr(ier,30)
         ifp=ifp+nbuf
         do i=1,nbuf
           data(i-1) = buf(i)
         enddo
      call prterr(ier,100)
      call ftclos(hdu,ier)
      call prterr(ier,113)

      end subroutine readfits_1dd


      SUBROUTINE readfits_3dd(data,filename,nx,ny,nz)

      implicit none
      integer*8 nx,ny,nz
      integer ier, naxes(3), j, i, hdu, blksize,nfnd
      integer nbuf, ifp, jfp,it
      character*80 comment, record
      character*120 filename
      real*8 data(0:nx-1,0:ny-1,0:nz-1),g
      real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf


      ier = 0
      hdu = 1
      nv  = 0.0
      naxes(1) = nx
      naxes(2) = ny
      naxes(3) = nz
      call ftopen(hdu,filename,0,blksize,ier)
      call prterr(ier,10)
      call ftgknj(hdu,'NAXIS',1,3,naxes,nfnd,ier)
      call prterr(ier,20)
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
      do it=1, nz
      do j=1,ny
         call ftgpvd(hdu,g,ifp,nbuf,nv,buf,anyf,ier)
         call prterr(ier,30)
         ifp=ifp+nbuf
         do i=1,nbuf
           data(i-1,j-1,it-1) = buf(i)
         enddo
       enddo
      enddo
      call prterr(ier,100)
      call ftclos(hdu,ier)
      call prterr(ier,113)

      end subroutine readfits_3dd

      SUBROUTINE readfits_3dl(data,filename,nx,ny,nz)

      implicit none
      integer*8 nx,ny,nz
      integer ier, naxes(3), j, i, hdu, blksize,nfnd
      integer nbuf, ifp, jfp,it
      character*80 comment, record
      character*120 filename
      integer*8 data(0:nx-1,0:ny-1,0:nz-1),g
      integer*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf


      ier = 0
      hdu = 1
      nv  = 0.0
      naxes(1) = nx
      naxes(2) = ny
      naxes(3) = nz
      call ftopen(hdu,filename,0,blksize,ier)
      call prterr(ier,10)
      call ftgknj(hdu,'NAXIS',1,3,naxes,nfnd,ier)
      call prterr(ier,20)
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
      do it=1, nz
      do j=1,ny
         call ftgpvk(hdu,g,ifp,nbuf,nv,buf,anyf,ier)
         call prterr(ier,30)
         ifp=ifp+nbuf
         do i=1,nbuf
           data(i-1,j-1,it-1) = buf(i)
         enddo
       enddo
      enddo
      call prterr(ier,100)
      call ftclos(hdu,ier)
      call prterr(ier,113)

      end subroutine readfits_3dl



      SUBROUTINE readfits_2dd(data,filename,nx,ny)

      implicit none
      integer*8 nx,ny,nz
      integer ier, naxes(2), j, i, hdu, blksize,nfnd
      integer nbuf, ifp, jfp,it
      character*80 comment, record
      character*120 filename
      real*8 data(0:nx-1,0:ny-1),g
      real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf


      ier = 0
      hdu = 1
      nv  = 0.0
      naxes(1) = nx
      naxes(2) = ny
      call ftopen(hdu,filename,0,blksize,ier)
      call prterr(ier,10)
      call ftgknj(hdu,'NAXIS',1,2,naxes,nfnd,ier)
      call prterr(ier,20)
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
      do j=1,ny
         call ftgpvd(hdu,g,ifp,nbuf,nv,buf,anyf,ier)
         call prterr(ier,30)
         ifp=ifp+nbuf
         do i=1,nbuf
           data(i-1,j-1) = buf(i)
         enddo
       enddo
      call prterr(ier,100)
      call ftclos(hdu,ier)
      call prterr(ier,113)

      end subroutine readfits_2dd

      SUBROUTINE readfits_2dl(data,filename,nx,ny)

      implicit none
      integer*8 nx,ny,nz
      integer ier, naxes(2), j, i, hdu, blksize,nfnd
      integer nbuf, ifp, jfp,it
      character*80 comment, record
      character*120 filename
      integer*8 data(0:nx-1,0:ny-1),g
      integer*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf


      ier = 0
      hdu = 1
      nv  = 0.0
      naxes(1) = nx
      naxes(2) = ny
      call ftopen(hdu,filename,0,blksize,ier)
      call prterr(ier,10)
      call ftgknj(hdu,'NAXIS',1,2,naxes,nfnd,ier)
      call prterr(ier,20)
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
      do j=1,ny
         call ftgpvk(hdu,g,ifp,nbuf,nv,buf,anyf,ier)
         call prterr(ier,30)
         ifp=ifp+nbuf
         do i=1,nbuf
           data(i-1,j-1) = buf(i)
         enddo
       enddo
      call prterr(ier,100)
      call ftclos(hdu,ier)
      call prterr(ier,113)

      end subroutine readfits_2dl

      SUBROUTINE readfits_4dd(data,filename,nx,ny,nz,nt)

      implicit none
      integer*8 nx,ny,nz,nt
      integer ier, naxes(4), j, i, hdu, blksize,nfnd
      integer nbuf, ifp, jfp,it, iitt
      character*80 comment, record
      character*120 filename
      real*8 data(0:nx-1,0:ny-1,0:nz-1,0:nt-1),g
      real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf


      ier = 0
      hdu = 1
      nv  = 0.0
      naxes(1) = nx
      naxes(2) = ny
      naxes(3) = nz
      naxes(4) = nt
      call ftopen(hdu,filename,0,blksize,ier)
      call prterr(ier,10)
      call ftgknj(hdu,'NAXIS',1,4,naxes,nfnd,ier)
      call prterr(ier,20)
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
      do iitt = 1, nt
      do it=1, nz
      do j=1,ny
         call ftgpvd(hdu,g,ifp,nbuf,nv,buf,anyf,ier)
         call prterr(ier,30)
         ifp=ifp+nbuf
         do i=1,nbuf
           data(i-1,j-1,it-1,iitt-1) = buf(i)
         enddo
       enddo
      enddo
      enddo
      call prterr(ier,100)
      call ftclos(hdu,ier)
      call prterr(ier,113)

      end subroutine readfits_4dd

      SUBROUTINE readfits_4dl(data,filename,nx,ny,nz,nt)

      implicit none
      integer*8 nx,ny,nz,nt
      integer ier, naxes(4), j, i, hdu, blksize,nfnd
      integer nbuf, ifp, jfp,it, iitt
      character*80 comment, record
      character*120 filename
      integer*8 data(0:nx-1,0:ny-1,0:nz-1,0:nt-1),g
      integer*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf


      ier = 0
      hdu = 1
      nv  = 0.0
      naxes(1) = nx
      naxes(2) = ny
      naxes(3) = nz
      naxes(4) = nt
      call ftopen(hdu,filename,0,blksize,ier)
      call prterr(ier,10)
      call ftgknj(hdu,'NAXIS',1,4,naxes,nfnd,ier)
      call prterr(ier,20)
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
      do iitt = 1, nt
      do it=1, nz
      do j=1,ny
         call ftgpvk(hdu,g,ifp,nbuf,nv,buf,anyf,ier)
         call prterr(ier,30)
         ifp=ifp+nbuf
         do i=1,nbuf
           data(i-1,j-1,it-1,iitt-1) = buf(i)
         enddo
       enddo
      enddo
      enddo
      call prterr(ier,100)
      call ftclos(hdu,ier)
      call prterr(ier,113)

      end subroutine readfits_4dl


      SUBROUTINE writefits_1dd(data,filename,nx)

      implicit none
      integer*8 nx
      integer ier
      integer naxes(1), j, i, hdw, blksize,nfnd
      integer :: nbuf, ifp, jfp,it, btpix, naxis
      character*80 comment, record
      character*120 filename
      real*8 data(0:nx-1),g
      real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf, simp, extend
	! write(6,*)'we are in'
      ier = 0
      hdw = 1
      nv  = 0.0
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
!      write(6,*)'here 1'
      call prterr(ier,100)
 !      write(6,*)data(0:10,0,0)
       simp=.TRUE.
       extend=.TRUE.
       btpix=-32
       btpix=-64
       naxes(1)=nx
       naxis=1
       call ftinit(hdw,filename,blksize,ier)
       call prterr(ier,15)
!	Write(6,*)'here 2'
       call ftphpr(hdw,simp,btpix,naxis,naxes,0,1,extend,ier)
       call prterr(ier,33)
       call ftrdef(hdw,ier)
       call prterr(ier,34)
       call prterr(ier,132)
       call flush(6)
        do i=0,nx-1
	! write(6,*)i
         buf(i+1) = data(i)
        enddo
        call flush(6)
        call ftpprd(hdw,g,jfp,nx,buf,ier)
        call prterr(ier,70)
        jfp=jfp+nx
      call ftclos(hdw,ier)
      call prterr(ier,113)

      end subroutine writefits_1dd

      SUBROUTINE writefits_2dd(data,filename,nx,ny)

      implicit none
      integer*8 :: nx,ny
      integer :: ier, naxes(2), j, i, hdw, blksize,nfnd
      integer :: nbuf, ifp, jfp,it, btpix, naxis
      character*80 comment, record
      character*120 filename
      real*8 data(0:nx-1,0:ny-1),g
      real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf, simp, extend
!	write(6,*)'we are in'
      ier = 0
      hdw = 1
      nv  = 0.0
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
!      write(6,*)'here 1'
      call prterr(ier,100)
!      write(6,*)data(0:10,0,0)
       simp=.TRUE.
       extend=.TRUE.
       btpix=-32
       btpix=-64
       naxes(1)=nx
       naxes(2) = ny
       naxis=2
       call ftinit(hdw,filename,blksize,ier)
       call prterr(ier,15)
!	Write(6,*)'here 2'
       call ftphpr(hdw,simp,btpix,naxis,naxes,0,1,extend,ier)
       call prterr(ier,33)
       call ftrdef(hdw,ier)
       call prterr(ier,34)
       call prterr(ier,132)
       call flush(6)
	do j = 0, ny -1
        do i=0,nx-1
	! write(6,*)i
         buf(i+1) = data(i,j)
        enddo
        call flush(6)
        call ftpprd(hdw,g,jfp,nx,buf,ier)
        call prterr(ier,70)
        jfp=jfp+nx
	enddo
      call ftclos(hdw,ier)
      call prterr(ier,113)

      end subroutine writefits_2dd
      SUBROUTINE writefits_3dd(data,filename,nx,ny,nz)

      implicit none
      integer*8 :: nx,ny, nz
      integer :: ier, naxes(3), j, i,k, hdw, blksize,nfnd
      integer :: nbuf, ifp, jfp,it, btpix, naxis
      character*80 comment, record
      character*120 filename
      real*8 data(0:nx-1,0:ny-1,0:nz-1),g
      real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf, simp, extend
!	write(6,*)'we are in'
      ier = 0
      hdw = 1
      nv  = 0.0
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
!      write(6,*)'here 1'
      call prterr(ier,100)
!      write(6,*)data(0:10,0,0)
       simp=.TRUE.
       extend=.TRUE.
       btpix=-32
       btpix=-64
       naxes(1)=nx
       naxes(2) = ny
       naxes(3) = nz
       naxis=3
       call ftinit(hdw,filename,blksize,ier)
       call prterr(ier,15)
!	Write(6,*)'here 2'
       call ftphpr(hdw,simp,btpix,naxis,naxes,0,1,extend,ier)
       call prterr(ier,33)
       call ftrdef(hdw,ier)
       call prterr(ier,34)
       call prterr(ier,132)
       call flush(6)
	do k = 0, nz -1 
	do j = 0, ny -1
        do i=0,nx-1
	! write(6,*)i
         buf(i+1) = data(i,j,k)
        enddo
        call flush(6)
        call ftpprd(hdw,g,jfp,nx,buf,ier)
        call prterr(ier,70)
        jfp=jfp+nx
	enddo
	enddo
      call ftclos(hdw,ier)
      call prterr(ier,113)

      end subroutine writefits_3dd

      SUBROUTINE writefits_4dd(data,filename,nx,ny,nz,nt)

      implicit none
      integer*8 :: nx,ny, nz,nt
      integer :: ier, naxes(4), j, i,k,kk, hdw, blksize,nfnd
      integer :: nbuf, ifp, jfp,it, btpix, naxis
      character*80 comment, record
      character*120 filename
      real*8 data(0:nx-1,0:ny-1,0:nz-1,0:nt-1),g
      real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf, simp, extend
!	write(6,*)'we are in'
      ier = 0
      hdw = 1
      nv  = 0.0
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
!      write(6,*)'here 1'
      call prterr(ier,100)
!      write(6,*)data(0:10,0,0)
       simp=.TRUE.
       extend=.TRUE.
       btpix=-32
       btpix=-64
       naxes(1)=nx
       naxes(2) = ny
       naxes(3) = nz
       naxes(4) = nt
       naxis=4
       call ftinit(hdw,filename,blksize,ier)
       call prterr(ier,15)
!	Write(6,*)'here 2'
       call ftphpr(hdw,simp,btpix,naxis,naxes,0,1,extend,ier)
       call prterr(ier,33)
       call ftrdef(hdw,ier)
       call prterr(ier,34)
       call prterr(ier,132)
       call flush(6)
	do kk = 0, nt -1
	do k = 0, nz -1 
	do j = 0, ny -1
        do i=0,nx-1
	! write(6,*)i
         buf(i+1) = data(i,j,k,kk)
        enddo
        call flush(6)
        call ftpprd(hdw,g,jfp,nx,buf,ier)
        call prterr(ier,70)
        jfp=jfp+nx
	enddo
	enddo
	enddo
      call ftclos(hdw,ier)
      call prterr(ier,113)

      end subroutine writefits_4dd

      SUBROUTINE writefits_5dd(data,filename,nx,ny,nz,nt,n5)

      implicit none
      integer*8 :: nx,ny, nz,nt,n5
      integer :: ier, naxes(5), j,jj, i,k,kk, hdw, blksize,nfnd
      integer :: nbuf, ifp, jfp,it, btpix, naxis
      character*80 comment, record
      character*120 filename
      real*8 data(0:nx-1,0:ny-1,0:nz-1,0:nt-1,0:n5-1),g
      real*8 nv, buf(nx)
      real timedif, dtime, timearray(2)
      logical anyf, simp, extend
!	write(6,*)'we are in'
      ier = 0
      hdw = 1
      nv  = 0.0
      nbuf = nx
      ifp=1
      jfp=1
      g = 0.0
!      write(6,*)'here 1'
      call prterr(ier,100)
!      write(6,*)data(0:10,0,0)
       simp=.TRUE.
       extend=.TRUE.
       btpix=-32
       btpix=-64
       naxes(1)=nx
       naxes(2) = ny
       naxes(3) = nz
       naxes(4) = nt
       naxes(5) = n5
       naxis=5
       call ftinit(hdw,filename,blksize,ier)
       call prterr(ier,15)
!	Write(6,*)'here 2'
       call ftphpr(hdw,simp,btpix,naxis,naxes,0,1,extend,ier)
       call prterr(ier,33)
       call ftrdef(hdw,ier)
       call prterr(ier,34)
       call prterr(ier,132)
       call flush(6)
	do jj = 0, n5 -1
	do kk = 0, nt -1
	do k = 0, nz -1 
	do j = 0, ny -1
        do i=0,nx-1
	! write(6,*)i
         buf(i+1) = data(i,j,k,kk,jj)
        enddo
        call flush(6)
        call ftpprd(hdw,g,jfp,nx,buf,ier)
        call prterr(ier,70)
        jfp=jfp+nx
	enddo
	enddo
	enddo
	enddo
      call ftclos(hdw,ier)
      call prterr(ier,113)

      end subroutine writefits_5dd

      SUBROUTINE writefits_4dl(data,filename,nx,ny,nz,nt)
	! actually writes doubles - need to fix this later, but I'm just
	! trying to write some output to debug some arrays right now - NF 9/20/2010
      implicit none
      integer*8 :: nx,ny, nz,nt
      integer :: ier, naxes(4), j, i,k,kk, hdw, blksize,nfnd
      integer :: nbuf, ifp, jfp,it, btpix, naxis
      character*80 comment, record
      character*120 filename
!      real*8 g
      Integer*8 data(0:nx-1,0:ny-1,0:nz-1,0:nt-1), buf(nx), g
      real*8 nv
      real timedif, dtime, timearray(2)
      logical anyf, simp, extend
!	write(6,*)'we are in'
      ier = 0
      hdw = 1
      nv  = 0.0
      nbuf = nx
      ifp=1
      jfp=1
      g = int(0.0)
!      write(6,*)'here 1'
      call prterr(ier,100)
!      write(6,*)data(0:10,0,0)
       simp=.TRUE.
       extend=.TRUE.
       btpix=-32
       btpix=-64
   !    btpix=-8
       naxes(1)=nx
       naxes(2) = ny
       naxes(3) = nz
       naxes(4) = nt
       naxis=4
	write(6,*)naxes
!	write(6,*)data
       call ftinit(hdw,filename,blksize,ier)
       call prterr(ier,15)
!	Write(6,*)'here 2'
       call ftphpr(hdw,simp,btpix,naxis,naxes,0,1,extend,ier)
       call prterr(ier,33)
       call ftrdef(hdw,ier)
       call prterr(ier,34)
       call prterr(ier,132)
       call flush(6)
	do kk = 0, nt -1
	do k = 0, nz -1 
	do j = 0, ny -1
        do i=0,nx-1
	! write(6,*)i
         buf(i+1) = data(i,j,k,kk)
        enddo
        call flush(6)
        call ftpprk(hdw,g,jfp,nx,buf,ier)
        call prterr(ier,70)
        jfp=jfp+nx
	enddo
	enddo
	enddo
      call ftclos(hdw,ier)
      call prterr(ier,113)

      end subroutine writefits_4dl
      SUBROUTINE writefits_3dl(data,filename,nx,ny,nz)
	! actually writes doubles - need to fix this later, but I'm just
	! trying to write some output to debug some arrays right now - NF 9/20/2010
      implicit none
      integer*8 :: nx,ny, nz
      integer :: ier, naxes(3), j, i,k,kk, hdw, blksize,nfnd
      integer :: nbuf, ifp, jfp,it, btpix, naxis
      character*80 comment, record
      character*120 filename
!      real*8 g
      Integer*8 data(0:nx-1,0:ny-1,0:nz-1), buf(nx), g
      real*8 nv
      real timedif, dtime, timearray(2)
      logical anyf, simp, extend
!	write(6,*)'we are in'
      ier = 0
      hdw = 1
      nv  = 0.0
      nbuf = nx
      ifp=1
      jfp=1
      g = int(0.0)
!      write(6,*)'here 1'
      call prterr(ier,100)
!      write(6,*)data(0:10,0,0)
       simp=.TRUE.
       extend=.TRUE.
       btpix=-32
       btpix=-64
   !    btpix=-8
       naxes(1)=nx
       naxes(2) = ny
       naxes(3) = nz
       naxis=3
	write(6,*)naxes
!	write(6,*)data
       call ftinit(hdw,filename,blksize,ier)
       call prterr(ier,15)
!	Write(6,*)'here 2'
       call ftphpr(hdw,simp,btpix,naxis,naxes,0,1,extend,ier)
       call prterr(ier,33)
       call ftrdef(hdw,ier)
       call prterr(ier,34)
       call prterr(ier,132)
       call flush(6)
!	do kk = 0, nt -1
	do k = 0, nz -1 
	do j = 0, ny -1
        do i=0,nx-1
	! write(6,*)i
         buf(i+1) = data(i,j,k)
        enddo
        call flush(6)
        call ftpprk(hdw,g,jfp,nx,buf,ier)
        call prterr(ier,70)
        jfp=jfp+nx
	enddo
	enddo
!	enddo
      call ftclos(hdw,ier)
      call prterr(ier,113)

      end subroutine writefits_3dl
      SUBROUTINE writefits_2dl(data,filename,nx,ny)
	! actually writes doubles - need to fix this later, but I'm just
	! trying to write some output to debug some arrays right now - NF 9/20/2010
      implicit none
      integer*8 :: nx,ny
      integer :: ier, naxes(2), j, i,k,kk, hdw, blksize,nfnd
      integer :: nbuf, ifp, jfp,it, btpix, naxis
      character*80 comment, record
      character*120 filename
!      real*8 g
      Integer*8 data(0:nx-1,0:ny-1), buf(nx), g
      real*8 nv
      real timedif, dtime, timearray(2)
      logical anyf, simp, extend
!	write(6,*)'we are in'
      ier = 0
      hdw = 1
      nv  = 0.0
      nbuf = nx
      ifp=1
      jfp=1
      g = int(0.0)
!      write(6,*)'here 1'
      call prterr(ier,100)
!      write(6,*)data(0:10,0,0)
       simp=.TRUE.
       extend=.TRUE.
       btpix=-32
       btpix=-64
   !    btpix=-8
       naxes(1)=nx
       naxes(2) = ny
       naxis=2
	write(6,*)naxes
!	write(6,*)data
       call ftinit(hdw,filename,blksize,ier)
       call prterr(ier,15)
!	Write(6,*)'here 2'
       call ftphpr(hdw,simp,btpix,naxis,naxes,0,1,extend,ier)
       call prterr(ier,33)
       call ftrdef(hdw,ier)
       call prterr(ier,34)
       call prterr(ier,132)
       call flush(6)
!	do kk = 0, nt -1

	do j = 0, ny -1
        do i=0,nx-1
	! write(6,*)i
         buf(i+1) = data(i,j)
        enddo
        call flush(6)
        call ftpprk(hdw,g,jfp,nx,buf,ier)
        call prterr(ier,70)
        jfp=jfp+nx
	enddo
	
!	enddo
      call ftclos(hdw,ier)
      call prterr(ier,113)

      end subroutine writefits_2dl
End Module Fitslib_NF
