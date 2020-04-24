c-----------------------------------------------------------------------
	program myspleval
c-----------------------------------------------------------------------
c
c  Leah Ziegler - Nov. 2010
c    questions: email to lziegler@ucsd.edu
c
c  Fortran 77 code made on MacBook Pro and compiled with gfortran. 
c     Untested on other platforms.
c
c
c    Input a parameter file of cubic B-spline model and output the model
c     evaluated either at specified regular intervals, or at a specific
c     list of times from a file.
c
c    IN:  parameter file (e.g. padm2m.2b.pf) formatted as:
c		nknot, tstart, tend
c		beta(1)
c		beta(2)
c	        .
c	        .
c		.
c		beta(n)
c
c	 where nknot is the number of spline knot points, tstart and
c        tend are the start and end time of the model, and beta() are
c        the spline coefficient values
c
c        If an input files supplies the times to evaluate than note:
c           NOTE:  format of times list is one age per line 
c	    NOTE:  times don't have to be evenly spaced, and don't have to 
c		   span the full 2 Myr of PADM2M, HOWEVER they must
c		   be sorted into sequential order...youngest to oldest
c	    NOTE:  the first line of the output file is always the model
c                  at time=tstart (0 Ma for padm2m)
c		   and the last is always the model at time=tend (2 Ma for padm2m)
c		   as specified in the model parameter file
c
c   OUT:  outfile is fort.21 and is formatted as:
c	time (in Myr ago) ADM (in x10^22 Am^2) derivative 2nd_deriv 


c-----------------------------------------------------------------------
c  Defining some variables
c-----------------------------------------------------------------------
c
c********* If parameter values are changed here they must also be changed in
c********* subroutines splev, dsplev, asplev, rowsd, rowsdd and rowsb
c
      parameter (ndmax=1000000,nbmax=3000,nsub=1000)
      implicit double precision (a-h,o-z)
c
c     ndmax   maximum number of data allowed in time series
c     nbmax    2 more than the maximum number of knots allowed
c    
      common /splin1/npts,nknot,data(ndmax),time(ndmax),t(nbmax),splin
     $(ndmax),dsplin(ndmax),asplin(ndmax) 
c
c     npts    number of data pairs
c     nknot   number of (equally spaced) knots along length of data series. 
c             This will be incremented by 2 to complete the B-spline basis.
c     data(i) value of data at time(i) - LEGACY, set to 0 here or DELETE
c     time(i) ith value in time array (note that these are not required to 
c             be equally spaced)
c     t(j)    position of jth knot along time array
c	splin(i) value of penalized spline at time(i)
c	dsplin(i) derivative of splin(i) at time(i)
c	asplin(i) 2nd derivative of splin(i) at time(i)
c
c
c
c******* unit 15 must be opened in calling program 
c        10 is used to print arrays for debugging, when iflag=1
c
c   set iflag=1 in calling prog to print arrays for debugging
c      common /data/deltat,iflag

      data nunit1/5/,nunit2/6/,nb/4/
      external rowsb
      character*100 timein,paramin,linein
      double precision tstart,tend,tmptime,betatmp,beta(nbmax),tsamp
      double precision splvadm,dsplvadm,asplvadm
      double precision b(4*nbmax*nbmax)
      integer nptsout



c-----------------------------------------------------------------------
c  Read in user-defined sampling interval (tsamp) or time array (timein) and 
c  spline parameters (from paramin)
c-----------------------------------------------------------------------

	print *,' Enter name of file containing model params: '
	read(5,'(a)')paramin
	print *, ' Enter time sampling interval for eval. of spline'
	print *, ' or 0 to read list of times from file: '
	read(5,*)tsamp
	if(tsamp.eq.0)then
	   print *,'Enter file name of times: '
	   read(5,'(a)')timein
	   print *,'Enter length of file: '
	   read(5,*)nintimes
	end if

	open(12,file=paramin,status='old')
	read(12,'(a)')linein
	call MYSPLIT(linein,nknot,tstart,tend)
	do j=1,nknot
	   read(12,*)betatmp
	   beta(j)=betatmp
	end do
	close(12)

c-----------------------------------------------------------------------
c  set data vec to 0 (because i'm adapting to some legacy code instead of
c  changing the legacy subroutines which would be better of me)
c-----------------------------------------------------------------------
	 	do j=1,ndmax
		   data(j)=0.
		end do


c-----------------------------------------------------------------------
c  set up time vector putting time min and max on ends of user-given
c  time vector 
c-----------------------------------------------------------------------

             if(tsamp.eq.0)then
                open(13,file=timein,status='old')
		time(1)=tstart
                do i=2,nintimes+1
                   read(13,*)tmptime
                   time(i)=tmptime
		end do
                close(13) 
                nptsout=nintimes+2
                time(nptsout)=tend
             else
                 npts=int((tend-tstart)/tsamp) +1
                 do 195 i=1,npts-1
195                 time(i)=tstart +tsamp*float(i-1)
                 time(npts)=tend
                 nptsout=npts
             end if
             call splev(beta,b,nptsout)
             call dsplev(beta,b,nptsout)
             call asplev(beta,b,nptsout)
             do j=1,nptsout
		splvadm=splin(j)*.258597                
		dsplvadm=dsplin(j)*.258597                
		asplvadm=asplin(j)*.258597                
                write(21,111)time(j),splvadm,dsplvadm,asplvadm
111		  format(4g18.8)
	     end do




      end

c__________________________________________________________________________
c__________________________________________________________________________
c
c  SUBROUTINES
c__________________________________________________________________________
c__________________________________________________________________________

c-----------------------------------------------------------------------
      subroutine MYSPLIT(linein,nknot,tstart,tend)
c-----------------------------------------------------------------------


      implicit none
      character (len=100) :: linein,numstring
      character *1 delim
      integer k,next,length,first,last,nknot
      double precision tstart,tend
      parameter (length=100,delim=',')

         k=0
         first=1
         do 10 next =1,length
            if (linein(next:next).eq.delim) then
               k=k+1
               last=next-1
               if (k.eq.1) then
                  numstring=linein(first:last)
                  read(numstring,*)nknot
		  print *,'nknot = ',nknot
               else if (k.eq.2) then
                  numstring=linein(first:last)
                  read(numstring,410)tstart
410		  format(g15.7)
		  print *,'tstart = ',tstart
               else if (k.eq.3) then
                  numstring=linein(first:last)
                  read(numstring,411)tend
411		  format(g15.7)
		  print *,'tend = ',tend
               endif
               first=next+1
            endif
10       continue
      return
      end subroutine

c__________________________________________________________________________
      subroutine rowsb(g,mg,mt,nb1,ndata,jt)
c__________________________________________________________________________
c  subroutine to compute g matrix for banded sequential least
c  squares qr solution of spline fit to data
c  jt is first non-zero column for the given block
c
c*************calls array, basis **************
c
      implicit double precision (a-h,o-z)
      parameter(nbmax=3000,ndmax=1000000)
      common /splin1/npts,nknot,data(ndmax),time(ndmax),t(nbmax),splin
     $(ndmax),dsplin(ndmax),asplin(ndmax) 
      common /data/deltat,iflag
      dimension g(mg,*)
      data nunit2/6/
c
c first set up t the knot array, knots are equally spaced along data.
c
       twhole=time(npts)-time(1)
       deltat = twhole/float(nknot-3)
      if(ndata.eq.0)then
       write(nunit2,100)deltat
100    format(' knots at intervals of ',g15.7)
       t(1)=time(1)-deltat
       t(2)=time(1)
       do 5 j=3,nknot
5      t(j)=t(1)+(j-1)*deltat
       t(nknot-1)=time(npts)
      endif
c
c     establish jt and knot position at which to start
c
      nd=ndata
      do 8 k=3,nknot
        if(time(ndata+1).le.t(k))then
         jt=k-2
         tt=t(k)
         go to 20
         endif
8     continue
c  construct g matrix
c
20    do 10 j=1,mt
      ndata=ndata+1
        do 25 l=1,nb1-1
25        g(j,l)=basis(t(l+jt-1),time(ndata),deltat)
        g(j,nb1)=data(ndata)
       if(time(ndata+1).gt.tt)goto 30
10    continue
30    mt=ndata-nd
c      if(iflag.eq.1)call array(g,mg,mt,nb1,10)
      return
      end
c
c__________________________________________________________________________
      double precision function basis(tl,timej,deltat)
c
c  Evaluates the cubic spline basis functions for knot tl at timej, for a 
c  knot spacing deltat
c
      implicit double precision (a-h,o-z)
      tt=(timej-tl)/deltat
      if(-2.0.le.tt.and.-1.0.gt.tt)then
        tt=tt+2.d0
        basis= tt*tt*tt
      else if (-1.0.le.tt.and.0.0.gt.tt)then
        tsq=tt*tt
        basis= 4.0 - 6.0*tt*tt - 3.0*tsq*tt
      else if (0.0.le.tt.and.1.0.gt.tt) then
        tsq=tt*tt
        basis= 4.0 - 6.0*tt*tt +3.0*tt*tsq
      else if(1.0.le.tt.and.2.0.ge.tt)then
        tt=2.d0-tt
        basis= tt*tt*tt
      else
        basis=0.d0
      endif
      return
      end

c__________________________________________________________________________
      subroutine splev(beta,b,nin)
c__________________________________________________________________________
c  spline evaluation routine
c
c  computes value of spline for array time.
c
c*********calls rowsb
c
      implicit double precision (a-h,o-z)
      parameter (ndmax=1000000,nbmax=3000,nsub=1000)
      dimension b(nsub,*), beta(nbmax)
      common /splin1/npts,nknot,data(ndmax),time(ndmax),t(nbmax),splin
     $(ndmax),dsplin(ndmax),asplin(ndmax) 
      data nb/4/
c  set up spline function bases for time array
      npts1=npts
      npts= nin
      mb=nsub
      nb1=nb+1
      ndata=0
      do 5000 k=1,100000
      mt=min(nsub,nin-ndata)
      nd=ndata+1
      call rowsb(b,mb,mt,nb1,ndata,jt)
      do 10 j=nd,ndata
      splin(j)=0.d0
      do 15 l=1,nb
15     splin(j)=splin(j) + b(j-nd+1,l)*beta(l+jt-1)
10    continue
      if(ndata.ge.nin)go to 6000
5000  continue
6000  npts=npts1
      return
      end
c
c******** subroutines for derivative calculations******************
c
c
c________________________________________________________________________
      subroutine rowsd(g,mg,mt,nb1,ndata,jt)
c  subroutine to compute g matrix for evaluation of first derivative of 
c  spline
c  jt is first non-zero column for the given block
c
c*********** calls dbasis
c
      implicit double precision (a-h,o-z)
      parameter(nbmax=3000,ndmax=1000000)
      common /splin1/npts,nknot,data(ndmax),time(ndmax),t(nbmax),splin
     $(ndmax),dsplin(ndmax),asplin(ndmax) 
      common /data/deltat,iflag
      dimension g(mg,*)
      data nunit2/6/
c     first set up t the knot array
       twhole=time(npts)-time(1)
       deltat = twhole/float(nknot-3)
      if(ndata.eq.0.and.mt.ne.3)then
       write(nunit2,100)deltat
100    format(' knots at intervals of ',g15.7)
       t(1)=time(1)-deltat
       t(2)=time(1)
       do 5 j=3,nknot
5      t(j)=t(1)+(j-1)*deltat
       t(nknot-1)=time(npts)
      endif
c
c     establish jt and knot position at which to start
      nd=ndata
      do 8 k=3,nknot
        if(time(ndata+1).le.t(k))then
         jt=k-2
         tt=t(k)
         go to 20
         endif
8     continue
c     construct g matrix
c
20    do 10 j=1,mt
      ndata=ndata+1
        do 25 l=1,nb1-1
25        g(j,l)=dbasis(t(l+jt-1),time(ndata),deltat)/deltat
        g(j,nb1)=data(ndata)
       if(time(ndata+1).gt.tt)goto 30
10    continue
30    mt=ndata-nd
c      if(iflag.eq.1)call array(g,mg,mt,nb1,10)
      return
      end
c
c___________________________________________________________________________
      double precision function dbasis(tl,timej,deltat)
c
c  finds derivative of cubic spline basis function for knot tl at timej,
c  with equal spacing deltat between knots
c
c**********calls no others
c
      implicit double precision (a-h,o-z)
      tt=(timej-tl)/deltat
      if(-2.0.le.tt.and.-1.0.gt.tt)then
        tt=tt+2.d0
        dbasis= 3*tt*tt
      else if (-1.0.le.tt.and.0.0.gt.tt)then
        tsq=tt*tt
        dbasis=  - 12.0*tt - 9.0*tsq
      else if (0.0.le.tt.and.1.0.gt.tt) then
        tsq=tt*tt
        dbasis=  - 12.0*tt +9.0*tsq
      else if(1.0.le.tt.and.2.0.ge.tt)then
        tt=2.d0-tt
        dbasis= -3.0*tt*tt
      else
        dbasis=0.d0
      endif
      return
      end
c
c___________________________________________________________________________
      subroutine dsplev(beta,b,nin)
c     computes value of spline derivative for array time.
c
c *******calls rowsd
c
      implicit double precision (a-h,o-z)
      parameter (ndmax=1000000,nbmax=3000,nsub=1000)
      dimension data(ndmax),b(nsub,*),
     $time(ndmax),beta(nbmax),t(nbmax),splin(ndmax),dsplin(ndmax),
     $asplin(ndmax)
      common /splin1/npts,nknot,data,time,t,splin,dsplin,asplin
      data nb/4/
c     set up spline function bases for time array
      npts1=npts
      npts= nin
      mb=nsub
      nb1=nb+1
      ndata=0
      do 5000 k=1,100000
      mt=min(nsub,nin-ndata)
      nd=ndata+1
      call rowsd(b,mb,mt,nb1,ndata,jt)
      do 10 j=nd,ndata
      dsplin(j)=0.d0
      do 15 l=1,nb
15     dsplin(j)=dsplin(j) + b(j-nd+1,l)*beta(l+jt-1)
10    continue
      if(ndata.ge.nin)go to 6000
5000  continue
6000  npts=npts1
      return
      end
c
c____________________________________________________________________________
      subroutine rowsdd(g,mg,mt,nb1,ndata,jt)
c  subroutine to compute g matrix for evaluation of second derivative of
c  spline fit to data
c    
c  jt is first non-zero column for the given block
c
c********** calls abasis
c
      implicit double precision (a-h,o-z)
      parameter(nbmax=3000,ndmax=1000000)
      common /splin1/npts,nknot,data(ndmax),time(ndmax),t(nbmax),splin
     $(ndmax),dsplin(ndmax),asplin(ndmax) 
      common /data/deltat,iflag
      dimension g(mg,*)
      data nunit2/6/
c     first set up t the knot array
       twhole=time(npts)-time(1)
       deltat = twhole/float(nknot-3)
      if(ndata.eq.0.and.mt.ne.3)then
       write(nunit2,100)deltat
100    format(' knots at intervals of ',g15.7)
       t(1)=time(1)-deltat
       t(2)=time(1)
       do 5 j=3,nknot
5      t(j)=t(1)+(j-1)*deltat
       t(nknot-1)=time(npts)
      endif
c
c     establish jt and knot position at which to start
      nd=ndata
      do 8 k=3,nknot
        if(time(ndata+1).le.t(k))then
         jt=k-2
         tt=t(k)
         go to 20
         endif
8     continue
c     construct g matrix
c
20    do 10 j=1,mt
      ndata=ndata+1
        do 25 l=1,nb1-1
25        g(j,l)=(abasis(t(l+jt-1),time(ndata),deltat))/(deltat**2)
        g(j,nb1)=data(ndata)
       if(time(ndata+1).gt.tt)goto 30
10    continue
30    mt=ndata-nd
c      if(iflag.eq.1)call array(g,mg,mt,nb1,10)
      return
      end
c
c___________________________________________________________________________
      double precision function abasis(tl,timej,deltat)
c  find 2nd derivative of spline basis functions for knot tl, at timej, 
c  equally spaced knots deltat
      implicit double precision (a-h,o-z)
      tt=(timej-tl)/deltat
      if(-2.0.le.tt.and.-1.0.gt.tt)then
        tt=tt+2.d0
        abasis= 6*tt
      else if (-1.0.le.tt.and.0.0.gt.tt)then
        abasis=  - 12.0 - 18.0*tt
      else if (0.0.le.tt.and.1.0.gt.tt) then
        abasis=  - 12.0 + 18.0*tt
      else if(1.0.le.tt.and.2.0.ge.tt)then
        tt=tt-2.d0
        abasis= -6.0*tt
      else
        abasis=0.d0
      endif
      return
      end
c
c___________________________________________________________________________
      subroutine asplev(beta,b,nin)
c  computes value of spline 2nd derivative for array time.
c
c**********calls rowsdd
c
      implicit double precision (a-h,o-z)
      parameter (ndmax=1000000,nbmax=3000,nsub=1000)
      dimension data(ndmax),b(nsub,*),
     $time(ndmax),beta(nbmax),t(nbmax),splin(ndmax),dsplin(ndmax),
     $asplin(ndmax)
      common /splin1/npts,nknot,data,time,t,splin,dsplin,asplin
      data nb/4/
c     set up spline function bases for time array
      npts1=npts
      npts= nin
      mb=nsub
      nb1=nb+1
      ndata=0
      do 5000 k=1,100000
      mt=min(nsub,nin-ndata)
      nd=ndata+1
      call rowsdd(b,mb,mt,nb1,ndata,jt)
      do 10 j=nd,ndata
      asplin(j)=0.d0
      do 15 l=1,nb
15      asplin(j)=asplin(j) + b(j-nd+1,l)*beta(l+jt-1)
10    continue
      if(ndata.ge.nin)go to 6000
5000  continue
6000  npts=npts1
      return
      end
