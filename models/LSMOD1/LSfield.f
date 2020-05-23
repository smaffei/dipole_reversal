      program LSfield
c
c     A program to calculate surface field prediction time series of
c     declination (D), inclination (I) and intensity (F) 
c     from time dependent coefficients in model file LSMOD.1 or LSMOD.2.
c     Output file format: 4 columns:
c     Age (ka BP)   D  I  F

c    Copyright © 2018, Monika Korte, Helmholtz Cenntre Potsdam GFZ German 
c    Research Centre f or Geosciecnes, Potsdam, Germany

c    Licensed under the Apache License, Version 2.0 (the "License");
c    you may not use this file except in compliance with the License.
c    You may obtain a copy of the License at
c    http://www.apache.org/licenses/LICENSE-2.0

c    Unless required by applicable law or agreed to in writing, software
c    distributed under the License is distributed on an "AS IS" BASIS,
c    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
c    See the License for the specific language governing permissions and
c    imitations under the License.

c    Citation: 
c    Korte, Monika; Brown, Maxwell; Gunnarson, Sydney (2018): LSMOD.1 - Global 
c    paleomagnetic field model for 50 -- 30 ka BP. V. 1. GFZ Data Services. 
c    http://doi.org/10.5880/GFZ.2.3.2018.008

c     December 2018, 
c     based on code for Holocene models,
c     version for CALS10k.1b June 2011
c     version for CALS10k.2  Feb. 2016
c
c     Uses pieces of code from the example program for evaluating 
c     time-dependent field model GUFM
c     by Jeremy Bloxham & Andrew Jackson,
c     which uses code by: David Gubbins, Kathy Whaler, David Barraclough,
c                         Rick O'Connell, and Carl de Boor
c
c---------------------------------------------------------------------      
c     CALLS:    interv   - calculates which knot lies immediately left
c                          of the current time point
c
c               bspline  - calculates B-splines at current time point
c---------------------------------------------------------------------

      implicit none

      integer lmax,nspl,n,np,nl,jord,nsplt
      real*8 gt,spl,tknts,g,gd,p,dp,dx,dy,dz,fac,dum
      character*15 infile
      character*25 outfile

      parameter (lmax=10)
      parameter (nsplt=403)

      parameter (n=lmax*(lmax+2))
      parameter (np=n*nsplt)
      parameter (nl=(lmax+1)*(lmax+2)/2)

      dimension gt(n,nsplt)
      dimension spl(nsplt),tknts(nsplt+4)
      
      dimension g(n),gd(n)
      
      dimension p(nl),dp(nl)
      dimension dx(lmax*(lmax+2)),dy(lmax*(lmax+2)),dz(lmax*(lmax+2))

      integer*4 it1,it2
      integer lm,nm,k,j,i,nleft,model
      real*8 alat,alon,alt,time,theta,phi,rad,sinth,costh,sd,cd
      real*8 x,y,z,h,f,ainc,d,tstartin,tendin,ki

      data jord/4/
      data fac/1.74532925e-2/

      write(*,*) 'Program to write time series of D, I, F'
      write(*,*) 'from model LSMOD.1 or LSMOD.2 to a file'
      write(*,*) '             ------------              '
      write(*,*) 'Give model number:'
      write(*,*)  '1 -- LSMOD.1'
      write(*,*) '2 -- LSMOD.2'
      read(*,*) model
      write(*,*) 'Give output file name:'
      read(*,*) outfile
      write(*,*) 'Give latitude (decimal degrees):'
      read(*,*) alat
      write(*,*) 'Give longitude (decimal degrees):'
      read(*,*) alon

      if (model.eq.1) then
         infile='LSMOD.1'
      else if (model.eq.2) then
         infile='LSMOD.2'
      else
         write(*,*) 'Error: Invalid model number!'
         stop
      endif

c********************************************************************
c     read model
      open(7,file=infile)
      read(7,*) tstartin,tendin
      read(7,*) lm,nm,nspl
      read(7,*) (tknts(i),i=1,nspl+4)
      read(7,*) gt
      close(7)


      it1=int(tstartin)
      it2=int(tendin)

      open(11,file=outfile) 
      write(11,*) 'Age(ka BP)   D(deg.)   I(deg.)  F(microT)'

      do i=it1,it2,20
      time = float(i)
      alt=0.0
      theta = (90.0-alat)*fac
      phi   = alon*fac
           
c-----
c     transform coordinates to geocentric colatitude,longitude,radius

      call coords(alt,theta,rad,sd,cd)      
      sinth=sin(theta)
      costh=cos(theta)
      call plmbar(p,dp,costh,lmax)

c-----
c     calculate main field coefficients at time time
c
      
10    call interv(tknts,time,nspl,nleft)
      call bspline(tknts,time,nspl,jord,nleft,spl(nleft-3))
      
      do  k=1,n
       g(k)=0.0
       do j=1,4
        g(k) = g(k) + spl(j+nleft-4)*gt(k,j+nleft-4)
       enddo 
      enddo 
          

c     calculate main field elements X,Y,Z,H,F,I(ainc) and D

      call magfdz(p,dp,theta,phi,rad,lmax,g,dx,dy,dz,x,y,z,h,f,
     >ainc,d,sd,cd)
     
c     convert inclination and declination to degrees     
      ainc=ainc/fac
      d=d/fac

      
      ki=(1950.0-int(i))/1000.
      write(11,6200) ki,d,ainc,f/1000.
      end do
99    continue
6200  format(f11.2,2f10.2,f10.1)
      close(11)
c********************************************************************
      stop
      end
      
c--------------------------------------------------------------------------      
      
      subroutine interv(tknts,time,nspl,nleft)
      implicit real*8 (a-h,o-z)
      
      dimension tknts(nspl+4)
     
      if(time.lt.tknts(4).or.time.gt.tknts(nspl+1)) return
      
      do 200 n=5,nspl+1
       if(time.le.tknts(n)) then
        nleft=n-1
        goto 210
       endif
200   continue
210   continue


      return
      end

c-------------------------------------------------------------------

       subroutine bspline(tknts,t,nspl,jorder,nleft,spl)
 
c calculate splines of order jorder where 1 <= jorder <= 4
       implicit real*8 (a-h,o-z)
       dimension tknts(nspl+4)
       dimension spl(4)
       
       dimension deltal(4),deltar(4)
       
       spl(1)=1.0
      
       do 200 j=1,jorder-1
       
       deltar(j) = tknts(nleft+j) - t
       deltal(j) = t - tknts(nleft+1-j)
       saved=0.0
       
       do 100 i=1,j
        term = spl(i)/(deltar(i)+deltal(j+1-i))
        spl(i) = saved + deltar(i)*term
        saved = deltal(j+1-i)*term
100    continue

       spl(j+1) = saved
       
200    continue

       
 
       return
       end
        
c-----------------------------------------------------------------

      subroutine coords(h,theta,r,sd,cd)
      implicit real*8 (a-h,o-z)

      pi=3.14159265
      b1=40680925.
      b2=40408585.
      theta=pi/2-theta
      clat=cos(theta)
      slat=sin(theta)
      one=b1*clat*clat
      two=b2*slat*slat
      three=one+two
      four=sqrt(three)
      r=sqrt(h*(h+2.*four)+(b1*one+b2*two)/three)
      cd=(h+four)/r
      sd=(b1-b2)/four*slat*clat/r
      sinth=slat*cd-clat*sd
      costh=clat*cd+slat*sd
      theta=pi/2.-atan2(sinth,costh)
      return
      end

c-----------------------------------------------------------------

      subroutine magfdz(p,dp,theta,phi,r,lmax,g,dx,dy,dz,
     >x,y,z,h,f,i,d,
     >sd,cd)

c
c***************************************************************
c
c     j bloxham  8 nov 1982 & 11 oct 1983
c
c     modified version of dg13.sv.fort:magfd & jb62.sv.progs:magfds
c
c     gives field components at radius r
c
c***************************************************************
c
c
c     this version 16 jan 87
c
c     saves dx dy dz in computation
c
cc======================================================================

      implicit real*8 (a-h,o-z)
      dimension g(lmax*(lmax+2))
      dimension dx(lmax*(lmax+2)),dy(lmax*(lmax+2)),dz(lmax*(lmax+2))
      dimension p((lmax+1)*(lmax+2)/2),dp((lmax+1)*(lmax+2)/2)
      real*8 i


      b=6371.2/r
      x=0.
      y=0.
      z=0.
      sinth=sin(theta)
      if(abs(sinth).lt.1.e-10) sinth=1.e-10

      do 20 l=1,lmax

      l1=l+1
      bb=b**(l+2)
      k=l*l
      k1=(l*l1)/2+1

      dx(k)=dp(k1)*bb
      dy(k)=0.
      dz(k)=-p(k1)*l1*bb
      x=x+g(k)*dx(k)
      z=z+g(k)*dz(k)

      do 20 m=1,l

      t=float(m)*phi
      k=l*l+2*m-1
      k1=(l*l1)/2+m+1
      sint=sin(t)
      cost=cos(t)

      dxd = dp(k1)*bb
      dx(k) = dxd*cost
      dx(k+1) = dxd*sint
      x = x + (g(k)*dx(k)) + (g(k+1)*dx(k+1))

      dxy = m*p(k1)*bb/sinth
      dy(k) = dxy*sint
      dy(k+1) = -dxy*cost
      y = y + (g(k)*dy(k)) + (g(k+1)*dy(k+1))

      dzd = -l1*p(k1)*bb
      dz(k) = dzd*cost
      dz(k+1) = dzd*sint
      z = z + (g(k)*dz(k)) + (g(k+1)*dz(k+1))

20    continue
      
      xs = x
      x = x*cd + z*sd
      z = z*cd - xs*sd
   
      do 50 k=1,lmax*(lmax+2)
      dxk = dx(k)
      dzk = dz(k)
      dx(k) = dxk*cd + dzk*sd
      dz(k) = dzk*cd - dxk*sd
50    continue
    
      
      h=sqrt(x*x+y*y)
      f=sqrt(h*h+z*z)
      i=asin(z/f)
      d=atan2(y,x)

      return
      end
       
c---------------------------------------------------------------------

      subroutine plmbar(p,dp,z,lmax)
c
c  evaluates normalized associated legendre function p(l,m) as function of
c   z=cos(colatitude) using recurrence relation starting with p(l,l) 
c   and then increasing l keeping m fixed.  normalization is: 
c   integral(y(l,m)*y(l,m))=4.*pi, where y(l,m) = p(l,m)*exp(i*m*longitude),
c   which is incorporated into the recurrence relation. p(k) contains p(l,m)
c   with k=(l+1)*l/2+m+1; i.e. m increments through range 0 to l before 
c   incrementing l. routine is stable in single and double precision to
c   l,m = 511 at least; timing proportional to lmax**2
c   r.j.o'connell 7 sept. 1989

c   a.jackson 19 october 1989  code added at end:
c   (2) derivatives added and stored in dp(k)
c       using same arrangement as for p(k)
c
      implicit real*8(a-h,o-z)
      dimension p(*),dp(*)
c     --dimension of p, dp must be (lmax+1)*(lmax+2)/2 in calling program
      if (lmax.lt.0.or.abs(z).gt.1.d0) stop 'bad arguments'
c       --case for p(l,0) 
        pm2=1.d0
        p(1)=1.d0
        dp(1)=0.d0
        if (lmax .eq. 0) return
        pm1=z
        p(2)=dsqrt(3.d0)*pm1
        k=2
        do 4 l=2,lmax
          k=k+l
          plm=(dfloat(2*l-1)*z*pm1-dfloat(l-1)*pm2)/dfloat(l)
          p(k)=dsqrt(dfloat(2*l+1))*plm
          pm2=pm1
4         pm1=plm
c       --case for m > 0
        pmm = 1.d0
        sintsq = (1.d0-z)*(1.d0+z)
        fnum = -1.0d0
        fden = 0.0d0
        kstart = 1
        do 20 m =1 ,lmax
c         --case for p(m,m) 
          kstart = kstart+m+1
          fnum = fnum+2.0d0
          fden = fden+2.0d0
          pmm = pmm*sintsq*fnum/fden
          pm2 = dsqrt(dfloat(4*m+2)*pmm)
          p(kstart) = pm2
          if (m .eq. lmax) goto 100
c         --case for p(m+1,m)
          pm1=z*dsqrt(dfloat(2*m+3))*pm2
          k = kstart+m+1
          p(k) = pm1
c         --case for p(l,m) with l > m+1
          if (m .lt. (lmax-1)) then
           do 10 l = m+2,lmax
            k = k+l
            f1=dsqrt(dfloat((2*l+1)*(2*l-1))/dfloat((l+m)*(l-m)))
            f2=dsqrt(dfloat((2*l+1)*(l-m-1)*(l+m-1))
     &              /dfloat((2*l-3)*(l+m)*(l-m)))
            plm=z*f1*pm1-f2*pm2
            p(k) = plm
            pm2 = pm1
10          pm1 = plm
          endif
20        continue

100     continue

c       Gauss-Schmidt normalisation:
        k=1
        do 30 l=1,lmax
        fac=1.d0/dsqrt(dfloat(2*l+1))
        do 30 m=0,l
        k=k+1
        p(k)=p(k)*fac
30      continue

c       now find derivatives of p(z) wrt theta, where z=cos(theta)
        dp(2)=-p(3)
        dp(3)=p(2)
        k=3
        do 200 l=2,lmax
        
          k=k+1
c         treat m=0 and m=l separately
          dp(k)=-dsqrt(dfloat(l*(l+1))/2.d0)*p(k+1)
          dp(k+l)=dsqrt(dfloat(l)/2.d0)*p(k+l-1)
          do 300 m=1,l-1
            k=k+1
            fac1=dsqrt( dfloat( (l-m)*(l+m+1) ) )
            fac2=dsqrt( dfloat( (l+m)*(l-m+1) ) )
            if(m.eq.1)fac2=fac2*dsqrt(2.d0)
            dp(k)=0.5d0*( fac2*p(k-1) - fac1*p(k+1) )
300       continue
          k=k+1

200     continue
        return
        end

c****************************************      


