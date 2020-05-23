      program LScoefs
c
c     A program to calculate Schmidt normalised SHA coefficients for a particular 
c     time from time dependent coefficients in model file LSMOD.1 or LSMOD.2.
c     Output of main field coefficients into file coefs.dat 
c     Output file format: 4 columns:
c     spherical harmonic degree l, order m, g_l^m, h_l^m
c
c    Copyright © 2018, Monika Korte, Helmholtz Cenntre Potsdam GFZ German 
c    Research Centre f or Geosciecnes, Potsdam, Germany

c    Licensed under the Apache License, Version 2.0 (the "License");
c    you may not use this file except in compliance with the License.
c    You may obtain a copy of the License at
c    http://www.apache.org/licenses/LICENSE-2.0

c   Unless required by applicable law or agreed to in writing, software
c   distributed under the License is distributed on an "AS IS" BASIS,
c   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
c   See the License for the specific language governing permissions and
c   imitations under the License.

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

      integer lmax,n,nsplt,jord
      parameter (lmax=10)
      parameter (nsplt=403)
      parameter (n=lmax*(lmax+2))

      real*8 g(n),dg(n),gt(n,nsplt),dgt(n,nsplt),h0
      real*8 time,tknts(nsplt+4),spl(nsplt),age
      real*8 tstartin,tendin
      integer lm,nm,i,k,j,l,m,nspl,nleft,itime,model
      character*50 outfile,infile

      data jord/4/

      write(*,*) 'Program to write SH coefficients for one epoch'
      write(*,*) 'from model LSMOD.1 or LSMOD.2 to file coefs.dat'
      write(*,*) '           --------------------                '
      write(*,*) 'Give input age (between 49.95 and 30 ka BP):'
      read(*,*) age
      write(*,*) 'Give model number:'
      write(*,*) '1 -- LSMOD.1'
      write(*,*) '2 -- LSMOD.2'
      read(*,*) model

      if (model.eq.1) then
         infile='LSMOD.1'
      else if (model.eq.2) then
         infile='LSMOD.2'
      else
         write(*,*) 'Error: Invalid model number!'
         stop
      endif

      if ((age.gt.49.95).or.(age.lt.30)) then
         write(*,*) 'Error: Age outside valid range!'
         write(*,*) 'Age must be between 49.95 and 30 ka.'
         stop
      endif

      time=-(age*1000.)+1950.
      outfile='coefs.dat'

c     read model
      open(7,file=infile)

      read(7,*) tstartin,tendin
      read(7,*) lm,nm,nspl,(tknts(i),i=1,nspl+4)
      read(7,*) gt
      close(7)


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
        

c     output coeffients:
      open(7,file=outfile)
      write(7,*) '  SH_degree    SH_order         g(nT)          h(nT)'
      h0=0.0
 5000 format(2I3,2F20.10)
 5001 format(2I3,2F20.10)
      do l=1,lmax
         m=0
         k=l*l
         write(7,*) l,m,g(k),h0
         do m=1,l
            k=l*l+2*m-1
            write(7,*) l,m,g(k),g(k+1)
         end do
      end do
      close(7)

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
        
c------------------------------------------------------------------

