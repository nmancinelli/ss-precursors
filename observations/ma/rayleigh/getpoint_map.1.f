c** getpoint and plot them
      character*255 filename
      dimension phl(15,41252),gvl(15,41252)
      dimension xval(26,41252)
      dimension xval_prem(26)
      dimension BB(9000,30),t(40),freq(9000)
      dimension Bi(9000,30)
      dimension cref(15),uref(13)
      data pi/3.14159265358979d0/

      rad=pi/180.d0

      print *,'input long and lat (geographic)'
      read(*,*) p1,t1

      call blks(nblk)
      r0=6371.
      tt=t1
      pp=p1
      if (pp.lt.0) pp=pp+360.
      tt=geocen((90.-tt)*rad)/rad
      call fblk(r0,tt,pp,ib)

c** spline interpolation
      filename='xval.prem'
      open(8,file=filename,form='unformatted')
      read(8) nblk,dum
      nblk=41252
      do i=1,22
        read(8) idum,xval_prem(i),dum,dum
      enddo
      close(8)
      filename='xval.map'
      open(8,file=filename,form='unformatted')
      read(8) nblk,dum
      nblk=41252
      do i=1,22
        read(8) idum,(xval(i,j),j=1,nblk),dum,dum
        do j=1,nblk
          xval(i,j)=xval(i,j)+xval_prem(i)
        enddo
      enddo
      close(8)
      do i=1,25
        t(i)=(i-1)*2.0-2.
      enddo
      do i=1,21
        jj=0
        do ix=1,15
          x=(ix-1)*2.5+5.
          jj=jj+1
          if (i.eq.1) freq(jj)=x
          call bspline(t,i,x,bval,bpval,bp2val,bp3val,bival)
          BB(jj,i)=bval
          Bi(jj,i)=bival
        enddo
      enddo
      filename='point.1'
        open(8,file=filename)
        do j=1,jj
	  sum=0.
	  sump=0.
	  do i=1,21
	    sum=sum+BB(j,i)*xval(i,ib)
	    sump=sump+Bi(j,i)*xval(i,ib)
	  enddo
          sum=1./sum
          sump=sump+xval(22,ib)
          sump=sump/freq(j)
	  sump=1./sump
          if (freq(j).lt.7.5) sum=-99.
          if (freq(j).gt.35) sump=-99.
	  write(8,*) freq(j),sump,sum
	enddo
      end


      real function geocen(arg)
c input:
c   arg    = geographic colatitude (radians)
c output:
c   geocen = geocentric colatitude (radians)
c (n.b. fac=(1-f)**2)
      data pi2,fac/1.570796326794895,0.993305621334896/
      geocen=pi2-atan(fac*cos(arg)/amax1(1.e-30,sin(arg)))
      return
      end

      subroutine blks(nblk)
c*** routine to establish block parameters
      common/blocks/r1(100),r2(100),hsz(1800),bsz,
     +    nshell,mshell,mlat(1801),nlat
      character*256 filename
      data rad/57.29578/
      bsize=1.
      nlat=nint(180./bsize)
      bsize=180./nlat
      bsz=bsize
      nshell=1
      ib=0
      do 10 n=1,nshell
        rad1=6371.0
	rad2=6371.0
        tmax=0.
        r1(n)=rad1
        r2(n)=rad2
        mlat(1)=0
        do 30 ii=1,nlat                  
        tmin=tmax
        tmax=tmin+bsize
        th=0.5*(tmin+tmax)/rad
        s1=sin(th)
        mlon=max(nint(360./bsize*s1),1)
        hsize=360./mlon
        ib=ib+mlon
        if(n.eq.1) then
           hsz(ii)=hsize
           mlat(ii+1)=ib
        end if
   30   continue
      if(n.eq.1) mshell=ib
   10 continue
      nblk=ib
      close(12)
      return
      end

      subroutine fblk(r0,t0,p0,ib)
c*** routine to find block number given radial index and colat and long
      common/blocks/r1(100),r2(100),hsz(1800),bsz,
     +    nshell,mshell,mlat(1801),nlat
      do 10 i=1,nshell
      if(r0.ge.r1(i).and.r0.le.r2(i)) then
         n=i
         goto 20
      end if
   10 continue
      ib=-1
      return
   20 it=t0/bsz+1
      if(it.ge.1.and.it.le.nlat) then
        j1=(n-1)*mshell+mlat(it)+1
        ib=p0/hsz(it)+j1
      else
        ib=-1
      end if
      return
      end

      subroutine bspline(t,i,x,bval,bpval,bp2val,bp3val,bival)
c** calculate the bspline at point x with
c** knots t(i),t(i+1),t(i+2),t(i+3),t(i+4)
c** i should be specified before calling bspline
      dimension t(40)
      Del=t(i+1)-t(i)
      if (x.le.t(i).or.x.ge.t(i+4)) then
        bval=0.
        bpval=0.
        bp2val=0.
        bp3val=0.
        bival=0.
        if (x.ge.t(i+4)) bival=6.
      endif
      if (x.ge.t(i).and.x.lt.t(i+1)) then
        ex=(x-t(i))/Del
        bval=ex**3
        bpval=3.*ex**2
        bp2val=6.*ex
        bp3val=6.
        bival=0.25*ex**4
      endif
      if (x.ge.t(i+1).and.x.lt.t(i+2)) then
        ex=(x-t(i+1))/Del
        bval=-3.*ex*(ex*(ex-1.)-1.)+1.
        bpval=-3.*(ex*(3.*ex-2.)-1.)
        bp2val=-18.*ex+6.
        bp3val=-18.
        bival=-0.75*ex**4+ex**3+1.5*ex**2+ex+0.25
      endif
      if (x.ge.t(i+2).and.x.lt.t(i+3)) then
        ex=(x-t(i+2))/Del
        bval=3.*ex*ex*(ex-2.)+4.
        bpval=3.*ex*(3.*ex-4.)
        bp2val=18.*ex-12.
        bp3val=18.
        bival=0.75*ex**4-2*ex**3+4*ex+3.
      endif
      if (x.ge.t(i+3).and.x.lt.t(i+4)) then
        ex=(x-t(i+3))/Del
        bval=-(ex-1.)**3
        bpval=-3.*(ex-1.)**2
        bp2val=-6.*ex+6.
        bp3val=-6.
        bival=-0.25*ex**4+ex**3-1.5*ex**2+ex+5.75
      endif
      bval=bval/6.
      bpval=bpval/6./Del
      bp2val=bp2val/6./Del/Del
      bp3val=bp3val/6./Del/Del/Del
      bival=bival/6.*Del
      end
