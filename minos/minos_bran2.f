c**** program minos_bran ****
c  this program uses one input and two output files
c  the input file contains the model 
c  the output files are 1) a model listing + a summary of mode properties
c  and 2) a file for the eigenfunctions 
c  the program then asks for some control info described below
c
c structure of model file
c  card 1  :   title (80 chars)                (20a4 format)
c  card 2  :   ifanis,tref,ifdeck              (unformatted)
c           ifanis=1 for anisotropic model, 0 for isotropic
c           tref=ref period(secs) of model for dispersion correction.
c           if tref is .le. 0. no correction is made.
c           ifdeck=1 for card deck model, 0 for polynomial model.
c             *** card deck model ***
c  card 3  :   n,nic,noc                       (unformatted)
c           n=no of levels,nic=index of solid side of icb,noc=index of
c           fluid side of mcb.note that n must be .le. 223.
c  card 4-n:  r,rho,vpv,vsv,qkappa,qshear,vph,vsh,eta (f8.0,3f9.2,2f9.1,2f9.2,
c                                                           f9.5)
c           if isotropic vp=vpv,vs=vsv and vph,vsh,eta are unspecified.
c           if the q model is not specified then no disp. correction is made.
c           s.i. units are used,e.g. radius in metres and velocities in m/s.
c             *** polynomial model ***
c  card 3  :   nreg,nic,noc,rx                 (unformatted)
c           nreg is the number of regions in the model,nic and noc as before
c           and rx is the normalising radius for the polynomials.
c           rx is given in kms and is usually 6371.
c  card 4  :   nlay,r1,r2                      (unformatted)
c           nlay is the number of levels to be used in the region extending
c           from radius r1 to r2 (in kms).
c  card 5-9(iso) or card 5-12(ani) : coefs     (5f9.5)
c           5 sets of coefficients are required for each region of an isotropic
c           model and 8 sets for an anisotropic model.each polynomial can  be
c           up to a quartic in r/rx ( 5 coefficients) and are ordered thusly:
c           rho,vpv,vsv,qkappa,qshear,vph,vsh,eta. the coeffs are given in the
c           usual mixed seismological units (rho in g/cc,vel in km/s etc.)
c           conversion to s.i. is done by the program.
c           cards 4-9(iso) or 4-12(ani) are repeated for each region of the
c           model
c
c control parameters  (from screen)
c   eps,wgrav          
c           eps controls the accuracy of the runge-kutta integration scheme.
c           the relative accuracy of an eigen frequency will be 2-3 x eps.
c           it also controls the precision with which a root is found and
c           the minimum relative separation of two roots with the same angular 
c           order.
c           wgrav is the frequency in millihertz above which gravitational terms
c           are neglected-this gives about a factor of 3 increase in speed.
c   jcom       
c           jcom=1 radial modes, 2 toroidal modes, 3 spheroidal modes,
c           4 inner core toroidal modes. 
c   lmin,lmax,wmin,wmax,nmin,nmax            
c           lmin - lmax defines the range of angular orders to be computed.
c           if jcom=1 this is ignored. wmin - wmax defines the frequency range
c           to be computed (in millihertz)
c           nmin-nmax specifies the branch numbers to be computed -- nmin=0
c           is the fundamental mode
c                          
c  model listing
c    this is an ascii file which lists the model and mode properties
c    i.e. phase velocity in km/s,frequency,period,group velocity in km/s,
c    q and a parameter which is the ratio of kinetic to potential energy
c    minus 1 which should be small ( of order eps )if the eigenfunction is 
c    accurate. (you will probably see some degradation
c    in this parameter for strongly exponential modes such as stoneley modes ). 
c  
c  eigenfunction file
c****** file name "none" will suppress calculation of eigenfunctions
c    this is a fixed record length binary file with an entry for each mode
c    written as
c      write(ioeig) (abuf(i),i=1,nvec)
c    where nvec is 5+6*npts for spheroidal modes and 5+2*npts for toroidal
c    and radial modes. the first five words of abuf are n,l,frequecy,q and
c    group velocity. the rest of abuf contains w(1..npts) and wp(1..npts)
c    for toroidal modes and u(1..npts),up(1..npts),v(1..npts),vp(1..npts),
c    p(1..npts) and pp(1..npts) for spheroidal modes. these are as in 
c    woodhouse and dahlen 1978 except that w,wp,v and vp must be divided
c    by sqrt(l*(l+1)). the normalisation is such that 
c       frequency**2 times integral(rho*w*w*r*r) is 1 for toroidal modes
c       frequency**2 times integral(rho*(u*u+l*(l+1)*v*v)*r*r) is 1 for 
c    spheroidal modes
c    the model has been normalised such that a density of 5515 mg/m**3 = 1 ;
c    pi*G=1 where G is the grav constant (6.6723e-11) and the radius of the
c    earth (rn=6371000 m) is 1. these normalizations result in 
c      acceleration normalisation = pi*G*rhobar*rn
c      velocity normalisation     = rn*sqrt(pi*G*rhobar)= vn
c      frequency normalization    = vn/rn
c
      character*256 filnam
      print *,'input model file:'
      read(5,100) filnam
  100 format(a256)
      open(7,file=filnam,status='old',form='formatted',iostat=iret)
      print *,'output file:'
      read(5,100) filnam
      open(8,file=filnam,form='formatted',iostat=iret)
      call model(7,8)
      close(7)
      print *,'eigenfunction file (output):'
      read(5,100) filnam
      ifreq=1
      if(filnam(1:4).eq.'none') ifreq=0
      open(3,file=filnam,form='unformatted',iostat=iret)
      call timer
      call wbran(8,3,ifreq)
      call timer
      close(8)
      close(3)
      stop
      end 

      subroutine timer
      save oldtim
      real tarray(2)
      data if1/1/
      if(if1.eq.0) goto 5
      if1=0
      oldtim=etime(tarray)
      return
    5 continue
      tnow=etime(tarray)
      diftim=tnow-oldtim
      oldtim=tnow
      print *,'elapsed seconds =',diftim
      return
      end

      subroutine wbran(iout,ioeig,ifreq)
c*** chases mode branches
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      common/stdep/fick,ls
      dimension we(2),de(2),ke(2)
      character*2 ichar(4)
      data ichar/' s',' t',' s',' c'/
      cmhz=pi/500.d0
      modtot=0
      print *,'enter eps and wgrav'
      read(*,*) eps,wgrav
      eps=max(eps,1.d-12)
      wgrav=wgrav*cmhz
      write(iout,100) eps,eps,wgrav
  100 format(/,'integration precision =',g12.4,'  root precision =',
     +   g12.4,'  gravity cut off =',g12.4,' rad/s',///,6x,'mode',
     +   8x,'phs vel',7x,'w(mhz)',10x,'t(secs)',6x,'grp vel(km/s)',
     +   8x,'q',13x,'raylquo',/)
      call steps(eps)
      print *,'enter jcom (1=rad;2=tor;3=sph;4=ictor)'
      read(*,*) jcom
      if(jcom.lt.1.or.jcom.gt.4) jcom=3
      print *,'enter lmin,lmax,wmin,wmax,nmin,nmax'
      read(*,*) lmin,lmax,wmin,wmax,normin,normax
      wmin=max(wmin,0.2d0)
      wmin=wmin*cmhz
      wmax=wmax*cmhz
      if(lmin.le.0) lmin=1
      if(jcom.eq.1) then
        lmin=0
        lmax=0
      end if
      normin=max(normin,0)
      normax=max(normax,normin)
      do 50 nor=normin,normax
      modbr=0
      wlast=wmin
      wtry=min(wmin+0.5d0*cmhz,wmax)
      bm=0.d0
      do 10 l=lmin,lmax
      if(wlast.ge.wmax) goto 45
      if(nor.eq.0.and.l.eq.1) goto 10
      if(nor.eq.1.and.l.eq.1) goto 10
      nord=nor
      nup=nor
      ndn=nor-1
      if(l.eq.1) then
        nup=ndn
        ndn=ndn-1
      end if
      knsw=1
      maxo=5
      fl=l
      fl1=fl+1.d0
      fl2=fl+fl1
      fl3=fl*fl1
      sfl3=sqrt(fl3)
      we(1)=wlast
      call detqn(we(1),ke(1),de(1),0)
      we(2)=wmax
      if(wtry.ne.0.d0) we(2)=min(wtry,wmax)
      call detqn(we(2),ke(2),de(2),0)
      if(ke(2).lt.nup) then
         we(2)=wmax
         call detqn(we(2),ke(2),de(2),0)
         if(ke(2).lt.nup) goto 45
      end if
c*** bracket this mode ***
      ktry2=0
   15 if(ke(1).eq.ndn.and.ke(2).eq.nup) goto 40
      if(abs(we(1)-we(2))/we(1).lt.eps) then
         print *,'unable to isolate mode : ',nord,l
c*** come here if there is a problem in the counter -- what occasionally happens is that
c*** the counter will skip over two closely spaced zero crossings and will be too low by 2
c*** -- this will cause a problem if ke(1) is wrong because then we(1) will be too high
   16    print *,we(1),we(2),ke(1),ke(2),ndn,nup
         ktry2=ktry2+1
         if(ktry2.gt.20) goto 10
         we(1)=wlast-ktry2*1.d-3
         if(wtry.gt.0.d0) then
           we(2)=min(wtry+ktry2*1.d-3,wmax)
         else
           we(2)=wmax
         end if
         call detqn(we(1),ke(1),de(1),0)
         if(ke(1).gt.ndn) goto 16
         call detqn(we(2),ke(2),de(2),0)
         if(ke(2).lt.nup) goto 16
         goto 15
      end if
      wx=0.5d0*(we(1)+we(2))
      call detqn(wx,kx,dx,0)
      if(kx.le.ndn) then
        we(1)=wx
        ke(1)=kx
        de(1)=dx
      else
        we(2)=wx
        ke(2)=kx
        de(2)=dx
      end if
      goto 15
   40 knsw=0
      maxo=8
c*** find root ***
      call rbrent(we(1),de(1),we(2),de(2),eps,wc)
      tcom=2.d0*pi/wc
      wmhz=1000.d0/tcom
      cvel=wc*rn/(l+0.5)/1000.
      if(ifreq.eq.1) then
c*** come here if eigenfunction computed
        call detqn(wc,knt,fb,ifreq)
        wdiff=(wc-wray*wn)/wc
        gcom=vn*cg/1000.d0
        qmod=0.d0
        if(qinv.gt.0.d0) then
           qmod=1.d0/qinv
c*** note that this is physical dispersion correction to group velocity
           gcom=gcom*(1.d0+qinv/pi)
           cg=cg*(1.d0+qinv/pi)
        end if
        write(iout,200) nord,ichar(jcom),l,cvel,wmhz,tcom,gcom,qmod,
     +      wdiff,ls,fb,fick,wray*wn/cmhz
  200   format(i5,a2,i5,6g16.7,i5,3g16.7)
        if(abs(wdiff).gt.1.d-4.and.ls.lt.n) then
          print *,'rejecting mode : ' 
          print 200, nord,ichar(jcom),l,cvel,wmhz,tcom,gcom,qmod,
     +      wdiff,ls,fb,fick,wray*wn/cmhz
        else
         call modout(wc,qmod,gcom,ioeig)
        end if
      else
c*** come here if no eigenfunction -- group vel is estimated from last two frequencies
        cg=cvel
        if(bm.ne.0.d0) cg=(wc-bm)/wn
        gcom=vn*cg/1000.d0
        write(iout,200) nord,ichar(jcom),l,cvel,wmhz,tcom
        bm=wc
      end if
      wlast=wc
      if(l.eq.lmin) wmin=max(wc-.001,wmin)
      if(jcom.eq.1) then
        wtry=wc+0.01
      else
        wtry=wc+2.d0*cg*wn
      end if
      modbr=modbr+1
   10 continue
   45 modtot=modtot+modbr
      if(modbr.eq.0) return
      print *,nor,modtot,modbr
   50 continue
      return
      end

      subroutine rbrent(a,fa,b,fb,eps,ans)
c*** finds root using quadratic interpolation (brent's method from numerical recipes)
      parameter (itmax=100)
      implicit real*8(a-h,o-z)
c*** now refine root
      ans=0.d0
      if(fa*fb.ge.0.d0) return
      c=b
      fc=fb
      do iter=1,itmax
         if(fb*fc.gt.0.d0) then
            c=a
            fc=fa
            d=b-a
            e=d
         end if
         if(abs(fc).lt.abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         end if
         xm=0.5d0*(c-b)
         t=2.d0*abs(eps*b)
c**** here is convergence
         if(abs(xm).lt.t) then
            ans=b
            return
         end if
         if(abs(e).gt.t.and.abs(fa).gt.abs(fb)) then
c*** attempt inverse quadratic interpolation
           s=fb/fa
           if(a.eq.c) then
              p=2.d0*xm*s
              q=1.d0-s
           else
              q=fa/fc
              r=fb/fc
              p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
              q=(q-1.d0)*(r-1.d0)*(s-1.d0)
           end if
           if(p.gt.0.d0) q=-q
           p=abs(p)
           if(2.d0*p.lt.min(3.d0*xm*q-abs(t*q),abs(e*q))) then
c*** accept interpolation
             e=d
             d=p/q
           else
c*** interpolation failed -- use bisection
             d=xm
             e=d
           end if
         else
c*** bounds decreasing too slowly -- use bisection
             d=xm
             e=d
         end if
         a=b
         fa=fb
         if(abs(d).gt.t) then
             b=b+d
         else
             b=b+sign(t,xm)
         end if
         call detqn(b,knt,fb,0)
      enddo
      print *,'max interations exceeded'
      stop
      end

      subroutine modout(wcom,qmod,gcom,ioeig)
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/eifx/a(14,mk),dum(mk)
      common/bufX/nn,ll,ww,qq,gc,buf(6*mk)
      common/stdep/fick,ls
      dimension abuf(6*mk+5)
      equivalence (nn,abuf)
      if(ls.gt.n) return
      nn=nord
      ll=l
      ww=wcom
      qq=qmod
      gc=gcom

      print *,w,wsq,cg,fl1,n,nsl
      if (ll.eq.30) then
        do i=n,n-30,-1
          print *,a(1,i)
        enddo
      endif

      if(jcom.eq.3) then
        nvec=6*n+5
c*** vector is U,U',V,V',P,P' -- last two will be zero if kg=0
        do i=1,n
          buf(i)=a(1,i)
          buf(i+n)=a(2,i)
          buf(i+n+n)=a(3,i)
          buf(i+3*n)=a(4,i)
          buf(i+4*n)=a(5,i)
          buf(i+5*n)=a(6,i)
        enddo
      else
        nvec=2*n+5
        if(jcom.eq.2) then
          do i=1,noc
            a(1,i)=0.d0
            a(2,i)=0.d0
          enddo
        end if
        do i=1,n
          buf(i)=a(1,i)
          buf(i+n)=a(2,i)
        enddo
      end if
      write(ioeig) (abuf(i),i=1,nvec)
      return
      end

      subroutine model(iin,iout)
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      character*4 ititle(20)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/vpv(mk),vph(mk),vsv(mk),vsh(mk),eta(mk),wrk(mk*10)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      data bigg,tau/6.6723d-11,1.d3/,rhobar/5515.d0/
      pi=3.14159265358979d0
      read(iin,100) (ititle(i),i=1,20)
  100 format(20a4)
      read(iin,*) ifanis,tref,ifdeck
      if(ifdeck.eq.0) go to 1000
c*** card deck model ***
      read(iin,*) n,nic,noc
      if(ifanis.ne.0) then
c           read(iin,*) (r(i),rho(i),vpv(i),vph(i),vsv(i),vsh(i),
c     +     eta(i),qshear(i),qkappa(i),i=1,n)
c** for the STW105.txt
           read(iin,*) (r(i),rho(i),vpv(i),vsv(i),
     +    qkappa(i),qshear(i),vph(i),vsh(i),eta(i),i=1,n)
      else
           read(iin,*) (r(i),rho(i),vpv(i),vsv(i),
     +     qkappa(i),qshear(i),i=1,n)
      endif
      go to 2000
c*** polynomial model ***
 1000 read(iin,*) nreg,nic,noc,rx
      rx=rx*tau
      n=0
      knt=0
      jj=5
      if(ifanis.ne.0) jj=8
      do nn=1,nreg
        read(iin,*) nlay,r1,r2
        r1=r1*tau
        r2=r2*tau
        dr=(r2-r1)/float(nlay-1)
        do i=1,nlay
          n=n+1
          r(n)=r1+dr*float(i-1)
        enddo
        do j=1,jj
          read(iin,110) (wrk(i),i=1,5)
  110     format(5f9.5)
          if(ifanis.eq.2) then
            do i=1,nlay
              ind=knt+i
              rt=r(ind)/rx
              val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
              if(j.eq.1) rho(ind)=val*tau
              if(j.eq.2) vph(ind)=val*tau
              if(j.eq.3) vsv(ind)=val*tau
              if(j.eq.4) qkappa(ind)=val
              if(j.eq.5) qshear(ind)=val
              if(j.eq.6) vpv(ind)=sqrt(val)*vph(ind)
              if(j.eq.7) vsh(ind)=sqrt(val)*vsv(ind)
              if(j.eq.8) eta(ind)=val
            enddo
          else
            do i=1,nlay
              ind=knt+i
              rt=r(ind)/rx
              val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
              if(j.eq.1) rho(ind)=val*tau
              if(j.eq.2) vpv(ind)=val*tau
              if(j.eq.3) vsv(ind)=val*tau
              if(j.eq.4) qkappa(ind)=val
              if(j.eq.5) qshear(ind)=val
              if(j.eq.6) vph(ind)=val*tau
              if(j.eq.7) vsh(ind)=val*tau
              if(j.eq.8) eta(ind)=val
            enddo
          end if
        enddo
        knt=knt+nlay
      enddo
 2000 if(ifanis.eq.0) then
      do i=1,n
        vph(i)=vpv(i)
        vsh(i)=vsv(i)
        eta(i)=1.d0
      enddo
      end if
      if(iout.ge.0) then
c*** write out model ***
        write(iout,900) (ititle(k),k=1,20),tref
  900   format(1x,20a4,' ref per =',f6.1,' secs',///,2x,'level',
     1   4x,'radius',8x,'rho',9x,'vpv',9x,'vph',9x,'vsv',
     2   9x,'vsh',9x,'eta',9x,'qmu ',8x,'qkap',/)
        write(iout,905) (i,r(i),rho(i),vpv(i),vph(i),vsv(i),vsh(i),
     1   eta(i),qshear(i),qkappa(i),i=1,n)
  905   format(3x,i3,f12.1,5f12.2,f12.5,2f12.2)
      end if
c*** normalise and spline ***
      rn=r(n)
      gn=pi*bigg*rhobar*rn
      vn2=gn*rn
      vn=sqrt(vn2)
      wn=vn/rn
      do i=1,n
        r(i)=r(i)/rn
        if(i.gt.1.and.dabs(r(i)-r(i-1)).lt.1.d-7) r(i)=r(i-1)
        if(qshear(i).gt.0.d0) qshear(i)=1.d0/qshear(i)
        if(qkappa(i).gt.0.d0) qkappa(i)=1.d0/qkappa(i)
        rho(i)=rho(i)/rhobar
        acon(i)=rho(i)*vph(i)*vph(i)/vn2
        ccon(i)=rho(i)*vpv(i)*vpv(i)/vn2
        lcon(i)=rho(i)*vsv(i)*vsv(i)/vn2
        ncon(i)=rho(i)*vsh(i)*vsh(i)/vn2
        fcon(i)=eta(i)*(acon(i)-2.d0*lcon(i))
        fmu(i)=(acon(i)+ccon(i)-2.d0*fcon(i)+5.d0*ncon(i)
     +   +6.d0*lcon(i))/15.d0
        flam(i)=(4.d0*(acon(i)+fcon(i)-ncon(i))+ccon(i))/9.d0
     +    -2.d0*fmu(i)/3.d0
        rat=4.d0*fmu(i)/(3.d0*(flam(i)+2.d0*fmu(i)))
        xlam(i)=((1.d0-rat)*qkappa(i)-.5d0*rat*qshear(i))/(1.d0
     +    -1.5d0*rat)
        xa2(i)=(1.d0-rat)*qkappa(i)+rat*qshear(i)
      enddo
      call drspln(1,n,r,rho,qro,wrk)
      call grav(g,rho,qro,r,n)
      call drspln(1,n,r,g,qg,wrk)
      call drspln(1,n,r,fcon,fspl,wrk)
      call drspln(1,n,r,lcon,lspl,wrk)
      if(ifanis.ne.0) then
        call drspln(1,n,r,acon,aspl,wrk)
        call drspln(1,n,r,ccon,cspl,wrk)
        call drspln(1,n,r,ncon,nspl,wrk)
      end if
      nsl=n+1
   65 nsl=nsl-1
      if(vsv(nsl).le.0.d0) go to 65
      nicp1=nic+1
      nocp1=noc+1
      nslp1=nsl+1
      tref=0.5d0*tref/pi
      return
      end

      subroutine detqn(wdim,knt,det,ifeif)
c**** supevises the integration of the equations,it returns the value
c**** of the secular determinant as det and the count of zero crossings.
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),inorm(mk),jjj(mk)
      common/upx/up(14,mk),iup(mk)
      common/stdep/fick,ls
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension ass(14),vf(mk),zi(4),f(14)
      iback=0
      w=wdim/wn
      wsq=w*w
      iexp=0
      kount=0
      kg=0
      fct=0.d0
      if(tref.gt.0.d0) fct=2.d0*dlog(tref*wdim)/pi
      goto (2,3,1,3),jcom
    1 if(wdim.le.wgrav) kg=1
      nvefm=2+kg*3
      nvesm=5+kg*9
      call sdepth(wdim,ls)
      if(ls.le.2) then
c*** would like to have first node be below turning point -- constant should be = V_s/a=6.d-4 in following
c*** not really adequate for low l modes -- maybe because solint not good enough?
c        r10=4.5d-4*(fl+.5d0)/wdim
c*** this works better in practice  -- possibly because there are integration issues when r10 is very small
        fc=1.40*fl-1.40
        r10=max(fc,0.5)/6371.d0
        if(r10.lt.r(2)) then
          r(1)=r10
          g(1)=rho(1)*r(1)*1.333333333333333d0
          ls=1
        end if
      end if
      if(ls.le.nic) then
c*** propagate through inner core ***
        call spsm(ls,nvesm,ass)
        call sprpmn(ls,nic,ass,vf,nvesm,iexp)
        call sfbm(ass,kg,iback)
      end if
      if(ls.le.noc) then
c*** propagate through outer core ***
        is=max(ls,nicp1)
        if(is.eq.ls) call fpsm(ls,nvefm,ass)
        call fprpmn(is,noc,ass,vf,nvefm,iexp)
        call fsbm(ass,kg,iback)
      end if
      is=max(ls,nocp1)
      if(is.eq.ls) call spsm(ls,nvesm,ass)
c*** propagate through mantle ***
      call sprpmn(is,nsl,ass,vf,nvesm,iexp)
      if(nsl.eq.n) then
        dnorm=a(1,nsl)*a(1,nsl)
        do i=2,nvesm
          dnorm=dnorm+a(i,nsl)*a(i,nsl)
        enddo
        det=a(5,nsl)/sqrt(dnorm)
      else
        call sfbm(ass,kg,iback)
c*** propagate through ocean ***
        call fprpmn(nslp1,n,ass,vf,nvefm,iexp)
        if(kg.eq.0) then
           det=a(2,n)/sqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
        else
          det=a(5,n)/sqrt(a(1,n)**2+a(2,n)**2+a(3,n)**2+
     +   a(4,n)**2+a(5,n)**2)
        end if
      end if
      if(ls.gt.noc) det=-det
      if(knsw.eq.1) then
        if(ls.gt.noc) kount=kount-2
        irem=mod(kount,2)
        if(irem.eq.0.and.det.lt.0.d0) kount=kount+1
        if(irem.ne.0.and.det.gt.0.d0) kount=kount+1
        knt=kount
      end if
      if(ifeif.eq.0) then
         r(1)=0.d0
         g(1)=0.d0
         return
      end if
cc***
c      do i=ls,n
c        do j=1,14
c          f(j)=a(j,i)
c        enddo
c      aaa=4.d0*r(i)/(2.d0*l+1.d0)
c      a1=(f(6)/f(1))
cc      a2=(f(7)/f(2))/aaa
cc      a3=(f(8)/f(3))/aaa
cc      a4=(f(9)/f(4))/aaa
cc      a5=(f(10)/f(5))/aaa
c      a2=a1+f(13)*f(14)/(f(1)*f(2))
c      a3=a1+f(13)*f(13)/(f(1)*f(3))
c      a4=a1-f(14)*f(14)/(f(1)*f(4))
c      a5=a1-(f(12)*f(14)+f(11)*f(13))/(f(1)*f(5))
cc      a3=a2+f(12)*f(13)/(f(2)*f(3))
cc      a4=a2-f(11)*f(14)/(f(2)*f(4))
cc      a5=a2-f(11)*f(12)/(f(2)*f(5))
c      print 666,a2*f(2)/f(7),a3*f(3)/f(8),a4*f(4)/f(9),a5*f(5)/f(10)
cc      print 666,a1*aaa,a2*aaa,a3*aaa,a4*aaa,sss,sss2
c  666 format(6g14.6)
c      enddo
c***

c*** eigenfunction calculation for spheroidal modes ***
c*** following is for a potential problem mode
c*** we have found it best to use the woodhouse 1988 algorithm toe get the eigen function down to
c*** either the cmb or icb for such modes
      fcut=1.212*l-18.18
      fdim=500.*wdim/pi
c*** above this frequency, the regular algorithm always works -- if below, copy upgoing solution
      if(fdim.lt.fcut.and.abs(det).gt..02.and.ls.le.noc) then
c**** catch high l stoneley modes which have negligible surface displacement
c*** beyond about l=500, the mode eigenfunctions are not very precise -- but shouldn't have any
c*** surface displacement.
        if(ls.gt.nic.and.l.gt.500) then
           cg=0.
           wray=0.
           qinv=0.
           ls=n+1
           r(1)=0.d0
           g(1)=0.d0
           return
        end if
c*** following is a conservative calculation of the frequency below which the algorithm can be imprecies
c*** -- but the IC mode has no significant surface displacement -- if kg=0. most modes are accurately computed
        if(ls.lt.nic) then
c          fcut2=2.222*l-244
          fcut2=2.667*l-333
          fdim=500.*wdim/pi
          if(fdim.lt.fcut2) then
            cg=0.
            wray=0.
            qinv=0.
            ls=n+1
            r(1)=0.d0
            g(1)=0.d0
            return
          end if
        end if
c*** copy original upgoing solution
        do i=ls,n
          iup(i)=inorm(i)
          do j=1,nvesm
            up(j,i)=a(j,i)
          enddo
        enddo
      end if
c*** now do regular eigenfunction case
      imtch=0
      iback=1
      jexp=0
      nbakf=1+kg*3
      nbaks=4+kg*10
      do i=1,nbaks
        ass(i)=0.d0
      enddo
      if(n.ne.nsl) then
        if(kg.eq.0) then
          ass(1)=dsign(1.d0,a(1,n))
        else
          asi1=a(3,n)*a(3,n)
          asi2=a(4,n)*a(4,n)
          if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
          if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
        end if
        call fprpmn(n,nslp1,ass,vf,nbakf,jexp)
        call fsbm(ass,kg,iback)
      else
        asi1=a(3,n)*a(3,n)
        asi2=a(4,n)*a(4,n)
        if(kg.ne.0) then
          asi1=asi1+a(12,n)*a(12,n)
          asi2=asi2+a(11,n)*a(11,n)
        end if
        if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
        if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
      end if
      nto=max(ls,nocp1)
      call sprpmn(nsl,nto,ass,vf,nbaks,jexp)
      if(nto.eq.ls) goto 90
      call sfbm(ass,kg,iback)
      nto=max(ls,nicp1)
      call fprpmn(noc,nto,ass,vf,nbakf,jexp)
      if(nto.eq.ls) goto 90
      call fsbm(ass,kg,iback)
      nto=max(ls,1)
      call sprpmn(nic,nto,ass,vf,nbaks,jexp)
  90  if(fdim.lt.fcut.and.abs(det).gt..02.and.ls.le.noc)  then
        call remedy(ls,imtch)
      end if
      call eifout(ls,fick,imtch)
      r(1)=0.d0
      g(1)=0.d0
      return
c
c*** radial modes ***
    2 ls=2
c need to start close to the center to get the integrals correct
      r10=5.d-5
      if(r10.lt.r(2)) then
        r(1)=r10
        g(1)=rho(1)*r(1)*1.333333333333333d0
        ls=1
      end if
      call rps(ls,ass)
      call rprop(ls,n,ass)
      det=a(2,n)/sqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
      knt=kount-1
      if(ifeif.eq.0) return
      call radint(ls,ass)
      call erprop(ls,n,ass)
      rnrm=1.d0/(w*sqrt(ass(3)))
      qinv=ass(4)/(wsq*ass(3))
      wray=sqrt(ass(5)/ass(3))
      cg=0.d0
      do i=ls,n
        ff=fcon(i)*(1.d0+xlam(i)*fct)
        cc=ccon(i)*(1.d0+xa2(i)*fct)
        a(2,i)=(a(2,i)-2.d0*ff*a(1,i)/r(i))/cc
        a(1,i)=a(1,i)*rnrm
        a(2,i)=a(2,i)*rnrm
      enddo
      a(1,1)=0.d0
      r(1)=0.d0
      g(1)=0.d0
      return
c
c*** toroidal modes ***
    3 if(jcom.eq.2) then
        ls=nocp1
        n2=nsl
        ass(1)=1.d0
        ass(2)=0.d0
      else
        ls=2
        r10=5.d-5
        if(r10.lt.r(2)) then
           r(1)=r10
           ls=1
        end if
        n2=nic
      end if
      q=0.d0
      call startl(ls,n2,fmu,ls,q)
      if(ls.ne.nocp1) call tps(ls,ass)
      if(ifeif.ne.1) then
        call tprop(ls,n2,ass)
        det=ass(2)/dsqrt(ass(1)*ass(1)+ass(2)*ass(2))
        r(1)=0.d0
        if(knsw.ne.1) return
        if(jcom.eq.4) then
           knt=kount-1
           return
        end if
        knt=kount-2
        if(l.eq.1) return
        knt=knt+1
        irem=mod(knt,2)
        if(irem.eq.0.and.det.lt.0.d0) return
        if(irem.ne.0.and.det.gt.0.d0) return
        knt=knt+1
        return
      end if
      do ind=3,6
        ass(ind)=0.d0
      enddo
      call etprop(ls,n2,ass)
      det=ass(2)/dsqrt(ass(1)*ass(1)+ass(2)*ass(2))
      do i=ls,n2
        a(2,i)=a(1,i)/r(i)+a(2,i)/(lcon(i)*(1.d0+qshear(i)*fct))
      enddo
      rnrm=1.d0/(w*dsqrt(ass(3)))
      cg=(fl+0.5d0)*ass(4)/(w*ass(3))
      qinv=ass(5)/(wsq*ass(3))
      wray=dsqrt(ass(6)/ass(3))
      if(ls.ne.1) then
        do i=1,ls-1
          a(1,i)=0.d0
          a(2,i)=0.d0
        enddo
      end if
      do i=ls,n2
          a(1,i)=a(1,i)*rnrm
          a(2,i)=a(2,i)*rnrm
      enddo
      r(1)=0.d0
      a(1,1)=0.
      return
      end

      subroutine sprpmn(jf,jl,f,h,nvesm,iexp)
c*** propagate a minor vector in a solid region from level jf to jl ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension f(*),h(nvesm,*),s(14),fp(14),rne(6)
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 call norm(f,iexp,nvesm)
      if(iback.ne.0) then
        inorm(i)=inorm(i)+iexp
        if(kg.ne.0) then
          t1=f(4)+f(8)
          t2=t1+f(4)
          t1=t1+f(8)
          t3=f(8)-f(4)
          rne(1)=ar(6,i)*f(10)-ar(14,i)*f(9)+ar(13,i)*t3
     1          -ar(1,i)*f(7)-ar(7,i)*f(6)+ar(8,i)*f(5)
     2          +ar(12,i)*f(3)-ar(2,i)*f(2)+ar(3,i)*f(1)
          rne(2)=ar(6,i)*f(13)+ar(14,i)*t2+ar(13,i)*f(12)
     1          -ar(1,i)*f(11)-ar(9,i)*f(6)-ar(7,i)*f(5)
     2          +ar(11,i)*f(3)-ar(4,i)*f(2)-ar(2,i)*f(1)
          rne(3)=ar(6,i)*f(14)-ar(7,i)*t1-ar(8,i)*f(12)
     1          +ar(13,i)*f(11)-ar(9,i)*f(9)+ar(14,i)*f(7)
     2          +ar(10,i)*f(3)+ar(11,i)*f(2)+ar(12,i)*f(1)
          rne(4)=ar(14,i)*f(14)+ar(7,i)*f(13)+ar(12,i)*f(12)
     1          -ar(2,i)*f(11)-ar(9,i)*f(10)-ar(11,i)*t3
     2          +ar(4,i)*f(7)+ar(10,i)*f(5)+ar(5,i)*f(1)
          rne(5)=ar(13,i)*f(14)+ar(8,i)*f(13)-ar(12,i)*t2
     1          -ar(3,i)*f(11)+ar(7,i)*f(10)-ar(11,i)*f(9)
     2          -ar(2,i)*f(7)+ar(10,i)*f(6)+ar(5,i)*f(2)
          rne(6)=ar(1,i)*f(14)+ar(13,i)*f(13)-ar(2,i)*t1
     1          -ar(3,i)*f(12)+ar(14,i)*f(10)-ar(4,i)*f(9)
     2          -ar(11,i)*f(6)-ar(12,i)*f(5)+ar(5,i)*f(3)
        else
          rne(1)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
          rne(2)=-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
          rne(3)=-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
          rne(4)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
        end if
        do jj=1,6
          ar(jj,i)=rne(jj)
        enddo
      else
        inorm(i)=iexp
        do j=1,nvesm
          ar(j,i)=f(j)
        enddo
      end if
      if(i.eq.jl) return
      i=i+jud
      x=y
      y=r(i)
      if(y.eq.x) goto 5
      iq=min(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=sqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=.5d0*(alfsq+betasq+delsq)-fl3/(x*x)
      qt=sqrt(abs(fksq))+sqrt(abs(fksq-delsq))+2.d0/r(iq)
      q=qt+float(kg)*sfl3/x
      del=float(jud)*step(maxo)/q
      dxs=0.d0
   15   do j=1,nvesm
          s(j)=f(j)
        enddo
        y=x+del
   30   if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q)
        dxs=dx
        do ni=1,in
          z=x+c(ni)
          call derms(iq,z,f,h(1,ni),qff,qll,qaa)
          call rkdot(f,s,h,nvesm,ni)
        enddo
        if(knsw.eq.1) then
          call zknt(s,h,f,x,y,1,isplit)
          if(isplit.eq.1) then
            if(abs(dx).lt.1.d-6) goto 40
            do j=1,nvesm
              f(j)=s(j)
            enddo
            dx=dx*0.5d0
            y=x+dx
            goto 30
          end if
        end if
   40   x=y
        if(y.ne.r(i)) goto 15
      go to 5
      end

      subroutine fprpmn(jf,jl,f,h,nvefm,iexp)
c*** propagate the minor vector in a fluid region from level jf to jl ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension f(*),h(nvefm,1),s(5),fp(5)
      if(nvefm.eq.1) then
        do i=jl,jf
          inorm(i)=inorm(i)+iexp
          do j=1,2
            ar(j,i)=ar(j,i)*f(1)
          enddo
        enddo
        return
      end if
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 call norm(f,iexp,nvefm)
      if(iback.ne.0) then
        inorm(i)=inorm(i)+iexp
        rne2   =-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
        ar(1,i)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
        rne3   =-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
        ar(4,i)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
        ar(2,i)=rne2
        ar(3,i)=rne3
      else
        inorm(i)=iexp
        do j=1,nvefm
          ar(j,i)=f(j)
        enddo
      end if
      if(i.eq.jl) return
      i=i+jud
      x=y
      y=r(i)
      if(y.eq.x) goto 5
      iq=min(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=sqrt(abs(alfsq-fl3/(x*x)))+1.d0/r(iq)+float(kg)*sfl3/x
      del=float(jud)*step(maxo)/q
      dxs=0.d0
   15   do j=1,nvefm
          s(j)=f(j)
        enddo
        y=x+del
   30   if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q)
        dxs=dx
        do ni=1,in
          z=x+c(ni)
          call dermf(iq,z,f,h(1,ni),qff)
          call rkdot(f,s,h,nvefm,ni)
        enddo
        if(knsw.eq.1) then
          call zknt(s,h,f,x,y,0,isplit)
          if(isplit.eq.1) then
            if(abs(dx).lt.1.d-6) goto 40
            do j=1,nvefm
              f(j)=s(j)
            enddo
            dx=dx*0.5d0
            y=x+dx
            goto 30
          end if
        end if
   40   x=y
        if(y.ne.r(i)) go to 15
      go to 5
      end

      subroutine derms(iq,z,f,fp,qff,qll,qaa)
c*** calculates minor vector derivative (fp) in a solid ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 nn,ll,lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension f(*),fp(*)
      t=z-r(iq)
      if(t.eq.0.d0) then
        ro=rho(iq)
        gr=g(iq)
        ff=fcon(iq)*qff
        ll=lcon(iq)*qll
        nn=ncon(iq)*qll
        cc=ccon(iq)*qaa
        aa=acon(iq)*qaa
      else
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
        ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
        if(ifanis.eq.0) then
          nn=ll
          cc=ff+ll+ll
          aa=cc
        else
          nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
          aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
        end if
      end if
      zr=1.d0/z
      sfl3z=sfl3*zr
      rogr=ro*gr
      c11=1.d0/cc
      c22=1.d0/ll
      dmg=aa-nn-ff*ff*c11
      zdmg=zr*dmg
      t11=-2.d0*ff*zr*c11+zr
      t12=sfl3z*ff*c11
      t21=-sfl3z
      t22=zr+zr
      s22=-ro*wsq
      s11=s22+4.d0*zr*(zdmg-rogr)
      s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn)
      s12=sfl3z*(rogr-zdmg-zdmg)
      if(kg.eq.0) then
        s11=s11+4.d0*ro*ro
        if(iback.eq.1) then
          fp(1)=t22*f(1)-t21*f(2)-c22*f(3)
          fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
          fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
          fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
        else
          b11=t11+t22
          b33=t11-t22
          fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
          fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
          fp(3)=s22*f(1)-2.d0*t12*f(2)+b33*f(3)+c11*f(5)
          fp(4)=-s11*f(1)+2.d0*t21*f(2)-b33*f(4)-c22*f(5)
          fp(5)=-2.d0*s12*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
        endif
      else
        t31=-4.d0*ro
        t33=-fl*zr
        s13=-fl1*zr*ro
        s23=ro*sfl3z
        if(iback.eq.1) then
          b11=t22+t33
          b22=t11+t33
          b33=t11+t22
          b55=t22-t33
          b66=t11-t33
          b99=t11-t22
          t4=f(4)+f(8)
          t5=t4+f(8)
          t4=t4+f(4)
          fp(1)=b11*f(1)-t21*f(2)-t31*f(3)-4.d0*f(5)+c22*f(7)
          fp(2)=-t12*f(1)+b22*f(2)-4.d0*f(6)+c11*f(11)
          fp(3)=b33*f(3)-c22*f(9)+c11*f(12)
          fp(4)=-s23*f(1)+s13*f(2)+t31*f(6)
          fp(5)=s13*f(3)+b55*f(5)-t21*f(6)-c22*f(10)
          fp(6)=s23*f(3)-t12*f(5)+b66*f(6)-c11*f(13)
          fp(7)=s22*f(1)-s12*f(2)-b55*f(7)+t31*f(9)+4.d0*f(10)+t12*f(11)
          fp(8)=s23*f(1)-s12*f(3)-t21*f(9)+t12*f(12)
          fp(9)=s23*f(2)-s22*f(3)-t12*t5+b99*f(9)-c11*f(14)
          fp(10)=s23*(f(4)-f(8))-s22*f(5)+s12*f(6)+s13*f(9)-b11*f(10)
     1      +t12*f(13)
          fp(11)=-s12*f(1)+s11*f(2)-t4*t31+t21*f(7)-b66*f(11)+4.d0*f(13)
          fp(12)=-s13*f(1)+s11*f(3)+t21*t5-t31*f(5)-b99*f(12)+c22*f(14)
          fp(13)=-t4*s13+s12*f(5)-s11*f(6)+t21*f(10)-s23*f(12)-b22*f(13)
          fp(14)=s12*t5-s13*f(7)-s11*f(9)+t31*f(10)-s23*f(11)+s22*f(12)
     1      -b33*f(14)
       else
          b11=t11+t22-t33
          b33=t11-t22-t33
          b44=t22-t11-t33
          b55=-t11-t22-t33
          b32=-t12-t12
          b42=t21+t21
          b52=-s12-s12
          b313=-s23-s23
          b414=s13+s13
          b914=t31+t31
          fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
          fp(2)=s12*f(1)-t33*f(2)-t21*f(3)+t12*f(4)-s13*f(13)-s23*f(14)
          fp(3)=s22*f(1)+b32*f(2)+b33*f(3)+c11*f(5)+b313*f(13)
          fp(4)=-s11*f(1)+b42*f(2)+b44*f(4)-c22*f(5)+b414*f(14)
          fp(5)=b52*f(2)+s11*f(3)-s22*f(4)+b55*f(5)-b313*f(11)
     +        +b414*f(12)
          fp(6)=4.d0*f(1)-b55*f(6)+c22*f(8)-c11*f(9)
          fp(7)=4.d0*f(2)+s12*f(6)+t33*f(7)-t21*f(8)+t12*f(9)-t31*f(13)
          fp(8)=4.d0*f(3)+s22*f(6)+b32*f(7)-b44*f(8)+c11*f(10)
          fp(9)=4.d0*f(4)-s11*f(6)+b42*f(7)-b33*f(9)-c22*f(10)
     +         +b914*f(14)
          fp(10)=4.d0*f(5)+b52*f(7)+s11*f(8)-s22*f(9)-b11*f(10)
     +          +b914*f(12)
          fp(11)=-t31*f(2)+s13*f(7)+s23*f(9)-t11*f(11)+t21*f(12)
     +      -s11*f(13)+s12*f(14)
          fp(12)=t31*f(3)+s23*f(7)-s13*f(8)+t12*f(11)-t22*f(12)
     +      +s12*f(13)-s22*f(14)
          fp(13)=s23*f(6)-c11*f(11)+t11*f(13)-t12*f(14)
          fp(14)=-t31*f(1)+s13*f(6)-c22*f(12)-t21*f(13)+t22*f(14)
        end if
      end if
      return
      end

      subroutine dermf(iq,z,f,fp,qff)
c*** calculates minor vector derivative (fp) in a fluid ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension f(*),fp(*)
      t=z-r(iq)
      if(t.eq.0.d0) then
        ro=rho(iq)
        flu=fcon(iq)*qff
        gr=g(iq)
      else
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        flu=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      end if
      t21=-4.d0*ro
      zr=1.d0/z
      t12=fl3*zr*zr/wsq
      t11=gr*t12-zr
      s11=ro*(gr*gr*t12-wsq)+t21*gr*zr
      c11=-t12/ro+1.d0/flu
      if(kg.eq.0) then
        fp(1)=t11*f(1)+c11*f(2)
        fp(2)=(s11-t21*ro)*f(1)-t11*f(2)
      else
        t22=-fl*zr
        s22=ro*t12
        b11=t11+t22
        s12=ro*b11
        if(iback.eq.1) then
          fp(1)=t22*f(1)-t21*f(2)-4.d0*f(3)
          fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
          fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
          fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
        else
          b33=t11-t22
          fp(1)=b11*f(1)+4.d0*f(3)-c11*f(4)
          fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
          fp(3)=s22*f(1)-(t12+t12)*f(2)+b33*f(3)+c11*f(5)
          fp(4)=-s11*f(1)+(t21+t21)*f(2)-b33*f(4)-4.d0*f(5)
          fp(5)=-(s12+s12)*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
        end if
      end if
      return
      end

      subroutine eifout(lsmin,fick,imtch)
c*** massages spheroidal mode eigenfunctions before output ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 ll,lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),inorm(mk),jjj(mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension zi(4)
      jfun=(kg+2)*2
c*** take out normalization here
      imax=0
      lsmin=max(lsmin,1)
      do iq=lsmin,n
        imax=max(inorm(iq),imax)
      enddo
      izer=0
      lsnew=lsmin
      do iq=lsmin,n
        iexp=inorm(iq)-imax
c        al=0.d0
c        if(iexp.ge.-80) al=2.d0**iexp
        al=2.d0**iexp
        if(al.eq.0.d0.and.izer.eq.0) lsnew=iq
        if(al.ne.0.d0) izer=1
        do j=1,jfun
          a(j,iq)=a(j,iq)*al
        enddo
c        print 222,iq,inorm(iq),(a(j,iq),j=1,jfun)
  222 format(i5,i10,6g14.6)
      enddo
      lsmin=lsnew
c      if(inorm(nsl)-imax.lt.-80) then
cc**** no point in outputing a mode with zero surface displacement
c          cg=0.
c          wray=0.
c          qinv=0.
c          lsmin=n+1
c          return
c      end if
c*** do integrals here
      call quod(lsmin,zi,fick,imtch)
      if(zi(1).eq.0.d0) then
         lsmin=n+1
         return
      end if
      fick=fick/zi(1)
      cg=(fl+0.5d0)*zi(2)/(w*sfl3*zi(1))
      wray=sqrt(zi(4)/zi(1))
      qinv=zi(3)/(wsq*zi(1))
      rnorm=1.d0/(w*sqrt(zi(1)))
c******
      lsm1=max(1,lsmin-1)
      do i=1,lsm1
        do j=1,6
          a(j,i)=0.d0
        enddo
      enddo
      i1=min(nic,max(2,lsmin))
      i2=nic
    5 if(i1.ne.i2) then
        do iq=i1,i2
          ff=fcon(iq)*(1.d0+xlam(iq)*fct)
          ll=lcon(iq)*(1.d0+qshear(iq)*fct)
          zr=1.d0/r(iq)
          sfl3z=sfl3*zr
          d=1.d0/(ccon(iq)*(1.d0+xa2(iq)*fct))
          v=a(2,iq)
          if(kg.eq.0) then
            a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(3,iq)
            a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(4,iq)/ll
            a(5,iq)=0.
            a(6,iq)=0.
          else
            a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(4,iq)
            a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(5,iq)/ll
            a(5,iq)=a(3,iq)
            a(6,iq)=4.d0*(a(6,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
          end if
          a(3,iq)=v
        enddo
      end if
      if(i2.eq.nsl) goto 25
      i1=min(nsl,max(lsmin,nocp1))
      i2=nsl
      goto 5
   25 i1=min(noc,max(lsmin,nicp1))
      i2=noc
   30 if(i1.ne.i2) then
        do iq=i1,i2
          zr=1.d0/r(iq)
          sfl3z=sfl3*zr
          ffi=1.d0/(flam(iq)*(1.d0+xlam(iq)*fct))
c*** p=y4=rR
          if(kg.eq.0) then
            p=a(2,iq)
            a(5,iq)=0.
            a(6,iq)=0.
          else
            p=a(3,iq)
            a(5,iq)=a(2,iq)
            a(6,iq)=4.d0*(a(4,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
          end if
          a(3,iq)=sfl3z*(g(iq)*a(1,iq)-p/rho(iq)+a(5,iq))/wsq
          a(2,iq)=sfl3z*a(3,iq)-a(1,iq)*zr+p*ffi
        a(4,iq)=sfl3z*(a(1,iq)+p*(qro(1,iq)/(rho(iq)**2)+g(iq)*ffi)/wsq)
        enddo
      end if
      if(n.eq.nsl.or.i2.eq.n) goto 55
      i1=nslp1
      i2=n
      goto 30
   55 continue
      do iq=lsmin,n
        zr=1.d0/r(iq)
        do j=1,jfun,2
          a(j,iq)=a(j,iq)*zr
          a(j+1,iq)=(a(j+1,iq)-a(j,iq))*zr
        enddo
      enddo
      do iq=lsmin,n
        do j=1,jfun
          a(j,iq)=a(j,iq)*rnorm
        enddo
      enddo
      if(lsmin.gt.2.or.l.gt.2) return
      if(l.eq.1) then
c*** for l=1, U and V tend to constants and the derivative is zero at the center -- this is quadratic approx
        a(1,1)=a(1,2)-.5d0*a(2,2)*r(2)
        a(2,1)=0.d0
        a(3,1)=a(3,2)-.5d0*a(4,2)*r(2)
        a(4,1)=0.d0
        if(kg.ne.0) a(6,1)=1.5d0*a(5,2)/r(2)-.5d0*a(6,2)
      else
c*** for l=2, U and V go to zero but derivs tend to constants
        a(2,1)=1.5d0*a(1,2)/r(2)-.5d0*a(2,2)
        a(4,1)=1.5d0*a(3,2)/r(2)-.5d0*a(4,2)
      end if
      return
      end

      subroutine sdepth(wdim,ls)
c*** finds starting level,ls, for a given l and w ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      data aw,bw,dw/-2.d-3,2.25d-3,1.28d-3/
      q=0.d0
      w=wdim/wn
      wsoc=aw+dw*fl
      if(wdim.gt.wsoc) goto 10
      call startl(nocp1,nsl,fmu,ls,q)
      if(ls.eq.nsl) ls=ls-1
      if(ls.gt.nocp1) return
   10 wsic=aw+bw*fl
      if(wdim.gt.wsic) goto 20
      call startl(nicp1,noc,flam,ls,q)
      if(ls.eq.noc) ls=ls-1
      if(ls.gt.nicp1) return
   20 call startl(2,nic,fmu,ls,q)
      if(ls.eq.nic) ls=ls-1
      return
      end

      subroutine sfbm(ass,kg,iback)
c*** convert minor vector at a solid/fluid boundary ***
      implicit real*8(a-h,o-z)
      dimension ass(14),as(14)
      do j=1,14
        as(j)=ass(j)
        ass(j)=0.d0
      enddo
      if(iback.ne.1) then
        if(kg.eq.0) then
          ass(1)=as(3)
          ass(2)=as(5)
        else
          ass(1)=as(8)
          ass(2)=-as(12)
          ass(3)=as(3)
          ass(4)=-as(10)
          ass(5)=as(5)
        end if
      else
        if(kg.eq.0) then
          ass(1)=-as(3)
        else
          ass(1)=as(7)
          ass(2)=-as(9)
          ass(3)=-as(10)
          ass(4)=-as(14)
        end if
      end if
      return
      end

      subroutine fsbm(ass,kg,iback)
c*** convert minor vector at a fluid/solid boundary ***
      implicit real*8(a-h,o-z)
      dimension ass(14),as(14)
      do j=1,14
        as(j)=ass(j)
        ass(j)=0.d0
      end do
      if(iback.ne.1) then
        if(kg.eq.0) then
          ass(1)=as(1)
          ass(4)=-as(2)
        else
          ass(6)=as(1)
          ass(14)=as(2)
          ass(1)=as(3)
          ass(9)=as(4)
          ass(4)=-as(5)
        end if
      else
        if(kg.eq.0) then
          ass(1)=-as(1)
        else
          ass(1)=-as(1)
          ass(3)=-as(2)
          ass(5)=-as(3)
          ass(12)=as(4)
        end if
      end if
      return
      end

      subroutine zknt(s,sp,f,x,y,ifsol,isplit)
c*** given minor vector and derivs,constructs mode count ***
      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension s(*),sp(*),f(*)
      isplit=0
      if(ifsol.eq.0.and.kg.eq.0) then
        y1=s(2)
        y2=f(2)
        y2p=s(2)+sp(2)*(y-x)
        t1=s(1)
        t2=f(1)
      else
        y1=s(5)
        y2=f(5)
        y2p=s(5)+sp(5)*(y-x)
        t1=s(3)-s(4)
        t2=f(3)-f(4)
      end if
      if(y1*y2.le.0.d0) goto 5
      if(y2*y2p.lt.0.d0) isplit=1
      return
    5 if(t1*t2.le.0.d0) then
        isplit=1
        return
      end if
      if(y2.eq.0.d0) then
        tes=-y1*t1
      else
        tes=y2*t1-y1*t2
      endif
      if(tes.lt.0.d0) kount=kount+1
      if(tes.gt.0.d0) kount=kount-1
      return
      end

      subroutine baylis(q)
c    baylis returns the coefficients for rks integration.
c    see e. baylis shanks(1966 a. m. s.) and references therein for the
c    coefficients. the eight runge-kutta-shanks formulae are (1-1) (2-2)
c    (3-3) (4-4) (5-5) (6-6) (7-7) (8-10). for orders greater than 4 the
c    formulae are approximate rather than exact so incurring less roundoff.
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),maxo,i
      ds=q*dabs(dx)
      do 10 j=1,maxo-1
      if(ds.gt.step(j)) go to 10
      i=j
      go to 15
   10 continue
      i=maxo
   15 c(1)=0.d0
      go to (1,2,3,4,5,6,7,8),i
    1 b(1)=dx
      return
    2 c(2)=dx
      b(1)=dx
      b(2)=.5d0*dx
      b(3)=1.d0
      return
    3 c(2)=.5d0*dx
      c(3)=dx
      b(1)=c(2)
      b(2)=-dx
      b(3)=-2.d0
      b(4)=.16666666666667d0*dx
      b(5)=4.d0
      b(6)=1.d0
      return
    4 c(2)=.01d0*dx
      c(3)=.6d0*dx
      c(4)=dx
      b(1)=c(2)
      b( 2)=-.17461224489790d+02*dx
      b( 3)=-.10343618513324d+01
      b( 4)= .59691275167780d+02*dx
      b( 5)=-.10140620414448d+01
      b( 6)= .30814908546230d-01
      b( 7)=-.25555555555556d+01*dx
      b( 8)=-.11165449632656d+01
      b( 9)=-.22568165070006d+00
      b(10)=-.49077733860351d-01
      return
    5 c( 2)= 1.1111111111111d-04*dx
      c( 3)= 3.0d-01*dx
      c( 4)= 7.5d-01*dx
      c( 5)= dx
      b( 1)=c(2)
      b( 2)=-.40470000000000d+03*dx
      b( 3)=-.10007412898443d+01
      b( 4)= .25301250000000d+04*dx
      b( 5)=-.10004446420631d+01
      b( 6)= .74107010523195d-03
      b( 7)=-.11494333333333d+05*dx
      b( 8)=-.10004929965491d+01
      b( 9)= .52629261224803d-03
      b(10)=-.12029545422812d-03
      b(11)= .92592592592593d-01*dx
      b(12)= .00000000000000d+00
      b(13)= .47619047619048d+01
      b(14)= .42666666666667d+01
      b(15)= .77142857142857d+00
      return
    6 c(2)=3.3333333333333d-03*dx
      c(3)=.2d0*dx
      c(4)=.6d0*dx
      c(5)=9.3333333333333d-01*dx
      c(6)=dx
      b( 1)=c(2)
      b( 2)=-.58000000000000d+01*dx
      b( 3)=-.10344827586207d+01
      b( 4)= .64600000000000d+02*dx
      b( 5)=-.10216718266254d+01
      b( 6)= .30959752321982d-01
      b( 7)=-.62975802469136d+03*dx
      b( 8)=-.10226149961576d+01
      b( 9)= .24906685695466d-01
      b(10)=-.37737402568887d-02
      b(11)=-.54275714285714d+04*dx
      b(12)=-.10225567867765d+01
      b(13)= .25375487829097d-01
      b(14)=-.31321559234596d-02
      b(15)= .12921040478749d-03
      b(16)= .53571428571429d-01*dx
      b(17)= .00000000000000d+00
      b(18)= .61868686868687d+01
      b(19)= .77777777777778d+01
      b(20)= .40909090909091d+01
      b(21)=-.38888888888889d+00
      return
    7 c(2)=5.2083333333333d-03*dx
      c(3)=1.6666666666667d-01*dx
      c(4)=.5d0*dx
      c(5)=dx
      c(6)=8.3333333333333d-01*dx
      c(7)=dx
      b( 1)=c(2)
      b( 2)=-.25000000000000d+01*dx
      b( 3)=-.10666666666667d+01
      b( 4)= .26166666666667d+02*dx
      b( 5)=-.10421204027121d+01
      b( 6)= .61228682966918d-01
      b( 7)=-.64500000000000d+03*dx
      b( 8)=-.10450612653163d+01
      b( 9)= .51262815703925d-01
      b(10)=-.77519379844961d-02
      b(11)=-.93549382716049d+02*dx
      b(12)=-.10450293206756d+01
      b(13)= .48394546673620d-01
      b(14)=-.11877268228307d-01
      b(15)=-.39590894094358d-03
      b(16)= .35111904761905d+03*dx
      b(17)=-.10446476812124d+01
      b(18)= .52479782656724d-01
      b(19)=-.71200922221468d-02
      b(20)=-.61029361904114d-03
      b(21)= .27463212856852d-02
      b(22)= .46666666666667d-01*dx
      b(23)= .57857142857143d+01
      b(24)= .78571428571429d+01
      b(25)= .00000000000000d+00
      b(26)= b(23)
      b(27)= .10000000000000d+01
      return
    8 c(2)=.14814814814815d0*dx
      c(3)=.22222222222222d0*dx
      c(4)=.33333333333333d0*dx
      c(5)= .5d0*dx
      c(6)=.66666666666667d0*dx
      c(7)=.16666666666667d0*dx
      c(8)=dx
      c(9)=.83333333333333d0*dx
      c(10)=dx
      b( 1)=c(2)
      b( 2)= .55555555555556d-01*dx
      b( 3)= .30000000000000d+01
      b( 4)= .83333333333333d-01*dx
      b( 5)= .00000000000000d+00
      b( 6)= .30000000000000d+01
      b( 7)= .12500000000000d+00*dx
      b( 8)= .00000000000000d+00
      b( 9)= .00000000000000d+00
      b(10)= .30000000000000d+01
      b(11)= .24074074074074d+00*dx
      b(12)= .00000000000000d+00
      b(13)=-.20769230769231d+01
      b(14)= .32307692307692d+01
      b(15)= .61538461538461d+00
      b(16)= .90046296296295d-01*dx
      b(17)= .00000000000000d+00
      b(18)=-.13881748071980d+00
      b(19)= .24832904884319d+01
      b(20)=-.21182519280206d+01
      b(21)= .62467866323908d+00
      b(22)=-.11550000000000d+02*dx
      b(23)=-.35064935064935d+00
      b(24)= .50389610389610d+01
      b(25)=-.28398268398268d+01
      b(26)= .52813852813853d+00
      b(27)=-.34632034632035d+01
      b(28)=-.44097222222222d+00*dx
      b(29)=-.14173228346457d+00
      b(30)= .53385826771654d+01
      b(31)=-.35905511811023d+01
      b(32)= .70866141732284d-01
      b(33)=-.45354330708661d+01
      b(34)=-.31496062992126d-01
      b(35)= .18060975609756d+01*dx
      b(36)=-.54692775151925d-01
      b(37)= .47967589466576d+01
      b(38)=-.22795408507765d+01
      b(39)= .48615800135044d-01
      b(40)=-.34031060094530d+01
      b(41)=-.40513166779204d-01
      b(42)= .48615800135044d+00
      b(43)= .48809523809524d-01*dx
      b(44)= .65853658536585d+00
      b(45)= .66341463414634d+01
      b(46)= .52682926829268d+01
      i=10
      return
      end

      subroutine grav(g,rho,qro,r,n)
c*** given rho and spline coeffs,computes gravity ***
      implicit real*8(a-h,o-z)
      dimension g(*),rho(*),qro(3,1),r(*)
      g(1)=0.d0
      do i=2,n
        im1=i-1
        del=r(i)-r(im1)
        rn2=r(im1)*r(im1)
        trn=2.d0*r(im1)
        c1=rho(im1)*rn2
        c2=(qro(1,im1)*rn2+trn*rho(im1))*0.5d0
        c3=(qro(2,im1)*rn2+trn*qro(1,im1)+rho(im1))/3.d0
        c4=(qro(3,im1)*rn2+trn*qro(2,im1)+qro(1,im1))*.25d0
        c5=(trn*qro(3,im1)+qro(2,im1))*0.2d0
        g(i)=(g(im1)*rn2+4.d0*del*(c1+del*(c2+del*(c3+del*(c4
     +   +del*(c5+del*qro(3,im1)/6.d0))))))/(r(i)*r(i))
      enddo
      return
      end

      subroutine startl(jf,jl,v,ls,q)
c*** finds start level between jf and jl using velocityv and ang. ord. l.
c*** upon entry q is the value of the exponent at r(jf) or at the turning
c*** point(q=0) depending on previous calls to startl. upon exit q is the
c*** value of the exponent at the starting level ls.
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      save rrlog, vertno
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension rrlog(mk),p(mk),v(*)
      data ifirst/1/
      if(ifirst.eq.1) then
        ifirst=0
        vertno=-dlog(eps)
        do i=3,n
          rrlog(i)=.5d0*dlog(r(i)/r(i-1))
        enddo
      end if
      do j=jf,jl
        pp=fl3-wsq*r(j)*r(j)*rho(j)/v(j)
        if(pp.le.0.d0) goto 15
        p(j)=dsqrt(pp)
      enddo
   15 p(j)=0.d0
   20 k=j
      j=j-1
      if(j.le.jf) go to 25
      q=q+rrlog(k)*(p(j)+p(k))
      if(q.lt.vertno) go to 20
      ls=j
      return
   25 ls=jf
      return
      end

c      subroutine startl(jf,jl,v,ls,q)
cc*** finds start level between jf and jl using velocity v and ang. ord. l.
cc*** upon entry q is the value of the exponent at r(jf) or at the turning
cc*** point(q=0) depending on previous calls to startl. upon exit q is the
cc*** value of the exponent at the starting level ls.
c      implicit real*8(a-h,o-z)
c      save
c      parameter (mk=950)
c      real*8 lcon,ncon,lspl,nspl
c      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
c     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
c     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
c     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
c      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
c     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
c      dimension p(mk),v(*)
c      vertno=-dlog(eps)
c      print *,jf,jl,q
c      do j=jf,jl
c        pp=fl3-wsq*r(j)*r(j)*rho(j)/v(j)
c        if(pp.le.0.d0) goto 15
c        p(j)=dsqrt(pp)
c      enddo
c   15 p(j)=0.d0
c   20 k=j
c      j=j-1
c      if(j.le.jf) go to 25
c      q=q+(p(j)+p(k))*0.5d0*dlog(r(k)/r(j))
c      if(q.lt.vertno) go to 20
c      ls=j
c      print *,ls
c      return
c   25 ls=jf
c      print *,ls
c      return
c      end

      subroutine steps(eps)
c*** computes 8 dimensionless step sizes for rks integration
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      ps=dlog(eps)
      fac=1.d0
      do n=1,8
        fn=n+1
        fac=fac*fn
        x=(dlog(fac)+ps)/fn
        x=exp(x)
        s=x
        do i=1,n
          s=x*dexp(-s/fn)
        enddo
        step(n)=s
      enddo
      return
      end

      subroutine drspln(i1,i2,x,y,q,f)
      implicit real*8(a-h,o-z)
c   rspln computes cubic spline interpolation coefficients
c   for y(x) between grid points i1 and i2 saving them in q.  the
c   interpolation is continuous with continuous first and second
c   derivitives.  it agrees exactly with y at grid points and with the
c   three point first derivitives at both end points (i1 and i2).
c   x must be monotonic but if two successive values of x are equal
c   a discontinuity is assumed and seperate interpolation is done on
c   each strictly monotonic segment.  the arrays must be dimensioned at
c   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
c   for rspln.
c                                                     -rpb
      dimension x(*),y(*),q(3,*),f(3,*)
      j1=i1
   10 if(j1.ge.i2) return
c   search for discontinuities.
      do i=j1+1,i2
        if(x(i).eq.x(i-1)) then
           j2=i-1
           goto 15
        end if
      enddo
      j2=i2
   15 if(j2.eq.j1+1) then
c  only two points.  use linear interpolation.
        q(1,j1)=(y(j2)-y(j1))/(x(j2)-x(j1))
        q(1,j2)=q(1,j1)
        q(2,j1)=0.d0
        q(2,j2)=0.d0
        q(3,j1)=0.d0
        q(3,j2)=0.d0
      else
c  more than two points.  do spline interpolation.
        a0=0.d0
        h=x(j1+1)-x(j1)
        h2=x(j1+2)-x(j1)
        y0=h*h2*(h2-h)
        h=h*h
        h2=h2*h2
c  calculate derivitive at near end.
        b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
        b1=b0
c  explicitly reduce banded matrix to an upper banded matrix.
        j22=j2-2
        do i=j1,j22
          h=x(i+1)-x(i)
          y0=y(i+1)-y(i)
          h2=h*h
          ha=h-a0
          h2a=h-2.d0*a0
          h3a=2.d0*h-3.d0*a0
          h2b=h2*b0
          q(1,i)=h2/ha
          q(2,i)=-ha/(h2a*h2)
          q(3,i)=-h*h2a/h3a
          f(1,i)=(y0-h*b0)/(h*ha)
          f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
          f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
          a0=q(3,i)
          b0=f(3,i)
        enddo
c  take care of last two rows.
        i=j2-1
        h=x(i+1)-x(i)
        y0=y(i+1)-y(i)
        h2=h*h
        ha=h-a0
        h2a=h*ha
        h2b=h2*b0-y0*(2.d0*h-a0)
        q(1,i)=h2/ha
        f(1,i)=(y0-h*b0)/h2a
        ha=x(j22)-x(i+1)
        y0=-h*ha*(ha+h)
        ha=ha*ha
c  calculate derivitive at far end.
        y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j22)*h2)/y0
        q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
        q(2,i)=f(1,i)-q(1,i)*q(3,i)
c  solve upper banded matrix by reverse iteration.
        do j=j1,j22
          k=i-1
          q(1,i)=f(3,k)-q(3,k)*q(2,i)
          q(3,k)=f(2,k)-q(2,k)*q(1,i)
          q(2,k)=f(1,k)-q(1,k)*q(3,k)
          i=k
        enddo
        q(1,i)=b1
c  fill in the last point with a linear extrapolation.
        q(1,j2)=y0
        q(2,j2)=0.d0
        q(3,j2)=0.d0
      end if
c  see if this discontinuity is the last.
      j1=j2+1
      goto 10
      end

      subroutine rprop(jf,jl,f)
c*** propagates soln ,f, for radial modes from jf to jl ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl,nn
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),dum(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension h(2,10),s(2),f(2)
      y=r(jf)
      vy=sqrt((flam(jf)+2.d0*fmu(jf))/rho(jf))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.eq.jl) return
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      vx=vy
      vy=sqrt((flam(i)+2.d0*fmu(i))/rho(i))
      q=max(w/vx+1.d0/x,w/vy+1.d0/y)
      del=step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q)
        dxs=dx
        s(1)=f(1)
        s(2)=f(2)
        do ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
          ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
          if(ifanis.eq.0) then
           nn=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
           cc=ff+nn+nn
           aa=cc
          else
           nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
           cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
           aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
          end if
          z=1.d0/z
          a21=-ro*wsq+4.d0*z*(z*(aa-nn-ff*ff/cc)-ro*gr)
          h(1,ni)=(f(2)-2.d0*ff*z*f(1))/cc
          h(2,ni)=a21*f(1)+2.d0*z*f(2)*(ff/cc-1.d0)
          call rkdot(f,s,h,2,ni)
        enddo
        if(knsw.eq.1) then
          if(s(2)*f(2).le.0.d0) then
            if(f(2).eq.0.d0) then
              tes=-s(2)*s(1)
            else
              tes=f(2)*s(1)-s(2)*f(1)
            endif
            if(tes.lt.0.d0) kount=kount+1
            if(tes.gt.0.d0) kount=kount-1
          end if
        end if
        x=y
        if(y.ne.r(i)) go to 15
      go to 10
      end

      subroutine erprop(jf,jl,f)
c*** propagates soln ,f, for radial modes from jf to jl ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl,nn
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),dum(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension h(5,10),s(5),f(5)
      fot=4.d0/3.d0
      y=r(jf)
      vy=sqrt((flam(jf)+2.d0*fmu(jf))/rho(jf))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.eq.jl) return
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      vx=vy
      vy=sqrt((flam(i)+2.d0*fmu(i))/rho(i))
      q=max(w/vx+1.d0/x,w/vy+1.d0/y)
      del=step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q)
        dxs=dx
        do nf=1,5
          s(nf)=f(nf)
        enddo
        do ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
          ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
          if(ifanis.eq.0) then
           nn=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
           cc=ff+nn+nn
           aa=cc
          else
           nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
           cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
           aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
          end if
          h(3,ni)=ro*f(1)*f(1)*z*z
          z=1.d0/z
          a21=-ro*wsq+4.d0*z*(z*(aa-nn-ff*ff/cc)-ro*gr)
          h(1,ni)=(f(2)-2.d0*ff*z*f(1))/cc
          h(2,ni)=a21*f(1)+2.d0*z*f(2)*(ff/cc-1.d0)
          qp=h(1,ni)/z
          qrka=(4.d0*(aa+ff-nn)+cc)*qkappa(iq)/9.d0
          qrmu=(aa+cc-2.d0*ff+5.d0*nn+6.d0*ll)*qshear(iq)/15.d0
          t1=(qp+2.d0*f(1))**2
          t2=fot*(qp-f(1))**2
          h(4,ni)=t1*qrka+t2*qrmu
          h(5,ni)=qp*(cc*qp+4.d0*ff*f(1))+4.d0*f(1)*f(1)*(aa-nn-ro*gr/z)
          call rkdot(f,s,h,5,ni)
        enddo
        x=y
        if(y.ne.r(i)) go to 15
      go to 10
      end

      subroutine radint(ls,f)
c*** this does starting values for integrals close to center of
c*** earth where U is proportional to r.
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),dum(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension f(*)
      f(1)=a(1,ls)
      f(2)=a(2,ls)
      fnn=r(ls)/5.d0
      fmm=r(ls)/3.d0
      qff=1.d0+xlam(ls)*fct
      qll=1.d0+qshear(ls)*fct
      qaa=1.d0+xa2(ls)*fct
      ff=fcon(ls)*qff
      ll=lcon(ls)*qll
      if(ifanis.eq.0) then
        nn=ll
        cc=ff+nn+nn
        aa=cc
      else
        nn=ncon(ls)*qll
        cc=ccon(ls)*qaa
        aa=acon(ls)*qaa
      end if
      xx=4.d0*(aa+ff-nn)+cc
      f1sq=f(1)*f(1)
      f(3)=rho(ls)*f1sq*r(ls)*r(ls)*fnn
      f(4)=f1sq*xx*qkappa(ls)*fmm
      f(5)=(xx*fmm-4.d0*rho(ls)*g(ls)*r(ls)*fnn)*f1sq
      return
      end


      subroutine tprop(jf,jl,f)
c*** propagates f from jf to jl - toroidal modes ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),dum(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension h(2,10),s(2),f(2)
      fl3m2=fl3-2.d0
      y=r(jf)
      qy=1.d0/y+sqrt(abs(rho(jf)*wsq/fmu(jf)-fl3/(y*y)))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.eq.jl) return
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qll=1.d0+qshear(iq)*fct
      qx=qy
      qy=1.d0/y+sqrt(abs(rho(i)*wsq/fmu(i)-fl3/(y*y)))
      q=max(qx,qy)
      del=step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q)
        dxs=dx
        s(1)=f(1)
        s(2)=f(2)
        do ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
          nn=ll
          if(ifanis.ne.0) nn=(ncon(iq)+
     +        t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          z=1.d0/z
          h(1,ni)=z*f(1)+f(2)/ll
          h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
          call rkdot(f,s,h,2,ni)
        enddo
        if(knsw.eq.1) then
          if(s(2)*f(2).le.0.d0) then
            if(f(2).eq.0.d0) then
              tes=-s(2)*s(1)
            else
              tes=f(2)*s(1)-s(2)*f(1)
            endif
            if(tes.lt.0.d0) kount=kount+1
            if(tes.gt.0.d0) kount=kount-1
          end if
        end if
        x=y
        if(y.ne.r(i)) goto 15
      go to 10
      end

      subroutine etprop(jf,jl,f)
c*** propagates f from jf to jl - toroidal modes ***
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,mk),dum(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension h(6,10),s(6),f(6)
      fl3m2=fl3-2.d0
      y=r(jf)
      qy=1.d0/y+sqrt(abs(rho(jf)*wsq/fmu(jf)-fl3/(y*y)))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.eq.jl) return
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qll=1.d0+qshear(iq)*fct
      qx=qy
      qy=1.d0/y+sqrt(abs(rho(i)*wsq/fmu(i)-fl3/(y*y)))
      q=max(qx,qy)
      del=step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q)
        dxs=dx
        do nf=1,6
          s(nf)=f(nf)
        enddo
        do ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
          nn=ll
          if(ifanis.ne.0) nn=(ncon(iq)+
     +        t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          h(3,ni)=ro*f(1)*f(1)*z*z
          z=1.d0/z
          h(1,ni)=z*f(1)+f(2)/ll
          h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
          t1=(h(1,ni)/z-f(1))**2
          t2=fl3m2*f(1)*f(1)
          h(4,ni)=nn*f(1)*f(1)
          h(5,ni)=(t1+t2)*ll*qshear(iq)
          h(6,ni)=ll*t1+t2*nn
          call rkdot(f,s,h,6,ni)
        enddo
        x=y
        if(y.ne.r(i)) goto 15
      go to 10
      end

      subroutine rkdot(f,s,h,nvec,ni)
c*** performs dot product with rks coefficients ***
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension s(*),f(*),h(nvec,*)
      goto (1,2,3,4,5,6,7,8,9,10),ni
    1 do 21 j=1,nvec
   21 f(j)=s(j)+b(1)*h(j,1)
      return
    2 do 22 j=1,nvec
   22 f(j)=s(j)+b(2)*(h(j,1)+b(3)*h(j,2))
      return
    3 do 23 j=1,nvec
   23 f(j)=s(j)+b(4)*(h(j,1)+b(5)*h(j,2)+b(6)*h(j,3))
      return
    4 do 24 j=1,nvec
   24 f(j)=s(j)+b(7)*(h(j,1)+b(8)*h(j,2)+b(9)*h(j,3)+b(10)*h(j,4))
      return
    5 do 25 j=1,nvec
   25 f(j)=s(j)+b(11)*(h(j,1)+b(12)*h(j,2)+b(13)*h(j,3)+b(14)*h(j,4)+
     +b(15)*h(j,5))
      return
    6 do 26 j=1,nvec
   26 f(j)=s(j)+b(16)*(h(j,1)+b(17)*h(j,2)+b(18)*h(j,3)+b(19)*h(j,4)+
     +b(20)*h(j,5)+b(21)*h(j,6))
      return
    7 do 27 j=1,nvec
   27 f(j)=s(j)+b(22)*(h(j,1)+b(23)*h(j,3)+b(24)*h(j,4)+b(25)*h(j,5)+
     +b(26)*h(j,6)+b(27)*h(j,7))
      return
    8 do 28 j=1,nvec
   28 f(j)=s(j)+b(28)*(h(j,1)+b(29)*h(j,3)+b(30)*h(j,4)+b(31)*h(j,5)+
     +b(32)*h(j,6)+b(33)*h(j,7)+b(34)*h(j,8))
      return
    9 do 29 j=1,nvec
   29 f(j)=s(j)+b(35)*(h(j,1)+b(36)*h(j,3)+b(37)*h(j,4)+b(38)*h(j,5)+
     +b(39)*h(j,6)+b(40)*h(j,7)+b(41)*h(j,8)+b(42)*h(j,9))
      return
   10 do 30 j=1,nvec
   30 f(j)=s(j)+b(43)*(h(j,1)+h(j,10)+b(44)*(h(j,4)+h(j,6))+
     +b(45)*h(j,5)+b(46)*(h(j,7)+h(j,9)))
      return
      end

      subroutine fpsm(ls,nvefm,ass)
c*** spheroidal mode start solution in a fluid region using sph. bessel fns.
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension ass(*)
      x=r(ls)
      fla=flam(ls)*(1.d0+xlam(ls)*fct)
      vpsq=fla/rho(ls)
      xi=g(ls)/x
      qsq=(wsq+float(kg)*4.d0*rho(ls)+xi-fl3*xi*xi/wsq)/vpsq
      zsq=qsq*x*x
      call bfs(l,zsq,eps,fp)
      if(kg.ne.0) then
        u=(fl-fp)/qsq
        c1=fl*g(ls)-wsq*x
        c2=fl2*c1*0.25d0/x-rho(ls)*fl
        ass(1)=-x*fl*vpsq-c1*u
        ass(2)=-x*fl*fla
        ass(3)=-fl*fl2*vpsq*0.25d0-u*c2
        ass(4)=x*fla*c1
        ass(5)=-x*fla*c2
      else
        ass(1)=-(fl3*xi/wsq+fp)/qsq
        ass(2)=x*fla
      end if
      sum=ass(1)*ass(1)
      do i=2,nvefm
        sum=sum+ass(i)*ass(i)
      enddo
      sum=1.d0/dsqrt(sum)
      if(ass(nvefm).lt.0.d0) sum=-sum
      do i=1,nvefm
        ass(i)=ass(i)*sum
      enddo
      return
      end

      subroutine spsm(ls,nvesm,ass)
c*** spheroidal mode start solution in a solid region using sph. bessel fns.
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(6,2),e(15),ass(*)
      x=r(ls)
      ro=rho(ls)
      fu=fmu(ls)*(1.d0+qshear(ls)*fct)
      flu=flam(ls)*(1.d0+xlam(ls)*fct)+2.d0*fu
      vssq=fu/ro
      vpsq=flu/ro
      zeta=4.d0*ro
      xi=g(ls)/x
      alfsq=(wsq+float(kg)*zeta+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vpsq*vssq))
      fksq=.5d0*(alfsq+betasq+delsq)
      qsq=fksq-delsq
      do k=1,2
        if(k.eq.1) then
          zsq=qsq*x*x
          b=xi/(vssq*(betasq-qsq))
        else
          zsq=fksq*x*x
          b=-flu/(fu*fl3*b)
        end if
        call bfs(l,zsq,eps,fp)
        a(1,k)=fl3*b+fp
        a(2,k)=1.d0+b+b*fp
        a(3,k)=-zsq
        a(4,k)=b*a(3,k)
        a(5,k)=1.d0
        a(6,k)=fl1-fl3*b
      enddo
      jj=3+2*kg
      kk=jj+1
      ll=0
      do i=1,jj
        i1=i+1
        do j=i1,kk
          ll=ll+1
          e(ll)=a(i,1)*a(j,2)-a(j,1)*a(i,2)
        enddo
      enddo
      if(kg.eq.0) then
        ass(1)=x*x*e(1)
        ass(2)=fu*x*sfl3*(2.d0*e(1)-e(5))
        ass(3)=fu*x*(e(3)-2.d0*e(1))
        ass(4)=x*(flu*e(4)+4.d0*fu*e(1))
        ass(5)=fu*(flu*(e(6)+2.d0*e(4))+4.d0*fu*(fl3*(e(5)-e(1))
     +     -e(3)+2.d0*e(1)))
      else
        c0=wsq-xi*fl
        c1=ro*fl+0.25d0*fl2*c0
        c2=2.d0*fu/x
        c3=c2*(fl-1.d0)
        ass(6)=x*x*(c0*e(1)-zeta*(fl*e(8)-e(4)))
        ass(14)=flu*(fl*e(6)-e(2))
        ass(13)=fu*sfl3*(fl*e(7)-e(3))
        ass(1)=x*(c1*e(1)-ro*(fl*e(9)-e(5)))
        ass(7)=x*flu*(c0*e(2)-zeta*fl*e(11))/sfl3+c2*sfl3*ass(6)
        ass(8)=x*fu*(c0*e(3)-zeta*fl*e(13))-c2*ass(6)
        ass(12)=(flu*fl*e(10)+2.d0*(ass(14)+sfl3*ass(13)))*fu/x
        ass(2)=flu*(c1*e(2)-ro*fl*e(12))/sfl3+c2*sfl3*ass(1)
        ass(3)=fu*(c1*e(3)-ro*fl*e(14))-c2*ass(1)
        ass(9)=(x*c0*ass(14)+sfl3*ass(7)-c3*fl*ass(6))/fl
        ass(11)=(sfl3*ass(12)+c3*(sfl3*ass(14)-fl*ass(13)))/fl
        ass(4)=(c1*ass(14)+sfl3*ass(2)-c3*fl*ass(1))/fl
        ass(10)=(x*c0*ass(11)-c3*(sfl3*ass(9)+fl*ass(7)))/sfl3
        ass(5)=(c1*ass(11)-c3*(sfl3*ass(4)+fl*ass(2)))/sfl3
      end if
      sum=ass(1)*ass(1)
      do i=2,nvesm
        sum=sum+ass(i)*ass(i)
      enddo
      sum=1.d0/dsqrt(sum)
      if(ass(5).lt.0.d0) sum=-sum
      do i=1,nvesm
        ass(i)=ass(i)*sum
      enddo
      return
      end

      subroutine rps(i,a)
c*** radial mode start soln using sph bessel fns.
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(2)
      fla=flam(i)*(1.d0+xlam(i)*fct)
      sig=fla+2.d0*fmu(i)*(1.d0+qshear(i)*fct)
      zsq=r(i)*r(i)*rho(i)*(wsq+4.d0*g(i)/r(i))/sig
      call bfs(1,zsq,eps,fp)
      a(1)=r(i)
      a(2)=sig*fp+2.d0*fla
      return
      end

      subroutine tps(i,a)
c*** toroidal mode start soln using sph bessel fns.
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(2)
      fu=fmu(i)*(1.d0+qshear(i)*fct)
      zsq=r(i)*r(i)*wsq*rho(i)/fu
      call bfs(l,zsq,eps,fp)
      a(1)=r(i)
      a(2)=fu*(fp-1.d0)
      return
      end

      subroutine bfs(l,xsq,eps,fp)
c  this routine calculates spherical bessel function of the ist kind.
c  fp is equivalent to (r*dj/dr)/j
c  where r is radius and j is the sbf of order l and argument x=k*r
c  the technique employs the continued fraction approach
c  described in w. lentz's article in applied optics, vol.15, #3, 1976
      implicit real*8(a-h,o-z)
      real*8 numer,nu
c*** positive argument uses continued fraction
      if(xsq.gt.0.d0) then
        x=sqrt(xsq)
        lp1=l+1
        rx=2.0d0/x
        nu=lp1-0.5d0
        rj=nu*rx
        rx=-rx
        denom=(nu+1.d0)*rx
        numer=denom+1.0d0/rj
        rj=rj*numer/denom
        nm1=1
    5     nm1=nm1+1
          rx=-rx
          a3=(nu+nm1)*rx
          denom=a3+1.d0/denom
          numer=a3+1.d0/numer
          ratio=numer/denom
          rj=rj*ratio
          if(abs(abs(ratio)-1.d0).gt.eps) goto 5
        fp=rj*x-lp1
      else
c  series solution
        f=1.d0
        fp=l
        a=1.d0
        b=l+l+1.d0
        c=2.d0
        d=l+2.d0
   15     a=-a*xsq/(c*(b+c))
          f=f+a
          fp=fp+a*d
          if(abs(a*d).lt.eps) goto 20
          c=c+2.d0
          d=d+2.d0
          goto 15
   20   fp=fp/f
      end if
      return
      end

      subroutine quod(ls,zi,fick,imtch)
c*** does integrals using RK -- don't have to worry about model knot spacing
      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension f(10),zi(*),s(10),h(140)
      jj=2*(kg+2)
      nint=jj+4
      fick=0.d0
      do i=1,nint
        f(i)=0.d0
      enddo
      if(ls.gt.nocp1) goto 25
      if(ls.gt.nicp1) goto 15
c*** propagate through inner core ***
c*** should try and get better estimate of start solutions -- this uses power series
c      if(ls.le.2) call solint(ls,f)
      call solint(ls,f)
      call esprop(ls,nic,f,h,nint)
      fick=f(jj+1)
      do i=1,4
        f(i+jj-2)=f(i+jj)
      enddo
      if(imtch.eq.nic) goto 100
   15 is=max(ls,nicp1) 
c*** propagate through outer core ***
      call efprop(is,noc,f,h,nint-2) 
      if(imtch.eq.noc) goto 100
      do i=1,4
        f(nint+1-i)=f(nint-1-i)
      enddo
   25 is=max(ls,nocp1)      
c*** propagate through mantle ***
      call esprop(is,nsl,f,h,nint)   
      if(nsl.ne.n) then
c*** propagate through ocean ***
        do i=1,4
          f(i+jj-2)=f(i+jj)
        enddo
        call efprop(nslp1,n,f,h,nint-2)
        do i=1,4
          zi(i)=f(i+jj-2)
        enddo
      else                
        do i=1,4
          zi(i)=f(i+jj)
        enddo
      end if
      return
  100 continue
c*** integrate to CMB or ICB from surface
      do i=1,nint-2
        s(i)=f(i)
      enddo
      do i=1,nint
        f(i)=0.d0
      enddo
c*** prop to ocean floor
      if(n.ne.nsl) then 
         call efprop(n,nslp1,f,h,nint-2)
         do i=1,4
            f(i+jj)=f(i+jj-2)
         enddo
      end if
c*** prop to cmb
      call esprop(nsl,nocp1,f,h,nint)   
      do i=1,4
        f(i+jj-2)=f(i+jj)
      enddo
c*** prop to icb if necessary
      if(imtch.eq.nic) call efprop(noc,nicp1,f,h,nint-2)
      do i=1,4
        zi(i)=s(i+jj-2)-f(i+jj-2)
      enddo
      return
      end

      subroutine solint(ls,f)
c*** does starting solutions for spheroidal mode integrals -- assume below turning point
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      dimension f(*)
      data d1,d2/.111111111111111d0,0.066666666666667d0/
      jfun=2*(kg+2)
c****
      p=fl+fl
      fnn=r(ls)/(p+1)
      fmm=r(ls)/(p-1)
      fpp=r(ls)/(p+2)
      fqq=r(ls)/(p+3)
c****
      qff=1.d0+xlam(ls)*fct
      qll=1.d0+qshear(ls)*fct
      qaa=1.d0+xa2(ls)*fct
      ro=rho(ls)
      gr=g(ls)
      ff=fcon(ls)*qff
      ll=lcon(ls)*qll
      if(ifanis.eq.0)then
        nn=ll
        cc=ff+ll+ll
        aa=cc
      else
        nn=ncon(ls)*qll
        cc=ccon(ls)*qaa
        aa=acon(ls)*qaa
      endif
      qrka=d1*(4.d0*(aa+ff-nn)+cc)*qkappa(ls)
      qrmu=d2*(aa+cc-2.d0*ff+5.d0*nn+6.d0*ll)*qshear(ls)
      z=r(ls)
      zr=1.d0/z
      do i=1,jfun
        f(i)=ar(i,ls)
      enddo
      f1=f(1)*zr
      f2=f(2)*zr
      f1sq=f1*f1
      f2sq=f2*f2
      g1=fl*f1
      g2=fl*f2
      e1=g1+f1-sfl3*f2
      e2=g2+sfl3*f1-2.d0*f2
      e3=2.d0*g1-4.d0*f1+sfl3*f2
c*** this integrand goes like r**(2l)
      f(jfun+1)=z*z*ro*(f1sq+f2sq)*fnn
c*** this integrand goes like r**(2l-2)
      f(jfun+3)=(qrka*e1*e1+qrmu*((fl3-2)*f2sq+e2*e2+e3*e3/3.d0))*fmm
      if(kg.eq.0) then 
c*** f(6) integrand goes like r**(2l-2)+ a part that goes like r**(2l)
        f(6)=(sfl3*(ll*f1sq+aa*f2sq)+f2*((2.d0*(nn-aa-ll)+ff)
     +     *f1-ff*g1)+ll*g2*f1)*fmm+f1*f2*z*ro*gr*fnn
c*** f(8) integrand goes like r**(2l-2)+ a part that goes like r**(2l)
        f(8)=((fl3*ll+4.d0*(aa-nn-ff)+cc)*f1sq+
     +    (4.d0*ll-nn-nn+fl3*aa)*f2sq+cc*g1*g1+
     +    ll*g2*g2 +2.d0*f1*(sfl3*(2.d0*(nn-aa-ll)+ff)*f2+
     +    (ff+ff-cc)*g1+sfl3*ll*g2)-2.d0*f2*(sfl3*ff*g1+(ll+ll)*g2))*fmm
     +    +2.d0*(2.d0*z*(z*ro-gr)*f1sq+f1*f2*sfl3*z*gr)*fnn*ro
      else
        f3=f(3)*zr
        g3=(fl+2.d0)*f3
c*** f(8) integrand goes like r**(2l-2)+ a part that goes like r**(2l)
c*** and a part that goes like r**(2l+1) and r**(2l+2)
        f(8)=(sfl3*(ll*f1sq+aa*f2sq)+f2*((2.d0*(nn-aa-ll)+ff)
     +     *f1-ff*g1)+ll*g2*f1)*fmm+f1*f2*z*ro*gr*fnn
     +   +f2*z*ro*f3*fpp+sfl3/(fl+0.5d0)*.25d0*f3*(g3+fl*f3)*fqq
c*** f(10) integrand goes like r**(2l-2)+ a part that goes like r**(2l)
c*** and a part that goes like r**(2l+1) and r**(2l+2)
        f(10)=((fl3*ll+4.d0*(aa-nn-ff)+cc)*f1sq+
     +    (4.d0*ll-nn-nn+fl3*aa)*f2sq+cc*g1*g1+
     +    ll*g2*g2 +2.d0*f1*(sfl3*(2.d0*(nn-aa-ll)+ff)*f2+
     +    (ff+ff-cc)*g1+sfl3*ll*g2)-2.d0*f2*(sfl3*ff*g1+(ll+ll)*g2))*fmm
     +    +2.d0*(2.d0*z*(z*ro-gr)*f1sq+f1*f2*sfl3*z*gr)*fnn*ro
     +    + (0.25d0*(fl*f3+g3)**2)*fqq+
     +    2.d0*ro*z*(f3*sfl3*f2+ f1*(g3-f3))*fpp
      endif
      return
      end

      subroutine efprop(jf,jl,f,h,nint)
c    fprop propagates the fundamental matrix f from jf to jl (a fluid region)
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension f(*),s(8),h(nint,*)
      d=fl3/wsq
      kk=kg+1
      jj=2*kk
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 if(i.eq.jl) return
      sum=0.
      do j=1,jj
        f(j)=ar(j,i)
        sum=sum+abs(f(j))
      enddo
      i=i+jud
      x=y
      y=r(i)
      if(y.eq.x) goto 5
      if(sum.eq.0.d0) goto 5
      iq=min(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-d*xi*xi)*rho(i)/flam(i)
      q=max(sfl3/x,sqrt(abs(alfsq-fl3/(x*x)))+1.d0/r(iq))
      del=jud*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q)
      dxs=dx
      do j=1,nint
        s(j)=f(j)
      enddo
      do ni=1,in
        z=x+c(ni)
        t=z-r(iq)
        zr=1.d0/z
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        gr=(g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq))))*zr
        qrka=ff*qkappa(iq)
        rogr=ro*gr
        t11=(gr*d-1.d0)*zr
        s11=-ro*(wsq+4.d0*gr-gr*gr*d)
        c11=-d*zr*zr/ro+1.d0/ff
        if(kg.eq.0) then
          s11=s11+4.d0*ro*ro
          h(1,ni)=t11*f(1)+c11*f(2)
          h(2,ni)=s11*f(1)-t11*f(2)
          f2=sfl3*(gr*f(1)-zr*f(2)/ro)/(wsq*z)
          f1=f(1)/z
          f1sq=f1*f1
          e1=h(1,ni)+f1-sfl3*f2
          rf=f1+f1-sfl3*f2
          h(3,ni)=z*z*ro*(f1sq+f2*f2)
          h(4,ni)=f2*(z*z*rogr*f1-e1*ff)
          h(5,ni)=qrka*e1*e1
          h(6,ni)=ff*e1*e1+2.d0*ro*z*z*f1*(2.d0*f1*ro-gr*rf)
          call rkdot(f,s,h,nint,ni)
        else
          t12=d*zr*zr
          t21=-4.d0*ro
          t22=-fl*zr
          s22=ro*t12
          s12=ro*(t11+t22)
          h(1,ni)=t11*f(1)+t12*f(2)+c11*f(3)
          h(2,ni)=t21*f(1)+t22*f(2)+4.d0*f(4)
          h(3,ni)=s11*f(1)+s12*f(2)-t11*f(3)-t21*f(4)
          h(4,ni)=s12*f(1)+s22*f(2)-t12*f(3)-t22*f(4)
          f2=sfl3*(gr*f(1)+zr*f(2)-zr*f(3)/ro)/(wsq*z)
          f1=f(1)/z
          f3=f(2)/z
          f1sq=f1*f1
          f2sq=f2*f2
          f3sq=f3*f3
          g3=h(2,ni)
          e1=h(1,ni)+f1-sfl3*f2
          rf=f1+f1-sfl3*f2
          h(5,ni)=z*z*ro*(f1sq+f2sq)
          h(6,ni)=f2*(z*z*rogr*f1-e1*ff)+z*ro*f2*f3
     +      +sfl3/(fl+0.5d0)*0.25*f3*(g3+fl*f3)
          h(7,ni)=qrka*e1*e1
          h(8,ni)=ff*e1*e1+2.d0*ro*z*z*f1*(2.d0*f1*ro-gr*rf)
     +             +0.25d0*(fl*fl*f3sq+g3*g3)+
     +        2.d0*f3*(z*ro*sfl3*f2+fl*.25d0*g3)+2.d0*f1*z*ro*(g3-f3)
          call rkdot(f,s,h,nint,ni)
        end if
      enddo
      x=y
      if(y.ne.r(i)) go to 15
      goto 5
      end

      subroutine esprop(jf,jl,f,h,nint)
c    sprop propagates the fundamental matrix f from jf to jl (a solid region)
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension f(*),s(10),h(nint,*)
      data d1,d2/.111111111111111d0,0.066666666666667d0/
      jj=2*(kg+2)
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 if(i.eq.jl) return
      sum=0.
      do j=1,jj
        f(j)=ar(j,i)
        sum=sum+abs(f(j))
      enddo
      i=i+jud
      x=y
      y=r(i)
      if(x.eq.y) goto 5
      if(sum.eq.0.) goto 5
      iq=min(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=sqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=.5d0*(alfsq+betasq+delsq)
      al=fl3/(x*x)
      aq=fksq-delsq-al
      qs=sqrt(abs(fksq-al))+1.d0/r(iq)
      qf=sqrt(abs(aq))+1.d0/r(iq)
      q=max(sfl3/x,qs,qf)
      del=jud*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q)
      dxs=dx
      do j=1,nint
        s(j)=f(j)
      enddo
      do 40 ni=1,in
      z=x+c(ni)
      t=z-r(iq)
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      if(ifanis.eq.0) then
        nn=ll
        cc=ff+ll+ll
        aa=cc
      else
        nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
        cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
        aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
      end if
      gg=2.d0*(nn-aa-ll)+ff
      qrka=d1*(4.d0*(aa+ff-nn)+cc)*qkappa(iq)
      qrmu=d2*(aa+cc-2.d0*ff+5.d0*nn+6.d0*ll)*qshear(iq)
      zr=1.d0/z
      sfl3z=sfl3*zr
      rogr=ro*gr
      c11=1.d0/cc
      c22=1.d0/ll
      dmg=aa-nn-ff*ff*c11
      zdmg=zr*dmg
      t11=-2.d0*ff*zr*c11+zr
      t12=sfl3z*ff*c11
      t21=-sfl3z
      t22=zr+zr
      s11=-ro*wsq+4.d0*zr*(zdmg-rogr)
      s22=-ro*wsq+zr*zr*(fl3*(dmg+nn)-nn-nn)
      s12=sfl3z*(rogr-zdmg-zdmg)
      if(kg.eq.0) then
        s11=s11+4.d0*ro*ro
        h(1,ni)=t11*f(1)+t12*f(2)+c11*f(3)
        h(2,ni)=t21*f(1)+t22*f(2)+c22*f(4)
        h(3,ni)=s11*f(1)+s12*f(2)-t11*f(3)-t21*f(4)
        h(4,ni)=s12*f(1)+s22*f(2)-t12*f(3)-t22*f(4)
        f1=f(1)*zr
        f2=f(2)*zr
        f1sq=f1*f1
        f2sq=f2*f2
        g1=h(1,ni)
        g2=h(2,ni)
        h(5,ni)=z*z*ro*(f1sq+f2sq)
        h(6,ni)=sfl3*(ll*f1sq+aa*f2sq)+f2*((z*rogr+gg)*f1-ff*g1)
     +         +ll*g2*f1
        e1=g1+f1-sfl3*f2
        e2=g2+sfl3*f1-2.d0*f2
        e3=2.d0*g1-4.d0*f1+sfl3*f2
        rf=f1+f1-sfl3*f2
        h(7,ni)=qrka*e1*e1+qrmu*((fl3-2)*f2sq+e2*e2+e3*e3/3.d0)
        h(8,ni)=(fl3*ll+4.d0*(z*ro*(z*ro-gr)+aa-nn-ff)+cc)*f1sq+
     +    (4.d0*ll-nn-nn+fl3*aa)*f2sq+cc*g1*g1+ll*g2*g2 +
     +     2.d0*f1*(sfl3*(z*rogr+gg)*f2+(ff+ff-cc)*g1+sfl3*ll*g2)-
     +       2.d0*f2*(sfl3*ff*g1+(ll+ll)*g2)
        call rkdot(f,s,h,nint,ni)
      else
        t31=-4.d0*ro
        t33=-fl*zr
        s13=-fl1*zr*ro
        s23=ro*sfl3z
        h(1,ni)=t11*f(1)+t12*f(2)+c11*f(4)
        h(2,ni)=t21*f(1)+t22*f(2)+c22*f(5)
        h(3,ni)=t31*f(1)+t33*f(3)+4.d0*f(6)
        h(4,ni)=s11*f(1)+s12*f(2)+s13*f(3)-t11*f(4)-t21*f(5)-t31*f(6)
        h(5,ni)=s12*f(1)+s22*f(2)+s23*f(3)-t12*f(4)-t22*f(5)
        h(6,ni)=s13*f(1)+s23*f(2)-t33*f(6)
        f1=f(1)*zr
        f2=f(2)*zr
        f3=f(3)*zr
        f1sq=f1*f1
        f2sq=f2*f2
        g1=h(1,ni)
        g2=h(2,ni)
        g3=h(3,ni)
        h(7,ni)=z*z*ro*(f1sq+f2sq)
        h(8,ni)=sfl3*(ll*f1sq+aa*f2sq)+f2*((z*rogr+gg)*f1-ff*g1)+
     +     ll*g2*f1
     +   +f2*z*ro*f3+sfl3/(fl+0.5d0)*.25d0*f3*(g3+fl*f3)
        e1=g1+f1-sfl3*f2
        e2=g2+sfl3*f1-2.d0*f2
        e3=2.d0*g1-4.d0*f1+sfl3*f2
        h(9,ni)=qrka*e1*e1+qrmu*((fl3-2)*f2sq+e2*e2+e3*e3/3.d0)
        h(10,ni)=(fl3*ll+4.d0*(z*ro*(z*ro-gr)+aa-nn-ff)+cc)*f1sq+
     +    (4.d0*ll-nn-nn+fl3*aa)*f2sq+cc*g1*g1+ll*g2*g2 +
     +     2.d0*f1*(sfl3*(z*rogr+gg)*f2+(ff+ff-cc)*g1+sfl3*ll*g2)-
     +         2.d0*f2*(sfl3*ff*g1+(ll+ll)*g2)
     +    + 0.25d0*(fl*f3+g3)**2+2.d0*z*ro*(f3*sfl3*f2+f1*(g3-f3))
        call rkdot(f,s,h,nint,ni)
      end if
   40 continue
      x=y
      if(y.ne.r(i)) go to 15
      goto 5
      end

      subroutine norm(f,iexp,nvec)
      implicit real*8(a-h,o-z)
      dimension f(*)
      data econst/1048576.d0/
      size=abs(f(1))
      do j=2,nvec
        size=max(size,abs(f(j)))
      enddo
   20 if(size.lt.1024.d0) return
        do j=1,nvec
          f(j)=f(j)/econst
        enddo
        size=size/econst
        iexp=iexp+20
        goto 20
      end

      subroutine sprp_jhw(jf,jl,f,h,nvesm,iexp)
c*** propagate a minor vector in a solid region from level jf to jl ***
C*** THIS VERSION IS USED ONLY IN THE DOWN DIRECTION and makes the
c*** eigenvector using the woodhouse 1988 algorithm. 
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      common/stdep/fick,ls
      dimension f(*),h(nvesm,*),s(14),fp(14),rne(6)
      data econst/1048576.d0/
      jfun=2*(kg+2)
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 call norm(f,iexp,nvesm)
      inorm(i)=inorm(i)+iexp
c*** these are diagonal elements of X so rne(1)=(+-)b1**2, etc
c*** we shall choose the largest to determine which row of X to use
c*** this will be one of the tractions unless we are very close to the free surface
      if(kg.eq.0) then
        dum1=ar(1,i)*f(5)-ar(5,i)*f(1)
        dum2=ar(3,i)*f(4)-ar(4,i)*f(3)
        rne(1)=ar(1,i)*f(3)-ar(3,i)*f(1)
        rne(2)=ar(4,i)*f(1)-ar(1,i)*f(4)
        rne(3)=ar(5,i)*f(4)-ar(4,i)*f(5)
        rne(4)=ar(3,i)*f(5)-ar(5,i)*f(3)
      else
        dum1=ar(12,i)*f(14)-ar(14,i)*f(12)
        dum2=ar(7,i)*f(2)-ar(2,i)*f(7)
        dum3=ar(11,i)*f(13)-ar(13,i)*f(11)
c*** don't bother with b3 or b6 as always small
        rne(1)=ar(6,i)*f(3)+ar(3,i)*f(6)-ar(1,i)*f(8)-ar(8,i)*f(1)+
     +      2.d0*ar(13,i)*f(13)
        rne(2)=ar(1,i)*f(9)+ar(9,i)*f(1)-ar(6,i)*f(4)-ar(4,i)*f(6)+
     +      2.d0*ar(14,i)*f(14)
        rne(3)=ar(4,i)*f(10)+ar(10,i)*f(4)-ar(9,i)*f(5)-ar(5,i)*f(9)+
     +      2.d0*ar(11,i)*f(11)
        rne(4)=ar(8,i)*f(5)+ar(5,i)*f(8)-ar(3,i)*f(10)-ar(10,i)*f(3)+
     +      2.d0*ar(12,i)*f(12)
      end if
      bmax=0.d0
      do j=1,4
        if(abs(rne(j)).gt.bmax) then
            bmax=abs(rne(j))
            jmax=j
        end if
      end do
      fnr=1.d0/sqrt(bmax)
c**choose best row to construct eigenvector
      if(kg.eq.0) then
        if(jmax.eq.1) then
          rne(2)=ar(2,i)*f(1)-ar(1,i)*f(2)
          rne(3)=0.5d0*(dum1+dum2)
          rne(4)=ar(2,i)*f(3)-ar(3,i)*f(2)
        else if(jmax.eq.2) then
          rne(1)=ar(2,i)*f(1)-ar(1,i)*f(2)
          rne(3)=ar(4,i)*f(2)-ar(2,i)*f(4)
          rne(4)=0.5d0*(dum1-dum2)
        else if(jmax.eq.3) then
          rne(1)=0.5d0*(dum1+dum2)
          rne(2)=ar(4,i)*f(2)-ar(2,i)*f(4)
          rne(4)=ar(2,i)*f(5)-ar(5,i)*f(2)
        else if(jmax.eq.4) then
          rne(1)=ar(2,i)*f(3)-ar(3,i)*f(2)
          rne(2)=0.5d0*(dum1-dum2)
          rne(3)=ar(2,i)*f(5)-ar(5,i)*f(2)
        end if
      else
        if(jmax.eq.1) then
          rne(2)=ar(1,i)*f(7)+ar(7,i)*f(1)-ar(6,i)*f(2)-ar(2,i)*f(6)-
     +           ar(14,i)*f(13)-ar(13,i)*f(14)
          rne(3)=ar(6,i)*f(12)+ar(12,i)*f(6)-ar(13,i)*f(7)-ar(7,i)*
     +     f(13)+ar(14,i)*f(8)+ar(8,i)*f(14)
          rne(4)=dum1+dum2-2.d0*ar(13,i)*f(11)+ar(6,i)*f(5)-ar(1,i)*
     +     f(10)+ar(8,i)*f(4)-ar(3,i)*f(9)
          rne(5)=ar(7,i)*f(3)+ar(3,i)*f(7)-ar(2,i)*f(8)-ar(8,i)*f(2)+
     +           ar(12,i)*f(13)+ar(13,i)*f(12)
          rne(6)=ar(14,i)*f(3)+ar(3,i)*f(14)-ar(13,i)*f(2)-ar(2,i)*
     +     f(13)+ar(12,i)*f(1)+ar(1,i)*f(12)
        else if(jmax.eq.2) then
          rne(1)=ar(1,i)*f(7)+ar(7,i)*f(1)-ar(6,i)*f(2)-ar(2,i)*f(6)-
     +           ar(14,i)*f(13)-ar(13,i)*f(14)
          rne(3)=ar(6,i)*f(11)+ar(11,i)*f(6)-ar(14,i)*f(7)-ar(7,i)*
     +     f(14)-ar(13,i)*f(9)-ar(9,i)*f(13)
          rne(4)=ar(9,i)*f(2)+ar(2,i)*f(9)-ar(4,i)*f(7)-ar(7,i)*f(4)+
     +           ar(14,i)*f(11)+ar(11,i)*f(14)
          rne(5)=dum2+dum3-2.d0*ar(14,i)*f(12)+ar(6,i)*f(5)-ar(1,i)*
     +     f(10)-ar(4,i)*f(8)+ar(9,i)*f(3)
          rne(6)=ar(1,i)*f(11)+ar(11,i)*f(1)-ar(14,i)*f(2)-ar(2,i)*
     +     f(14)-ar(13,i)*f(4)-ar(4,i)*f(13)
        else if(jmax.eq.3) then
          rne(4)=rne(3)
          rne(1)=-dum1-dum2-2.d0*ar(11,i)*f(13)+ar(5,i)*f(6)-ar(10,i)*
     +      f(1)+ar(4,i)*f(8)-ar(9,i)*f(3)
          rne(2)=ar(9,i)*f(2)+ar(2,i)*f(9)-ar(4,i)*f(7)-ar(7,i)*f(4)+
     +           ar(14,i)*f(11)+ar(11,i)*f(14)
          rne(3)=ar(7,i)*f(11)+ar(11,i)*f(7)-ar(9,i)*f(12)-ar(12,i)*
     +     f(9)+ar(14,i)*f(10)+ar(10,i)*f(14)
          rne(5)=ar(7,i)*f(5)+ar(5,i)*f(7)-ar(12,i)*f(11)-ar(11,i)*
     +     f(12)-ar(2,i)*f(10)-ar(10,i)*f(2)
          rne(6)=ar(2,i)*f(11)+ar(11,i)*f(2)-ar(4,i)*f(12)-ar(12,i)*
     +     f(4)+ar(14,i)*f(5)+ar(5,i)*f(14)
        else if(jmax.eq.4) then
          rne(5)=rne(4)
          rne(1)=ar(7,i)*f(3)+ar(3,i)*f(7)-ar(2,i)*f(8)-ar(8,i)*f(2)+
     +           ar(12,i)*f(13)+ar(13,i)*f(12)
          rne(2)=-dum3-dum2-2.d0*ar(12,i)*f(14)+ar(5,i)*f(6)-ar(10,i)
     +     *f(1)-ar(8,i)*f(4)+ar(3,i)*f(9)
          rne(3)=ar(7,i)*f(12)+ar(12,i)*f(7)+ar(8,i)*f(11)+ar(11,i)*
     +      f(8)+ar(13,i)*f(10)+ar(10,i)*f(13)
          rne(4)=ar(7,i)*f(5)+ar(5,i)*f(7)-ar(12,i)*f(11)-ar(11,i)*
     +     f(12)-ar(2,i)*f(10)-ar(10,i)*f(2)
          rne(6)=ar(2,i)*f(12)+ar(12,i)*f(2)+ar(3,i)*f(11)+ar(11,i)*
     +      f(3)+ar(13,i)*f(5)+ar(5,i)*f(13)
        end if
      end if
      do j=1,jfun
        ar(j,i)=rne(j)*fnr
      enddo
c****hmmmm -- technically took sqrt
      inorm(i)=inorm(i)/2
c*** fix sign
      if(i.eq.jf.and.jf.ne.n) then
          jm=2+kg
          if(ar(jm+1,jf)*ar(jm,jf+1).lt.0.d0) then
             do jj=1,jfun
               ar(jj,jf)=-ar(jj,jf)
             enddo
          end if
      end if
      if(i.lt.jf) then
         sum=0.d0
         sum2=0.d0
         idiff=inorm(i)-inorm(i+1)
         dnrm=1.d0
         if(idiff.ne.0) dnrm=2.d0**idiff
         do jj=1,jfun
             sum=sum+(ar(jj,i+1)-ar(jj,i)*dnrm)**2
             sum2=sum2+(ar(jj,i+1)+ar(jj,i)*dnrm)**2
         enddo
         if(sum2.lt.sum) then
           do jj=1,jfun
              ar(jj,i)=-ar(jj,i)
           enddo
         endif
      end if
      if(i.eq.jl) return
      i=i+jud
      x=y
      y=r(i)
      if(y.eq.x) goto 5
      iq=min(i,i-jud)
      zs=min(x,y)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=sqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=.5d0*(alfsq+betasq+delsq)-fl3/(x*x)
      qt=sqrt(abs(fksq))+sqrt(abs(fksq-delsq))+2.d0/zs
      q=qt+float(kg)*sfl3/x
      del=float(jud)*step(maxo)/q
      dxs=0.d0
   20   do j=1,nvesm
          s(j)=f(j)
        enddo
        y=x+del
        if(jud*(y-r(i)).gt.0.d0) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q)
        dxs=dx
        do ni=1,in
          z=x+c(ni)
          call derms(iq,z,f,h(1,ni),qff,qll,qaa)
          call rkdot(f,s,h,nvesm,ni)
        enddo
        x=y
      if(y.ne.r(i)) goto 20
      goto 5
      end

      subroutine fprp_jhw(jf,jl,f,h,nvefm,iexp)
c*** propagate the minor vector in a fluid region from level jf to jl ***
C*** THIS VERSION IS USED ONLY IN THE DOWN DIRECTION and makes the
c*** eigenvector using the woodhouse 1988 algorithm. This technique seems
c*** to work best for problem modes.
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      real*8 lcon,ncon,lspl,nspl
      common r(mk),fmu(mk),flam(mk),qshear(mk),qkappa(mk),
     + xa2(mk),xlam(mk),rho(mk),qro(3,mk),g(mk),qg(3,mk),
     + fcon(mk),fspl(3,mk),lcon(mk),lspl(3,mk),ncon(mk),
     + nspl(3,mk),ccon(mk),cspl(3,mk),acon(mk),aspl(3,mk)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/shanks/b(46),c(10),dx,step(8),maxo,in
      dimension f(*),h(nvefm,*),s(5),fp(5),rne(4)
      jfun=2*(kg+1)
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 call norm(f,iexp,nvefm)
      if(kg.eq.1) then
        inorm(i)=inorm(i)+iexp
c*** these are diagonal elements of X so rne(1)=(+-)b1**2, etc
c*** we shall choose the largest to determine which row of X to use
        dum1=ar(1,i)*f(5)-ar(5,i)*f(1)
        dum2=ar(3,i)*f(4)-ar(4,i)*f(3)
        rne(1)=ar(1,i)*f(3)-ar(3,i)*f(1)
        rne(2)=ar(4,i)*f(1)-ar(1,i)*f(4)
        rne(3)=ar(5,i)*f(4)-ar(4,i)*f(5)
        rne(4)=ar(3,i)*f(5)-ar(5,i)*f(3)
        bmax=0.d0
        do j=1,jfun
          if(abs(rne(j)).gt.bmax) then
              bmax=abs(rne(j))
              jmax=j
          end if
        end do
        fnr=1.d0/sqrt(bmax)
c**choose best row to construct eigenvector -- this will nearly always be jmax=3
        if(jmax.eq.1) then
          rne(2)=ar(2,i)*f(1)-ar(1,i)*f(2)
          rne(3)=0.5d0*(dum1+dum2)
          rne(4)=ar(2,i)*f(3)-ar(3,i)*f(2)
        else if(jmax.eq.2) then
          rne(1)=ar(2,i)*f(1)-ar(1,i)*f(2)
          rne(3)=ar(4,i)*f(2)-ar(2,i)*f(4)
          rne(4)=0.5d0*(dum1-dum2)
        else if(jmax.eq.3) then
          rne(1)=0.5d0*(dum1+dum2)
          rne(2)=ar(4,i)*f(2)-ar(2,i)*f(4)
          rne(4)=ar(2,i)*f(5)-ar(5,i)*f(2)
        else if(jmax.eq.4) then
          rne(1)=ar(2,i)*f(3)-ar(3,i)*f(2)
          rne(2)=0.5d0*(dum1-dum2)
          rne(3)=ar(2,i)*f(5)-ar(5,i)*f(2)
        end if
        do j=1,4
          ar(j,i)=rne(j)*fnr
        enddo
        ar(5,i)=0.
        ar(6,i)=0.
c****hmmmm -- technically took sqrt
        inorm(i)=inorm(i)/2
c*** fix sign
        if(i.eq.jf.and.jf.ne.n) then
          if(ar(3,jf)*ar(4,jf+1).lt.0.d0) then
              do jj=1,4
                ar(jj,jf)=-ar(jj,jf)
             enddo
          end if
        end if
        if(i.lt.jf) then
           sum=0.d0
           sum2=0.d0
           idiff=inorm(i)-inorm(i+1)
           dnrm=1.d0
           if(idiff.ne.0) dnrm=2.d0**idiff
           do jj=1,4
               sum=sum+(ar(jj,i+1)-ar(jj,i)*dnrm)**2
               sum2=sum2+(ar(jj,i+1)+ar(jj,i)*dnrm)**2
           enddo
           if(sum2.lt.sum) then
             do jj=1,4
                ar(jj,i)=-ar(jj,i)
             enddo
           endif
        end if
c*** this is for kg=0 -- just integrate across the outer core or ocean
      else
        inorm(i)=iexp
        ar(1,i)=f(1)
        ar(2,i)=f(2)
        ar(3,i)=0.
        ar(4,i)=0.
      end if
      if(i.eq.jl) return
      i=i+jud
      x=y
      y=r(i)
      if(y.eq.x) goto 5
      iq=min(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=sqrt(abs(alfsq-fl3/(x*x)))+1.d0/r(iq)+float(kg)*sfl3/x
      del=float(jud)*step(maxo)/q
      dxs=0.d0
   15   y=x+del
        if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q)
        dxs=dx
        do j=1,nvefm
          s(j)=f(j)
        enddo
        do ni=1,in
          z=x+c(ni)
          call dermf(iq,z,f,h(1,ni),qff)
          call rkdot(f,s,h,nvefm,ni)
        enddo
        x=y
        if(y.ne.r(i)) go to 15
      goto 5
      end

      subroutine remedy(ls,imtch)
c    obtains the eigenfunction of an awkward spheroidal mode by using the technique described in
c    j woodhouse, 1988 (in seismological algorithms, ed doornbos, academic press, pp321--370)
c*** we have found it best to use the original eigenfunction below the match level, and this eigenfunction
c*** above it.
c*** this works extremely well for almost all IC modes when kg=0, can get to about l=500 f0r cmb stoneleys
      implicit real*8(a-h,o-z)
      parameter (mk=950)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,mk),inorm(mk),jjj(mk)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/upx/up(14,mk),iup(mk)
      dimension ass(14),h(140),asav(6)
c      print *,'in remedy  ',nord,l 
      nvefm=2+kg*3
      nvesm=5+kg*9
      jfun=2*(kg+2)
      jfl=jfun-2
c*** copy original upgoing solution in the mantle
      do i=nocp1,n
         inorm(i)=iup(i)
         do j=1,nvesm
           ar(j,i)=up(j,i)
         enddo
      enddo
c*** need to set iback=2 so that derms behaves correctly
      iback=2
      ar(5,n)=0.d0
      do i=1,nvesm
        ass(i)=0.d0
      enddo
      if(kg.eq.1.and.nsl.eq.n) then
         ass(6)=1.d0
      else
         ass(1)=1.d0
      end if
      iexp=0
      if(nsl.ne.n) then
c*** propagate through ocean if there is one 
        call fprp_jhw(n,nslp1,ass,h,nvefm,iexp)
        call fsbm(ass,kg,iback)
      end if
c*** propagate through mantle
      nto=max(ls,nocp1)
      call sprp_jhw(nsl,nto,ass,h,nvesm,iexp)
      if(nsl.ne.n) then
         if(kg.eq.0) then
           con=ar(2,nslp1)/ar(3,nsl)
           do i=nto,nsl
               do j=1,4
                 ar(j,i)=ar(j,i)*con
               enddo
           enddo
         end if
         idiff=inorm(nslp1)-inorm(nsl)
         do i=nto,nsl
             inorm(i)=inorm(i)+idiff
         enddo
      end if
c*** the next never happens
      if(nto.gt.nocp1) return
c next case for CMB stoneley where upgoing solution should be fine
      if(ls.gt.nic) then
c*** now we have  fluid vector (ass) at noc so match to upgoing solution
        if(kg.eq.0) then
          ass(1)=ar(1,nocp1)
          ass(2)=ar(3,nocp1)
        else
          ass(1)=ar(1,nocp1)
          ass(2)=ar(3,nocp1)
          ass(3)=ar(4,nocp1)
          ass(4)=ar(6,nocp1)
        end if
        imtch=noc
c*** now do match
        sum1=0.d0
        sum2=0.d0
        do j=1,jfl
          sum1=sum1+ass(j)*ar(j,noc)
          sum2=sum2+ar(j,noc)**2
        enddo
        c=sum1/sum2
        idiff=inorm(nocp1)-inorm(noc)
        do i=ls,noc
          inorm(i)=inorm(i)+idiff
          do j=1,jfl
            ar(j,i)=c*ar(j,i)
          enddo
        enddo
      else
c*** here for IC mode
        imtch=nic
c*** save old solution just above icb to do the match
        do j=1,jfl
          asav(j)=ar(j,nicp1)
        enddo
c***  propagate from cmb to icb -- copy original upgoing solution first
         do i=nicp1,noc
            inorm(i)=iup(i)
            do j=1,nvesm
              ar(j,i)=up(j,i)
            enddo
         enddo
        iexp=0
        call sfbm(ass,kg,iback)
        call fprp_jhw(noc,nicp1,ass,h,nvefm,iexp)
        idiff=inorm(nocp1)-inorm(noc)
        do i=nicp1,noc
          inorm(i)=inorm(i)+idiff
        enddo
c*** now do match
        sum1=0.d0
        sum2=0.d0
        do j=1,jfl
          sum1=sum1+ar(j,nicp1)*asav(j)
          sum2=sum2+asav(j)**2
        enddo
        c=sum1/sum2
        idiff=inorm(nicp1)-inorm(nic)
        do i=ls,nic
          inorm(i)=inorm(i)+idiff
          do jj=1,jfun
            ar(jj,i)=c*ar(jj,i)
          enddo
        enddo
      end if
      iback=0
      return
      end
