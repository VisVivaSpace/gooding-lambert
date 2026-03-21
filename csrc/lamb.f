C----=------------------------------------------------------------------ 
      subroutine vlamb (gm,r1,r2,th,tdelt,n,vr11,vt11,vr12,vt12,vr21,
     &    vt21,vr22,vt22)
C----=------------------------------------------------------------------ 
      implicit none
      integer i,n,m
      real*8 gm,r1,r2,th,tdelt,vr11,vt11,vr12,vt12,vr21,vt21,vr22,vt22
      real*8 thr2,dr,r1r2,r1r2th,csq,c,s,gms,qsqfm1,q,rho,sig,t,x1,x2
      real*8 x,unused,qzminx,qzplx,zplqx,vt1,vr1,vr2,vt2
      real*8 PI,TWOPI
      parameter (PI = 3.141592653589793D0, TWOPI = 2D0*PI)

      m = th/TWOPI
      thr2 = th/2D0 - m*PI
      dr = r1 -r2
      r1r2 = r1*r2
      r1r2th = 4D0*r1r2*dsin(thr2)**2
      csq = dr**2 + r1r2th
      c = dsqrt(csq)
      s = (r1+r2+c)/2D0
      gms = dsqrt(gm*s/2D0)
      qsqfm1 = c/s
      q = dsqrt(r1r2)*dcos(thr2)/s
      
      if (c.ne.0D0) then
        rho = dr/c
        sig = r1r2th/csq
      else
        rho = 0D0
        sig = 1D0
      end if

      t = 4D0*gms*tdelt/s**2

      call xlamb(m,q,qsqfm1,t,n,x1,x2)
C        proceed for either a single soultion or a pair

      do 100 i = 1,n
        if (i.eq.1) then
          x=x1
        else
          x=x2
        end if

        call tlamb (m,q,qsqfm1,x,-1,unused,qzminx,qzplx,zplqx)
        vt2 = gms*zplqx*dsqrt(sig)
        vr1 = gms*(qzminx-qzplx*rho)/r1
        vt1 = vt2/r1
        vr2 = -gms*(qzminx+qzplx*rho)/r2
        vt2 = vt2/r2
        if (i.eq.1) then
          vr11 = vr1
          vt11 = vt1 
          vr12 = vr2
          vt12 = vt2 
        else
          vr21 = vr1
          vt21 = vt1 
          vr22 = vr2
          vt22 = vt2 
        end if
 100  CONTINUE

      return
      END
C----=------------------------------------------------------------------ 


C----=------------------------------------------------------------------ 
      subroutine xlamb (m,q,qsqfm1,tin,n,x,xpl)
C----=------------------------------------------------------------------ 
      implicit none
      integer i,n,m
      real*8 q,qsqfm1,tin,x,xpl
      real*8 thr2,t0,dt,d2t,d3t,tdiff,w,xm,tmin,xmold,xtest,tdiffm,d2t2
      real*8 t,tdiff0
      real*8 PI,TOL,C0,C1,C2,C3,C41,C42,D8RT
      parameter (PI = 3.141592653589793D0, TOL=1D-10, C0=1.7D0)
      parameter (C1=0.5D0, C2=0.03D0, C3 = 0.15D0, C41=1D0, C42=0.24D0)
      D8RT(X) = dsqrt(dsqrt(dsqrt(x)))
      
      thr2 = datan2(qsqfm1,2D0*q)/PI

      if (m.eq.0) then
C         single-rev starter from T (at X = 0) & bilinear (usually)
        n=1
        call tlamb(m,q,qsqfm1,0D0,0,t0,dt,d2t,d3t)
        tdiff = tin-t0

        if (tdiff.le.0D0) then 
          x = t0*tdiff/(-4D0*tin)
C           (-4 is the value of dt, for x =0)
        else
          x = -tdiff/(tdiff+4D0)
          w = x+c0*dsqrt(2D0*(1D0-thr2))
          if (w.lt.0D0)
     &       x = x-dsqrt(D8RT(-w))*(x+dsqrt(tdiff/(tdiff+1.5D0*t0)))
          w = 4D0/(4D0+tdiff)
          x = x*(1D0+x*(c1*w-c2*x*dsqrt(w))) 
        end if
      else
C      with multirevs, first get t(min) as basis for starter
        xm = 1D0/(1.5D0*(M+5D-1)*PI)
        if (thr2.lt.5D-1) xm = D8RT(2D0*thr2)*xm
        if (thr2.gt.5D-1) xm = (2D0 - D8RT(2D0-2D0*thr2))*XM

C        (starter for tmin)
        do 110 i = 1,16
          call tlamb (m,q,qsqfm1,xm,3,tmin,dt,d2t,d3t)
          if (d2t.eq.0D0) go to 120
          xmold = xm
          xm = xm - dt*d2t/(d2t*d2t-dt*d3t/2D0)
          xtest = dabs(xmold/xm-1D0)
          if (xtest.le.TOL) go to 120
 110    CONTINUE
      
        n = -1
        return
C         (break off & exit if tmin not located - should never happen)
C         now proceed from t(min) to full starter

 120    CONTINUE

        tdiffm = tin - tmin
        if (tdiffm.lt.0D0) then
          n = 0 
          return
C           (exit if no solution with this m)
        else if (tdiffm.eq.0D0) then
          x = xm
          n = 1
          return
C           (exit if unique solution already from x(tmin))
        else
          n = 3
          if (d2t.eq.0D0) d2t = 6D0*M*PI
          x = dsqrt(tdiffm/(d2t/2d0+tdiffm/(1D0-xm)**2))
          w = xm+x
          w = w*4D0/(4D0+tdiffm)+(1D0-w)**2
          x = x*(1D0-(1D0+M+c41*(thr2-0.5D0))/(1D0+C3*M)*
     &      x*(c1*w+c2*x*dsqrt(w)))+xm
          d2t2 = d2t/2D0
          if (x.ge.1D0) then
            n=1
            go to 150
          end if
C          (no finite solution with X > XM)
        end if 
C
      end if

C          (now have a starter, so proceed by halley)

 130  CONTINUE
     
      do 140 i =1,3
      call tlamb (m,q,qsqfm1,x,2,t,dt,d2t,d3t)
      t = tin-t
      if (dt.ne.0D0) x = x+t*dt/(dt*dt+t*d2t/2D0)
  
 140  CONTINUE
 
      if (n.ne.3) return
C       (exit if only one solution, normally when m = 0)

      n=2
      xpl = x
C       (second multi-rev starter)

 150  call tlamb(m,q,qsqfm1,0D0,0,t0,dt,d2t,d3t)
      tdiff0 = t0-tmin
      tdiff = tin-t0
      if (tdiff.le.0) then
        x = xm-dsqrt(tdiffm/(d2t2-tdiffm*(d2t2/tdiff0
     &                - 1D0/XM**2)))
      else
        x = -tdiff/(tdiff+4D0)
        w = x+c0*dsqrt(2D0*(1D0-thr2)) 
        if (w.lt.0D0) x = 
     &   x-dsqrt(D8RT(-w))*(x+dsqrt(tdiff/(tdiff+1.5D0*t0)))
        w = 4D0/(4D0+tdiff)
        x = x*(1D0+(1D0+M+c42*(thr2-0.5D0))/(1D0+c3*m)*
     &             x*(c1*w-c2*x*dsqrt(w)))
        if (x.le.-1D0) then
          n = n-1
C           (no finite solution with X < XM)
          if (n.eq.1) x = xpl
        end if
      end if

      go to 130

      END
C----=------------------------------------------------------------------ 



C----=------------------------------------------------------------------ 
      subroutine tlamb (m,q,qsqfm1,x,n,t,dt,d2t,d3t)
C----=------------------------------------------------------------------ 
      implicit none
      integer i,n,m,p
      real*8 q,qsqfm1,x,t,dt,d2t,d3t
      real*8 qsq,xsq,y,z,qx,a,b,aa,bb,g,f,fg1,term,fg1sq,twoi1,told
      real*8 qz,qz2,u0i,u1i,u2i,u3i,tq,tqsum,ttmold,tterm,tqterm,u
      real*8 PI,SW
      logical LM1,L1,L2,L3
      parameter (PI = 3.141592653589793D0, SW = 0.4D0)
 
      LM1 = n.eq.-1
      L1 = n.ge.1
      L2 = n.ge.2
      L3 = n.eq.3
      qsq = q*q
      xsq = x*x
      u = (1D0-x)*(1D0+x)
      if (.not.LM1) then
C         (needed if series, and otherwise useful when z=0)
        dt = 0D0
        d2t = 0D0
        d3t = 0D0
      end if
      if (LM1.or.m.gt.0.or.x.lt.0D0.or.dabs(u).gt.SW) then
C         direct computation (not series)
        y = dsqrt(dabs(u)) 
        z = dsqrt(qsqfm1+qsq*xsq) 
        qx = q*x
        if (qx.le.0D0) then
          a = z-qx
          b = q*z-x
        end if
        if (qx.lt.0D0.and.LM1) then
          aa = qsqfm1/a
          bb = qsqfm1*(qsq*u-xsq)/b
        end if
        if (qx.eq.0D0.and.LM1.or.qx.gt.0D0) then
          aa = z+qx
          bb = q*z+x
        end if
        if (qx.gt.0D0) then
          a = qsqfm1/aa
          b = qsqfm1*(qsq*u-xsq)/bb
        end if
        if (.not.LM1) then
          if (qx*u.ge.0D0) then 
            g = x*z +q*u
          else
            g = (xsq-qsq*u)/(x*z-q*u)
          end if
          f = a*y
          if (x.le.1D0) then 
            t = m*PI + datan2(f,g)
          else
            if (f.gt.SW) then
              t = dlog(f+g)
            else
              fg1 = f/(g+1D0)
              term = 2D0*fg1
              fg1sq = fg1*fg1
              t = term
              twoi1 = 1D0

 200          twoi1 = twoi1 + 2D0
              term = term*fg1sq
              told = t
              t = t+term/twoi1
              if (t.ne.told) go to 200
C        (continue looping for inverse tanh)

            end if
          end if
          t = 2D0*(t/y +b)/u
          if (L1.and.z.ne.0D0) then
            qz = q/z
            qz2 = qz*qz
            qz = qz*qz2
            dt = (3D0*x*t -4D0*(a+qx*qsqfm1)/z)/u
            if (L2) d2t = (3D0*t+5D0*x*dt+4D0*qz*qsqfm1)/u
            if (L3) d3t = (8D0*dt+7D0*x*d2t-12D0*qz*qz2*x*qsqfm1)/u
          end if
        else
          dt = b
          d2t = bb
          d3t = aa
        end if
      else
C   compute by series
        u0i = 1D0
        if (L1) u1i = 1D0
        if (L2) u2i = 1D0
        if (L3) u3i = 1D0
        term = 4D0
        tq = q*qsqfm1
        i = 0
        if (q.lt.5D-1) tqsum = 1D0-q*qsq
        if (q.ge.5D-1) tqsum = (1D0/(1D0+q)+q)*qsqfm1
        ttmold = term / 3D0
        t = ttmold*tqsum

C         (start of loop)
 210    i = i+1
        p = i
        u0i = u0i*u
        if (L1.and.i.gt.1) u1i = u1i*u
        if (L2.and.i.gt.2) u2i = u2i*u
        if (L3.and.i.gt.3) u3i = u3i*u
        term = term*(p-0.5D0)/P
        tq = tq*qsq
        tqsum = tqsum+tq
        told = t
        tterm = term/(2D0*p+3D0) 
        tqterm = tterm*tqsum
        t = t-u0i*((1.5D0*p+0.25D0)*TQTERM/(p*p-0.25D0)
     &         -ttmold*tq)
        ttmold = tterm
        tqterm = tqterm*p
        if (L1) dt = dt+tqterm*u1i
        if (L2) d2t = d2t+tqterm*u2i*(p-1D0)
        if (L3) d3t = d3t+tqterm*u3i*(p-1D0)*(p-2D0)
        if (i.lt.n.or.t.ne.told) go to 210
C         (end of loop)
    
        if (L3) d3t = 8D0*x*(1.5D0*d2t-xsq*d3t)
        if (L2) d2t = 2D0*(2D0*xsq*d2t-dt)
        if (L1) dt = -2D0*x*dt
        t = t/xsq
      end if

      return
      END
C----=------------------------------------------------------------------ 


C----=------------------------------------------------------------------ 
C     END OF FILE
C----=------------------------------------------------------------------ 
