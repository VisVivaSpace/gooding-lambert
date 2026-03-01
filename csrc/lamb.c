#define PI 3.141592653589793
#define TWOPI 6.283185307179586
#define SW 0.4
#define TOL 1E-10
#define C0 1.7
#define C1 0.5
#define C2 0.03
#define C3 0.15
#define C41 1.0
#define C42 0.24
#define sign(x) (( x > 0 ) - ( x < 0 ))
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* forward declarations */
void tlamb(int m, double q, double qsqfm1, double x, int n,
           double *t, double *dt, double *d2t, double *d3t);
int xlamb(int m, double q, double qsqfm1, double tin, double *x, double *xpl);
int vlamb(double gm, double r1, double r2, double th, int nrev, double tdelt,
          double v1[2], double v2[2]);

void tlamb (int m, double q, double qsqfm1, double x, int n, double *t, double *dt, double *d2t, double *d3t) {
    
	double qsq,xsq,u,y,z,qx,a,b,aa,bb,g,f,term,tterm,tqterm,twoi1,qz,qz2,u0i,u1i,u2i,u3i,tq,p;
    double ttmold = 0.0;
    double told = 0.0;
    double tqsum = 0.0;
    int i;
	
    double M = (double) m;
    
    if (fabs(x)==1.0){
        return; /* was: return -1; (bug in void function) */
    }
    qsq = q*q;
	
	xsq = x*x;
	u = (1.0-x)*(1.0+x);

	if (n!=-1) {
		/* needed if series and otherwise useful when z==0 */
		*dt = 0.0; *d2t =0.0; *d3t = 0.0;
    } 
	if ((n==-1)||(m>0)||(x<0.0)||(u>SW)||(-u>SW)){
		/*direct computation (not series) */
		y=sqrt(fabs(u));
		z=sqrt(qsqfm1+qsq*xsq);
		qx = q*x;
		if (qx <= 0.0) {
			a = z-qx;
			b = q*z-x;
		}
		if ((qx < 0.0)&&(n==-1)) {
			aa = qsqfm1/a;
			bb = qsqfm1*(qsq*u-xsq)/b;
		}
		if ((qx==0.0)&&(n==-1)||(qx>0.0)) {
			aa = z+qx;
			bb = q*z+x;
		}
		if (qx > 0.0) {
			a = qsqfm1/aa;
			b = qsqfm1*(qsq*u-xsq)/bb;
		}
		if (n!=-  1) {
			if ((qx*u)>=0.0) {
				g = x*z+q*u;
			} else {
				g = (xsq-qsq*u)/(x*z-q*u);
			}
			f = a*y;
			if (x <=1.0) {
				*t = M*PI+atan2(f,g);
			} else {
				if (f > SW) {
					*t = log(f+g);
				} else {
					double fg1 = f/(g+1.0);
					term = 2.0*fg1;
					double fg1sq = fg1*fg1;
					*t = term;
					twoi1 = 1.0;
					do {
						twoi1 = twoi1+2.0;
						term = term*fg1sq;
						told = *t;
						*t = (*t)+term/twoi1;
					} while ((*t)!=told);
				}
			}
			*t = 2.0*((*t)/y+b)/u;
			if ((n>=1)&&(z!=0.0)) {
				qz = q/z;
				qz2 = qz*qz;
				qz = qz*qz2;
				*dt = (3.0*x*(*t)-4.0*(a+qx*qsqfm1)/z)/u;
				if (n>=2) {
					*d2t = (3.0*(*t)+5.0*x*(*dt)+4.0*qz*qsqfm1)/u;
				}
				if (n>=3) {
					*d3t = (8.0*(*dt)+7.0*x*(*d2t)-12.0*qz*qz2*x*qsqfm1)/u;
				}
			}
		} else {
			*dt = b;
			*d2t = bb;
			*d3t = aa;
		}
        
	} else {
		/*compute by series */
        u0i = 1.0;
        
		if (n >=1) u1i = 1.0;
		if (n >=2) u2i = 1.0;
		if (n >=3) u3i = 1.0;
		term = 4.0;
		tq = q*qsqfm1;
		i = 0;
		if (q < 0.5) tqsum = 1.0-q*qsq;
		if (q >= 0.5) tqsum = (1.0/(1.0+q)+q)*qsqfm1;
		ttmold = term / 3.0;
		*t = ttmold*tqsum; 
		
		do {
			i=i+1;
			p= (double) i;
			u0i = u0i*u;
			if ((n >=1)&&(i>1)) u1i = u1i*u;
			if ((n >=2)&&(i>2)) u2i = u2i*u;
			if ((n >=3)&&(i>3)) u3i = u3i*u;
			term = term*(p-0.5)/p;
			tq = tq*qsq;
			tqsum = tqsum+tq;
			told = *t;
			tterm = term/(2.0*p+3.0);
			tqterm = tterm*tqsum;
			*t = (*t)-u0i*((1.5*p+0.25)*tqterm/(p*p-0.25)-ttmold*tq);
			ttmold = tterm;
			tqterm = tqterm*p;
			if (n >= 1) *dt = (*dt)+tqterm*u1i;
			if (n >= 2) *d2t = (*d2t)+tqterm*u2i*(p-1.0);
			if (n >= 3) *d3t = (*d3t)+tqterm*u3i*(p-1.0)*(p-2.0);
		} while ((i<n)||((*t)!=told));
        
	    if (n >= 3) *d3t = 8.0*x*(1.5*(*d2t)-xsq*(*d3t));
	    if (n >= 2) *d2t = 2.0*(2.0*xsq*(*d2t)-(*dt));	
        if (n >= 1) *dt = -2.0*x*(*dt);	
        *t = (*t)/xsq;
    }

}

int xlamb (int m, double q, double qsqfm1, double tin, double *x, double *xpl) {
    
    double t,t0,dt,d2t,d3t,tdiff,tdiff0,w,xm,tdiffm,d2t2,tmin;
    int count;
    int n = 0;
    
    double thr2 = atan2(qsqfm1,2.0*q)/PI;
    double M = (double) m;
     
    if (m==0) {
        /*single-rev starter from T (at X = 0) & bilinear (usually) */
        n=1;
        tlamb(m,q,qsqfm1,0.0,0,&t0,&dt,&d2t,&d3t);
        tdiff = tin-t0; 
        if(tdiff <= 0.0) {
            *x = t0*tdiff/(-4.0*tin);  /* dt is -4 for x = 0 */
        } else {
            *x = -tdiff/(tdiff+4.0);
            w = (*x)+C0*sqrt(2.0*(1.0-thr2));
            if (w < 0.0) 
                *x = (*x)-sqrt(sqrt(sqrt(sqrt(-w))))*((*x)+sqrt(tdiff/(tdiff+1.5*t0))); 
            
            w = 4.0/(4.0+tdiff);
            *x = (*x)*(1.0+(*x)*(C1*w-C2*(*x)*sqrt(w)));
        }
    } else {
        /* with multirevs, first get t(min) as basis for starter */
        xm = 1.0/(1.5*(M+0.5)*PI);
        if (thr2 < 0.5) {
            xm = sqrt(sqrt(sqrt(2.0*thr2)))*xm;
        } else if (thr2 > 0.5) {
            xm = (2.0-sqrt(sqrt(sqrt(2.0-2.0*thr2))))*xm;
        } 
        
     
        double tmin,xmold;
        for (count=0;count<16;count++) {
            tlamb(m,q,qsqfm1,xm,3,&tmin,&dt,&d2t,&d3t);
            if (d2t==0.0) break;
            xmold = xm;
            xm = xm - dt*d2t/(d2t*d2t-dt*d3t/2.0);
            if (fabs(xmold/xm-1.0)<=TOL) break;
        }
        /* check if used max iterations on loop */
        if (count>=15) return -1; 
        
        tdiffm = tin-tmin;
        if (tdiffm == 0) {
            /* unique solution */
            *x = xm;
            return 1; 
        } else if (tdiffm < 0) {
            /*no solution */
            return 0; 
        } else {
        
          n = 3;
        
          if (d2t==0.0) d2t = 6.0*M*PI;
          *x = sqrt(tdiffm/(d2t/2.0+tdiffm/((1.0-xm)*(1.0-xm))));
          w = xm+(*x);
          w = w*4.0/(4.0+tdiffm)+(1.0-w)*(1.0-w);
          *x = (*x)*(1.0-(1.0+M+C41*(thr2-0.5))/(1.0+C3*M)*(*x)*(C1*w+C2*(*x)*sqrt(w)))+xm;
          d2t2 = d2t/2.0;
        
          if ((*x) >= 1.0) {
              n = 1;
              goto JUMP_POINT2;
          }
        } 
    }

JUMP_POINT1:
    for (count=0;count<3;count++) {
        /* this loop is run twice... not sure why --NS */
        tlamb(m,q,qsqfm1,(*x),2,&t,&dt,&d2t,&d3t);
        t = tin-t;
        if (dt) {
            *x = (*x)+t*dt/(dt*dt+t*d2t/2.0);
        }
    } 
    if (n != 3) return n;
    
    n = 2;
    *xpl = (*x);
    
JUMP_POINT2:
    
    tlamb(m,q,qsqfm1,0.0,0,&t0,&dt,&d2t,&d3t);
    tdiff0 = t0-tmin;
    tdiff = tin-t0;

    if(tdiff <= 0.0) {
        *x = xm-sqrt(tdiffm/(d2t2-tdiffm*(d2t2/tdiff0-1.0/(xm*xm))));
    } else {
        *x = -tdiff/(tdiff+4.0);
        w = (*x)+C0*sqrt(2.0*(1.0-thr2));
        if (w < 0.0) *x = (*x)-sqrt(sqrt(sqrt(sqrt(-w))))*((*x)+sqrt(tdiff/(tdiff+1.5*t0)));
        w = 4.0/(4.0+tdiff);
        *x = (*x)*(1.0+(1.0+M+C42*(thr2-0.5))/(1.0+C3*M)*((*x)*(C1*w-C2*(*x)*sqrt(w))));
        if ((*x) <= -1.0) {
            n = n-1;
            if (n==1) *x = *xpl;
        }
    }
    goto JUMP_POINT1;
}

int vlamb (double gm, double r1, double r2, double th, int nrev, double tdelt, 
           double v1[2], double v2[2]) 
{
    
    /* let nrev be negative for short period, pos for long period */
    /* let th be positive for prograde, neg for retrograde */

    int m;
    m= abs(nrev);
    
    double thr2 = th/2.0;
    double dr = r1-r2;
    double r1r2th = 4.0*r1*r2*pow(sin(thr2),2);
    double csq = dr*dr + r1r2th;
    double c = sqrt(csq);
    double s = (r1+r2+c)/2.0;
    double gms = sqrt(gm*s/2.0);
    double qsqfm1 = c/s;
    double q = sqrt(r1*r2)*cos(thr2)/s;
    
    double rho, sig;  
    if (c) {
        rho = dr/c;
        sig = r1r2th/csq;
    } else {
        rho = 0.0;
        sig = 1.0;
    }
    
    double t = 4.0*gms*tdelt/(s*s);
    
    double x1,x2;
    
    int code = xlamb(m,q,qsqfm1,t,&x1,&x2);
    
    double x,unused,qzminx,qzplx,zplqx;
     
    if (code > 0 ) {
        if (nrev>0) {
            x = x2;
        } else {
            x = x1;
        }
        
        tlamb(m,q,qsqfm1,x,-1,&unused,&qzminx,&qzplx,&zplqx);    
        v1[0] = gms*(qzminx-qzplx*rho)/r1;
        v1[1] = gms*zplqx*sqrt(sig); 
        v2[0] = -gms*(qzminx+qzplx*rho)/r2;
        v2[1] = v1[1]/r2;
        v1[1] = v1[1]/r1;
    }
    
    return code;
}

int vlamb2 (double gm, double r1, double r2, double th, int nrev, double tdelt, 
           double v11[2], double v12[2], double v21[2], double v22[2]) 
{
    
    /* let th be positive for prograde, neg for retrograde */
    
    int m = abs(nrev);
    
    double thr2;
    thr2 = th/2.0;
    
    double dr;
    dr = r1-r2;
    
    double r1r2th;
    r1r2th = 4.0*r1*r2*pow(sin(thr2),2);
    
    double csq;
    csq = dr*dr + r1r2th;
    
    double c;
    c = sqrt(csq);
    
    double s;
    s = (r1+r2+c)/2.0;
    
    double gms;
    gms = sqrt(gm*s/2.0);
    
    double qsqfm1;
    qsqfm1 = c/s;
    
    double q;
    q = sqrt(r1*r2)*cos(thr2)/s; 
    double rho, sig;  
    if (c) {
        rho = dr/c;
        sig = r1r2th/csq;
    } else {
        rho = 0.0;
        sig = 1.0;
    }
    
    double t;
    t = 4.0*gms*tdelt/(s*s);
    
    double x1,x2;
    int code = xlamb(m,q,qsqfm1,t,&x1,&x2);
    
    double unused,qzminx,qzplx,zplqx;
    
    if (code > 0) {
        tlamb(m,q,qsqfm1,x1,-1,&unused,&qzminx,&qzplx,&zplqx);    
        
        v11[0] = gms*(qzminx-qzplx*rho)/r1;
        v11[1] = gms*zplqx*sqrt(sig); 
        v12[0] = -gms*(qzminx+qzplx*rho)/r2;
        v12[1] = v11[1]/r2;
        v11[1] = v11[1]/r1;
        if (code > 1) {
            tlamb(m,q,qsqfm1,x2,-1,&unused,&qzminx,&qzplx,&zplqx);    
            
            v21[0] = gms*(qzminx-qzplx*rho)/r1;
            v21[1] = gms*zplqx*sqrt(sig); 
            v22[0] = -gms*(qzminx+qzplx*rho)/r2;
            v22[1] = v21[0]/r2;
            v21[1] = v21[0]/r1;
        }
    }
    
    return code;
}

int lambert (double gm, double r1[3], double r2[3], int nrev, double tdelt, 
           double v1[3], double v2[3]) 
{
    /* let nrev be negative for short period, pos for long period */
    /* let dt be positive for prograde, neg for retrograde */
     
    double th,rad1,rad2,va1[2],va2[2];
    
    rad1 = sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
    rad2 = sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);
    
    th = acos((r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2])/(rad1*rad2));
    if (th<0) {
        th = th + TWOPI;
    }

    int code = vlamb (gm,rad1,rad2,th,nrev,tdelt,va1,va2);    

    if (code > 0 ) {
        
        /* coordinate bases */
        
        double x1[3],x2[3],y1[3],y2[3],z[3];
        
        
        x1[0] = r1[0]/rad1;
        x1[1] = r1[1]/rad1;
        x1[2] = r1[2]/rad1;
        x2[0] = r2[0]/rad2;
        x2[1] = r2[1]/rad2;
        x2[2] = r2[2]/rad2;
        
        z[0] = x1[1]*x2[2]-x1[2]*x2[1];
        z[1] = x1[2]*x2[0]-x1[0]*x2[2];
        z[2] = x1[0]*x2[1]-x1[1]*x2[0];
        double zm = sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
    
        if (zm < TOL) {
            z[0] = 0.0; z[1] = 0.0; z[2] = 1.0;
        } else {
            z[0] = z[0]/zm;
            z[1] = z[1]/zm;
            z[2] = z[2]/zm; 
        }

        
        y1[0] = z[1]*x1[2]-z[2]*x1[1];
        y1[1] = z[2]*x1[0]-z[0]*x1[2];
        y1[2] = z[0]*x1[1]-z[1]*x1[0];
        y2[0] = z[1]*x2[2]-z[2]*x2[1];
        y2[1] = z[2]*x2[0]-z[0]*x2[2];
        y2[2] = z[0]*x2[1]-z[1]*x2[0];

        
        /* convert velocity vectors */
        v1[0] = va1[0]*x1[0]+va1[1]*y1[0];
        v1[1] = va1[0]*x1[1]+va1[1]*y1[1];
        v1[2] = va1[0]*x1[2]+va1[1]*y1[2];
        v2[0] = va2[0]*x2[0]+va2[1]*y2[0];
        v2[1] = va2[0]*x2[1]+va2[1]*y2[1];
        v2[2] = va2[0]*x2[2]+va2[1]*y2[2]; 
    }
    
    return code;
}


int vlamb_test (double gm, double r1, double r2, double th, int nrev, double tdelt, 
                double v1[2], double v2[2]) 
{
    double vm1sq,vm2sq,en1,en2,am1,am2;
    int ok=1;
     
    int m = abs(nrev);
    double M = (double)m;
    
    vm1sq = v1[0]*v1[0]+v1[1]*v1[1];
    vm2sq = v2[0]*v2[0]+v2[1]*v2[1];
    
    en1 = 0.5*vm1sq-gm/r1;
    en2 = 0.5*vm2sq-gm/r2;
    
    am1 = r1*v1[1];
    am2 = r2*v2[1];
    
    double sma,ecc,ea1,ea2,per,tfp1,tfp2;
    
    sma = -gm/2/en1;
    per = TWOPI*sqrt(sma*sma*sma/gm);
    ecc = sqrt(1-am1*am1/gm/sma);
    ea1 = sign(v1[1])*acos((1-r1/sma)/ecc);
    ea2 = sign(v2[1])*acos((1-r2/sma)/ecc);
    double mm = sqrt(gm/sma/sma/sma);
    tfp1 = (ea1-ecc*sin(ea1))/mm;
    tfp2 = (ea2-ecc*sin(ea2))/mm;
    if (tfp2<tfp1) {
        tfp2 = tfp2+per;
    }
   
    
    ok = ok&&(am1==am2);
   
    if (!ok) {
    printf("Ang. Mom. Error: %6.3g \n",(am1-am2));
    }
    
    ok = ok&&(fabs(en1-en2)<TOL);
    
    if (!ok) {
    printf("Energy Error: %6.3g \n",(en1-en2));
    }
    
    ok = ok&&(fabs(1-((tfp2-tfp1+M*per)/tdelt))<TOL);
    
    if (!ok) {
    printf("TOF Error: %6.3g \n",(1-((tfp2-tfp1+M*per)/tdelt)));
    }
  
    return ok;
}