/*
This code runs the patch simulation whose data is seen in figure 3 of the paper.

The simulation consists of a sphere composed of 128 patches uniformly distributed on
the surface of the sphere. The patches move with the sphere as a solid body.
The sphere bounces on a curved surface with radius of curvature 11cm. 
At each time step we calculate the force acting on the sphere due to all the interactions between
charged patches and image charges. Charge is added to each patch upon contact with the surface. Charging
is limited by the local electric field. This is the sum of fields arising from all the charged patches and all their
corresponding image charges in the surface. 

gcc -o d.ex patch_simulation.c -lm -O3
*/


//Import libraries
# include<GL/glut.h> 
# include<math.h> 
# include<stdio.h>
# include<stdlib.h>
# include<assert.h>


//Define simulation parameters
# define RA0 5.0E-3
# define RHO0 955
# define R 0.11

# define GAMMA 2.80
# define FRQ 40.0
# define OMEGA 2.0*M_PI*(FRQ)

# define F40 40.0
# define A40 5.00e-5
# define O40 2.0*M_PI*(F40)

# define F1 1.0
# define A1 4.0e-3
# define O1 2.0*M_PI*(F1)

# define DT 1.E-4
# define G 9.81

# define K 1000.0
# define KW 1000.0
# define RES 0.80
# define FRIC 0.4
# define VC 0.001  

# define E0 8.85e-12
# define DQ 0.8125e-11
# define QMAX 6.5e-11
# define EMAX 3.0e+6
# define ER 10.0

# define RDT 1.0/DT
# define DT2 DT*DT 

# define GAP (8*2.* RA0)
# define SIZEX 4.E-2
# define SIZEY 6.0E-2
# define SIZEZ 4.E-2
# define ASPEC_YX SIZEY/(SIZEX)
# define ASPEC_ZX SIZEZ/(SIZEX)
# define SCALE  180

# define HALFX -0.5*SIZEX
# define HALFY -0.5*SIZEY
# define HALFZ -0.5*SIZEZ

# define HSIZEX 0.5*SIZEX
# define HSIZEZ 0.5*SIZEZ

# define TMAX 600.0

# define IOUT 1000

# define NP  128
# define NSX 25
# define NSZ 25
# define NS (NSX*NSZ)

# define AS (SIZEX*SIZEZ/NS)
# define RP (SIZEX/(1.77*NSX))

//define variables
double x,y,z;
double xold,yold,zold;
double vx,vy,vz;
double fx,fy,fz;
double ax,ay,az;
double axold,ayold,azold;
double wx,wy,wz;
double tx,ty,tz;
double ys0,vs0;
double m,rm,r;
double ii,ri;
double mu,muw;
double a,t,gam,trestart;
double theta,fh1,fh40;
double top,bottom;

//define arrays
double xp[NP],yp[NP],zp[NP],rp[NP];
double xpstart[NP],ypstart[NP],zpstart[NP];
double qp[NP],cp[NP],d0[NP],dq;
double xs[NS],ys[NS],zs[NS],qs[NS],cs[NS];

//particle type
int ptype[NP];

int ipmax;

//collision detection with surface
int col,colold;

int run;

long iseed=-1234567899; 

float ran2(long * );
void setup(void);
void force(void);
void fwall(void);
void move(void);
void data(void);


int main (int argc, char** argv)
{
 int i,it;
 double factor,qwtot,qbtot;
 char fname[100];
 
 FILE *fp;
 
 run=atoi(argv[1]);
 
 iseed=iseed+17*run;
 
 a=(GAMMA*G)/(OMEGA*OMEGA);
 t=-M_PI/(2.0*OMEGA);
 
 ran2(&iseed);
 setup(); 
 
  while(t<TMAX)
   {
     for(it=0;it<IOUT;it++)
        {
        force();
        fwall();
        move();
        t+=DT;
        }

     printf("run %d time %f\n",run,t);
     
     data();       
    }
}

//Calculate charge and dipole moments and prints to file
void data(void)
 {
  int i,is;
  int psat,ssat;
  double qpx,qpy,qpz,qptot;
  double qwtot,qbtot;
  char fname[200];
  FILE *fp;
 
      qpx=0.0;
      qpy=0.0;
      qpz=0.0;
      qptot=0.0;
      for(i=0;i<NP;i++)
        {
        qpx+=qp[i]*xp[i];
        qpy+=qp[i]*yp[i];
        qpz+=qp[i]*zp[i];
        qptot+=qp[i];
        }
      
         
      sprintf(fname,"eg0_%d_er%.1f_dt%.1e_dq%.1e_r%.2f_f%.2f_1Hz%.1e_40Hz%.2e_g%.2f_%d.dat",
                        NP,ER,DT,DQ,RES,FRIC,A1,A40,GAMMA,run);
      fp=fopen(fname,"a");
      fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e\n",
                 t,qptot,qpx,qpy,qpz,x,y,z,xp[ipmax],yp[ipmax],zp[ipmax]);
      fclose(fp);    
      
 }



//calculates all forces
void force(void)
 {    
  int i,j,is; 
  double d,d2,dx,dy,dz;
  double fc,fcx,fcy,fcz;
  double alpha;
  
     alpha=(ER-1.0)/(ER+1.0);
     
     fx=fh40; 
     fz=fh40+fh1;  
     fy=0.0;   
     
     tx=0.0;
     ty=0.0;
     tz=0.0;
     
     if(col == 0)
     {
     
     for(i=0;i<NP;i++)
      for(j=0;j<NP;j++)
      {     
      dx=x+xp[i]-(x+xp[j]); 
      dy=y+yp[i]-(-y-yp[j]+2.0*ys0);
      dz=z+zp[i]-(z+zp[j]);

      d2=dx*dx+dy*dy+dz*dz;
      d=sqrt(d2);
     
      fc=((qp[j]/AS)/(2.0*E0))*(1.0-d/sqrt(RP*RP+d2))*(qp[i])*alpha;
     
      fcx=-fc*dx/d;
      fcy=-fc*dy/d;
      fcz=-fc*dz/d;

      fx+=fcx;
      fy+=fcy;
      fz+=fcz;
      tx+=-(fcy*zp[i]-fcz*yp[i]);
      ty+=-(fcz*xp[i]-fcx*zp[i]);
      tz+=-(fcx*yp[i]-fcy*xp[i]);
      }   
     

 
     for(i=0;i<NP;i++)
      for(is=0;is<NS;is++)
      {     
      dx=x+xp[i]-xs[is];
      dy=y+yp[i]-ys[is];
      dz=z+zp[i]-zs[is];

      d2=dx*dx+dy*dy+dz*dz;
      d=sqrt(d2);
     
      fc=((qs[is]/AS)/(2.0*E0))*(1.0-d/sqrt(RP*RP+d2))*(qp[i])*(1.0-alpha);
     
      fcx=-fc*dx/d;
      fcy=-fc*dy/d;
      fcz=-fc*dz/d;

      fx+=fcx;
      fy+=fcy;
      fz+=fcz;
      tx+=-(fcy*zp[i]-fcz*yp[i]);
      ty+=-(fcz*xp[i]-fcx*zp[i]);
      tz+=-(fcx*yp[i]-fcy*xp[i]);
      }
      
     }
     
 }

// Collisions with the surface
void fwall(void)
{
 int i,imin;
 int is,isx,isz;
 double d,d2,d2min;
 double dx,dy,dz;
 double nx,ny,nz;
 double vrx,vry,vrz,vn;
 double vrtotx,vrtoty,vrtotz;
 double vrtanx,vrtany,vrtanz;
 double f1,f2,fn;
 double vtx,vtz,vtm;
 double ntx,nty,ntz;
 double ftx,fty,ftz;
 double qsum,qmax;
 double xg,yg,zg;
 double ec,etot;
 double epx,epy,epz;
 double esx,esy,esz;
 double alpha;

 
 FILE *fp;
 
      d=x-r+HSIZEX;

      if(d<0) 
      {
       vrx=vx;     
       fn=-K*d-muw*vrx;
       fx+=fn;
      }
      
      d=HSIZEX-x-r;

      if(d<0) 
      {
       vrx=vx;     
       fn=K*d-muw*vrx;
       fx+=fn;
      }
  
      d=z-r+HSIZEZ;

      if(d<0.0) 
      {
       vrz=vz;     
       fn=-K*d-muw*vrz;
       fz+=fn;
      }
      
      d=HSIZEZ-z-r;

      if(d<0.0) 
      {
       vrz=vz;     
       fn=K*d-muw*vrz;
       fz+=fn;
      }
      
      d=GAP+ys0-y-r;

      if(d<0.0) 
      {      
       vry=vy-vs0;   
       fn=K*d-muw*vry; 
       fy+=fn;
      }
  
  
      dx=-x;
      dy=R+ys0-y;
      dz=-z;
      
      d=sqrt(dx*dx+dy*dy+dz*dz);

      if(d+r>R) 
      {
       col=1;
        
       nx=dx/d;
       ny=dy/d;
       nz=dz/d;
       
       f1=K*(d+r-R);
       
       vrx=vx;
       vry=vy-vs0;
       vrz=vz;
       
       vn=vrx*nx+vry*ny+vrz*nz;
       
       f2=-muw*vn;
       
       fn=f1+f2;
       
       fx+=fn*nx;
       fy+=fn*ny;
       fz+=fn*nz;
       
       
       dx=-r*nx;
       dy=-r*ny;
       dz=-r*nz;
	          
       vrtotx=vrx+wy*dz-wz*dy;
       vrtoty=vry+wz*dx-wx*dz;
       vrtotz=vrz+wx*dy-wy*dx;
	     
       vrtanx=vrtotx-vn*nx;
       vrtany=vrtoty-vn*ny;
       vrtanz=vrtotz-vn*nz;
       
       vtm=sqrt(vrtanx*vrtanx+vrtany*vrtany+vrtanz*vrtanz);
      
       if(vtm > VC)
       {
       ntx=vrtanx/vtm;
       nty=vrtany/vtm;
       ntz=vrtanz/vtm;
       
       ftx=-(FRIC*fn)*ntx;
       fty=-(FRIC*fn)*nty;
       ftz=-(FRIC*fn)*ntz;
            	     
       fx+=ftx;
       fy+=fty;
       fz+=ftz;
	     
       tx-=fty*dz-ftz*dy;
       ty-=ftz*dx-ftx*dz;
       tz-=ftx*dy-fty*dx;      
       }
       else
       {  
       ntx=vrtanx/VC;
       nty=vrtany/VC;
       ntz=vrtanz/VC;
       
       ftx=-(FRIC*fn)*ntx;
       fty=-(FRIC*fn)*nty;
       ftz=-(FRIC*fn)*ntz;
            	     
       fx+=ftx;
       fy+=fty;
       fz+=ftz;
	     
       tx-=fty*dz-ftz*dy;
       ty-=ftz*dx-ftx*dz;
       tz-=ftx*dy-fty*dx;      
       }

      }
      else if(d+r<R)
      {
      col=0;
      }

    if(col==0 && colold==1)
     {

       d2min=1.e10;
     
       for(i=0;i<NP;i++)
         {
         d2=yp[i];
         if(d2<d2min)
          {
          d2min=d2;
          imin=i;
          }
         }
       
        isx=((x+HSIZEX)/SIZEX)*NSX+0.0;
        isz=((z+HSIZEZ)/SIZEZ)*NSZ+0.0;
        is=isx+NSX*isz;
      
        cp[imin]=cp[imin]+1.0;
        cs[is]=cs[is]+1.0;
        
        alpha=(ER-1.0)/(ER+1.0);
        
        xg=(x+xp[imin]);   
        yg=(y+yp[imin]);      
        zg=(z+zp[imin]);       
      
        epx=0.0;
        epy=0.0;
        epz=0.0;
        for(i=0;i<NP;i++)
          {     
          dx=x+xp[i]-xg; 
          dy=y+yp[i]-yg; 
          dz=z+zp[i]-zg; 

          d2=dx*dx+dy*dy+dz*dz;
          d=sqrt(d2);
     
          ec=((qp[i]/AS)/(2.0*E0))*(1.0-d/sqrt(RP*RP+d2));

          if(d == 0.0)
           {
           epy+=-ec*(1.0+alpha);
           }
          else
           {
           epx+=-ec*(1.0-alpha)*dx/d;
           epy+=-ec*(1.0+alpha)*dy/d;
           epz+=-ec*(1.0-alpha)*dz/d;
           } 
          }        
                        
          esy=-((qs[is]/AS)/(2.0*E0))*(1.0-alpha);
          
          etot=sqrt(epx*epx+(epy+esy)*(epy+esy)+epz*epz);        
        
          if(etot<EMAX)
           {
           qp[imin]=qp[imin]+DQ;     
           }       


          if(etot>EMAX)
           {
           qp[imin]=qp[imin]-DQ;  
           if(qp[imin]<0.0) qp[imin]=0.0;                   
           } 
      }
     
    colold=col;

}

//update timestep
void move(void)
 {
 int i,is;
 double xnew,ynew,znew;
 double axnew,aynew,aznew;
 double thetax,thetay,thetaz;
 double x1,y1,z1;
 double x2,y2,z2;
    
       xnew=2.0*x-xold+DT2*fx*rm;
       vx=RDT*(xnew-x);
       xold=x;
       x=xnew;

       ynew=2.0*y-yold+DT2*(fy*rm-G);
       vy=RDT*(ynew-y);
       yold=y;
       y=ynew;

       znew=2.0*z-zold+DT2*fz*rm;
       vz=RDT*(znew-z);
       zold=z;
       z=znew;
       
       axnew=2.0*ax-axold+DT2*tx*ri;
       wx=RDT*(axnew-ax);
       axold=ax;
       ax=axnew;

       aynew=2.0*ay-ayold+DT2*ty*ri;
       wy=RDT*(aynew-ay);
       ayold=ay;
       ay=aynew;

       aznew=2.0*az-azold+DT2*tz*ri;
       wz=RDT*(aznew-az);
       azold=az;
       az=aznew;
       
       thetax=ax-axold;
       thetay=ay-ayold;
       thetaz=az-azold;
       
       for(i=0;i<NP;i++)
         {
         x1=xp[i];
         y1=yp[i]*cos(thetax)-zp[i]*sin(thetax);
         z1=yp[i]*sin(thetax)+zp[i]*cos(thetax);
         
         x2=x1*cos(thetay)+z1*sin(thetay);
         y2=y1;
         z2=-x1*sin(thetay)+z1*cos(thetay);
         
         xp[i]=x2*cos(thetaz)-y2*sin(thetaz);
         yp[i]=x2*sin(thetaz)+y2*cos(thetaz);
         zp[i]=z2;
         }
            
       theta = OMEGA*(t+DT);
       ys0=a*sin(theta) + a;
       vs0=a*OMEGA*cos(theta);
       
       fh1=m*A1*O1*O1*sin(O1*(t+DT));
       fh40=m*A40*O40*O40*sin(O40*(t+DT));        
                 
       for(is=0;is<NS;is++)
         {
         ys[is]=ys0+R-sqrt(R*R-xs[is]*xs[is]-zs[is]*zs[is]);
         }        
       
 }

//setup simulation
void setup(void)
  {
  int i,isx,isz,is,fret;
  double th,ph,xpmax,factor;
  char fname[100];
  FILE *fp;
  
      r=RA0;
      m=(4.0/3.0)*M_PI*r*r*r*RHO0;
      rm=1.0/m;   
      
      ii=(2.0/5.0)*m*r*r;
      ri=1.0/ii;
                   
      factor=sqrt(1.0+log(RES)*log(RES)/(M_PI*M_PI));
      mu=-(2.0/M_PI)*log(RES)*sqrt(K*m*m/(m+m))/factor;
      muw=-(2.0/M_PI)*log(RES)*sqrt(KW*m)/factor;
      
      x=0.0*RA0;   
      y=1.0*RA0;   
      z=0.0*SIZEZ;
     
      vx=0.0; 
      vy=0.001*ran2(&iseed); 
      vz=0.0; 

      xold=x-vx*DT;
      yold=y-vy*DT;
      zold=z-vz*DT;
      
      ax=0.0;
      ay=0.0; 
      az=0.0;
     
      wx=0.0*(ran2(&iseed)-0.5);
      wy=0.0*(ran2(&iseed)-0.5);
      wz=0.0*(ran2(&iseed)-0.5);

      axold=ax-wx*DT;
      ayold=ay-wy*DT;
      azold=az-wz*DT;
      
           
      sprintf(fname,"config.dat");
      fp=fopen(fname,"r");

      for(i=0;i<NP;i++) fret=fscanf(fp,"%lf %lf %lf %lf\n",&xpstart[i],&ypstart[i],&zpstart[i],&rp[i]);
      
      for(i=0;i<NP;i++)
        {
        xp[i]=xpstart[i];
        yp[i]=ypstart[i];
        zp[i]=zpstart[i];
        }

      dq=DQ;
      xpmax=-10;
      for(i=0;i<NP;i++) 
        {
        qp[i]=0.0; 
        cp[i]=0.0;
        if(xp[i] > 0)
          {
          ptype[i]=0;
          if(xp[i]>xpmax)
            {
            xpmax=xp[i];
            ipmax=i;
            }
          }
        else
          {
          ptype[i]=1;
          }
        }
      
      for(i=0;i<NP;i++) d0[i]=sqrt((xp[i]-xp[0])*(xp[i]-xp[0])
                                  +(yp[i]-yp[0])*(yp[i]-yp[0])
                                  +(zp[i]-zp[0])*(zp[i]-zp[0]));
      
      
      theta = OMEGA*t;
      ys0=a*sin(theta) + a;
      vs0=a*OMEGA*cos(theta);
      
      fh1=m*A1*O1*O1*sin(O1*t);
      fh40=m*A40*O40*O40*sin(O40*t);        
         
      for(isz=0;isz<NSZ;isz++)
       for(isx=0;isx<NSX;isx++)
         {
         is=isx+NSX*isz;
         xs[is]=(isx+0.5)*(SIZEX/NSX)-HSIZEX;
         zs[is]=(isz+0.5)*(SIZEZ/NSZ)-HSIZEZ;
         ys[is]=ys0+R-sqrt(R*R-xs[is]*xs[is]-zs[is]*zs[is]);
         qs[is]=0.0;
         cs[is]=0.0;
         }
      
      qs[0]=0.0;
      
      colold=0;
      
  }


//random number generation
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
