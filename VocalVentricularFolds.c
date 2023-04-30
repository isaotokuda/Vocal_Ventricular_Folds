/************************************************************************
   Program for simulating the vocal-ventricular fold oscillations
   written by Isao Tokuda (2022)

   Physical Units are in [cm], [msec], [g] 

  Reference: Rintaro Miyazaki, Tomoki Yoshitani, Mayuka Kanaya,
    Shigehiro Miyachi, Akihisa Kaneko, Yuki Kinoshita, 
    Kanta Nakamura, Takeshi Nishimura, Isao T. Tokuda
    ``Ventricular fold oscillations lower the vocal pitch in rhesus macaques,'' 
    Submitted to Journal of Experimental Biology
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h> 

#define QQupp      0.55          //Q-parameter to control frequency of the ventricular folds
#define QQlow      1.00          //Q-parameter to control frequency of the vocal folds

#define DIM       (2*8+1)

#define P_s        0.014         // subglottal pressure

#define M_1       (0.125/4.0)    // 1st mass [g]
#define M_2       (0.025/4.0)    // 2nd mass [g]
#define M_3       (0.125/4.0)    // 3rd mass [g]
#define M_4       (0.025/4.0)    // 4th mass [g]

#define K_1       (0.080)        // Stiffness
#define K_2       (0.008)
#define K_3       (0.080)
#define K_4       (0.008)
#define K_c34     (0.025)
#define K_c12     (0.025)

#define rho        0.001130          // atmospheric density
#define length    (1.4000/2.0)       // length of the vocal folds in cm 
#define dt        (1000.0/(40.0e3))  // integration step size in ms

double XX[DIM];
double kx[DIM];
double rx[DIM];
double xx[DIM];
double P_1,P_2,P_3,P_4;
double a_1;
double a_2;
double a_3;
double a_4;
double a_03;
double a_04;
double a_min;
double a_min1;
double a_min2;

/* Left Vocal Fold Parameters */ 
double m_1l,m_2l,m_3l,m_4l,k_1l,k_2l,k_c12l,k_3l,k_4l,k_c34l;
double c_1l,c_2l,c_3l,c_4l,r_1l,r_2l,r_3l,r_4l;

/* Right Vocal Fold Parameters */ 
double m_1r,m_2r,m_3r,m_4r,k_1r,k_2r,k_c12r,k_3r,k_4r,k_c34r;
double c_1r,c_2r,c_3r,c_4r,r_1r,r_2r,r_3r,r_4r;

int main()
{     
  int    i,j;
  double adduction;
  double U_g;
  double Theta();
  void   Runge_kutta_gill(),ParamQ(),Start();
  FILE   *fp;

  adduction = 0.100/4.0;
  a_03 = a_04 = adduction;

  Start();
  ParamQ();
  for(i=0;i<100000;i++) Runge_kutta_gill();

  fp=fopen("Waveform.txt","w");
  XX[16]=0;
  for(i=0;i<800;i++){
    for(j=0;j<2;j++) Runge_kutta_gill();
    U_g=sqrt(2.0*P_s/rho)*a_min*Theta(a_min);
    if(a_min1<0) a_min1=0;
    if(a_min2<0) a_min2=0;
    fprintf(fp,"%lf %lf %lf %lf\n",XX[16],a_min1,a_min2,U_g);
  }
  fclose(fp);

  return(0);
}     

void ParamQ()
{     
  double dratio;
  dratio = 0.15;
 
  /* left parameters */ 
  m_1l   = M_1/QQlow;
  m_2l   = M_2/QQlow;
  k_1l   = K_1*QQlow;
  k_2l   = K_2*QQlow;
  k_c12l = K_c12*QQlow;
  c_1l   = 3.0*k_1l;
  c_2l   = 3.0*k_2l;

  m_3l   = M_3/QQupp;
  m_4l   = M_4/QQupp;
  k_3l   = K_3*QQupp;
  k_4l   = K_4*QQupp;
  k_c34l = K_c34*QQupp;
  c_3l   = 3.0*k_3l;
  c_4l   = 3.0*k_4l;

  /* right parameters */ 
  m_1r   = m_1l;
  m_2r   = m_2l;
  k_1r   = k_1l;
  k_2r   = k_2l;
  k_c12r = k_c12l;
  c_1r   = c_1l;
  c_2r   = c_2l;

  m_3r   = m_3l;
  m_4r   = m_4l;
  k_3r   = k_3l;
  k_4r   = k_4l;
  k_c34r = k_c34l;
  c_3r   = c_3l;
  c_4r   = c_4l;

  r_1l = (2.0*dratio*sqrt(m_1l*k_1l));
  r_2l = (2.0*dratio*sqrt(m_2l*k_2l));
  r_3l = (2.0*dratio*sqrt(m_3l*k_3l));
  r_4l = (2.0*dratio*sqrt(m_4l*k_4l));
  r_1r = (2.0*dratio*sqrt(m_1r*k_1r));
  r_2r = (2.0*dratio*sqrt(m_2r*k_2r));
  r_3r = (2.0*dratio*sqrt(m_3r*k_3r));
  r_4r = (2.0*dratio*sqrt(m_4r*k_4r));
}

void Start()
{     
  int j;
  XX[0+0]=0.10;
  XX[0+2]=0.10;
  XX[0+4]=0.10;
  XX[0+6]=0.10;
  XX[0+1]=XX[0+3]=XX[0+5]=XX[0+7]=0.0; 
  XX[8+0]=0.10;
  XX[8+2]=0.10;
  XX[8+4]=0.10;
  XX[8+6]=0.10;
  XX[8+1]=XX[8+3]=XX[8+5]=XX[8+7]=0.0; 
  XX[16]=0.0; 
  for(j=0;j<DIM;j++) xx[j]=0.0;
}

void Params()
{
  double Theta();
  double a_01,a_02;
  a_01 = 1.00*0.050/4.0;
  a_02 = 1.00*0.050/4.0;

  a_1 = a_01+length*(XX[0+0]+XX[8+0]);
  a_2 = a_02+length*(XX[0+2]+XX[8+2]);
  a_3 = a_03+length*(XX[0+4]+XX[8+4]);
  a_4 = a_04+length*(XX[0+6]+XX[8+6]);

  if     (a_1<=  a_2) a_min1=a_1;
  else                a_min1=a_2;
  if     (a_3<=  a_4) a_min2=a_3;
  else                a_min2=a_4;

  if     (a_1<=  a_2) a_min=a_1;
  else                a_min=a_2;
  if     (a_3<=a_min) a_min=a_3;
  if     (a_4<=a_min) a_min=a_4;

  if(a_min>0){
    P_1 = P_s*(1.0-Theta(a_min)*pow(a_min/a_1,2.0));
    P_2 = P_s*(1.0-Theta(a_min)*pow(a_min/a_2,2.0));
    P_3 = P_s*(1.0-Theta(a_min)*pow(a_min/a_3,2.0));
    P_4 = P_s*(1.0-Theta(a_min)*pow(a_min/a_4,2.0));
    if(a_min==a_1) P_1=P_2=P_3=P_4=0.0;
    if(a_min==a_2) P_2=P_3=P_4=0.0;
    if(a_min==a_3) P_3=P_4=0.0;
    if(a_min==a_4) P_4=0.0;
  }
  else{
    if     (a_1<0){ P_1=P_2=P_3=P_4=0.0; }
    else if(a_2<0){ 
      P_2=P_3=P_4=0.0;
      if(a_1<0) P_1=0;
      else      P_1=P_s;
    }
    else if(a_3<0){ 
      P_3=P_4=0.0;
      if(a_1<0){ P_1=P_2=0; }
      else     { 
        P_1=P_s;
        if(a_2<0) P_2=0;
        else      P_2=P_s;
      }
    }
    else if(a_4<0){ 
      P_4=0.0;
      if(a_1<0){ P_1=P_2=P_3=0; }
      else     { 
        P_1=P_s;
        if(a_2<0){ P_2=P_3=0; }
        else     { 
          P_2=P_s;
          if(a_3<0) P_3=0;
          else      P_3=P_s;
        }
      }
    }
  }
}

void Runge_kutta_gill()
{
  int    i,j;
  double root;
  double J();
  void   Params();
  root = pow(2.0,-0.50);
  Params();
  for(j=0;j<DIM;j++) kx[j]=dt*J(j);  
  for(j=0;j<DIM;j++){
    rx[j] =0.50*(kx[j]-2.0*xx[j]); 
    XX[j] +=rx[j]; 
    xx[j]+=(3.0*rx[j]-0.50*kx[j]);
  }

  Params();
  for(j=0;j<DIM;j++) kx[j]=dt*J(j);  
  for(j=0;j<DIM;j++){
    rx[j] =(1.0-root)*(kx[j]-xx[j]); 
    XX[j] +=rx[j]; 
    xx[j]+=(3.0*rx[j]-(1.0-root)*kx[j]);
  }

  Params();
  for(j=0;j<DIM;j++) kx[j]=dt*J(j);  
  for(j=0;j<DIM;j++){
    rx[j] =(1.0+root)*(kx[j]-xx[j]); 
    XX[j] +=rx[j]; 
    xx[j]+=(3.0*rx[j]-(1.0+root)*kx[j]);
  }

  Params();
  for(j=0;j<DIM;j++) kx[j]=dt*J(j);  
  for(j=0;j<DIM;j++){
    rx[j] =(kx[j]-2.0*xx[j])/6.0; 
    XX[j] +=rx[j];  
    xx[j]+=(3.0*rx[j]-0.50*kx[j]);
  }
}

double J(j)
{
  double y;
  double d_1,d_2,d_3,d_4;
  d_1 = 0.250/2.0;
  d_2 = 0.050/2.0;
  d_3 = 0.250/2.0;
  d_4 = 0.050/2.0;

  double Theta();
  /* right vocal fold */
  if     (j==(0+0+0)) y = XX[0+0+1];
  else if(j==(0+0+1))
    y = (1.0/m_1r)*(P_1*length*d_1-r_1r*XX[0+0+1]-k_1r*XX[0+0+0]
        -Theta(-a_1)*c_1r*0.5*(a_1/length)
        -k_c12r*(XX[0+0+0]-XX[0+0+2]));
  else if(j==(0+0+2)) y = XX[0+0+3];
  else if(j==(0+0+3)) 
    y = (1.0/m_2r)*(P_2*length*d_2-r_2r*XX[0+0+3]-k_2r*XX[0+0+2]
        -Theta(-a_2)*c_2r*0.5*(a_2/length)
        -k_c12r*(XX[0+0+2]-XX[0+0+0]));

  else if(j==(0+4+0)) y = XX[0+4+1];
  else if(j==(0+4+1))
    y = (1.0/m_3r)*(P_3*length*d_3-r_3r*XX[0+4+1]-k_3r*XX[0+4+0]
        -Theta(-a_3)*c_3r*0.5*(a_3/length)
        -k_c34r*(XX[0+4+0]-XX[0+4+2]));
  else if(j==(0+4+2)) y = XX[0+4+3];
  else if(j==(0+4+3)) 
    y = (1.0/m_4r)*(P_4*length*d_4-r_4r*XX[0+4+3]-k_4r*XX[0+4+2]
        -Theta(-a_4)*c_4r*0.5*(a_4/length)
	-k_c34r*(XX[0+4+2]-XX[0+4+0]));

  /* left vocal fold */
  else if(j==(8+0+0)) y = XX[8+0+1];
  else if(j==(8+0+1))
    y = (1.0/m_1l)*(P_1*length*d_1-r_1l*XX[8+0+1]-k_1l*XX[8+0+0]
        -Theta(-a_1)*c_1l*0.5*(a_1/length)
        -k_c12l*(XX[8+0+0]-XX[8+0+2]));
  else if(j==(8+0+2)) y = XX[8+0+3];
  else if(j==(8+0+3)) 
    y = (1.0/m_2l)*(P_2*length*d_2-r_2l*XX[8+0+3]-k_2l*XX[8+0+2]
        -Theta(-a_2)*c_2l*0.5*(a_2/length)
        -k_c12l*(XX[8+0+2]-XX[8+0+0]));

  else if(j==(8+4+0)) y = XX[8+4+1];
  else if(j==(8+4+1))
    y = (1.0/m_3l)*(P_3*length*d_3-r_3l*XX[8+4+1]-k_3l*XX[8+4+0]
        -Theta(-a_3)*c_3l*0.5*(a_3/length)
        -k_c34l*(XX[8+4+0]-XX[8+4+2]));
  else if(j==(8+4+2)) y = XX[8+4+3];
  else if(j==(8+4+3)) 
    y = (1.0/m_4l)*(P_4*length*d_4-r_4l*XX[8+4+3]-k_4l*XX[8+4+2]
        -Theta(-a_4)*c_4l*0.5*(a_4/length)
	-k_c34l*(XX[8+4+2]-XX[8+4+0]));

  else if(j==16) y = 1.0;
  return (y);
}

double Theta(x)
  double x;
{
  double m;
  m = 1000.0;   
  if (x < 0.0) return(0.0);
  else         return(tanh(m*x));
}



