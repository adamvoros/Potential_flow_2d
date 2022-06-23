#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include "nrutil.h"
#include <stdio.h>
FILE  *  initgr(int i,int j,int k,int l,char * s);
void closegr(FILE  *g1);
void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	    int itmax, int *iter, double *err);
double pythag(double a, double b);

using namespace std;


const double v=1.;
unsigned long NN,M;

double Lxmin,Lxmax,Lymin,Lymax;
double *node_x,*node_y,*b,Lx,Ly,*psi;
int *A,*B,*C,meg[4],**ipszi;
unsigned long *ija;
double *sa;


bool belulvan(double x,double y,int n,double *z){
  double ax=node_x[A[n]],ay=node_y[A[n]],bx=node_x[B[n]],by=node_y[B[n]],cx=node_x[C[n]],cy=node_y[C[n]],apsi=psi[A[n]],bpsi=psi[B[n]],cpsi=psi[C[n]],u,v,t;
  t=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
  u=((bpsi-apsi)*(cy-ay)-(by-ay)*(cpsi-apsi))/t;
  v=((bx-ax)*(cpsi-apsi)-(bpsi-apsi)*(cx-ax))/t;
  *z=apsi+u*(x-ax)+v*(y-ay);
  bool b1=((x-ax)*(by-ay)>=(y-ay)*(bx-ax))&&((x-bx)*(cy-by)>=(y-by)*(cx-bx))&&((x-cx)*(ay-cy)>=(y-cy)*(ax-cx));
  bool b2=((x-ax)*(by-ay)<=(y-ay)*(bx-ax))&&((x-bx)*(cy-by)<=(y-by)*(cx-bx))&&((x-cx)*(ay-cy)<=(y-cy)*(ax-cx));
  if(t>0){
    return b2;
  }else{
    return b1;
  }
  
}

bool belulvan_phi(double x,double y,int n){
  double ax=node_x[A[n]],ay=node_y[A[n]],bx=node_x[B[n]],by=node_y[B[n]],cx=node_x[C[n]],cy=node_y[C[n]],t;
  t=(bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
  bool b1=((x-ax)*(by-ay)>=(y-ay)*(bx-ax))&&((x-bx)*(cy-by)>=(y-by)*(cx-bx))&&((x-cx)*(ay-cy)>=(y-cy)*(ax-cx));
  bool b2=((x-ax)*(by-ay)<=(y-ay)*(bx-ax))&&((x-bx)*(cy-by)<=(y-by)*(cx-bx))&&((x-cx)*(ay-cy)<=(y-cy)*(ax-cx));
  if(t>0){
    return b2;
  }else{
    return b1;
  }
  
}

void interpol_phi(double *phi,double **xpsi,double h){
  for(int n=0;n<=NN;n++){
    //ciklus a háromszögekre
    double xmin,xmax,ymin,ymax;
    xmin=node_x[A[n]];
    xmax=node_x[A[n]];
    ymin=node_y[A[n]];
    ymax=node_y[A[n]];
    if(xmin>node_x[B[n]]) xmin=node_x[B[n]];
    if(xmin>node_x[C[n]]) xmin=node_x[C[n]];
    if(xmax<node_x[B[n]]) xmax=node_x[B[n]];
    if(xmax<node_x[C[n]]) xmax=node_x[C[n]];
    if(ymin>node_y[B[n]]) ymin=node_y[B[n]];
    if(ymin>node_y[C[n]]) ymin=node_y[C[n]];
    if(ymax<node_y[B[n]]) ymax=node_y[B[n]];
    if(ymax<node_y[C[n]]) ymax=node_y[C[n]];
    int imin=(int)((xmin-Lxmin)/h);
    int imax=(int)((xmax-Lxmin)/h);
    int jmin=(int)((ymin-Lymin)/h);
    int jmax=(int)((ymax-Lymin)/h);
    for(int i=imin;i<=imax;i++){
      double x=i*h+Lxmin;
      for(int j=jmin;j<=jmax;j++){
	double y=j*h+Lymin;
	if(belulvan_phi(x,y,n)) xpsi[i][j]=phi[n];
      }
    }
  }
}


void interpol(double **xpsi,double h){
  for(int n=0;n<=NN;n++){
    //ciklus a háromszögekre
    double xmin,xmax,ymin,ymax;
    xmin=node_x[A[n]];
    xmax=node_x[A[n]];
    ymin=node_y[A[n]];
    ymax=node_y[A[n]];
    if(xmin>node_x[B[n]]) xmin=node_x[B[n]];
    if(xmin>node_x[C[n]]) xmin=node_x[C[n]];
    if(xmax<node_x[B[n]]) xmax=node_x[B[n]];
    if(xmax<node_x[C[n]]) xmax=node_x[C[n]];
    if(ymin>node_y[B[n]]) ymin=node_y[B[n]];
    if(ymin>node_y[C[n]]) ymin=node_y[C[n]];
    if(ymax<node_y[B[n]]) ymax=node_y[B[n]];
    if(ymax<node_y[C[n]]) ymax=node_y[C[n]];
    int imin=(int)((xmin-Lxmin)/h);
    int imax=(int)((xmax-Lxmin)/h);
    int jmin=(int)((ymin-Lymin)/h);
    int jmax=(int)((ymax-Lymin)/h);
    for(int i=imin;i<=imax;i++){
      double x=i*h+Lxmin;
      for(int j=jmin;j<=jmax;j++){
	double y=j*h+Lymin,z;
	if(belulvan(x,y,n,&z)) xpsi[i][j]=z;
      }
    }
  }
}


int teszt0(int i){
  //Az i-edik (már bejárt) háromszög olyan élszomszédját keressük, 
  //melynek a harmadik csúcsában psi még ismeretlen
  //ha nincs ilyen, akkor -1 a visszatérési érték, 
  //egyébként a háromszög sorszáma
  int j,k,n;
  for(j=1;j<=3;j++){
    k=ipszi[i][j];//a szomszédos háromszög sorszáma
    n=-1;
    if(k==-1) continue; //nincs szomszédos háromszög
    if(ipszi[k][0]<0) continue;//már bejártuk a szomszédos háromszöget
    n=k;
    for(int m=1;m<=3;m++){
      if(ipszi[k][m]==i||ipszi[k][m]<0){
	continue;
      }else{
	int s=ipszi[k][m];//cout<<j<<"  "<<k<<"  "<<s<<"  eddig jó\n";
	if(ipszi[s][0]<0){
//két bejárt szomszéd háromszög van (i és s)
	  ipszi[k][0]=-1;
//ezért a k-adik háromszög minden csúcsában ismert a psi
	  n=-1;
	  break;
	}
      }
    }
    if(n>-1) return n; 
  }
  return n;
}

int teszt(int i){
  //Az i-edik (már bejárt) háromszög olyan élszomszédját keressük, 
  //melynek a harmadik csúcsában psi még ismeretlen
  //ha nincs ilyen, akkor -1 a visszatérési érték, 
  //egyébként a háromszög sorszáma
  int j,k;
  for(j=1;j<=3;j++){
    k=ipszi[i][j];//a szomszédos háromszög sorszáma
    if(k==-1) continue; //nincs szomszédos háromszög
    if(psi[A[k]]>-1.e9&&psi[B[k]]>-1.e9&&psi[C[k]]>-1.e9) continue;
    //már bejártuk a szomszédos háromszöget
    return k; 
  }
  return -1;
}
  


double hossz(int i, int j){
  //az i-edik és a j-edik pont távolsága
  return pythag(node_x[i]-node_x[j],node_y[i]-node_y[j]);
}

void center(int n,double *x,double *y){
  //az n-edik háromszög köré írható kör középpontjának koordinátáit adja vissza
  double x1,x2,x3,y1,y2,y3;
  x1=node_x[A[n]];
  x2=node_x[B[n]];
  x3=node_x[C[n]];
  y1=node_y[A[n]];
  y2=node_y[B[n]];
  y3=node_y[C[n]];
  *x=.5*(x2+x1)+.5*((x3-x1)*(x3-x2)+(y3-y1)*(y3-y2))/((x3-x1)*(y2-y1)-(y3-y1)*(x2-x1))*(y2-y1); 
  *y=.5*(y2+y1)-.5*((x3-x1)*(x3-x2)+(y3-y1)*(y3-y2))/((x3-x1)*(y2-y1)-(y3-y1)*(x2-x1))*(x2-x1);

}

void ritka_matrix_indexelese(){
    //az ija vektor megadasa
  int nsz[3],xa,xb,xc,ia,ib,n;
    unsigned long ii=NN+1;
  ifstream be("kep.1.neigh");
  be>>n;be>>n;
  be>>n;
  be>>n;if(n>-1) meg[1]=n-1;
  be>>n;if(n>-1) meg[2]=n-1;
  be>>n;if(n>-1) meg[3]=n-1;
  ipszi[0][0]=0;
  ipszi[0][1]=meg[1];
  ipszi[0][2]=meg[2];
  ipszi[0][3]=meg[3];
    for(int n=1;n<=NN;n++){
      //szomszedok(n,nsz);
	ija[n]=0;
	 be>>nsz[0];be>>nsz[1];be>>nsz[2];be>>nsz[3];
	 ipszi[n][0]=n;
	 ipszi[n][1]=nsz[1];if(ipszi[n][1]>=0) ipszi[n][1]--;
	 ipszi[n][2]=nsz[2];if(ipszi[n][2]>=0) ipszi[n][2]--;
	 ipszi[n][3]=nsz[3];if(ipszi[n][3]>=0) ipszi[n][3]--;
	 for(int i=1;i<=3;i++){nsz[i]--;if(nsz[i]==0) nsz[i]--;}
	 xa=nsz[1];ia=1;
	 if(xa>nsz[2]){ xa=nsz[2];ia=2;}
	 if(xa>nsz[3]){ xa=nsz[3];ia=3;}
	 ib=(ia+1)%3+1;
	 xb=nsz[ib];
	 if(xb>nsz[6-ia-ib]){ ib=6-ia-ib;xb=nsz[ib];}
	 xc=nsz[6-ia-ib];
	 nsz[1]=xa;
	 nsz[2]=xb;
	 nsz[3]=xc;
	for(int j=1;j<=3;j++){
	    if(nsz[j]>-1){
		ii++;
		if(ija[n]==0) ija[n]=ii;
		ija[ii]=nsz[j];
	    }
	}
    }
    ija[NN+1]=ii+1;
    //cout<<ija[1]<<"  "<<NN+2<<endl;
    be.close();
}





void ritka_matrix_elemei(){
//az sa vektor megadasa

  int n;
  double xc,yc,xp,yp,s,a;

  for(n=1;n<ija[NN+1];n++) sa[n]=0.;

  for(n=1;n<=NN;n++) b[n]=0.;
//A b vektor a lineáris egyenletrendszer jobboldalát tartalmazza.

  for(n=1;n<=NN;n++){
    int i1,i2,j1,j2;
    i1=ija[n];
    i2=ija[n+1];
    center(n,&xp,&yp);
    for(int i=i1;i<i2;i++){
      center(ija[i],&xc,&yc);
      s=pythag(xc-xp,yc-yp);
      j1=A[ija[i]];
      if(j1!=A[n]&&j1!=B[n]&&j1!=C[n]){
	j1=B[ija[i]];j2=C[ija[i]];
      }else{
	j2=B[ija[i]];
	if(j2!=A[n]&&j2!=B[n]&&j2!=C[n]) j2=C[ija[i]];
      }
      a=hossz(j1,j2);
      sa[n]-=a/s;
      sa[i]=a/s;
    }
//külön kezelendők a 0. háromszög szomszédai. 
//A diagonális elemhez kapunk járulékot.   
    for(int i=1;i<=3;i++){
      if(meg[i]==n){
	center(0,&xc,&yc);
	s=pythag(xc-xp,yc-yp);
	j1=A[0];
      if(j1!=A[n]&&j1!=B[n]&&j1!=C[n]){
	j1=B[0];j2=C[0];
      }else{
	j2=B[0];
	if(j2!=A[n]&&j2!=B[n]&&j2!=C[n]) j2=C[0];
      }
      a=hossz(j1,j2);
      sa[n]-=a/s;
      }
    }

    if(i2-i1<3){
//ha a háromszög a határon/sarokban van, csak az az él ad járulékot, 
//amelyik a bal szélre (x=2) vagy a jobb szélre (x=22) illeszkedik.
      j1=-1;j2=-1;
      if(fabs(node_x[A[n]]-Lxmin)<1.e-5||fabs(node_x[A[n]]-Lxmax)<1.e-5) j1=A[n];
      if(fabs(node_x[B[n]]-Lxmin)<1.e-5||fabs(node_x[B[n]]-Lxmax)<1.e-5){ 
	if(j1<0){j1=B[n];}else{j2=B[n];}
      }
      if(fabs(node_x[C[n]]-Lxmin)<1.e-5||fabs(node_x[C[n]]-Lxmax)<1.e-5){ 
	if(j1<0){j1=C[n];}else{if(j2<0) j2=C[n];}
      }
      if(j1>0&&j2>0){ 
	b[n]=v*hossz(j1,j2);
	if(fabs(node_x[j1]-Lxmax)<1.e-5) b[n]*=-1.;
	//a jobbszél járuléka ellentétes előjelű
      }
    }
  }


}





int main(){
  int n;
 
  ifstream be("kep.1.ele");
  be>>NN;
  NN-=1;
  be>>n;be>>n;
  A=ivector(0,NN);
  B=ivector(0,NN);
  C=ivector(0,NN);
  for(int i=0;i<=NN;i++) be>>n>>A[i]>>B[i]>>C[i];
  be.close();

  ipszi=imatrix(0,NN,0,3);

  

  be.open("kep.1.node");
  be>>M;
  be>>n;be>>n;be>>n;
  node_x=dvector(1,M);
  node_y=dvector(1,M);
  Lxmin=1.e10;
  Lxmax=-1.e10;
  Lymin=1.e10;
  Lymax=-1.e10;
  for(int i=1;i<=M;i++){
    be>>n>>node_x[i]>>node_y[i]>>n;
    if(Lxmin>node_x[i]) Lxmin=node_x[i];
    if(Lxmax<node_x[i]) Lxmax=node_x[i];
    if(Lymin>node_y[i]) Lymin=node_y[i];
    if(Lymax<node_y[i]) Lymax=node_y[i];
  }
  be.close();
  Lx=Lxmax-Lxmin;
  Ly=Lymax-Lymin;
  
  double *b0,*x;
  b=dvector(1,NN);
  b0=dvector(1,NN);

  x=dvector(1,NN);
  int iter;
  double err; 
  ija=lvector(1,4*NN); 
  sa=dvector(1,4*NN);
  ritka_matrix_indexelese();
  ritka_matrix_elemei();

  for(int j=1;j<=NN;j++) b0[j]=b[j];

  linbcg(NN, b, x, 1, 1.e-6, 10000, &iter, &err);

  //Az áramlási függvény meghatározása következik
  psi=dvector(1,M);int k;
  //v_x=dpsi/dy;v_y=-dpsi/dx;
  for(int i=1;i<=M;i++) psi[i]=-1.e10;
  int *ell;ell=ivector(1,M);
  for(int i=1;i<=M;i++) ell[i]=0;
  ell[1]=A[0];
  ell[2]=B[0];
  ell[3]=C[0];

  int jj=3,j=0;srand(1234);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //Az áramlási függvény meghatározása a 0. háromszög csúcsaiban
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  //Az A csúcsban legyen psi=0
  psi[A[0]]=0.;
  //B csúcs
  int volt=0;
  double x1,x2,y1,y2,z;
  for(int i=1;i<=3;i++){
    k=ipszi[j][i];
    if(k>-1){
      if(A[k]!=A[j]&&B[k]!=A[j]&&C[k]!=A[j]) continue;
      if(A[k]!=B[j]&&B[k]!=B[j]&&C[k]!=B[j]) continue;
      center(j,&x1,&y1);
      center(k,&x2,&y2);
      z=pythag(x1-x2,y1-y2);
      psi[B[j]]=psi[A[j]]+hossz(A[j],B[j])/z*x[k];
      volt++;
    } 	
  }
  if(volt==0){
    //ha nincs az AB éllel szomszédos háromszög
    if(fabs(node_x[A[j]]-2.)<1.e-5&&fabs(node_x[B[j]]-2.)<1.e-5){
      if(node_y[B[j]]>node_y[A[j]]){
	psi[B[j]]=psi[A[j]]-hossz(A[j],B[j])*v;
      }else{ 
	psi[B[j]]=psi[A[j]]+hossz(A[j],B[j])*v;
      }
    }
    
    if(fabs(node_x[A[j]]-22.)<1.e-5&&fabs(node_x[B[j]]-22.)<1.e-5){
      if(node_y[B[j]]>node_y[A[j]]){
	psi[B[j]]=psi[A[j]]+hossz(A[j],B[j])*v;
      }else{ 
	psi[B[j]]=psi[A[j]]-hossz(A[j],B[j])*v;
      }
    }
  }
    
    //C csúcs
    volt=0;
    for(int i=1;i<=3;i++){
      k=ipszi[j][i];
      if(k>-1){
	if(A[k]!=A[j]&&B[k]!=A[j]&&C[k]!=A[j]) continue;
	if(A[k]!=C[j]&&B[k]!=C[j]&&C[k]!=C[j]) continue;
	center(j,&x1,&y1);
	center(k,&x2,&y2);
	z=pythag(x1-x2,y1-y2);
	psi[C[j]]=psi[A[j]]+hossz(A[j],C[j])/z*x[k];
	volt++;
      } 	
    }
    if(volt==0){
      //ha nincs az AC éllel szomszédos háromszög
      if(fabs(node_x[A[j]]-2.)<1.e-5&&fabs(node_x[C[j]]-2.)<1.e-5){
	if(node_y[C[j]]>node_y[A[j]]){
	  psi[C[j]]=psi[A[j]]-hossz(A[j],C[j])*v;
	}else{ 
	  psi[C[j]]=psi[A[j]]+hossz(A[j],C[j])*v;
	}
      }
      
      if(fabs(node_x[A[j]]-22.)<1.e-5&&fabs(node_x[C[j]]-22.)<1.e-5){
	if(node_y[C[j]]>node_y[A[j]]){
	  psi[C[j]]=psi[A[j]]+hossz(A[j],C[j])*v;
	}else{ 
	  psi[C[j]]=psi[A[j]]-hossz(A[j],C[j])*v;
	}
      }
    }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //Az áramlási függvény meghatározása a többi csúcsban 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while(jj<M){
    //olyan háromszöget keresünk, aminek két csúcsában psi már ismert,
    //a harmadikban még nem.
    //a j-edik háromszöget már bejártuk (kezdetben j=0)
    
    k=teszt(j);
    while(k<0){
      int m=(int)(3.*(double)rand()/(double)RAND_MAX)+1;
      if(m==4) m=3;
      k=ipszi[j][m];
      while(k<0){
	m%=3;m++;
	k=ipszi[j][m];
      }
      j=k;
      k=teszt(j);
    }
      
    jj++;
    //A k-adik háromszögnek a j-edikhez nem tartozó csúcsában keressük psi-t

    int m,j1,j2;
//m lesz a j-edik háromszöghöz nem tartozó csúcs
//j1 és j2 a j-edik háromszöghöz tartozó két csúcs
 
    m=A[k];j1=B[k];j2=C[k];
    if(m==A[j]||m==B[j]||m==C[j]){ m=B[k];j1=C[k];j2=A[k];}
    if(m==A[j]||m==B[j]||m==C[j]){ m=C[k];j1=A[k];j2=B[k];}

    //Meghatározzuk psi növekményét, ami a j1-ből m-be jutáskor fellép 
    //A k-adik háromszög szomszédai ipszi[k][1], ipszi[k][2], ipszi[k][3]
    //Ezek közül az egyik j-vel egyenlő. 
    //A másik kettő lehet mindkettő -1, és lehet egyik vagy mindkettő >=0  
    n=ipszi[k][1];
    if(n==j||n==-1) n=ipszi[k][2];
    if(n==j||n==-1) n=ipszi[k][3];

    if(n==j||n==-1){
      psi[m]=psi[j1];
      //két oldalon nincs szomszéd, így az m-j1 él a peremen van
      //csak akkor nem nulla a járulék, ha a balszélen vagy a jobbszélen van
      if(fabs(node_x[m]-Lxmin)<1.e-5&&fabs(node_x[j1]-Lxmin)<1.e-5){
	
	if(node_y[m]>node_y[j1]){
	  psi[m]=psi[j1]-hossz(m,j1)*v;
	}else{ 
	  psi[m]=psi[j1]+hossz(m,j1)*v;
	}
      }

      if(fabs(node_x[m]-Lxmax)<1.e-5&&fabs(node_x[j1]-Lxmax)<1.e-5){
	if(node_y[m]>node_y[j1]){
	  psi[m]=psi[j1]+hossz(m,j1)*v;
	}else{ 
	  psi[m]=psi[j1]-hossz(m,j1)*v;
	}
      }  
    }else{
      //az n-edik háromszög az egyik oldali szomszéd
      //először azt kell tisztázni, hogy a j1 vagy a j2 tartozik-e hozzá
      int q;
      if(A[n]==j1||B[n]==j1||C[n]==j1){
	for(q=ija[k];q<ija[k+1];q++) if(ija[q]==n) break;
	psi[m]=psi[j1]+sa[q]*(x[n]-x[k])*(((node_x[m]-node_x[j1])*(node_y[j2]-node_y[j1])>(node_y[m]-node_y[j1])*(node_x[j2]-node_x[j1]))?1:-1);
      }else{
	for(q=ija[k];q<ija[k+1];q++) if(ija[q]==n) break;
	psi[m]=psi[j2]+sa[q]*(x[n]-x[k])*(((node_x[m]-node_x[j2])*(node_y[j1]-node_y[j2])>(node_y[m]-node_y[j2])*(node_x[j1]-node_x[j2]))?1:-1);
      }
    }

    ipszi[k][0]=-1;
    j=k;//miután a k-adik háromszög csúcsaiban psi-t meghatároztuk, 
    //ez lesz a következő j. 
    ell[jj]=m;
     for(int i=1;i<jj;i++) if(ell[i]==m) cout<<"A(z) "<<m<<"-edik csúcs ismételten szerepel!\n";
    //cout<<j<<"  "<<m<<"  "<<psi[m]<<endl;
    }
  /* for(int i=0;i<=NN;i++){ 
    cout<<i<<"  "<<ipszi[i][0]<<"  "<<ipszi[i][1]<<"  "<<ipszi[i][2]<<"  "<<ipszi[i][3];
    if(ipszi[i][0]>-1){
      cout<<"  "<<node_x[A[i]]<<"  "<<node_y[A[i]]<<"  "<<psi[A[i]]<<endl;
    }else{
      cout<<endl;
    }
    }*/
  //cout<<"---------------------------------------------\n";
  //for(int i=1;i<=M;i++) cout<<i<<"  "<<psi[i]<<endl;
   


    FILE * graphics,*g2; 
    char cim[40];
    sprintf(cim,"Potential flow");
    graphics=initgr(1200,600,0,0,cim);
    fprintf(graphics,"set term x11\n");
    fprintf(graphics,"set xrange[%f:%f]\n",Lxmin,Lxmax);
    fprintf(graphics,"set yrange[%f:%f]\n",Lymin,Lymax);
    fprintf(graphics,"set contour\n");
    fprintf(graphics,"unset surface\n");
    fprintf(graphics,"set view 0,0,1.5,1\n");
   
    
  
    double h=.1,**xpsi;
    int Nx=(int)(Lx/h),Ny=(int)(Ly/h),nvonal=17;
    xpsi=dmatrix(0,Nx,0,Ny);

    double xpsimin=1.e10,xpsimax=-1.e10;

    for(int i=1;i<=M;i++){
      if(xpsimin>psi[i]) xpsimin=psi[i];
      if(xpsimax<psi[i]) xpsimax=psi[i];
    }
    //fprintf(graphics,"set zrange[%f:%f]\n",xpsimin,xpsimax);

    fprintf(graphics,"set cntrparam levels inc %f, %f, %f \n",xpsimin-(xpsimax-xpsimin)/nvonal,(xpsimax-xpsimin)/nvonal,xpsimax);

    /*for(int i=0;i<=NN;i++){
      if(xpsimin>x[i]) xpsimin=x[i];
      if(xpsimax<x[i]) xpsimax=x[i];
      }*/

    for(int i=0;i<=Nx;i++){
      for(int j=0;j<=Ny;j++){
	xpsi[i][j]=xpsimin-(xpsimax-xpsimin)/nvonal;
      }
    }
    interpol(xpsi,h);
    //interpol_phi(x,xpsi,h);
 


    fprintf(graphics,"splot '-' u 1:2:3 w l lc 2\n");
    for(int i=0;i<=Nx;i++){
      for(int j=0;j<=Ny;j++){
	x1=i*h+Lxmin;
	y1=j*h+Lymin;
	fprintf(graphics,"%f %f %f\n",x1,y1,xpsi[i][j]);
      }
      fprintf(graphics," \n");
    }
    fprintf(graphics,"e \n");

    

    closegr(graphics);

    free_dvector(psi,1,M);
    free_ivector(A,0,NN);
    free_ivector(B,0,NN);
    free_ivector(C,0,NN);
    free_lvector(ija,1,4*NN);
    free_dvector(sa,1,4*NN);
    free_dvector(b,1,NN);
    free_dvector(b0,1,NN);
    free_dvector(x,1,NN);
    free_imatrix(ipszi,0,NN,0,3);
    free_dmatrix(xpsi,0,Nx,0,Ny);
}
