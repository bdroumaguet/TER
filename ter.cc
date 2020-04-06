#include <iostream>
#include "Dense"
#include "Sparse"
#include <complex>
#include <chrono>
#include "Timescheme.h"
using namespace Eigen;
using namespace std;


int main()
{
  Vector2d F1,F2;

  TimeScheme* tim(0);
  fonction* fct(0);
  Flux* flx(0);
  tim = new TimeScheme();
  fct = new fonction();
  flx = new Flux();
  int T=300;
  int N=100;
  double dx(0);
  dx=1./N;
  Matrix<double,Dynamic,Dynamic> U, U_1,U0,UTsur2;
  U.resize(N,2);
  U_1.resize(N,2);
  U0.resize(N,2);
  UTsur2.resize(N,2);
  Vector2d Ug,Ud; Ug[0]=6; Ug[1]=0; Ud[0]=2; Ud[1]=0;
  Vector2d Unext;
  string file_name("water_height200_ordre1T300.txt");
  string file_name2("water_flow200_ordre1T300.txt");




  tim->Initialize(N,Ug,Ud,fct,flx);

  U0=tim->Get_U0();
  U=U0;

  Vector2d sol_ex;
  double h0(Ug[0]-Ud[0]);
  double err(0);
  double t;
  ofstream file1;
  file1.open(file_name,ios::out);
  ofstream file2;
  file2.open(file_name2,ios::out);


  for (int i=0 ; i<T; i++)
  {


      for (int j=0 ; j<N ; j++)
        {
          cout << t << endl;
          Unext.resize(2);

          tim->Advance(dx,U,j);
          Unext=tim->Get_Unext();
          U_1(j,0)=Unext[0];
          U_1(j,1)=Unext[1];
          if ((t>2.0)&&(t<2.06))
          {
            cout << "ouiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii" << endl;
            double x = j*dx*100-50;
          /*  if (x<-t*sqrt(9.81*h0))
              {
                sol_ex = Ug;
              }
            if (x>2*t*sqrt(9.81*h0))
            {
              sol_ex=Ud;
            }
            if ((x<2*t*sqrt(9.81*h0)) && (x>-t*sqrt(9.81*h0)))
            {
              sol_ex[0]=4.0/(9.0*9.81)*pow(sqrt(9.81*h0)-x/(2.0*t),2)+Ud[0];
              sol_ex[1]=2.0/3.0*(x*1.0/t+sqrt(9.81*h0));
            }*/
            file1 << x << " " << U0(j,0) << " " << sol_ex[0]  << " " << Unext[0] << endl;
            file2 << x << " " << U0(j,1) << " " << sol_ex[1]/sol_ex[0] << " " << Unext[1]/Unext[0] << endl;
            //UTsur2(j,0)=U_1(j,0);
            //UTsur2(j,1)=U_1(j,1);

          }

        }
        t += tim->Get_dt()*100;
        U=U_1;


  }
  file2.close();
  file1.close();
delete tim;
delete fct;
delete flx;

TimeScheme* tim2(0);
fonction* fct2(0);
Flux* flx2(0);

tim2 = new TimeScheme();
fct2 = new fonction();
flx2 = new Flux();

double dx2(0);
dx2=1./(2*N);
Matrix<double,Dynamic,Dynamic> U2, U_12,U02;
U2.resize(2*N,2);
U_12.resize(2*N,2);
U02.resize(2*N,2);

tim2->Initialize(2*N,Ug,Ud,fct2,flx2);
U02=tim2->Get_U0();
U2=U02;
err=0;
t=0;
for (int i=0 ; i<T; i++)
{
  cout << t << endl;
    for (int j=0 ; j<2*N ; j++)
      {

        Unext.resize(2);
        tim2->Advance(dx2,U2,j);
        Unext=tim2->Get_Unext();
        U_12(j,0)=Unext[0];
        U_12(j,1)=Unext[1];
        if ((t>2.0)&&(t<2.02))
        {
          cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << endl;
          if (sqrt(pow((int(j*1./2.))-j*1./2.,2))<0.01)
          {
            err+=pow(sqrt(pow(U_12(j,0)-U_1(int(j/2),0),2)),2)*dx;
          }

        }
      }
      t+= tim2->Get_dt()*100;
      U2=U_12;
}
err=sqrt(err);
cout << "erreur= "<< err << endl;


delete tim2;
delete fct2;
delete flx2;
return 0;
}
