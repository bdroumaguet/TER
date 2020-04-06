#ifndef _TIME_SCHEME_CPP

#include "Timescheme.h"


using namespace Eigen;
using namespace std;


TimeScheme::TimeScheme():_fct(0),_flx(0)
{}

TimeScheme::~TimeScheme(){}

void TimeScheme::Initialize(int N, Vector2d Uleft, Vector2d Uright, fonction* fct, Flux* flx)
{
  _fct=fct;
  _flx = flx;
_N=N;
_U0.resize(_N,2);
_Ul=Uleft;
_Ur=Uright;
for (int i=0;i<_N;i++)
{
  if (i<N/2.)
{
    _U0(i,0)=Uleft[0];
    _U0(i,1)=Uleft[1];
  }
  else
  {
      _U0(i,0)=Uright[0];
      _U0(i,1)=Uright[1];
    }
}
}

void TimeScheme::Advance(double dx, Matrix<double,Dynamic,Dynamic> U, int i)
{
  Vector2d F_uiplus1demi, F_uimoins1demi,Ui,Uiplus1,Uimoins1;

  Ui[0]=U(i,0); Ui[1]=U(i,1);
  // clacul de max uin, u
  double maxb(_fct->max_b(U));
  _dt=0.5*dx/maxb;
  double Dmmhi, Dmmui;
  double himoins1demi,hiplus1demi;



  Vector2d Uiplus1demi,Uimoins1demi;
    if (i==0)
    {
      Uiplus1[0]=U(i+1,0); Uiplus1[1]=U(i+1,1);
     Dmmhi=_fct->minmod((U(i,0)-_U0(0,0))/dx,(U(i+1,0)-U(i,0))/dx);// ordre 2
      Dmmui=_fct->minmod((U(i,1)-_U0(0,1))/dx,(U(i+1,1)-U(i,1))/dx);// ordre 2



      hiplus1demi=U(i,0)+0.5*dx*Dmmhi;// ordre 2
      Uiplus1demi[0]=hiplus1demi; Uiplus1demi[1]=U(i,1)+0.5*dx*Dmmui*hiplus1demi/U(i,0);// ordre 2
      himoins1demi=U(i,0)-0.5*dx*Dmmhi;// ordre 2
      Uimoins1demi[0]=himoins1demi; Uimoins1demi[1]=U(i,1)-0.5*dx*Dmmui*himoins1demi/U(i,0);// ordre 2



      _flx->Rusanov(Uiplus1,Uiplus1demi,_fct); //ordre 2
    //  _flx->Rusanov(Uiplus1,Ui,_fct);
    //  F_uiplus1demi=_flx->GetFlux();
    _flx->Rusanov(Uimoins1demi,_Ul,_fct); //ordre2
    //_flx->Rusanov(Ui,_Ul,_fct);
      //F_uimoins1demi=_flx->GetFlux();
    }
    else if (i==_N-1)
    {

            Uimoins1[0]=U(i-1,0); Uimoins1[1]=U(i-1,1);

     Dmmhi=_fct->minmod((U(i,0)-U(i-1,0))/dx,(_U0(_N-1,0)-U(i,0))/dx);// ordre 2*
      Dmmui=_fct->minmod((U(i,1)-U(i-1,1))/dx,(_U0(_N-1,1)-U(i,1))/dx);// ordre 2*

      hiplus1demi=U(i,0)+0.5*dx*Dmmhi;// ordre 2*
      Uiplus1demi[0]=hiplus1demi; Uiplus1demi[1]=U(i,1)+0.5*dx*Dmmui*hiplus1demi/U(i,0);// ordre 2*
      himoins1demi=U(i,0)-0.5*dx*Dmmhi; // ordre 2
      Uimoins1demi[0]=himoins1demi; Uimoins1demi[1]=U(i,1)-0.5*dx*Dmmui*himoins1demi/U(i,0); // ordre 2
    _flx->Rusanov(_Ur,Uiplus1demi,_fct); // ordre 2
      //  _flx->Rusanov(_Ur,Ui,_fct);
      //F_uiplus1demi=_flx->GetFlux();
      _flx->Rusanov(Uimoins1demi,Uimoins1,_fct); // ordre 2
    //  _flx->Rusanov(Ui,Uimoins1,_fct);
    //  F_uimoins1demi=_flx->GetFlux();

    }
    else
    {      Uiplus1[0]=U(i+1,0); Uiplus1[1]=U(i+1,1);
          Uimoins1[0]=U(i-1,0); Uimoins1[1]=U(i-1,1);

     Dmmhi=_fct->minmod((U(i,0)-U(i-1,0))/dx,(U(i+1,0)-U(i,0))/dx); //ordre2
      Dmmui=_fct->minmod((U(i,1)-U(i-1,1))/dx,(U(i+1,1)-U(i,1))/dx);//ordre2

      hiplus1demi=U(i,0)+0.5*dx*Dmmhi;//ordre2
      Uiplus1demi[0]=hiplus1demi; Uiplus1demi[1]=U(i,1)+0.5*dx*Dmmui*hiplus1demi/U(i,0);//ordre2
      himoins1demi=U(i,0)-0.5*dx*Dmmhi;//ordre2
      Uimoins1demi[0]=himoins1demi; Uimoins1demi[1]=U(i,1)-0.5*dx*Dmmui*himoins1demi/U(i,0);//ordre2
     _flx->Rusanov(Uiplus1,Uiplus1demi,_fct);//ordre2
    //  _flx->Rusanov(Uiplus1,Ui,_fct); //ordre1
      //F_uiplus1demi=_flx->GetFlux();
      _flx->Rusanov(Uimoins1demi,Uimoins1,_fct);//ordre2
    //  _flx->Rusanov(Ui,Uimoins1,_fct); //ordre1
    //  F_uimoins1demi=_flx->GetFlux();
    }

    _Unext[0]=U(i,0)-_dt/dx*(F_uiplus1demi[0]-F_uimoins1demi[0]);
    _Unext[1]=U(i,1)-_dt/dx*(F_uiplus1demi[1]-F_uimoins1demi[1]);
}



#define _TIME_SCHEME_CPP
#endif
