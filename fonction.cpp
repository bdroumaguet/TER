#include "fonction.h"




using namespace Eigen;
using namespace std;




fonction :: fonction() // constructeur par défaut
{}
fonction ::  ~fonction() // destructeur par défaut
{}

  void  fonction::f( double q,double h)
  {
    _vect[0]=q;
    if (h<pow(10,-16))
    {
      _vect[1]=0;
    }
    else
    {_vect[1]=q*q*1./(h*1.) + 9.81*h*h/2.;}
  }

  double fonction::max_b(MatrixXd U)
  {
    double m1(0),m2(0),m(0);
    for (int i=0; i<U.rows();i++)
      {
          m1 = max(sqrt(pow(U(i,1)/U(i,0)+sqrt(9.81*U(i,0)),2)),sqrt(pow(U(i,1)/U(i,0)-sqrt(9.81*U(i,0)),2)));
        m = max(m1,m);}
    return m;
  }

  double fonction::max_lambda(Vector2d u1, Vector2d u2)
  {
    double m1(0),m2(0),m(0);
    if (sqrt(pow(u1[0],2))>pow(10,-16))
    {m1 = max(sqrt(pow(u1[1]/u1[0]+sqrt(9.81*u1[0]),2)),sqrt(pow(u1[1]/u1[0]-sqrt(9.81*u1[0]),2)));}
    if (sqrt(pow(u2[0],2))>pow(10,-16))
    {m2 = max(sqrt(pow(u2[1]/u2[0]+sqrt(9.81*u2[0]),2)),sqrt(pow(u2[1]/u2[0]-sqrt(9.81*u2[0]),2)));}
    m = max(m1,m2);
    return m;
  }

  double fonction::min_lambda(Vector2d u1, Vector2d u2)
  {
    double m1(0),m2(0),m(0);
    if (sqrt(pow(u1[0],2))>pow(10,-16))
    {m1 = min(sqrt(pow(u1[1]/u1[0]+sqrt(9.81*u1[0]),2)),sqrt(pow(u1[1]/u1[0]-sqrt(9.81*u1[0]),2)));}
    if (sqrt(pow(u2[0],2))>pow(10,-16))
    {m2 = min(sqrt(pow(u2[1]/u2[0]+sqrt(9.81*u2[0]),2)),sqrt(pow(u2[1]/u2[0]-sqrt(9.81*u2[0]),2)));}
    m = min(m1,m2);
    return m;
  }

  double fonction::minmod(double x,double y)
  {
    double minmod_xy;
    if ((x>0)&&(y>0))
    {
      minmod_xy=min(x,y);
    }
    else if ((x<0)&&(y<0))
    {
      minmod_xy=max(x,y);
    }
    else
    {
      minmod_xy=0;
    }
    return minmod_xy;
  }
