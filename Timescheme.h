#ifndef _TIME_SCHEME_H

#include "Flux.h"

class TimeScheme
{

protected:
  int _N ;
  double _dt;
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> _U0;
  Eigen::Vector2d _Unext, _Ul, _Ur;

private:
  fonction* _fct;
  Flux* _flx;


public:
  TimeScheme();
  ~TimeScheme();
  void Initialize(int N, Eigen::Vector2d Uleft, Eigen::Vector2d Uright, fonction* fct, Flux* flx);
  const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> & Get_U0(){return _U0;};
  void Advance(double dx,Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> U, int i);
  const Eigen::Vector2d & Get_Unext(){return _Unext;};
  const double Get_dt(){return _dt;};
};


#define _TIME_SCHEME_H
#endif
