#ifndef _FLUX_H


#include "fonction.h"

class Flux
{
private:

  Eigen::Vector2d _F_uv;

public:
  Flux();
  ~Flux();

  void Rusanov(Eigen::Vector2d u1, Eigen::Vector2d u2, fonction* fct);
  const Eigen::Vector2d & GetFlux(){return _F_uv;};
};


#define _Flux_H
#endif
