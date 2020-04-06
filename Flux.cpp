#include "Flux.h"




using namespace Eigen;
using namespace std;




Flux :: Flux() // constructeur par défaut
{}
Flux ::  ~Flux() // destructeur par défaut
{}

  void Flux::Rusanov(Vector2d u1, Vector2d u2, fonction* fct)
  {
    double maxlambda(fct->max_lambda(u1,u2));
    Vector2d f_u1,f_u2;
    fct->f(u1[1],u1[0]);
    f_u1=fct->Get_f();


    fct->f(u2[1],u2[0]);
    f_u2=fct->Get_f();
    _F_uv=0.5*(f_u1+f_u2)-maxlambda*0.5*(u1-u2);
  }
