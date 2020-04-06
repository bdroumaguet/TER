#include "Dense"
#include "Sparse"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>


class fonction
{
protected:
  Eigen::Vector2d _vect;
public:
  fonction(); // constructeur par défaut
  virtual ~fonction(); // destructeur par défaut

  void f( double q, double h);
  Eigen::Vector2d & Get_f() {return _vect;};
  double max_b(Eigen::MatrixXd U);
  double max_lambda(Eigen::Vector2d u1,Eigen::Vector2d u2) ;
  double min_lambda(Eigen::Vector2d u1,Eigen::Vector2d u2) ;
  double minmod(double x,double y);
};
