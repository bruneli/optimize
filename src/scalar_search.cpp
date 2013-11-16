#include <cmath>
#include "optimize/scalar_search.hpp"

int optimize::bracket(const optimize::t_func & func, double x0, double h, 
                      double & a, double & b, double c, unsigned int max_iter)
{
  optimize::function function(func);
  return optimize::bracket(&function, x0, h, a, b, c, max_iter);
}

int optimize::bracket(optimize::t_boofunc & func, double x0, double h,
                      double & a, double & b, double c, unsigned int max_iter)
{
  optimize::function function(func);
  return optimize::bracket(&function, x0, h, a, b, c, max_iter);
}

int optimize::bracket(const optimize::t_func_obj & func, void * obj,
                      double x0, double h, double & a, double & b,
                      double c, unsigned int max_iter)
{
  optimize::function function(func, obj);
  return optimize::bracket(&function, x0, h, a, b, c, max_iter);
}

int optimize::bracket(optimize::function * const func,
                      double x0, double h, double & a, double & b,
                      double c, unsigned int max_iter)
{
  h = std::fabs(h);
  boost::scoped_array<double> x1(new double[1]);
  boost::scoped_array<double> x2(new double[1]);
  x1[0] = x0;
  double f1 = func->eval_ubd(x1);
  x2[0] = x0 + h;
  double f2 = func->eval_ubd(x2);

  // Check downhill direction
  if (f2 > f1) {
    h = -h;
    x2[0] = x0 + h;
    f2 = func->eval_ubd(x2);
    if (f2 > f1) {
      // [x0 - h, x0 + h] is bracketing a (local) minimum
      a = x2[0];
      b = x0 - h;
      return 0;
    }
  }

  // Iterate till max iterations is reached or a local minimum is found
  unsigned int iter = 0;
  bool cond = true;
  while (cond) {
    h = c*h; // Extend step at each iteration
    x1[0] = x2[0]; f1 = f2;
    x2[0] += h;
    f2 = func->eval_ubd(x2);
    cond = ((iter++) < max_iter) && (f2 < f1);
  }

  if (x1[0] < x2[0]) {
    a = x1[0];
    b = x2[0];
  } else {
    a = x2[0];
    b = x1[0];
  }
  return (iter >= max_iter);
}

double optimize::minimize_1d_golden(const optimize::t_func & func, 
                                    double a, double b,
                                    double & fmin, double tol)
{
  optimize::function function(func);
  return optimize::minimize_1d_golden(&function, a, b, fmin, tol);
}

double optimize::minimize_1d_golden(optimize::t_boofunc & func, 
                                    double a, double b,
                                    double & fmin, double tol)
{
  optimize::function function(func);
  return optimize::minimize_1d_golden(&function, a, b, fmin, tol);
}

double optimize::minimize_1d_golden(const optimize::t_func_obj & func,
                                    void * obj, 
                                    double a, double b,
                                    double & fmin, double tol)
{
  optimize::function function(func, obj);
  return optimize::minimize_1d_golden(&function, a, b, fmin, tol);
}

double optimize::minimize_1d_golden(optimize::function * const func,
                                    double a, double b,
                                    double & fmin, double tol)
{
  // The stopping rule is defined as:
  //    abs(b - a)*pow(r,n_iter) = tol
  // => n_iter = log(tol/abs(b - a))/log(r)
  // with r = (-1 + sqrt(5))/2 the Golden ratio
  size_t n_iter = size_t(-2.078087 * std::log(tol / std::fabs(b - a))) + 1;
  double r = (-1. + std::sqrt(5.))/2.;
  double c = 1. - r;

  // First telescoping
  boost::scoped_array<double> x1(new double[1]);
  x1[0] = r*a + c*b;
  double f1 = func->eval_ubd(x1);
  boost::scoped_array<double> x2(new double[1]);
  x2[0] = c*a + r*b;
  double f2 = func->eval_ubd(x2);

  // Iterations loop
  for (size_t i = 0; i < n_iter; i++) {
    if (f1 > f2) {
      a = x1[0];
      x1[0] = x2[0];
      f1 = f2;
      x2[0] = c*a + r*b;
      f2 = func->eval_ubd(x2);
    } else {
      b = x2[0];
      x2[0] = x1[0];
      f2 = f1;
      x1[0] = r*a + c*b;
      f1 = func->eval_ubd(x1);
    }
  }

  fmin = (f1 < f2) ? f1 : f2;
  if (f1 < f2) return x1[0];
  return x2[0];
}
