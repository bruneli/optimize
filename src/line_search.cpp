#include <cmath>
#include "optimize/line_search.hpp"

using namespace boost::numeric::ublas;

double optimize::line_search_wolfe(const optimize::function * const fptr,
                                   const vector<double>& xk,
                                   const vector<double>& pk, double& fstar,
                                   double c1, double c2, double c3)
{
  // Define alpha0
  double alpha0 = 0.;
  double fa0 = fptr->eval_ubd(xk);
  vector<double> grad0(fptr->get_n());
  fptr->get_approx_gradient(xk, fa0, grad0, true);
  double dfa0 = inner_prod(pk, grad0);
  double f0  = fa0;
  double df0 = dfa0;

  // Define alpha1
  double alpha1 = 1.; // Standard choice for quasi-Newton methods
  double fa1 = 0.;
  vector<double> x1(fptr->get_n());
  vector<double> grad1(fptr->get_n());
  double dfa1 = 0.;

  double alpha_star = alpha1; // candidate step length

  // Iterative loop: try to bracket a local minimum
  unsigned int niter = 10;
  for (unsigned int i = 0; i < niter; i++) {
    x1 = xk + alpha1 * pk;
    fa1 = fptr->eval_ubd(x1);

    // Test sufficient decrease condition
    if ((fa1 > f0 + c1 * alpha1 * df0) ||
        ((i > 1) && (fa1 >= fa0))) {
      // If not satisfied, find minimum between alpha0, alpha1
      alpha_star = line_search_zoom(fptr, xk, pk, alpha0, fa0, dfa0, 
                                    alpha1, fa1, f0, df0, fstar, c1, c2);
      break;
    }

    // Test the curvature condition
    fptr->get_approx_gradient(x1, fa1, grad1);
    dfa1 = inner_prod(pk, grad1);
    if (std::fabs(dfa1) <= -c2 * df0) {
      // alpha1 is satisfying both conditions, stop iterations
      alpha_star = alpha1;
      fstar = fa1;
      break;
    }

    // Check if high value has a positive derivative
    if (dfa1 >= 0.) {
      // If that is the case, find minimum between alpha0, alpha1
      alpha_star = line_search_zoom(fptr, xk, pk, alpha0, fa0, dfa0, 
                                    alpha1, fa1, f0, df0, fstar, c1, c2);
      break;
    }

    // If none of the conditions above are satisfyed, extend the interval
    alpha_star = c3 * alpha1;
    alpha0 = alpha1;
    alpha1 = alpha_star;
    fa0 = fa1;
    dfa0 = dfa1;
  }

  return alpha_star;
}

double optimize::line_search_zoom(const optimize::function * const fptr,
                                  const vector<double>& xk,
                                  const vector<double>& pk,
                                  double alpha_lo, double flo, double dflo,
                                  double alpha_hi, double fhi, 
                                  double f0, double df0, double& fstar,
                                  double c1, double c2)
{
  double alpha_star = alpha_lo;
  fstar = flo;
  vector<double> x_star(fptr->get_n());
  vector<double> grad_star(fptr->get_n());
  double dfstar = 0.;
  double alpha3 = alpha_lo;
  double f3 = flo;
  bool condition = true;
  unsigned int iter  = 0;
  unsigned int niter = 10;
  while (condition) {
    // Estimate alpha_star
    // try third order polynom (except at first iteration)
    if (iter)
      alpha_star = get_pol3_minimum(alpha_lo, flo, dflo, alpha_hi, fhi,
                                    alpha3, f3);
    // if it fails or too close to bounds, try second order polynom
    if (alpha_star < alpha_lo + 0.1 * (alpha_hi - alpha_lo) ||
        alpha_star > alpha_lo + 0.9 * (alpha_hi - alpha_lo)) {
      alpha_star = get_pol2_minimum(alpha_lo, flo, dflo, alpha_hi, fhi);
      // if it fails or too close to bounds, do bisection
      if (alpha_star < alpha_lo + 0.1 * (alpha_hi - alpha_lo) ||
          alpha_star > alpha_lo + 0.9 * (alpha_hi - alpha_lo))
        alpha_star = 0.5 * (alpha_lo + alpha_hi);
    }
    x_star = xk + alpha_star * pk;
    fstar = fptr->eval_ubd(x_star);

    // Test the strong Wolfe conditions
    // Test sufficient decrease condition
    if ((fstar > f0 + c1 * alpha_star * df0) || (fstar >= flo)) {
      // minimum is in the interval [alpha_lo, alpha_star]
      alpha3 = alpha_hi;
      f3 = fhi;
      alpha_hi = alpha_star;
      fhi = fstar;
    } else {
      fptr->get_approx_gradient(x_star, fstar, grad_star, true);
      dfstar = inner_prod(pk, grad_star);

      // Test the curvature condition
      if (std::fabs(dfstar) <= -c2 * df0) {
        // alpha_star is satisfying both conditions, stop iterations
        condition = false;
      } else if (dfstar * (alpha_hi - alpha_lo) >= 0.) {
        // minimum is in the interval [alpha_lo, alpha_star]
        alpha3 = alpha_hi;
        f3 = fhi;
        alpha_hi = alpha_star;
        fhi = fstar;
      } else {
        // minimum is in the interval [alpha_star, alpha_hi]
        alpha3 = alpha_lo;
        f3 = flo;
        alpha_lo = alpha_star;
        flo = fstar;
        dflo = dfstar;
      }
    } // End while

    if (iter++ > niter) condition = false;
  }

  return alpha_star;
}

double optimize::get_pol3_minimum(double a, double fa, double dfa,
                                  double b, double fb, double c, double fc)
{
  if (a == b || a == c || b == c) return -1.;
  // f(x) = c3 * (x - a)**3 + c2 * (x - a)**2 + c1 * (x - a) + c0
  double c0 = fa;
  double c1 = dfa;
  // With (b, fb) and (c, fc), solve 2 equations with 2 unknowns (c3, c2)
  double m12 = (b - a) * (b - a);
  double m11 = m12 * (b - a);
  double m22 = (c - a) * (c - a);
  double m21 = m22 * (c - a);
  double det = m11 * m22 - m12 * m21;
  if (det == 0.) return -1.;
  double dfb = fb - c1 * (b - a) - c0;
  double dfc = fc - c1 * (c - a) - c0;
  double c3 = (m22 * dfb - m12 * dfc) / det;
  double c2 = (m11 * dfc - m21 * dfb) / det;
  // Find the minimum
  // df(x)/dx = 3 * c3 * (x - a)**2 + 2 * c2 * (x - a) + c1 = 0
  if (c3 == 0.) return -1.;
  det = c2 * c2 - 3 * c1 * c3;
  if (det < 0.) return -1.;
  return a + (std::sqrt(det) - c2) / (3 * c3);
}

double optimize::get_pol2_minimum(double a, double fa, double dfa,
                                  double b, double fb)
{
  // f(x) = c2 * (x - a)**2 + c1 * (x - a) + c0
  double c0 = fa;
  double c1 = dfa;
  double d = b - a;
  if (d == 0.) return -1.;
  double c2 = (fb - c1 * d - c0) / (d * d);
  if (c2 <= 0.) return -1.;
  // x_min = a - c1 / (2 * c2)
  return (a - c1 / (2. * c2));
}
