/*!\file   line_search.hpp
 * \brief  Functions used to perform minimization along a line
 * \author Renaud Bruneliere
 * \date   04.11.2013
 *
 * The header contains different functions which are used to find an approximate
 * step length along a line. These functions are used in quasi-Newton methods.
 */

#ifndef LINE_SEARCH_HPP
#define LINE_SEARCH_HPP

#include "optimize/function.hpp"

namespace optimize {

//! Find a step length along a line satisfying the strong wolfe conditions.
/*!
  \param function pointer to a function
  \param xk       starting position
  \param pk       line (descent) direction
  \param fstar    value of f correspoding to the step length (output)
  \param eps      increment to determine the function gradient
  \param c1       wolfe parameter for the first condition (decrease)
  \param c2       wolfe parameter for the second condition (curvature)
  \param c3       extension parameter
  \return         step length alpha (in unit of pk)
 */
double line_search_wolfe(
    const function * const fptr,
    const boost::numeric::ublas::vector<double>& xk,
    const boost::numeric::ublas::vector<double>& pk,
    double& fstar, double c1=1e-4, double c2=0.9, double c3=2.);

//! Find step length satisfying wolfe conditions in interval alpha_lo, alpha_hi
double line_search_zoom(
    const function * const fptr,
    const boost::numeric::ublas::vector<double>& xk,
    const boost::numeric::ublas::vector<double>& pk,
    double alpha_lo, double flo, double dflo, double alpha_hi, double fhi,
    double f0, double df0, double& fstar, double c1=1e-4, double c2=0.9);

//! Find the minimum between 3 points via a third order polynom
double get_pol3_minimum(double a, double fa, double dfa, 
                        double b, double fb,
                        double c, double fc);

//! Find the minimum between 2 points via a second order polynom
double get_pol2_minimum(double a, double fa, double dfa,
                        double b, double fb);

} // namespace optimize

#endif // LINE_SEARCH_HPP
