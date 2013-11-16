#define BOOST_TEST_MODULE LineSearchTests
#include <cmath>

#include <boost/test/unit_test.hpp>

#include "optimize/function.hpp"
#include "optimize/line_search.hpp"

using namespace boost::numeric::ublas;
using namespace optimize;

BOOST_AUTO_TEST_CASE(Pol3_Minimum_Test)
{
  // a == b || a == c || b == c => -1
  BOOST_CHECK_EQUAL( get_pol3_minimum(0.,1.,-0.1,0.,1.,1.,1.), -1. );
  BOOST_CHECK_EQUAL( get_pol3_minimum(0.,1.,-0.1,1.,1.,0.,1.), -1. );
  BOOST_CHECK_EQUAL( get_pol3_minimum(0.,1.,-0.1,1.,1.,1.,1.), -1. );

  // first order polynom => c3 = 0 => -1
  BOOST_CHECK_EQUAL( get_pol3_minimum(0.,1.,-0.1,1.,0.9,2.,0.8), -1. );

  // second order polynom => c3 = 0 => -1
  BOOST_CHECK_EQUAL( get_pol3_minimum(0.,1.,-4.,1.,1.,2.,9.), -1. );

  // third order polynom
  double sol = (5. - std::sqrt(13))/3.;
  BOOST_CHECK_CLOSE( get_pol3_minimum(0.,1.,-4.,1.,1.,2.,5.), sol, 0.0001 );
}

BOOST_AUTO_TEST_CASE(Pol2_Minimum_Test)
{
  // a == b => -1
  BOOST_CHECK_EQUAL( get_pol2_minimum(0.,1.,-0.1,0.,1.), -1. );

  // first order polynom => -1
  BOOST_CHECK_EQUAL( get_pol2_minimum(0.,1.,-0.1,1.,0.9), -1. );

  // second order polynom, f(a) == f(b) => solution = a + (a + b)/2
  BOOST_CHECK_CLOSE( get_pol2_minimum(0.1,1.,-0.1,0.7,1.), 0.4, 0.0001 );

  // second order polynom, f(b) < f(a) + df(a) * (b - a) => no minimum => -1
  BOOST_CHECK_EQUAL( get_pol2_minimum(0.,1.,-0.4,1.,0.), -1. );
}

double fpol2(double* x)
{
  return 4. * std::pow(x[0] - 0.5, 2);
}

BOOST_AUTO_TEST_CASE(Line_Search_Zoom_Test)
{
  vector<double> x0(1);
  x0[0] = 0.;
  vector<double> p0(1);
  p0[0] = 1.; // Newton search direction = -f'(x)/f''(x)
  function func(&fpol2, 1);
  double fstar = 0.;

  // Parabola between 0 and 1 with a minimum at 0.5 
  BOOST_CHECK_CLOSE( line_search_zoom(&func,x0,p0,0.,1.,-4.,
                                      1.,1.,1.,-4.,fstar), 0.5, 0.0001 );
  BOOST_CHECK_CLOSE( fstar, 0., 0.0001 );
}

BOOST_AUTO_TEST_CASE(Line_Search_Wolfe_Test)
{
  vector<double> x0(1);
  x0[0] = 0.;  
  vector<double> p0(1);
  p0[0] = 0.5 - x0[0]; // Newton search direction = -f'(x)/f''(x)
  function func(&fpol2, 1);
  double fstar = 0.;

  // Parabola between 0 and 1 with a minimum at 0.5 
  BOOST_CHECK_CLOSE( line_search_wolfe(&func,x0,p0,fstar), 1., 0.0001 );
  BOOST_CHECK_CLOSE( fstar, 0., 0.0001 );
}

// EOF
