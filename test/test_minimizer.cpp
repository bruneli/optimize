#define BOOST_TEST_MODULE MinimizerTests
#include <cmath>
#include <stdexcept>

#include <boost/bind.hpp>
#include <boost/scoped_array.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "optimize/minimizer_nm.hpp"
#include "optimize/minimizer_bfgs.hpp"

using namespace optimize;
using namespace boost::numeric::ublas;

double quad_xyz(double * const x)
{
  double val = std::pow(x[0] - 0.5, 2);
  val += 10. * std::pow(x[1] - 2., 2);
  val += 4. * std::pow(x[2] + 1., 2);
  return val;
}

void quad_xyz_gradient(const vector<double>& x, vector<double>& grad)
{
  grad[0] = 2. * (x[0] - 0.5);
  grad[1] = 20. * (x[1] - 2.);
  grad[2] = 8. * (x[2] + 1.);
}

class quad_func : public function
{
 public:
  quad_func(size_t n, double * const x0, double * const c);

  // Overload the function::eval
  virtual double eval(double * const x) const;

  // Static function used as a wrapper to access the class method eval
  static double evaluate(double * const x, void * obj);

  // Compute exact gradient
  void grad(const vector<double>& x, vector<double>& grad) const;

 private:
  boost::scoped_array<double> x0_;
  boost::scoped_array<double> c_;
};

quad_func::quad_func(size_t n, double * const x0, double * const c) :
  function(n)
{
  x0_.reset(new double[this->get_n()]);
  c_.reset(new double[this->get_n()]);
  for (size_t i = 0; i < this->get_n(); i++) {
    x0_[i] = x0[i];
    c_[i] = c[i];
  }
}

double quad_func::eval(double * const x) const
{
  double val = 0.;
  for (size_t i = 0; i < this->get_n(); i++) 
    val += c_[i] * std::pow(x[i] - x0_[i], 2);
  return val;
}

double quad_func::evaluate(double * const x, void * obj)
{
  quad_func * qf = (quad_func*) obj;
  return qf->eval(x);
}

void quad_func::grad(const vector<double>& x, vector<double>& grad) const
{
  for (size_t i = 0; i < this->get_n(); i++)
    grad[i] = 2. * c_[i] * (x[i] - x0_[i]);
}

BOOST_AUTO_TEST_CASE(Base_Minimizer_Test)
{
  size_t n = 3;
  boost::scoped_array<double> x(new double[n]);
  for (size_t i = 0; i < n; i++)
    x[i] = 2.;
  double x0[3] = {0.5, 2., -1.};
  double c[3]  = {1., 10., 4.};

  // Show the various ways to construct a minimizer with an objective function
  // minimizer without any objective function
  minimizer_bfgs m0(n);
  // minimizer built from an optimize::function
  quad_func f1(n, (double*)x0, (double*)c);
  minimizer_bfgs m1(&f1);
  // minimizer pointing to a C-like function address
  minimizer_bfgs m2(&quad_xyz, n);
  // minimizer pointing to a static wrapper function and its class instance
  minimizer_bfgs m3(&quad_func::evaluate, &f1, n);
  // minimizer linked to a boost::function which is binding a class method to
  // its object
  boost::function<double (double * const)> feval;
  feval = boost::bind(&quad_func::eval, &f1, _1);
  minimizer_bfgs m4(feval, n);

  // check number of variables
  BOOST_CHECK_EQUAL( m0.get_function()->get_n(), n );
  BOOST_CHECK_EQUAL( m1.get_function()->get_n(), n );
  BOOST_CHECK_EQUAL( m2.get_function()->get_n(), n );
  BOOST_CHECK_EQUAL( m3.get_function()->get_n(), n );
  BOOST_CHECK_EQUAL( m4.get_function()->get_n(), n );

  // check function associated to a minimizer
  BOOST_CHECK_THROW( m0.get_function()->eval(x), std::runtime_error );
  BOOST_CHECK_EQUAL( m1.get_function(), &f1 );
  BOOST_CHECK_EQUAL( m2.get_function()->eval(x), quad_xyz(x.get()) );
  BOOST_CHECK_EQUAL( m3.get_function()->eval(x), f1.eval(x.get()) );
  BOOST_CHECK_EQUAL( m4.get_function()->eval(x), f1.eval(x.get()) );
}

BOOST_AUTO_TEST_CASE(Minimizer_BFGS_Test)
{
  size_t n = 3;
  boost::scoped_array<double> x(new double[n]);
  for (size_t i = 0; i < n; i++)
    x[i] = 3.;
  double x0[3] = {0.5, 2., -1.};
  double c[3]  = {1., 10., 4.};

  // Show the various ways to construct a minimizer with an objective function
  // minimizer without any objective function
  minimizer_bfgs m0(n);
  // minimizer built from an optimize::function
  quad_func f1(n, (double*)x0, (double*)c);
  minimizer_bfgs m1(&f1);
  // minimizer pointing to a C-like function address
  minimizer_bfgs m2(&quad_xyz, n);
  // minimizer pointing to a static wrapper function and its class instance
  minimizer_bfgs m3(&quad_func::evaluate, &f1, n);
  // minimizer linked to a boost::function which is binding a class method to
  // its object
  boost::function<double (double * const)> feval;
  feval = boost::bind(&quad_func::eval, &f1, _1);
  minimizer_bfgs m4(feval, n);

  // check number of variables
  BOOST_CHECK_EQUAL( m0.get_function()->get_n(), n );
  BOOST_CHECK_EQUAL( m1.get_function()->get_n(), n );
  BOOST_CHECK_EQUAL( m2.get_function()->get_n(), n );
  BOOST_CHECK_EQUAL( m3.get_function()->get_n(), n );
  BOOST_CHECK_EQUAL( m4.get_function()->get_n(), n );

  // check function associated to a minimizer
  BOOST_CHECK_THROW( m0.get_function()->eval(x), std::runtime_error );
  BOOST_CHECK_EQUAL( m1.get_function(), &f1 );
  BOOST_CHECK_EQUAL( m2.get_function()->eval(x), quad_xyz(x.get()) );
  BOOST_CHECK_EQUAL( m3.get_function()->eval(x), f1.eval(x.get()) );
  BOOST_CHECK_EQUAL( m4.get_function()->eval(x), f1.eval(x.get()) );

  // minimize
  BOOST_CHECK_THROW( m0.minimize(), std::runtime_error );
  m0.minimize(&quad_xyz, x.get());
  m1.minimize(x.get());
  m2.minimize(x.get());
  m3.minimize(x.get());
  m4.minimize(x.get());
  for (size_t i = 0; i < n; i++) {
    BOOST_CHECK_CLOSE( m0.get_opt_var(i).value(), x0[i], 0.01 );
    BOOST_CHECK_CLOSE( m1.get_opt_var(i).value(), x0[i], 0.01 );
    BOOST_CHECK_CLOSE( m2.get_opt_var(i).value(), x0[i], 0.01 );
    BOOST_CHECK_CLOSE( m3.get_opt_var(i).value(), x0[i], 0.01 );
    BOOST_CHECK_CLOSE( m4.get_opt_var(i).value(), x0[i], 0.01 );
  }  
}

// EOF
