#define BOOST_TEST_MODULE Function_Variable_Tests
#include <cmath>
#include <stdexcept>

#include <boost/bind.hpp>
#include <boost/scoped_array.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "optimize/variable.hpp"
#include "optimize/function.hpp"

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

BOOST_AUTO_TEST_CASE(Variable_Test)
{
  variable var1;
  variable var2(2.);
  variable var3;

  // value just after constructor
  BOOST_CHECK_EQUAL( var1.value(), 0. );
  BOOST_CHECK_EQUAL( var2.value(), 2. );
  BOOST_CHECK_EQUAL( var3.value(), 0. );

  // set_value_ptr
  double x3 = 5.;
  var3.set_value_ptr(&x3);
  BOOST_CHECK_EQUAL( var3.value(), x3 );
  x3 = 6.;
  BOOST_CHECK_EQUAL( var3.value(), x3 );

  // bound_type after constructor
  BOOST_CHECK_EQUAL( var1.bound_type(), UNBOUNDED );

  // set_lower_bound
  BOOST_CHECK_EQUAL( var2.bound_type(), UNBOUNDED );
  BOOST_CHECK( var2.value() == var2.value_ubd() );
  var2.set_lower_bound(0.);
  BOOST_CHECK_EQUAL( var2.lower_bound(), 0. );
  BOOST_CHECK_EQUAL( var2.bound_type(), LOWER_BOUND );
  BOOST_CHECK( var2.value() != var2.value_ubd() );

  // set_upper_bound
  var2.set_upper_bound(4.);
  BOOST_CHECK_EQUAL( var2.upper_bound(), 4. );
  BOOST_CHECK_EQUAL( var2.bound_type(), BOUNDED );
  BOOST_CHECK_EQUAL( var3.bound_type(), UNBOUNDED );
  BOOST_CHECK( var3.value() == var3.value_ubd() );
  var3.set_upper_bound(10.);
  BOOST_CHECK_EQUAL( var3.upper_bound(), 10. );
  BOOST_CHECK_EQUAL( var3.bound_type(), UPPER_BOUND );
  BOOST_CHECK_EQUAL( var3.value(), x3 );
  BOOST_CHECK( var3.value() != var3.value_ubd() );

  // check throw when setting improper bounds
  BOOST_CHECK_THROW( var1.set_bound_type((VarBd)5), std::runtime_error );
  BOOST_CHECK_EQUAL( var1.bound_type(), UNBOUNDED );
  BOOST_CHECK_THROW( var3.set_lower_bound(12.), std::runtime_error );
  BOOST_CHECK_EQUAL( var3.bound_type(), UPPER_BOUND );
  var3.set_lower_bound(0.);
  BOOST_CHECK_EQUAL( var3.bound_type(), BOUNDED );
  BOOST_CHECK_THROW( var3.set_upper_bound(-2.), std::runtime_error );
  BOOST_CHECK_EQUAL( var3.upper_bound(), 10. );

  // setting/getting value for an unbounded variable
  var1.set_value(123.);
  BOOST_CHECK_EQUAL( var1.value(), 123. );
  BOOST_CHECK_EQUAL( var1.value_ubd(), 123. );
  var1.set_value_ubd(146.);
  BOOST_CHECK_EQUAL( var1.value(), 146. );
  BOOST_CHECK_EQUAL( var1.value_ubd(), 146. );
}

BOOST_AUTO_TEST_CASE(Function_Test)
{
  size_t n = 3;
  boost::scoped_array<double> x(new double[n]);
  for (size_t i = 0; i < n; i++)
    x[i] = 2.;
  double x0[3] = {0.5, 2., -1.};
  double c[3]  = {1., 10., 4.};

  // Show the various ways to construct objective functions
  // f0 is an empty function
  function f0(n);
  // f1 is built from the derived class quad_func from function
  quad_func f1(n, (double*)x0, (double*)c);
  // f2 is pointing to a C-like function address
  function f2(&quad_xyz, n);
  // f3 is pointing to a static wrapper function and its class instance
  function f3(&quad_func::evaluate, &f1, n);
  // f4 is linked to a boost::function which is binding a class method to its 
  // object
  boost::function<double (double * const)> feval;
  feval = boost::bind(&quad_func::eval, &f1, _1);
  function f4(feval, n);

  // get_n just after constructor
  BOOST_CHECK_EQUAL( f0.get_n(), n );
  BOOST_CHECK_EQUAL( f1.get_n(), n );
  BOOST_CHECK_EQUAL( f2.get_n(), n );
  BOOST_CHECK_EQUAL( f3.get_n(), n );
  BOOST_CHECK_EQUAL( f4.get_n(), n );

  // access variable and check integrity of variables
  BOOST_CHECK_THROW( f0.get_var(n), std::out_of_range );
  for (size_t i = 0; i < n; i++) {
    BOOST_CHECK_SMALL( f0.get_var(i).value(), 0.000001 );
    BOOST_CHECK_EQUAL( f0.get_var(i).bound_type(), UNBOUNDED );
    BOOST_CHECK_EQUAL( f0.get_var(i).value(), f0.get_var(i).value_ubd() );
  }

  // change number of variables
  f0.set_n(10);
  BOOST_CHECK_EQUAL( f0.get_n(), 10 );
  BOOST_CHECK_SMALL( f0.get_var(9).value(), 0.000001 );
  BOOST_CHECK_THROW( f0.get_var(10), std::out_of_range );

  // evaluate the function
  double val = 0.;
  for (size_t i = 0; i < n; i++) 
    val += c[i] * std::pow(x[i] - x0[i], 2);
  BOOST_CHECK_THROW( f0.eval(x.get()), std::runtime_error );
  BOOST_CHECK_CLOSE( f1.eval(x.get()), val, 0.0001 );
  BOOST_CHECK_EQUAL( f2.eval(x), quad_xyz(x.get()) );
  BOOST_CHECK_EQUAL( f3.eval(x), f1.eval(x.get()) );
  BOOST_CHECK_EQUAL( f4.eval(x), f1.eval(x.get()) );
  BOOST_CHECK_THROW( f0.eval_ubd(x.get()), std::runtime_error );
  BOOST_CHECK_CLOSE( f1.eval_ubd(x), val, 0.0001 );
  BOOST_CHECK_EQUAL( f2.eval_ubd(x), quad_xyz(x.get()) );
  BOOST_CHECK_EQUAL( f3.eval_ubd(x), f1.eval(x.get()) );
  BOOST_CHECK_EQUAL( f4.eval_ubd(x), f1.eval(x.get()) );


  // working with vector
  vector<double> vx(n);
  for (size_t i = 0; i < n; i++) 
    vx[i] = x[i];
  BOOST_CHECK_THROW( f0.eval(vx), std::runtime_error );
  BOOST_CHECK_EQUAL( f2.eval(vx), quad_xyz(x.get()) );
  BOOST_CHECK_EQUAL( f3.eval(vx), f1.eval(x.get()) );
  BOOST_CHECK_EQUAL( f4.eval(vx), f1.eval(x.get()) );

  // Evaluate the gradient
  vector<double> grad_xyz(n);
  vector<double> grad_f1(n);
  quad_xyz_gradient(vx, grad_xyz); // Exact gradient
  f1.grad(vx, grad_f1);            // Exact gradient
  vector<double> grad1(n);
  vector<double> grad2(n);
  vector<double> grad3(n);
  vector<double> grad4(n);
  BOOST_CHECK_THROW( f0.get_approx_gradient(vx, grad1), std::runtime_error );
  f1.get_approx_gradient(vx, grad1);
  f2.get_approx_gradient(vx, grad2);
  f3.get_approx_gradient(vx, grad3);
  f4.get_approx_gradient(vx, grad4);
  for (size_t i = 0; i < n; i++) {
    BOOST_CHECK_CLOSE( grad1[i], grad_f1[i], 0.0001 );
    BOOST_CHECK_CLOSE( grad2[i], grad_xyz[i], 0.0001 );
    BOOST_CHECK_CLOSE( grad3[i], grad_f1[i], 0.0001 );
    BOOST_CHECK_CLOSE( grad4[i], grad_f1[i], 0.0001 );
  }

  // set all values simultaneously
  vx[0] = x[0] = 135.; vx[1] = x[1] = 124.; vx[2] = x[2] = 89.;
  f1.set_values(x.get()); // pointer
  f2.set_values(x);       // boost::scoped_array
  f3.set_values(vx);      // boost::numeric::ublas::vector
  for (size_t i = 0; i < n; i++) {
    BOOST_CHECK_EQUAL( f1.get_var(i).value(), x[i] );
    BOOST_CHECK_EQUAL( f2.get_var(i).value(), x[i] );
    BOOST_CHECK_EQUAL( f3.get_var(i).value(), x[i] );
  }
}

// EOF
