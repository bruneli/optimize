#include <stdexcept>
#include <boost/throw_exception.hpp>
#include "optimize/function.hpp"

using namespace boost::numeric::ublas;

optimize::function::function(size_t n) : n_(n), 
  func_(NULL), boofunc_(NULL), func_obj_(NULL), obj_(NULL)
{ 
  init();
}

optimize::function::function(const optimize::t_func & func, size_t n) : n_(n), 
  func_(&func), boofunc_(NULL), func_obj_(NULL), obj_(NULL)
{ 
  init();
}

optimize::function::function(optimize::t_boofunc & boofunc, size_t n) : 
  n_(n), func_(NULL), boofunc_(&boofunc), func_obj_(NULL), obj_(NULL)
{ 
  init();
}

optimize::function::function(const optimize::t_func_obj & func_obj,
                             void * const obj, size_t n) : n_(n), 
  func_(NULL), boofunc_(NULL), func_obj_(&func_obj), obj_(obj)
{ 
  init();
}

double optimize::function::eval(double * const x) const
{
  if (func_obj_ && obj_)
    return (*func_obj_)(x, obj_);
  else if (func_)
    return (*func_)(x);
  else if (boofunc_)
    return (*boofunc_)(x);
  BOOST_THROW_EXCEPTION( std::runtime_error("No function is found") );
  return 0.;
}

double optimize::function::eval(const boost::scoped_array<double>& x) const
{
  return this->eval(x.get());
}

double optimize::function::eval(const vector<double>& x) const
{
  return this->eval((double*)&x[0]);
}

double optimize::function::eval_ubd(double * const x) const
{
  if (!is_bounded_) return this->eval(x);

  // Apply a non-linear transformation if needed
  for (size_t i = 0; i < n_; i++)
    vars_[i].set_value_ubd(x[i]);
  return this->eval(vx_.get());
}

double optimize::function::eval_ubd(const boost::scoped_array<double>& x) const
{
  return this->eval_ubd(x.get());
}

double optimize::function::eval_ubd(const vector<double>& x) const
{
  return this->eval_ubd((double*)&x[0]);
}

void optimize::function::get_approx_gradient(const vector<double>& x,
                                             vector<double>& grad,
                                             bool check_bounds) const
{
  double f = (check_bounds) ? this->eval_ubd(x) : this->eval(x);
  get_approx_gradient(x, f, grad);
}

void optimize::function::get_approx_gradient(const vector<double>& x0,
                                             double f0,
                                             vector<double>& grad,
                                             bool check_bounds) const
{
  vector<double> x1(x0);
  double f1 = 0.;
  for (size_t i = 0; i < n_; i++) {
    x1[i] += eps_;
    f1 = (check_bounds) ? this->eval_ubd(x1) : this->eval(x1);
    grad[i] = (f1 - f0)/eps_;
    x1[i] = x0[i];
  }
}

void optimize::function::set_n(size_t n)
{
  n_ = n;
  init(); // Reset array size and variable settings
}

optimize::variable& optimize::function::get_var(size_t i) const
{
  if (i >= n_)
    BOOST_THROW_EXCEPTION( std::out_of_range("i exceeds variables range") );
  return vars_[i];
}

void optimize::function::set_values(double * const x)
{
  for (size_t i = 0; i < n_; i++)
    vars_[i].set_value(x[i]);
}

void optimize::function::set_values(const boost::scoped_array<double>& x)
{
  this->set_values(x.get());
}

void optimize::function::set_values(const vector<double>& x)
{
  this->set_values((double*)&x[0]);
}

void optimize::function::init()
{
  // Set proper array size
  vars_.reset(new optimize::variable[n_]);
  vx_.reset(new double[n_]);
  // Set default values
  for (unsigned int i = 0; i < n_; i++) {
    vars_[i].set_value_ptr(&(vx_[i]));
    vars_[i].set_function_ptr(this);
  }
  eps_ = 1.e-8;
  is_bounded_ = false;
}
