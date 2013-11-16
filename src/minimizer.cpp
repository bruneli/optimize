#include <stdexcept>
#include <boost/throw_exception.hpp>
#include "optimize/minimizer.hpp"

optimize::minimizer::minimizer(size_t n, double * x0)
{
  func_.reset(new optimize::function(n));
  fptr_ = func_.get();
  optimize::minimizer::init(x0);
}

optimize::minimizer::minimizer(optimize::function * const fptr, double * x0)
{
  // Initialize func_ even if unused in case minimize(func) is called later
  func_.reset(new optimize::function());
  fptr_ = fptr;
  optimize::minimizer::init(x0);
}

optimize::minimizer::minimizer(const optimize::t_func & func, 
                               size_t n, double * x0)
{
  func_.reset(new optimize::function(func, n));
  fptr_ = func_.get();
  optimize::minimizer::init(x0);
}

optimize::minimizer::minimizer(optimize::t_boofunc & func,
                               size_t n, double * x0)
{
  func_.reset(new optimize::function(func, n));
  fptr_ = func_.get();
  optimize::minimizer::init(x0);
}

optimize::minimizer::minimizer(const optimize::t_func_obj & func,
                               void * obj, size_t n, double * x0)
{
  func_.reset(new optimize::function(func, obj, n));
  fptr_ = func_.get();
  optimize::minimizer::init(x0);
}

int optimize::minimizer::minimize(const optimize::t_func & func, double * x0)
{
  func_->set_func_ptr(func);
  if (fptr_ != func_.get()) fptr_ = func_.get();
  return this->minimize(x0);
}

int optimize::minimizer::minimize(optimize::t_boofunc & func, double * x0)
{
  func_->set_boost_func(func);
  if (fptr_ != func_.get()) fptr_ = func_.get();
  return this->minimize(x0);
}

int optimize::minimizer::minimize(const optimize::t_func_obj & func,
                                  void * obj, double * x0)
{
  func_->set_func_obj_ptr(func);
  func_->set_obj_ptr(obj);
  if (fptr_ != func_.get()) fptr_ = func_.get();
  return this->minimize(x0);
}

optimize::variable& optimize::minimizer::get_opt_var(size_t i)
{
  if (i >= fptr_->get_n())
    BOOST_THROW_EXCEPTION( std::out_of_range("i exceeds variables range") );
  return vopt_[i];
}

void optimize::minimizer::store_opt_values(double * const xopt)
{
  size_t n = fptr_->get_n();
  vopt_.reset(new optimize::variable[n]);
  for (size_t i = 0; i < n; i++) {
    vopt_[i].copy_bounds(fptr_->get_var(i));
    vopt_[i].set_value_ubd(xopt[i]);
    fptr_->get_var(i).set_value_ubd(xopt[i]);
  }
}

void optimize::minimizer::init(double * x0)
{
  // Set default values
  if (x0) fptr_->set_values(x0);
  vopt_.reset(new optimize::variable[fptr_->get_n()]);
}
