#include <limits>
#include <cmath>
#include <string>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include "optimize/variable.hpp"
#include "optimize/function.hpp"

optimize::variable::variable(double val) : val_(val), ptr_(&val_),
  bd_type_(optimize::UNBOUNDED),
  lower_bd_(std::numeric_limits<double>::min()),
  upper_bd_(std::numeric_limits<double>::max()),
  fptr_(NULL)
{
  bd_to_ubd_ = &optimize::variable::bd_to_ubd_none;
  ubd_to_bd_ = &optimize::variable::bd_to_ubd_none;
}

double optimize::variable::value_ubd() const
{
  return (*bd_to_ubd_)(*ptr_, lower_bd_, upper_bd_);
}

void optimize::variable::set_value_ubd(double val)
{
  *ptr_ = (*ubd_to_bd_)(val, lower_bd_, upper_bd_);
}

void optimize::variable::set_bound_type(VarBd type)
{
  switch (type) {
  case optimize::UNBOUNDED:
    bd_type_ = type;
    bd_to_ubd_ = &optimize::variable::bd_to_ubd_none;
    ubd_to_bd_ = &optimize::variable::bd_to_ubd_none;
    break;
  case optimize::LOWER_BOUND:
    bd_type_ = type;
    bd_to_ubd_ = &optimize::variable::bd_to_ubd_lower;
    ubd_to_bd_ = &optimize::variable::bd_to_ubd_lower;
    break;
  case optimize::UPPER_BOUND:
    bd_type_ = type;
    bd_to_ubd_ = &optimize::variable::bd_to_ubd_upper;
    ubd_to_bd_ = &optimize::variable::bd_to_ubd_upper;
    break;
  case optimize::BOUNDED:
    bd_type_ = type;
    bd_to_ubd_ = &optimize::variable::bd_to_ubd_both;
    ubd_to_bd_ = &optimize::variable::bd_to_ubd_both;
    break;
  default:
    BOOST_THROW_EXCEPTION( std::runtime_error("Bound type is not recognized") );
  }
  if (fptr_ && bd_type_ != optimize::UNBOUNDED) fptr_->set_is_bounded(true);
}

void optimize::variable::set_lower_bound(double lower_bd)
{
  if (bd_type_ > 1 && lower_bd > upper_bd_) {
    std::string msg = "Lower bound cannot be larger than the upper bound";
    BOOST_THROW_EXCEPTION( std::runtime_error(msg) );
    return;
  }
  lower_bd_ = lower_bd;
  if (bd_type_ == optimize::UNBOUNDED)
    set_bound_type(optimize::LOWER_BOUND);
  else if (bd_type_ == optimize::UPPER_BOUND)
    set_bound_type(optimize::BOUNDED);
}

void optimize::variable::set_upper_bound(double upper_bd)
{
  if ((bd_type_ == optimize::LOWER_BOUND || bd_type_ == optimize::BOUNDED) &&
      upper_bd < lower_bd_) {
    std::string msg = "Upper bound cannot be smaller than the lower bound";
    BOOST_THROW_EXCEPTION( std::runtime_error(msg) );
    return;
  }
  upper_bd_ = upper_bd;
  if (bd_type_ == optimize::UNBOUNDED)
    set_bound_type(optimize::UPPER_BOUND);
  else if (bd_type_ == optimize::LOWER_BOUND)
    set_bound_type(optimize::BOUNDED);
}

void optimize::variable::copy_bounds(const variable& var)
{
  // Calls set_lower_bound and set_upper_bound before set_bound_type
  set_lower_bound(var.lower_bound());
  set_upper_bound(var.upper_bound());
  set_bound_type(var.bound_type());
}

double optimize::variable::bd_to_ubd_none(double val, 
                                          double lower, double upper)
{
  return val;
}

double optimize::variable::bd_to_ubd_lower(double val, 
                                           double lower, double upper)
{
  return std::sqrt(std::pow(val - lower + 1., 2) - 1.);
}

double optimize::variable::ubd_to_bd_lower(double val, 
                                           double lower, double upper)
{
  return lower - 1. + std::sqrt(val*val + 1.);
}

double optimize::variable::bd_to_ubd_upper(double val, 
                                           double lower, double upper)
{
  return std::sqrt(std::pow(upper - val + 1., 2) - 1.);
}

double optimize::variable::ubd_to_bd_upper(double val, 
                                           double lower, double upper)
{
  return upper + 1. - std::sqrt(val*val + 1.);
}

double optimize::variable::bd_to_ubd_both(double val, 
                                          double lower, double upper)
{
  if (upper - lower == 0.) return upper;
  return std::asin(2.*(val - lower)/(upper - lower) - 1.);
}

double optimize::variable::ubd_to_bd_both(double val, 
                                          double lower, double upper)
{
  return lower + (upper - lower)/2.*(std::sin(val) + 1.);
}
