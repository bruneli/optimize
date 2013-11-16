/*!\file   function.hpp
 * \brief  Define a class to represent the objective function.
 * \author Renaud Bruneliere
 * \date   13.11.2013
 *
 * The header contains the base class used to define a real-valued objective
 * function.
 */

#ifndef OPTIMIZE_FUNCTION_HPP
#define OPTIMIZE_FUNCTION_HPP

#include <boost/function.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "optimize/variable.hpp"

namespace optimize {

//! t_func = C-like function pointer with a double array as argument
typedef double (*t_func)(double * const);
//! t_func_obj = function pointer with a double array and a void* as arguments
/*
  This typedef can be used to write a static wrapper function to callback
  a non-static class method via a pointer to the object.
*/
typedef double (*t_func_obj)(double * const, void *);
//! t_boofunc = boost::function<double (double * const)> typedef
typedef boost::function<double (double * const)> t_boofunc;

//! t_vars = dynamic array of function variables
typedef boost::scoped_array<variable> t_vars;

//! Base class defining a real-valued objective function
class function
{
 public:
  //! Default constructor.
  function(size_t n = 0);

  //! Constructor with function (and possibly the number of variables) specified
  /*!
    \param func     C-like function pointer
    \param boofunc  boost::function
    \param func_obj static (wrapper) function pointer
    \param obj      pointer to a class instance passed to func_obj
    \param n        number of variables
  */
  // @{
  function(const t_func & func, size_t n = 0);
  function(t_boofunc & boofunc, size_t n = 0);
  function(const t_func_obj & func_obj, void * const obj, size_t n = 0);
  // @}

  //! Destructor
  virtual ~function() {}

  //! Evaluate the function for a given set of variables
  /*!
    This method can be overloaded.
  */
  // @{
  virtual double eval(double * const x) const;
  double eval(const boost::scoped_array<double>& x) const;
  double eval(const boost::numeric::ublas::vector<double>& x) const;
  // @}

  //! Evaluate the function for a given set of unbounded variables
  // @{
  double eval_ubd(double * const x) const;
  double eval_ubd(const boost::scoped_array<double>& x) const;
  double eval_ubd(const boost::numeric::ublas::vector<double>& x) const;
  // @}

  //! Estimate via finite difference the gradient grad of f in x
  // @{
  void get_approx_gradient(const boost::numeric::ublas::vector<double>& x,
                           boost::numeric::ublas::vector<double>& grad,
                           bool check_bounds = false) const;
  void get_approx_gradient(const boost::numeric::ublas::vector<double>& x0,
                           double f0,
                           boost::numeric::ublas::vector<double>& grad,
                           bool check_bounds = false) const;
  // @}

  //! Get the number of variables
  size_t get_n() const { return n_; }
  //! Set the number of variables (be careful, it resets all variable settings)
  void set_n(size_t n);
  //! Get access to a variable
  variable& get_var(size_t i) const;
  //! Get a pointer to the array of variables
  variable * const get_vars() const { return vars_.get(); }
  //! Set all variables values
  // @{
  void set_values(double * const x);
  void set_values(const boost::scoped_array<double>& x);
  void set_values(const boost::numeric::ublas::vector<double>& x);
  // @}
  //! Get the step length
  double get_eps() const { return eps_; }
  //! Set the step length
  void set_eps(double eps) { eps_ = eps; }
  //! Set the C-like function pointer
  void set_func_ptr(const t_func & func) { func_ = &func; }
  //! Set the boost function
  void set_boost_func(t_boofunc & boofunc) { boofunc_ = &boofunc; }
  //! Set the static function pointer
  void set_func_obj_ptr(const t_func_obj & func_obj) { func_obj_ = &func_obj; }
  //! Set the object pointer called by func_obj
  void set_obj_ptr(void * obj) { obj_ = obj; }
  //! Returns bound flag telling if any variable is constrained via bounds
  bool is_bounded() const { return is_bounded_; }
  //! Set bound flag
  void set_is_bounded(bool is_bounded) { is_bounded_ = is_bounded; }

 private:
  size_t n_;    ///< Number of variables
  t_vars vars_; ///< Array of variables
  boost::scoped_array<double> vx_; ///< Array of variables values
  double eps_; ///< Step length used to compute finite-difference derivatives
  // Three different types of functions can be minimized
  const t_func * func_;         ///< Pointer to the C-like function pointer
  t_boofunc * boofunc_;         ///< Pointer to the boost::function
  const t_func_obj * func_obj_; ///< Pointer to the function pointer + obj
  void * obj_;                  ///< Pointer to an object
  bool is_bounded_; ///< Internal check of parameter bounds

  // Private methods
  void init();
};

} // namespace optimize

#endif // OPTIMIZE_FUNCTION_HPP
