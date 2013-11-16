/*!\file   minimizer.hpp
 * \brief  Abstract base class for minimization algorithms
 * \author Renaud Bruneliere
 * \date   01.10.2013
 *
 * The header contains the abstract base class for minimization algorithms,
 * and typedef for array of parameters and function pointers.
 */

#ifndef OPTIMIZE_MINIMIZER_HPP
#define OPTIMIZE_MINIMIZER_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include "optimize/function.hpp"
#include "optimize/variable.hpp"

namespace optimize {

//! Base class for minimization algorithms
class minimizer
{
 public:
  //! Constructor with number of variables to minimize.
  /*!
    \param n  number of variables
    \param x0 starting coordinate values
  */
  minimizer(size_t n, double * x0 = NULL);

  //! Constructor with function and starting coordinate values.
  /*!
    \param fptr pointer to an optimize::function
    \param func C-like function pointer or boost::function
    \param n    number of variables
    \param x0   starting coordinate values
    \param obj  pointer to a class instance passed to func_obj
  */
  // @{
  minimizer(function * const fptr, double * x0 = NULL);
  minimizer(const t_func & func, size_t n, double * x0 = NULL);
  minimizer(t_boofunc & func, size_t n, double * x0 = NULL);
  minimizer(const t_func_obj & func, void * obj, size_t n, double * x0 = NULL);
  // @}

  //! Destructor
  virtual ~minimizer() {}

  //! Minimize using function previously set.
  virtual int minimize(double * x0 = NULL) = 0;

  //! Minimize specifying a new function.
  /*!
    \param func new function pointer (to minimize)
    \param obj  object pointer
  */
  // @{
  int minimize(const t_func & func, double * x0 = NULL);
  int minimize(t_boofunc & func, double * x0 = NULL);
  int minimize(const t_func_obj & func, void * obj, double * x0 = NULL);
  // @}

  //! Get access to the objective function pointer
  function * const get_function() const { return fptr_; }
  //! Get access to a variable at the local minimum
  variable& get_opt_var(size_t i);
  //! Get a pointer to the array of variables at the local minimum
  variable * const get_opt_vars() const { return vopt_.get(); }

 protected:
  void store_opt_values(double * const xopt);

 private:
  function * fptr_;                  ///< Objective function pointer
  boost::scoped_ptr<function> func_; ///< Objective function (optional)
  t_vars vopt_; ///< Coordinates at the local minimum

  // Private methods
  void init(double * x0 = NULL);
};

} // namespace optimize

#endif // OPTIMIZE_MINIMIZER_HPP
