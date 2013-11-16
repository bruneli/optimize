/*!\file   minimizer_bfgs.hpp
 * \brief  Broyden–Fletcher–Goldfarb–Shanno algorithm 
 * \author Renaud Bruneliere
 * \date   28.10.2013
 *
 * The class contains an implementation of the Broyden–Fletcher–Goldfarb–Shanno
 * algorithm, an approximate Newton's method, used to find a multi-dimensional
 * local minimum. The algorithm is similar to the gradient descent method with
 * an Hessian matrix approximated.
 */

#ifndef OPTIMIZE_MINIMIZER_BFGS_HPP
#define OPTIMIZE_MINIMIZER_BFGS_HPP

#include <boost/scoped_ptr.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include "optimize/minimizer.hpp"

namespace optimize {

//! Class implementing the Broyden–Fletcher–Goldfarb–Shanno algorithm
class minimizer_bfgs : public minimizer
{
 public:
  //! Constructor with number of parameters.
  /*!
    \param n  number of variables
    \param x0 starting coordinate values
  */
  minimizer_bfgs(size_t n, double * x0 = NULL);

  //! Constructor with function and starting coordinate values.
  /*!
    \param fptr pointer to an optimize::function 
    \param func C-like function pointer or boost::function
    \param n    number of variables
    \param p0   starting coordinate values
    \param obj  pointer to a class instance passed to func_obj
  */
  // @{
  minimizer_bfgs(function * const fptr, double * x0 = NULL);
  minimizer_bfgs(const t_func & func, size_t n, double * x0 = NULL);
  minimizer_bfgs(t_boofunc & func, size_t n, double * x0 = NULL);
  minimizer_bfgs(const t_func_obj & func, void * obj, size_t n,
                 double * x0 = NULL);
  // @}

  //! Destructor
  virtual ~minimizer_bfgs() {}

  // This is to avoid name hiding of all minimizer::minimize methods
  using minimizer::minimize;

  //! Minimize using parameters and function previously set.
  virtual int minimize(double * x0 = NULL);

  //! Get value of the function at minimum
  double get_fmin() const { return fmin_; }

  //! Get current approximated inverse Hessian matrix value
  double get_inv_hessian(size_t i, size_t j) const { return (*invh_)(i, j); }
  //! Get pointer to the approximated inverse Hessian matrix
  boost::numeric::ublas::matrix<double> * const get_inv_hessian() \
    { return invh_.get(); }

  //! Set maximum number of iterations
  void set_max_iterations(unsigned int max_iter) { max_iter_ = max_iter; }

 private:
  // Class members
  ///< Inverse of the Hessian matrix
  boost::scoped_ptr< boost::numeric::ublas::matrix<double> > invh_;
  double fmin_; ///< Value of the function at minimum

  // Members used to stop iterations (error tolerance, number of iterations)
  double tol_;            ///< Error tolerance
  unsigned int max_iter_; ///< Maximum number of iterations

  // Private methods
  void init();
};

} // namespace optimize

#endif // OPTIMIZE_MINIMIZER_BFGS_HPP
