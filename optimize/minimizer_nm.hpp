/*!\file   minimizer_nm.hpp
 * \brief  Implementation of the Nelder-Mead method
 * \author Renaud Bruneliere
 * \date   01.10.2013
 *
 * The header contains an implementation of the Nelder-Mead method,
 * also named the downhill simplex method, to find a multi-dimensional local
 * minimum. The Nelder-Mead method is an iterative gradient-less minimization 
 * algorithm.
 */

#ifndef OPTIMIZE_MINIMIZER_NM_HPP
#define OPTIMIZE_MINIMIZER_NM_HPP

#include "optimize/minimizer.hpp"

namespace optimize {

//! Class implementing the Nelder-Mead method
class minimizer_nm : public minimizer
{
 public:
  //! Constructor with number of parameters.
  /*!
    \param n  number of variables
    \param x0 starting coordinate values
    \param is_mt turn on/off the multithreading flag
  */
  minimizer_nm(size_t n, double * x0 = NULL, bool is_mt = false);

  //! Constructor with function and starting values for the parameters.
  /*!
    \param fptr pointer to an optimize::function 
    \param func C-like function pointer or boost::function
    \param n    number of variables
    \param p0   starting coordinate values
    \param is_mt turn on/off the multithreading flag
    \param obj  pointer to a class instance passed to func_obj
  */
  // @{
  minimizer_nm(function * const fptr,
               double * x0 = NULL, bool is_mt = false);
  minimizer_nm(const t_func & func, size_t n,
               double * x0 = NULL, bool is_mt = false);
  minimizer_nm(t_boofunc & func, size_t n,
               double * x0 = NULL, bool is_mt = false);
  minimizer_nm(const t_func_obj & func, void * obj, size_t n,
               double * x0 = NULL, bool is_mt = false);
  // @}

  //! Destructor
  virtual ~minimizer_nm() {}

  // This is to avoid name hiding of all minimizer::minimize methods
  using minimizer::minimize;

  //! Minimize using parameters and function previously set.
  virtual int minimize(double * x0 = NULL);

  //! Get value of the function at minimum
  double get_fmin() const { return f_[imin_]; }

  //! Is multithreading used ?
  bool is_mt() const { return is_mt_; }
  //! Turn on/off multithreading flag
  void set_is_mt(bool is_mt) { is_mt_ = is_mt; }

  //! Set reflexion coefficient
  void set_reflexion_coef(double c_reflexion) { c_reflexion_ = c_reflexion; }
  //! Set expansion coefficient
  void set_expansion_coef(double c_expansion) { c_expansion_ = c_expansion; }
  //! Set contraction coefficient
  void set_contraction_coef(double c_contraction) { 
    c_contraction_ = c_contraction; 
  }
  //! Set shrinkage coefficient
  void set_shrinkage_coef(double c_shrinkage) { c_shrinkage_ = c_shrinkage; }

  //! Set tolerance error
  void set_tolerance_error(double tol) { tol_ = tol; }
  //! Set maximum number of iterations
  void set_max_iterations(unsigned int max_iter) { max_iter_ = max_iter; }

  //! Set parameter starting step
  void set_param_step(size_t i, double step) { params_step_[i] = step; }
  //! Set relative delta starting step
  void set_rel_delta(double rel_delta) { rel_delta_ = rel_delta; }
  //! Set absolute delta starting step
  void set_abs_delta(double abs_delta) { abs_delta_ = abs_delta; }

 private:
  // Class members
  bool is_mt_; ///< Use multithreading during the minimization process
  //! Simplex, array of n+1 vertices
  boost::scoped_array< boost::scoped_array<double> > x_;
  boost::scoped_array<double> f_;   ///< Array of n+1 function values
  size_t imin_; ///< Index of mimimum vertex
  size_t imax_; ///< Index of maximum vertex
  boost::scoped_array<double> barycenter_; ///< Simplex barycenter

  // Members used to define the four different simplex moves
  double c_reflexion_;   ///< Reflexion move coefficient (2 by default)
  double c_expansion_;   ///< Expansion move coefficient (3 by default)
  double c_contraction_; ///< Contraction move coefficient (0.5 by default)
  double c_shrinkage_;   ///< Shrinkage move coefficient (0.5 by default)

  // Members used to stop iterations (error tolerance, number of iterations)
  double tol_;            ///< Error tolerance
  unsigned int max_iter_; ///< Maximum number of iterations

  // Members used to build the starting simplex
  /// Array of step size used to construct the starting simplex
  boost::scoped_array<double> params_step_;
  double rel_delta_; ///< Simplex vtx x[k] = (1 + rel_delta)*x[0] (x[0] != 0.)
  double abs_delta_; ///< Simplex vtx x[x] = abs_delta (x[0] == 0.)

  // Private methods
  void init();
  void copy(const boost::scoped_array<double> & x, double f, size_t idx);
  void find_min_max();
  void generate_starting_simplex(bool is_mt = false);
  void generate_starting_vtx(size_t ivtx);
  void get_move_vector(boost::scoped_array<double> & d, double & norm);
  void move(const boost::scoped_array<double> & x0,
            const boost::scoped_array<double> & d, double coef,
            boost::scoped_array<double> & xnew, double & fnew);
  void move_mt(double * xold, double * d, double coef,
               double * xnew, double & fnew);
  void shrink();
  void shrink_mt();
  void shrink_vtx(size_t ivtx);
};

} // namespace optimize

#endif // OPTIMIZE_MINIMIZER_NM_HPP
