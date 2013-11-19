/*!\file   pybindings.hpp
 * \brief  List of functions which can be called in python
 * \author Renaud Bruneliere
 * \date   07.11.2013
 *
 * The header contains a list of functions which can be called from python to
 * perform a minimization with C++ classes.
 */

#ifndef OPTIMIZE_PYBINDINGS_HPP
#define OPTIMIZE_PYBINDINGS_HPP

#include <exception>
#include <string>
#include <boost/python.hpp>
#include "optimize/function.hpp"

namespace optimize {

//! Call a python function and minimize it
/*!
  \param cb callable python object
  \param x0 tuple containing starting coordinates
  \param method algorithm name: Nelder-Mead, BFGS
  \return tuple with coordinates of variables at their local minimum
*/
boost::python::tuple minimize(PyObject * cb,
                              boost::python::tuple& x0,
                              std::string method = "Nelder-Mead");

//! Wrapper to a real-valued objective function
class pyf : public function
{
 public:
  //! Constructor
  /*!
    \param cb callable python object
    \param n  number of variables
  */
  pyf(PyObject * cb, size_t n) : pycb_(cb), function(n) {}

  //! Destructor
  virtual ~pyf() {}

  //! Evaluate the function for a given set of variables
  virtual double eval(double * const x) const;

 private:
  PyObject* pycb_; ///< Callable python function
};

//! Translate the C++ exception to a Python exception
void translate(std::runtime_error const& e);

} // namespace optimize

#endif // OPTIMIZE_PYBINDINGS_HPP
