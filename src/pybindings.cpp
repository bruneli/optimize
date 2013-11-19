#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/range.hpp>
#include <boost/foreach.hpp>

#include "optimize/pybindings.hpp"
#include "optimize/minimizer_nm.hpp"
#include "optimize/minimizer_bfgs.hpp"

namespace bp = boost::python;

boost::python::tuple optimize::minimize(PyObject * cb, bp::tuple& x0,
                                        std::string method)
{
  optimize::pyf pyf(cb, bp::len(x0));
  boost::scoped_array<double> x(new double[pyf.get_n()]);
  for (size_t i = 0; i < pyf.get_n(); i++)
    x[i] = bp::extract<float>(x0[i]);
  boost::scoped_ptr< optimize::minimizer > pymin;
  if (method == "Nelder-Mead")
    pymin.reset(new optimize::minimizer_nm(&pyf, x.get()));
  else if (method == "BFGS")
    pymin.reset(new optimize::minimizer_bfgs(&pyf, x.get()));
  else {
    std::string msg = "method " + method + " is not recognized.";
    throw std::runtime_error(msg);
  }
  pymin->minimize();
  bp::list l;
  for (size_t i = 0; i < pyf.get_n(); i++)
    l.append(pymin->get_opt_var(i).value());
  return bp::tuple(l);
}

double optimize::pyf::eval(double * const x) const
{
  bp::list lx;
  BOOST_FOREACH(double& val, boost::make_iterator_range(x, x + this->get_n()))
    lx.append(val);
  return bp::call<double>(pycb_, lx);
}

/*
double optimize::pyf::eval(double * const x) const
{
  npy_intp size = this->get_n();
  PyObject * pyObj = PyArray_SimpleNewFromData( 1, &size, NPY_DOUBLE, x );
  boost::python::handle<> handle( pyObj );
  boost::python::numeric::array arr( handle );
  return bp::call<double>(pycb_, arr);
}
*/

//void optimize::translate(optimize::pyopt_exception const& e)
void optimize::translate(std::runtime_error const& e)
{
  PyErr_SetString(PyExc_RuntimeError, e.what());
}

// Expose classes and methods to Python
BOOST_PYTHON_MODULE(liboptimize) 
{
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  boost::python::register_exception_translator<std::runtime_error>(&optimize::translate);
  def("minimize", &optimize::minimize, 
      ( "fun", "x0", bp::arg("method")="Nelder-Mead") );
}
