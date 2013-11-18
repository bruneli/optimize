#include <iostream>

#include <boost/scoped_array.hpp>

#include "optimize/pybindings.hpp"
#include "optimize/minimizer_nm.hpp"

namespace bp = boost::python;

boost::python::tuple optimize::minimize(PyObject * cb, bp::tuple& x0)
{
  optimize::pyf pyf(cb, bp::len(x0));
  boost::scoped_array<double> x(new double[pyf.get_n()]);
  for (size_t i = 0; i < pyf.get_n(); i++)
    x[i] = bp::extract<float>(x0[i]);
  optimize::minimizer_nm pymin(&pyf, x.get());
  pymin.minimize();
  bp::list l;
  for (size_t i = 0; i < pyf.get_n(); i++)
    l.append(pymin.get_opt_var(i).value());
  return bp::tuple(l);
}

double optimize::pyf::eval(double * const x) const
{
  bp::list lx;
  for (size_t i = 0; i < this->get_n(); i++)
    lx.append(x[i]);
  return bp::call<double>(pycb_, lx);
}

// Expose classes and methods to Python
BOOST_PYTHON_MODULE(liboptimize) 
{
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  def("minimize", &optimize::minimize);
}
