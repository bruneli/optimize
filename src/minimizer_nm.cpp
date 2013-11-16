#include "optimize/minimizer_nm.hpp"
#include <boost/thread.hpp>

optimize::minimizer_nm::minimizer_nm(size_t n, double * x0, bool is_mt) :
  optimize::minimizer(n, x0), is_mt_(is_mt)
{
  optimize::minimizer_nm::init();
}

optimize::minimizer_nm::minimizer_nm(optimize::function * const fptr, 
                                     double * x0, bool is_mt) :
  optimize::minimizer(fptr, x0), is_mt_(is_mt)
{
  optimize::minimizer_nm::init();
}

optimize::minimizer_nm::minimizer_nm(const optimize::t_func & func, 
                                     size_t n, double * x0, bool is_mt) :
  optimize::minimizer(func, n, x0), is_mt_(is_mt)
{
  optimize::minimizer_nm::init();
}

optimize::minimizer_nm::minimizer_nm(optimize::t_boofunc & func, 
                                     size_t n, double * x0, bool is_mt) :
  optimize::minimizer(func, n, x0), is_mt_(is_mt)
{
  optimize::minimizer_nm::init();
}

optimize::minimizer_nm::minimizer_nm(const optimize::t_func_obj & func,
                                     void * obj, size_t n, double * x0, 
                                     bool is_mt) :
  optimize::minimizer(func, obj, n, x0), is_mt_(is_mt)
{
  optimize::minimizer_nm::init();
}

int optimize::minimizer_nm::minimize(double * x0)
{
  size_t n = this->get_function()->get_n();

  // Update starting coordinates if specified
  if (x0) this->get_function()->set_values(x0);

  // Generate starting simplex, compute values of func at the vertices,
  // find min and max indices, find the barycenter
  generate_starting_simplex(is_mt_);

  boost::scoped_array<double> xnew(new double[n]);
  double fnew;
  boost::scoped_array<double> d(new double[n]);
  double norm_sq;

  // Main loop
  unsigned int iter = 0;
  for ( ; iter < max_iter_; iter++) {
    // Get the move vector and check for convergence
    get_move_vector(d, norm_sq);
    if (norm_sq/n < tol_*tol_) break;

    // Try reflexion : xnew = x_[imax_] + c_reflexion_*d
    move(x_[imax_], d, c_reflexion_, xnew, fnew);
    if (fnew <= f_[imin_]) { // Accept reflexion
      copy(xnew, fnew, imax_);
      // Try expansion : xnew = x_[imax_] + c_expansion_*d
      move(x_[imax_], d, c_expansion_, xnew, fnew);
      if (fnew <= f_[imin_]) // Accept expansion
        copy(xnew, fnew, imax_);
      find_min_max();
    } else {
      if (fnew <= f_[imax_]) { // Accept reflexion
        copy(xnew, fnew, imax_);
        find_min_max();
      } else {
        // Try contraction : xnew = x_[imax_] + c_contraction_*d
        move(x_[imax_], d, c_contraction_, xnew, fnew);
        if (fnew <= f_[imax_]) { // Accept contraction
          copy(xnew, fnew, imax_);
          find_min_max();
        } else { // Shrink the simplex
          if (is_mt_)
            shrink_mt();
          else
            shrink();
        }         
      }
    }
  } // End loop over iterations

  // Store minimum coordinates into vopt_
  optimize::minimizer::store_opt_values(x_[imin_].get());

  int ierr = (iter + 1 >= max_iter_);   
  return ierr;
}

void optimize::minimizer_nm::init()
{
  // Retrieve number of parameters
  size_t n = this->get_function()->get_n();

  // Set default values
  c_reflexion_ = 2.;
  c_expansion_ = 1.;
  c_contraction_ = c_shrinkage_ = 0.5;
  tol_ = 1.e-5;
  max_iter_ = n * 200;
  rel_delta_ = 0.05;
  abs_delta_ = 0.00025;

  // Set proper array size
  x_.reset(new boost::scoped_array<double>[n + 1]);
  f_.reset(new double[n + 1]);
  barycenter_.reset(new double[n]);
  for (size_t i = 0; i < n + 1; i++)
    x_[i].reset(new double[n]);
  params_step_.reset(new double[n]);
  for (size_t i = 0; i < n; i++)
    params_step_[i] = 0.;
}

void optimize::minimizer_nm::copy(const boost::scoped_array<double> & x,
                                  double f, size_t idx)
{
  size_t n = this->get_function()->get_n();
  for (size_t i = 0; i < n; i++) {
    barycenter_[i] += (x[i] - x_[idx][i])/n;
    x_[idx][i] = x[i];
  }
  f_[idx] = f;
}

void optimize::minimizer_nm::find_min_max()
{
  for (size_t ivtx = 0; ivtx < this->get_function()->get_n() + 1; ivtx++) {
    if (ivtx) {
      if (f_[ivtx] > f_[imax_]) imax_ = ivtx;
      if (f_[ivtx] < f_[imin_]) imin_ = ivtx;
    } else
      imax_ = imin_ = 0;
  }
}

void optimize::minimizer_nm::generate_starting_simplex(bool is_mt)
{
  size_t n = this->get_function()->get_n();

  boost::thread_group g;
  for (size_t ivtx = 0; ivtx < n + 1; ivtx++) {
    if (is_mt)
      g.create_thread(
        boost::bind(&optimize::minimizer_nm::generate_starting_vtx, 
                    this, ivtx));
    else
      generate_starting_vtx(ivtx);
  }
  if (is_mt) g.join_all();

  double wgt = 1./n;
  for (size_t ivtx = 0; ivtx < n + 1; ivtx++) {
    // Evaluate the barycenter
    for (size_t i = 0; i < n; i++) {
      if (ivtx)
        barycenter_[i] += wgt*(x_[ivtx][i]);
      else
        barycenter_[i] = wgt*(x_[ivtx][i]);
    }
    // Find the min and max indices
    if (ivtx) {
      if (f_[ivtx] > f_[imax_]) imax_ = ivtx;
      if (f_[ivtx] < f_[imin_]) imin_ = ivtx;
    } else
      imax_ = imin_ = 0;
  }
}

void optimize::minimizer_nm::generate_starting_vtx(size_t ivtx)
{
  // Define starting vertex
  for (size_t i = 0; i < this->get_function()->get_n(); i++) {
    x_[ivtx][i] = this->get_function()->get_var(i).value_ubd();
    if (ivtx && i == ivtx - 1) {
      if (params_step_[i] != 0.)
        x_[ivtx][i] += params_step_[i];
      else if (x_[ivtx][i] == 0.)
        x_[ivtx][i] = abs_delta_;
      else
        x_[ivtx][i] += rel_delta_*(x_[ivtx][i]);
    }
  }
  // Evaluate function at vertex
  f_[ivtx] = this->get_function()->eval_ubd(x_[ivtx]);
}

void optimize::minimizer_nm::get_move_vector(boost::scoped_array<double> & d,
                                             double & norm_sq)
{
  size_t n = this->get_function()->get_n();
  double wgt = 1. + 1./n;
  norm_sq = 0.;
  for (size_t i = 0; i < n; i++) {
    d[i] = barycenter_[i] - wgt*(x_[imax_][i]);
    norm_sq += d[i]*d[i];
  }
}

void optimize::minimizer_nm::move(const boost::scoped_array<double> & x0,
                                  const boost::scoped_array<double> & d, 
                                  double coef,
                                  boost::scoped_array<double> & xnew,
                                  double & fnew)
{
  for (size_t i = 0; i < this->get_function()->get_n(); i++)
    xnew[i] = x0[i] + coef*d[i];
  fnew = this->get_function()->eval_ubd(xnew);
}

void optimize::minimizer_nm::move_mt(double * xold, double * d, 
                                     double coef,
                                     double * xnew, double & fnew)
{
  try {
    for (size_t i = 0; i < this->get_function()->get_n(); i++)
      xnew[i] = xold[i] + coef*d[i];
    fnew = this->get_function()->eval_ubd(xnew);
    boost::this_thread::interruption_point();
  } catch(boost::thread_interrupted const&) {} // catch and ignore interruption
}

void optimize::minimizer_nm::shrink()
{
  size_t n = this->get_function()->get_n();
  size_t imin0 = imin_;
  imin_ = imax_ = imin0;
  double wgt = 1./n;
  for (size_t i = 0; i < n; i++)
    barycenter_[i] = wgt*(x_[imin0][i]);
  for (size_t ivtx = 0; ivtx < n + 1; ivtx++) {
    if (ivtx == imin0) continue;
    // Create a new vertex closer to current minimum
    for (size_t i = 0; i < n; i++) {
      x_[ivtx][i] = x_[imin0][i] + c_shrinkage_*(x_[ivtx][i] - x_[imin0][i]);
      barycenter_[i] += wgt*(x_[ivtx][i]);
    }
    // Evaluate function at vertex
    f_[ivtx] = this->get_function()->eval_ubd(x_[ivtx]);
    // Find min and max indices
    if (f_[ivtx] > f_[imax_]) imax_ = ivtx;
    if (f_[ivtx] < f_[imin_]) imin_ = ivtx;
  }
}

void optimize::minimizer_nm::shrink_mt()
{
  size_t n = this->get_function()->get_n();
  // Parallel computing of new vertices and function values
  boost::thread_group g;
  for (size_t ivtx = 0; ivtx < n + 1; ivtx++) {
    if (ivtx == imin_) continue;
    g.create_thread(boost::bind(&optimize::minimizer_nm::shrink_vtx, 
                                this, ivtx));
  }
  g.join_all();
  // Find new min, max, barycenter
  // Assume it is quicker to do that outside threads (avoids locks)
  size_t imin0 = imin_;
  imin_ = imax_ = imin0;
  double wgt = 1./n;
  for (size_t i = 0; i < n; i++)
    barycenter_[i] = wgt*(x_[imin0][i]);
  for (size_t ivtx = 0; ivtx < n + 1; ivtx++) {
    if (ivtx == imin0) continue;
    for (size_t i = 0; i < n; i++)
      barycenter_[i] += wgt*(x_[ivtx][i]);
    // Find min and max indices
    if (f_[ivtx] > f_[imax_]) imax_ = ivtx;
    if (f_[ivtx] < f_[imin_]) imin_ = ivtx;
  }
}

void optimize::minimizer_nm::shrink_vtx(size_t ivtx)
{
  // Create a new vertex closer to current minimum
  for (size_t i = 0; i < this->get_function()->get_n(); i++)
    x_[ivtx][i] = x_[imin_][i] + c_shrinkage_*(x_[ivtx][i] - x_[imin_][i]);
  // Evaluate function at vertex
  f_[ivtx] = this->get_function()->eval_ubd(x_[ivtx]);
}
