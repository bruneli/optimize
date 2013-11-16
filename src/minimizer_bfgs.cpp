#include "optimize/minimizer_bfgs.hpp"
#include "optimize/line_search.hpp"

using namespace boost::numeric::ublas;

optimize::minimizer_bfgs::minimizer_bfgs(size_t n, double * x0) :
  optimize::minimizer(n, x0)
{
  optimize::minimizer_bfgs::init();
}

optimize::minimizer_bfgs::minimizer_bfgs(optimize::function * const fptr, 
                                         double * x0) :
  optimize::minimizer(fptr, x0)
{
  optimize::minimizer_bfgs::init();
}

optimize::minimizer_bfgs::minimizer_bfgs(const optimize::t_func & func, 
                                         size_t n, double * x0) :
  optimize::minimizer(func, n, x0)
{
  optimize::minimizer_bfgs::init();
}

optimize::minimizer_bfgs::minimizer_bfgs(optimize::t_boofunc & func, 
                                         size_t n, double * x0) :
  optimize::minimizer(func, n, x0)
{
  optimize::minimizer_bfgs::init();
}

optimize::minimizer_bfgs::minimizer_bfgs(const optimize::t_func_obj & func,
                                         void * obj, size_t n, double * x0) :
  optimize::minimizer(func, obj, n, x0)
{
  optimize::minimizer_bfgs::init();
}

int optimize::minimizer_bfgs::minimize(double * x0)
{
  size_t n = this->get_function()->get_n();

  // Update starting coordinates if specified
  if (x0) this->get_function()->set_values(x0);

  // Declare objects
  vector<double> xk(n);
  vector<double> gradfk(n);
  vector<double> xkpp(n);
  vector<double> gradfkpp(n);
  vector<double> sk(n); // := xkpp - xk
  vector<double> yk(n); // := gradfkpp - gradfk
  for (size_t i = 0; i < n; i++)
    xk[i] = this->get_function()->get_var(i).value();
  this->get_function()->get_approx_gradient(xk, gradfk, true);
  double norm = norm_2(gradfk);
  vector<double> pk(n); // search direction
  double alpha_k = 1.;
  double fkpp = 0.;
  matrix<double> m1(n, n);
  matrix<double> m2(n, n);

  // Starting inverse hessian candidate
  *invh_ = identity_matrix<double>(n);

  // Main loop
  unsigned int iter = 0;
  while ((norm > tol_) && (iter < max_iter_)) {
    // (quasi-)Newton search direction
    pk = (-1.) * prod((*invh_), gradfk);

    // Get a step length along pk satisfying the wolfe conditions to define the
    // next position
    alpha_k = optimize::line_search_wolfe(this->get_function(), xk, pk, fkpp);
    xkpp = xk + alpha_k * pk;
    this->get_function()->get_approx_gradient(xkpp, fkpp, gradfkpp, true);
    norm = norm_2(gradfkpp);
    if (norm < tol_) break;

    // New vs old difference in position and gradient
    sk = xkpp - xk;
    yk = gradfkpp - gradfk;
    double rhok = inner_prod(yk, sk);
    rhok = (rhok != 0.) ? 1./rhok : 1000.;

    // Update approximate inverse Hessian via the DFP formula
    // Hkpp = (I - rhok*sk*yk^T) * Hk * (I - rhok*yk*sk^T) + rhok*sk*sk^T
    m1 = identity_matrix<double>(n);
    m1 -= rhok * outer_prod(sk, yk); // rank-1 update
    m2 = identity_matrix<double>(n);
    m2 -= rhok * outer_prod(yk, sk); // rank-1 update
    *invh_ = prod((*invh_), m2); // matrices product
    *invh_ = prod(m1, (*invh_)); // matrices product
    *invh_ += rhok * outer_prod(sk, sk); // rank-1 update
    xk     = xkpp;     // Update position
    gradfk = gradfkpp; // Update gradient
    iter++;
  } // End while 

  // Store minimum coordinates into vopt_
  optimize::minimizer::store_opt_values((double*)&xkpp[0]);
  fmin_ = fkpp;

  int ierr = (iter + 1 >= max_iter_);   
  return ierr;
}

void optimize::minimizer_bfgs::init()
{
  // Retrieve number of parameters
  size_t n = this->get_function()->get_n();

  // Set default values
  fmin_ = -999.;
  tol_  = 1.e-5;
  max_iter_ = n * 200;

  // Define matrix size
  invh_.reset(new matrix<double>(n, n));
}
