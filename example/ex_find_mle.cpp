/*!\file   ex_find_mle.cpp
 * \brief  Example on how to find the maximum likelihood estimates
 * \author Renaud Bruneliere
 * \date   14.10.2013
 *
 * In this example, the Nelder-Mead method is used to find the minimum of the
 * log-likelihood function and then perform a likelihood ratio test.
 */

#include <iostream>
#include <cmath>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/scoped_array.hpp>
#include <boost/math/distributions/poisson.hpp>

#include "optimize/minimizer_bfgs.hpp"
#include "optimize/minimizer_nm.hpp"

//! Class defining the likelihood function
/*!
  The likelihood function is constructed as the product of two poisson
  distributions with mean mu_on and mu_off respectively.
 */
class likelihood
{
 public:
  //! Constructor
  likelihood(unsigned int n_on, unsigned int n_off,
             double mu, double s, double b, double tau) : 
    n_on_(n_on), n_off_(n_off), mu_(mu), s_(s), b_(b), tau_(tau) { }

  //! Evaluate the negative log-likelihood function value
  double eval(double mu_on, double mu_off) const
  {
    boost::math::poisson pmf_off(mu_off);
    boost::math::poisson pmf_on(mu_on);
    return -2. * log(pdf(pmf_on, n_on_) * pdf(pmf_off, n_off_));
  }

  //! non-static class function used with boost::function and boost::bind
  double eval_unconstrained(double * const params)
  {
    double mu_off = tau_ * std::fabs(params[1]);
    double mu_on  = std::fabs(params[0]) * s_ + std::fabs(params[1]);
    return eval(mu_on, mu_off);
  }

  //! static wrapper function to callback the member function
  static double eval_constrained(double * const params, void * obj)
  {
    likelihood * l = (likelihood*) obj;
    double mu_off = l->tau() * std::fabs(params[0]);
    double mu_on  = l->mu() * l->s() + std::fabs(params[0]);
    return l->eval(mu_on, mu_off);
  }

  //! Getters
  // {@
  unsigned int n_on()  const { return n_on_; }
  unsigned int n_off() const { return n_on_; }
  double mu() const  { return mu_; }
  double s() const   { return s_; }
  double b() const   { return b_; }
  double tau() const { return tau_; }
  // @}

 private:
  unsigned int n_on_;
  unsigned int n_off_;
  double mu_;
  double s_;
  double b_;
  double tau_;
};

int main(int argc, char** argv)
{
  double mu = 0.; // Test a background only hypothesis => mu=0
  double b  = 1.; // Expected number of background events in signal region

  // Define the likelihood function
  likelihood l(4, 5, mu, 3., b, 5.);

  // A first minimizer is used to find the unconstrained maximum likelihood
  // estimates, i.e. the simultaneous minumum of mu and b.
  boost::scoped_array<double> p0(new double[2]);
  p0[0] = mu;
  p0[1] = b;
  optimize::minimizer_bfgs unconstrained_mle(2, p0.get());
  //unconstrained_mle.get_function()->get_var(0).set_lower_bound(0.);
  //unconstrained_mle.get_function()->get_var(0).set_value(1.);
  //unconstrained_mle.get_function()->get_var(1).set_lower_bound(0.);
  // Minimize with the help of a boost::function
  boost::function<double (double * const)> feval;
  feval = boost::bind(&likelihood::eval_unconstrained, &l, _1);
  unconstrained_mle.minimize(feval);
  std::cout << "unconstrained mle fmin = " << unconstrained_mle.get_fmin()
            << ", mu = " << unconstrained_mle.get_opt_var(0).value()
            << ", b = " << unconstrained_mle.get_opt_var(1).value()
            << std::endl;

  // A second minimizer is used to find the constrained maximum likelihood
  // estimates, i.e. the minumum of b with mu fixed to 0.
  p0.reset(new double[1]);
  p0[0] = b;
  optimize::minimizer_bfgs constrained_mle(1, p0.get());
  //constrained_mle.get_function()->get_var(0).set_lower_bound(0.);
  // Minimize via a static wrapper function
  constrained_mle.minimize(&likelihood::eval_constrained, &l);
  std::cout << "constrained mle   fmin = " << constrained_mle.get_fmin()
            << ", b = " << constrained_mle.get_opt_var(0).value() << std::endl;

  // Evaluate log-likelihood ratio and extract Z-score assuming Wilks theorem
  double lambda = constrained_mle.get_fmin() - unconstrained_mle.get_fmin();
  double zscore = sqrt(lambda);
  std::cout << "Z-score = " << zscore << std::endl;
}
