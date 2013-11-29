#include <iostream>
#include <cmath>

// The boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/timer.hpp>
#include <boost/thread.hpp>
#include <boost/scoped_array.hpp>
#include <boost/multi_array.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>

#include "optimize/minimizer_nm.hpp"
#include "optimize/minimizer_bfgs.hpp"

//! Class defining the likelihood function
/*!
  The likelihood function is constructed as the product of multiple poisson
  and normal distributions.
 */
class likelihood
{
 public:
  //! Constructor
  likelihood(unsigned int nr, unsigned int np, unsigned int ns);

  ~likelihood() { }

  //! Compute number of expected events for a region
  double nexp(unsigned int ir, double * const params) const;
  unsigned int nobs(unsigned int ir) const { return nobs_[ir]; }

  //! non-static class function used with boost::function and boost::bind
  double eval(double * const params) const;

  unsigned int nr() const { return nr_; } ///< Number of regions
  unsigned int np() const { return np_; } ///< Number of processes
  unsigned int ns() const { return ns_; } ///< Number of systematics

  void set_nobs(unsigned int ir, unsigned int val) { nobs_[ir] = val; }
  void set_nexp(unsigned int ir, unsigned int ip, double val) 
  { 
    nexp_[ir][ip] = val;
  }
  void set_dsyst(unsigned int ir, unsigned int ip, unsigned int is, double val) 
  { 
    dsyst_[ir][ip][is] = val;
  }

 private:
  unsigned int nr_;
  unsigned int np_;
  unsigned int ns_;
  double mu_;
  boost::scoped_array<unsigned int> nobs_;
  boost::multi_array<double, 2> nexp_;
  boost::multi_array<double, 3> dsyst_;
};

likelihood::likelihood(unsigned int nr, unsigned int np, unsigned int ns) :
   nr_(nr), np_(np), ns_(ns)
{
  nobs_.reset(new unsigned int[nr_]);
  nexp_.resize(boost::extents[nr_][np_]);
  dsyst_.resize(boost::extents[nr_][np_][ns_]);
}

double likelihood::nexp(unsigned int ir, double * const params) const
{
  double nexp = 0.;
  for (unsigned int ip = 0; ip < np_; ip++) {
    double delta = 1.;
    for (unsigned int is = 0; is < ns_; is++)
      delta += dsyst_[ir][ip][is] * params[np_+is];
    nexp += params[ip] * nexp_[ir][ip] * delta;
  }
  //return (nexp > 0.001) ? nexp : 0.001;
  return (std::fabs(nexp) > 0.000001) ? std::fabs(nexp) : 0.000001;
}

double likelihood::eval(double * const params) const
{
  double lf = 1.;
  for (unsigned int ir = 0; ir < nr_; ir++) {
    double nexp = this->nexp(ir, params);
    if (nexp > 100. * nobs_[ir]) nexp = 100. * nobs_[ir];
    boost::math::poisson pmf(nexp);
    lf *= pdf(pmf, nobs_[ir]);
  }
  for (unsigned int is = 0; is < ns_; is++) {
    double loc = params[np_+is];
    //if (isnan(loc)) loc = 0.;
    if (loc > 20.) loc = 20.;
    if (loc < -20.) loc = -20.;
    boost::math::normal pdf_syst(loc);
    lf *= pdf(pdf_syst, 0.);
  }
  //for (unsigned int i = 0; i < 100000; i++)
  //  double x = std::sqrt(i + 1.) / (i + 1.);
  return -2. * log(lf);
}

int main(int argc, char* argv[])
{
  // Number of threads
  unsigned int n_threads = boost::thread::hardware_concurrency();
  std::cout << "Number of available threads " << n_threads << std::endl;

  likelihood l(5, 5, 9);
  double nexp[5][5] = {{1.,  1.6,  0.6, 0.2,  0.0},
                       {0., 15.4,  2.3, 0.0,  0.0},
                       {0.,  2.5,  4.1, 0.0,  0.0},
                       {0.,  3.3,  0.0, 1.1,  0.0},
                       {0.,  9.5, 10.0, 0.7, 39.2}};
  std::cout << "n(regions) = " << l.nr()
            << ", n(processes) = " << l.np()
            << ", n(systematics) = " << l.ns() << std::endl;
  for (unsigned int ir = 0; ir < l.nr(); ir++)
    for (unsigned int ip = 0; ip < l.np(); ip++) {
      l.set_nexp(ir, ip, nexp[ir][ip]);
      for (unsigned int is = 0; is < l.ns(); is++) {
        double dsyst = 0.;
        if (is < 2 && ip != 4) // Jet/MET
          dsyst = (0.1 - 0.05 * is) * (1. + ir * 0.05);
        else if (is < 7 && (is - 2 == ip)) // Theory
          dsyst = (ip == 4) ? 1. : 0.5;
        else if (is == 7 && (ir == 1 || ir == 2)) // b-tagging
          dsyst = (-1. + (ir - 1) * 2.) * 0.1;
        else if (is == 8) // Photon id efficiency
          dsyst = 0.05;
        l.set_dsyst(ir, ip, is, dsyst);
      }      
    }

  // Setup minimizer
  unsigned int nparams = l.np()+l.ns();
  boost::scoped_array<double> p0(new double[nparams]);
  boost::scoped_array<double> popt(new double[nparams]);
  for (unsigned int i = 0; i < nparams; i++)
    p0[i] = (i < l.np()) ? 1. : 0.;
  boost::function<double (double * const)> feval;
  feval = boost::bind(&likelihood::eval, &l, _1);
  optimize::minimizer_nm   mle_nm(feval, l.np()+l.ns());
  optimize::minimizer_bfgs mle_bfgs(feval, l.np()+l.ns());

  // Setup random generation of poisson distributed numbers
  boost::mt19937 gen;
  typedef boost::poisson_distribution<unsigned int> t_poisson;
  typedef boost::scoped_ptr< t_poisson > t_ptr_poisson;
  typedef boost::variate_generator< boost::mt19937&, t_poisson > t_gen;
  typedef boost::scoped_ptr< t_gen > t_ptr_gen;
  boost::scoped_array< t_ptr_poisson > pdist(new t_ptr_poisson[l.nr()]);
  boost::scoped_array< t_ptr_gen > rndgen(new t_ptr_gen[l.ns()]);
  for (unsigned int ir = 0; ir < l.nr(); ir++) {
    std::cout << "region" << ir << " nexp = " 
              << l.nexp(ir, p0.get()) << std::endl;
    pdist[ir].reset(new t_poisson(l.nexp(ir, p0.get())));
    rndgen[ir].reset(new t_gen(gen, *(pdist[ir])));
  }

  boost::timer time_monitor;
  double dt[3] = {0., 0., 0.};
  unsigned int ntoys = 10;
  for (unsigned int iexp = 0; iexp < ntoys; iexp++) {
    std::cout << "Starting experiment " << iexp << std::endl;

    // Generate pseudo-experiment random numbers
    for (unsigned int ir = 0; ir < l.nr(); ir++) {
      l.set_nobs(ir, (*(rndgen[ir]))());
      std::cout << "  region" << ir
                << " nobs = " << l.nobs(ir)
                << " nexp = " << l.nexp(ir, p0.get()) << std::endl;
    }

    std::cout << "  Perform a minimization " << std::endl;
    
    // Minimize with Nelder-Mean without mt
    mle_nm.set_is_mt(false);
    time_monitor.restart();
    mle_nm.minimize(p0.get());
    dt[0] = time_monitor.elapsed();
    std::cout << "  min = (";
    for (size_t i = 0; i < nparams; i++)
       popt[i] = mle_nm.get_opt_var(i).value();
    for (size_t i = 0; i < nparams; i++) {
      if (i < l.nr())
        std::cout << l.nexp(i, popt.get());
      else
        std::cout << mle_nm.get_opt_var(i).value();
      if (i < nparams - 1) std::cout << ", ";
    }
    std::cout << ") fmin = " << mle_nm.get_fmin() << std::endl;

    // Minimize with Nelder-Mean with mt
    mle_nm.set_is_mt(true);
    time_monitor.restart();
    mle_nm.minimize(p0.get());
    dt[1] = time_monitor.elapsed();
    std::cout << "  min = (";
    for (size_t i = 0; i < nparams; i++)
       popt[i] = mle_nm.get_opt_var(i).value();
    for (size_t i = 0; i < nparams; i++) {
      if (i < l.nr())
        std::cout << l.nexp(i, popt.get());
      else
        std::cout << mle_nm.get_opt_var(i).value();
      if (i < nparams - 1) std::cout << ", ";
    }
    std::cout << ") fmin = " << mle_nm.get_fmin() << std::endl;

    // Minimize with BFGS
    time_monitor.restart();
    mle_bfgs.minimize(p0.get());
    dt[2] = time_monitor.elapsed();
    std::cout << "  min = (";
    for (size_t i = 0; i < nparams; i++)
       popt[i] = mle_bfgs.get_opt_var(i).value();
    for (size_t i = 0; i < nparams; i++) {
      if (i < l.nr())
        std::cout << l.nexp(i, popt.get());
      else
        std::cout << mle_bfgs.get_opt_var(i).value();
      if (i < nparams - 1) std::cout << ", ";
    }
    std::cout << ") fmin = " << mle_bfgs.get_fmin() << std::endl;

    std::cout << "  timing " << dt[0] << " " << dt[1] << " " << dt[2] 
              << std::endl;
    std::cout << "  ... done." << std::endl;
  }

  return 0;
}
