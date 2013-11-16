/*!\file   variable.hpp
 * \brief  Define a class for a function variable.
 * \author Renaud Bruneliere
 * \date   13.11.2013
 */

#ifndef OPTIMIZE_PARAM_HPP
#define OPTIMIZE_PARAM_HPP

namespace optimize {

class function;

//! Enum defining the different possible variable bounds
enum VarBd { UNBOUNDED = 0, LOWER_BOUND, UPPER_BOUND, BOUNDED };

//! Class defining a real-valued variable
/*!
  Many optimization algorithms act only on unconstrained objective functions.
  For very simple constrains on variables (lower/upper bound(s) independent of
  the other variables), a 'trick' consists into the transformation of a bounded
  variable into an unbounded variable via a non-linear transformation.
  To handle this type of simple constrains, the object variable has attributes
  defining its type of boundaries and possibly their values.
  Also, as all functions used in this package require double* as input values,
  the variable object does not always own its value but can point it toward an 
  array. 
*/
class variable
{
 public:
  //! Constructor with possibly a variable value.
  variable(double val = 0.);

  //! Destructor
  ~variable() {}

  //! Set the pointer to its value
  void set_value_ptr(double * ptr) { ptr_ = ptr; }
  //! Get the variable value
  double value() const { return *ptr_; }
  //! Get the variable value in an unconstrained representation
  /*!
   * Transform the variable value from a bounded to an unbounded 
   * representation via non-linear formula.
   * When the variableeter is one-sided, the sqrt function is used:
   *   new = sqrt(pow(value - lower + 1, 2) - 1)
   *   new = sqrt(pow(upper - value + 1, 2) - 1)
   * When the variableeter is two-sided, the arcsin function is used:
   *   new = asin(2*(value - lower)/(upper - lower) - 1)
   */
  double value_ubd() const;
  //! Set the variable value
  void set_value(double val) { *ptr_ = val; }
  //! Set the variable value using an unconstrained representation as input
  void set_value_ubd(double val);
  //! Get the bound type
  VarBd bound_type() const { return bd_type_; }
  //! Set the bound type
  void set_bound_type(VarBd type);
  //! Get the lower bound
  double lower_bound() const { return lower_bd_; }
  //! Set the lower bound (be cautious, it affects the bound type)
  void set_lower_bound(double lower_bd);
  //! Get the upper bound
  double upper_bound() const { return upper_bd_; }
  //! Set the upper bound (be cautious, it affects the bound type)
  void set_upper_bound(double upper_bd);
  //! Copy bounds from another variable (safe)
  void copy_bounds(const variable& var);
  //! Set the pointer to its function
  void set_function_ptr(function * const fptr) { fptr_ = fptr; } 

  //! Default bounded <-> unbounded function, no transformation is done 
  static double bd_to_ubd_none(double val, double lower, double upper);
  //! Transform a one-sided lower bounded variable to an unbounded variable
  static double bd_to_ubd_lower(double val, double lower, double upper);
  //! Transform an unbounded variable to a one-sided lower bounded variable
  static double ubd_to_bd_lower(double val, double lower, double upper);
  //! Transform a one-sided upper bounded variable to an unbounded variable
  static double bd_to_ubd_upper(double val, double lower, double upper);
  //! Transform an unbounded variable to a one-sided upper bounded variable
  static double ubd_to_bd_upper(double val, double lower, double upper);
  //! Transform a two-sided bounded variable to an unbounded variable
  static double bd_to_ubd_both(double val, double lower, double upper);
  //! Transform an unbounded variable to a two-sided bounded variable
  static double ubd_to_bd_both(double val, double lower, double upper);

 private:
  double val_;      ///< Parameter value (optional)
  double * ptr_;    ///< Pointer to its value
  VarBd bd_type_;   ///< Type of bounds
  double lower_bd_; ///< Parameter lower bound
  double upper_bd_; ///< Parameter upper bound
  //! Bounded to unbounded function pointer
  double (*bd_to_ubd_)(double, double, double);
  //! Unbounded to bounded function pointer
  double (*ubd_to_bd_)(double, double, double);
  function * fptr_;  ///< Pointer to its function (optional)
};

} // namespace optimize

#endif // OPTIMIZE_PARAM_HPP
