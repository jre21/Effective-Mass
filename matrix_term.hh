#ifndef _MATRIX_TERM_HH
#define _MATRIX_TERM_HH

#include <gsl/gsl_matrix.h>

// Base class for matrix terms.  It needs to be able to generate a
// matrix based on the scale granularity parameters.
class matrix_term
{
public:
  ~matrix_term() {on_delete();}

  // Return a matrix based on arguments and internal state.
  // Deallocating this matrix is the caller's responsibility.
  // Min and max set minimum and maximum scales in multiples of the
  // inverse effective Bohr radius.  Num basis states will be used
  // scaling between these values in a geometric progression.
  virtual gsl_matrix_complex
  *matrix(double min, double max, unsigned int num) = 0;
protected:
  // Inheriting classes should override this with clean-up code.
  virtual void on_delete() {}
};

// the basic types of matrix terms
class crystal_term   : public matrix_term {};
class potential_term : public matrix_term {};
class overlap_term   : public matrix_term {};

#endif // _MATRIX_TERM_HH
