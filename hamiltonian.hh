#ifndef _HAMILTONIAN_HH
#define _HAMILTONIAN_HH 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// classes defining the matrix terms used
#ifndef _MATRIX_TERM_HH
class crystal_term;
class potential_term;
class overlap_term;
#endif

// class to construct and modify the hamiltonian matrix and determine
// eigenvalues (might implement eigenvectors later)
class hamiltonian
{
public:
  hamiltonian();
  hamiltonian(crystal_term *c, potential_term *p, overlap_term *o);
  hamiltonian(crystal_term *c, potential_term *p, overlap_term *o,
	      double min, double max, size_t num);
  ~hamiltonian();

  // Set hamiltonian terms.  Note that these pointers should not be
  // reused between hamiltonian instances, because they will be
  // deleted when no longer needed.
  void set_crystal(crystal_term *c);
  void set_potential(potential_term *p);
  void set_overlap(overlap_term *o);

  // Override default granularity of basis states.  Min and max set
  // minimum and maximum scales in multiples of the inverse effective
  // Bohr radius.  Num basis states will be used scaling between these
  // values in a geometric progression.
  void set_granularity(double min, double max, size_t num);

  // Get the n'th eigenvalue, recalculating all eigenvalues if
  // necessary.
  double get_eval(int n);
protected:
  // generate eigenvalues
  void gen_evals();
  // Measure degree by which matrix fails to be hermitian, and return
  // 0 if the failure is within an allowed tolerance attributable to
  // roundoff error.
  double nonhermiticity(gsl_matrix_complex *m);

  // pointers to the current term objects
  crystal_term *crystal;
  potential_term *potential;
  overlap_term *overlap;

  // constants defining basis granularity
  double basis_min;
  double basis_max;
  size_t basis_num;

  // true if eigenvalues are generated from current parameters
  int clean_evals;

  // calculated eigenvalues
  gsl_vector *evals;
};

#endif // _HAMILTONIAN_HH
