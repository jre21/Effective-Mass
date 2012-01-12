#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "hamiltonian.hh"
#include "matrix_term.hh"

// default basis granularity
#define BASIS_MIN 1e-2
#define BASIS_MAX 1e2
#define BASIS_NUM 30


hamiltonian::hamiltonian()
{
  crystal = NULL;
  impurity = NULL;
  overlap = NULL;
  basis_min = BASIS_MIN;
  basis_max = BASIS_MAX;
  basis_num = BASIS_NUM;
  inv_radius = dielectric = 1.0;
  clean_evals = 0;
  evals = NULL;
}

hamiltonian::hamiltonian(crystal_term *c, impurity_term *p, overlap_term *o)
{
  crystal = NULL;
  impurity = NULL;
  overlap = NULL;
  set_crystal(c);
  set_impurity(p);
  set_overlap(o);
  basis_min = BASIS_MIN;
  basis_max = BASIS_MAX;
  basis_num = BASIS_NUM;
  inv_radius = dielectric = 1.0;
  clean_evals = 0;
  evals = NULL;
}

hamiltonian::hamiltonian(crystal_term *c, impurity_term *p, overlap_term *o,
			 double min, double max, size_t num)
{
  crystal = NULL;
  impurity = NULL;
  overlap = NULL;
  set_crystal(c);
  set_impurity(p);
  set_overlap(o);
  basis_min = min;
  basis_max = max;
  basis_num = num;
  inv_radius = dielectric = 1.0;
  clean_evals = 0;
  evals = NULL;
}

hamiltonian::~hamiltonian()
{
  if(crystal) delete crystal;
  if(impurity) delete impurity;
  if(overlap) delete overlap;
  if(evals) delete evals;
}

void hamiltonian::set_crystal(crystal_term *c)
{
  // load crystal term
  crystal = c;
  // ensure old eigenvalues don't get used
  clean_evals = 0;
  // fetch dielectric constant
  dielectric = c->get_dielectric_constant();
  // set dielectric constant and inv_radius from crystal term and
  // propagate to other terms
  inv_radius = c->get_inv_radius();
  if(impurity) impurity->set_inv_radius(inv_radius);
  if(overlap) overlap->set_inv_radius(inv_radius);
  dielectric = c->get_dielectric_constant();
  if(impurity)
    impurity->set_dielectric_constant(inv_radius);
}

void hamiltonian::set_impurity(impurity_term *p)
{
  // load impurity term
  impurity = p;
  // ensure old eigenvalues don't get used
  clean_evals = 0;
  // pass in dielectric constant and inv_radius
  p->set_dielectric_constant(dielectric);
  p->set_inv_radius(inv_radius);
}

void hamiltonian::set_overlap(overlap_term *o)
{
  // load overlap term
  overlap = o;
  // ensure old eigenvalues don't get used
  clean_evals = 0;
  // pass in inv_radius
  o->set_inv_radius(inv_radius);
}

void hamiltonian::set_granularity(double min, double max, size_t num)
{
  // free eigenvalues if their size needs to be changed
  if(evals && basis_num != num)
    {
      gsl_vector_free(evals);
      evals = NULL;
    }
  // load granularity terms
  basis_min = min;
  basis_max = max;
  basis_num = num;
  // ensure old eigenvalues don't get used
  clean_evals = 0;
}

double hamiltonian::get_eval(int n)
{
  if(!clean_evals) gen_evals();
  return gsl_vector_get(evals, n);
}

void hamiltonian::gen_evals()
{
  // ensure Hamiltonian is fully built
  if(!crystal || !impurity || !overlap)
    {
      printf("Error: eigenvalues requested from incomplete Hamiltonian\n");
      exit(-1);
    }

  // build matrices
  gsl_matrix_complex *c =
    crystal->matrix(basis_min, basis_max, basis_num);
  gsl_matrix_complex *p =
    impurity->matrix(basis_min, basis_max, basis_num);
  gsl_matrix_complex *o =
    overlap->matrix(basis_min, basis_max, basis_num);
  gsl_matrix_complex_add(c,p);
  gsl_matrix_complex_free(p);
  p = NULL;

  // abort if passed nonhermitian matrices
  double norm = nonhermiticity(o);
  if(norm)
    {
      printf("Nonhermitian overlap matrix: %g\n", norm);
      exit(-1);
    }
  norm = nonhermiticity(c);
  if(norm)
    {
      printf("Nonhermitian hamiltonian matrix: %g\n", norm);
      exit(-1);
    }

  // generate eigenvalues
  if(!evals) evals = gsl_vector_alloc(c->size1);
  gsl_eigen_genherm_workspace *w = gsl_eigen_genherm_alloc(c->size1);
  gsl_eigen_genherm(c,o,evals,w);
  gsl_eigen_genherm_free(w);

  // sort eigenvalues (c now contains garbage, so we lose nothing by
  // reordering it)
  gsl_eigen_genhermv_sort(evals,c,GSL_EIGEN_SORT_VAL_ASC);
  // mark eigenvalues as usable
  clean_evals = 1;

  // clean up remaining local variables
  gsl_matrix_complex_free(c);
  gsl_matrix_complex_free(o);
}

// Calculate nonhermiticity index associated with each element m_ij
// (for i <= j) as abs(m_ij-conjugate(m_ji))^2 / abs(m_ij+m_ji)^2.
// Sum over indices which exceed THRESHOLD and return this value as
// the nonhermiticity index of the matrix.
double hamiltonian::nonhermiticity(gsl_matrix_complex *m)
{
#define THRESHOLD 1e-20
  double norm = 0.0, addend;
  for(size_t i = 0; i < m->size1; i++)
    for(size_t j = i; j < m->size1; j++)
      {
	addend = gsl_complex_abs2
	  (
	   gsl_complex_sub
	   (
	    gsl_matrix_complex_get(m,i,j),
	    gsl_complex_conjugate(gsl_matrix_complex_get(m,j,i))
	   )
	  )
	  / gsl_complex_abs2
	  (
	   gsl_complex_add
	   (
	    gsl_matrix_complex_get(m,i,j),
	    gsl_matrix_complex_get(m,j,i)
	   )
	  );
	if(addend > THRESHOLD) norm += addend;
      }
  return norm;
}
