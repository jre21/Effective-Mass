#include <stdio.h>
#include <stdlib.h>

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
  this->crystal = 0;
  this->potential = 0;
  this->overlap = 0;
  this->basis_min = BASIS_MIN;
  this->basis_max = BASIS_MAX;
  this->basis_num = BASIS_NUM;
  this->evals = 0;
}

hamiltonian::hamiltonian(crystal_term *c, potential_term *p, overlap_term *o)
{
  this->crystal = c;
  this->potential = p;
  this->overlap = o;
  this->basis_min = BASIS_MIN;
  this->basis_max = BASIS_MAX;
  this->basis_num = BASIS_NUM;
  this->evals = 0;
}

hamiltonian::hamiltonian(crystal_term *c, potential_term *p, overlap_term *o,
			 double min, double max, unsigned int num)
{
  this->crystal = c;
  this->potential = p;
  this->overlap = o;
  this->basis_min = min;
  this->basis_max = max;
  this->basis_num = num;
  this->evals = 0;
}

void hamiltonian::set_crystal(crystal_term *c)
{
  this->crystal = c;
  this->clean_evals = 0;
}

void hamiltonian::set_potential(potential_term *p)
{
  this->potential = p;
  this->clean_evals = 0;
}

void hamiltonian::set_overlap(overlap_term *o)
{
  this->overlap = o;
  this->clean_evals = 0;
}

void hamiltonian::set_granularity(double min, double max, unsigned int num)
{
  this->basis_min = min;
  this->basis_max = max;
  this->basis_num = num;
  this->clean_evals = 0;
  if(this->evals)
    {
      gsl_vector_free(this->evals);
      this->evals = 0;
    }
}

double hamiltonian::get_eval(int n)
{
  if(!clean_evals) gen_evals();
  return gsl_vector_get(this->evals, n);
}

void hamiltonian::gen_evals()
{
  // build matrices
  gsl_matrix_complex *c =
    this->crystal->matrix(basis_min, basis_max, basis_num);
  gsl_matrix_complex *p =
    this->potential->matrix(basis_min, basis_max, basis_num);
  gsl_matrix_complex *o =
    this->overlap->matrix(basis_min, basis_max, basis_num);
  gsl_matrix_complex_add(c,p);
  gsl_matrix_complex_free(p);
  p = 0;

  // check for non-hermiticity
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
  if(!this->evals) this->evals = gsl_vector_alloc(c->size1);
  gsl_eigen_genherm_workspace *w = gsl_eigen_genherm_alloc(c->size1);
  gsl_eigen_genherm(c,o,evals,w);
  gsl_eigen_genherm_free(w);

  // sort eigenvalues (c now contains garbage, so we lose nothing by
  // reordering it)
  gsl_eigen_genhermv_sort(this->evals,c,GSL_EIGEN_SORT_VAL_ASC);
  // mark eigenvalues as usable
  this->clean_evals=1;

  // clean up remaining local variables
  gsl_matrix_complex_free(c);
  gsl_matrix_complex_free(o);
}

// stub for now
double hamiltonian::nonhermiticity(gsl_matrix_complex *m)
{
  return 0;
}
