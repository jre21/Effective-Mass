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
  this->crystal = NULL;
  this->potential = NULL;
  this->overlap = NULL;
  this->basis_min = BASIS_MIN;
  this->basis_max = BASIS_MAX;
  this->basis_num = BASIS_NUM;
  this->inv_radius = this->dielectric = 1.0;
  this->clean_evals = 0;
  this->evals = NULL;
}

hamiltonian::hamiltonian(crystal_term *c, potential_term *p, overlap_term *o)
{
  this->set_crystal(c);
  this->set_potential(p);
  this->set_overlap(o);
  this->basis_min = BASIS_MIN;
  this->basis_max = BASIS_MAX;
  this->basis_num = BASIS_NUM;
  this->inv_radius = this->dielectric = 1.0;
  this->clean_evals = 0;
  this->evals = NULL;
}

hamiltonian::hamiltonian(crystal_term *c, potential_term *p, overlap_term *o,
			 double min, double max, size_t num)
{
  this->set_crystal(c);
  this->set_potential(p);
  this->set_overlap(o);
  this->basis_min = min;
  this->basis_max = max;
  this->basis_num = num;
  this->inv_radius = this->dielectric = 1.0;
  this->clean_evals = 0;
  this->evals = NULL;
}

hamiltonian::~hamiltonian()
{
  if(this->crystal) delete this->crystal;
  if(this->potential) delete this->potential;
  if(this->overlap) delete this->overlap;
  if(this->evals) delete this->evals;
}

void hamiltonian::set_crystal(crystal_term *c)
{
  // load crystal term
  this->crystal = c;
  // ensure old eigenvalues don't get used
  this->clean_evals = 0;
  // pass in dielectric constant
  c->set_dielectric_constant(this->dielectric);
  // set inv_radius from crystal term and propagate to other terms
  this->inv_radius = c->get_inv_radius();
  if(this->potential) this->potential->set_inv_radius(this->inv_radius);
  if(this->overlap) this->overlap->set_inv_radius(this->inv_radius);
}

void hamiltonian::set_potential(potential_term *p)
{
  // load potential term
  this->potential = p;
  // ensure old eigenvalues don't get used
  this->clean_evals = 0;
  // pass in dielectric constant and inv_radius
  p->set_dielectric_constant(this->dielectric);
  p->set_inv_radius(this->inv_radius);
}

void hamiltonian::set_overlap(overlap_term *o)
{
  // load overlap term
  this->overlap = o;
  // ensure old eigenvalues don't get used
  this->clean_evals = 0;
  // pass in inv_radius
  o->set_inv_radius(this->inv_radius);
}

void hamiltonian::set_granularity(double min, double max, size_t num)
{
  // free eigenvalues if their size needs to be changed
  if(this->evals && this->basis_num != num)
    {
      gsl_vector_free(this->evals);
      this->evals = NULL;
    }
  // load granularity terms
  this->basis_min = min;
  this->basis_max = max;
  this->basis_num = num;
  // ensure old eigenvalues don't get used
  this->clean_evals = 0;
}

double hamiltonian::set_dielectric_constant(double k)
{
  // load dielectric constant and propagate to the terms that depend
  // on it
  this->dielectric = k;
  if(this->crystal) this->crystal->set_dielectric_constant(k);
  if(this->potential) this->potential->set_dielectric_constant(k);
  return k;
}

double hamiltonian::get_eval(int n)
{
  if(!this->clean_evals) this->gen_evals();
  return gsl_vector_get(this->evals, n);
}

void hamiltonian::gen_evals()
{
  // ensure Hamiltonian is fully built
  if(!this->crystal || !this->potential || !this->overlap)
    {
      printf("Error: eigenvalues requested from incomplete Hamiltonian\n");
      exit(-1);
    }

  // build matrices
  gsl_matrix_complex *c =
    this->crystal->matrix(this->basis_min, this->basis_max, this->basis_num);
  gsl_matrix_complex *p =
    this->potential->matrix(this->basis_min, this->basis_max, this->basis_num);
  gsl_matrix_complex *o =
    this->overlap->matrix(this->basis_min, this->basis_max, this->basis_num);
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
  if(!this->evals) this->evals = gsl_vector_alloc(c->size1);
  gsl_eigen_genherm_workspace *w = gsl_eigen_genherm_alloc(c->size1);
  gsl_eigen_genherm(c,o,this->evals,w);
  gsl_eigen_genherm_free(w);

  // sort eigenvalues (c now contains garbage, so we lose nothing by
  // reordering it)
  gsl_eigen_genhermv_sort(this->evals,c,GSL_EIGEN_SORT_VAL_ASC);
  // mark eigenvalues as usable
  this->clean_evals = 1;

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
