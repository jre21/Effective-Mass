#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "hamiltonian.hh"

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
}

double hamiltonian::get_eval(int n)
{
  if(!clean_evals) gen_evals();
  return gsl_vector_get(this->evals, n);
}

// stub for now
void gen_eivals()
{
}
