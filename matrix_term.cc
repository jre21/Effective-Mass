#include <math.h>

#include <gsl/gsl_matrix.h>

#include "matrix_term.hh"

// size of matrix blocks.  This is the number of angular momentum
// states (four spin-3/2 states plus two spin-1/2 states) multiplied
// by the number of spherical harmonics used (currently one l=0 plus
// five l=2).  This should only be modified if the mathematica-
// generated matrix blocks are changed to include more spherical
// harmonics.
int matrix_term_block_size = 36;

// ########################### matrix_term ###########################
matrix_term::matrix_term()
{
  this->inv_radius = this->dielectric = 1.0;
}

matrix_term::~matrix_term()
{
  this->on_delete();
}

gsl_matrix_complex *matrix_term::matrix(double min, double max, size_t num)
{
  gsl_matrix_complex *work, *output =
    gsl_matrix_complex_alloc(matrix_term_block_size * num, matrix_term_block_size * num);
  gsl_matrix_complex_view view;

  // rescale min and max
  min *= this->inv_radius;
  max *= this->inv_radius;

  for(size_t i = 0; i < num; i++)
    for(size_t j = 0; j < num; j++)
      {
	// copy contents of each block into output
	work = this->matrix_block
	  (
	   min * pow(max/min,(double)i/(num-1)),
	   min * pow(max/min,(double)j/(num-1))
	  );
	view = gsl_matrix_complex_submatrix
	  (output, matrix_term_block_size * i, matrix_term_block_size * j, matrix_term_block_size, matrix_term_block_size);
	gsl_matrix_complex_memcpy(&view.matrix, work);
	gsl_matrix_complex_free(work);
      }
  return output;
}

double matrix_term::get_inv_radius()
{
  return this->inv_radius;
}

double matrix_term::set_inv_radius(double r)
{
  this->on_set_inv_radius(r);
  return this->inv_radius = r;
}

double matrix_term::get_dielectric_constant()
{
  return this->dielectric;
}

double matrix_term::set_dielectric_constant(double k)
{
  this->on_set_dielectric_constant(k);
  return this->dielectric = k;
}

// hooks provided for use by child classes
void matrix_term::on_set_inv_radius(double r) {}
void matrix_term::on_set_dielectric_constant(double k) {}
void matrix_term::on_delete() {}
