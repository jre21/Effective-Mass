#include <math.h>
#include <stdio.h>

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
  inv_radius = dielectric = 1.0;
}

matrix_term::~matrix_term()
{
  on_delete();
}

gsl_matrix_complex *matrix_term::matrix(double min, double max, size_t num)
{
  gsl_matrix_complex *work, *output =
    gsl_matrix_complex_alloc(matrix_term_block_size * num, matrix_term_block_size * num);
  gsl_matrix_complex_view view;

  // rescale min and max
  min *= inv_radius;
  max *= inv_radius;

  for(size_t i = 0; i < num; i++)
    for(size_t j = 0; j < num; j++)
      {
	// copy contents of each block into output
	work = matrix_block
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
  return inv_radius;
}

double matrix_term::set_inv_radius(double r)
{
  on_set_inv_radius(r);
  return inv_radius = r;
}

double matrix_term::get_dielectric_constant()
{
  return dielectric;
}

double matrix_term::set_dielectric_constant(double k)
{
  on_set_dielectric_constant(k);
  return dielectric = k;
}

// hooks provided for use by child classes
void matrix_term::on_set_inv_radius(double r) {}
void matrix_term::on_set_dielectric_constant(double k) {}
void matrix_term::on_delete() {}

// ########################### crystal_term ##########################
double crystal_term::get_parameter(crystal_parameters_t param)
{
  double ret = _get_parameter(param);
  if(1.0 / 0.0 == ret)
    {
      // mismatch between crystal term type and parameter
      printf("Error: crystal_term::get_parameter(): invalid parameter: %s\n",
	     crystal_parameters_to_string(param));
      exit(-1);
    }
  return ret;
}

double crystal_term::set_parameter(crystal_parameters_t param, double val)
{
  double ret = _set_parameter(param, val);
  if(1.0 / 0.0 == ret)
    {
      // mismatch between crystal term type and parameter
      printf("Error: crystal_term::set_parameter(): invalid parameter: %s\n",
	     crystal_parameters_to_string(param));
      exit(-1);
    }
  return ret;
}

double crystal_term::_get_parameter(crystal_parameters_t param)
{
  return 1.0 / 0.0;
}

double crystal_term::_set_parameter(crystal_parameters_t param, double val)
{
  return 1.0 / 0.0;
}

// ########################## impurity_term ##########################
double impurity_term::get_parameter(impurity_parameters_t param)
{
  double ret = _get_parameter(param);
  if(1.0 / 0.0 == ret)
    {
      // mismatch between impurity term type and parameter
      printf("Error: impurity_term::get_parameter(): invalid parameter: %s\n",
	     impurity_parameters_to_string(param));
      exit(-1);
    }
  return ret;
}

double impurity_term::set_parameter(impurity_parameters_t param, double val)
{
  double ret = _set_parameter(param, val);
  if(1.0 / 0.0 == ret)
    {
      // mismatch between impurity term type and parameter
      printf("Error: impurity_term::set_parameter(): invalid parameter: %s\n",
	     impurity_parameters_to_string(param));
      exit(-1);
    }
  return ret;
}

double impurity_term::_get_parameter(impurity_parameters_t param)
{
  return 1.0 / 0.0;
}

double impurity_term::_set_parameter(impurity_parameters_t param, double val)
{
  return 1.0 / 0.0;
}

