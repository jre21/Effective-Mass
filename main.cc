#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>

#include "hamiltonian.hh"
#include "matrix_term.hh"

int main(int argc, char **argv)
{
  hamiltonian ham
    (
     new exp_wz(AlN),
     new exp_coulomb(),
     new exp_overlap()
    );
  printf("AlN: Energy vs crystal field splitting\n");
  ham.set_crystal_parameter(_d1,105);
  for(double i = ham.get_crystal_parameter(_d1); i < 385;
      i = ham.get_crystal_parameter(_d1))
    {
      printf("%.1f: %f\n", i, ham.get_eval(0));
      i = ham.get_crystal_parameter(_d1);
      ham.set_crystal_parameter(_d1, i-1);
    }

  return 0;
}
