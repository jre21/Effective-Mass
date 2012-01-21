#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>

#include "hamiltonian.hh"
#include "matrix_term.hh"

int main(int argc, char **argv)
{
  hamiltonian ham
    (
     new exp_zb(GaN),
     new exp_coulomb(),
     new exp_overlap()
    );
  printf("%g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(H, H));
  printf("%g\n", ham.get_eval(0));
  ham.set_impurity(new exp_LCZ(Ga, Zn));
  printf("%g\n", ham.get_eval(0));
  ham.set_impurity(new gauss_LCZ(Ga, Zn));
  printf("%g\n", ham.get_eval(0));
  ham.set_impurity(new gauss_HGH(Ga, Zn));
  printf("%g\n", ham.get_eval(0));
  return 0;
}
