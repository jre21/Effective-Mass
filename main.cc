#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>

#include "hamiltonian.hh"
#include "matrix_term.hh"

int main(int argc, char **argv)
{
  hamiltonian ham
    (
     new exp_wz(GaN),
     new exp_coulomb(),
     new exp_overlap()
    );
  printf("GaN: %g\n", ham.get_eval(0));
  ham.set_crystal(new exp_wz(AlN));
  printf("AlN: %g\n", ham.get_eval(0));
  ham.set_crystal(new exp_wz(InN));
  printf("InN: %g\n", ham.get_eval(0));
  ham.set_crystal(new exp_gwz(ZnGeN2));
  printf("AlN: %g\n", ham.get_eval(0));
  ham.set_crystal(new exp_gwz(ZnSnN2));
  printf("AlN: %g\n", ham.get_eval(0));
  return 0;
}
