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
  for(double i = ham.get_crystal_parameter(_d1); i < 386;
      i = ham.get_crystal_parameter(_d1))
    {
      printf("%.1f: %f\n", i, ham.get_eval(0));
      i = ham.get_crystal_parameter(_d1);
      ham.set_crystal_parameter(_d1, i+5);
    }

  ham.set_crystal(new exp_wz(AlN));
  ham.set_impurity(new exp_LCZ(Al, Mg));
  printf("Mg_Al in AlN: %f\n", ham.get_eval(0));

  ham.set_crystal(new exp_wz(GaN));
  ham.set_impurity(new exp_LCZ(Ga, Mg));
  printf("Mg_Ga in GaN: %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga, Zn));
  printf("Zn_Ga in GaN (3d as core): %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga_val, Zn));
  printf("Zn_Ga in GaN (3d as valance): %f\n", ham.get_eval(0));

  ham.set_crystal(new exp_wz(InN));
  ham.set_impurity(new exp_LCZ(Ga, Zn));
  printf("Zn_Ga in InN (3d as core): %f\n", ham.get_eval(0));

  ham.set_impurity(new exp_LCZ(Ga_val, Zn));
  printf("Zn_Ga in InN (3d as valance): %f\n", ham.get_eval(0));

  return 0;
}
