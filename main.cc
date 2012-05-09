#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>

#include "hamiltonian.hh"
#include "matrix_term.hh"

int main(int argc, char **argv)
{
  hamiltonian ham
    (
     new exp_zb(AlAs),
     new exp_wang(AlAs, Be),
     new exp_overlap()
    );
  printf("%.1f\n", ham.get_eval(0));
  double level = ham.get_impurity_parameter(_V) - 4000;
  for(int i = 0; i >= -6500; i -= 500)
    {
      ham.set_impurity_parameter(_V, level - i);
      printf("%d: %.1f\n", i, ham.get_eval(0));
    }

  ham.set_crystal(new exp_wz(AlN));
  printf("\nWang+Chen potentials\n");

  ham.set_impurity(new exp_wang(AlN, Be));
  printf("Be_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 2500);
  printf("Be_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) + 4000);
  printf("Be_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(AlN, Mg));
  printf("Mg_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 2500);
  printf("Mg_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) + 4000);
  printf("Mg_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(AlN, Ca));
  printf("Ca_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 2500);
  printf("Ca_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) + 4000);
  printf("Ca_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(AlN, Zn));
  printf("Zn_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 2500);
  printf("Zn_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) + 4000);
  printf("Zn_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(AlN, Cd));
  printf("Cd_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 2500);
  printf("Cd_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) + 4000);
  printf("Cd_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(AlN, C));
  printf("C_N: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(AlN, Si));
  printf("Si_N: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(AlN, Ge));
  printf("Ge_N: %.1f\n", ham.get_eval(0));

  return 0;
}
