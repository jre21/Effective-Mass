#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>

#include "hamiltonian.hh"
#include "matrix_term.hh"

int main(int argc, char **argv)
{
  hamiltonian ham
    (
     new exp_zb(3.76, 0.90, 1.42 ,300, 10.06),
     new exp_wang(GaN, Be),
     new exp_overlap()
    );
  double level = ham.get_impurity_parameter(_V);
  for(int i = -2500; i >= -6500; i -= 500)
    {
      ham.set_impurity_parameter(_V, level - i);
      printf("%d: %.1f\n", i, ham.get_eval(0));
    }

  ham.set_crystal(new exp_wz(AlN));
  printf("\nWang+Chen potentials\n");

  ham.set_impurity(new exp_wang(GaN, Be));
  printf("Be_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 1200);
  printf("Be_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(GaN, Mg));
  printf("Mg_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 1200);
  printf("Mg_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(GaN, Ca));
  printf("Ca_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 1200);
  printf("Ca_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(GaN, Zn));
  printf("Zn_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 1200);
  printf("Zn_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(GaN, Cd));
  printf("Cd_III: %.1f\n", ham.get_eval(0));
  ham.set_impurity_parameter(_V, ham.get_impurity_parameter(_V) - 1200);
  printf("Cd_III: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(GaN, C));
  printf("C_N: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(GaN, Si));
  printf("Si_N: %.1f\n", ham.get_eval(0));

  ham.set_impurity(new exp_wang(GaN, Ge));
  printf("Ge_N: %.1f\n", ham.get_eval(0));

  return 0;
}
