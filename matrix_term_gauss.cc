#include <math.h>
#include <stdio.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_erf.h>

#include "matrix_term.hh"
#include "defs.hh"

// defined in matrix_term.cc
extern int matrix_term_block_size;

// ############################# gauss_zb ############################
gauss_zb::gauss_zb
(double _g1, double _g2, double _g3, double _d0, double _dielectric)
{
  g1 = _g1;
  g2 = _g2;
  g3 = _g3;
  d0 = _d0 * RYD_PER_MEV;
  dielectric = _dielectric;
  inv_radius = 1.0 / _g1 / _dielectric;
}

gauss_zb::gauss_zb(crystals_t c)
{
  switch(c)
    {
    case GaN:
      g1 = 2.46;
      g2 = 0.65;
      g3 = 0.98;
      d0 = 19 * RYD_PER_MEV;
      dielectric = 10.13;
      break;
    case AlN:
      g1 = 1.40;
      g2 = 0.35;
      g3 = 0.59;
      d0 = 19 * RYD_PER_MEV;
      dielectric = 9.67;
      break;
    default:
      printf("Error: gauss_zb instantiated with unknown crystal");
      exit(-1);
    }
  inv_radius = 1.0 / g1 / dielectric;

}

gsl_matrix_complex *gauss_zb::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_zb_def.hh"
  return output;
}


// ############################# gauss_wz ############################
gauss_wz::gauss_wz(double _A1, double _A2, double _A3, double _A4,
		   double _A5, double _A6, double _d1, double _d2,
		   double _d3, double _dielectric)
{
  A1 = _A1;
  A2 = _A2;
  A3 = _A3;
  A4 = _A4;
  A5 = _A5;
  A6 = _A6;
  d1 = _d1 * RYD_PER_MEV;
  d2 = _d2 * RYD_PER_MEV;
  d3 = _d3 * RYD_PER_MEV;
  dielectric = _dielectric;
  inv_radius = 1.0 / (_A2 + _A4) / _dielectric;
}

gauss_wz::gauss_wz(crystals_t c)
{
  switch(c)
    {
    case GaN:
      A1 = 5.98;
      A2 = 0.58;
      A3 = -5.44;
      A4 = 2.46;
      A5 = 2.53;
      A6 = 1.55;
      d1 = -12 * RYD_PER_MEV;
      d2 = -3.9 * RYD_PER_MEV;
      d3 = -5.4 * RYD_PER_MEV;
      dielectric = 10.13;
      break;
    case AlN:
      A1 = 4.05;
      A2 = 0.28;
      A3 = -3.71;
      A4 = 1.71;
      A5 = 1.90;
      A6 = 1.05;
      d1 = 245 * RYD_PER_MEV;
      d2 = -6.2 * RYD_PER_MEV;
      d3 = -7.5 * RYD_PER_MEV;
      dielectric = 9.67;
      break;
    case InN:
      A1 = 15.72;
      A2 = 0.7;
      A3 = -15.23;
      A4 = 5.79;
      A5 = 5.9;
      A6 = 4.58;
      d1 = 14.8 * RYD_PER_MEV;
      d2 = 3.4 * RYD_PER_MEV;
      d3 = -2.2 * RYD_PER_MEV;
      dielectric = 14;
      break;
    default:
      printf("Error: gauss_wz instantiated with unknown crystal");
      exit(-1);
    }
  inv_radius = 1.0 / (A2 + A4) / dielectric;
}


gsl_matrix_complex *gauss_wz::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_wz_def.hh"
  return output;
}


// ############################ gauss_gwz ############################
gauss_gwz::gauss_gwz(double _A1, double _A2, double _A3, double _B1,
		     double _B2, double _B3, double _C1, double _C2,
		     double _C3, double _D1, double _D2, double _D3,
		     double _d1c, double _d2c, double _d1so,
		     double _d2so, double _d3so, double _dielectric)
{
  A1 = _A1;
  A2 = _A2;
  A3 = _A3;
  B1 = _B1;
  B2 = _B2;
  B3 = _B3;
  C1 = _C1;
  C2 = _C2;
  C3 = _C3;
  D1 = _D1;
  D2 = _D2;
  D3 = _D3;
  d1c = _d1c * RYD_PER_MEV;
  d2c = _d2c * RYD_PER_MEV;
  d1so = _d1so * RYD_PER_MEV;
  d2so = _d2so * RYD_PER_MEV;
  d3so = _d3so * RYD_PER_MEV;
  dielectric = _dielectric;
  inv_radius = 1.0 / (_B1 + _B2) / _dielectric;
}

gauss_gwz::gauss_gwz(crystals_t c)
{
  switch(c)
    {
    case ZnGeN2:
      A1 = 6.82;
      A2 = -6.44;
      A3 = -0.01;
      B1 = 0.51;
      B2 = 2.19;
      B3 = -0.09;
      C1 = 0.02;
      C2 = 0.05;
      C3 = -2.30;
      D1 = -4.6;
      D2 = -2.76;
      D3 = -2.76;
      d1c = -115 * RYD_PER_MEV;
      d2c = 14 * RYD_PER_MEV;
      d1so = 0 * RYD_PER_MEV;
      d2so = 0 * RYD_PER_MEV;
      d3so = 0 * RYD_PER_MEV;
      dielectric = 9.7;
      break;
    case ZnSnN2:
      A1 = 8.57;
      A2 = -8.10;
      A3 = -0.02;
      B1 = 0.53;
      B2 = 3.12;
      B3 = -0.11;
      C1 = 0.03;
      C2 = 0.05;
      C3 = -3.18;
      D1 = -6.36;
      D2 = -4.62;
      D3 = -4.62;
      d1c = -82 * RYD_PER_MEV;
      d2c = 94 * RYD_PER_MEV;
      d1so = 0 * RYD_PER_MEV;
      d2so = 0 * RYD_PER_MEV;
      d3so = 0 * RYD_PER_MEV;
      dielectric = 12.71;
      break;
    default:
      printf("Error: gauss_gwz instantiated with unknown crystal");
      exit(-1);
    }
  inv_radius = 1.0 / (B1 + B2) / dielectric;
}

gsl_matrix_complex *gauss_gwz::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_gwz_def.hh"
  return output;
}


// ########################### gauss_coulomb #########################
gauss_coulomb::gauss_coulomb()
{
  dielectric = inv_radius = 1.0;
}

gsl_matrix_complex *gauss_coulomb::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_coulomb_def.hh"
  return output;
}


// ############################ gauss_wang ###########################
gauss_wang::gauss_wang(double _V, double _ra, double _rb, double _r1)
{
  V = _V;
  ra = _ra;
  rb = _rb;
  r1 = _r1;
  dielectric = inv_radius = 1.0;
}

gsl_matrix_complex *gauss_wang::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_well_def.hh"
  gsl_matrix_complex *out2 = output;
  output = gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_dielectric_def.hh"
  gsl_matrix_complex_add(output, out2);
  gsl_matrix_complex_free(out2);
  return output;
}


// ############################ gauss_LCZ ############################
gauss_LCZ::gauss_LCZ(elements_t _host, elements_t _impurity)
{
  host = new gauss_LCZ_atom(_host);
  impurity = new gauss_LCZ_atom(_impurity);
  coulomb = new gauss_coulomb();
}

gsl_matrix_complex *gauss_LCZ::matrix(double min, double max, size_t num)
{
  gsl_matrix_complex *h = host->matrix(min, max, num),
    *i = impurity->matrix(min, max, num),
    *c = coulomb->matrix(min, max, num);
  gsl_matrix_complex_sub(h,i);
  gsl_matrix_complex_add(h,c);
  gsl_matrix_complex_free(i);
  gsl_matrix_complex_free(c);
  return h;
}

// no-op because blocks are calculated in the lower-level terms
gsl_matrix_complex *gauss_LCZ::matrix_block(double a1, double a2)
{
  return NULL;
}

void gauss_LCZ::on_set_inv_radius(double r)
{
  host->set_inv_radius(r);
  impurity->set_inv_radius(r);
  coulomb->set_inv_radius(r);
}

void gauss_LCZ::on_set_dielectric_constant(double k)
{
  host->set_dielectric_constant(k);
  impurity->set_dielectric_constant(k);
  coulomb->set_dielectric_constant(k);
}

void gauss_LCZ::on_delete()
{
  delete host;
  delete impurity;
  delete coulomb;
}


// ######################### gauss_LCZ_atom ##########################
gauss_LCZ_atom::gauss_LCZ_atom(elements_t atom)
{
  set_d_core(atom);
}

void gauss_LCZ_atom::set_d_core(elements_t atom)
{
  switch(atom)
    {
    case Li:
      Zc = 2;
      C1 = 2.9990;
      C2 = 1.0472;
      C3 = 1.9686;
      break;
    case Be:
      Zc = 2;
      C1 = 3.0169;
      C2 = 1.4229;
      C3 = 2.9456;
      break;
    case B:
      Zc = 2;
      C1 = 3.0315;
      C2 = 1.8014;
      C3 = 3.9668;
      break;
    case C:
      Zc = 2;
      C1 = 3.0731;
      C2 = 2.1922;
      C3 = 5.0062;
      break;
    case N:
      Zc = 2;
      C1 = 3.0487;
      C2 = 2.5638;
      C3 = 6.0734;
      break;
    case O:
      Zc = 2;
      C1 = 3.0500;
      C2 = 2.9600;
      C3 = 6.8700;
      break;
    case F:
      Zc = 2;
      C1 = 3.0744;
      C2 = 3.3609;
      C3 = 7.7131;
      break;
    case Ne:
      Zc = 2;
      C1 = 3.0765;
      C2 =  3.7456;
      C3 = 8.7366;
      break;
    case Na:
      Zc = 10;
      C1 = 8.6691;
      C2 = 1.354;
      C3 = 2.3326;
      break;
    case Mg:
      Zc = 10;
      C1 = 8.0022;
      C2 = 1.4094;
      C3 = 2.6122;
      break;
    case Al:
      Zc = 10;
      C1 = 8.0223;
      C2 = 1.5274;
      C3 = 3.0200;
      break;
    case Si:
      Zc = 10;
      C1 = 8.0929;
      C2 = 1.6490;
      C3 = 3.4389;
      break;
    case P:
      Zc = 10;
      C1 = 8.3029;
      C2 = 1.8115;
      C3 = 3.8615;
      break;
    case S:
      Zc = 10;
      C1 = 8.4635;
      C2 = 1.9634;
      C3 = 4.3209;
      break;
    case Cl:
      Zc = 10;
      C1 = 8.6365;
      C2 = 2.1205;
      C3 = 4.8235;
      break;
    case Ar:
      Zc = 10;
      C1 = 8.8225;
      C2 = 2.2825;
      C3 = 5.4118;
      break;
    case K:
      Zc = 18;
      C1 = 10.615;
      C2 = 0.8133;
      C3 = 1.5698;
      break;
    case Ca:
      Zc = 18;
      C1 = 11.154;
      C2 = 0.8898;
      C3 = 1.8822;
      break;
    case Sc:
      Zc = 18;
      C1 = 11.706;
      C2 = 0.9809;
      C3 = 2.1414;
      break;
    case Ti:
      Zc = 18;
      C1 = 12.210;
      C2 = 1.0726;
      C3 = 2.3685;
      break;
    case V:
      Zc = 18;
      C1 = 12.655;
      C2 = 1.1619;
      C3 = 2.5793;
      break;
    case Cr:
      Zc = 18;
      C1 = 12.781;
      C2 = 1.2449;
      C3 = 2.7634;
      break;
    case Mn:
      Zc = 18;
      C1 = 13.257;
      C2 = 1.3321;
      C3 = 2.9744;
      break;
    case Fe:
      Zc = 18;
      C1 = 13.661;
      C2 = 1.4247;
      C3 = 3.1611;
      break;
    case Co:
      Zc = 18;
      C1 = 13.761;
      C2 = 1.4997;
      C3 = 3.3435;
      break;
    case Ni:
      Zc = 18;
      C1 = 14.134;
      C2 = 1.5925;
      C3 = 3.5240;
      break;
    case Cu:
      Zc = 18;
      C1 = 14.230;
      C2 = 1.6794;
      C3 = 3.6701;
      break;
    case Zn:
      Zc = 18;
      C1 = 14.557;
      C2 = 1.7601;
      C3 = 3.8772;
      break;
    case Ga:
      Zc = 28;
      C1 = 14.934;
      C2 = 1.8452;
      C3 = 3.3476;
      break;
    case Ge:
      Zc = 28;
      C1 = 14.888;
      C2 = 1.8923;
      C3 = 3.5107;
      break;
    case As:
      Zc = 28;
      C1 = 14.954;
      C2 = 1.9544;
      C3 = 3.6166;
      break;
    case Se:
      Zc = 28;
      C1 = 15.171;
      C2 = 2.0278;
      C3 = 3.8322;
      break;
    case Br:
      Zc = 28;
      C1 = 15.224;
      C2 = 2.0910;
      C3 = 4.0707;
      break;
    default:
      printf("Error: LCZ pseudopotential not available for %s\n",
	     elements_to_string(atom));
    }
}

gsl_matrix_complex *gauss_LCZ_atom::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_lcz_def.hh"
  return output;
}


// ############################ gauss_HGH ############################
gauss_HGH::gauss_HGH(elements_t _host, elements_t _impurity)
{
  host = new gauss_HGH_atom(_host);
  impurity = new gauss_HGH_atom(_impurity);
  coulomb = new gauss_coulomb();
}

gsl_matrix_complex *gauss_HGH::matrix(double min, double max, size_t num)
{
  gsl_matrix_complex *h = host->matrix(min, max, num),
    *i = impurity->matrix(min, max, num),
    *c = coulomb->matrix(min, max, num);
  gsl_matrix_complex_sub(h,i);
  gsl_matrix_complex_add(h,c);
  gsl_matrix_complex_free(i);
  gsl_matrix_complex_free(c);
  return h;
}

// no-op because blocks are calculated in the lower-level terms
gsl_matrix_complex *gauss_HGH::matrix_block(double a1, double a2)
{
  return NULL;
}

void gauss_HGH::on_set_inv_radius(double r)
{
  host->set_inv_radius(r);
  impurity->set_inv_radius(r);
  coulomb->set_inv_radius(r);
}

void gauss_HGH::on_set_dielectric_constant(double k)
{
  host->set_dielectric_constant(k);
  impurity->set_dielectric_constant(k);
  coulomb->set_dielectric_constant(k);
}

void gauss_HGH::on_delete()
{
  delete host;
  delete impurity;
  delete coulomb;
}


// ######################### gauss_HGH_atom ##########################
gauss_HGH_atom::gauss_HGH_atom(elements_t atom)
{
  set_semicore(atom);
}

void gauss_HGH_atom::set_semicore(elements_t atom)
{
  switch(atom)
    {
    case H:
      Zion = 1;
      rloc = 0.2;
      C1 = -4.180237;
      C2 = .725075;
      C3 = 0;
      C4 = 0;
      r0 = 0;
      h011 = 0;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case He:
      Zion = 2;
      rloc = 0.2;
      C1 = -9.112023;
      C2 = 1.698368;
      C3 = 0;
      C4 = 0;
      r0 = 0;
      h011 = 0;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Li:
      Zion = 3;
      rloc = .4;
      C1 = -14.034868;
      C2 = 9.553476;
      C3 = -1.766488;
      C4 = 0.08437;
      r0 = 0;
      h011 = 0;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Be:
      Zion = 4;
      rloc = .325;
      C1 = -24.015041;
      C2 = 17.204014;
      C3 = -3.32639;
      C4 = .165419;
      r0 = 0;
      h011 = 0;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case B:
      Zion = 3;
      rloc = .43393;
      C1 = -5.578642;
      C2 = .804251;
      C3 = 0;
      C4 = 0;
      r0 = .373843;
      h011 = 6.233928;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case C:
      Zion = 4;
      rloc = .34883;
      C1 = -8.513771;
      C2 = 1.228432;
      C3 = 0;
      C4 = 0;
      r0 = 0;
      h011 = .304553;
      h022 = 9.522842;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case N:
      Zion = 5;
      rloc = .289179;
      C1 = -12.234820;
      C2 = 1.766407;
      C3 = 0;
      C4 = 0;
      r0 = .256605;
      h011 = 13.552243;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case O:
      Zion = 6;
      rloc = .247621;
      C1 = -16.580318;
      C2 = 2.395701;
      C3 = 0;
      C4 = 0;
      r0 = .221786;
      h011 = 18.266917;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case F:
      Zion = 7;
      rloc = .218525;
      C1 = -21.307361;
      C2 = 3.072869;
      C3 = 0;
      C4 = 0;
      r0 = .195567;
      h011 = 23.584942;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Ne:
      Zion = 8;
      rloc = .19;
      C1 = -27.692852;
      C2 = 4.005906;
      C3 = 0;
      C4 = 0;
      r0 = .179488;
      h011 = 28.506098;
      h022 = -1.076245;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Na:
      Zion = 9;
      rloc = .246318;
      C1 = -7.545593;
      C2 = 1.125997;
      C3 = 0;
      C4 = 0;
      r0 = .141251;
      h011 = 36.556987;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Mg:
      Zion = 10;
      rloc = .210950;
      C1 = -19.419008;
      C2 = 2.871331;
      C3 = 0;
      C4 = 0;
      r0 = .141547;
      h011 = 40.316626;
      h022 = 0;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Al:
      Zion = 3;
      rloc = .45;
      C1 = -8.491351;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .460104;
      h011 = 5.08834;
      h022 = 2.6797;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Si:
      Zion = 4;
      rloc = .44;
      C1 = -7.336103;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .422738;
      h011 = 5.906928;
      h022 = 3.258196;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case P:
      Zion = 5;
      rloc = .43;
      C1 = -6.65422;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .389803;
      h011 = 6.842136;
      h022 = 3.856693;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case S:
      Zion = 6;
      rloc = .42;
      C1 = -6.554492;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .361757;
      h011 = 7.905303;
      h022 = 4.471698;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Cl:
      Zion = 7;
      rloc = .41;
      C1 = -6.864754;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .338208;
      h011 = 9.06224;
      h022 = 5.065682;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Ar:
      Zion = 8;
      rloc = .4;
      C1 = -7.1;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .317381;
      h011 = 10.249487;
      h022 = 5.602516;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case K:
      Zion = 9;
      rloc = .4;
      C1 = -4.989348;
      C2 = -.756048;
      C3 = 0;
      C4 = 0;
      r0 = .294826;
      h011 = 11.238705;
      h022 = 7.067779;
      h033 = 0;
      r2 = 0;
      h211 = 0;
      h222 = 0;
      h233 = 0;
      break;
    case Ca:
      Zion = 10;
      rloc = .39;
      C1 = -4.928146;
      C2 = -1.232854;
      C3 = 0;
      C4 = 0;
      r0 = .281909;
      h011 = 12.35234;
      h022 = 7.657455;
      h033 = 0;
      r2 = .90433;
      h211 = .016806;
      h222 = 0;
      h233 = 0;
      break;
    case Sc:
      Zion = 11;
      rloc = .385;
      C1 = 7.425036;
      C2 = -.489852;
      C3 = 0;
      C4 = 0;
      r0 = .359707;
      h011 = 6.119585;
      h022 = -2.563453;
      h033 = 0;
      r2 = .252945;
      h211 = -8.020892;
      h222 = 0;
      h233 = 0;
      break;
    case Ti:
      Zion = 12;
      rloc = .38;
      C1 = 7.548789;
      C2 = -.588377;
      C3 = 0;
      C4 = 0;
      r0 = .334235;
      h011 = 6.92574;
      h022 = -3.142005;
      h033 = 0;
      r2 = .242947;
      h211 = -9.125896;
      h222 = 0;
      h233 = 0;
      break;
    case V:
      Zion = 13;
      rloc = .375;
      C1 = 4.941291;
      C2 = -.096443;
      C3 = 0;
      C4 = 0;
      r0 = .326651;
      h011 = 7.65939;
      h022 = -3.892229;
      h033 = 0;
      r2 = .240792;
      h211 = -8.828518;
      h222 = 0;
      h233 = 0;
      break;
    case Cr:
      Zion = 14;
      rloc = .37;
      C1 = 5.113362;
      C2 = -.646819;
      C3 = 0;
      C4 = 0;
      r0 = .306011;
      h011 = 8.617835;
      h022 = -4.137695;
      h033 = 0;
      r2 = .219577;
      h211 = -11.157868;
      h222 = 0;
      h233 = 0;
      break;
    case Mn:
      Zion = 15;
      rloc = .365;
      C1 = 6.748683;
      C2 = -.576569;
      C3 = 0;
      C4 = 0;
      r0 = .280753;
      h011 = 9.379532;
      h022 = -5.57528;
      h033 = 0;
      r2 = .221422;
      h211 = -12.115385;
      h222 = 0;
      h233 = 0;
      break;
    case Fe:
      Zion = 16;
      rloc = .36;
      C1 = 5.392507;
      C2 = -.030066;
      C3 = 0;
      C4 = 0;
      r0 = .269268;
      h011 = 10.193723;
      h022 = -6.834982;
      h033 = 0;
      r2 = .223021;
      h211 = -12.026941;
      h222 = 0;
      h233 = 0;
      break;
    case Co:
      Zion = 17;
      rloc = .355;
      C1 = 3.418391;
      C2 = .482078;
      C3 = 0;
      C4 = 0;
      r0 = .25914;
      h011 = 11.195226;
      h022 = -7.845825;
      h033 = 0;
      r2 = .221665;
      h211 = -12.075354;
      h222 = 0;
      h233 = 0;
      break;
    case Ni:
      Zion = 18;
      rloc = .35;
      C1 = 3.610311;
      C2 = .449638;
      C3 = 0;
      C4 = 0;
      r0 = .245105;
      h011 = 12.161131;
      h022 = -9.078929;
      h033 = 0;
      r2 = .21495;
      h211 = -13.395062;
      h222 = 0;
      h233 = 0;
      break;
    case Cu:
      Zion = 11;
      rloc = .53;
      C1 = 0;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .423734;
      h011 = 3.88805;
      h022 = 3.276584;
      h033 = 2.290091;
      r2 = .266143;
      h211 = -12.676957;
      h222 = 0;
      h233 = 0;
      break;
    case Zn:
      Zion = 12;
      rloc = .51;
      C1 = 0;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .400866;
      h011 = 4.27871;
      h022 = 3.627342;
      h033 = 2.849567;
      r2 = .252151;
      h211 = -14.338368;
      h222 = 0;
      h233 = 0;
      break;
    case Ga:
      Zion = 13;
      rloc = .49;
      C1 = 0;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .384713;
      h011 = 4.831779;
      h022 = 4.238168;
      h033 = 2.833238;
      r2 = .240803;
      h211 = -15.795675;
      h222 = 0;
      h233 = 0;
      break;
    case Ge:
      Zion = 4;
      rloc = .54;
      C1 = 0;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .493743;
      h011 = 3.826891;
      h022 = 1.100231;
      h033 = -1.344218;
      r2 = .788369;
      h211 = .191205;
      h222 = 0;
      h233 = 0;
      break;
    case As:
      Zion = 5;
      rloc = .52;
      C1 = 0;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .4564;
      h011 = 4.560761;
      h022 = 1.692389;
      h033 = -1.373804;
      r2 = .685283;
      h211 = .312373;
      h222 = 0;
      h233 = 0;
      break;
    case Se:
      Zion = 6;
      rloc = .51;
      C1 = 0;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .432531;
      h011 = 5.145131;
      h022 = 2.052009;
      h033 = -1.369203;
      r2 = .61342;
      h211 = .434829;
      h222 = 0;
      h233 = 0;
      break;
    case Br:
      Zion = 7;
      rloc = .5;
      C1 = 0;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .428207;
      h011 = 5.398837;
      h022 = 1.820292;
      h033 = -1.323974;
      r2 = .557847;
      h211 = .555903;
      h222 = 0;
      h233 = 0;
      break;
    case Kr:
      Zion = 8;
      rloc = .5;
      C1 = 0;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .410759;
      h011 = 5.911194;
      h022 = 1.967372;
      h033 = -1.458069;
      r2 = .51712;
      h211 = .629228;
      h222 = 0;
      h233 = 0;
      break;
    case Rb:
      Zion = 9;
      rloc = .49;
      C1 = 4.504151;
      C2 = -.741018;
      C3 = 0;
      C4 = 0;
      r0 = .282301;
      h011 = 9.536329;
      h022 = 9.486634;
      h033 = 0;
      r2 = .514895;
      h211 = .449376;
      h222 = 0;
      h233 = 0;
      break;
    case Sr:
      Zion = 10;
      rloc = .48;
      C1 = 5.571455;
      C2 = -1.079963;
      C3 = 0;
      C4 = 0;
      r0 = .275441;
      h011 = 9.995135;
      h022 = 9.336679;
      h033 = 0;
      r2 = .502045;
      h211 = .43728;
      h222 = 0;
      h233 = 0;
      break;
    case Y:
      Zion = 11;
      rloc = .475;
      C1 = 13.217914;
      C2 = 1.353178;
      C3 = 0;
      C4 = 0;
      r0 = .41436;
      h011 = 2.522621;
      h022 = -4.363769;
      h033 = 0;
      r2 = .513304;
      h211 = -1.571003;
      h222 = .627616;
      h233 = 0;
      break;
    case Zr:
      Zion = 12;
      rloc = .47;
      C1 = 15.782342;
      C2 = .433648;
      C3 = 0;
      C4 = 0;
      r0 = .39654;
      h011 = 2.571767;
      h022 = -4.714509;
      h033 = 0;
      r2 = .520496;
      h211 = -1.548402;
      h222 = .826127;
      h233 = 0;
      break;
    case Nb:
      Zion = 13;
      rloc = .46;
      C1 = 13.505394;
      C2 = .752434;
      C3 = 0;
      C4 = 0;
      r0 = .393708;
      h011 = 3.222025;
      h022 = -4.599342;
      h033 = 0;
      r2 = .513644;
      h211 = -1.489848;
      h222 = .823817;
      h233 = 0;
      break;
    case Mo:
      Zion = 14;
      rloc = .43;
      C1 = 16.237452;
      C2 = 1.496536;
      C3 = 0;
      C4 = 0;
      r0 = .376255;
      h011 = 3.362426;
      h022 = -5.289276;
      h033 = 0;
      r2 = .525828;
      h211 = -1.543211;
      h222 = 1.074388;
      h233 = 0;
      break;
    case Tc:
      Zion = 15;
      rloc = .43;
      C1 = 14.910011;
      C2 = 1.046381;
      C3 = 0;
      C4 = 0;
      r0 = .369721;
      h011 = 3.917408;
      h022 = -5.268399;
      h033 = 0;
      r2 = .510487;
      h211 = -1.586709;
      h222 = 1.132307;
      h233 = 0;
      break;
    case Ru:
      Zion = 16;
      rloc = .43;
      C1 = 13.582571;
      C2 = .596227;
      C3 = 0;
      C4 = 0;
      r0 = .364084;
      h011 = 4.480632;
      h022 = -5.268679;
      h033 = 0;
      r2 = .49585;
      h211 = -1.59787;
      h222 = 1.165495;
      h233 = 0;
      break;
    case Rh:
      Zion = 17;
      rloc = .42;
      C1 = 15.225012;
      C2 = .415911;
      C3 = 0;
      C4 = 0;
      r0 = .350052;
      h011 = 4.715292;
      h022 = -5.805525;
      h033 = 0;
      r2 = .49695;
      h211 = -1.685594;
      h222 = 1.387707;
      h233 = 0;
      break;
    case Pd:
      Zion = 18;
      rloc = .41;
      C1 = 15.720259;
      C2 = .140765;
      C3 = 0;
      C4 = 0;
      r0 = .342151;
      h011 = 5.177686;
      h022 = -5.852819;
      h033 = 0;
      r2 = .494916;
      h211 = -1.608273;
      h222 = 1.446609;
      h233 = 0;
      break;
    case Ag:
      Zion = 11;
      rloc = .57;
      C1 = 1.017053;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .4989;
      h011 = 2.990284;
      h022 = 3.912395;
      h033 = 2.205847;
      r2 = .38766;
      h211 = -3.420076;
      h222 = -1.019949;
      h233 = 0;
      break;
    case Cd:
      Zion = 12;
      rloc = .55;
      C1 = 2.382713;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .491505;
      h011 = 3.207932;
      h022 = 4.140963;
      h033 = 1.584234;
      r2 = .377874;
      h211 = -4.190072;
      h222 = -.770156;
      h233 = 0;
      break;
    case In:
      Zion = 13;
      rloc = .53;
      C1 = 2.395404;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .474081;
      h011 = 3.554411;
      h022 = 4.754135;
      h033 = 1.565040;
      r2 = .360488;
      h211 = -4.566414;
      h222 = -.773785;
      h233 = 0;
      break;
    case Sn:
      Zion = 4;
      rloc = .605;
      C1 = 4.610912;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .663544;
      h011 = 1.648791;
      h022 = -.141974;
      h033 = -.576546;
      r2 = .944459;
      h211 = .225115;
      h222 = 0;
      h233 = 0;
      break;
    case Sb:
      Zion = 5;
      rloc = .59;
      C1 = 6.680228;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .597684;
      h011 = 1.951477;
      h022 = .037537;
      h033 = -.786631;
      r2 = .856557;
      h211 = .300103;
      h222 = 0;
      h233 = 0;
      break;
    case Te:
      Zion = 6;
      rloc = .575;
      C1 = 9.387085;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .556456;
      h011 = 2.04689;
      h022 = -.029333;
      h033 = -.881119;
      r2 = .805101;
      h211 = .317411;
      h222 = 0;
      h233 = 0;
      break;
    case I:
      Zion = 7;
      rloc = .56;
      C1 = 14.661825;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .55283;
      h011 = 1.338054;
      h022 = -.834851;
      h033 = -.467438;
      r2 = .794325;
      h211 = .224345;
      h222 = 0;
      h233 = 0;
      break;
    case Xe:
      Zion = 8;
      rloc = .56;
      C1 = 12.73428;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .507371;
      h011 = 2.236451;
      h022 = -.403551;
      h033 = -1.132507;
      r2 = .729821;
      h211 = .280131;
      h222 = 0;
      h233 = 0;
      break;
    case Cs:
      Zion = 9;
      rloc = .54;
      C1 = 35.234438;
      C2 = -3.31807;
      C3 = 0;
      C4 = 0;
      r0 = .456821;
      h011 = -.282378;
      h022 = -2.780956;
      h033 = 0;
      r2 = .761462;
      h211 = .183754;
      h222 = 0;
      h233 = 0;
      break;
    case Ba:
      Zion = 10;
      rloc = .54;
      C1 = 24.478653;
      C2 = -2.50085;
      C3 = 0;
      C4 = 0;
      r0 = .514776;
      h011 = 1.046729;
      h022 = -.977217;
      h033 = 0;
      r2 = .665403;
      h211 = .378419;
      h222 = 0;
      h233 = 0;
      break;
    case La:
      Zion = 11;
      rloc = .535;
      C1 = 19.909308;
      C2 = -1.47483;
      C3 = 0;
      C4 = 0;
      r0 = .551775;
      h011 = 1.293272;
      h022 = -1.121819;
      h033 = 0;
      r2 = .626672;
      h211 = .328377;
      h222 = 0;
      h233 = 0;
      break;
    case Ce:
      Zion = 12;
      rloc = .535;
      C1 = 18.84747;
      C2 = -.765636;
      C3 = 0;
      C4 = 0;
      r0 = .52179;
      h011 = 1.321616;
      h022 = -1.700444;
      h033 = 0;
      r2 = .703593;
      h211 = .074241;
      h222 = 0;
      h233 = 0;
      break;
    case Pr:
      Zion = 13;
      rloc = .532083;
      C1 = 18.424739;
      C2 = -.657669;
      C3 = 0;
      C4 = 0;
      r0 = .52685;
      h011 = 1.012621;
      h022 = -1.717982;
      h033 = 0;
      r2 = .74761;
      h211 = .017571;
      h222 = 0;
      h233 = 0;
      break;
    case Nd:
      Zion = 14;
      rloc = .529167;
      C1 = 17.81503;
      C2 = -.594798;
      C3 = 0;
      C4 = 0;
      r0 = .503;
      h011 = 1.52911;
      h022 = -2.153732;
      h033 = 0;
      r2 = .32529;
      h211 = -.54324;
      h222 = 0;
      h233 = 0;
      break;
    case Pm:
      Zion = 15;
      rloc = .52625;
      C1 = 18.251723;
      C2 = -.492107;
      C3 = 0;
      C4 = 0;
      r0 = .489879;
      h011 = 1.308978;
      h022 = -2.507751;
      h033 = 0;
      r2 = 0.473709;
      h211 = -.429952;
      h222 = 0;
      h233 = 0;
      break;
    case Sm:
      Zion = 16;
      rloc = .523333;
      C1 = 17.206792;
      C2 = -.532803;
      C3 = 0;
      C4 = 0;
      r0 = .479677;
      h011 = 1.723635;
      h022 = -2.659367;
      h033 = 0;
      r2 = .47084;
      h211 = -.41063;
      h222 = 0;
      h233 = 0;
      break;
    case Eu:
      Zion = 17;
      rloc = .520417;
      C1 = 17.373516;
      C2 = -.648468;
      C3 = 0;
      C4 = 0;
      r0 = .469043;
      h011 = 1.763638;
      h022 = -2.916932;
      h033 = 0;
      r2 = .490038;
      h211 = -.42612;
      h222 = 0;
      h233 = 0;
      break;
    case Gd:
      Zion = 18;
      rloc = .5175;
      C1 = 17.512556;
      C2 = -.719534;
      C3 = 0;
      C4 = 0;
      r0 = .462014;
      h011 = 1.551856;
      h022 = -3.068703;
      h033 = 0;
      r2 = .482368;
      h211 = -.562601;
      h222 = 0;
      h233 = 0;
      break;
    case Tb:
      Zion = 19;
      rloc = .514583;
      C1 = 17.603616;
      C2 = -.82808;
      C3 = 0;
      C4 = 0;
      r0 = .448694;
      h011 = 1.718481;
      h022 = -3.435239;
      h033 = 0;
      r2 = .482809;
      h211 = -.625802;
      h222 = 0;
      h233 = 0;
      break;
    case Dy:
      Zion = 20;
      rloc = .511667;
      C1 = 16.994331;
      C2 = -.955298;
      C3 = 0;
      C4 = 0;
      r0 = .44059;
      h011 = 1.94032;
      h022 = -3.559798;
      h033 = 0;
      r2 = .467229;
      h211 = -.668924;
      h222 = 0;
      h233 = 0;
      break;
    case Ho:
      Zion = 21;
      rloc = .50875;
      C1 = 16.78157;
      C2 = -1.173514;
      C3 = 0;
      C4 = 0;
      r0 = .432212;
      h011 = 2.052797;
      h022 = -3.674534;
      h033 = 0;
      r2 = .447131;
      h211 = -.742863;
      h222 = 0;
      h233 = 0;
      break;
    case Er:
      Zion = 22;
      rloc = .505833;
      C1 = 17.105293;
      C2 = -1.430953;
      C3 = 0;
      C4 = 0;
      r0 = .419948;
      h011 = 2.144503;
      h022 = -3.984203;
      h033 = 0;
      r2 = .418385;
      h211 = -.999006;
      h222 = 0;
      h233 = 0;
      break;
    case Tm:
      Zion = 23;
      rloc = .502917;
      C1 = 17.247293;
      C2 = -1.627697;
      C3 = 0;
      C4 = 0;
      r0 = .413373;
      h011 = 1.947196;
      h022 = -4.121556;
      h033 = 0;
      r2 = .39287;
      h211 = -1.353308;
      h222 = 0;
      h233 = 0;
      break;
    case Yb:
      Zion = 24;
      rloc = .5;
      C1 = 17.357144;
      C2 = -1.773916;
      C3 = 0;
      C4 = 0;
      r0 = .402309;
      h011 = 2.120771;
      h022 = -4.80299;
      h033 = 0;
      r2 = .444025;
      h211 = -.889967;
      h222 = 0;
      h233 = 0;
      break;
    case Lu:
      Zion = 25;
      rloc = .497;
      C1 = 17.037053;
      C2 = -1.66161;
      C3 = 0;
      C4 = 0;
      r0 = .391206;
      h011 = 2.184678;
      h022 = -5.432346;
      h033 = 0;
      r2 = .436518;
      h211 = -1.173245;
      h222 = 0;
      h233 = 0;
      break;
    case Hf:
      Zion = 12;
      rloc = .56;
      C1 = 5.134801;
      C2 = .529191;
      C3 = 0;
      C4 = 0;
      r0 = .42281;
      h011 = 2.564442;
      h022 = -6.013732;
      h033 = 1.10976;
      r2 = .426388;
      h211 = 1.459363;
      h222 = -5.282764;
      h233 = 0;
      break;
    case Ta:
      Zion = 13;
      rloc = .55;
      C1 = 4.546236;
      C2 = .779422;
      C3 = 0;
      C4 = 0;
      r0 = .421853;
      h011 = 2.708136;
      h022 = -5.790959;
      h033 = .947663;
      r2 = .410994;
      h211 = 1.348495;
      h222 = -5.386947;
      h233 = 0;
      break;
    case W:
      Zion = 14;
      rloc = .54;
      C1 = 4.800251;
      C2 = .901544;
      C3 = 0;
      C4 = 0;
      r0 = .41857;
      h011 = 2.692204;
      h022 = -6.022637;
      h033 = 1.218316;
      r2 = .399602;
      h211 = 1.177436;
      h222 = -5.553621;
      h233 = 0;
      break;
    case Re:
      Zion = 15;
      rloc = .53;
      C1 = 5.59266;
      C2 = .943957;
      C3 = .403252;
      C4 = 2.76072;
      r0 = -6.396415;
      h011 = .868732;
      h022 = 0;
      h033 = 0;
      r2 = .390395;
      h211 = .875251;
      h222 = -5.672543;
      h233 = 0;
      break;
    case Os:
      Zion = 16;
      rloc = .52;
      C1 = 5.613073;
      C2 = .921955;
      C3 = 0;
      C4 = 0;
      r0 = .410578;
      h011 = 2.785758;
      h022 = -6.69213;
      h033 = 2.247034;
      r2 = .380252;
      h211 = .880133;
      h222 = -5.732892;
      h233 = 0;
      break;
    case Ir:
      Zion = 17;
      rloc = .51;
      C1 = 4.904509;
      C2 = 1.313786;
      C3 = 0;
      C4 = 0;
      r0 = .404469;
      h011 = 3.243278;
      h022 = -7.315509;
      h033 = 2.956978;
      r2 = .376428;
      h211 = .754315;
      h222 = -5.87558;
      h233 = 0;
      break;
    case Pt:
      Zion = 18;
      rloc = .5;
      C1 = 5.445832;
      C2 = 1.156382;
      C3 = 0;
      C4 = 0;
      r0 = .409942;
      h011 = 2.994366;
      h022 = -7.448772;
      h033 = 4.243095;
      r2 = .367964;
      h211 = .632067;
      h222 = -5.755431;
      h233 = 0;
      break;
    case Au:
      Zion = 11;
      rloc = .59;
      C1 = 11.604428;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .52118;
      h011 = 2.538614;
      h022 = 2.701113;
      h033 = 0;
      r2 = .440706;
      h211 = -4.71907;
      h222 = -1.650429;
      h233 = 0;
      break;
    case Hg:
      Zion = 12;
      rloc = .57;
      C1 = 2.134572;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .521802;
      h011 = 3.29392;
      h022 = 4.661001;
      h033 = 0;
      r2 = .401894;
      h211 = -1.669886;
      h222 = -2.473265;
      h233 = 0;
      break;
    case Tl:
      Zion = 13;
      rloc = .55;
      C1 = 7.301886;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .502423;
      h011 = 3.32656;
      h022 = 4.34139;
      h033 = 0;
      r2 = .393185;
      h211 = -3.200652;
      h222 = -3.008296;
      h233 = 0;
      break;
    case Pb:
      Zion = 4;
      rloc = .6175;
      C1 = .753143;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .705259;
      h011 = 1.979927;
      h022 = -.16496;
      h033 = -.80606;
      r2 = .971939;
      h211 = .374967;
      h222 = 0;
      h233 = 0;
      break;
    case Bi:
      Zion = 5;
      rloc = .605;
      C1 = 6.679437;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .678858;
      h011 = 1.377634;
      h022 = -.513697;
      h033 = -.471028;
      r2 = .934683;
      h211 = .378476;
      h222 = 0;
      h233 = 0;
      break;
    case Po:
      Zion = 6;
      rloc = .5925;
      C1 = 10.411731;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .64795;
      h011 = 1.144203;
      h022 = -.735851;
      h033 = -.339386;
      r2 = .880468;
      h211 = .433232;
      h222 = 0;
      h233 = 0;
      break;
    case At:
      Zion = 7;
      rloc = .58;
      C1 = 13.520411;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .627827;
      h011 = .945557;
      h022 = -.965903;
      h033 = -.190429;
      r2 = .838365;
      h211 = .468948;
      h222 = 0;
      h233 = 0;
      break;
    case Rn:
      Zion = 8;
      rloc = .57;
      C1 = 14.629185;
      C2 = 0;
      C3 = 0;
      C4 = 0;
      r0 = .615182;
      h011 = .981832;
      h022 = -1.038963;
      h033 = -.120456;
      r2 = .788337;
      h211 = .557746;
      h222 = 0;
      h233 = 0;
      break;
    default:
      printf("Error: HGH pseudopotential not available for %s\n",
	     elements_to_string(atom));
    }
}

gsl_matrix_complex *gauss_HGH_atom::matrix_block(double a1, double a2)
{
  double h012 = - 0.5 * sqrt(3.0/5) * h022,
    h013 = .5 * sqrt(5.0/21) * h033,
    h023 = - 5.0 / 3 * sqrt(1.0/7) * h033,
    h212 = - 1.0 / 6 * sqrt(7) * h222,
    h213 = 1.5 * sqrt(7/143) * h233,
    h223 = - 9.0 / sqrt(143) * h233;
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_hgh_def.hh"
  return output;
}


// ########################## gauss_overlap ##########################
gsl_matrix_complex *gauss_overlap::matrix_block(double a1, double a2)
{
  gsl_matrix_complex *output =
    gsl_matrix_complex_alloc(matrix_term_block_size, matrix_term_block_size);
  #include "gauss_overlap_def.hh"
  return output;
}
