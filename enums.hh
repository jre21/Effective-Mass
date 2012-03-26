#ifndef _ENUMS_HH
#define _ENUMS_HH 1

enum elements_t
  {
    H,He,
    Li,Be, B,C,N,O,F,Ne,
    Na,Mg, Al,Si,P,S,Cl,Ar,
    K,Ca,  Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn, Ga,Ge,As,Se,Br,Kr,
    Rb,Sr, Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd, In,Sn,Sb,Te,I,Xe,
    Cs,Ba, La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,
           Lu,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg, Tl,Pb,Bi,Po,At,Rn,
    Fr,Ra, Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,
           Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn
  };

const char *elements_to_string(elements_t e);

enum crystals_t
  {
    GaN, AlN, InN, ZnGeN2, ZnSnN2
  };

const char *crystals_to_string(crystals_t c);

enum crystal_parameters_t
  {
    g1, g2, g3, d0,
    A1, A2, A3, A4, A5, A6, d1, d2, d3,
    B1, B2, B3, C1, C2, C3, D1, D2, D3, d1c, d2c, d1so, d2so, d3so
  };

enum impurity_parameters_t
  {
    dielectric_ratio
  };

#endif // _ENUMS_HH
