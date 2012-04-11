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
           Lr,Rf,Db,Sg,Bh,Hs,Mt,Ds,Rg,Cn,
    // elements with 3d band treated as valance
    Ga_val, Ge_val, As_val, Se_val, Br_val
  };

const char *elements_to_string(elements_t e);

enum crystals_t
  {
    GaN, AlN, InN, ZnGeN2, ZnSnN2
  };

const char *crystals_to_string(crystals_t c);

enum crystal_parameters_t
  {
    _g1, _g2, _g3, _d0,
    _A1, _A2, _A3, _A4, _A5, _A6, _d1, _d2, _d3,
    _B1, _B2, _B3, _C1, _C2, _C3, _D1, _D2, _D3,
    _d1c, _d2c, _d1so, _d2so, _d3so,
    _dielectric
  };

const char *crystal_parameters_to_string(crystal_parameters_t c);

enum impurity_parameters_t
  {
    _dielectric_ratio
  };

const char *impurity_parameters_to_string(impurity_parameters_t c);

#endif // _ENUMS_HH
