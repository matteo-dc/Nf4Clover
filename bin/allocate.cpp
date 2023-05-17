#include "allocate.hpp"
#include "Zq.hpp"
#include "aliases.hpp"
#include "global.hpp"
#include "operations.hpp"
#include "sigmas.hpp"
#include "vertices.hpp"

// read file
void allocate_vec_internal(double &t, const vector<int> sizes, int isize) { t = 0.0; }

void allocate_vec_internal(oper_t &o, const vector<int> sizes, int isize) {}

void oper_t::clear_all()
{
  sigma.clear();
  sigma.shrink_to_fit();

  jG.clear();
  jG.shrink_to_fit();

  jZq.clear();
  jZq.shrink_to_fit();

  jZ.clear();
  jZ.shrink_to_fit();

  jZVoverZA.clear();
  jZVoverZA.shrink_to_fit();
  jZPoverZS.clear();
  jZPoverZS.shrink_to_fit();
  jZAoverZV.clear();
  jZAoverZV.shrink_to_fit();

  bval.clear();
  bval.shrink_to_fit();
  bsea.clear();
  bsea.shrink_to_fit();

  stepfunc.clear();
  stepfunc.shrink_to_fit();
}

void oper_t::allocate_val()
{
  // allocate effective masses
  allocate_vec(eff_mass, {njacks, _nm, _nm});
  allocate_vec(eff_mass_sea, {njacks});
}

void oper_t::check_allocation()
{
  printf("sigma: \t\t %zu/%d %zu/%d %zu/%d %zu/%d %zu/%d \n", 
         sigma.size(), _linmoms,
         sigma[0].size(), sigma::nproj,
         sigma[0][0].size(), sigma::nins,
         sigma[0][0][0].size(), njacks,
         sigma[0][0][0][0].size(), _nmr);
         
  printf("jG: \t\t %zu/%d %zu/%d %zu/%d %zu/%d %zu/%d \n", 
         jG.size(), _bilmoms, 
         jG[0].size(), gbil::nins,
         jG[0][0].size(), nbil,
         jG[0][0][0].size(), njacks,
         jG[0][0][0][0].size(), _nmr);
}

void oper_t::allocate()
{
  clear_all();

  printf("\n--- Allocating memory --- \n");

  // allocate sigmas
  allocate_vec(sigma, {_linmoms, sigma::nproj, sigma::nins, njacks, _nmr});
  // allocate pr_bil
  allocate_vec(jG, {_bilmoms, gbil::nins, nbil, njacks, _nmr, _nmr});

  // allocate Zq
  allocate_vec(jZq, {_linmoms, njacks, _nmr});
  // allocate Zbil
  allocate_vec(jZ, {_bilmoms, nbil, njacks, _nmr, _nmr});

  // allocate ZV/ZA
  allocate_vec(jZVoverZA, {_bilmoms, 1, njacks, _nmr, _nmr});
  // allocate ZP/ZS
  allocate_vec(jZPoverZS, {_bilmoms, 1, njacks, _nmr, _nmr});
  // allocate ZA/ZV
  allocate_vec(jZAoverZV, {_bilmoms, 1, njacks, _nmr, _nmr});

  // allocate bval/bsea
  allocate_vec(bval, {_linmoms, njacks, nbil + 1});
  allocate_vec(bsea, {_linmoms, njacks, nbil + 1});

  // allocate stepfunc
  allocate_vec(stepfunc, {njacks, nbil + 1 + 3}); // nbil+zq+zv/za+zs/zp+za/zv

  check_allocation();
  printf("------------------------ \n");
}
