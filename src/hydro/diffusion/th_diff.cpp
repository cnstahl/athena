//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena++ headers
// TODO: clean includes
#include "diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"
#include "../../eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroDiffusion::ThDiff
//  \brief Adds thermal diffusion as flux to the hydro flux

void HydroDiffusion::ThDiff(const AthenaArray<Real> &prim,
    const AthenaArray<Real> &cons, AthenaArray<Real> *diflx)
{
  AthenaArray<Real> &x1flux=diflx[X1DIR];
  AthenaArray<Real> &x2flux=diflx[X2DIR];
  AthenaArray<Real> &x3flux=diflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real denf, dTdx, dTdy, dTdz;

  // TODO: support EOS

  // step-1. calculate the flux due to thermal difusion
  // i-direction
  // set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) {
      if(pmb_->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }
  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
      // compute and store fluxes
      for (int i=is; i<=ie+1; ++i){
        // if viscosity is set, fluxes are updated there. Otherwise they need
        // to be updated here
        if (!hydro_diffusion_defined) x1flux(IEN,k,j,i) = 0;
		    Real chi1 = chiiso1();
        denf = 0.5*(prim(IDN,k,j,i)+prim(IDN,k,j,i-1));
        dTdx = (prim(IPR,k,j,i)/prim(IDN,k,j,i) - prim(IPR,k,j,i-1)/
                  prim(IDN,k,j,i-1))/pco_->dx1v(i-1);
        // TODO: minus sign?
        x1flux(IEN,k,j,i) += -denf * chi1 * dTdx;
      }
    }
  }

  // j-direction
  // set the loop limits
  il=is, iu=ie, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx3 == 1) // 2D
      il=is-1, iu=ie+1, kl=ks, ku=ke;
    else // 3D
      il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
  }
  if(pmb_->block_size.nx2 > 1) { //2D or 3D
    for (int k=kl; k<=ku; ++k){
      for (int j=js; j<=je+1; ++j){
        // compute fluxes
        for (int i=il; i<=iu; ++i){
          if (!hydro_diffusion_defined) x1flux(IEN,k,j,i) = 0;
          dTdy = (prim(IPR,k,j,i)/prim(IDN,k,j,i)-prim(IPR,k,j-1,i)/
                    prim(IDN,k,j-1,i))/pco_->h2v(i)/pco_->dx2v(j-1);
		      Real chi1 = chiiso1();
          denf = 0.5*(prim(IDN,k,j+1,i)+prim(IDN,k,j,i));
          x2flux(IEN,k,j,i) += -denf * chi1 * dTdy;
        }
      }
    }
  }

  // k-direction
  // set the loop limits
  il=is, iu=ie, jl=js, ju=je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if(pmb_->block_size.nx2 > 1) // 2D or 3D
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;
    else // 1D
      il=is-1, iu=ie+1;
  }
  if(pmb_->block_size.nx3 > 1) { //3D
    for (int k=ks; k<=ke+1; ++k){
      for (int j=jl; j<=ju; ++j){
        // compute fluxes
        for (int i=il; i<=iu; ++i){
          if (!hydro_diffusion_defined) x1flux(IEN,k,j,i) = 0;
          dTdz = (prim(IPR,k,j,i)/prim(IDN,k,j,i)-prim(IPR,k-1,j,i)/
                   prim(IDN,k-1,j,i))/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j);
	        Real chi1 = chiiso1();
          denf = 0.5*(prim(IDN,k+1,j,i)+prim(IDN,k,j,i));
          x3flux(IEN,k,j,i) += -denf * chi1 * dTdz;

        }
      }
    }
  }
  
  return;
}

//-------------------------------------------------------------------------------------
// Get the coefficient chiiso1
Real HydroDiffusion::chiiso1(void)
{
  return (chiiso_);
}
