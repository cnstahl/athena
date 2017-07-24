//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file thermal_diff.cpp
//  \brief Problem generator to test thermal diffusion
//
// Sets up two different problems:
//   - iprob=1: Decaying Gaussian profile
//   - iprob=2: 
//========================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Gaussian profile
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Read problem parameters
  int iprob = pin->GetInteger("problem","iprob");
  Real chi = pin->GetReal("problem","chiiso");

//--- iprob=1.  Gaussian profile

  if (iprob == 1) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IDN,k,j,i) = 2.0;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
      phydro->u(IEN,k,j,i) = exp(-(SQR(pcoord->x1v(i)))/(4.0*chi*2/5));
    }}}
  }

//--- iprob=2. Sound wave

  if (iprob == 2) {
    Real a = 0.05;
  }

  return;
}
