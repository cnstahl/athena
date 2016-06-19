//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
//======================================================================================
//! \file pointmass.cpp
//  \brief Adds source terms due to point mass AT ORIGIN
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"

// this class header
#include "hydro_srcterms.hpp"

//Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
//
//[JMSHI
//--------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::ShearingBoxSourceTerms(const Real dt,
//  const AthenaArray<Real> *flux, const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
//  \brief Shearing Box source terms
//
//  Detailed description starts here.
//  We add shearing box source term via operator splitting method. The source terms are
//  added after the fluxes are computed in each level of the integration (in
//  FluxDivergence) to give predictions of the conservative variables for either the next
//  level or the final update.
//  Currently, it is hard-wired to Van-Leer integration and in 2d x-y plane.
//  Note:
//  (1) when update CONS, should use area and vol to be consistent with other coordinates?
//  (2)

void HydroSourceTerms::ShearingBoxSourceTerms(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{
  if (Omega_0_==0.0 || qshear_==0.0 ) {
      std::cout << "[ShearingBoxSourceTerms]: Omega_0 or qshear not stated " << std::endl;
	  return;
  }
  Real phic,phil,phir;

  MeshBlock *pmb = pmy_hydro_->pmy_block;

  if (pmb->block_size.nx3 > 1) {
    std::cout << "[ShearingBoxSourceTerms]: not compatible to 3D yet!!" << std::endl;
	return;
  }
  else if (pmb->block_size.nx2 > 1) {
	int k = pmb->ks;
	if (ShBoxCoord_== 1) {
//#pragma omp parallel for schedule(static)
      for (int j=pmb->js; j<=pmb->je; ++j) {
//#pragma simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den = prim(IDN,k,j,i);
          phic = UnstratifiedDisk(pmb->pcoord->x1v(i),
		  						pmb->pcoord->x2v(j),
		  						pmb->pcoord->x3v(k));
          phil = UnstratifiedDisk(pmb->pcoord->x1f(i),
		  						pmb->pcoord->x2v(j),
		  						pmb->pcoord->x3v(k));
          phir = UnstratifiedDisk(pmb->pcoord->x1f(i+1),
		  						pmb->pcoord->x2v(j),
		  						pmb->pcoord->x3v(k));
// 1) S_M = -rho*grad(Phi); S_E = -rho*v*grad(Phi)
//    dM1/dt = 2q\rho\Omega^2 x
//    dE /dt = 2q\Omega^2 (\rho v_x)
// 2) Coriolis forces:
//    dM1/dt = 2\Omega(\rho v_y)
//    dM2/dt = -2\Omega(\rho v_x)
//

//		  cons(IM1,k,j,i) -= dt*((phir-phil)*den/pmb->pcoord->dx1v(i) -
//			                   2.0*Omega_0_*den*prim(IVY,k,j,i));
		  cons(IM1,k,j,i) += dt*(2.0*qshear_*Omega_0_*Omega_0_*den*pmb->pcoord->x1v(i) +
		  	                   2.0*Omega_0_*den*prim(IVY,k,j,i));
		  cons(IM2,k,j,i) -= dt*2.0*Omega_0_*den*prim(IVX,k,j,i);
          if (NON_BAROTROPIC_EOS) {
		  	  cons(IEN,k,j,i) -= dt*(flux[X1DIR](IDN,k,j,i)*(phic-phil) +
			 					     flux[X1DIR](IDN,k,j,i+1)*(phir-phic))
			                        /pmb->pcoord->dx1v(i);
		  }
        }
      }
	} else if (ShBoxCoord_ == 2) {
//#pragma omp parallel for schedule(static)
        for (int j=pmb->js; j<=pmb->je; ++j) {
//#pragma simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real den = prim(IDN,k,j,i);
            phic = UnstratifiedDisk(pmb->pcoord->x1v(i),
		      					pmb->pcoord->x2v(j),
		      					pmb->pcoord->x3v(k));
            phil = UnstratifiedDisk(pmb->pcoord->x1f(i),
		      					pmb->pcoord->x2v(j),
		      					pmb->pcoord->x3v(k));
            phir = UnstratifiedDisk(pmb->pcoord->x1f(i+1),
		      					pmb->pcoord->x2v(j),
		      					pmb->pcoord->x3v(k));
		    cons(IM1,k,j,i) += dt*(2.0*qshear_*Omega_0_*Omega_0_*den*pmb->pcoord->x1v(i) +
		                         2.0*Omega_0_*den*prim(IVZ,k,j,i));
		    cons(IM3,k,j,i) -= dt*2.0*Omega_0_*den*prim(IVX,k,j,i);
            if (NON_BAROTROPIC_EOS) {
		        cons(IEN,k,j,i) -= dt*(flux[X1DIR](IDN,k,j,i)*(phic-phil) +
		      				     flux[X1DIR](IDN,k,j,i+1)*(phir-phic))
		                          /pmb->pcoord->dx1v(i);
		  }
        }
      }
    }
  }
  else {
	std::cout << "[ShearingBoxSourceTerms]: not compatible to 1D !!" << std::endl;
	return;
  }


  return;
}


//--------------------------------------------------------------------------------------
//! \fn Real HydroSourceTerms::UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
//  \brief shearing box tidal potential
//
//  Detailed description starts here.
//

Real HydroSourceTerms::UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{

  Real phi=0.0;
  phi -= qshear_*Omega_0_*Omega_0_*x1*x1;
  return phi;
}
//JMSHI]
//--------------------------------------------------------------------------------------
