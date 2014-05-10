//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>
#include <stdio.h>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../rsolvers/riemann_solver.hpp"
#include "fluid_integrator.hpp"

#define SUM_ON 1

//======================================================================================
/*! \file van_leer2.cpp
 *  \brief van-Leer (MUSCL-Hancock) second-order integrator
 *====================================================================================*/

//--------------------------------------------------------------------------------------
///*! \fn  void FluidIntegrator::PredictVanLeer2
// *  \brief predictor step for 2nd order VL integrator */

void FluidIntegrator::Predict(Block *pb)
{
  int is = pb->is; int js = pb->js; int ks = pb->ks;
  int ie = pb->ie; int je = pb->je; int ke = pb->ke;

  Real sum=0.0; Real sum_2=0.0;
  int ndata = (pb->block_size.nx1)*(pb->block_size.nx2)*(pb->block_size.nx3)*(NVAR);
 
  AthenaArray<Real> u = pb->pfluid->u;
  AthenaArray<Real> w = pb->pfluid->w;
  AthenaArray<Real> u1 = pb->pfluid->u1;
  AthenaArray<Real> w1 = pb->pfluid->w1;

  AthenaArray<Real> wl = pb->pfluid->wl_;
  AthenaArray<Real> wr = pb->pfluid->wr_;
  AthenaArray<Real> flx = pb->pfluid->flx_;
 
#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
    sum += pb->pfluid->u(i);
    sum_2 += pb->pfluid->w(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif

//--------------------------------------------------------------------------------------
// i-direction 

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je; ++j){

    lr_states_->PiecewiseLinear(k,j,is,ie+1,1,w,wl,wr);

    flux_->HLLC(is,ie+1,wl,wr,flx);

    for (int n=0; n<NVAR; ++n){
#pragma simd
      for (int i=is; i<=ie; ++i){
        Real& ui  = u (n,k,j,i);
        Real& u1i = u1(n,k,j,i);
        Real& flxi   = flx(n,  i);
        Real& flxip1 = flx(n,i+1);
        Real& dx = pb->dx1f(i);
 
        u1i = ui - 0.5*pb->pfluid->dt*(flxip1 - flxi)/dx;
      }
    }

  }}

#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
      sum += pb->pfluid->u(i);
      sum_2 += u1(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif


//--------------------------------------------------------------------------------------
// j-direction

  for (int k=ks; k<=ke; ++k){
  for (int j=js; j<=je+1; ++j){

    lr_states_->PiecewiseLinear(k,j,is,ie,2,w,wl,wr);

    flux_->HLLC(is,ie,wl,wr,flx); 

    Real dtodxjm1 = 0.5*pb->pfluid->dt/(pb->dx2f(j-1));
    Real dtodxj   = 0.5*pb->pfluid->dt/(pb->dx2f(j));
    for (int n=0; n<NVAR; ++n){
#pragma simd
      for (int i=is; i<=ie; ++i){
        Real& u1jm1 = u1(n,k,j-1,i);
        Real& u1j   = u1(n,k,j  ,i);
        Real& flxj  = flx(n,i);
  
        u1jm1 -= dtodxjm1*flxj;
        u1j   += dtodxj  *flxj;
      }
    }

  }}


#if SUM_ON>0
  /**** output to force calcs ****/
  for (int i=0; i<ndata; ++i){
    sum += pb->pfluid->u(i);
    sum_2 += u1(i);
  }
  printf("sum = %e sum_2 = %e\n", sum, sum_2);
#endif



//--------------------------------------------------------------------------------------
// k-direction 

  for (int k=ks; k<=ke+1; ++k){
  for (int j=js; j<=je; ++j){

    lr_states_->PiecewiseLinear(k,j,is,ie,3,w,wl,wr);

    flux_->HLLC(is,ie,wl,wr,flx);

    Real dtodxkm1 = 0.5*pb->pfluid->dt/(pb->dx3f(k-1));
    Real dtodxk   = 0.5*pb->pfluid->dt/(pb->dx3f(k));
    for (int n=0; n<NVAR; ++n){
#pragma simd
      for (int i=is; i<=ie; ++i){
        Real& u1km1 = u1(n,k-1,j,i);
        Real& u1k   = u1(n,k  ,j,i);
        Real& flxk = flx(n,i);

        u1km1 -= dtodxkm1*flxk;
        u1k   += dtodxk  *flxk;
      }
    }

  }}

  return;
}
