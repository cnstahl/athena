// Local Lax-Friedrichs Riemann solver for relativistic hydrodynamics

// Primary header
#include "../../hydro_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../hydro.hpp"                       // Hydro
#include "../../../eos/eos.hpp"                     // HydroEqnOfState
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../mesh.hpp"                     // MeshBlock
#include "../../../../coordinates/coordinates.hpp"  // Coordinates

// Declarations
static void NonTransformingLLF(const int k, const int j, const int i,
    const AthenaArray<Real> &prim_l, const AthenaArray<Real> &prim_r,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi, Hydro *pmy_hydro,
    AthenaArray<Real> &flux);

//--------------------------------------------------------------------------------------

// Riemann solver
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields (not used)
//   prim_l, prim_r: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   prim_l, prim_r overwritten
//   implements LLF algorithm similar to that of fluxcalc() in step_ch.c in Harm
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB)
void HydroIntegrator::RiemannSolver(const int k, const int j, const int il,
    const int iu, const int ivx, const AthenaArray<Real> &bb, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &flux)
{
  // Transform primitives to locally flat coordinates if in GR
  bool pole_top = false, pole_bottom = false;
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_hydro->pmy_block->pcoord->PrimToLocal1(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
      case IVY:
        pmy_hydro->pmy_block->pcoord->IdentifyPoles(&pole_top, &pole_bottom);
        pmy_hydro->pmy_block->pcoord->PrimToLocal2(k, j, il+(pole_top?1:0),
            iu-(pole_bottom?1:0), bb, prim_l, prim_r, bb_normal_);
        break;
      case IVZ:
        pmy_hydro->pmy_block->pcoord->PrimToLocal3(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
    }

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_hydro->peos->GetGamma();

  // Go through each non-pole interface
  #pragma simd
  for (int i = il+(pole_top?1:0); i <= iu-(pole_bottom?1:0); ++i)
  {
    // Extract left primitives
    const Real &rho_l = prim_l(IDN,i);
    const Real &pgas_l = prim_l(IEN,i);
    Real u_l[4];
    if (GENERAL_RELATIVITY)
    {
      u_l[1] = prim_l(ivx,i);
      u_l[2] = prim_l(ivy,i);
      u_l[3] = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 + SQR(u_l[1]) + SQR(u_l[2]) + SQR(u_l[3]));
    }
    else  // SR
    {
      const Real &vx_l = prim_l(ivx,i);
      const Real &vy_l = prim_l(ivy,i);
      const Real &vz_l = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 / (1.0 - SQR(vx_l) - SQR(vy_l) - SQR(vz_l)));
      u_l[1] = u_l[0] * vx_l;
      u_l[2] = u_l[0] * vy_l;
      u_l[3] = u_l[0] * vz_l;
    }

    // Extract right primitives
    const Real &rho_r = prim_r(IDN,i);
    const Real &pgas_r = prim_r(IEN,i);
    Real u_r[4];
    if (GENERAL_RELATIVITY)
    {
      u_r[1] = prim_r(ivx,i);
      u_r[2] = prim_r(ivy,i);
      u_r[3] = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 + SQR(u_r[1]) + SQR(u_r[2]) + SQR(u_r[3]));
    }
    else  // special relativity
    {
      const Real &vx_r = prim_r(ivx,i);
      const Real &vy_r = prim_r(ivy,i);
      const Real &vz_r = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 / (1.0 - SQR(vx_r) - SQR(vy_r) - SQR(vz_r)));
      u_r[1] = u_r[0] * vx_r;
      u_r[2] = u_r[0] * vy_r;
      u_r[3] = u_r[0] * vz_r;
    }

    // Calculate wavespeeds in left state (MB 23)
    Real lambda_p_l, lambda_m_l;
    Real wgas_l = rho_l + gamma_adi/(gamma_adi-1.0) * pgas_l;
    pmy_hydro->peos->SoundSpeedsSR(wgas_l, pgas_l, u_l[1]/u_l[0], SQR(u_l[0]),
        &lambda_p_l, &lambda_m_l);

    // Calculate wavespeeds in right state (MB 23)
    Real lambda_p_r, lambda_m_r;
    Real wgas_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r;
    pmy_hydro->peos->SoundSpeedsSR(wgas_r, pgas_r, u_r[1]/u_r[0], SQR(u_r[0]),
        &lambda_p_r, &lambda_m_r);

    // Calculate extremal wavespeed
    Real lambda_l = std::min(lambda_m_l, lambda_m_r);
    Real lambda_r = std::max(lambda_p_l, lambda_p_r);
    Real lambda = std::max(lambda_r, -lambda_l);

    // Calculate conserved quantities in L region (MB 3)
    Real cons_l[NWAVE];
    cons_l[IDN] = rho_l * u_l[0];
    cons_l[IEN] = wgas_l * u_l[0] * u_l[0] - pgas_l;
    cons_l[ivx] = wgas_l * u_l[1] * u_l[0];
    cons_l[ivy] = wgas_l * u_l[2] * u_l[0];
    cons_l[ivz] = wgas_l * u_l[3] * u_l[0];

    // Calculate fluxes in L region (MB 2,3)
    Real flux_l[NWAVE];
    flux_l[IDN] = rho_l * u_l[1];
    flux_l[IEN] = wgas_l * u_l[0] * u_l[1];
    flux_l[ivx] = wgas_l * u_l[1] * u_l[1] + pgas_l;
    flux_l[ivy] = wgas_l * u_l[2] * u_l[1];
    flux_l[ivz] = wgas_l * u_l[3] * u_l[1];

    // Calculate conserved quantities in R region (MB 3)
    Real cons_r[NWAVE];
    cons_r[IDN] = rho_r * u_r[0];
    cons_r[IEN] = wgas_r * u_r[0] * u_r[0] - pgas_r;
    cons_r[ivx] = wgas_r * u_r[1] * u_r[0];
    cons_r[ivy] = wgas_r * u_r[2] * u_r[0];
    cons_r[ivz] = wgas_r * u_r[3] * u_r[0];

    // Calculate fluxes in R region (MB 2,3)
    Real flux_r[NWAVE];
    flux_r[IDN] = rho_r * u_r[1];
    flux_r[IEN] = wgas_r * u_r[0] * u_r[1];
    flux_r[ivx] = wgas_r * u_r[1] * u_r[1] + pgas_r;
    flux_r[ivy] = wgas_r * u_r[2] * u_r[1];
    flux_r[ivz] = wgas_r * u_r[3] * u_r[1];

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n)
      flux(n,i) = 0.5 * (flux_l[n] + flux_r[n] - lambda * (cons_r[n] - cons_l[n]));

    // Set conserved quantities in GR
    if (GENERAL_RELATIVITY)
      for (int n = 0; n < NWAVE; ++n)
        cons_(n,i) = 0.5 * (cons_r[n] + cons_l[n] + (flux_l[n] - flux_r[n]) / lambda);
  }

  // Transform non-pole fluxes to global coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal1(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
      case IVY:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal2(k, j, il+(pole_top?1:0),
            iu-(pole_bottom?1:0), cons_, bb_normal_, flux);
        break;
      case IVZ:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal3(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
    }

  // Calculate pole fluxes if necessary
  if (GENERAL_RELATIVITY and (pole_top or pole_bottom))
    pmy_hydro->pmy_block->pcoord->Face2Metric(k, j, il, iu, g_, gi_);
  if (GENERAL_RELATIVITY and pole_top)
    NonTransformingLLF(k, j, il, prim_l, prim_r, g_, gi_, pmy_hydro, flux);
  if (GENERAL_RELATIVITY and pole_bottom)
    NonTransformingLLF(k, j, iu, prim_l, prim_r, g_, gi_, pmy_hydro, flux);
  return;
}

//--------------------------------------------------------------------------------------

// Fallback Riemann solver for poles
// Inputs:
//   k,j,i: x3-, x2-, and x1-indices
//   prim_l, prim_r: left and right primitive states
//   g,gi: 1D arrays of metric covariant and contravariant coefficients
//   pmy_hydro: pointer to Hydro
// Outputs:
//   flux: fluxes across interface
// Notes:
//   implements LLF algorithm similar to that of fluxcalc() in step_ch.c in Harm
//   adapted from RiemannSolver() in llf_rel_no_transform.cpp
static void NonTransformingLLF(const int k, const int j, const int i,
    const AthenaArray<Real> &prim_l, const AthenaArray<Real> &prim_r,
    const AthenaArray<Real> &g, const AthenaArray<Real> &gi, Hydro *pmy_hydro,
    AthenaArray<Real> &flux)
{
  // Fix variable direction
  int ivx = IVY;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_hydro->peos->GetGamma();

  // Extract metric
  const Real
      &g_00 = g(I00,i), &g_01 = g(I01,i), &g_02 = g(I02,i), &g_03 = g(I03,i),
      &g_10 = g(I01,i), &g_11 = g(I11,i), &g_12 = g(I12,i), &g_13 = g(I13,i),
      &g_20 = g(I02,i), &g_21 = g(I12,i), &g_22 = g(I22,i), &g_23 = g(I23,i),
      &g_30 = g(I03,i), &g_31 = g(I13,i), &g_32 = g(I23,i), &g_33 = g(I33,i);
  const Real
      &g00 = gi(I00,i), &g01 = gi(I01,i), &g02 = gi(I02,i), &g03 = gi(I03,i),
      &g10 = gi(I01,i), &g11 = gi(I11,i), &g12 = gi(I12,i), &g13 = gi(I13,i),
      &g20 = gi(I02,i), &g21 = gi(I12,i), &g22 = gi(I22,i), &g23 = gi(I23,i),
      &g30 = gi(I03,i), &g31 = gi(I13,i), &g32 = gi(I23,i), &g33 = gi(I33,i);
  Real alpha = std::sqrt(-1.0/g00);

  // Extract left primitives
  const Real &rho_l = prim_l(IDN,i);
  const Real &pgas_l = prim_l(IEN,i);
  const Real &uu1_l = prim_l(IVX,i);
  const Real &uu2_l = prim_l(IVY,i);
  const Real &uu3_l = prim_l(IVZ,i);

  // Extract right primitives
  const Real &rho_r = prim_r(IDN,i);
  const Real &pgas_r = prim_r(IEN,i);
  const Real &uu1_r = prim_r(IVX,i);
  const Real &uu2_r = prim_r(IVY,i);
  const Real &uu3_r = prim_r(IVZ,i);

  // Calculate 4-velocity in left state
  Real ucon_l[4], ucov_l[4];
  Real tmp = g_11*SQR(uu1_l) + 2.0*g_12*uu1_l*uu2_l + 2.0*g_13*uu1_l*uu3_l
           + g_22*SQR(uu2_l) + 2.0*g_23*uu2_l*uu3_l
           + g_33*SQR(uu3_l);
  Real gamma_l = std::sqrt(1.0 + tmp);
  ucon_l[0] = gamma_l / alpha;
  ucon_l[1] = uu1_l - alpha * gamma_l * g01;
  ucon_l[2] = uu2_l - alpha * gamma_l * g02;
  ucon_l[3] = uu3_l - alpha * gamma_l * g03;
  ucov_l[0] = g_00*ucon_l[0] + g_01*ucon_l[1] + g_02*ucon_l[2] + g_03*ucon_l[3];
  ucov_l[1] = g_10*ucon_l[0] + g_11*ucon_l[1] + g_12*ucon_l[2] + g_13*ucon_l[3];
  ucov_l[2] = g_20*ucon_l[0] + g_21*ucon_l[1] + g_22*ucon_l[2] + g_23*ucon_l[3];
  ucov_l[3] = g_30*ucon_l[0] + g_31*ucon_l[1] + g_32*ucon_l[2] + g_33*ucon_l[3];

  // Calculate 4-velocity in right state
  Real ucon_r[4], ucov_r[4];
  tmp = g_11*SQR(uu1_r) + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
      + g_22*SQR(uu2_r) + 2.0*g_23*uu2_r*uu3_r
      + g_33*SQR(uu3_r);
  Real gamma_r = std::sqrt(1.0 + tmp);
  ucon_r[0] = gamma_r / alpha;
  ucon_r[1] = uu1_r - alpha * gamma_r * g01;
  ucon_r[2] = uu2_r - alpha * gamma_r * g02;
  ucon_r[3] = uu3_r - alpha * gamma_r * g03;
  ucov_r[0] = g_00*ucon_r[0] + g_01*ucon_r[1] + g_02*ucon_r[2] + g_03*ucon_r[3];
  ucov_r[1] = g_10*ucon_r[0] + g_11*ucon_r[1] + g_12*ucon_r[2] + g_13*ucon_r[3];
  ucov_r[2] = g_20*ucon_r[0] + g_21*ucon_r[1] + g_22*ucon_r[2] + g_23*ucon_r[3];
  ucov_r[3] = g_30*ucon_r[0] + g_31*ucon_r[1] + g_32*ucon_r[2] + g_33*ucon_r[3];

  // Calculate wavespeeds in left state
  Real lambda_p_l, lambda_m_l;
  Real wgas_l = rho_l + gamma_adi/(gamma_adi-1.0) * pgas_l;
  pmy_hydro->peos->SoundSpeedsGR(wgas_l, pgas_l, ucon_l[0], ucon_l[ivx], g00, g02,
      g22, &lambda_p_l, &lambda_m_l);

  // Calculate wavespeeds in right state
  Real lambda_p_r, lambda_m_r;
  Real wgas_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r;
  pmy_hydro->peos->SoundSpeedsGR(wgas_r, pgas_r, ucon_r[0], ucon_r[ivx], g00, g02,
      g22, &lambda_p_r, &lambda_m_r);

  // Calculate extremal wavespeed
  Real lambda_l = std::min(lambda_m_l, lambda_m_r);
  Real lambda_r = std::max(lambda_p_l, lambda_p_r);
  Real lambda = std::max(lambda_r, -lambda_l);

  // Calculate conserved quantities in L region (rho u^0 and T^0_\mu)
  Real cons_l[NWAVE];
  cons_l[IDN] = rho_l * ucon_l[0];
  cons_l[IEN] = wgas_l * ucon_l[0] * ucov_l[0] + pgas_l;
  cons_l[IVX] = wgas_l * ucon_l[0] * ucov_l[1];
  cons_l[IVY] = wgas_l * ucon_l[0] * ucov_l[2];
  cons_l[IVZ] = wgas_l * ucon_l[0] * ucov_l[3];

  // Calculate fluxes in L region (rho u^i and T^i_\mu, where i = ivx)
  Real flux_l[NWAVE];
  flux_l[IDN] = rho_l * ucon_l[ivx];
  flux_l[IEN] = wgas_l * ucon_l[ivx] * ucov_l[0];
  flux_l[IVX] = wgas_l * ucon_l[ivx] * ucov_l[1];
  flux_l[IVY] = wgas_l * ucon_l[ivx] * ucov_l[2];
  flux_l[IVZ] = wgas_l * ucon_l[ivx] * ucov_l[3];
  flux_l[ivx] += pgas_l;

  // Calculate conserved quantities in R region (rho u^0 and T^0_\mu)
  Real cons_r[NWAVE];
  cons_r[IDN] = rho_r * ucon_r[0];
  cons_r[IEN] = wgas_r * ucon_r[0] * ucov_r[0] + pgas_r;
  cons_r[IVX] = wgas_r * ucon_r[0] * ucov_r[1];
  cons_r[IVY] = wgas_r * ucon_r[0] * ucov_r[2];
  cons_r[IVZ] = wgas_r * ucon_r[0] * ucov_r[3];

  // Calculate fluxes in R region (rho u^i and T^i_\mu, where i = ivx)
  Real flux_r[NWAVE];
  flux_r[IDN] = rho_r * ucon_r[ivx];
  flux_r[IEN] = wgas_r * ucon_r[ivx] * ucov_r[0];
  flux_r[IVX] = wgas_r * ucon_r[ivx] * ucov_r[1];
  flux_r[IVY] = wgas_r * ucon_r[ivx] * ucov_r[2];
  flux_r[IVZ] = wgas_r * ucon_r[ivx] * ucov_r[3];
  flux_r[ivx] += pgas_r;

  // Set fluxes
  for (int n = 0; n < NWAVE; ++n)
    flux(n,i) = 0.5 * (flux_l[n] + flux_r[n] - lambda * (cons_r[n] - cons_l[n]));
  return;
}