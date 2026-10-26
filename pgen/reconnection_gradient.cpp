//========================================================================================
// Athena++ astrophysical MHD code
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconnection_gradient.cpp
//! \brief 2D resistive-MHD reconnection with a prescribed resistivity profile.
//!
//! Reproduces the setup of Lin, Liu & Li, "Fast Magnetic Reconnection induced by
//! Resistivity Gradients in 2D MHD" (arXiv:2109.07526):
//!   - force-free current sheet along y centred at x = xc:
//!         B = B0 tanh((x-xc)/lambda) yhat + B0 sech((x-xc)/lambda) zhat
//!   - uniform rho, uniform p = beta*B0^2/2  (so V_A = B0/sqrt(rho) = 1)
//!   - Ohmic resistivity eta(x,y) (see iprofile below), imposed through the
//!     Athena++ user field-diffusivity hook
//!   - zero-gradient ("outflow") boundaries on all sides (set in the input file)
//!
//! Resistivity profiles (<problem>/iprofile):
//!   0: hyperbolic tangent in y (Eq. 6):  eta = 0.5*eta1*(1 - tanh((y-yc)/l_eta)) + eta2
//!   1: exponential in radius (Eq. 12):   eta = eta1*exp(-r/l_eta) + eta2,
//!                                        r = sqrt((x-xc)^2 + (y-yc)^2)
//!   2: uniform:                          eta = eta1 + eta2
//!
//! NOTE: Athena++ only evaluates field diffusion when <problem>/eta_ohm != 0, so the
//! input file must set eta_ohm to any non-zero value (it is otherwise ignored here).
//!
//! User output variables (use `variable = uov` in an <outputN> block):
//!   Jz   out-of-plane current density
//!   eta  resistivity
//!   Ez   out-of-plane electric field  (= -(v x B)_z + eta*Jz); R = Ez/(B0*VA) at x-line
//!   Az   z-component of the vector potential on cell centres is NOT provided; use the
//!        analysis scripts to integrate it from the field if flux contours are needed.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
int iprofile;
Real eta1, eta2, leta, xc, yc;

Real EtaProfile(Real x, Real y) {
  switch (iprofile) {
    case 0:
      return 0.5*eta1*(1.0 - std::tanh((y - yc)/leta)) + eta2;
    case 1: {
      Real r = std::sqrt(SQR(x - xc) + SQR(y - yc));
      return eta1*std::exp(-r/leta) + eta2;
    }
    default:
      return eta1 + eta2;
  }
}

void ResistivityProfile(FieldDiffusion *pfdif, MeshBlock *pmb,
                        const AthenaArray<Real> &w, const AthenaArray<Real> &bmag,
                        const int is, const int ie, const int js, const int je,
                        const int ks, const int ke) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        pfdif->etaB(FieldDiffusion::DiffProcess::ohmic, k, j, i) =
            EtaProfile(pmb->pcoord->x1v(i), pmb->pcoord->x2v(j));
      }
    }
  }
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief read the resistivity parameters and enroll the resistivity function

void Mesh::InitUserMeshData(ParameterInput *pin) {
  iprofile = pin->GetOrAddInteger("problem", "iprofile", 0);
  eta1     = pin->GetReal("problem", "eta1");
  eta2     = pin->GetReal("problem", "eta2");
  leta     = pin->GetOrAddReal("problem", "l_eta", 0.05);
  xc       = pin->GetOrAddReal("problem", "xc", 1.0);
  yc       = pin->GetOrAddReal("problem", "yc", 1.0);

  if (pin->GetOrAddReal("problem", "eta_ohm", 0.0) == 0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in reconnection_gradient.cpp InitUserMeshData" << std::endl
        << "<problem>/eta_ohm must be non-zero to switch on Ohmic diffusion "
        << "(its value is not used)" << std::endl;
    ATHENA_ERROR(msg);
  }
  EnrollFieldDiffusivity(ResistivityProfile);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief force-free current sheet with an optional localized flux perturbation

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const Real B0     = pin->GetOrAddReal("problem", "B0", 1.0);
  const Real lambda = pin->GetOrAddReal("problem", "lambda", 0.04);
  const Real beta   = pin->GetOrAddReal("problem", "beta", 0.1);
  const Real d0     = pin->GetOrAddReal("problem", "d0", 1.0);
  // Flux perturbation dAz = pert_amp*B0*lambda*exp(-(r/pert_width)^2), r from (xc, yc).
  // The paper does not specify a perturbation; the default is small and mimics a seed.
  const Real pert_amp   = pin->GetOrAddReal("problem", "pert_amp", 0.0);
  const Real pert_width = pin->GetOrAddReal("problem", "pert_width", 0.2);
  const Real pgas = 0.5*beta*B0*B0;

  auto Az = [&](Real x, Real y) {
    return pert_amp*B0*lambda*std::exp(-(SQR(x - xc) + SQR(y - yc))/SQR(pert_width));
  };

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = d0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
      }
    }
  }

  // Divergence-free field from the vector potential (perturbation) plus the equilibrium
  // By(x), Bz(x). Bx = dAz/dy, By = B0 tanh - dAz/dx (face-centred finite differences).
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          Real xf = pcoord->x1f(i);
          pfield->b.x1f(k,j,i) =
              (Az(xf, pcoord->x2f(j+1)) - Az(xf, pcoord->x2f(j)))/pcoord->dx2f(j);
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real yf = pcoord->x2f(j);
          pfield->b.x2f(k,j,i) = B0*std::tanh((pcoord->x1v(i) - xc)/lambda)
              - (Az(pcoord->x1f(i+1), yf) - Az(pcoord->x1f(i), yf))/pcoord->dx1f(i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = B0/std::cosh((pcoord->x1v(i) - xc)/lambda);
        }
      }
    }
    pfield->CalculateCellCenteredField(pfield->b, pfield->bcc, pcoord,
                                       is, ie, js, je, ks, ke);
  }

  // Total energy: thermal + magnetic (v = 0).
  if (NON_BAROTROPIC_EOS) {
    const Real igm1 = 1.0/(peos->GetGamma() - 1.0);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          phydro->u(IEN,k,j,i) = pgas*igm1;
          if (MAGNETIC_FIELDS_ENABLED) {
            phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->bcc(IB1,k,j,i))
                                       + SQR(pfield->bcc(IB2,k,j,i))
                                       + SQR(pfield->bcc(IB3,k,j,i)));
          }
        }
      }
    }
  }

  AllocateUserOutputVariables(3);
  SetUserOutputVariableName(0, "Jz");
  SetUserOutputVariableName(1, "eta");
  SetUserOutputVariableName(2, "Ez");
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//! \brief fill Jz, eta and Ez (cell-centred, second-order differences)

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        // Jz = dBy/dx - dBx/dy with one-sided stencils at the block edge
        const int im = std::max(i-1, is-NGHOST), ip = std::min(i+1, ie+NGHOST);
        const int jm = std::max(j-1, js-NGHOST), jp = std::min(j+1, je+NGHOST);
        Real dbydx = (pfield->bcc(IB2,k,j,ip) - pfield->bcc(IB2,k,j,im))
                   / (pcoord->x1v(ip) - pcoord->x1v(im));
        Real dbxdy = (pfield->bcc(IB1,k,jp,i) - pfield->bcc(IB1,k,jm,i))
                   / (pcoord->x2v(jp) - pcoord->x2v(jm));
        Real jz  = dbydx - dbxdy;
        Real eta = EtaProfile(pcoord->x1v(i), pcoord->x2v(j));
        Real vx = phydro->w(IVX,k,j,i), vy = phydro->w(IVY,k,j,i);
        Real bx = pfield->bcc(IB1,k,j,i), by = pfield->bcc(IB2,k,j,i);
        user_out_var(0,k,j,i) = jz;
        user_out_var(1,k,j,i) = eta;
        user_out_var(2,k,j,i) = -(vx*by - vy*bx) + eta*jz;
      }
    }
  }
}
