// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"

template <typename propagator_t, typename propagator_options_t>
Acts::Result<Acts::LinearizedTrack> Acts::
    HelicalTrackLinearizer<propagator_t, propagator_options_t>::linearizeTrack(
        const BoundTrackParameters& params, const Vector4D& linPoint,
        const Acts::GeometryContext& gctx,
        const Acts::MagneticFieldContext& mctx, State& state) const {
  Vector3D linPointPos = VectorHelpers::position(linPoint);

  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(linPointPos);

  // Create propagator options
  auto logger = getDefaultLogger("HelTrkLinProp", Logging::INFO);
  propagator_options_t pOptions(gctx, mctx, LoggerWrapper{*logger});
  pOptions.direction = backward;

  const BoundTrackParameters* endParams = nullptr;
  // Do the propagation to linPointPos
  auto result = m_cfg.propagator->propagate(params, *perigeeSurface, pOptions);
  if (result.ok()) {
    endParams = (*result).endParameters.get();

  } else {
    return result.error();
  }

  BoundVector paramsAtPCA = endParams->parameters();
  Vector4D positionAtPCA = Vector4D::Zero();
  {
    auto pos = endParams->position(gctx);
    positionAtPCA[ePos0] = pos[ePos0];
    positionAtPCA[ePos1] = pos[ePos1];
    positionAtPCA[ePos2] = pos[ePos2];
    positionAtPCA[eTime] = endParams->time();
  }
  BoundSymMatrix parCovarianceAtPCA = *(endParams->covariance());

  if (endParams->covariance()->determinant() == 0) {
    // Use the original parameters
    paramsAtPCA = params.parameters();
    auto pos = endParams->position(gctx);
    positionAtPCA[ePos0] = pos[ePos0];
    positionAtPCA[ePos1] = pos[ePos1];
    positionAtPCA[ePos2] = pos[ePos2];
    parCovarianceAtPCA = *(params.covariance());
  }

  // phiV and functions
  double phiV = paramsAtPCA(BoundIndices::eBoundPhi);
  double sinPhiV = std::sin(phiV);
  double cosPhiV = std::cos(phiV);

  // theta and functions
  double theta = paramsAtPCA(BoundIndices::eBoundTheta);
  const double sinTh = std::sin(theta);
  const double tanTh = std::tan(theta);

  // q over p
  double qOvP = paramsAtPCA(BoundIndices::eBoundQOverP);
  double sgnH = (qOvP < 0.) ? -1 : 1;

  Vector3D momentumAtPCA(phiV, theta, qOvP);

  // get B-field z-component at current position
  double Bz = m_cfg.bField.getField(VectorHelpers::position(positionAtPCA),
                                    state.fieldCache)[eZ];
  double rho;
  // Curvature is infinite w/o b field
  if (Bz == 0. || std::abs(qOvP) < m_cfg.minQoP) {
    rho = m_cfg.maxRho;
  } else {
    rho = sinTh * (1. / qOvP) / Bz;
  }

  const double lambda = M_PI/2. - theta;
  const double cosLambda = std::cos(lambda);
  const double Tx = std::cos(lambda)*cosPhiV;
  const double Ty = cosLambda * sinPhiV;
  const double Tz = std::sin(lambda);

  const double qpBz = qOvP * Bz;
  const double sqrtTz = std::sqrt(1-Tz*Tz);

  double newX = positionAtPCA(0) - linPointPos.x() + Ty/(qpBz);
  double newY = positionAtPCA(1) - linPointPos.y() - Tx/(qpBz);
  double newS2 = newX * newX + newY * newY;
  double newS = std::sqrt(newS2);

  // Eq. 5.34 in Ref(1) (see .hpp)
  double X = positionAtPCA(0) - linPointPos.x() + rho * sinPhiV;
  double Y = positionAtPCA(1) - linPointPos.y() - rho * cosPhiV;
  const double S2 = (X * X + Y * Y);
  const double S = std::sqrt(S2);

  int sgnX = (newX < 0.) ? -1 : 1;
  int sgnY = (newY < 0.) ? -1 : 1;

  double phiAtPCA;
  if (std::abs(newX) > std::abs(newY)) {
    phiAtPCA = sgnH * sgnX * std::acos(-sgnH * newY / newS);
  } else {
    phiAtPCA = std::asin(sgnH * newX / newS);
    if ((sgnH * sgnY) > 0) {
      phiAtPCA = sgnH * sgnX * M_PI - phiAtPCA;
    }
  }

  double d0 = sqrtTz/(qOvP*Bz) - sgnH * newS;
  double z0 = positionAtPCA[eZ] + linPointPos.z() + Tz/(qOvP*Bz) * (phiV - phiAtPCA);
  /// F(V, p_i) at PCA in Billoir paper
  /// (see FullBilloirVertexFitter.hpp for paper reference,
  /// Page 140, Eq. (2) )
  BoundVector predParamsAtPCA;

  // Eq. 5.33 in Ref(1) (see .hpp)
  predParamsAtPCA[0] = d0;
  predParamsAtPCA[1] = z0;
  predParamsAtPCA[2] = phiAtPCA;
  predParamsAtPCA[3] = theta;
  predParamsAtPCA[4] = qOvP;
  predParamsAtPCA[5] = 0.;

  double dPhi = phiAtPCA - phiV;

  FreeToBoundMatrix freeToBoundJacobian{FreeToBoundMatrix::Zero()};

  freeToBoundJacobian(0,0) = -sgnH * newX / S;
  freeToBoundJacobian(0,1) = -sgnH * neWY / S;
  freeToBoundJacobian(0,2) = 0.;
  freeToBoundJacobian(0,3) = 0.;
  freeToBoundJacobian(0,4) = sgnH * newY / (qpBz * newS);
  freeToBoundJacobian(0,5) = - sgnH * newX / (qpBz * newS);;
  freeToBoundJacobian(0,6) = - Tz/(qpBz * sqrtTz);
  freeToBoundJacobian(0,7) = - 1/(qpBz * qOvP) * (sgnH*(Tx*newY - Ty*newX)/newS + sqrtTz);

  freeToBoundJacobian(1,0) = Tz * newX / (qpBz * newS2);
  freeToBoundJacobian(1,1) = -Tz * newY / (qpBz * newS2);
  freeToBoundJacobian(1,2) = 1.;
  freeToBoundJacobian(1,3) = 0.;
  freeToBoundJacobian(1,4) = Tz * newX / (qpBz * qpBz * newS2);
  freeToBoundJacobian(1,5) = Tz / (qpBz * qpBz * newS2) * (newS2 * qOvP / std::sqrt(1-Ty*Ty-Tz*Tz) + newY);
  freeToBoundJacobian(1,6) = 1./(qpBz) * ( Ty*Tz*Tz/((1-Tz*Tz)*(1-Tz*Tz)*std::sqrt(1-Ty*Ty-Tz*Tz)) + dPhi);
  freeToBoundJacobian(1,7) = -Tz / (qpBz * qpBz * newS2) * ( (Tx*newX + Ty*newY)/qOvP + Bz * newS2 * dPhi);

  freeToBoundJacobian(2,0) = - newY / newS2;
  freeToBoundJacobian(2,1) = newX / newS2;
  freeToBoundJacobian(2,2) = 0.;
  freeToBoundJacobian(2,3) = 0.;
  freeToBoundJacobian(2,4) = -newX / (qpBz * newS2);
  freeToBoundJacobian(2,5) = -newY / (qpBz * newS2)
  freeToBoundJacobian(2,6) = 0;
  freeToBoundJacobian(2,7) = (Tx*newX + Ty*newY)/(qpBz * qOvP * newS2);


  // TODO: add other derivatives...

// Fill position jacobian (D_k matrix), Eq. 5.36 in Ref(1)
  ActsMatrix<BoundScalar, eBoundSize, 4> positionJacobian;
  positionJacobian.setZero();
  // // First row
  positionJacobian(0, 0) = -sgnH * X / S;
  positionJacobian(0, 1) = -sgnH * Y / S;

  const double S2tanTh = S2 * tanTh;

  // Second row
  positionJacobian(1, 0) = rho * Y / S2tanTh;
  positionJacobian(1, 1) = -rho * X / S2tanTh;
  positionJacobian(1, 2) = 1.;

  // Third row
  positionJacobian(2, 0) = -Y / S2;
  positionJacobian(2, 1) = X / S2;

  // TODO: include timing in track linearization
  // Last row
  positionJacobian(5, 3) = 1;

  // Fill momentum jacobian (E_k matrix), Eq. 5.37 in Ref(1)
  ActsMatrixD<eBoundSize, 3> momentumJacobian;
  momentumJacobian.setZero();

  double R = X * cosPhiV + Y * sinPhiV;
  double Q = X * sinPhiV - Y * cosPhiV;

  // First row
  momentumJacobian(0, 0) = -sgnH * rho * R / S;

  double qOvSred = 1 - sgnH * Q / S;

  momentumJacobian(0, 1) = qOvSred * rho / tanTh;
  momentumJacobian(0, 2) = -qOvSred * rho / qOvP;

  const double rhoOverS2 = rho / S2;

  // Second row
  momentumJacobian(1, 0) = (1 - rhoOverS2 * Q) * rho / tanTh;
  momentumJacobian(1, 1) = (dPhi + rho * R / (S2tanTh * tanTh)) * rho;
  momentumJacobian(1, 2) = (dPhi - rhoOverS2 * R) * rho / (qOvP * tanTh);

  // Third row
  momentumJacobian(2, 0) = rhoOverS2 * Q;
  momentumJacobian(2, 1) = -rho * R / S2tanTh;
  momentumJacobian(2, 2) = rhoOverS2 * R / qOvP;

  // Last two rows:
  momentumJacobian(3, 1) = 1.;
  momentumJacobian(4, 2) = 1.;

  // const term F(V_0, p_0) in Talyor expansion
  BoundVector constTerm = predParamsAtPCA - positionJacobian * positionAtPCA -
                          momentumJacobian * momentumAtPCA;

  // The parameter weight
  ActsSymMatrixD<5> parWeight =
      (parCovarianceAtPCA.block<5, 5>(0, 0)).inverse();

  BoundSymMatrix weightAtPCA{BoundSymMatrix::Identity()};
  weightAtPCA.block<5, 5>(0, 0) = parWeight;

  return LinearizedTrack(paramsAtPCA, parCovarianceAtPCA, weightAtPCA, linPoint,
                         positionJacobian, momentumJacobian, positionAtPCA,
                         momentumAtPCA, constTerm);
}
