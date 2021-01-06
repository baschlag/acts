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
  double tanPhiV = std::tan(phiV);

  // theta and functions
  double theta = paramsAtPCA(BoundIndices::eBoundTheta);
  const double sinTh = std::sin(theta);
  const double cosTh = std::cos(theta);
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

  const double Tx = sinTh*cosPhiV;
  const double Ty = sinTh * sinPhiV;
  const double Tz = std::cos(theta);

  const double qpBz = qOvP * Bz;
  const double sqrtTzTz = std::sqrt(1-Tz*Tz);

  double deltaX = positionAtPCA(0) - linPointPos.x();
  double deltaY = positionAtPCA(1) - linPointPos.y();

  double newX = deltaX + Ty/(qpBz);
  double newY = deltaY - Tx/(qpBz);
  double newS2 = newX * newX + newY * newY;
  double newS = std::sqrt(newS2);

  // Eq. 5.34 in Ref(1) (see .hpp)
  double X = positionAtPCA(0) - linPointPos.x() + rho * sinPhiV;
  double Y = positionAtPCA(1) - linPointPos.y() - rho * cosPhiV;
  const double S2 = (X * X + Y * Y);
  const double S = std::sqrt(S2);

  int sgnX = (newX < 0.) ? -1 : 1;
  int sgnY = (newY < 0.) ? -1 : 1;

  // phi at the point of closest approach (perigee)
  double phiP;
  if (std::abs(newX) > std::abs(newY)) {
    phiP = sgnH * sgnX * std::acos(-sgnH * newY / newS);
  } else {
    phiP = std::asin(sgnH * newX / newS);
    if ((sgnH * sgnY) > 0) {
      phiP = sgnH * sgnX * M_PI - phiP;
    }
  } 
  
  double d0 = sqrtTzTz/(qOvP*Bz) - sgnH * newS;
  double z0 = positionAtPCA[eZ] + linPointPos.z() + Tz/(qOvP*Bz) * (phiV - phiP);
  /// F(V, p_i) at PCA in Billoir paper
  /// (see FullBilloirVertexFitter.hpp for paper reference,
  /// Page 140, Eq. (2) )
  BoundVector predParamsAtPCA;

  // Eq. 5.33 in Ref(1) (see .hpp)
  predParamsAtPCA[0] = d0;
  predParamsAtPCA[1] = z0;
  predParamsAtPCA[2] = phiP;
  predParamsAtPCA[3] = theta;
  predParamsAtPCA[4] = qOvP;
  predParamsAtPCA[5] = 0.;

  double dPhi = phiP - phiV;
  double hOvS = sgnH/newS;  

  FreeToBoundMatrix freeToBoundJacobian{FreeToBoundMatrix::Zero()};

  double Tt = std::sqrt(Tx*Tx + Ty*Ty);

  
  double dRhodTxNew = cosTh /(qpBz * Tt) * Tx*Tz/(Tz*Tz + Tt*Tt);
  double dRhodTyNew = cosTh /(qpBz * Tt) * Ty*Tz/(Tz*Tz + Tt*Tt);
  double dRhodTzNew = -cosTh /qpBz * Tt/(Tz*Tz + Tt*Tt);
  double dRhodQPNew = -sinTh/(qpBz*qOvP);

  double dRhodTx = Tz /(qpBz * Tt) * Tx*Tz/(Tz*Tz + Tt*Tt);
  double dRhodTy = Tz /(qpBz * Tt) * Ty*Tz/(Tz*Tz + Tt*Tt);
  double dRhodTz = -Tz /qpBz * Tt/(Tz*Tz + Tt*Tt);
  double dRhodQP = -sinTh/(qpBz*qOvP);

  //std::cout << "old/new dRhodTx: " << dRhodTx << "," << dRhodTxNew << std::endl;
  //astd::cout << "old/new dRhodTy: " << dRhodTy << "," << dRhodTyNew << std::endl;
  //std::cout << "old/new dRhodTz: " << dRhodTz << "," << dRhodTzNew << std::endl;

  double dPhidTx = -Ty/(Tx*Tx + Ty*Ty);
  double dPhidTy = Tx/(Tx*Tx + Ty*Ty);

  double myR = newY * sinPhiV + newX * cosPhiV;
  double myQ = newX * sinPhiV - newY * cosPhiV;

  double Tt3 = Tt*Tt*Tt;

  // d0 derivatives
  freeToBoundJacobian(0,0) = -sgnH * newX / S; // checked
  freeToBoundJacobian(0,1) = -sgnH * newY / S; // checked
  freeToBoundJacobian(0,2) = 0.; // checked
  freeToBoundJacobian(0,3) = 0.; // checked
  freeToBoundJacobian(0,4) = dRhodTx -hOvS *(dRhodTx*myQ + rho * dPhidTx*myR);
  freeToBoundJacobian(0,5) = dRhodTy -hOvS *(dRhodTy*myQ + rho * dPhidTy*myR);
  freeToBoundJacobian(0,6) = dRhodTz * (1- hOvS*myQ);
  freeToBoundJacobian(0,7) = dRhodQP * (1- hOvS*myQ);

  // phi derivatives
  freeToBoundJacobian(2,0) = - newY / newS2; // checked
  freeToBoundJacobian(2,1) = newX / newS2; // checked
  freeToBoundJacobian(2,2) = 0.; // checked
  freeToBoundJacobian(2,3) = 0.; // checked
  freeToBoundJacobian(2,4) = 1/newS2 * (-dRhodTx* myR + rho * dPhidTx * myQ);
  freeToBoundJacobian(2,5) = 1/newS2 * (-dRhodTy* myR + rho * dPhidTy * myQ);
  freeToBoundJacobian(2,6) = -1/newS2 * dRhodTz * myR;
  freeToBoundJacobian(2,7) = -1/newS2 * dRhodQP * myR;

  // z0 derivatives
  freeToBoundJacobian(1,0) = Tz * newY / (qpBz * newS2); // checked
  freeToBoundJacobian(1,1) = -Tz * newX / (qpBz * newS2); // checked
  freeToBoundJacobian(1,2) = 1.; // checked
  freeToBoundJacobian(1,3) = 0.; // checked
  freeToBoundJacobian(1,4) = -(dRhodTx * Tz/Tt - rho * Tz * Tx / Tt3) * dPhi + rho * Tz/Tt *(dPhidTx - freeToBoundJacobian(2,4));
  freeToBoundJacobian(1,5) = -(dRhodTy * Tz/Tt - rho * Tz * Ty / Tt3) * dPhi + rho * Tz/Tt *(dPhidTy - freeToBoundJacobian(2,5));
  freeToBoundJacobian(1,6) = -(dRhodTz * Tz/Tt + rho/Tt) * dPhi - rho * Tz/Tt * freeToBoundJacobian(2,6);
  freeToBoundJacobian(1,7) = -dRhodQP * Tz/Tt * dPhi -rho * Tz/Tt * freeToBoundJacobian(2,7);

  // theta derivatives
  freeToBoundJacobian(3,4) = Tx*Tz/(Tz*Tz + Tt * Tt) * 1/Tt;
  freeToBoundJacobian(3,5) = Ty*Tz/(Tz*Tz + Tt * Tt) * 1/Tt;
  freeToBoundJacobian(3,6) = -Tt/(Tz*Tz + Tt * Tt);
  // q/p
  freeToBoundJacobian(4,7) = 1.;
  // time
  freeToBoundJacobian(5,3) = 1.;

  // Calculate boundToFree matrix for covariance transformation
  BoundToFreeMatrix boundToFreeJacobian{BoundToFreeMatrix::Zero()};
  double rhoMinusD0 = rho - d0;
  double dRhodTheta = cosTh/qpBz;
  double sinPhiP = std::sin(phiP);
  double cosPhiP = std::cos(phiP);
  double deltaSinPhi = sinPhiP - sinPhiV;
  double deltaCosPhi = cosPhiP - cosPhiV;

  // x-component derivatives
  boundToFreeJacobian(0,0) = -sinPhiP;
  boundToFreeJacobian(0,2) = cosPhiP * rhoMinusD0;
  boundToFreeJacobian(0,3) = dRhodTheta * deltaSinPhi;
  boundToFreeJacobian(0,4) = dRhodQP * deltaSinPhi;

  // y-component derivatives
  boundToFreeJacobian(1,0) = cosPhiP;
  boundToFreeJacobian(1,2) = sinPhiP * rhoMinusD0;
  boundToFreeJacobian(1,3) = -dRhodTheta * deltaCosPhi;
  boundToFreeJacobian(1,4) = -dRhodQP * deltaCosPhi;

  // z-component derivatives
  boundToFreeJacobian(2,1) = 1.;
  boundToFreeJacobian(2,2) = dRhodTheta;
  boundToFreeJacobian(2,3) = -sinTh/qpBz * dPhi;
  boundToFreeJacobian(2,4) = -cosTh/(qpBz*qOvP) * dPhi;

  // t-component derivatives
  boundToFreeJacobian(3,5) = 1.;

  double xi = deltaY - deltaX*tanPhiV;
  double dPhidPhiP = (deltaX/cosPhiV + sinPhiV) * xi / std::sqrt(rho*rho - cosPhiV*cosPhiV * xi * xi);
  
  // Tx-component derivatives
  boundToFreeJacobian(4,2) = -sinPhiP*sinTh * dPhidPhiP;
  boundToFreeJacobian(4,3) = cosPhiP*cosTh;

  // Ty-component derivatives
  boundToFreeJacobian(5,2) = sinTh*cosPhiP * dPhidPhiP;
  boundToFreeJacobian(5,3) = cosTh*sinPhiP;
  
  // Tz-component derivatives
  boundToFreeJacobian(6,3) = -sinTh;

  // q/p-component derivatives
  boundToFreeJacobian(7,4) = 1.;

  // Calculate free covariance matrix
  FreeSymMatrix freeCovarianceAtPCA;
  freeCovarianceAtPCA = boundToFreeJacobian * parCovarianceAtPCA * boundToFreeJacobian.transpose();

  // Calculate free weight matrix
  FreeSymMatrix freeWeightMatrix = freeCovarianceAtPCA.inverse();

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

  FreeVector freeVector;
  freeVector << positionAtPCA[ePos0], positionAtPCA[ePos1], positionAtPCA[ePos2], 0.,// TODO positionAtPCA[eTime],
                Tx, Ty, Tz, qOvP;
  BoundVector newConstTerm = predParamsAtPCA - freeToBoundJacobian * freeVector;

  // The parameter weight
  ActsSymMatrixD<5> parWeight =
      (parCovarianceAtPCA.block<5, 5>(0, 0)).inverse();

  BoundSymMatrix weightAtPCA{BoundSymMatrix::Identity()};
  weightAtPCA.block<5, 5>(0, 0) = parWeight;

  return LinearizedTrack(paramsAtPCA, freeVector, parCovarianceAtPCA, freeCovarianceAtPCA, weightAtPCA,
                        freeWeightMatrix, linPoint, freeToBoundJacobian,
                         positionJacobian, momentumJacobian, positionAtPCA,
                         momentumAtPCA, constTerm, newConstTerm);
}
