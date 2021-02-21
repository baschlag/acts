// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

template <typename propagator_t, typename propagator_options_t>
Acts::Result<Acts::LinearizedTrack> Acts::
    HelicalTrackLinearizer<propagator_t, propagator_options_t>::linearizeTrack(
        const BoundTrackParameters& params, const Vector4D& linPoint,
        const Acts::GeometryContext& gctx,
        const Acts::MagneticFieldContext& mctx, State& state) const {
  Vector3D linPointPos = VectorHelpers::position(linPoint);

  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(linPointPos);

      const std::shared_ptr<PerigeeSurface> perigeeSurface2 =
      Surface::makeShared<PerigeeSurface>(linPointPos + Vector3D(1e-3,1e-3,1e-3));

  // std::cout << "linPointPos: " << linPointPos.transpose() << std::endl;

  // Create propagator options
  auto logger = getDefaultLogger("HelTrkLinProp", Logging::INFO);
  propagator_options_t pOptions(gctx, mctx, LoggerWrapper{*logger});
  pOptions.direction = backward;
  // std::cout << "here 1 " << std::endl;
  // STEP 1: propagate to pergiee at linpoint (should be inside beam pipe)

  const BoundTrackParameters* endParams = nullptr;
  // Do the propagation to linPointPos
  auto result = m_cfg.propagator->propagate(params, *perigeeSurface, pOptions);
  if (result.ok()) {
    endParams = (*result).endParameters.get();
  } else {
    return result.error();
  }
  pOptions.direction = forward;
  // std::cout << "here 2 " << std::endl;
  // STEP 2: propagate beam pipe like surface

  double radius(10.0), halfZ(1e5);
  Translation3D translation{0., 0., 0.};
  auto pTransform = Transform3D(translation);
  auto cylinderSurface =
      Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);

  const BoundTrackParameters* endParamsBeampipe = nullptr;
  auto resultBeampipe = m_cfg.propagator->propagate(params, *perigeeSurface2, pOptions);
  //auto resultBeampipe = m_cfg.propagator->propagate(params, *cylinderSurface, pOptions);
  if (resultBeampipe.ok()) {
    endParamsBeampipe = (*resultBeampipe).endParameters.get();
  } else {
    return resultBeampipe.error();
  }
  // std::cout << "here 3 " << std::endl;
  // STEP 3: convert to free parameters at beam pipe surface

  using PropagatorOptions = PropagatorOptions<ActionList<>, AbortList<TestAborter>>;
  PropagatorOptions pOptionsAbort(gctx, mctx, LoggerWrapper{*logger});

  const FreeTrackParameters* endFreeParams = nullptr;
  auto freeresult = m_cfg.propagator->template propagate<FreeTrackParameters>(*endParamsBeampipe, pOptionsAbort);
  if (freeresult.ok()) {
    endFreeParams = (*freeresult).endParameters.get();
  } else {
    return freeresult.error();
  }

  // std::cout << "here 4 " << std::endl;
  // STEP 4: propagate back to perigee at linpoint, going from FreeParameters to BoundParameters --> retrieve jacobian
  pOptions.direction = backward;
  FreeToBoundMatrix freeToBoundJacobian;
  const BoundTrackParameters* backBoundParams = nullptr;
  // Do the propagation to linPointPos
  auto backBoundResult = m_cfg.propagator->propagate(*endFreeParams,*perigeeSurface, pOptions);
  if (backBoundResult.ok()) {
    backBoundParams = (*backBoundResult).endParameters.get();
    auto jacVariant = *((*backBoundResult).transportJacobian);
    freeToBoundJacobian = std::get<1>(jacVariant);
  } else {
    return result.error();
  }

  // std::cout << "here 5 " << std::endl;
  ActsMatrix<BoundScalar, eBoundSize, 4> newPositionJacobian;
  newPositionJacobian = freeToBoundJacobian.block<6,4>(0,0);

  // std::cout << "new positionJacobian: \n" << newPositionJacobian << std::endl;

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
  double th = paramsAtPCA(BoundIndices::eBoundTheta);
  const double sinTh = std::sin(th);
  const double cosTh = std::cos(th);
  const double tanTh = std::tan(th);

  // q over p
  double qOvP = paramsAtPCA(BoundIndices::eBoundQOverP);
  double sgnH = (qOvP < 0.) ? -1 : 1;

  Vector3D momentumAtPCA(phiV, th, qOvP);

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

  // Eq. 5.34 in Ref(1) (see .hpp)
  double X = positionAtPCA(0) - linPointPos.x() + rho * sinPhiV;
  double Y = positionAtPCA(1) - linPointPos.y() - rho * cosPhiV;
  const double S2 = (X * X + Y * Y);
  const double S = std::sqrt(S2);

  /// F(V, p_i) at PCA in Billoir paper
  /// (see FullBilloirVertexFitter.hpp for paper reference,
  /// Page 140, Eq. (2) )
  BoundVector predParamsAtPCA;

  int sgnX = (X < 0.) ? -1 : 1;
  int sgnY = (Y < 0.) ? -1 : 1;

  double phiAtPCA;
  if (std::abs(X) > std::abs(Y)) {
    phiAtPCA = sgnH * sgnX * std::acos(-sgnH * Y / S);
  } else {
    phiAtPCA = std::asin(sgnH * X / S);
    if ((sgnH * sgnY) > 0) {
      phiAtPCA = sgnH * sgnX * M_PI - phiAtPCA;
    }
  }

  // std::cout << "phi: " << phiAtPCA << ", " << phiV << std::endl;

  // Eq. 5.33 in Ref(1) (see .hpp)
  // predParamsAtPCA[0] = rho - sgnH * S;
  // predParamsAtPCA[1] =
  //     positionAtPCA[eZ] - linPointPos.z() + rho * (phiV - phiAtPCA) / tanTh;
  // predParamsAtPCA[2] = phiAtPCA;
  // predParamsAtPCA[3] = th;
  // predParamsAtPCA[4] = qOvP;
  // predParamsAtPCA[5] = 0.;

  // std::cout << "predParamsAtPCA: " <<  predParamsAtPCA << std::endl;
  // std::cout << "paramsAtPCA: " <<  paramsAtPCA << std::endl;

  // Fill position jacobian (D_k matrix), Eq. 5.36 in Ref(1)
  ActsMatrix<BoundScalar, eBoundSize, 4> positionJacobian;
  positionJacobian.setZero();
  // First row
  // positionJacobian(0, 0) = -sgnH * X / S;
  // positionJacobian(0, 1) = -sgnH * Y / S;

  // const double S2tanTh = S2 * tanTh;

  // // Second row
  // positionJacobian(1, 0) = rho * Y / S2tanTh;
  // positionJacobian(1, 1) = -rho * X / S2tanTh;
  // positionJacobian(1, 2) = 1.;

  // // Third row
  // positionJacobian(2, 0) = -Y / S2;
  // positionJacobian(2, 1) = X / S2;

  // // TODO: include timing in track linearization
  // // Last row
  // positionJacobian(5, 3) = 1;

  // std::cout << "old pos jacobian:\n" << positionJacobian << std::endl;

  // Fill momentum jacobian (E_k matrix), Eq. 5.37 in Ref(1)
  ActsMatrixD<eBoundSize, 3> momentumJacobian;
  momentumJacobian.setZero();

  double R = X * cosPhiV + Y * sinPhiV;
  double Q = X * sinPhiV - Y * cosPhiV;
  double dPhi = phiAtPCA - phiV;

  // // First row
  // momentumJacobian(0, 0) = -sgnH * rho * R / S;

  // double qOvSred = 1 - sgnH * Q / S;

  // momentumJacobian(0, 1) = qOvSred * rho / tanTh;
  // momentumJacobian(0, 2) = -qOvSred * rho / qOvP;

  // const double rhoOverS2 = rho / S2;

  // // Second row
  // momentumJacobian(1, 0) = (1 - rhoOverS2 * Q) * rho / tanTh;
  // momentumJacobian(1, 1) = (dPhi + rho * R / (S2tanTh * tanTh)) * rho;
  // momentumJacobian(1, 2) = (dPhi - rhoOverS2 * R) * rho / (qOvP * tanTh);

  // // Third row
  // momentumJacobian(2, 0) = rhoOverS2 * Q;
  // momentumJacobian(2, 1) = -rho * R / S2tanTh;
  // momentumJacobian(2, 2) = rhoOverS2 * R / qOvP;

  // // Last two rows:
  // momentumJacobian(3, 1) = 1.;
  // momentumJacobian(4, 2) = 1.;



  

  ActsMatrixD<4, 3> myTransform = ActsMatrixD<4, 3>::Zero();
  myTransform(0,0) = - sinTh * sinPhiV;
  myTransform(1,0) = sinTh * cosPhiV;
  myTransform(0,1) = cosTh * cosPhiV;
  myTransform(1,1) = cosTh * sinPhiV;
  myTransform(2,1) = -sinTh;
  myTransform(3,2) = 1.;

  ActsMatrix<BoundScalar, eBoundSize, 4> newMomentumJacobian;
  newMomentumJacobian = freeToBoundJacobian.block<6,4>(0,4);

  ActsMatrixD<6, 3> transformedMomJac = newMomentumJacobian * myTransform;

  // std::cout << "newMomentumJacobian: \n" << transformedMomJac << std::endl;

  // std::cout << "old mom jacobian:\n" << momentumJacobian << std::endl;


  positionJacobian = newPositionJacobian;
  momentumJacobian = transformedMomJac;
  predParamsAtPCA = paramsAtPCA;


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
