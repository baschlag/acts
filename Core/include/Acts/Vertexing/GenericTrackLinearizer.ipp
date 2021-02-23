// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

template <typename propagator_t, typename propagator_options_t>
Acts::Result<Acts::LinearizedTrack> Acts::
    GenericTrackLinearizer<propagator_t, propagator_options_t>::linearizeTrack(
        const BoundTrackParameters& params, const Vector4D& linPoint,
        const Acts::GeometryContext& gctx,
        const Acts::MagneticFieldContext& mctx, State& /*unused*/) const {
  Vector3D linPointPos = VectorHelpers::position(linPoint);

  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(linPointPos);

      const std::shared_ptr<PerigeeSurface> perigeeSurface2 =
      Surface::makeShared<PerigeeSurface>(linPointPos + Vector3D(1e-2,1e-2,1e-2));

  // std::cout << "linPointPos: " << linPointPos.transpose() << std::endl;

  // Create propagator options
  auto logger = getDefaultLogger("HelTrkLinProp", Logging::INFO);
  propagator_options_t pOptions(gctx, mctx, LoggerWrapper{*logger});
  pOptions.direction = forward;
 
  // STEP 1: propagate to perigee surfac close to linPointPos

  // double radius(10.0), halfZ(1e5);
  // Translation3D translation{0., 0., 0.};
  // auto pTransform = Transform3D(translation);
  // auto cylinderSurface =
  //     Surface::makeShared<CylinderSurface>(pTransform, radius, halfZ);

  const BoundTrackParameters* endParamsBeampipe = nullptr;
  auto resultBeampipe = m_cfg.propagator->propagate(params, *perigeeSurface2, pOptions);
  //auto resultBeampipe = m_cfg.propagator->propagate(params, *cylinderSurface, pOptions);
  if (resultBeampipe.ok()) {
    endParamsBeampipe = (*resultBeampipe).endParameters.get();
  } else {
    return resultBeampipe.error();
  }
  
  // STEP 2: convert to free parameters

  using PropagatorOptions = PropagatorOptions<ActionList<>, AbortList<NoPropagationAborter>>;
  PropagatorOptions pOptionsAbort(gctx, mctx, LoggerWrapper{*logger});

  const FreeTrackParameters* endFreeParams = nullptr;
  auto freeresult = m_cfg.propagator->template propagate<FreeTrackParameters>(*endParamsBeampipe, pOptionsAbort);
  if (freeresult.ok()) {
    endFreeParams = (*freeresult).endParameters.get();
  } else {
    return freeresult.error();
  }

  // STEP 3: propagate back to perigee at linpoint, going from FreeParameters to BoundParameters --> retrieve jacobian
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
    return backBoundResult.error();
  }

  BoundVector paramsAtPCA = backBoundParams->parameters();
  Vector4D positionAtPCA = Vector4D::Zero();
  {
    auto pos = backBoundParams->position(gctx);
    positionAtPCA[ePos0] = pos[ePos0];
    positionAtPCA[ePos1] = pos[ePos1];
    positionAtPCA[ePos2] = pos[ePos2];
    positionAtPCA[eTime] = backBoundParams->time();
  }
  BoundSymMatrix parCovarianceAtPCA = *(backBoundParams->covariance());

  if (backBoundParams->covariance()->determinant() == 0) {
    // Use the original parameters
    paramsAtPCA = params.parameters();
    auto pos = backBoundParams->position(gctx);
    positionAtPCA[ePos0] = pos[ePos0];
    positionAtPCA[ePos1] = pos[ePos1];
    positionAtPCA[ePos2] = pos[ePos2];
    parCovarianceAtPCA = *(params.covariance());
  }

  // phiV and functions
  double phiV = paramsAtPCA(BoundIndices::eBoundPhi);
  double sinPhiV = std::sin(phiV);
  double cosPhiV = std::cos(phiV);
  double qOvP = paramsAtPCA(BoundIndices::eBoundQOverP);

  // theta and functions
  double th = paramsAtPCA(BoundIndices::eBoundTheta);
  const double sinTh = std::sin(th);
  const double cosTh = std::cos(th);

  Vector3D momentumAtPCA(phiV, th, qOvP);

  ActsMatrix<BoundScalar, eBoundSize, 4> positionJacobian;
  positionJacobian = freeToBoundJacobian.block<6,4>(0,0);

  ActsMatrixD<4, 3> freeToBoundMomentumTransform = ActsMatrixD<4, 3>::Zero();
  freeToBoundMomentumTransform(0,0) = - sinTh * sinPhiV;
  freeToBoundMomentumTransform(1,0) = sinTh * cosPhiV;
  freeToBoundMomentumTransform(0,1) = cosTh * cosPhiV;
  freeToBoundMomentumTransform(1,1) = cosTh * sinPhiV;
  freeToBoundMomentumTransform(2,1) = -sinTh;
  freeToBoundMomentumTransform(3,2) = 1.;

  ActsMatrix<BoundScalar, eBoundSize, 4> blockMomentumJacobian;
  blockMomentumJacobian = freeToBoundJacobian.block<6,4>(0,4);

  ActsMatrixD<6, 3> momentumJacobian = blockMomentumJacobian * freeToBoundMomentumTransform;

  momentumJacobian(0,0) = 0.;
  momentumJacobian(1,1) = 0.;

  // const term F(V_0, p_0) in Talyor expansion
  BoundVector constTerm = paramsAtPCA - positionJacobian * positionAtPCA -
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
