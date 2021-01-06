// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t>
void Acts::KalmanVertexTrackUpdater::update(TrackAtVertex<input_track_t>& track,
                                            const Vertex<input_track_t>& vtx) {
  const Vector3D vtxPos = vtx.fullPosition().template head<3>();

  // Get the linearized track
  const LinearizedTrack& linTrack = track.linearizedState;

  // Check if linearized state exists
  if (linTrack.covarianceAtPCA.determinant() == 0.) {
    // Track has no linearized state, returning w/o update
    return;
  }

  // Retrieve linTrack information
  const ActsMatrixD<5, 3> posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  const ActsMatrixD<5, 3> momJac = linTrack.momentumJacobian.block<5, 3>(0, 0);
  const ActsVectorD<5> trkParams = linTrack.parametersAtPCA.head<5>();
  const ActsSymMatrixD<5> trkParamWeight =
      linTrack.weightAtPCA.block<5, 5>(0, 0);

  const Vector4D fullVtxPos = vtx.fullPosition();
  const ActsVectorD<6> fullTrkParams = linTrack.parametersAtPCA;
  const FreeToBoundMatrix& freeToBoundJacobian = linTrack.freeToBoundJacobian;
  const ActsMatrixD<6, 4> fullPosJac = freeToBoundJacobian.block<6,4>(0,0);
  const ActsMatrixD<6, 4> fullMomJac = freeToBoundJacobian.block<6,4>(0,4);
  ActsSymMatrixD<6> fullTrkParamWeight = linTrack.weightAtPCA;

  // TODO REMOVE... discard time for now
  fullTrkParamWeight.row(5).setZero();
  fullTrkParamWeight.col(5).setZero();

  // Calculate S matrix
  ActsSymMatrixD<4> sMatFull =
      (fullMomJac.transpose() * (fullTrkParamWeight * fullMomJac)).inverse();

  // Calculate S matrix
  ActsSymMatrixD<3> sMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  const ActsVectorD<6> fullResidual = linTrack.newConstantTerm;
  const ActsVectorD<5> residual = linTrack.constantTerm.head<5>();

  std::cout << "fullMomJac: " << fullMomJac << std::endl;
  std::cout << "fullPosJac: " << fullPosJac << std::endl;
  std::cout << "old PosJac: " << posJac << std::endl;
  std::cout << "fullTrkParamWeight: " << fullTrkParamWeight << std::endl;
  std::cout << "fullVtxPos: " << fullVtxPos << std::endl;
  std::cout << "fullTrkParams: " << fullTrkParams << std::endl;
  std::cout << "sMatFull w/o INV: " << fullMomJac.transpose() * (fullTrkParamWeight * fullMomJac) << std::endl;
  std::cout << "sMatFull: " << sMatFull << std::endl;

  // Refit track momentum
  Vector4D fullNewTrkMomentum = sMatFull * fullMomJac.transpose() * fullTrkParamWeight *
                            (fullTrkParams - fullResidual - fullPosJac * fullVtxPos);

  std::cout << "residual: " << residual << std::endl;
  std::cout << "fullResidual: " << fullResidual << std::endl;

  // Refit track momentum
  Vector3D newTrkMomentum = sMat * momJac.transpose() * trkParamWeight *
                            (trkParams - residual - posJac * vtxPos);

  // Refit track parameters
  BoundVector newFullTrkParams(BoundVector::Zero());

  // Refit track parameters
  BoundVector newTrkParams(BoundVector::Zero());

  // Get phi and theta and correct for possible periodicity changes
  auto correctedPhiTheta =
      Acts::detail::ensureThetaBounds(newTrkMomentum(0), newTrkMomentum(1));

  double Tx = fullNewTrkMomentum(0);
  double Ty = fullNewTrkMomentum(1);
  double Tz = fullNewTrkMomentum(2);
  double qOverP = fullNewTrkMomentum(3); 

  std::cout << "Tx,Ty,Tz,q/p: " << Tx << "," << Ty << "," << Tz << "," << qOverP << std::endl; 
 
  newFullTrkParams(BoundIndices::eBoundPhi) = std::atan2(Ty, Tx) ;     // phi
  newFullTrkParams(BoundIndices::eBoundTheta) = std::atan2(std::sqrt(Tx*Tx + Ty * Ty), Tz);  // theta
  newFullTrkParams(BoundIndices::eBoundQOverP) = qOverP;      // qOverP

  newTrkParams(BoundIndices::eBoundPhi) = correctedPhiTheta.first;     // phi
  newTrkParams(BoundIndices::eBoundTheta) = correctedPhiTheta.second;  // theta
  newTrkParams(BoundIndices::eBoundQOverP) = newTrkMomentum(2);        // qOverP

  std::cout << "COMPARE:" << std::endl;
  std::cout << newTrkParams(BoundIndices::eBoundPhi) << " = " << newFullTrkParams(BoundIndices::eBoundPhi) << std::endl;
  std::cout << newTrkParams(BoundIndices::eBoundTheta) << " = " << newFullTrkParams(BoundIndices::eBoundTheta) << std::endl;
  std::cout << newTrkParams(BoundIndices::eBoundQOverP) << " = " << newFullTrkParams(BoundIndices::eBoundQOverP) << std::endl;
  std::cout << std::endl;

  // Vertex covariance and weight matrices
  const ActsSymMatrixD<4> fullVtxCov = vtx.fullCovariance();
  const ActsSymMatrixD<4> fullVtxWeight = fullVtxCov.inverse();

  // Vertex covariance and weight matrices
  const ActsSymMatrixD<3> vtxCov =
      vtx.fullCovariance().template block<3, 3>(0, 0);
  const ActsSymMatrixD<3> vtxWeight = vtxCov.inverse();

  // New track covariance matrix
  ActsSymMatrixD<4> fullNewTrkCov =
      -fullVtxCov * fullPosJac.transpose() * fullTrkParamWeight * fullMomJac * sMatFull;

  // New track covariance matrix
  ActsSymMatrixD<3> newTrkCov =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * sMat;

  KalmanVertexUpdater::MatrixCache matrixCache;

  // Now determine the smoothed chi2 of the track in the following
  KalmanVertexUpdater::updatePosition<input_track_t>(
      vtx, linTrack, track.trackWeight, -1, matrixCache);

  // Corresponding weight matrix
  const ActsSymMatrixD<3>& reducedVtxWeight = matrixCache.newVertexWeight;

  // Difference in positions
  Vector3D posDiff = vtx.position() - matrixCache.newVertexPos;

  // Get smoothed params
  ActsVectorD<5> smParams =
      trkParams - (residual + posJac * vtx.fullPosition().template head<3>() +
                   momJac * newTrkMomentum);

  // New chi2 to be set later
  double chi2 = posDiff.dot(reducedVtxWeight * posDiff) +
                smParams.dot(trkParamWeight * smParams);

  // Not yet 4d ready. This can be removed together will all head<> statements,
  // once time is consistently introduced to vertexing
  ActsMatrixD<4, 3> newFullTrkCov(ActsMatrixD<4, 3>::Zero());
  newFullTrkCov.block<3, 3>(0, 0) = newTrkCov;

  SymMatrix4D vtxFullWeight(SymMatrix4D::Zero());
  vtxFullWeight.block<3, 3>(0, 0) = vtxWeight;

  SymMatrix4D vtxFullCov(SymMatrix4D::Zero());
  vtxFullCov.block<3, 3>(0, 0) = vtxCov;

  const Acts::BoundMatrix fullPerTrackCov = detail::createFullTrackCovariance(
      sMat, newFullTrkCov, vtxFullWeight, vtxFullCov, newTrkParams);

  // Create new refitted parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  BoundTrackParameters refittedPerigee = BoundTrackParameters(
      perigeeSurface, newTrkParams, std::move(fullPerTrackCov));

  // Set new properties
  track.fittedParams = refittedPerigee;
  track.chi2Track = chi2;
  track.ndf = 2 * track.trackWeight;

  return;
}

inline Acts::BoundMatrix
Acts::KalmanVertexTrackUpdater::detail::createFullTrackCovariance(
    const SymMatrix3D& sMat, const ActsMatrixD<4, 3>& newTrkCov,
    const SymMatrix4D& vtxWeight, const SymMatrix4D& vtxCov,
    const BoundVector& newTrkParams) {
  // Now new momentum covariance
  ActsSymMatrixD<3> momCov =
      sMat + (newTrkCov.block<3, 3>(0, 0)).transpose() *
                 (vtxWeight.block<3, 3>(0, 0) * newTrkCov.block<3, 3>(0, 0));

  // Full (x,y,z,phi, theta, q/p) covariance matrix
  // To be made 7d again after switching to (x,y,z,phi, theta, q/p, t)
  ActsSymMatrixD<6> fullTrkCov(ActsSymMatrixD<6>::Zero());

  fullTrkCov.block<3, 3>(0, 0) = vtxCov.block<3, 3>(0, 0);
  fullTrkCov.block<3, 3>(0, 3) = newTrkCov.block<3, 3>(0, 0);
  fullTrkCov.block<3, 3>(3, 0) = (newTrkCov.block<3, 3>(0, 0)).transpose();
  fullTrkCov.block<3, 3>(3, 3) = momCov;

  // Combined track jacobian
  ActsMatrixD<5, 6> trkJac(ActsMatrixD<5, 6>::Zero());

  // First row
  trkJac(0, 0) = -std::sin(newTrkParams[2]);
  trkJac(0, 1) = std::cos(newTrkParams[2]);

  double tanTheta = std::tan(newTrkParams[3]);

  // Second row
  trkJac(1, 0) = -trkJac(0, 1) / tanTheta;
  trkJac(1, 1) = trkJac(0, 0) / tanTheta;

  trkJac.block<4, 4>(1, 2) = ActsSymMatrixD<4>::Identity();

  // Full perigee track covariance
  BoundMatrix fullPerTrackCov(BoundMatrix::Identity());
  fullPerTrackCov.block<5, 5>(0, 0) =
      (trkJac * (fullTrkCov * trkJac.transpose()));

  return fullPerTrackCov;
}
