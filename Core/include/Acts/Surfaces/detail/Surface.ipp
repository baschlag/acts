// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline Vector3D Surface::center(const GeometryContext& gctx) const {
  // fast access via tranform matrix (and not translation())
  auto tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 3), tMatrix(1, 3), tMatrix(2, 3));
}

inline Vector3D Surface::normal(const GeometryContext& gctx,
                                const Vector3D& /*unused*/) const {
  return normal(gctx, s_origin2D);
}

inline const Transform3D& Surface::transform(
    const GeometryContext& gctx) const {
  if (m_associatedDetElement != nullptr) {
    return m_associatedDetElement->transform(gctx);
  }
  return m_transform;
}

inline bool Surface::insideBounds(const Vector2D& lposition,
                                  const BoundaryCheck& bcheck) const {
  return bounds().inside(lposition, bcheck);
}

inline RotationMatrix3D Surface::referenceFrame(
    const GeometryContext& gctx, const Vector3D& /*unused*/,
    const Vector3D& /*unused*/) const {
  return transform(gctx).matrix().block<3, 3>(0, 0);
}

inline void Surface::initJacobianToGlobal(const GeometryContext& gctx,
                                          BoundToFreeMatrix& jacobian,
                                          const Vector3D& position,
                                          const Vector3D& direction,
                                          const BoundVector& /*pars*/) const {
  // The trigonometry required to convert the direction to spherical
  // coordinates and then compute the sines and cosines again can be
  // surprisingly expensive from a performance point of view.
  //
  // Here, we can avoid it because the direction is by definition a unit
  // vector, with the following coordinate conversions...
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)

  // ...which we can invert to directly get the sines and cosines:
  const double cos_theta = z;
  const double sin_theta = sqrt(x * x + y * y);
  const double inv_sin_theta = 1. / sin_theta;
  const double cos_phi = x * inv_sin_theta;
  const double sin_phi = y * inv_sin_theta;
  // retrieve the reference frame
  const auto rframe = referenceFrame(gctx, position, direction);
  // the local error components - given by reference frame
  jacobian.topLeftCorner<3, 2>() = rframe.topLeftCorner<3, 2>();
  // the time component
  jacobian(3, eBoundTime) = 1;
  // the momentum components
  jacobian(4, eBoundPhi) = (-sin_theta) * sin_phi;
  jacobian(4, eBoundTheta) = cos_theta * cos_phi;
  jacobian(5, eBoundPhi) = sin_theta * cos_phi;
  jacobian(5, eBoundTheta) = cos_theta * sin_phi;
  jacobian(6, eBoundTheta) = (-sin_theta);
  jacobian(7, eBoundQOverP) = 1;
}

inline void Surface::initJacobianToLocal(const GeometryContext& gctx,
                                         FreeToBoundMatrix& jacobian,
                                         const Vector3D& position,
                                         const Vector3D& direction) const {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // The measurement frame of the surface
  RotationMatrix3D rframeT =
      referenceFrame(gctx, position, direction).transpose();
  // given by the refernece frame
  jacobian.block<2, 3>(0, 0) = rframeT.block<2, 3>(0, 0);
  // Time component
  jacobian(eBoundTime, 3) = 1;
  // Directional and momentum elements for reference frame surface
  jacobian(eBoundPhi, 4) = -sinPhi * invSinTheta;
  jacobian(eBoundPhi, 5) = cosPhi * invSinTheta;
  jacobian(eBoundTheta, 4) = cosPhi * cosTheta;
  jacobian(eBoundTheta, 5) = sinPhi * cosTheta;
  jacobian(eBoundTheta, 6) = -sinTheta;
  jacobian(eBoundQOverP, 7) = 1;
}

inline FreeRowVector Surface::freeToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global position
  const auto position = parameters.head<3>();
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The measurement frame of the surface
  const RotationMatrix3D rframe = referenceFrame(gctx, position, direction);
  // The measurement frame z axis
  const Vector3D refZAxis = rframe.col(2);
  // Cosine of angle between momentum direction and measurement frame z axis
  const double dz = refZAxis.dot(direction);
  // Initialize the derivative
  FreeRowVector freeToPath = FreeRowVector::Zero();
  freeToPath.head<3>() = -1.0 * refZAxis.transpose() / dz;
  return freeToPath;
}

inline const FreeRowVector Surface::derivativeFactors(
    const GeometryContext& /*unused*/, const Vector3D& /*unused*/,
    const Vector3D& dir, const RotationMatrix3D& rft,
    const FreeMatrix& jac) const {
  // Create the normal and scale it with the projection onto the direction
  ActsRowVectorD<3> norm_vec = rft.template block<1, 3>(2, 0);
  norm_vec /= (norm_vec * dir);
  // calculate the s factors
  return (norm_vec * jac.topLeftCorner<3, FreeParsDim>());
}
inline const BoundRowVector Surface::derivativeFactors(
    const GeometryContext& /*unused*/, const Vector3D& /*unused*/,
    const Vector3D& dir, const RotationMatrix3D& rft,
    const BoundToFreeMatrix& jac) const {
  // Create the normal and scale it with the projection onto the direction
  ActsRowVectorD<3> norm_vec = rft.template block<1, 3>(2, 0);
  norm_vec /= (norm_vec * dir);
  // calculate the s factors
  return (norm_vec * jac.topLeftCorner<3, BoundParsDim>());
}

inline const DetectorElementBase* Surface::associatedDetectorElement() const {
  return m_associatedDetElement;
}

inline const Layer* Surface::associatedLayer() const {
  return (m_associatedLayer);
}

inline const ISurfaceMaterial* Surface::surfaceMaterial() const {
  return m_surfaceMaterial.get();
}

inline const std::shared_ptr<const ISurfaceMaterial>&
Surface::surfaceMaterialSharedPtr() const {
  return m_surfaceMaterial;
}

inline void Surface::assignSurfaceMaterial(
    std::shared_ptr<const ISurfaceMaterial> material) {
  m_surfaceMaterial = std::move(material);
}

inline void Surface::associateLayer(const Layer& lay) {
  m_associatedLayer = (&lay);
}
