// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

namespace Acts {

/// @class GenericTrackLinearizer
/// Linearizes the measurement equation (dependance of track
/// parameters on the vertex position and track momentum at vertex)
/// at the vicinity of the user-provided linearization point.
///
/// The measurement equation is linearized in the following way:
///
/// q_k= A_k (x_k - x_0k) + B_k (p_k - p_0k) + c_k
///
/// where q_k are the parameters at perigee nearest to the lin point,
/// x_k is the position of the vertex, p_k the track momentum at the vertex,
/// and c_k is the constant term of expansion. A_k and B_k are matrices
/// of derivatives, denoted hereafter as "positionJacobian" and
/// "momentumJacobian" respectively.
///
/// Ref.(1) - CERN-THESIS-2010-027, Giacinto Piacquadio (Freiburg U.)
///
/// @tparam propagator_t Propagator type
/// @tparam propagator_options_t Propagator options type
template <typename propagator_t,
          typename propagator_options_t = PropagatorOptions<>>
class GenericTrackLinearizer {
 public:
  using Propagator_t = propagator_t;
  using BField_t = typename Propagator_t::Stepper::BField;

  /// @struct State struct
  struct State {
    /// @brief The state constructor
    ///
    /// @param mctx The magnetic field context
    State(const Acts::MagneticFieldContext& mctx) : fieldCache(mctx) {}
    /// Magnetic field cache
    typename BField_t::Cache fieldCache;
  };

  /// @brief Configuration struct
  struct Config {
    /// @ Config constructor if magnetic field is present
    ///
    /// @param bIn The magnetic field
    /// @param prop The propagator
    Config(const BField_t& bIn, std::shared_ptr<Propagator_t> prop)
        : bField(bIn), propagator(std::move(prop)) {}

    /// @brief Config constructor if BField_t == NullBField (no B-Field
    /// provided)
    ///
    /// @param prop The propagator
    template <typename T = BField_t,
              std::enable_if_t<std::is_same<T, NullBField>::value, int> = 0>
    Config(std::shared_ptr<Propagator_t> prop) : propagator(std::move(prop)) {}

    // The magnetic field
    BField_t bField;
    // The propagator
    std::shared_ptr<Propagator_t> propagator;
  };

  /// @brief Constructor
  ///
  /// @param config Configuration object
  GenericTrackLinearizer(const Config& config) : m_cfg(config) {}

  /// @brief Function that linearizes BoundTrackParameters at
  /// given linearization point
  ///
  /// @param params Parameters to linearize
  /// @param linPoint Linearization point
  /// @param gctx The geometry context
  /// @param mctx The magnetic field context
  /// @param state The state object
  ///
  /// @return Linearized track
  Result<LinearizedTrack> linearizeTrack(const BoundTrackParameters& params,
                                         const Vector4D& linPoint,
                                         const Acts::GeometryContext& gctx,
                                         const Acts::MagneticFieldContext& mctx,
                                         State& state) const;

 
struct NoPropagationAborter {
  NoPropagationAborter() = default;
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& /*state*/,
                  const stepper_t& /*unused*/) const {
    return true;
}
};

 private:
  /// Configuration object
  const Config m_cfg;
};

}  // namespace Acts

#include "GenericTrackLinearizer.ipp"
