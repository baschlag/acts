// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/FsmwMode1dFinder.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitterConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "Acts/Vertexing/ZScanVertexFinder.hpp"

namespace Acts {

///
/// @tparam vfitter_t Vertex fitter type
/// @tparam sfinder_t Seed finder type
template <typename vfitter_t, typename sfinder_t>
class AntiKtVertexFinder {

 public:
  using InputTrack_t = typename vfitter_t::InputTrack_t;

  /// @struct Config Configuration struct
  struct Config {

  	double maxZsearchWindow = 4_mm;
  	
  	double radiusParameter = 500_um;

  };

  /// @struct State State struct
  struct State {};

  /// @brief Constructor used if InputTrack_t type == BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param logger The logging instance
  template <
      typename T = InputTrack_t,
      std::enable_if_t<std::is_same<T, BoundTrackParameters>::value, int> = 0>
  AntiKtVertexFinder(Config& cfg)
      : m_cfg(std::move(cfg)),
        m_extractParameters([](T params) { return params; }){}

  /// @brief Constructor for user-defined InputTrack_t type =!
  /// BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from InputTrack_t
  /// object
  /// @param logger The logging instance
  AntiKtVertexFinder(Config& cfg,
                        std::function<BoundTrackParameters(InputTrack_t)> func)
      : m_cfg(std::move(cfg)),
        m_extractParameters(func) {}

  /// @brief Finds vertices corresponding to input trackVector
  ///
  /// @param trackVector Input tracks
  /// @param vertexingOptions Vertexing options
  /// @param state State for fulfilling interfaces
  ///
  /// @return Collection of vertices found by finder
  Result<std::vector<Vertex<InputTrack_t>>> find(
      const std::vector<const InputTrack_t*>& trackVector,
      const VertexingOptions<InputTrack_t>& vertexingOptions) const;

 private:

 	struct PseudoTrack
 	{
 		// If breaking condition is reached, declare PseudoTrack as vertex
 		bool isVertex = false;
 		// If track is merged with another track, it becomes invalid
 		bool isValid = true;
 		double sumPt = 0;
 		std::vector<double> zValues;
 		double zMean = 0;
 		// Sum of squared sigma_z0 values
 		double sumSigZ0squared = 0;
 		std::vector<const InputTrack_t*> tracks;
 	};


  /// Configuration object
  const Config m_cfg;

  /// @brief Function to extract track parameters,
  /// InputTrack_t objects are BoundTrackParameters by default, function to be
  /// overwritten to return BoundTrackParameters for other InputTrack_t objects.
  ///
  /// @param InputTrack_t object to extract track parameters from
  std::function<BoundTrackParameters(InputTrack_t)> m_extractParameters;
};

}  // namespace Acts

#include "AntiKtVertexFinder.ipp"
