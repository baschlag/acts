// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <set>
#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {

/// @class TrackDensity
///
/// @brief Class to model tracks as 2D density functions based on
/// their d0 and z0 perigee parameters (mean value) and covariance
/// matrices (determining the width of the function)
template<typename input_track_t>
class TrackDensity {
 public:
  /// @brief Struct to store information for a single track
  struct TrackEntry {
    // Default constructor
    TrackEntry() = default;
    // Constructor initializing all members
    TrackEntry(const input_track_t* trk, double z, double c0in, double c1in, double c2in, double zMin,
               double zMax)
        : trkPtr(trk),
          z(z),
          c0(c0in),
          c1(c1in),
          c2(c2in),
          lowerBound(zMin),
          upperBound(zMax) {}

    const input_track_t* trkPtr;

    double z = 0;
    // Cached information for a single track
    // z-independent term in exponent
    double c0 = 0;
    // linear coefficient in exponent
    double c1 = 0;
    // quadratic coefficient in exponent
    double c2 = 0;
    // The lower bound
    double lowerBound = 0;
    // The upper bound
    double upperBound = 0;
  };

private:
  /// @class TrackDensityStore
  /// @brief Helper class to evaluate and store track density at specific position
  class TrackDensityStore {
   public:
    struct Cache
    {
      double deltaDensity = 0;
      double deltaFirstDerivative = 0;
      double deltaSecondDerivative = 0;
    };
    // Initialise at the z coordinate at which the density is to be evaluated
    TrackDensityStore(double z_coordinate) : m_z(z_coordinate) {}

    // Add the contribution of a single track to the density
    void addTrackToDensity(const TrackEntry& entry, Cache& cache);

    // TODO
    void updateValues(Cache& cache);

    // Retrieve the density and its derivatives
    inline double density() const { return m_density; }
    inline double firstDerivative() const { return m_firstDerivative; }
    inline double secondDerivative() const { return m_secondDerivative; }

   private:
    // Store density and derivatives for z position m_z
    double m_z;
    double m_density{0};
    double m_firstDerivative{0};
    double m_secondDerivative{0};
  };

public:
  /// @brief The Config struct
  struct Config {

    Config(double d0Sig = 3.5, double z0Sig = 12.)
    : d0MaxSignificance(d0Sig)
    , z0MaxSignificance(z0Sig)
    , d0SignificanceCut(d0Sig * d0Sig)
    , z0SignificanceCut(z0Sig * z0Sig)
    {}

    // Assumed shape of density function:
    // Gaussian (true) or parabolic (false)
    bool isGaussianShaped = true;

    // Maximum d0 impact parameter significance to use a track
    double d0MaxSignificance;
    // Maximum z0 impact parameter significance to use a track
    double z0MaxSignificance;
    // Correspondig cut values
    const double d0SignificanceCut;
    const double z0SignificanceCut;
  };

  /// @brief The State struct
  struct State {
    // Constructor with size track map
    State(unsigned int nTracks) { trackEntries.reserve(nTracks); }
    // Vector to cache track information
    std::vector<TrackEntry> trackEntries;

    std::map<std::pair<const input_track_t*, double>, typename TrackDensityStore::Cache> testMap;
  };

  /// Default constructor
  TrackDensity() = default;

  /// Constructor with config
  TrackDensity(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Calculates z position of global maximum with Gaussian width
  /// for density function.
  /// Strategy:
  /// The global maximum must be somewhere near a track.
  /// Since we can calculate the first and second derivatives, at each point we
  /// can determine a) whether the function is curved up (minimum) or down
  /// (maximum) b) the distance to nearest maximum, assuming either Newton
  /// (parabolic) or Gaussian local behavior.
  /// For each track where the second derivative is negative, find step to
  /// nearest maximum, take that step and then do one final refinement. The
  /// largest density encountered in this procedure (after checking all tracks)
  /// is considered the maximum.
  ///
  /// @param state The track density state
  /// TODO
  ///
  /// @return Pair of position of global maximum and Gaussian width
  std::pair<double, double> globalMaximumWithWidth(State& state,
    const std::vector<const input_track_t*>& trackList,
    const std::function<BoundParameters(input_track_t)>& m_extractParameters) const;

  /// @brief Calculates the z position of the global maximum
  ///
  /// @param state The track density state
  /// TODO
  ///
  /// @return z position of the global maximum
  double globalMaximum(State& state,
    const std::vector<const input_track_t*>& trackList,
    const std::function<BoundParameters(input_track_t)>& m_extractParameters) const;

 private:
  /// The configuration
  Config m_cfg;

  // TODO
  void addAllTracks(State& state, 
    const std::vector<const input_track_t*>& trackList,
    const std::function<BoundParameters(input_track_t)>& m_extractParameters) const;

  /// @brief Add a track to the set being considered
  ///
  /// @param state The track density state
  /// @param trk Track parameters.
  /// TODO
  void addTrack(State& state, const input_track_t* trk,
    const std::function<BoundParameters(input_track_t)>& m_extractParameters) const;

  /// @brief Evaluate the density function at the specified
  /// coordinate along the beamline
  ///
  /// @param state The track density state
  /// @param z z-position along the beamline
  ///
  /// @return The track density value
  double trackDensity(State& state, double z) const;

  /// @brief Evaluate the density function and its two first
  /// derivatives at the specified coordinate along the beamline
  ///
  /// @param state The track density state
  /// @param z z-position along the beamline
  /// @param[out] firstDerivative The first derivative
  /// @param[out] secondDerivative The second derivative
  ///
  /// @return The track density value
  double trackDensity(State& state, double z, double& firstDerivative,
                      double& secondDerivative) const;

  /// @brief Update the current maximum values
  ///
  /// @param newZ The new z value
  /// @param newValue The new value at z position
  /// @param newSecondDerivative The new second derivative
  /// @param[out] maxZ The max z value
  /// @param[out] maxValue The max value at z position
  /// @param[out] maxSecondDerivative The max second derivative
  void updateMaximum(double newZ, double newValue, double newSecondDerivative,
                     double& maxZ, double& maxValue,
                     double& maxSecondDerivative) const;

  /// @brief Calculates the step size
  ///
  /// @param y Position value
  /// @param dy First derivative
  /// @param ddy Second derivative
  ///
  /// @return The step size
  double stepSize(double y, double dy, double ddy) const;

};

#include "Acts/Vertexing/TrackDensity.ipp"
}  // namespace Acts