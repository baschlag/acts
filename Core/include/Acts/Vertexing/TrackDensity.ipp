// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <math.h>

template <typename input_track_t>
std::pair<double, double> Acts::TrackDensity<input_track_t>::globalMaximumWithWidth(
    State& state, const std::vector<const input_track_t*>& trackList,
    const std::function<BoundParameters(input_track_t)>& extractParameters) const {

  addAllTracks(state, trackList, extractParameters);

  double maxPosition = 0.;
  double maxDensity = 0.;
  double maxSecondDerivative = 0.;

  for (const auto& track : state.trackEntries) {
    double trialZ = track.z;
    double density = 0.;
    double firstDerivative = 0.;
    double secondDerivative = 0.;
    density = trackDensity(state, trialZ, firstDerivative, secondDerivative);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    updateMaximum(trialZ, density, secondDerivative, maxPosition, maxDensity,
                  maxSecondDerivative);
    trialZ += stepSize(density, firstDerivative, secondDerivative);
    density = trackDensity(state, trialZ, firstDerivative, secondDerivative);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    updateMaximum(trialZ, density, secondDerivative, maxPosition, maxDensity,
                  maxSecondDerivative);
    trialZ += stepSize(density, firstDerivative, secondDerivative);
    density = trackDensity(state, trialZ, firstDerivative, secondDerivative);
    if (secondDerivative >= 0. || density <= 0.) {
      continue;
    }
    updateMaximum(trialZ, density, secondDerivative, maxPosition, maxDensity,
                  maxSecondDerivative);
  }

  return std::make_pair(maxPosition,
                        std::sqrt(-(maxDensity / maxSecondDerivative)));
}

template <typename input_track_t>
double Acts::TrackDensity<input_track_t>::globalMaximum(State& state,
  const std::vector<const input_track_t*>& trackList,
  const std::function<BoundParameters(input_track_t)>& extractParameters) const {
  return globalMaximumWithWidth(state, trackList, extractParameters).first;
}

template <typename input_track_t>
void Acts::TrackDensity<input_track_t>::addAllTracks(State& state, 
    const std::vector<const input_track_t*>& trackList,
    const std::function<BoundParameters(input_track_t)>& extractParameters) const {
  for (auto trk : trackList) {
    addTrack(state, trk, extractParameters);
  }
}

template <typename input_track_t>
void Acts::TrackDensity<input_track_t>::addTrack(State& state, const input_track_t* trk,
  const std::function<BoundParameters(input_track_t)>& extractParameters) const {
  const BoundParameters& boundParams = extractParameters(*trk);

  // Get required track parameters
  const double d0 = boundParams.parameters()[ParID_t::eLOC_D0];
  const double z0 = boundParams.parameters()[ParID_t::eLOC_Z0];
  // Get track covariance
  const auto perigeeCov = *(boundParams.covariance());
  const double covDD = perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_D0);
  const double covZZ = perigeeCov(ParID_t::eLOC_Z0, ParID_t::eLOC_Z0);
  const double covDZ = perigeeCov(ParID_t::eLOC_D0, ParID_t::eLOC_Z0);
  const double covDeterminant = covDD * covZZ - covDZ * covDZ;

  // Do track selection based on track cov matrix and m_cfg.d0SignificanceCut
  if ((covDD <= 0) || (d0 * d0 / covDD > m_cfg.d0SignificanceCut) || (covZZ <= 0) ||
      (covDeterminant <= 0)) {
    return;
  }

  // Calculate track density quantities
  double constantTerm =
      -(d0 * d0 * covZZ + z0 * z0 * covDD + 2. * d0 * z0 * covDZ) /
      (2. * covDeterminant);
  const double linearTerm =
      (d0 * covDZ + z0 * covDD) /
      covDeterminant;  // minus signs and factors of 2 cancel...
  const double quadraticTerm = -covDD / (2. * covDeterminant);
  double discriminant =
      linearTerm * linearTerm -
      4. * quadraticTerm * (constantTerm + 2. * m_cfg.z0SignificanceCut);
  if (discriminant < 0) {
    return;
  }

  // Add the track to the current maps in the state
  discriminant = std::sqrt(discriminant);
  const double zMax = (-linearTerm - discriminant) / (2. * quadraticTerm);
  const double zMin = (-linearTerm + discriminant) / (2. * quadraticTerm);
  constantTerm -= std::log(2. * M_PI * std::sqrt(covDeterminant));

  state.trackEntries.emplace_back(trk, z0, constantTerm, linearTerm, quadraticTerm,
                                  zMin, zMax);
}

template <typename input_track_t>
double Acts::TrackDensity<input_track_t>::trackDensity(State& state, double z) const {
  double firstDerivative = 0;
  double secondDerivative = 0;
  return trackDensity(state, z, firstDerivative, secondDerivative);
}

template <typename input_track_t>
double Acts::TrackDensity<input_track_t>::trackDensity(State& state, double z,
                                        double& firstDerivative,
                                        double& secondDerivative) const {
  TrackDensityStore densityResult(z);
  for (const auto& trackEntry : state.trackEntries) {
    if(state.testMap.find(std::make_pair(trackEntry.trkPtr, z)) == state.testMap.end()){
      // Pair(track, trial z position) not in cache yet
      typename TrackDensityStore::Cache cache;
      densityResult.addTrackToDensity(trackEntry, cache);
      state.testMap[std::make_pair(trackEntry.trkPtr, z)] = cache;
    }
    else{
      // Pair(track, trial z position) found in cache, use cached values
      typename TrackDensityStore::Cache cache = state.testMap.at(std::make_pair(trackEntry.trkPtr, z));
      densityResult.updateValues(cache);
    }
    
    
  }
  firstDerivative = densityResult.firstDerivative();
  secondDerivative = densityResult.secondDerivative();

  return densityResult.density();
}

template <typename input_track_t>
void Acts::TrackDensity<input_track_t>::updateMaximum(double newZ, double newValue,
                                       double newSecondDerivative, double& maxZ,
                                       double& maxValue,
                                       double& maxSecondDerivative) const {
  if (newValue > maxValue) {
    maxZ = newZ;
    maxValue = newValue;
    maxSecondDerivative = newSecondDerivative;
  }
}

template <typename input_track_t>
double Acts::TrackDensity<input_track_t>::stepSize(double y, double dy, double ddy) const {
  return (m_cfg.isGaussianShaped ? (y * dy) / (dy * dy - y * ddy) : -dy / ddy);
}

template <typename input_track_t>
void Acts::TrackDensity<input_track_t>::TrackDensityStore::updateValues(Cache& cache){
    m_density += cache.deltaDensity;
    m_firstDerivative += cache.deltaFirstDerivative;
    m_secondDerivative += cache.deltaSecondDerivative;
}

template <typename input_track_t>
void Acts::TrackDensity<input_track_t>::TrackDensityStore::addTrackToDensity(
    const TrackEntry& entry, Cache& cache) {
  // Take track only if it's within bounds
  if (entry.lowerBound < m_z && m_z < entry.upperBound) {
    cache.deltaDensity = std::exp(entry.c0 + m_z * (entry.c1 + m_z * entry.c2));
    double qPrime = entry.c1 + 2. * m_z * entry.c2;
    cache.deltaFirstDerivative = cache.deltaDensity * qPrime;
    cache.deltaSecondDerivative = 2. * entry.c2 * cache.deltaDensity + qPrime * cache.deltaFirstDerivative;

    updateValues(cache);

  }
}
