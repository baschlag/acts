// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <bitset>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "Acts/EventData/TrackParametersBase.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

// forward declarations
template <typename source_link_t>
class MultiTrajectory;
class Surface;

namespace detail_lt {
  /// Either type T or const T depending on the boolean.
  template <typename T, bool select>
  using ConstIf = std::conditional_t<select, const T, T>;
  /// wrapper for a dynamic Eigen type that adds support for automatic growth
  ///
  /// \warning Assumes the underlying storage has a fixed number of rows
  template <typename Storage, size_t kSizeIncrement>
  struct GrowableColumns
  {

    /// Make sure storage for @p n additional columns is allocated. Will update
    /// the size of the container accordingly. The indices added by this call
    /// can safely be written to.
    /// @param n Number of columns to add, defaults to 1.
    /// @return View into the last allocated column
    auto
    addCol(size_t n = 1)
    {
      size_t index = m_size + (n - 1);
      while (capacity() <= index) {
        data.conservativeResize(Eigen::NoChange, data.cols() + kSizeIncrement);
      }
      m_size = index + 1;

      // @TODO: do this or not? If we assume this happens only when something is
      // written, the expectation is that everything is zero
      data.col(index).setZero();

      return data.col(index);
    }

    /// Writable access to a column w/o checking its existence first.
    auto
    col(size_t index)
    {
      return data.col(index);
    }

    /// Read-only access to a column w/o checking its existence first.
    auto
    col(size_t index) const
    {
      return data.col(index);
    }

    /// Return the current allocated storage capacity
    size_t
    capacity() const
    {
      return static_cast<size_t>(data.cols());
    }

    size_t
    size() const
    {
      return m_size;
    }

  private:
    Storage data;
    size_t  m_size{0};
  };

  /// Type construction helper for coefficients and associated covariances.
  template <size_t Size, bool ReadOnlyMaps = true>
  struct Types
  {
    enum {
      Flags         = Eigen::ColMajor | Eigen::AutoAlign,
      SizeIncrement = 8,
    };
    using Scalar = double;
    // single items
    using Coefficients    = Eigen::Matrix<Scalar, Size, 1, Flags>;
    using Covariance      = Eigen::Matrix<Scalar, Size, Size, Flags>;
    using CoefficientsMap = Eigen::Map<ConstIf<Coefficients, ReadOnlyMaps>>;
    using CovarianceMap   = Eigen::Map<ConstIf<Covariance, ReadOnlyMaps>>;
    // storage of multiple items in flat arrays
    using StorageCoefficients
        = GrowableColumns<Eigen::Array<Scalar, Size, Eigen::Dynamic, Flags>,
                          SizeIncrement>;
    using StorageCovariance
        = GrowableColumns<Eigen::
                              Array<Scalar, Size * Size, Eigen::Dynamic, Flags>,
                          SizeIncrement>;
  };

  struct IndexData
  {
    using IndexType = uint16_t;

    static constexpr IndexType kInvalid = UINT16_MAX;

    const Surface& surface;
    IndexType      iprevious  = kInvalid;
    IndexType      ipredicted = kInvalid;
    IndexType      ifiltered  = kInvalid;
    IndexType      ismoothed  = kInvalid;
    IndexType      ijacobian  = kInvalid;
    IndexType      iprojector = kInvalid;

    IndexType iuncalibrated         = kInvalid;
    IndexType icalibrated           = kInvalid;
    IndexType icalibratedsourcelink = kInvalid;
    IndexType measdim               = 0;
  };

  /// Proxy object to access a single point on the trajectory.
  ///
  /// @tparam source_link_t Type to link back to an original measurement
  /// @tparam N         Number of track parameters
  /// @tparam M         Maximum number of measurements
  /// @tparam ReadOnly  true for read-only access to underlying storage
  template <typename source_link_t, size_t N, size_t M, bool ReadOnly = true>
  class TrackStateProxy
  {
  public:
    using SourceLink            = source_link_t;
    using Parameters            = typename Types<N, ReadOnly>::CoefficientsMap;
    using Covariance            = typename Types<N, ReadOnly>::CovarianceMap;
    using Measurement           = typename Types<M, ReadOnly>::CoefficientsMap;
    using MeasurementCovariance = typename Types<M, ReadOnly>::CovarianceMap;

    // as opposed to the types above, this is an actual Matrix (rather than a
    // map)
    // @TODO: Does not copy flags, because this fails: can't have col major row
    // vector, but that's required for 1xN projection matrices below.
    constexpr static auto ProjectorFlags = Eigen::RowMajor | Eigen::AutoAlign;
    using Projector
        = Eigen::Matrix<typename Covariance::Scalar, M, N, ProjectorFlags>;
    using EffectiveProjector = Eigen::Matrix<typename Projector::Scalar,
                                             Eigen::Dynamic,
                                             Eigen::Dynamic,
                                             ProjectorFlags,
                                             M,
                                             N>;

    /// Index within the trajectory.
    size_t
    index() const
    {
      return m_istate;
    }

    /// Reference surface.
    const Surface&
    referenceSurface() const
    {
      return m_data.surface;
    }

    /// Track parameters vector.
    Parameters
    parameters() const;

    /// Track parameters covariance matrix.
    Covariance
    covariance() const;

    /// Track parameters vector.
    Parameters
    predicted() const;

    /// Track parameters covariance matrix.
    Covariance
    predictedCovariance() const;

    bool
    hasPredicted() const
    {
      return m_data.ipredicted != IndexData::kInvalid;
    }

    /// Track parameters vector.
    Parameters
    filtered() const;

    /// Track parameters covariance matrix.
    Covariance
    filteredCovariance() const;

    bool
    hasFiltered() const
    {
      return m_data.ifiltered != IndexData::kInvalid;
    }

    /// Track parameters vector.
    Parameters
    smoothed() const;

    /// Track parameters covariance matrix.
    Covariance
    smoothedCovariance() const;

    bool
    hasSmoothed() const
    {
      return m_data.ismoothed != IndexData::kInvalid;
    }

    /// Returns the jacobian associated to this track state
    Covariance
    jacobian() const;

    bool
    hasJacobian() const
    {
      return m_data.ijacobian != IndexData::kInvalid;
    }

    Projector
    projector() const;

    bool
    hasProjector() const
    {
      return m_data.iprojector != IndexData::kInvalid;
    }

    EffectiveProjector
    effectiveProjector() const
    {
      return projector().topLeftCorner(m_data.measdim, M);
    }

    bool
    hasUncalibrated() const
    {
      return m_data.iuncalibrated != IndexData::kInvalid;
    }

    /// Uncalibrated measurement in the form of a source link
    const SourceLink&
    uncalibrated() const;

    /// Check if the point has an associated measurement.
    bool
    hasCalibrated() const
    {
      return m_data.icalibrated != IndexData::kInvalid;
    }

    const SourceLink&
    calibratedSourceLink() const;

    /// Full measurement vector. Might contain additional zeroed dimensions.
    Measurement
    calibrated() const;

    /// Full measurement covariance matrix.
    MeasurementCovariance
    calibratedCovariance() const;

    /// Dynamic measurement vector with only the valid dimensions.
    auto
    effectiveCalibrated() const
    {
      return calibrated().head(m_data.measdim);
    }

    /// Dynamic measurement covariance matrix with only the valid dimensions.
    auto
    effectiveCalibratedCovariance() const
    {
      return calibratedCovariance().topLeftCorner(m_data.measdim,
                                                  m_data.measdim);
    }

  private:
    // Private since it can only be created by the trajectory.
    TrackStateProxy(ConstIf<MultiTrajectory<SourceLink>, ReadOnly>& trajectory,
                    size_t istate);

    ConstIf<MultiTrajectory<SourceLink>, ReadOnly>& m_traj;
    size_t    m_istate;
    IndexData m_data;

    friend class Acts::MultiTrajectory<SourceLink>;
  };

  // implement track state visitor concept
  template <typename T, typename TS>
  using call_operator_t = decltype(std::declval<T>()(std::declval<TS&>()));

  template <typename T, typename TS>
  constexpr bool VisitorConcept
      = concept::require<concept::exists<call_operator_t, T, TS>>;

}  // namespace detail_lt

/// Store a trajectory of track states with multiple components.
///
/// This container supports both simple, sequential trajectories as well
/// as combinatorial or multi-component trajectories. Each point can store
/// a parent point such that the trajectory forms a directed, acyclic graph
/// of sub-trajectories. From a set of endpoints, all possible sub-components
/// can be easily identified. Some functionality is provided to simplify
/// iterating over specific sub-components.
/// @tparam source_link_t Type to link back to an original measurement
template <typename source_link_t>
class MultiTrajectory
{
public:
  enum {
    ParametersSize     = NGlobalPars,
    MeasurementSizeMax = NGlobalPars,
  };
  using SourceLink           = source_link_t;
  using ConstTrackStateProxy = detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, true>;
  using TrackStateProxy = detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, false>;

  using ProjectorBitset = std::bitset<ParametersSize * MeasurementSizeMax>;

  /// Create an empty trajectory.
  MultiTrajectory() = default;

  /// Add a point without measurement and return its index.
  ///
  /// @param trackParameters  at the local point
  /// @param iprevious        index of the previous state, SIZE_MAX if first
  template <typename parameters_t>
  size_t
  addTrackState(const TrackState<SourceLink, parameters_t>& ts,
                size_t iprevious = SIZE_MAX);

  /// Access a read-only point on the trajectory by index.
  ConstTrackStateProxy
  getTrackState(size_t istate) const
  {
    return {*this, istate};
  }

  /// Access a writable point on the trajectory by index.
  TrackStateProxy
  getTrackState(size_t istate)
  {
    return {*this, istate};
  }

  /// Visit all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   non-modifying functor to be called with each point
  template <typename F>
  void
  visitBackwards(size_t iendpoint, F&& callable) const;
  /// Apply a function to all previous states starting at a given endpoint.
  ///
  /// @param iendpoint  index of the last state
  /// @param callable   modifying functor to be called with each point
  ///
  /// @warning If the trajectory contains multiple components with common
  ///          points, this can have an impact on the other components.
  template <typename F>
  void
  applyBackwards(size_t iendpoint, F&& callable);

private:
  /// index to map track states to the corresponding
  std::vector<detail_lt::IndexData>                                  m_index;
  typename detail_lt::Types<ParametersSize>::StorageCoefficients     m_params;
  typename detail_lt::Types<ParametersSize>::StorageCovariance       m_cov;
  typename detail_lt::Types<MeasurementSizeMax>::StorageCoefficients m_meas;
  typename detail_lt::Types<MeasurementSizeMax>::StorageCovariance   m_measCov;
  std::vector<SourceLink>      m_sourceLinks;
  std::vector<ProjectorBitset> m_projectors;

  friend class detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, true>;
  friend class detail_lt::
      TrackStateProxy<SourceLink, ParametersSize, MeasurementSizeMax, false>;
};

// implementations

namespace detail_lt {
  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline TrackStateProxy<SL, N, M, ReadOnly>::TrackStateProxy(
      ConstIf<MultiTrajectory<SL>, ReadOnly>& trajectory,
      size_t istate)
    : m_traj(trajectory), m_istate(istate), m_data(trajectory.m_index[istate])
  {
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::parameters() const -> Parameters
  {
    IndexData::IndexType idx;
    if (hasSmoothed()) {
      idx = m_data.ismoothed;
    } else if (hasFiltered()) {
      idx = m_data.ifiltered;
    } else {
      idx = m_data.ipredicted;
    }

    return Parameters(m_traj.m_params.data.col(idx).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::covariance() const -> Covariance
  {
    IndexData::IndexType idx;
    if (hasSmoothed()) {
      idx = m_data.ismoothed;
    } else if (hasFiltered()) {
      idx = m_data.ifiltered;
    } else {
      idx = m_data.ipredicted;
    }
    return Covariance(m_traj.m_cov.data.col(idx).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::predicted() const -> Parameters
  {
    assert(m_data.ipredicted != IndexData::kInvalid);
    return Parameters(m_traj.m_params.col(m_data.ipredicted).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::predictedCovariance() const -> Covariance
  {
    assert(m_data.ipredicted != IndexData::kInvalid);
    return Covariance(m_traj.m_cov.col(m_data.ipredicted).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::filtered() const -> Parameters
  {
    assert(m_data.ifiltered != IndexData::kInvalid);
    return Parameters(m_traj.m_params.col(m_data.ifiltered).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::filteredCovariance() const -> Covariance
  {
    assert(m_data.ifiltered != IndexData::kInvalid);
    return Covariance(m_traj.m_cov.col(m_data.ifiltered).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::smoothed() const -> Parameters
  {
    assert(m_data.ismoothed != IndexData::kInvalid);
    return Parameters(m_traj.m_params.col(m_data.ismoothed).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::smoothedCovariance() const -> Covariance
  {
    assert(m_data.ismoothed != IndexData::kInvalid);
    return Covariance(m_traj.m_cov.col(m_data.ismoothed).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::jacobian() const -> Covariance
  {
    assert(m_data.ijacobian != IndexData::kInvalid);
    return Covariance(m_traj.m_cov.col(m_data.ijacobian).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::projector() const -> Projector
  {
    assert(m_data.iprojector != IndexData::kInvalid);
    return bitsetToMatrix<Projector>(m_traj.m_projectors[m_data.iprojector]);
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::uncalibrated() const -> const SourceLink&
  {
    assert(m_data.iuncalibrated != IndexData::kInvalid);
    return m_traj.m_sourceLinks[m_data.iuncalibrated];
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::calibrated() const -> Measurement
  {
    assert(m_data.icalibrated != IndexData::kInvalid);
    return Measurement(m_traj.m_meas.col(m_data.icalibrated).data());
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::calibratedSourceLink() const
      -> const SourceLink&
  {
    assert(m_data.icalibratedsourcelink != IndexData::kInvalid);
    return m_traj.m_sourceLinks[m_data.icalibratedsourcelink];
  }

  template <typename SL, size_t N, size_t M, bool ReadOnly>
  inline auto
  TrackStateProxy<SL, N, M, ReadOnly>::calibratedCovariance() const
      -> MeasurementCovariance
  {
    assert(m_data.icalibrated != IndexData::kInvalid);
    return MeasurementCovariance(
        m_traj.m_measCov.col(m_data.icalibrated).data());
  }

}  // namespace detail_lt

template <typename SL>
template <typename parameters_t>
inline size_t
MultiTrajectory<SL>::addTrackState(const TrackState<SL, parameters_t>& ts,
                                   size_t iprevious)
{
  using CovMap =
      typename detail_lt::Types<ParametersSize, false>::CovarianceMap;

  detail_lt::IndexData p = {ts.referenceSurface()};
  if (iprevious != SIZE_MAX) {
    p.iprevious = static_cast<uint16_t>(iprevious);
  }

  if (ts.parameter.predicted) {
    const auto& predicted         = *ts.parameter.predicted;
    m_params.addCol()             = predicted.parameters();
    CovMap(m_cov.addCol().data()) = *predicted.covariance();
    p.ipredicted                  = m_params.size() - 1;
  }

  if (ts.parameter.filtered) {
    const auto& filtered          = *ts.parameter.filtered;
    m_params.addCol()             = filtered.parameters();
    CovMap(m_cov.addCol().data()) = *filtered.covariance();
    p.ifiltered                   = m_params.size() - 1;
  }

  if (ts.parameter.smoothed) {
    const auto& smoothed          = *ts.parameter.smoothed;
    m_params.addCol()             = smoothed.parameters();
    CovMap(m_cov.addCol().data()) = *smoothed.covariance();
    p.ismoothed                   = m_params.size() - 1;
  }

  // store jacobian
  if (ts.parameter.jacobian) {
    CovMap(m_cov.addCol().data()) = *ts.parameter.jacobian;
    p.ijacobian                   = m_cov.size() - 1;
  }

  // handle measurements
  if (ts.measurement.uncalibrated) {
    m_sourceLinks.push_back(*ts.measurement.uncalibrated);
    p.iuncalibrated = m_sourceLinks.size() - 1;
  }

  if (ts.measurement.calibrated) {
    auto meas    = m_meas.addCol();
    auto measCov = m_measCov.addCol();
    std::visit(
        [&meas, &measCov, &p, this](const auto& m) {
          using meas_t                         = std::decay_t<decltype(m)>;
          meas.template head<meas_t::size()>() = m.parameters();
          CovMap(measCov.data())
              .template topLeftCorner<meas_t::size(), meas_t::size()>()
              = m.covariance();

          // We can only store the projection if we have a calibrated
          // measurement. Place (possibly asymmetric) projector into
          // full size projector, padded with zeroes.
          // Convert to bitset before setting.
          typename TrackStateProxy::Projector fullProjector;
          fullProjector.setZero();
          fullProjector
              .template topLeftCorner<meas_t::size(), MeasurementSizeMax>()
              = m.projector();

          m_projectors.push_back(matrixToBitset(fullProjector));

          // we also need to set the measurement dimension
          p.measdim = meas_t::size();

          m_sourceLinks.push_back(m.sourceLink());
          p.icalibratedsourcelink = m_sourceLinks.size() - 1;
        },
        *ts.measurement.calibrated);
    p.icalibrated = m_meas.size() - 1;
    p.iprojector  = m_projectors.size() - 1;
  }

  m_index.push_back(std::move(p));
  return m_index.size() - 1;
}

template <typename SL>
template <typename F>
void
MultiTrajectory<SL>::visitBackwards(size_t iendpoint, F&& callable) const
{
  static_assert(detail_lt::VisitorConcept<F, ConstTrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");

  while (true) {
    callable(getTrackState(iendpoint));
    // this point has no parent and ends the trajectory
    if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
      break;
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}

template <typename SL>
template <typename F>
void
MultiTrajectory<SL>::applyBackwards(size_t iendpoint, F&& callable)
{
  static_assert(detail_lt::VisitorConcept<F, TrackStateProxy>,
                "Callable needs to satisfy VisitorConcept");
  while (true) {
    callable(getTrackState(iendpoint));
    // this point has no parent and ends the trajectory
    if (m_index[iendpoint].iprevious == detail_lt::IndexData::kInvalid) {
      break;
    }
    iendpoint = m_index[iendpoint].iprevious;
  }
}

}  // namespace Acts