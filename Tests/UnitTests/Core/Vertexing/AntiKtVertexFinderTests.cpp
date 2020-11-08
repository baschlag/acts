// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Tests/CommonHelpers/DataDirectory.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/FsmwMode1dFinder.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/AntiKtVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"

#include "VertexingDataHelper.hpp"

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Propagator = Propagator<EigenStepper<ConstantBField>>;
using Linearizer = HelicalTrackLinearizer<Propagator>;

// Create a test context
GeometryContext geoContext = GeometryContext();
MagneticFieldContext magFieldContext = MagneticFieldContext();

// Vertex x/y position distribution
std::uniform_real_distribution<> vXYDist(-0.1_mm, 0.1_mm);
// Vertex z position distribution
std::uniform_real_distribution<> vZDist(-20_mm, 20_mm);
// Track d0 distribution
std::uniform_real_distribution<> d0Dist(-0.01_mm, 0.01_mm);
// Track z0 distribution
std::uniform_real_distribution<> z0Dist(-0.2_mm, 0.2_mm);
// Track pT distribution
std::uniform_real_distribution<> pTDist(0.4_GeV, 10_GeV);
// Track phi distribution
std::uniform_real_distribution<> phiDist(-M_PI, M_PI);
// Track theta distribution
std::uniform_real_distribution<> thetaDist(1.0, M_PI - 1.0);
// Track charge helper distribution
std::uniform_real_distribution<> qDist(-1, 1);
// Track IP resolution distribution
std::uniform_real_distribution<> resIPDist(0., 100_um);
// Track angular distribution
std::uniform_real_distribution<> resAngDist(0., 0.1);
// Track q/p resolution distribution
std::uniform_real_distribution<> resQoPDist(-0.01, 0.01);
// Number of vertices per test event distribution
std::uniform_int_distribution<> nVertexDist(1, 6);
// Number of tracks per vertex distribution
std::uniform_int_distribution<> nTracksDist(5, 15);

///
/// @brief
///
BOOST_AUTO_TEST_CASE(antikt_finder_test) {

  // Set up RNG
  int mySeed = 31415;
  std::mt19937 gen(mySeed);


    using BilloirFitter =
        FullBilloirVertexFitter<BoundTrackParameters, Linearizer>;
    using ZScanSeedFinder = ZScanVertexFinder<BilloirFitter>;
    // Vertex Finder
    using VertexFinder = AntiKtVertexFinder<BilloirFitter, ZScanSeedFinder>;

    VertexFinder::Config cfg;
    VertexFinder finder(cfg);

    // Vector to be filled with all tracks in current event
    std::vector<std::unique_ptr<const BoundTrackParameters>> tracks;

    std::vector<const BoundTrackParameters*> tracksPtr;

    // Vector to be filled with truth vertices for later comparison
    std::vector<Vertex<BoundTrackParameters>> trueVertices;

    // start creating event with nVertices vertices
    unsigned int nVertices = nVertexDist(gen);
    for (unsigned int iVertex = 0; iVertex < nVertices; ++iVertex) {
      // Number of tracks
      unsigned int nTracks = nTracksDist(gen);

      // Create perigee surface
      std::shared_ptr<PerigeeSurface> perigeeSurface =
          Surface::makeShared<PerigeeSurface>(Vector3D(0., 0., 0.));

      // Create position of vertex and perigee surface
      double x = vXYDist(gen);
      double y = vXYDist(gen);
      double z = vZDist(gen);

      // True vertex
      Vertex<BoundTrackParameters> trueV(Vector3D(x, y, z));
      std::vector<TrackAtVertex<BoundTrackParameters>> tracksAtTrueVtx;

      // Calculate d0 and z0 corresponding to vertex position
      double d0_v = sqrt(x * x + y * y);
      double z0_v = z;

      // Construct random track emerging from vicinity of vertex position
      // Vector to store track objects used for vertex fit
      for (unsigned int iTrack = 0; iTrack < nTracks; iTrack++) {
        // Construct positive or negative charge randomly
        double q = qDist(gen) < 0 ? -1. : 1.;

        // Construct random track parameters
        BoundVector paramVec;
        double z0track = z0_v + z0Dist(gen);
        paramVec << d0_v + d0Dist(gen), z0track, phiDist(gen), thetaDist(gen),
            q / pTDist(gen), 0.;

        // Resolutions
        double res_d0 = resIPDist(gen);
        double res_z0 = resIPDist(gen);
        double res_ph = resAngDist(gen);
        double res_th = resAngDist(gen);
        double res_qp = resQoPDist(gen);

        // Fill vector of track objects with simple covariance matrix
        Covariance covMat;
        covMat << res_d0 * res_d0, 0., 0., 0., 0., 0., 0., res_z0 * res_z0, 0.,
            0., 0., 0., 0., 0., res_ph * res_ph, 0., 0., 0., 0., 0., 0.,
            res_th * res_th, 0., 0., 0., 0., 0., 0., res_qp * res_qp, 0., 0.,
            0., 0., 0., 0., 1.;
        auto params =
            BoundTrackParameters(perigeeSurface, paramVec, std::move(covMat));

        tracks.push_back(std::make_unique<BoundTrackParameters>(params));

        TrackAtVertex<BoundTrackParameters> trAtVt(0., params,
                                                   tracks.back().get());
        tracksAtTrueVtx.push_back(trAtVt);
      }

      trueV.setTracksAtVertex(tracksAtTrueVtx);
      trueVertices.push_back(trueV);

    }  // end loop over vertices

    // shuffle list of tracks
    std::shuffle(std::begin(tracks), std::end(tracks), gen);

    for (const auto& trk : tracks) {
      tracksPtr.push_back(trk.get());
    }

    VertexingOptions<BoundTrackParameters> vertexingOptions(geoContext,
                                                            magFieldContext);

    // find vertices
    auto res = finder.find(tracksPtr, vertexingOptions);

    BOOST_CHECK(res.ok());

    if (!res.ok()) {
      std::cout << res.error().message() << std::endl;
    }

    // Retrieve vertices found by vertex finder
    auto vertexCollection = *res;

    // if (debug) {
    //   std::cout << "########## RESULT: ########## Event " << iEvent
    //             << std::endl;
    //   std::cout << "Number of true vertices: " << nVertices << std::endl;
    //   std::cout << "Number of reco vertices: " << vertexCollection.size()
    //             << std::endl;

    //   int count = 1;
    //   std::cout << "----- True vertices -----" << std::endl;
    //   for (const auto& vertex : trueVertices) {
    //     Vector3D pos = vertex.position();
    //     std::cout << count << ". True Vertex:\t Position:"
    //               << "(" << pos[eX] << "," << pos[eY] << "," << pos[eZ] << ")"
    //               << std::endl;
    //     std::cout << "Number of tracks: " << vertex.tracks().size() << std::endl
    //               << std::endl;
    //     count++;
    //   }
    //   std::cout << "----- Reco vertices -----" << std::endl;
    //   count = 1;
    //   for (const auto& vertex : vertexCollection) {
    //     Vector3D pos = vertex.position();
    //     std::cout << count << ". Reco Vertex:\t Position:"
    //               << "(" << pos[eX] << "," << pos[eY] << "," << pos[eZ] << ")"
    //               << std::endl;
    //     std::cout << "Number of tracks: " << vertex.tracks().size() << std::endl
    //               << std::endl;
    //     count++;
    //   }
    // }

}

}  // namespace Test
}  // namespace Acts
