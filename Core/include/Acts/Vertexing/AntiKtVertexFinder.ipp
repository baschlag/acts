// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename vfitter_t, typename sfinder_t>
auto Acts::AntiKtVertexFinder<vfitter_t, sfinder_t>::find(
    const std::vector<const InputTrack_t*>& trackVector,
    const VertexingOptions<InputTrack_t>& vertexingOptions) const
    -> Result<std::vector<Vertex<InputTrack_t>>> {

  // List of vertices to be filled below
  std::vector<Vertex<InputTrack_t>> vertexCollection;

  auto metricDistance = [](const PseudoTrack& a, const PseudoTrack& b){
  	return std::min(1./a.sumPt, 1./b.sumPt) * std::abs(a.zMean - b.zMean);
  };

  std::vector<PseudoTrack> pseudoTracks;
  for(const auto& trk : trackVector){
  	const BoundTrackParameters& params = m_extractParameters(*trk);
  	PseudoTrack pseudoTrack;
  	pseudoTrack.sumPt = 1./std::abs(params.parameters()[4]);
  	pseudoTrack.zValues.push_back(params.parameters()[1]);
  	pseudoTrack.zMean = params.parameters()[1];
  	pseudoTrack.tracks.push_back(trk);
  }

  double minDistance = 1e5;
  int minTrk1 = 0;
  int minTrk2 = 0;
  for(unsigned int trkCount1 = 0; trkCount1 < pseudoTracks.size(); trkCount1++){
  	const PseudoTrack& params1 = pseudoTracks[trkCount1];
  	if(not params1.isValid){
  		continue;
  	}
  	double z01 = params1.zMean;
  	for(unsigned int trkCount2 = trkCount1 + 1; trkCount2 < pseudoTracks.size(); trkCount2++){
  		const PseudoTrack& params2 = pseudoTracks[trkCount2];
  		if(not params2.isValid){
  		  continue;
  		}
  		double z02 = params2.zMean;
  	
  		if(std::abs(z01 - z02) > m_cfg.maxZsearchWindow){
	  		continue;
	  	}
	  	double distance = metricDistance(params1, params2);
	  	
	  	if(distance < minDistance){
	  		minTrk1 = trkCount1;
	  		minTrk2 = trkCount2;
	  	}
  	} // trk 2 loop
  } // trk 1 loop

  // Merge trk1 and trk2 as one combined pseudotrack
  PseudoTrack& updatedTrack = pseudoTracks[minTrk1];
  updatedTrack.tracks.push_back(trackVector[minTrk2]);
  updatedTrack.sumPt += pseudoTracks[minTrk2].sumPt;

  for(auto entry : pseudoTracks[minTrk2].zValues){
  	updatedTrack.zValues.push_back(entry);
  }
  updatedTrack.zMean = std::accumulate(updatedTrack.zValues.begin(), updatedTrack.zValues.end(), 0.0) / updatedTrack.zValues.size();
  // Do not use old pseudotrack anymore
  pseudoTracks[minTrk2].isValid = false;
  
  return vertexCollection;
}


