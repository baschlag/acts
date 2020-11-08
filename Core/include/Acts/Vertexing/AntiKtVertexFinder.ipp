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

  auto metricDistance = [&](const PseudoTrack& a, const PseudoTrack& b){
  	double minPtInv = std::min(1./a.sumPt, 1./b.sumPt);
  	double zdiff = std::abs(a.zMean - b.zMean);
  	// Note: according to gaus it's sqrt(sum) but sqrt not needed
  	double combinedError = a.sumSigZ0squared + b.sumSigZ0squared;
  	double cutoff = m_cfg.radiusParameter;
  	return  minPtInv * zdiff * combinedError / cutoff;
  };
  std::vector<PseudoTrack> pseudoTracks;
  for(const auto& trk : trackVector){
  	const BoundTrackParameters& params = m_extractParameters(*trk);
  	PseudoTrack pseudoTrack;
  	pseudoTrack.sumPt = 1./std::abs(params.parameters()[4]);
  	pseudoTrack.zValues.push_back(params.parameters()[1]);
  	pseudoTrack.zMean = params.parameters()[1];
  	const auto& covariance = *params.covariance();
  	double sigZ0 = covariance(Acts::eBoundLoc1, Acts::eBoundLoc1);
  	pseudoTrack.sumSigZ0squared = sigZ0 * sigZ0;
  	pseudoTrack.tracks.push_back(trk);
  	pseudoTracks.push_back(pseudoTrack);
  }


  int nValidTracks = pseudoTracks.size();

  while(nValidTracks != 0){
  	 std::cout << "nValidTracks: " << nValidTracks << std::endl;
  	  // antikt diB calculation
  	  double minSumPtInv = 1e5;
  	  int minSumPtInvTrackIdx = 0;
  	  for(unsigned int trkCount = 0; trkCount < pseudoTracks.size(); trkCount++){
	  	const PseudoTrack& pseudoTrack = pseudoTracks[trkCount];
  	  	double sumPtInv = 1./pseudoTrack.sumPt;
  	  	if(sumPtInv < minSumPtInv){
  	  		minSumPtInv = sumPtInv;
  	  		minSumPtInvTrackIdx = trkCount;
  	  	}
  	  }

  	  // antikt dij calculation
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
		  		minDistance = distance;
		  		std::cout << "\t new min: " << std::endl;
		  	}
		  	std::cout << trkCount1 << "," << trkCount2 << ": " << distance << std::endl;
	  	} // trk 2 loop
	  } // trk 1 loop

	  // 
	  if(minSumPtInv < minDistance){
	  	std::cout << "CREATING VERTEX NOW!!" << std::endl;
	  	PseudoTrack& currentTrack = pseudoTracks[minSumPtInvTrackIdx];
	  	currentTrack.isVertex = true;
	  	currentTrack.isValid = false;
	  	nValidTracks--;
	  }
	  else{
	  	std::cout << "################### updating now... selected tracks:" << std::endl;
		  std::cout << "Track 1: " << minTrk1 << ", z/pt: " << pseudoTracks[minTrk1].zMean << "/" << pseudoTracks[minTrk1].sumPt << std::endl;
		  std::cout << "Track 2: " << minTrk2 << ", z/pt: " << pseudoTracks[minTrk2].zMean << "/" << pseudoTracks[minTrk2].sumPt << std::endl;

		  // Merge trk1 and trk2 as one combined pseudotrack
		  PseudoTrack& updatedTrack = pseudoTracks[minTrk1];
		  updatedTrack.tracks.push_back(trackVector[minTrk2]);
		  updatedTrack.sumPt += pseudoTracks[minTrk2].sumPt;
		  updatedTrack.sumSigZ0squared += pseudoTracks[minTrk2].sumSigZ0squared;

		  for(auto entry : pseudoTracks[minTrk2].zValues){
		  	updatedTrack.zValues.push_back(entry);
		  }
		  updatedTrack.zMean = std::accumulate(updatedTrack.zValues.begin(), updatedTrack.zValues.end(), 0.0) / updatedTrack.zValues.size();
		  // Do not use old pseudotrack anymore
		  pseudoTracks[minTrk2].isValid = false;
		  nValidTracks--;

		  std::cout << "after updates::" << std::endl;
		  std::cout << "Track 1: " << minTrk1 << ", z/pt: " << pseudoTracks[minTrk1].zMean << "/" << pseudoTracks[minTrk1].sumPt << std::endl;
		  std::cout << "Track 2: " << minTrk2 << ", z/pt: " << pseudoTracks[minTrk2].zMean << "/" << pseudoTracks[minTrk2].sumPt << std::endl;
		  std::cout << "ntracks in pseudotrack 1: " << pseudoTracks[minTrk1].tracks.size() << std::endl;
		  std::cout << "is still valid track1/2: " << pseudoTracks[minTrk1].isValid << "/" << pseudoTracks[minTrk2].isValid << std::endl;
		  
	  }

	  
	} 

	for(const auto& trk : pseudoTracks){
		if(trk.isVertex){
			std::cout << "Vertex found at pos: " << trk.zMean << " with ntracks: " << trk.tracks.size() << std::endl;
		}
	}


  return vertexCollection;
}


