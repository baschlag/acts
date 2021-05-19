#!/bin/bash

#Run pythia gen + fatras + truth tracking + vertexing

# do several round with small event numbers each; workaround to mitigate segfault problem with high statistics

# loop for different pileup
 for i in $(seq 1 15); 

 	do 

	pu=$(($i * 10))

	echo Running chain on PU $pu now; 

	# ttbar event generation
	# ./build/bin/ActsExamplePythia8  --events=100 --rnd-seed=$i  --output-dir=build/data/gen/ttbar_mu$pu --output-csv --output-root --gen-cms-energy-gev=14000 --gen-hard-process=Top:qqbar2ttbar=on --gen-npileup=$pu

	# #simulation
	# ./build/bin/ActsExampleFatrasGeneric --input-dir=build/data/gen/ttbar_mu$pu   --output-dir=build/data/sim_generic/ttbar_mu$pu  --output-csv  --select-eta=-2.5:2.5 --select-pt-gev=0.4: --fatras-pmin-gev 0.4 --remove-neutral  --bf-constant-tesla=0:0:2

	# # truth tracking
	# ./build/bin/ActsExampleTruthTracksGeneric --input-dir=build/data/sim_generic/ttbar_mu$pu --output-dir=build/data/reco_generic/ttbar_mu$pu --bf-constant-tesla=0:0:2 --digi-config-file Examples/Algorithms/Digitization/share/default-smearing-config-generic.json

	# cp build/data/gen/ttbar_mu$pu/particles.root build/data/reco_generic/ttbar_mu$pu

	# vertexing + performance
	./build/bin/ActsExampleVertexFinderTrackReaderPerformanceWriter  --bf-constant-tesla=0:0:2 --input-dir=build/data/reco_generic/ttbar_mu$pu --output-dir=build/data/vertexing/ttbar_mu$pu --select-eta=-2.5:2.5 --select-pt-gev=0.4:  --remove-neutral


done

cd build/data/vertexing/
hadd merged$r.root ttbar_mu*/*.root
mv *.root ..
# cd ..
# rm -rf gen  reco_generic  sim_generic  vertexing
# cd ../..




