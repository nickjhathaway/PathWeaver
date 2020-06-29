#pragma once
/*
 * HaploPathFinder.hpp
 *
 *  Created on: May 22, 2017
 *      Author: nick
 */
 // PathWeaver - A library for running local haplotype assembly
 // Copyright (C) 2012-2020 Nicholas Hathaway <nickjhathaway@gmail.com>,
 //
 // This file is part of PathWeaver.
 //
 // PathWeaver is free software: you can redistribute it and/or modify
 // it under the terms of the GNU General Public License as published by
 // the Free Software Foundation, either version 3 of the License, or
 // (at your option) any later version.
 //
 // PathWeaver is distributed in the hope that it will be useful,
 // but WITHOUT ANY WARRANTY; without even the implied warranty of
 // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 // GNU General Public License for more details.
 //
 // You should have received a copy of the GNU General Public License
 // along with PathWeaver.  If not, see <http://www.gnu.org/licenses/>.
 //
 //


#include <njhseq/GenomeUtils.h>
#include <SeekDeep/objects/IlluminaUtils/PairedReadProcessor.hpp>


namespace njhseq {

class HaploPathFinder {
public:
	struct PathFinderCorePars {


		PathFinderCorePars();
		bool needToRePair_ = false;

		//uint32_t connectorCutOff = 1;

		bool trimEnds = false;



		bool trimEndsByCoverage = false;

		bool verbose = false;
		bool debug = false;
		bool exitOnException_ = false;

		bool writeNodeCounts_{false};

		bool graphVerbose_ = false;
		bool graphDebug_ = false;
		bool debugWriteEdgeInfo_ = false;
		bool writeOutFinalDot_ = false;
		bool writeOutDotsByGroups_ = false;
		uint32_t writeOutDotsByGroupsSizeCutOff_{2}; //! needs to be at least this size /

		bool removeDuplicatedSequences_ = false;
		bool trimRemoveDuplicatedSequences_ = false;


		bool trimSeqsEdgesForLowEntropy_ = true;
		readVecTrimmer::TrimEdgesByLowEntropyPars seqEdgeEntropyTrimPars_;


		bool trimSeqsWithNs_ = false;
		uint32_t qualTrim_ = 2;
		bool trimOnQual_ = false;

		double qualCheckCutOff_  = 0.5;
		uint32_t qualCheck_ = 30;
		bool preFilterReadsOnEntropy_{true};
		double preFilterReadsOnEntropyCutOff_{0.75};
		bool markPreFilterInfo_{true};

		double nodeBridgingReadPercCutOff_ = 0.00;
		bool startDisentanglementConservative_ {false};
		bool adjustForHeadLessAddedAfterDisentanglement{false};
		bool throwAwayConservedAddNodesDuringDisentaglement{false};

		bool collapseNodesWithAllowableError_ = false;
		comparison errorsToAllow_;
		double collapseNodesWithAllowableErrorFreqMultiplier_ = 3.0;


		bool collapseOneBaseIndelsNodes_ = false;
		bool collapseOneBaseIndelsBeforeDisentanglement_ = false;

		bool collapsePossibleSimpleLoops_ = false;

		bool trimTipsOfLowEntropyNodes_{false};
		double trimTipsOfLowEntropyNodesCutOff_{1.5};


		bool trimEdgesNodeTipsWithLowEntropy_{false};
		//double trimEdgesNodeTipsWithLowEntropyCutOff{1.25};
		readVecTrimmer::TrimEdgesByLowEntropyPars trimEdgesOfEndNodesPars_;


		double oneBaseIndelError_ = 0.50;
		double twoBaseIndelError_ = 0.50;
		double largeBaseIndelError_ = 0.50;
		double homopolymerIndelCollapseFreqMultiplier_ = 2.0;

		bool initialBreakingSinglePaths_ = false;

		bool splitTailed_ = false;
		bool splitEnds_ = false;
		bool splitEndsOnce_ = false;
		uint32_t maxSplitEndNodeSize_{std::numeric_limits<uint32_t>::max()};
		bool useOptAfterOneEndSplit_ = true;

		bool keepCycles_ = false;

		bool addSeqOfSingleHeadAndTailSeqs_= false;

		bool collapseFinalDeterminedTandems_ = false;
		uint32_t finalTandemOccurenceCutOff_ = 2;
		bool useSmallerThanKlenExpandedPositions_ = false;

		bool addGroups_ = false;
		bool writeGroupInputNames_ = false;

		bool removeHeadlessTaillessAfterDisentaglement = false;
		bool removeHeadlessTaillessAlongTheWay = false;
		uint32_t headlessTailessLenCutOff = 0; //good default is same size of kmer length so should change with klen change

		bool collapseFinalClusters_ = false;
		bool doNotCollapseIdenticalContigs_ = false;

		bool writeAllPossibleHaps_ = false;
		bool writeEstimatedMajorHaps_ = false;
		bool writeOutFinalConnections_ = false;

		uint32_t numThreads = 1;


		uint32_t lenCutOff = 0;
		uint32_t taillessOrHeaddLesslenCutOff = 0;
		uint32_t trimEndsBy = 0;

		bool removePossibleOutliersWithMuscle = false;
		double outliersCutOff = .90;

		bool removePossibleOutliersKSim = false;
		uint32_t outlierKLen = 10;
		double outliersKSimCutOff = 0.40;
		bool keepGroupsForOutlierFilter = false; //keep the whole group if one contig of a group passes the filter
		bool filterOffOutlierInputSeqs = false; //remove recruited sequences if they only fall within the outliter sequences


		bool stitchPairs = false;
		bool reOrientPairsForStitching = false;
	  PairedReadProcessor::ProcessParams pairProcessorParams_;
		std::vector<PairedReadProcessor::ReadPairOverLapStatus> acceptableOverlapStatuses {
				PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2,
				PairedReadProcessor::ReadPairOverLapStatus::PERFECTOVERLAP };


		std::vector<seqInfo> inputSeqs;
		std::vector<seqInfo> trimSeqs;


		bool trimToInputSeqs = false;
		bool trimWithGlobalAln = true;
		bool trimToCircularGenome = true;
		readVecTrimmer::trimCircularGenomeToRefPars circularTrimPars_;

		bool splitToRecruit_ {false};
		bool addHeadTailSeqsToRecruit = false;

		Muscler::TrimWithMusclePars mtPars;

		uint32_t optKLenStart = 35;
		uint32_t optKLenStop = 45;
		uint32_t optKLenStep = 10;
		std::vector<uint32_t> kmerLengths{35, 45};
		std::vector<uint32_t> originalKmerLengths;
		std::vector<uint32_t> initialKmerLengths;



		uint32_t optKcutStart = 2;
		uint32_t optKcutStop = 5;
		uint32_t optKcutStep = 1;
		//std::vector<uint32_t> kmerKOcurrenceCutOffs{2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40}; kmerKOcurrenceCutOffs =
		std::vector<uint32_t> kmerKOcurrenceCutOffs{2,3,4,5};
		std::vector<uint32_t> initialKmerKOcurrenceCutOffs;
		bool autoDetermineKCutsOnNodeCount = false;
		bool autoDetermineKCutsOnTotalBaseCount = false;
		bool autoDetermineKCuts = false;
		bool forceDetermineKCuts = false;
		//std::set<double> autoDetermineKCutsPercentages{0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.20};
		std::set<double> autoDetermineKCutsPercentages{0.02,0.03,0.04,0.05};
		std::set<double> initialAutoDetermineKCutsPercentages;


		bool trimShortTips_ = false;
		bool trimAllShortTips_ = false;
		bool forceTipsByFreq_ = false;
		uint32_t optShortTipNumberStart = 2;
		uint32_t optShortTipNumberStop  = 10;
		uint32_t optShortTipNumberStep  = 4;
		//std::vector<uint32_t> shortTipNumbers{2, 6, 10};
		std::vector<uint32_t> shortTipNumbers{2, 10, 30};
		std::vector<uint32_t> initialShortTipNumbers;


		bool forceFullKcutOptimize_ = false;
		uint32_t kCutOptimizationAttempts = 5;
		bool forceFullShortTipOptimize_ = false;
		uint32_t shortTipOptimizationAttempts = 1;


//		uint32_t optimizeNodeSizeCutOff_ = 150;
		uint32_t optimizeNodeSizeCutOff_ = 100;

		bool useFullOptimalCount_ = false;

		bool optimizeOnIteration = false;
		bool optAfterFirstRecruit = false;

		bool keepOptimizedSubDirs = false;
		//double percentOfInputReadsUsedForOptimization = .80;
		double percentOfInputReadsUsedForOptimization = .95;
		double percentOfInputReadsUsedForOptimizationStepDown = .05;
		double percentageOfKmersUsedCutOff = 0.75;
		bool calcPercentUsedByKmerUsage = false;

		uint32_t estimatorKlen = 15;

		/**@todo update with new members*/
		Json::Value toJson() const;

		void setTandemRepeatHandlingOptions(seqSetUp & setUp);
		void setOptimizationParameters(seqSetUp & setUp);
		void setStitchParameters(seqSetUp & setUp);
		void setOutlierRemovalOptions(seqSetUp & setUp);
		void setQualityTrimAndFiltOpts(seqSetUp & setUp);
		void setFurtherSplittingOpts(seqSetUp & setUp);
		void setHeadlessTaillessNodesOpts(seqSetUp & setUp);
		void setNodeProcessingOpts(seqSetUp & setUp);
		void setDuplicatedSeqHandlingOpts(seqSetUp & setUp);


		void setFinalLengthCutOffs(seqSetUp & setUp);
		void setRunningOpts(seqSetUp & setUp);
		void setPossibleHapsOpts(seqSetUp & setUp);
		void setAddingGroupInfoOpts(seqSetUp & setUp);

		void setPostProcessTrimmingOpts(seqSetUp & setUp);
		void setPostProcessContigHandlingOpts(seqSetUp & setUp);


	};

	struct ExtractParams {
		ExtractParams();
		bfs::path inputDir = "";
		bfs::path extractionDir = "";
		bfs::path regionDir = "";
		bool overWriteDir = false;
		bfs::path metaDataFnp = "";
		bfs::path genomeDir = "";
		std::string primaryGenome = "";
		std::string selectedGenomes = "";
		bfs::path genomeFnp = "";

		bfs::path outputDir = "";

		bfs::path trimSeqBedFnp;

		void setGenomeFnp();

		BamExtractor::extractReadsWtihCrossRegionMappingPars bamExtractPars_;
		BamExtractor::writeUnMappedSeqsAndSmallAlnsWithFiltersPars unmappedPars_;

		GenomicRegion region_;

		Json::Value toJson() const;

		PathFinderCorePars pFinderPars_;

		void setGenomeOpts(seqSetUp & setUp, bool genomeRequired);


		void setBamExtractOpts(seqSetUp & setUp);



	};

	Json::Value toJson() const;

	HaploPathFinder(ExtractParams inPars);

	ExtractParams pars_;


	std::unique_ptr<MultiGenomeMapper> gMapper_;


};

}  // namespace njhseq





