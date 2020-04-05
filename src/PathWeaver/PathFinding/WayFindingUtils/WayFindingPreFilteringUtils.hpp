#pragma once

/*
 * WayFindingPreFilteringUtils.hpp
 *
 *  Created on: Oct 17, 2019
 *      Author: nicholashathaway
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

#include "PathWeaver/PathFinding/HaploPathFinder.hpp"
#include "PathWeaver/PathFinding/WayFindingUtils/TandemRepeatUtils.hpp"


namespace njhseq {


std::string getPossibleSampleNameFromFnp(const bfs::path & fnp);

std::string getPossibleSampleNameFromSeqName(const std::string & seqName, const std::string & defaultSamp = "samp1");



struct preprocessSeqsForWayFindingPars{
	SeqIOOptions filteredPairedOpts;
	SeqIOOptions filteredSingletOuts;
	SeqIOOptions filteredOff_pairedOpts;
	SeqIOOptions filteredOff_singletOuts;
	SeqIOOptions filteredOffDups_pairedOpts;
	SeqIOOptions filteredOffDups_singletOuts;

};

struct preprocessSeqsForWayFindingRes{
	bfs::path usedFilteredSinglesFnp;
	bfs::path usedFilteredPairedR1Fnp;
	bfs::path usedFilteredPairedR2Fnp;
	uint64_t maxInputSeqLen = 0;

	uint32_t filteredPairsDups_{0};
	uint32_t filteredSinglesDups_{0};
	uint32_t filteredR1Qual_{0};
	uint32_t filteredR2Qual_{0};
	uint32_t filteredSinglesQual_{0};
	uint32_t filteredR1LowEntropy_{0};
	uint32_t filteredR2LowEntropy_{0};
	uint32_t filteredSinglesLowEntropy_{0};


	Json::Value filteredInfo() const{
		Json::Value ret;
		ret["filteredPairsDups_"] = filteredPairsDups_;
		ret["filteredSinglesDups_"] = filteredSinglesDups_;
		ret["filteredR1Qual_"] = filteredR1Qual_;
		ret["filteredR2Qual_"] = filteredR2Qual_;
		ret["filteredSinglesQual_"] = filteredSinglesQual_;
		ret["filteredR1LowEntropy_"] = filteredR1LowEntropy_;
		ret["filteredR2LowEntropy_"] = filteredR2LowEntropy_;
		ret["filteredSinglesLowEntropy_"] = filteredSinglesLowEntropy_;
		return ret;
	}
};

preprocessSeqsForWayFindingRes preprocessSeqsForWayFinding(
		const preprocessSeqsForWayFindingPars & outPars,
		const BamExtractor::ExtractedFilesOpts & inOpts,
		const HaploPathFinder::PathFinderCorePars & extractionPars,
		const std::string & sampName);



struct repairSeqsAfterRecruitmentPars {
	bfs::path singlesFnp_;
	bfs::path r1Fnp_;
	bfs::path r2Fnp_;
};
struct repairSeqsAfterRecruitmentResults {
	uint32_t numberOfPairsRePaired_{0};
};

repairSeqsAfterRecruitmentResults repairSeqsAfterRecruitment(const repairSeqsAfterRecruitmentPars & rePairPars);


struct writeOutTandemsAndOptionallyStitchPars {
	writeOutTandemsAndOptionallyStitchPars(const OutOptions & tandemInfoOpts,
			const bfs::path & workingDir) :
			tandemInfoOpts_(tandemInfoOpts),
			workingDir_(workingDir){

	}

	SeqIOOptions usedFilteredPairedOpts;
	SeqIOOptions usedFilteredSinglesOpts;
	OutOptions tandemInfoOpts_;
	bfs::path workingDir_;

};

struct writeOutTandemsAndOptionallyStitchRes{
	std::unordered_map<std::string, std::vector<TandemRepeatPlusAdjustedSize>> previouslyDeterminedTandemsForSingles_;
	std::unordered_map<std::string, std::vector<TandemRepeatPlusAdjustedSize>> previouslyDeterminedTandemsForPairs_;
	std::vector<TandemRepeatPlusAdjustedSize> allTandems;

};

writeOutTandemsAndOptionallyStitchRes writeOutTandemsAndOptionallyStitch(
		const writeOutTandemsAndOptionallyStitchPars & runPars,
		const uint64_t maxInputSeqsLen,
		const HaploPathFinder::PathFinderCorePars & extractionPars
		);

writeOutTandemsAndOptionallyStitchRes writeOutTandemsAndOptionallyStitch(
		const writeOutTandemsAndOptionallyStitchPars & runPars,
		const uint64_t maxInputSeqsLen,
		const HaploPathFinder::PathFinderCorePars & extractionPars,
		const writeOutTandemsAndOptionallyStitchRes & previouslyDeterminedTandems
		);




}  // namespace njhseq
