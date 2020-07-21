#pragma once

/*
 * WayFinderingUtils.hpp
 *
 *  Created on: Jan 30, 2019
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



namespace njhseq {



struct OptimizationReconResult {

	struct Params {
		Params(const uint32_t klen, const uint32_t kcut,
				const uint32_t shortNumber);
		uint32_t klen_;
		uint32_t kcut_;
		uint32_t shortNumber_;

		Json::Value toJson() const;

	};

	struct Dirs {
		Dirs(const bfs::path & klenDir, const bfs::path & kcutDir,
				const bfs::path & shortTipDir);
		bfs::path klenDir_;
		bfs::path kcutDir_;
		bfs::path shortTipDir_;

		Json::Value toJson() const;

	};

	OptimizationReconResult(const Params runParams, const Dirs & runDirs,uint32_t kcutIterNumber, uint32_t shortTipIterNumber);

	Params runParams_;
	Dirs runDirs_;

	uint32_t readLengthAverage_ { 0 };
	uint32_t readLengthMedian_ { 0 };

	uint32_t kcutIterNumber_ { 0 };
	uint32_t shortTipIterNumber_ { 0 };

	uint32_t headlessCount_ { 0 };
	uint32_t taillessCount_ { 0 };
	uint32_t optimalCount_ { 0 };

	uint32_t headlessCountBelowLen_ { 0 };
	uint32_t taillessCountBelowLen_ { 0 };
	uint32_t optimalCountBelowLen_ { 0 };

	uint32_t internalNodesCountBelowLen_ { 0 };

	double percentOfInputUsed_ { 0.0 };
	uint32_t totalInputReads_ { 0 };

	uint32_t numberOfFinalFilteredSeqs_ { 0 };

	Json::Value toJson() const;


	//captured output instead of written
	std::vector<seqInfo> finalFilteredOutSeqs_;
	std::unordered_set<std::string> keepSeqNames_;
	std::shared_ptr<std::stringstream> finalDot_;

	//


	static std::function<
			bool(const OptimizationReconResult&, const OptimizationReconResult&)> sortFunc;

	static std::vector<OptimizationReconResult> getBestResults(
			const std::vector<OptimizationReconResult> & allResults,
			double percentageUsed,
			double percentagedUsedStepDown, bool useFullOptimalCount);

};

}  // namespace njhseq
