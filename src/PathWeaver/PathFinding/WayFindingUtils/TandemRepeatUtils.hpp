#pragma once

/*
 * TandemRepeatUtils.hpp
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

#include <njhseq/objects/helperObjects/motif.hpp>
#include <njhseq/objects/helperObjects/tandemRepeat.hpp>

namespace njhseq {

struct StartStopTanPos {
	StartStopTanPos(uint32_t start, uint32_t stop);
	uint32_t start_ = std::numeric_limits<uint32_t>::max();
	uint32_t stop_ = std::numeric_limits<uint32_t>::max();
	uint32_t length() const;

};



struct TandemRepeatPlusAdjustedSize {

	TandemRepeatPlusAdjustedSize(const TandemRepeat & repeat, uint32_t adjustedPartialSize);
	TandemRepeat repeat_;
	uint32_t adjustedPartialSize_;
};


std::vector<TandemRepeatPlusAdjustedSize> readInTandemsWithPartInfo(const bfs::path & tandemInfoFnp);




struct ProcessedDeterminedTandemRepeatsResults{

	VecStr finalTandemSeqs;
	std::unordered_map<std::string, uint32_t> maxTandemSizesForEachTandem;
	std::unordered_map<std::string, std::shared_ptr<TandemRepeatPlusAdjustedSize>> maxTandemEachTandem;
	std::unordered_map<std::string, std::set<std::string>> tandemsAlts;
	std::unordered_map<std::string, uint32_t> tandemsCounts;
	std::unordered_map<std::string, uint32_t> tandemsCountsAdjusted;
	std::unordered_map<std::string, uint32_t> tandemsCountsAdjustedAfterCollapse;



	std::unordered_map<std::string, motif> tandemMots;
	std::unordered_map<std::string, std::unordered_map<std::string, motif>> tandemsAltMots;

	void writeOutInfos(const bfs::path & currentKmerDirectory, bool overWrite = true);
};


ProcessedDeterminedTandemRepeatsResults processGlobalTandems(
		const std::vector<TandemRepeatPlusAdjustedSize> & allTandems,
		uint32_t currentKLen,
		const HaploPathFinder::PathFinderCorePars & extractionPars);


}  // namespace njhseq
