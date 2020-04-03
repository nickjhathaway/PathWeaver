#pragma once

/*
 * KmerExpansionHelper.hpp
 *
 *  Created on: Oct 18, 2019
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


#include "PathWeaver/PathFinding/WayFindingUtils/TandemRepeatUtils.hpp"

#include "PathWeaver/objects/dataContainers/graphs/KmerPathwayGraph.hpp"

namespace njhseq {




class KmerExpansionHelper {
public:


	class CondensedSeq {
	public:
		CondensedSeq(const std::string & seq);

	  std::vector<char> base_;
	  std::vector<uint32_t> count_;

	};

	KmerExpansionHelper(const std::unordered_map<std::string, motif>& tandemMots,
			const std::unordered_map<std::string, std::unordered_map<std::string, motif>> & tandemsAltMots,
			const uint32_t klength);

	KmerExpansionHelper(const std::unordered_map<std::string, motif>& tandemMots,
			const std::unordered_map<std::string, std::unordered_map<std::string, motif>> & tandemsAltMots,
			const uint32_t klength,
			bool keepPositionsSmallerThanKlen);

	std::unordered_map<std::string, motif> tandemMots_;
	std::unordered_map<std::string, std::unordered_map<std::string, motif>> tandemsAltMots_;
	const uint32_t klength_;
	bool keepPositionsSmallerThanKlen_{false};

	static uint32_t tandemMismatchCutOff; // as it stands with the current method this cannot be not 0, especially for say dinucleotide repeats, even allowing 1 wouldn't work
	static uint32_t tandemRepeatNumbercutOff;
	static void addToTandemLocs(const std::string & seq,
			const std::unordered_map<std::string, motif>& tandemMots,
			std::unordered_map<std::string, std::unordered_map<std::string, motif>> & tandemsAltMots,
			uint32_t klength,
			std::unordered_map<std::string, std::vector<StartStopTanPos>> & tandemLocations);

	static void addToGlobalTandems(const std::string & seq,
			uint32_t klength,
			std::vector<TandemRepeatPlusAdjustedSize> & currentTandems);

	static std::unordered_map<std::string, std::vector<StartStopTanPos>> processForTandems(
			const std::string & seq,
			const std::unordered_map<std::string, motif>& tandemMots,
			std::unordered_map<std::string, std::unordered_map<std::string, motif>> & tandemsAltMots,
			uint32_t klength);

	static std::vector<KmerPathwayGraph::StartEndPos> processTandemLocsIntoPositions(
			const std::string & seq,
			std::unordered_map<std::string, std::vector<StartStopTanPos>> tandemLocations);

	static std::vector<KmerPathwayGraph::StartEndPos> getPositionsForHomopolymers(
			const std::string & seq, const uint32_t klength);

	std::vector<KmerPathwayGraph::StartEndPos> getExpansionPositions(
			const std::string & seq);

};





}  // namespace njhseq
