/*
 * KmerExpansionHelper.cpp
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




#include "KmerExpansionHelper.hpp"

namespace njhseq {



KmerExpansionHelperDev::CondensedSeq::CondensedSeq(const std::string & seq){
  uint32_t currentCount = 1;
  uint32_t i = 1;
  for (; i < seq.length(); ++i) {
    if (seq[i] == seq[i - 1]) {
      ++currentCount;
    } else {
    	base_.push_back(seq[i - 1]);
    	count_.push_back(currentCount);
      currentCount = 1;
    }
  }
	base_.push_back(seq[i - 1]);
	count_.push_back(currentCount);
}


KmerExpansionHelperDev::KmerExpansionHelperDev(const std::unordered_map<std::string, motif>& tandemMots,
		const std::unordered_map<std::string, std::unordered_map<std::string, motif>> & tandemsAltMots,
		const uint32_t klength) :
		tandemMots_(tandemMots),
		tandemsAltMots_(tandemsAltMots),
		klength_(klength) {
}

KmerExpansionHelperDev::KmerExpansionHelperDev(const std::unordered_map<std::string, motif>& tandemMots,
		const std::unordered_map<std::string, std::unordered_map<std::string, motif>> & tandemsAltMots,
		const uint32_t klength,
		bool keepPositionsSmallerThanKlen) :
		tandemMots_(tandemMots),
		tandemsAltMots_(tandemsAltMots),
		klength_(klength),
		keepPositionsSmallerThanKlen_(keepPositionsSmallerThanKlen){
}


void KmerExpansionHelperDev::addToTandemLocs(
		const std::string & seq,
		const std::unordered_map<std::string, motif>& tandemMots,
		std::unordered_map<std::string, std::unordered_map<std::string, motif>> & tandemsAltMots,
		uint32_t klength,
		std::unordered_map<std::string, std::vector<StartStopTanPos>> & tandemLocations){

	for(const auto & t : tandemMots){
		auto minSize = std::max<uint32_t>(t.second.size(), ((klength > t.second.size() ? klength - t.second.size() : t.second.size())));

		if(tandemsAltMots[t.first].empty()){
			auto locs = t.second.findPositionsFull(seq, tandemMismatchCutOff);
			njh::sort(locs);
			if (!locs.empty()) {
				uint32_t length = 1;
				size_t start = locs.front();
				for (const auto pos : iter::range<uint32_t>(1, locs.size())) {
					if (locs[pos] == locs[pos - 1] + t.second.size()) {
						++length;
					} else {
						//add
						if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
							tandemLocations[t.first].emplace_back(start, start + length * t.second.size());
						}
						length = 1;
						start = locs[pos];
					}
				}
				if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
					tandemLocations[t.first].emplace_back(start, start + length * t.second.size());
				}
			}
		}else{
			std::map<std::string, std::vector<StartStopTanPos>> allAltTandemLocations;
			{
				//main tandem
				auto locs = t.second.findPositionsFull(seq, tandemMismatchCutOff);
				njh::sort(locs);
				if (!locs.empty()) {
					uint32_t length = 1;
					size_t start = locs.front();
					for (const auto pos : iter::range<uint32_t>(1, locs.size())) {
						if (locs[pos] == locs[pos - 1] + t.second.size()) {
							++length;
						} else {
							//add
							if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
								allAltTandemLocations[t.first].emplace_back(start, start + length * t.second.size());
							}
							length = 1;
							start = locs[pos];
						}
					}
					//add
					if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
						allAltTandemLocations[t.first].emplace_back(start, start + length * t.second.size());
					}
				}
			}
			for(const auto & altT : tandemsAltMots[t.first]){
				//alt tandems
				auto locs = altT.second.findPositionsFull(seq, tandemMismatchCutOff);
				njh::sort(locs);
				if (!locs.empty()) {
					uint32_t length = 1;
					size_t start = locs.front();
					for (const auto pos : iter::range<uint32_t>(1, locs.size())) {
						if (locs[pos] == locs[pos - 1] + altT.second.size()) {
							++length;
						} else {
							//add
							if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
								allAltTandemLocations[altT.first].emplace_back(start, start + length * altT.second.size());
							}
							length = 1;
							start = locs[pos];
						}
					}
					//add
					if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
						allAltTandemLocations[altT.first].emplace_back(start, start + length * altT.second.size());
					}
				}
			}
			std::string longestTandem = "";
			uint32_t longestTandemLen = 0;
			for(const auto & allAltTandemLocation : allAltTandemLocations){
				for(const auto & tsection : allAltTandemLocation.second){
					if(tsection.length() > longestTandemLen){
						longestTandem = allAltTandemLocation.first;
						longestTandemLen = tsection.length();
					}
				}
			}
			if("" != longestTandem){
				for(const auto & tsection : allAltTandemLocations[longestTandem]){
					tandemLocations[longestTandem].emplace_back(tsection);
				}
			}
		}
	}
}


void KmerExpansionHelperDev::addToGlobalTandems(const std::string & seq,
		uint32_t klength,
		std::vector<TandemRepeatPlusAdjustedSize> & currentTandems){
	auto firstTandems = aligner::findTandemRepeatsInSequence(seq,  2, -2, -7, klength * 2);

	for (const auto & t : firstTandems) {
		uint32_t adjustedSize = t.getSize();
		if (t.startPos_ > 3 && t.startPos_ < t.repeat_.size()) {
			std::string begSeq = seq.substr(0, t.startPos_);
			std::string tPortionSeq = t.repeat_.substr(
					t.repeat_.size() - t.startPos_);
			if (numberOfMismatches(begSeq, tPortionSeq) <= 0) {
				adjustedSize += t.startPos_;
			}
		}
		if (seq.size() - t.stopPos_ > 3
				&& seq.size() - t.stopPos_ < t.repeat_.size()) {
			std::string backSeq = seq.substr(t.stopPos_);
			std::string tPortionSeq = t.repeat_.substr(0,
					seq.size() - t.stopPos_);
			if (numberOfMismatches(backSeq, tPortionSeq) <= 0) {
				adjustedSize += seq.size() - t.stopPos_;
			}
		}
		currentTandems.emplace_back(t, adjustedSize);
	}
}

std::unordered_map<std::string, std::vector<StartStopTanPos>> KmerExpansionHelperDev::processForTandems(const std::string & seq,
		const std::unordered_map<std::string, motif>& tandemMots,
		std::unordered_map<std::string, std::unordered_map<std::string, motif>> & tandemsAltMots,
		uint32_t klength){
	//single
	std::unordered_map<std::string, std::vector<StartStopTanPos>> tandemLocations;
	{
		for (const auto & t : tandemMots) {
			auto minSize = std::max<uint32_t>(t.second.size(), ((klength > t.second.size() ? klength - t.second.size() : t.second.size())));
			if(tandemsAltMots[t.first].empty()){
				auto locs = t.second.findPositionsFull(seq, tandemMismatchCutOff);
				njh::sort(locs);
				if (!locs.empty()) {
					uint32_t length = 1;
					size_t start = locs.front();
					for (const auto pos : iter::range<uint32_t>(1, locs.size())) {
						if (locs[pos] == locs[pos - 1] + t.second.size()) {
							++length;
						} else {
							//add
							if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
								tandemLocations[t.first].emplace_back(start, start + length * t.second.size());
							}
							length = 1;
							start = locs[pos];
						}
					}
					//add
					if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
						tandemLocations[t.first].emplace_back(start, start + length * t.second.size());
					}
				}
			} else {
				std::map<std::string, std::vector<StartStopTanPos>> allAltTandemLocations;
				{
					//main tandem
					auto locs = t.second.findPositionsFull(seq, tandemMismatchCutOff);
					njh::sort(locs);
					if (!locs.empty()) {
						uint32_t length = 1;
						size_t start = locs.front();
						for (const auto pos : iter::range<uint32_t>(1, locs.size())) {
							if (locs[pos] == locs[pos - 1] + t.second.size()) {
								++length;
							} else {
								//add
								if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
									allAltTandemLocations[t.first].emplace_back(start, start + length * t.second.size());
								}
								length = 1;
								start = locs[pos];
							}
						}
						//add
						if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
							allAltTandemLocations[t.first].emplace_back(start, start + length * t.second.size());
						}
					}
				}
				for(const auto & altT : tandemsAltMots[t.first]){
					//alt tandems
					auto locs = altT.second.findPositionsFull(seq, tandemMismatchCutOff);
					njh::sort(locs);
					if (!locs.empty()) {
						uint32_t length = 1;
						size_t start = locs.front();
						for (const auto pos : iter::range<uint32_t>(1, locs.size())) {
							if (locs[pos] == locs[pos - 1] + altT.second.size()) {
								++length;
							} else {
								//add
								if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
									allAltTandemLocations[altT.first].emplace_back(start, start + length * altT.second.size());
								}
								length = 1;
								start = locs[pos];
							}
						}
						//add
						if(length >= tandemRepeatNumbercutOff && length * t.second.size() >= minSize){
							allAltTandemLocations[altT.first].emplace_back(start, start + length * altT.second.size());
						}
					}
				}
				std::string longestTandem = "";
				uint32_t longestTandemLen = 0;
				for(const auto & allAltTandemLocation : allAltTandemLocations){
					for(const auto & tsection : allAltTandemLocation.second){
						if(tsection.length() > longestTandemLen){
							longestTandem = allAltTandemLocation.first;
							longestTandemLen = tsection.length();
						}
					}
				}
				if("" != longestTandem){
					for(const auto & tsection : allAltTandemLocations[longestTandem]){
						tandemLocations[longestTandem].emplace_back(tsection);
					}
				}
			}
		}
		// this is no longer needed with the now additional post processing
//		std::cout << "tandemLocations.size(): " << tandemLocations.size() << std::endl;
//
//		if(tandemLocations.size() > 1){
//			auto tandemLocationsCopy = tandemLocations;
//			std::unordered_map<std::string, StartStopTanPos> startsStopsTands;
//			for(const auto & tand : tandemLocations){
//				auto currentStart = std::min_element(
//						tand.second.begin(),
//						tand.second.end(),
//						[](const StartStopTanPos & pos1,const StartStopTanPos & pos2) {
//							return pos1.start_ < pos2.start_;
//						})->start_;
//				uint32_t currentStop = std::max_element(
//						tand.second.begin(),
//						tand.second.end(),
//						[](const StartStopTanPos & pos1,const StartStopTanPos & pos2) {
//							return pos1.stop_ < pos2.stop_;
//				})->stop_;
//				startsStopsTands.emplace(tand.first, StartStopTanPos{currentStart, currentStop});
//			}
//			std::set<std::string> withinOther;
//			for(const auto & tand : startsStopsTands){
//				for(const auto & secondTand : startsStopsTands){
//					if(tand.first == secondTand.first){
//						continue;
//					}
//					if(tand.second.start_ >= secondTand.second.start_ &&
//							tand.second.stop_ <= secondTand.second.stop_){
//						withinOther.emplace(tand.first);
//					}
//				}
//			}
//			if(!withinOther.empty()){
//				tandemLocations.clear();
//				for(const auto & tand : tandemLocationsCopy){
//					if(!njh::in(tand.first, withinOther)){
//						tandemLocations.emplace(tand);
//					}
//				}
//			}
//		}
	}
	return tandemLocations;
}

std::vector<KmerPathwayGraphDev::StartEndPos> KmerExpansionHelperDev::processTandemLocsIntoPositions(const std::string & seq,
		 std::unordered_map<std::string, std::vector<StartStopTanPos>> tandemLocations){
	std::vector<KmerPathwayGraphDev::StartEndPos> ret;
	//walk out the tandems to get the partials
	for(auto & tloc : tandemLocations){
		std::string tandemStr = tloc.first;
		for(auto & loc : tloc.second){
			uint32_t startWalkBack = 0;
			while(loc.start_> 0 && startWalkBack < tandemStr.size() && seq[loc.start_ - 1] == tandemStr[tandemStr.size() - 1 - startWalkBack]){
				--(loc.start_);
				++startWalkBack;
			}
			uint32_t stopWalkForward = 0;
			while(loc.stop_  < seq.size() && seq[loc.stop_] == tandemStr[stopWalkForward]){
				++stopWalkForward;
				++(loc.stop_);
			}
		}
	}

	for(auto & tloc : tandemLocations){
		njh::sort(tloc.second, [](const StartStopTanPos & p1, const StartStopTanPos & p2){
			if(p1.start_ == p2.start_){
				return p1.stop_ < p2.stop_;
			}else{
				return p1.start_ < p2.start_;
			}
		});

		std::vector<KmerPathwayGraphDev::StartEndPos> locsForCurrent;
		for(const auto & loc : tloc.second){
			if(locsForCurrent.empty()){
				locsForCurrent.emplace_back(loc.start_, loc.stop_);
			}else{
				if(locsForCurrent.back().end_ + tloc.first.size() >=loc.start_){
					locsForCurrent.back().end_ = loc.stop_;
				}else{
					locsForCurrent.emplace_back(loc.start_, loc.stop_);
				}
			}
		}
		addOtherVec(ret, locsForCurrent);
	}
	return ret;
}

std::vector<KmerPathwayGraphDev::StartEndPos> KmerExpansionHelperDev::getPositionsForHomopolymers(const std::string & seq,
		const uint32_t klength){

	std::vector<KmerPathwayGraphDev::StartEndPos> positions;
	CondensedSeq testC(seq);
	uint32_t seqPosition = 0;
	for(const auto pos : iter::range(testC.base_.size())){
		if(testC.count_[pos] >= klength){
			positions.emplace_back(KmerPathwayGraphDev::StartEndPos(seqPosition, seqPosition + testC.count_[pos]));
		}
		seqPosition += testC.count_[pos];
	}
	return positions;
}

std::vector<KmerPathwayGraphDev::StartEndPos> KmerExpansionHelperDev::getExpansionPositions(const std::string & seq){
	auto tandemLocations = processForTandems(seq, tandemMots_, tandemsAltMots_, klength_);
	auto expansionLocs = processTandemLocsIntoPositions(seq, tandemLocations);
	auto homopolymersExpLocs = getPositionsForHomopolymers(seq, klength_);
	addOtherVec(expansionLocs, homopolymersExpLocs);
	auto mergedPositions = KmerPathwayGraphDev::StartEndPos::merge(expansionLocs);
	if(!keepPositionsSmallerThanKlen_){
		std::vector<KmerPathwayGraphDev::StartEndPos> keep;
		for(const auto & pos : mergedPositions){
			if(pos.len() >= klength_ - 1){
				keep.emplace_back(pos);
			}
		}
		mergedPositions = keep;
	}
	return mergedPositions;
}



uint32_t KmerExpansionHelperDev::tandemMismatchCutOff = 0;
uint32_t KmerExpansionHelperDev::tandemRepeatNumbercutOff = 1;


}  // namespace njhseq
