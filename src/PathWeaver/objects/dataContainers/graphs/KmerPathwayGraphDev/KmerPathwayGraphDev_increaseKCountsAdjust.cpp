/*
 * KmerPathwayGraphDev_increaseKCountsAdjust.cpp
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



#include "KmerPathwayGraphDev.hpp"


namespace njhseq {

void KmerPathwayGraphDev::increaseKCountsAdjust(const std::string & seq,
		const std::vector<StartEndPos> & processedPositions){
	if(0 == processedPositions.size()){
		increaseKCounts(seq);
	} else	if(1 == processedPositions.size()){
		bool startPass = false;
		bool stopPass = false;
		uint32_t start = processedPositions.front().start_;
		uint32_t stop = processedPositions.front().end_;
		if (start > 1) {
			if (start < seq.size() + 1 - klen_) {
				startPass = true;
			}
			for (auto pos : iter::range<uint32_t>(0, std::min<uint32_t>(start - 1, seq.size() + 1 - klen_))) {
				kCounts_[seq.substr(pos, klen_)] += 1;
				if(seq.substr(pos, klen_).size() < 35){
					//std::cout << seq.substr(pos, klen_) << std::endl;
				}
			}
		}
		if(seq.size() - stop > 1){
			if (stop > klen_) {
				stopPass = true;
			}
			uint32_t rangeStart = stop + 2 < klen_ ? 0 : stop + 2 - klen_;;
			uint32_t rangeStop = seq.size() - klen_ + 1 ;

			for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
				kCounts_[seq.substr(pos, klen_)] += 1;
				if(seq.substr(pos, klen_).size() < 35){
					//std::cout << seq.substr(pos, klen_) << std::endl;
				}
			}
		}
		if(stopPass && startPass){
			//add the kmers that then span the region
			uint32_t newKmerSize = stop - start + 2;
			if(newKmerSize >klen_){
				kCounts_[seq.substr(start - 1, newKmerSize)] += 1;
				if(seq.substr(start - 1, newKmerSize).size() < 35){
					//std::cout << seq.substr(start - 1, newKmerSize) << std::endl;
				}
			}
		}
	} else if(processedPositions.size() > 1){
		//front
		{
			bool startPass = false;
			uint32_t start = processedPositions.front().start_;
			uint32_t stop = processedPositions.front().end_;
			if (start > 1) {
				if (start < seq.size() + 1 - klen_) {
					startPass = true;
				}
				for (auto pos : iter::range<uint32_t>(0, std::min<uint32_t>(start - 1, seq.size() + 1 - klen_))) {
					kCounts_[seq.substr(pos, klen_)] += 1;
					if(seq.substr(pos, klen_).size() < 35){
						//std::cout << seq.substr(pos, klen_) << std::endl;
					}
				}
			}
			if(startPass){
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize >klen_){
					kCounts_[seq.substr(start - 1, newKmerSize)] += 1;

					if(seq.substr(start - 1, newKmerSize).size() < 35){
						//std::cout << seq.substr(start - 1, newKmerSize) << std::endl;
					}
				}
			}
		}

		//connect the positions
		for(const auto pos : iter::range<uint32_t>(processedPositions.size() - 1)){
			if(pos > 0){
				uint32_t start = processedPositions[pos].start_;
				uint32_t stop = processedPositions[pos].end_;
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize > klen_){
					kCounts_[seq.substr(start - 1, newKmerSize)] += 1;
					if(seq.substr(start - 1, newKmerSize).size() < 35){
						//std::cout << seq.substr(start - 1, newKmerSize) << std::endl;
					}
				}
			}
			uint32_t stop = processedPositions[pos].end_;
			uint32_t nextStart = processedPositions[pos + 1].start_;
			uint32_t rangeStart = stop + 2 < klen_ ? stop + 2 : stop + 2 - klen_;
			uint32_t rangeStop = std::min<uint32_t>(nextStart - 1 , seq.size() - klen_ + 1);

			for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
				kCounts_[seq.substr(pos, klen_)] += 1;
				if(seq.substr(pos, klen_).size() < 35){
					//for(const auto & pos : processedPositions){
						//std::cout << pos.start_ << ":" << pos.end_ << ": " << seq.substr(pos.start_, pos.end_ - pos.start_) << std::endl;
					//}
					//std::cout << seq << std::endl;
					//std::cout << seq.substr(pos, klen_) << std::endl;
					//std::cout << "stop: " << stop << std::endl;
					//std::cout << "nextStart: " << nextStart << std::endl;
					//std::cout << "rangeStart: " << rangeStart << std::endl;
					//std::cout << "rangeStop: " << rangeStop << std::endl;
					//std::cout << "seq.size(): " << seq.size() << std::endl;
					//exit(1);
				}
			}
		}
		//end
		{
			bool stopPass = false;
			uint32_t start = processedPositions.back().start_;
			uint32_t stop = processedPositions.back().end_;
			if(seq.size() - stop > 1){
				if (stop > klen_) {
					stopPass = true;
				}
				uint32_t rangeStart = stop + 2 < klen_ ? 0 : stop + 2 - klen_;;
				uint32_t rangeStop = seq.size() - klen_ + 1 ;
				for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
					kCounts_[seq.substr(pos, klen_)] += 1;
					if(seq.substr(pos, klen_).size() < 35){
//					//std::cout << seq.substr(pos, klen_) << std::endl;
					}
				}
			}
			if(stopPass){
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize > klen_){
					kCounts_[seq.substr(start - 1, newKmerSize)] += 1;
					if(seq.substr(start - 1, newKmerSize).size() < 35){
//						//std::cout << seq.substr(start - 1, newKmerSize) << std::endl;
					}
				}
			}
		}
	}
}

void KmerPathwayGraphDev::threadThroughSequenceAdjust(const PairedRead & pseq,
		const std::vector<StartEndPos> & adjustPositionsFirstMate,
		const std::vector<StartEndPos> & adjustPositionsSecondMate){
	if (pseq.seqBase_.seq_.size() < klen_) {
		threadThroughSequenceAdjust(pseq.mateSeqBase_, adjustPositionsSecondMate,
				pseq.seqBase_.name_);
	} else {
		uint32_t threadingSeqNameFirstMateIdx = readNames_.size();
		std::string threadingSeqNameFirstMate = pseq.seqBase_.name_;
		auto firstCounts = threadThroughSequenceAdjust(pseq.seqBase_,
				adjustPositionsFirstMate, threadingSeqNameFirstMate);

		threadThroughSequenceAdjustMate(pseq, adjustPositionsSecondMate,
				firstCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
	}

}

std::unordered_map<std::string, uint32_t> KmerPathwayGraphDev::threadThroughSequenceAdjustMate(const PairedRead & pseq,
		const std::vector<StartEndPos> & processedPositions,
		std::unordered_map<std::string, uint32_t> & firstMateCounts,
		const std::string & threadingSeqNameFirstMate,
		uint32_t threadingSeqNameFirstMateIdx){

	if(processedPositions.size() == 0){
		return threadThroughSequenceMate(pseq, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
	}
	uint32_t threadingSeqNameSecondMateIdx = readNames_.size();
	std::string threadingSeqNameSecondMate = pseq.mateSeqBase_.name_  + "_mate";
	readNames_.emplace_back(threadingSeqNameSecondMate);
	std::unordered_map<std::string, uint32_t> internalCount;

	if(processedPositions.size() == 1){
		bool startPass = false;
		bool stopPass = false;
		uint32_t start = processedPositions.front().start_;
		uint32_t stop = processedPositions.front().end_;
		if (start > 1) {
			if (start < pseq.mateSeqBase_.seq_.size() + 1 - klen_) {
				startPass = true;
			}
			for (auto pos : iter::range<uint32_t>(0, std::min<uint32_t>(start - 2, pseq.mateSeqBase_.seq_.size() + 1 - klen_))) {
				std::string firstKmer = pseq.mateSeqBase_.seq_.substr(pos, klen_);
				std::string nextKmer = pseq.mateSeqBase_.seq_.substr(pos + 1, klen_);
				if(kCounts_[firstKmer] > occurenceCutOff_ &&
						kCounts_[nextKmer] > occurenceCutOff_){
					threadThroughSequenceHelperMate(firstKmer,
							nodes_[nodePositions_.at(firstKmer)],
							nextKmer,
							nodes_[nodePositions_.at(nextKmer)],
							internalCount,
							threadingSeqNameSecondMate,
							threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
				}
			}
		}
		if(pseq.mateSeqBase_.seq_.size() - stop > 1){
			if (stop > klen_) {
				stopPass = true;
			}

			uint32_t rangeStart = stop + 2 < klen_ ? 0 : stop + 2 - klen_;;
			uint32_t rangeStop = pseq.mateSeqBase_.seq_.size() - klen_  ;

			for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
				std::string firstKmer = pseq.mateSeqBase_.seq_.substr(pos, klen_);
				std::string nextKmer = pseq.mateSeqBase_.seq_.substr(pos + 1, klen_);
				if(kCounts_[firstKmer] > occurenceCutOff_ &&
						kCounts_[nextKmer] > occurenceCutOff_){
					threadThroughSequenceHelperMate(firstKmer,
							nodes_[nodePositions_.at(firstKmer)],
							nextKmer,
							nodes_[nodePositions_.at(nextKmer)],
							internalCount,
							threadingSeqNameSecondMate,
							threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
				}
			}
		}
		if(stopPass && startPass){
			//add the kmers that then span the region
			uint32_t newKmerSize = stop - start + 2;
			if(newKmerSize >klen_){
				std::string expandedKmer = pseq.mateSeqBase_.seq_.substr(start - 1, newKmerSize);
				{
					//thread in
					std::string firstKmer = pseq.mateSeqBase_.seq_.substr(start - 2, klen_);
					if(kCounts_[firstKmer] > occurenceCutOff_ &&
							kCounts_[expandedKmer] > occurenceCutOff_){
						threadThroughSequenceHelperMate(firstKmer,
								nodes_[nodePositions_.at(firstKmer)],
								expandedKmer,
								nodes_[nodePositions_.at(expandedKmer)],
								internalCount,
								threadingSeqNameSecondMate,
								threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
					}
				}
				{
					//thread out
					std::string nextKmer = pseq.mateSeqBase_.seq_.substr(stop + 2 - klen_, klen_);
					if(kCounts_[expandedKmer] > occurenceCutOff_ &&
							kCounts_[nextKmer] > occurenceCutOff_){
						threadThroughSequenceHelperMate(expandedKmer,
								nodes_[nodePositions_.at(expandedKmer)],
								nextKmer,
								nodes_[nodePositions_.at(nextKmer)],
								internalCount,
								threadingSeqNameSecondMate,
								threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
					}
				}
			}
		}
	} else if(processedPositions.size() > 1){
		//front
		{
			bool startPass = false;
			uint32_t start = processedPositions.front().start_;
			uint32_t stop = processedPositions.front().end_;
			if (start > 1) {
				if (start < pseq.mateSeqBase_.seq_.size() + 1 - klen_) {
					startPass = true;
				}
				for (auto pos : iter::range<uint32_t>(0, std::min<uint32_t>(start - 2, pseq.mateSeqBase_.seq_.size() + 1 - klen_))) {
					std::string firstKmer = pseq.mateSeqBase_.seq_.substr(pos, klen_);
					std::string nextKmer = pseq.mateSeqBase_.seq_.substr(pos + 1, klen_);
					if(kCounts_[firstKmer] > occurenceCutOff_ &&
							kCounts_[nextKmer] > occurenceCutOff_){
						threadThroughSequenceHelperMate(firstKmer,
								nodes_[nodePositions_.at(firstKmer)],
								nextKmer,
								nodes_[nodePositions_.at(nextKmer)],
								internalCount,
								threadingSeqNameSecondMate,
								threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
					}
				}
			}
			if(startPass){
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize >klen_){
					std::string expandedKmer = pseq.mateSeqBase_.seq_.substr(start - 1, newKmerSize);
					{
						//thread in
						std::string firstKmer = pseq.mateSeqBase_.seq_.substr(start - 2, klen_);
						if(kCounts_[firstKmer] > occurenceCutOff_ &&
								kCounts_[expandedKmer] > occurenceCutOff_){
							threadThroughSequenceHelperMate(firstKmer,
									nodes_[nodePositions_.at(firstKmer)],
									expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									internalCount,
									threadingSeqNameSecondMate,
									threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
						}
					}
					{
						//thread out
						std::string nextKmer = pseq.mateSeqBase_.seq_.substr(stop + 2 - klen_, klen_);
						if(kCounts_[expandedKmer] > occurenceCutOff_ &&
								kCounts_[nextKmer] > occurenceCutOff_){
							threadThroughSequenceHelperMate(expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									nextKmer,
									nodes_[nodePositions_.at(nextKmer)],
									internalCount,
									threadingSeqNameSecondMate,
									threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
						}
					}
				}
			}
		}

		//connect the positions
		for(const auto pos : iter::range<uint32_t>(processedPositions.size() - 1)){
			if(pos > 0){
				uint32_t start = processedPositions[pos].start_;
				uint32_t stop = processedPositions[pos].end_;
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize > klen_){
					std::string expandedKmer = pseq.mateSeqBase_.seq_.substr(start - 1, newKmerSize);
					{
						//thread in
						std::string firstKmer = pseq.mateSeqBase_.seq_.substr(start - 2, klen_);
						if(kCounts_[firstKmer] > occurenceCutOff_ &&
								kCounts_[expandedKmer] > occurenceCutOff_){
							threadThroughSequenceHelperMate(firstKmer,
									nodes_[nodePositions_.at(firstKmer)],
									expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									internalCount,
									threadingSeqNameSecondMate,
									threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
						}
					}
					{
						//thread out
						std::string nextKmer = pseq.mateSeqBase_.seq_.substr(stop + 2 - klen_, klen_);
						if(kCounts_[expandedKmer] > occurenceCutOff_ &&
								kCounts_[nextKmer] > occurenceCutOff_){
							threadThroughSequenceHelperMate(expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									nextKmer,
									nodes_[nodePositions_.at(nextKmer)],
									internalCount,
									threadingSeqNameSecondMate,
									threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
						}
					}
				}
			}
			uint32_t stop = processedPositions[pos].end_;
			uint32_t nextStart = processedPositions[pos + 1].start_;
			uint32_t rangeStart = stop + 2 < klen_ ? stop + 2 : stop + 2 - klen_;
			uint32_t rangeStop = std::min<uint32_t>(nextStart - 2 , pseq.mateSeqBase_.seq_.size() - klen_);
			for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
				std::string firstKmer = pseq.mateSeqBase_.seq_.substr(pos, klen_);
				std::string nextKmer = pseq.mateSeqBase_.seq_.substr(pos + 1, klen_);
				if(kCounts_[firstKmer] > occurenceCutOff_ &&
						kCounts_[nextKmer] > occurenceCutOff_){
					threadThroughSequenceHelperMate(firstKmer,
							nodes_[nodePositions_.at(firstKmer)],
							nextKmer,
							nodes_[nodePositions_.at(nextKmer)],
							internalCount,
							threadingSeqNameSecondMate,
							threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
				}
			}
		}
		//end
		{
			bool stopPass = false;
			uint32_t start = processedPositions.back().start_;
			uint32_t stop = processedPositions.back().end_;
			if(pseq.mateSeqBase_.seq_.size() - stop > 1){
				if (stop > klen_) {
					stopPass = true;
				}
			}
			if(stopPass){
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize > klen_){
					std::string expandedKmer = pseq.mateSeqBase_.seq_.substr(start - 1, newKmerSize);
					{
						//thread in
						std::string firstKmer = pseq.mateSeqBase_.seq_.substr(start - 2, klen_);
						if(kCounts_[firstKmer] > occurenceCutOff_ &&
								kCounts_[expandedKmer] > occurenceCutOff_){
							threadThroughSequenceHelperMate(firstKmer,
									nodes_[nodePositions_.at(firstKmer)],
									expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									internalCount,
									threadingSeqNameSecondMate,
									threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
						}
					}
					{
						//thread out
						std::string nextKmer = pseq.mateSeqBase_.seq_.substr(stop + 2 - klen_, klen_);
						if(kCounts_[expandedKmer] > occurenceCutOff_ &&
								kCounts_[nextKmer] > occurenceCutOff_){
							threadThroughSequenceHelperMate(expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									nextKmer,
									nodes_[nodePositions_.at(nextKmer)],
									internalCount,
									threadingSeqNameSecondMate,
									threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
						}
					}
				}
			}
			if(pseq.mateSeqBase_.seq_.size() - stop > 1){
				if (stop > klen_) {
					stopPass = true;
				}
				uint32_t rangeStart = stop + 2 < klen_ ? 0 : stop + 2 - klen_;;
				uint32_t rangeStop = pseq.mateSeqBase_.seq_.size() - klen_;
				for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
					std::string firstKmer = pseq.mateSeqBase_.seq_.substr(pos, klen_);
					std::string nextKmer = pseq.mateSeqBase_.seq_.substr(pos + 1, klen_);
					if(kCounts_[firstKmer] > occurenceCutOff_ &&
							kCounts_[nextKmer] > occurenceCutOff_){
						threadThroughSequenceHelperMate(firstKmer,
								nodes_[nodePositions_.at(firstKmer)],
								nextKmer,
								nodes_[nodePositions_.at(nextKmer)],
								internalCount,
								threadingSeqNameSecondMate,
								threadingSeqNameSecondMateIdx, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
					}
				}
			}
		}
	}
	return internalCount;


}

std::unordered_map<std::string, uint32_t> KmerPathwayGraphDev::threadThroughSequenceAdjust(const seqInfo & seq,
		const std::vector<StartEndPos> & processedPositions,
		const std::string & threadingSeqName){

	if(processedPositions.size() == 0){
		return threadThroughSequence(seq, threadingSeqName);
	}
	uint32_t threadingSeqNameIdx = readNames_.size();
	readNames_.emplace_back(threadingSeqName);
	std::unordered_map<std::string, uint32_t> internalCount;

	if(processedPositions.size() == 1){
		bool startPass = false;
		bool stopPass = false;
		uint32_t start = processedPositions.front().start_;
		uint32_t stop = processedPositions.front().end_;
		if (start > 1) {
			if (start < seq.seq_.size() + 1 - klen_) {
				startPass = true;
			}
			for (auto pos : iter::range<uint32_t>(0, std::min<uint32_t>(start - 2, seq.seq_.size() + 1 - klen_))) {
				std::string firstKmer = seq.seq_.substr(pos, klen_);
				std::string nextKmer = seq.seq_.substr(pos + 1, klen_);
				if(kCounts_[firstKmer] > occurenceCutOff_ &&
						kCounts_[nextKmer] > occurenceCutOff_){
					threadThroughSequenceHelper(firstKmer,
							nodes_[nodePositions_.at(firstKmer)],
							nextKmer,
							nodes_[nodePositions_.at(nextKmer)],
							internalCount,
							threadingSeqName,
							threadingSeqNameIdx);
//					nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//					nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//					addEdge(firstKmer, nextKmer, 1, threadingSeqNameIdx);
				}
			}
		}
		if(seq.seq_.size() - stop > 1){
			if (stop > klen_) {
				stopPass = true;
			}

			uint32_t rangeStart = stop + 2 < klen_ ? 0 : stop + 2 - klen_;;
			uint32_t rangeStop = seq.seq_.size() - klen_  ;

			for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
				std::string firstKmer = seq.seq_.substr(pos, klen_);
				std::string nextKmer = seq.seq_.substr(pos + 1, klen_);
				if(kCounts_[firstKmer] > occurenceCutOff_ &&
						kCounts_[nextKmer] > occurenceCutOff_){
					threadThroughSequenceHelper(firstKmer,
							nodes_[nodePositions_.at(firstKmer)],
							nextKmer,
							nodes_[nodePositions_.at(nextKmer)],
							internalCount,
							threadingSeqName,
							threadingSeqNameIdx);
//					nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//					nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//					addEdge(firstKmer, nextKmer, 1, threadingSeqNameIdx);
				}
			}
		}
		if(stopPass && startPass){
			//add the kmers that then span the region
			uint32_t newKmerSize = stop - start + 2;
			if(newKmerSize >klen_){
				std::string expandedKmer = seq.seq_.substr(start - 1, newKmerSize);
				{
					//thread in
					std::string firstKmer = seq.seq_.substr(start - 2, klen_);
					if(kCounts_[firstKmer] > occurenceCutOff_ &&
							kCounts_[expandedKmer] > occurenceCutOff_){
						threadThroughSequenceHelper(firstKmer,
								nodes_[nodePositions_.at(firstKmer)],
								expandedKmer,
								nodes_[nodePositions_.at(expandedKmer)],
								internalCount,
								threadingSeqName,
								threadingSeqNameIdx);
//						nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//						nodes_[nodePositions_.at(expandedKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//						addEdge(firstKmer, expandedKmer, 1, threadingSeqNameIdx);
					}
				}
				{
					//thread out
					std::string nextKmer = seq.seq_.substr(stop + 2 - klen_, klen_);
					if(kCounts_[expandedKmer] > occurenceCutOff_ &&
							kCounts_[nextKmer] > occurenceCutOff_){
						threadThroughSequenceHelper(expandedKmer,
								nodes_[nodePositions_.at(expandedKmer)],
								nextKmer,
								nodes_[nodePositions_.at(nextKmer)],
								internalCount,
								threadingSeqName,
								threadingSeqNameIdx);
//						nodes_[nodePositions_.at(expandedKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//						nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//						addEdge(expandedKmer, nextKmer, 1, threadingSeqNameIdx);
					}
				}
			}
		}
	} else if(processedPositions.size() > 1){
		//front
		{
			bool startPass = false;
			uint32_t start = processedPositions.front().start_;
			uint32_t stop = processedPositions.front().end_;
			if (start > 1) {
				if (start < seq.seq_.size() + 1 - klen_) {
					startPass = true;
				}
				for (auto pos : iter::range<uint32_t>(0, std::min<uint32_t>(start - 2, seq.seq_.size() + 1 - klen_))) {
					std::string firstKmer = seq.seq_.substr(pos, klen_);
					std::string nextKmer = seq.seq_.substr(pos + 1, klen_);
					if(kCounts_[firstKmer] > occurenceCutOff_ &&
							kCounts_[nextKmer] > occurenceCutOff_){
						threadThroughSequenceHelper(firstKmer,
								nodes_[nodePositions_.at(firstKmer)],
								nextKmer,
								nodes_[nodePositions_.at(nextKmer)],
								internalCount,
								threadingSeqName,
								threadingSeqNameIdx);
//						nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//						nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//						addEdge(firstKmer, nextKmer, 1, threadingSeqNameIdx);
					}
				}
			}
			if(startPass){
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize >klen_){
					std::string expandedKmer = seq.seq_.substr(start - 1, newKmerSize);
					{
						//thread in
						std::string firstKmer = seq.seq_.substr(start - 2, klen_);
						if(kCounts_[firstKmer] > occurenceCutOff_ &&
								kCounts_[expandedKmer] > occurenceCutOff_){
							threadThroughSequenceHelper(firstKmer,
									nodes_[nodePositions_.at(firstKmer)],
									expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									internalCount,
									threadingSeqName,
									threadingSeqNameIdx);
//							nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							nodes_[nodePositions_.at(expandedKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							addEdge(firstKmer, expandedKmer, 1, threadingSeqNameIdx);
						}
					}
					{
						//thread out
						std::string nextKmer = seq.seq_.substr(stop + 2 - klen_, klen_);
						if(kCounts_[expandedKmer] > occurenceCutOff_ &&
								kCounts_[nextKmer] > occurenceCutOff_){
							threadThroughSequenceHelper(expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									nextKmer,
									nodes_[nodePositions_.at(nextKmer)],
									internalCount,
									threadingSeqName,
									threadingSeqNameIdx);
//							nodes_[nodePositions_.at(expandedKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							addEdge(expandedKmer, nextKmer, 1, threadingSeqNameIdx);
						}
					}
				}
			}
		}

		//connect the positions
		for(const auto pos : iter::range<uint32_t>(processedPositions.size() - 1)){
			if(pos > 0){
				uint32_t start = processedPositions[pos].start_;
				uint32_t stop = processedPositions[pos].end_;
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize > klen_){
					std::string expandedKmer = seq.seq_.substr(start - 1, newKmerSize);
					{
						//thread in
						std::string firstKmer = seq.seq_.substr(start - 2, klen_);
						if(kCounts_[firstKmer] > occurenceCutOff_ &&
								kCounts_[expandedKmer] > occurenceCutOff_){
							threadThroughSequenceHelper(firstKmer,
									nodes_[nodePositions_.at(firstKmer)],
									expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									internalCount,
									threadingSeqName,
									threadingSeqNameIdx);
//							error
//							nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							nodes_[nodePositions_.at(expandedKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							addEdge(firstKmer, expandedKmer, 1, threadingSeqNameIdx);
						}
					}
					{
						//thread out
						std::string nextKmer = seq.seq_.substr(stop + 2 - klen_, klen_);
						if(kCounts_[expandedKmer] > occurenceCutOff_ &&
								kCounts_[nextKmer] > occurenceCutOff_){
							threadThroughSequenceHelper(expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									nextKmer,
									nodes_[nodePositions_.at(nextKmer)],
									internalCount,
									threadingSeqName,
									threadingSeqNameIdx);
//							error
//							nodes_[nodePositions_.at(expandedKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							addEdge(expandedKmer, nextKmer, 1, threadingSeqNameIdx);
						}
					}
				}
			}
			uint32_t stop = processedPositions[pos].end_;
			uint32_t nextStart = processedPositions[pos + 1].start_;
			uint32_t rangeStart = stop + 2 < klen_ ? stop + 2 : stop + 2 - klen_;
			uint32_t rangeStop = std::min<uint32_t>(nextStart - 2 , seq.seq_.size() - klen_);
			for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
				std::string firstKmer = seq.seq_.substr(pos, klen_);
				std::string nextKmer = seq.seq_.substr(pos + 1, klen_);
				if(kCounts_[firstKmer] > occurenceCutOff_ &&
						kCounts_[nextKmer] > occurenceCutOff_){
					threadThroughSequenceHelper(firstKmer,
							nodes_[nodePositions_.at(firstKmer)],
							nextKmer,
							nodes_[nodePositions_.at(nextKmer)],
							internalCount,
							threadingSeqName,
							threadingSeqNameIdx);
//					error
//					nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//					nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//					addEdge(firstKmer, nextKmer, 1, threadingSeqNameIdx);
				}
			}
		}
		//end
		{
			bool stopPass = false;
			uint32_t start = processedPositions.back().start_;
			uint32_t stop = processedPositions.back().end_;
			if(seq.seq_.size() - stop > 1){
				if (stop > klen_) {
					stopPass = true;
				}
			}
			if(stopPass){
				//add the kmers that then span the region
				uint32_t newKmerSize = stop - start + 2;
				if(newKmerSize > klen_){
					std::string expandedKmer = seq.seq_.substr(start - 1, newKmerSize);
					{
						//thread in
						std::string firstKmer = seq.seq_.substr(start - 2, klen_);
						if(kCounts_[firstKmer] > occurenceCutOff_ &&
								kCounts_[expandedKmer] > occurenceCutOff_){
							threadThroughSequenceHelper(firstKmer,
									nodes_[nodePositions_.at(firstKmer)],
									expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									internalCount,
									threadingSeqName,
									threadingSeqNameIdx);
//							error
//							nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							nodes_[nodePositions_.at(expandedKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							addEdge(firstKmer, expandedKmer, 1, threadingSeqNameIdx);
						}
					}
					{
						//thread out
						std::string nextKmer = seq.seq_.substr(stop + 2 - klen_, klen_);
						if(kCounts_[expandedKmer] > occurenceCutOff_ &&
								kCounts_[nextKmer] > occurenceCutOff_){
							threadThroughSequenceHelper(expandedKmer,
									nodes_[nodePositions_.at(expandedKmer)],
									nextKmer,
									nodes_[nodePositions_.at(nextKmer)],
									internalCount,
									threadingSeqName,
									threadingSeqNameIdx);
//							error
//							nodes_[nodePositions_.at(expandedKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//							addEdge(expandedKmer, nextKmer, 1, threadingSeqNameIdx);
						}
					}
				}
			}
			if(seq.seq_.size() - stop > 1){
				if (stop > klen_) {
					stopPass = true;
				}
				uint32_t rangeStart = stop + 2 < klen_ ? 0 : stop + 2 - klen_;;
				uint32_t rangeStop = seq.seq_.size() - klen_;
				for (auto pos : iter::range<uint32_t>(rangeStart, rangeStop )) {
					std::string firstKmer = seq.seq_.substr(pos, klen_);
					std::string nextKmer = seq.seq_.substr(pos + 1, klen_);
					if(kCounts_[firstKmer] > occurenceCutOff_ &&
							kCounts_[nextKmer] > occurenceCutOff_){
						threadThroughSequenceHelper(firstKmer,
								nodes_[nodePositions_.at(firstKmer)],
								nextKmer,
								nodes_[nodePositions_.at(nextKmer)],
								internalCount,
								threadingSeqName,
								threadingSeqNameIdx);
//						error
//						nodes_[nodePositions_.at(firstKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//						nodes_[nodePositions_.at(nextKmer)]->inReadNamesIdx_.emplace(threadingSeqNameIdx);
//						addEdge(firstKmer, nextKmer, 1, threadingSeqNameIdx);
					}
				}
			}
		}
	}
	return internalCount;
}


} // namespace njhseq
