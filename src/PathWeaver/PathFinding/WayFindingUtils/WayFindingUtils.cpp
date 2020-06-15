/*
 * WayFindingUtils.cpp
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


#include "WayFindingUtils.hpp"

namespace njhseq {

OptimizationReconResult::Params::Params(const uint32_t klen,
		const uint32_t kcut, const uint32_t shortNumber) :
		klen_(klen), kcut_(kcut), shortNumber_(shortNumber) {

}
Json::Value OptimizationReconResult::Params::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["kcut_"] = njh::json::toJson(kcut_);
	ret["klen_"] = njh::json::toJson(klen_);
	ret["shortNumber_"] = njh::json::toJson(shortNumber_);
	return ret;
}

OptimizationReconResult::Dirs::Dirs(const bfs::path & klenDir,
		const bfs::path & kcutDir, const bfs::path & shortTipDir) :
		klenDir_(klenDir), kcutDir_(kcutDir), shortTipDir_(shortTipDir) {

}

Json::Value OptimizationReconResult::Dirs::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["klenDir_"] = njh::json::toJson(klenDir_);
	ret["kcutDir_"] = njh::json::toJson(kcutDir_);
	ret["shortTipDir_"] = njh::json::toJson(shortTipDir_);
	return ret;
}

OptimizationReconResult::OptimizationReconResult(const Params runParams,
		const Dirs & runDirs, uint32_t kcutIterNumber, uint32_t shortTipIterNumber) :
		runParams_(runParams), runDirs_(runDirs), kcutIterNumber_(kcutIterNumber), shortTipIterNumber_(
				shortTipIterNumber) {

}

Json::Value OptimizationReconResult::toJson() const {
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));
	ret["runParams_"] = njh::json::toJson(runParams_);
	ret["runDirs_"] = njh::json::toJson(runDirs_);
	ret["kcutIterNumber_"] = njh::json::toJson(kcutIterNumber_);
	ret["shortTipIterNumber_"] = njh::json::toJson(shortTipIterNumber_);
	ret["readLengthAverage_"] = njh::json::toJson(readLengthAverage_);
	ret["readLengthMedian_"] = njh::json::toJson(readLengthMedian_);

	ret["headlessCount_"] = njh::json::toJson(headlessCount_);
	ret["taillessCount_"] = njh::json::toJson(taillessCount_);
	ret["optimalCount_"] = njh::json::toJson(optimalCount_);

	ret["headlessCountBelowLen_"] = njh::json::toJson(headlessCountBelowLen_);
	ret["taillessCountBelowLen_"] = njh::json::toJson(taillessCountBelowLen_);
	ret["optimalCountBelowLen_"] = njh::json::toJson(optimalCountBelowLen_);

	ret["internalNodesCountBelowLen_"] = njh::json::toJson(internalNodesCountBelowLen_);

	ret["percentOfInputUsed_"] = njh::json::toJson(percentOfInputUsed_);
	ret["totalInputReads_"] = njh::json::toJson(totalInputReads_);
	ret["numberOfFinalFilteredSeqs_"] = njh::json::toJson(
			numberOfFinalFilteredSeqs_);

	return ret;
}


std::function<
			bool(const OptimizationReconResult&, const OptimizationReconResult&)> OptimizationReconResult::sortFunc =
			[](const OptimizationReconResult& optRes1,const OptimizationReconResult& optRes2) {
				if(optRes1.runParams_.kcut_ == optRes2.runParams_.kcut_) {
					if(optRes1.runParams_.shortNumber_ == optRes2.runParams_.shortNumber_) {
						return optRes1.runParams_.klen_ < optRes2.runParams_.klen_;
					} else {
						return optRes1.runParams_.shortNumber_ < optRes2.runParams_.shortNumber_;
					}
				} else {
					return optRes1.runParams_.kcut_ < optRes2.runParams_.kcut_;
				}
			};



std::vector<OptimizationReconResult> OptimizationReconResult::getBestResults(const std::vector<OptimizationReconResult> & allResults,
			double percentageUsed,
			double percentagedUsedStepDown,
			bool useFullOptimalCount){
		std::vector<OptimizationReconResult> bestOptResults;
		uint32_t optimalCount = std::numeric_limits<uint32_t>::max();
		uint32_t internalCount = std::numeric_limits<uint32_t>::max(); //this is used to break ties of optimal counts

		double percentUsed = percentageUsed;
		while(bestOptResults.empty() && percentUsed >= 0){
			for(const auto & optResTest : allResults){
				if(optResTest.percentOfInputUsed_ > percentUsed && 0 != optResTest.numberOfFinalFilteredSeqs_ ){
					uint32_t currentOptCount = useFullOptimalCount ? optResTest.optimalCount_ : optResTest.optimalCountBelowLen_;
					//uint32_t currentInternalCountBelowCount = useFullOptimalCount ? 1 : optResTest.internalNodesCountBelowLen_;
					uint32_t currentInternalCountBelowCount = optResTest.internalNodesCountBelowLen_;
					if(bestOptResults.empty()){
						bestOptResults.emplace_back(optResTest);
						optimalCount = currentOptCount;
						internalCount = currentInternalCountBelowCount;
					}else{
						if(currentOptCount < optimalCount){
							bestOptResults.clear();
							bestOptResults.emplace_back(optResTest);
							optimalCount = currentOptCount;
							internalCount = currentInternalCountBelowCount;
						}else if(currentOptCount == optimalCount){
							if(currentInternalCountBelowCount < internalCount){
								bestOptResults.clear();
								bestOptResults.emplace_back(optResTest);
								optimalCount = currentOptCount;
								internalCount = currentInternalCountBelowCount;
							}else{
								bestOptResults.emplace_back(optResTest);
							}
						}
					}
				}
			}
			if(percentUsed >0 && percentagedUsedStepDown > percentUsed){
				percentUsed = 0;
			}else{
				percentUsed -= percentagedUsedStepDown;
			}
		} //end while
		njh::sort(bestOptResults,OptimizationReconResult::sortFunc);
		return bestOptResults;
	}


}  // namespace njhseq
