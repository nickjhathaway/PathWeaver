/*
 * CoverageEstimator.cpp
 *
 *  Created on: Apr 23, 2019
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

#include "CoverageEstimator.hpp"
namespace njhseq {



CoverageEstimator::CovInfoPerPos::CovInfoPerPos() {
}

CoverageEstimator::CovInfoPerPos::CovInfoPerPos(uint32_t pos, uint32_t minPos,
		double sdCov, double avgCov) :
		pos_(pos), minPos_(minPos), sdCov_(sdCov), avgCov_(avgCov) {
}


CoverageEstimator::CoverageEstimatedResult CoverageEstimator::estimateCov(
		const seqInfo & seq,  KmerPathwayGraph & estimatingGraph) {
	CoverageEstimatedResult ret;
	MetaDataInName meta;
	if (MetaDataInName::nameHasMetaData(seq.name_)) {
		meta = MetaDataInName(seq.name_);
	}
	return estimateCov(seq.seq_, meta, estimatingGraph);
}

CoverageEstimator::CoverageEstimatedResult CoverageEstimator::estimateCov(const std::string & seq, const MetaDataInName & meta,  KmerPathwayGraph & estimatingGraph){
	CoverageEstimatedResult ret;
	if(seq.length() >= estimatingGraph.klen_){
		for(const auto pos : iter::range(seq.length() - estimatingGraph.klen_ + 1)){
//			//this is terrible, because it could fail
//			if(estimatingGraph.kCounts_.end() == estimatingGraph.kCounts_.find(seq.substr(pos, estimatingGraph.klen_))){
//				std::cout << "seq:     " << seq << std::endl;
//				std::cout << "seqsize: " << seq.size() << std::endl;
//				std::cout << "pos:     " << pos << std::endl;
//				std::cout << "klen:    " << estimatingGraph.klen_ << std::endl;
//				std::cout << "kmer:    " << seq.substr(pos, estimatingGraph.klen_) << std::endl;
//				//exit(1);
//			}
//			ret.allCounts_.emplace_back(estimatingGraph.kCounts_.at(seq.substr(pos, estimatingGraph.klen_)));
			ret.allCounts_.emplace_back(estimatingGraph.kCounts_[seq.substr(pos, estimatingGraph.klen_)]);
		}
		bool headTrimStatus = false;
		bool tailTrimStatus = false;
		if (meta.containsMeta("headTrimStatus") && meta.containsMeta("tailTrimStatus")) {
			headTrimStatus = meta.getMeta<bool>("headTrimStatus");
			tailTrimStatus = meta.getMeta<bool>("tailTrimStatus");
		} else if (meta.containsMeta("trimStatus")) {
			headTrimStatus = meta.getMeta<bool>("trimStatus");
			tailTrimStatus = meta.getMeta<bool>("trimStatus");
		}

		bool headless = false;
		if(meta.containsMeta("headless")){
			headless = meta.getMeta<bool>("headless");
		}
		bool tailless = false;
		if(meta.containsMeta("tailless")){
			tailless = meta.getMeta<bool>("tailless");
		}
		auto currentStartAdjust = (headless && !headTrimStatus) ? estimatingGraph.klen_ - 1 : 0;
		auto currentStopAdjust  = (tailless && !tailTrimStatus) ? estimatingGraph.klen_     : 0;
		if (ret.allCounts_.size() > currentStartAdjust + currentStopAdjust) {
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			for(const auto & pos : iter::range<uint64_t>(currentStartAdjust, ret.allCounts_.size() - currentStopAdjust)){
				CovInfoPerPos minSd;
				uint64_t startSurroundPos = pos + 1 > estimatingGraph.klen_ ? pos + 1 - estimatingGraph.klen_: 0;
				for(const auto & surPos : iter::range(startSurroundPos, pos + 1)){
					auto subSet = getSubVector(ret.allCounts_, surPos, std::min<uint32_t>(pos + 1,std::min<uint32_t>(estimatingGraph.klen_, ret.allCounts_.size() - surPos)));
					double currentSd = vectorStandardDeviationPop(subSet);
					if(currentSd < minSd.sdCov_){
						minSd = CovInfoPerPos(pos, surPos, currentSd, vectorMean(subSet));
					}
				}
				ret.coverageInfos_.emplace_back(minSd);
			}

			for(const auto & cInfo : ret.coverageInfos_){
				if(cInfo.avgCov_ < ret.minCov_.avgCov_){
					ret.minCov_ = cInfo;
				}
			}
			//meta.addMeta("estimatedPerBaseCoverage", minAvg, true);
		} else {
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			ret.minCov_.avgCov_ = vectorMean(ret.allCounts_);
			//meta.addMeta("estimatedPerBaseCoverage", vectorMean(allCounts), true);
		}
		//meta.resetMetaInName(seq.name_);
	}
	return ret;
}

}  // namespace njhseq
