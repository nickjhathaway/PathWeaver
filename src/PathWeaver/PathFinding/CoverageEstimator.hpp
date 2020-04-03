#pragma once

/*
 * CoverageEstimator.hpp
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


#include "PathWeaver/objects/dataContainers/graphs/KmerPathwayGraph.hpp"


namespace njhseq {

class CoverageEstimator {
public:
	struct CovInfoPerPos {
		CovInfoPerPos();
		CovInfoPerPos(uint32_t pos, uint32_t minPos, double sdCov, double avgCov);
		uint32_t pos_ { std::numeric_limits<uint32_t>::max() };
		uint32_t minPos_ { std::numeric_limits<uint32_t>::max() };
		double sdCov_ { std::numeric_limits<double>::max() };
		double avgCov_ { std::numeric_limits<double>::max() };
	};

	struct CoverageEstimatedResult {
		std::vector<uint32_t> allCounts_;
		std::vector<CovInfoPerPos> coverageInfos_;
		CovInfoPerPos minCov_;
	};

	static CoverageEstimatedResult estimateCov(const seqInfo & seq,
			 KmerPathwayGraph & estimatingGraph);

	static CoverageEstimatedResult estimateCov(const std::string & seq,
			const MetaDataInName & meta, KmerPathwayGraph & estimatingGraph);



};

}  // namespace njhseq
