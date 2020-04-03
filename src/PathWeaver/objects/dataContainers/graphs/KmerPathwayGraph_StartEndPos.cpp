/*
 * KmerPathwayGraph_StartEndPos.cpp
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

#include "KmerPathwayGraph.hpp"


namespace njhseq {




KmerPathwayGraph::StartEndPos::StartEndPos(uint32_t start, uint32_t end) :
		start_(start), end_(end) {
	if (start_ >= end_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "start_: " << start_
				<< "greater than or equal to end_: " << end_ << "\n";
		throw std::runtime_error { ss.str() };
	}
}

uint32_t KmerPathwayGraph::StartEndPos::len() const {
	return end_ - start_;
}

std::vector<KmerPathwayGraph::StartEndPos> KmerPathwayGraph::StartEndPos::merge(
		std::vector<StartEndPos> positions, uint32_t expand) {
	njh::sort(positions, [](const StartEndPos & p1, const StartEndPos & p2) {
		if(p1.start_ == p2.start_) {
			return p1.end_ < p2.end_;
		} else {
			return p1.start_ < p2.start_;
		}
	});
	std::vector<StartEndPos> ret;
	for (const auto & position : positions) {
		if (ret.empty()) {
			ret.emplace_back(position);
		} else {
			uint32_t previousEnd = ret.back().end_ + expand;
			uint32_t currentStart =
					position.start_ < expand ? 0 : position.start_ - expand;
			if (currentStart <= previousEnd) {
				ret.back().end_ = position.end_;
			} else {
				ret.emplace_back(position);
			}
		}
	}
	return ret;
}




} // namespace njhseq
