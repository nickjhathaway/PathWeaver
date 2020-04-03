/*
 * ContigsCompareGraph.cpp
 *
 *  Created on: Nov 5, 2019
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


#include "ContigsCompareGraph.hpp"

namespace njhseq {

ContigsCompareGraph::edge::edge(const std::shared_ptr<node> & head,
		const std::shared_ptr<node> & tail, uint32_t cnt, const std::vector<ConnectorInfo> & connectorInfos) :
		head_(head), tail_(tail), cnt_(cnt), connectorInfos_(connectorInfos) {
}



std::string ContigsCompareGraph::edge::createUid() const {
	return head_.lock()->uid_ + "_" + tail_.lock()->uid_;
}



}  // namespace njhseq
