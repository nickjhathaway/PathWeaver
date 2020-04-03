#pragma once

/*
 * KmerGraphDebugWriter.hpp
 *
 *  Created on: May 20, 2019
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

class KmerGraphDebugWriter {
public:

	KmerGraphDebugWriter(const bfs::path &currentDir,
			const KmerPathwayGraph & mainGraph, KmerPathwayGraph & covEstimatorGraph);

	uint32_t rectGraphCount_ = 0;
	bool writeEdgeInfo_ = false;
	bfs::path currentDir_;

	const KmerPathwayGraph & mainGraph_;
	KmerPathwayGraph & covEstimatorGraph_;

	void writeOutDotsAndSeqs(const std::string & nameStub);
};




}  // namespace njhseq
