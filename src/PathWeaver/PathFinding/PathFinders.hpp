#pragma once
/*
 * PathFinders.hpp
 *
 *  Created on: May 19, 2018
 *      Author: nick
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

#include <TwoBit.h>

#include <njhseq/BamToolsUtils.h>
#include <njhseq/objects/BioDataObject.h>
#include <njhseq/GenomeUtils.h>
#include <njhcpp/concurrency/LockableJsonLog.hpp>

#include "PathWeaver/objects/dataContainers/graphs.h"
#include "PathWeaver/PathFinding/HaploPathFinder.hpp"
#include "PathWeaver/PathFinding/WayFindingUtils.h"


namespace njhseq {



struct PathFinderFromSeqsRes{
	Json::Value log_;

	writeOutTandemsAndOptionallyStitchRes previousDetTandems_;
};

PathFinderFromSeqsRes PathFinderFromSeqs(
		const BamExtractor::ExtractedFilesOpts & inOpts,
		const bfs::path & workingDir,
		const std::string & sampName,
		const HaploPathFinder::PathFinderCorePars & pars,
		const std::unique_ptr<MultipleGroupMetaData> & meta,
		const PathFinderFromSeqsRes & previousRunRes = PathFinderFromSeqsRes{});



}  // namespace njhseq
