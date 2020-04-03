#pragma once
/*
 * WeaverRunner.h
 *
 *  Created on: Feb 4, 2015
 *      Author: nickhathaway
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

#include <njhcpp/progutils.h>
#include <njhseq.h>

namespace njhseq {

class WeaverRunner : public njh::progutils::ProgramRunner {
 public:
	WeaverRunner();

	static int BamExtractPathawaysFromRegion(const njh::progutils::CmdArgs & inputCommands);
	static int SeqsExtractPathaways(const njh::progutils::CmdArgs & inputCommands);
	static int ExtractPathWaysReadsFallingInMultipleRegions(const njh::progutils::CmdArgs & inputCommands);

};

} /* namespace njhseq */
