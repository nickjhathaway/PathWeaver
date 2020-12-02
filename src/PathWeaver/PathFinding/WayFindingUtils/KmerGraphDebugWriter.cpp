/*
 * KmerGraphDebugWriter.cpp
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


#include "KmerGraphDebugWriter.hpp"
namespace njhseq {

KmerGraphDebugWriter::KmerGraphDebugWriter(const bfs::path &currentDir,
		const KmerPathwayGraph & mainGraph,
					KmerPathwayGraph & covEstimatorGraph):currentDir_(currentDir), mainGraph_(mainGraph),
					covEstimatorGraph_(covEstimatorGraph){}


void KmerGraphDebugWriter::writeOutDotsAndSeqs(
		const std::string & nameStub) {
	OutputStream rectDotOut(
			njh::files::make_path(currentDir_,
					leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
							+ "-rect.dot"));
	mainGraph_.writeRectangleWithEstimatedCovDot(rectDotOut, covEstimatorGraph_);
	OutputStream rectDotNoLabelsOut(
			njh::files::make_path(currentDir_,
					leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
							+ "-rect-noLabels.dot"));
	mainGraph_.writeRectangleWithEstimatedCovDot(rectDotNoLabelsOut,
			covEstimatorGraph_, true);
	auto outSeqsOpts = SeqIOOptions::genFastaOut(
			njh::files::make_path(currentDir_,
					leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
							+ "-seqs.fasta"));
	SeqOutput::write(mainGraph_.nodesToSeqs(false), outSeqsOpts);
	auto outSeqsWithHeadsAndTailsOpts = SeqIOOptions::genFastaOut(
			njh::files::make_path(currentDir_,
					leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
							+ "-seqsWithHeadsTailsAddedOn.fasta"));
	SeqOutput::write(mainGraph_.nodesToSeqs(true), outSeqsWithHeadsAndTailsOpts);

	if(writeEdgeInfo_){
		OutputStream rectDotEdgeInfoOut(
				njh::files::make_path(currentDir_,
						leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
								+ "-rect-edgeInfo.tab.txt"));
		rectDotEdgeInfoOut << njh::conToStr(mainGraph_.writeEdgesWithNamesHeader(), "\t") << std::endl;
		mainGraph_.writeEdgesWithNames(rectDotEdgeInfoOut);
	}

	++rectGraphCount_;
}


KmerGraphDebugWriterDev::KmerGraphDebugWriterDev(const bfs::path &currentDir,
		const KmerPathwayGraphDev & mainGraph,
		KmerPathwayGraphDev & covEstimatorGraph):currentDir_(currentDir), mainGraph_(mainGraph),
					covEstimatorGraph_(covEstimatorGraph){}


void KmerGraphDebugWriterDev::writeOutDotsAndSeqs(
		const std::string & nameStub) {
	OutputStream rectDotOut(
			njh::files::make_path(currentDir_,
					leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
							+ "-rect.dot"));
	mainGraph_.writeRectangleWithEstimatedCovDot(rectDotOut, covEstimatorGraph_);
	OutputStream rectDotNoLabelsOut(
			njh::files::make_path(currentDir_,
					leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
							+ "-rect-noLabels.dot"));
	mainGraph_.writeRectangleWithEstimatedCovDot(rectDotNoLabelsOut,
			covEstimatorGraph_, true);
	auto outSeqsOpts = SeqIOOptions::genFastaOut(
			njh::files::make_path(currentDir_,
					leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
							+ "-seqs.fasta"));
	SeqOutput::write(mainGraph_.nodesToSeqs(false), outSeqsOpts);
	auto outSeqsWithHeadsAndTailsOpts = SeqIOOptions::genFastaOut(
			njh::files::make_path(currentDir_,
					leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
							+ "-seqsWithHeadsTailsAddedOn.fasta"));
	SeqOutput::write(mainGraph_.nodesToSeqs(true), outSeqsWithHeadsAndTailsOpts);

	if(writeEdgeInfo_){
		OutputStream rectDotEdgeInfoOut(
				njh::files::make_path(currentDir_,
						leftPadNumStr<uint32_t>(rectGraphCount_, 9999) + "-" + nameStub
								+ "-rect-edgeInfo.tab.txt"));
		rectDotEdgeInfoOut << njh::conToStr(mainGraph_.writeEdgesWithNamesHeader(), "\t") << std::endl;
		mainGraph_.writeEdgesWithNames(rectDotEdgeInfoOut);
	}

	++rectGraphCount_;
}


}  // namespace njhseq
