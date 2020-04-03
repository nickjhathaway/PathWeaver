/*
 * PathFinders.cpp
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

#include "PathFinders.hpp"
#include "PathWeaver/PathFinding/CoverageEstimator.hpp"

namespace njhseq {


PathFinderFromSeqsRes PathFinderFromSeqs(
		const BamExtractor::ExtractedFilesOpts & inOpts,
		const bfs::path & workingDir,
		const std::string & sampName,
		const HaploPathFinder::PathFinderCorePars & extractionPars,
		const std::unique_ptr<MultipleGroupMetaData> & meta,
		const PathFinderFromSeqsRes & previousRunRes) {
	PathFinderFromSeqsRes ret;
	//for optimization results
	std::vector<OptimizationReconResult> allOptRunResults;
	//create a directory for the working directory
	auto finalCurrentDir = njh::files::makeDir(workingDir, njh::files::MkdirPar(sampName));
	//start log file
	ret.log_["date"] = njh::json::toJson(njh::getCurrentDateFull());
	ret.log_["sample"] =  sampName;
	ret.log_["pairedOpts"] =              njh::json::toJson(inOpts.inPairs_.toJson());
	ret.log_["pairedMateUnmappedOpts"] =  njh::json::toJson(inOpts.inPairsMateUnmapped_.toJson());
	ret.log_["singleOpts"] =              njh::json::toJson(inOpts.inUnpaired_.toJson());
	ret.log_["HaploPathFinder-ExtractParams"] =  njh::json::toJson(extractionPars);
	ret.log_["lenCutOff"] =  extractionPars.lenCutOff;
	njh::stopWatch watch;
	watch.setLapName(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps(), 10000), ": Pre-filter Seqs") );
	if(extractionPars.needToRePair_){
		//re-pair reads, this is necessary because a mate could have recruited in a previous step
		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": Re-pairing Sequences" ) );
		repairSeqsAfterRecruitmentPars rePairPars;
		rePairPars.singlesFnp_ = inOpts.inUnpaired_.firstName_;
		rePairPars.r1Fnp_ = inOpts.inPairs_.firstName_;
		rePairPars.r2Fnp_ = inOpts.inPairs_.secondName_;
		auto rePairRes = repairSeqsAfterRecruitment(rePairPars);
		ret.log_["numberOfPairsRePaired"] = rePairRes.numberOfPairsRePaired_;
	}

	//generate output for processed filtered sequences
	preprocessSeqsForWayFindingPars preprocessPars;
	preprocessPars.filteredPairedOpts = SeqIOOptions::genPairedOut(njh::files::make_path(finalCurrentDir,         "filteredExtractedPairs"));
	preprocessPars.filteredSingletOuts = SeqIOOptions::genFastqOut(njh::files::make_path(finalCurrentDir,         "filteredSingles"));
	preprocessPars.filteredOff_pairedOpts = SeqIOOptions::genPairedOut(njh::files::make_path(finalCurrentDir,     "filteredOff_extractedPairs"));
	preprocessPars.filteredOff_singletOuts = SeqIOOptions::genFastqOut(njh::files::make_path(finalCurrentDir,     "filteredOff_extractedSingles"));
	preprocessPars.filteredOffDups_pairedOpts = SeqIOOptions::genPairedOut(njh::files::make_path(finalCurrentDir, "filteredOffDups_extractedPairs"));
	preprocessPars.filteredOffDups_singletOuts = SeqIOOptions::genFastqOut(njh::files::make_path(finalCurrentDir, "filteredOffDups_extractedSingles"));

	//pre-process sequences this includes (use to include removes seqs with Ns), optionally trimming bad quality, and optionally removing duplicated sequences



	auto preprocessResults = preprocessSeqsForWayFinding(preprocessPars, inOpts, extractionPars, sampName);
	ret.log_["preProcessFilteredInfo"] = preprocessResults.filteredInfo();
//	std::cout << preprocessResults.filteredInfo() << std::endl;
//	exit(1);

	//stitching and global tandem repeat determination
	watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": Processing Sequences For Tandems  - Finding Global Tandems and Stitching" ) );
	//generate the input options for next step
	// global tandem info file

	bfs::path tandemInfoFnp = njh::files::make_path(finalCurrentDir, "tandemInfos.tab.txt");
	writeOutTandemsAndOptionallyStitchPars findTanStitchPars{OutOptions(tandemInfoFnp), finalCurrentDir};

	findTanStitchPars.usedFilteredPairedOpts = SeqIOOptions::genPairedIn(preprocessResults.usedFilteredPairedR1Fnp, preprocessResults.usedFilteredPairedR2Fnp);
	findTanStitchPars.usedFilteredSinglesOpts = SeqIOOptions::genFastqIn(preprocessResults.usedFilteredSinglesFnp);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	ret.previousDetTandems_ = writeOutTandemsAndOptionallyStitch(findTanStitchPars, preprocessResults.maxInputSeqLen, extractionPars, previousRunRes.previousDetTandems_);
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//for breaking up Ns
	std::string patStr = "N+";
	njh::PatPosFinder pFinder(patStr);

//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(const auto currentKLen : extractionPars.kmerLengths){
//	while(currentKLen <= kmerLenStop){
	//while(currentKLen <= kmerLenStop && !foundTrimmed){
//		//
//		std::cout << "\tcurrentKLen: " << currentKLen << ", kmerLenStop: " << kmerLenStop << std::endl;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		bfs::path currentKmerDirectory = njh::files::makeDir(finalCurrentDir, njh::files::MkdirPar(sampName + njh::pasteAsStr("_klen-", leftPadNumStr(currentKLen, vectorMaximum(extractionPars.kmerLengths))) ));;
		Json::Value & klenLog = ret.log_[sampName + njh::pasteAsStr("_klen-", leftPadNumStr(currentKLen, vectorMaximum(extractionPars.kmerLengths)))];
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		watch.startNewLap(
				njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen, "-",
				" Pre-filtering Sequences") );
		SeqIOOptions pairedOpts = SeqIOOptions::genPairedOut(
				njh::files::make_path(currentKmerDirectory, "extractedPairs"));
		SeqIOOptions singletOuts = SeqIOOptions::genFastqOut(
				njh::files::make_path(currentKmerDirectory, "extractedSingles"));
		SeqIOOptions unused_pairedOpts = SeqIOOptions::genPairedOut(
				njh::files::make_path(currentKmerDirectory, "unused_extractedPairs"));
		SeqIOOptions unused_singletOuts = SeqIOOptions::genFastqOut(
				njh::files::make_path(currentKmerDirectory, "unused_extractedSingles"));
		SeqIOOptions pairedTandemsOpts = SeqIOOptions::genPairedOut(
					njh::files::make_path(currentKmerDirectory, "extractedPairsWithTandems"));
		SeqIOOptions pairedNoTandemsOpts = SeqIOOptions::genPairedOut(
					njh::files::make_path(currentKmerDirectory, "extractedPairsWithNoTandems"));
		SeqIOOptions singletTandemsOuts = SeqIOOptions::genFastqOut(
				njh::files::make_path(currentKmerDirectory, "extractedSinglesWithTandems"));
		SeqIOOptions singletNoTandemsOuts = SeqIOOptions::genFastqOut(
				njh::files::make_path(currentKmerDirectory, "extractedSinglesWithNoTandems"));
		bfs::path usedSinglesFnp;
		bfs::path usedPairedR1Fnp;
		bfs::path usedPairedR2Fnp;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		{
			SeqOutput writer(pairedOpts);
			SeqOutput singleWriter(singletOuts);
			SeqOutput unused_writer(unused_pairedOpts);
			SeqOutput unused_singleWriter(unused_singletOuts);
			if (findTanStitchPars.usedFilteredPairedOpts.inExists()) {
				PairedRead pSeq;
				SeqInput pairedReader(findTanStitchPars.usedFilteredPairedOpts);
				pairedReader.openIn();
				while (pairedReader.readNextRead(pSeq)) {

					bool firstSkipped = KmerPathwayGraph::skipInputSeqForKCount(pSeq.seqBase_.seq_, currentKLen);
					bool secondSkipped = KmerPathwayGraph::skipInputSeqForKCount(pSeq.mateSeqBase_.seq_, currentKLen);
					if(extractionPars.trimSeqsWithNs_){
						auto breaksFirst = readVecTrimmer::breakUpSeqOnPat(pSeq.seqBase_, pFinder);
						{
							uint32_t maxFirstBreakSize = 0;
							uint32_t maxPosition = 0;
							uint32_t pos = 0;
							for(const auto & b : breaksFirst){
								if(len(b.seqBase_) > maxFirstBreakSize){
									maxFirstBreakSize = len(b.seqBase_);
									maxPosition = pos;
								}
								++pos;
							}
							if(maxFirstBreakSize > currentKLen + 1){
								firstSkipped = false;
								pSeq.seqBase_ = breaksFirst[maxPosition].seqBase_;
							}
						}
						auto breaksSecond = readVecTrimmer::breakUpSeqOnPat(pSeq.mateSeqBase_, pFinder);
						{
							uint32_t maxSecondBreakSize = 0;
							uint32_t maxPosition = 0;
							uint32_t pos = 0;
							for(const auto & b : breaksSecond){
								if(len(b.seqBase_) > maxSecondBreakSize){
									maxSecondBreakSize = len(b.seqBase_);
									maxPosition = pos;
								}
								++pos;
							}
							if(maxSecondBreakSize > currentKLen + 1){
								secondSkipped = false;
								pSeq.mateSeqBase_ = breaksSecond[maxPosition].seqBase_;
							}
						}
					}
					if (!firstSkipped && !secondSkipped) {
						writer.openWrite(pSeq);
					}else if(!firstSkipped && secondSkipped){
						singleWriter.openWrite(pSeq.seqBase_);
						unused_singleWriter.openWrite(pSeq.mateSeqBase_);
					}else if(firstSkipped && !secondSkipped){
						singleWriter.openWrite(pSeq.mateSeqBase_);
						unused_singleWriter.openWrite(pSeq.seqBase_);
					}else{
						unused_writer.openWrite(pSeq);
					}
				}
			}
			if (findTanStitchPars.usedFilteredSinglesOpts.inExists()) {
				seqInfo seq;
				SeqInput singleReader(findTanStitchPars.usedFilteredSinglesOpts);
				singleReader.openIn();
				while (singleReader.readNextRead(seq)) {
					bool skip = KmerPathwayGraph::skipInputSeqForKCount(seq.seq_,
							currentKLen);

					auto breaksFirst = readVecTrimmer::breakUpSeqOnPat(seq, pFinder);
					{
						uint32_t maxFirstBreakSize = 0;
						uint32_t maxPosition = 0;
						uint32_t pos = 0;
						for(const auto & b : breaksFirst){
							if(len(b.seqBase_) > maxFirstBreakSize){
								maxFirstBreakSize = len(b.seqBase_);
								maxPosition = pos;
							}
							++pos;
						}
						if(maxFirstBreakSize > currentKLen + 1){
							skip = false;
							seq = breaksFirst[maxPosition].seqBase_;
						}
					}
					if (!skip) {
						singleWriter.openWrite(seq);
					}else{
						unused_singleWriter.openWrite(seq);
					}
				}
			}
			if(writer.outOpen()){
				usedPairedR1Fnp = writer.getPrimaryOutFnp();
				usedPairedR2Fnp = writer.getSecondaryOutFnp();
			}
			if(singleWriter.outOpen()){
				usedSinglesFnp = singleWriter.getPrimaryOutFnp();
			}
			writer.closeOut();
			singleWriter.closeOut();
			unused_writer.closeOut();
			unused_singleWriter.closeOut();
		}
		SeqInput pairedReader(
				SeqIOOptions::genPairedIn(pairedOpts.getPriamryOutName(),
						pairedOpts.getSecondaryOutName()));
		SeqInput singleReader(
				SeqIOOptions::genFastqIn(singletOuts.getPriamryOutName()));


		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen,
				"- Processing Sequences For Tandems - Pre-Processing Tandems" ) );
		//uint32_t klength = currentKLen;
		std::vector<TandemRepeatPlusAdjustedSize> allTandems = readInTandemsWithPartInfo(tandemInfoFnp);
		auto processedTandems = processGlobalTandems(allTandems, currentKLen, extractionPars);
		processedTandems.writeOutInfos(currentKmerDirectory, false);
		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen,
				"- Processing Sequences For Tandems - Finding Tandems In Seqs" ) );
		//files for processed tandems
		KmerExpansionHelper expansionHelperInitial(processedTandems.tandemMots,
				processedTandems.tandemsAltMots, currentKLen,
				extractionPars.useSmallerThanKlenExpandedPositions_);
		uint32_t seqsWithTandems = 0;
		SeqOutput tandemPairWriter(pairedTandemsOpts);
		SeqOutput tandemSingleWriter(singletTandemsOuts);
		SeqOutput noTandemPairWriter(pairedNoTandemsOpts);
		SeqOutput noTandemSingleWriter(singletNoTandemsOuts);
		std::vector<uint32_t> readLengths;
		uint32_t pairedInReadCount = 0;
		if (pairedReader.ioOptions_.inExists()) {
			PairedRead pSeq;
			pairedReader.reOpenIn();
			while (pairedReader.readNextRead(pSeq)) {
				++pairedInReadCount;
				readLengths.emplace_back(len(pSeq.seqBase_));
				readLengths.emplace_back(len(pSeq.mateSeqBase_));
				auto firstExpansionPositions = expansionHelperInitial.getExpansionPositions(pSeq.seqBase_.seq_);
				auto secondExpansionPositions = expansionHelperInitial.getExpansionPositions(pSeq.mateSeqBase_.seq_);

				bool foundAnyTandemsFirst = !firstExpansionPositions.empty();
				bool foundAnyTandemsSecond = !secondExpansionPositions.empty();
				if (foundAnyTandemsFirst || foundAnyTandemsSecond) {
					tandemPairWriter.openWrite(pSeq);
					++seqsWithTandems;
				} else {
					noTandemPairWriter.openWrite(pSeq);
				}
//				error
//				if(foundAnyTandemsFirst && foundAnyTandemsSecond){
//					tandemPairWriter.openWrite(pSeq);
//					++seqsWithTandems;
//				} else if(foundAnyTandemsFirst){
//					++seqsWithTandems;
//					tandemSingleWriter.openWrite(pSeq.seqBase_);
//					noTandemSingleWriter.openWrite(pSeq.mateSeqBase_);
//				} else if(foundAnyTandemsSecond){
//					++seqsWithTandems;
//					tandemSingleWriter.openWrite(pSeq.mateSeqBase_);
//					noTandemSingleWriter.openWrite(pSeq.seqBase_);
//				} else {
//					noTandemPairWriter.openWrite(pSeq);
//				}
			}
		}
		uint32_t singleReadCount = 0;
		if (singleReader.ioOptions_.inExists()) {
			seqInfo seq;
			singleReader.reOpenIn();
			while (singleReader.readNextRead(seq)) {
				++singleReadCount;
				readLengths.emplace_back(len(seq));
				auto expansionPositions = expansionHelperInitial.getExpansionPositions(seq.seq_);
				bool foundAnyTandems = !expansionPositions.empty();
				if(foundAnyTandems){
					++seqsWithTandems;
					tandemSingleWriter.openWrite(seq);
				} else {
					noTandemSingleWriter.openWrite(seq);
				}
			}
		}
		klenLog["seqsWithTandems"] = seqsWithTandems;
		klenLog["pairedInReadCount"] = njh::json::toJson(pairedInReadCount);
		klenLog["singleReadCount"] = njh::json::toJson(singleReadCount);
	//	std::cout << "Line: " << __LINE__ << std::endl;
	//	for(const auto & count : graph.kCounts_){
	//		if(count.first.size() != klength){
	//			std::cout << count.first << "\t" << count.second << std::endl;
	//		}
	//	}
	//	std::cout << "TACACGACATACACGACATACACGACATACACGACATACA: " << graph.kCounts_["TACACGACATACACGACATACACGACATACACGACATACA"] << std::endl;
		bfs::path tandemR1Fnp;
		bfs::path tandemR2Fnp;
		bfs::path tandemSingleFnp;
		//close out writers so their buffers get flushed
		if(tandemPairWriter.outOpen()){
			tandemR1Fnp = tandemPairWriter.getPrimaryOutFnp();
			tandemR2Fnp = tandemPairWriter.getSecondaryOutFnp();
			tandemPairWriter.closeOut();
		}
		if(tandemSingleWriter.outOpen()){
			tandemSingleFnp = tandemSingleWriter.getPrimaryOutFnp();
			tandemSingleWriter.closeOut();
		}
		bfs::path notandemR1Fnp;
		bfs::path notandemR2Fnp;
		bfs::path notandemSingleFnp;
		if(noTandemPairWriter.outOpen()){
			notandemR1Fnp = noTandemPairWriter.getPrimaryOutFnp();
			notandemR2Fnp = noTandemPairWriter.getSecondaryOutFnp();
			noTandemPairWriter.closeOut();
		}else if(!bfs::exists(tandemR1Fnp) && !bfs::exists(tandemR2Fnp) ){
			//this check is needed in case all paired reads ended up with tandems
			notandemR1Fnp = pairedReader.ioOptions_.firstName_;
			notandemR2Fnp = pairedReader.ioOptions_.secondName_;
		}
		if(noTandemSingleWriter.outOpen()){
			notandemSingleFnp = noTandemSingleWriter.getPrimaryOutFnp();
			noTandemSingleWriter.closeOut();
		}else if(!bfs::exists(tandemSingleFnp) ){
			//this check is needed in case all single reads ended up with tandems
			notandemSingleFnp = singleReader.ioOptions_.firstName_;
		}
		auto medReadLeng = vectorMedianRef(readLengths);
		klenLog["readLengthAverage"] = njh::json::toJson(vectorMean(readLengths));
		klenLog["readLengthMedian"] = njh::json::toJson(medReadLeng);
		KmerPathwayGraph debugEstimatingGraph(std::min(extractionPars.estimatorKlen, *std::min_element(extractionPars.kmerLengths.begin(),extractionPars.kmerLengths.end() )));
		if(extractionPars.debug|| extractionPars.writeOutFinalDot_ || extractionPars.calcPercentUsedByKmerUsage){
		//{
			watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": Debug Getting Estimate Counts") );
			//debugEstimatingGraph.setOccurenceCutOff(kmerOccurenceCutOff);
			debugEstimatingGraph.debug_ = extractionPars.graphDebug_;
			debugEstimatingGraph.verbose_ = extractionPars.graphVerbose_;
//			auto usedSinglesFnp = singletOuts.getPriamryOutName();
//			auto usedPairedR1Fnp = pairedOpts.getPriamryOutName();
//			auto usedPairedR2Fnp = pairedOpts.getSecondaryOutName();
			auto usedSinglesFnp =  njh::files::make_path(finalCurrentDir, "filteredSingles.fastq");
			auto usedPairedR1Fnp = njh::files::make_path(finalCurrentDir, "filteredExtractedPairs_R1.fastq");
			auto usedPairedR2Fnp = njh::files::make_path(finalCurrentDir, "filteredExtractedPairs_R2.fastq");
			if(bfs::exists(usedPairedR1Fnp)){
				SeqInput reader(SeqIOOptions::genPairedIn(usedPairedR1Fnp, usedPairedR2Fnp));
				reader.openIn();
				PairedRead seq;
				while(reader.readNextRead(seq)){
					kmerInfo firstInfo(seq.seqBase_.seq_, debugEstimatingGraph.klen_, false);
					debugEstimatingGraph.increaseKCounts(seq.seqBase_.seq_);
					if(seq.mateSeqBase_.seq_.size() > debugEstimatingGraph.klen_){
						for (auto pos : iter::range(seq.mateSeqBase_.seq_.size() - debugEstimatingGraph.klen_ + 1)) {
							auto mateSeqKmer = seq.mateSeqBase_.seq_.substr(pos, debugEstimatingGraph.klen_);
							debugEstimatingGraph.kCounts_[mateSeqKmer] += 1;
//							if(!njh::in(mateSeqKmer, firstInfo.kmers_)){
//								debugEstimatingGraph.kCounts_[mateSeqKmer] += 1;
//							}
						}
					}
				}
			}
			if(bfs::exists(usedSinglesFnp)){
				SeqInput reader(SeqIOOptions::genFastqIn(usedSinglesFnp));
				reader.openIn();
				seqInfo seq;
				while(reader.readNextRead(seq)){
					debugEstimatingGraph.increaseKCounts(seq.seq_);
				}
			}
		}

		//StitchResForKLen klenRes(currentKLen, topCurrentDir);
		uint32_t totalInputBases = 0;
		uint32_t totalInputSequences = 0;
		std::shared_ptr<KmerPathwayGraph> firstGraph=std::make_shared<KmerPathwayGraph>(currentKLen);
		firstGraph->setOccurenceCutOff(extractionPars.kmerKOcurrenceCutOffs.front());
		firstGraph->debug_ = extractionPars.graphDebug_;
		firstGraph->verbose_ = extractionPars.graphVerbose_;
		firstGraph->numThreads_ = extractionPars.numThreads;
		firstGraph->bridgingCutOff_ = extractionPars.nodeBridgingReadPercCutOff_;
		//firstGraph->throwAwayConservedAddNodesDuringDisentaglement_ = extractionPars.throwAwayConservedAddNodesDuringDisentaglement;

		firstGraph->homopolymerIndelCollapseFreqMultiplier_ = extractionPars.homopolymerIndelCollapseFreqMultiplier_;
		firstGraph->allowableErrorForHPIndexCollapse_.oneBaseIndel_   = extractionPars.oneBaseIndelError_;
		firstGraph->allowableErrorForHPIndexCollapse_.twoBaseIndel_   = extractionPars.twoBaseIndelError_;
		firstGraph->allowableErrorForHPIndexCollapse_.largeBaseIndel_ = extractionPars.largeBaseIndelError_;


		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen,
		"- Adding Counts With No Tandems"));

		// populate nodes
		if (bfs::exists(notandemR1Fnp) && bfs::exists(notandemR2Fnp) ) {
			PairedRead pSeq;
			SeqInput noTandemPairReader(
					SeqIOOptions::genPairedIn(notandemR1Fnp, notandemR2Fnp));
			noTandemPairReader.openIn();
			while (noTandemPairReader.readNextRead(pSeq)) {
				++totalInputSequences;
				/**@todo handle the increase of counts to handle k-mers that appear in both pairs*/
				firstGraph->increaseKCounts(pSeq.seqBase_.seq_);
				firstGraph->increaseKCounts(pSeq.mateSeqBase_.seq_);
				totalInputBases += pSeq.seqBase_.seq_.size();
				totalInputBases += pSeq.mateSeqBase_.seq_.size();
			}
		}
		if (bfs::exists(notandemSingleFnp) ) {
			seqInfo seq;
			SeqInput noTandemSingleReader(SeqIOOptions::genFastqIn(notandemSingleFnp));
			noTandemSingleReader.openIn();
			while (noTandemSingleReader.readNextRead(seq)) {
				++totalInputSequences;
				firstGraph->increaseKCounts(seq.seq_);
				totalInputBases += seq.seq_.size();
			}
		}
		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen,
		"- Getting Node Counts"));

		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen,
		"- Adding Counts With Tandems"));
		if (bfs::exists(tandemR1Fnp)
		 || bfs::exists(tandemSingleFnp)) {

			if(bfs::exists(tandemR1Fnp)){
				SeqInput tandemPairReader(SeqIOOptions::genPairedIn(tandemR1Fnp, tandemR2Fnp ));
				PairedRead pSeq;
				tandemPairReader.openIn();
				while(tandemPairReader.readNextRead(pSeq)){
					/**@todo handle the increase of counts to handle k-mers that appear in both pairs*/
					//first mate
					{
						auto firstExpansionPositions = expansionHelperInitial.getExpansionPositions(pSeq.seqBase_.seq_);
						firstGraph->increaseKCountsAdjust(pSeq.seqBase_.seq_, firstExpansionPositions);
					}
					//second mate
					{
						auto secondExpansionPositions = expansionHelperInitial.getExpansionPositions(pSeq.mateSeqBase_.seq_);
						firstGraph->increaseKCountsAdjust(pSeq.mateSeqBase_.seq_, secondExpansionPositions);
					}
				}
			}
			//singles
			if(bfs::exists(tandemSingleFnp)){
				SeqInput singleReader(SeqIOOptions::genFastqIn(tandemSingleFnp ));
				seqInfo seq;
				singleReader.openIn();
				while(singleReader.readNextRead(seq)){
					//single
					{
						auto expansionPositions = expansionHelperInitial.getExpansionPositions(seq.seq_);
						firstGraph->increaseKCountsAdjust(seq.seq_, expansionPositions);
					}
				}
			}
		}
		double medianNodeCount = 0;
		double meanNodeCount = 0;
		double medianNodeCountAbove2 = 0;
		double meanNodeCountAbove2 = 0;
		double medianNodeCountAbove50Percent = 0;
		double meanNodeCountAbove50Percent = 0;

		double nodeCountForAutoCutOff = 0;
		auto used_kmerKOcurrenceCutOffs = extractionPars.kmerKOcurrenceCutOffs;
		{
			std::vector<uint32_t> nodeCounts;
			std::vector<uint32_t> nodeCountsAbove2;
			for(const auto & k : firstGraph->kCounts_){
				nodeCounts.emplace_back(k.second);
				if(k.second > 2){
					nodeCountsAbove2.emplace_back(k.second);
				}
			}
			medianNodeCount= vectorMedianRef(nodeCounts);
			meanNodeCount= vectorMean(nodeCounts);
			medianNodeCountAbove2= vectorMedianRef(nodeCountsAbove2);
			meanNodeCountAbove2= vectorMean(nodeCountsAbove2);
			klenLog["medianNodeCount"] = medianNodeCount;
			klenLog["meanNodeCount"] = meanNodeCount;
			klenLog["medianNodeCountAbove2"] = medianNodeCountAbove2;
			klenLog["meanNodeCountAbove2"] = meanNodeCountAbove2;
			nodeCountForAutoCutOff = std::max(meanNodeCountAbove2, medianNodeCountAbove2);

			if (meanNodeCountAbove2 >= 30 || medianNodeCountAbove2 > 30){
				uint32_t cutOff = std::round(std::max(medianNodeCountAbove2, meanNodeCountAbove2) * 0.5);
				std::vector<uint32_t> nodeCountsAbove50Percent;
				for (const auto & k : firstGraph->kCounts_) {
					nodeCounts.emplace_back(k.second);
					if (k.second > cutOff) {
						nodeCountsAbove50Percent.emplace_back(k.second);
					}
				}
				medianNodeCountAbove50Percent = vectorMedianRef(nodeCountsAbove50Percent);
				meanNodeCountAbove50Percent = vectorMean(nodeCountsAbove50Percent);

				klenLog["nodeCountCutOff"] = cutOff;
				klenLog["medianNodeCountAbove50Percent"] = medianNodeCountAbove50Percent;
				klenLog["meanNodeCountAbove50Percent"] = meanNodeCountAbove50Percent;
				nodeCountForAutoCutOff = std::max(meanNodeCountAbove50Percent, medianNodeCountAbove50Percent);


			}
			klenLog["nodeCountForAutoCutOff"] = nodeCountForAutoCutOff;
//			std::cout << "klen: " << currentKLen << std::endl;
//			std::cout << "medianNodeCount       : " << medianNodeCount << std::endl;
//			std::cout << "meanNodeCount         : " << meanNodeCount << std::endl;
//			std::cout << "medianNodeCountAbove2 : " << medianNodeCountAbove2 << std::endl;
//			std::cout << "meanNodeCountAbove2   : " << meanNodeCountAbove2 << std::endl;
			if(extractionPars.autoDetermineKCutsOnNodeCount){
				if(extractionPars.kmerKOcurrenceCutOffs.size() > 1 && extractionPars.autoDetermineKCutsOnNodeCount && (nodeCountForAutoCutOff > 100 || extractionPars.forceDetermineKCuts)){
					used_kmerKOcurrenceCutOffs.clear();
					for(const auto & perc : extractionPars.autoDetermineKCutsPercentages){
						used_kmerKOcurrenceCutOffs.emplace_back(std::round(perc * nodeCountForAutoCutOff));
					}
				}
			} else if (extractionPars.autoDetermineKCutsOnTotalBaseCount && !extractionPars.inputSeqs.empty()) {
				uint64_t totalBases = 0;
				for(const auto & rl : readLengths){
					totalBases += rl;
				}
				double totalInputLen = 0;
				for(const auto & input : extractionPars.inputSeqs){
					totalInputLen += len(input);
				}
				double perBaseCoverage = totalBases/totalInputLen;
				klenLog["perBaseCoverageByTotalBaseCount"] = perBaseCoverage;
				if(extractionPars.kmerKOcurrenceCutOffs.size() > 1 && (perBaseCoverage > 100 || extractionPars.forceDetermineKCuts)){
					used_kmerKOcurrenceCutOffs.clear();
					for(const auto & perc : extractionPars.autoDetermineKCutsPercentages){
						used_kmerKOcurrenceCutOffs.emplace_back(std::round(perc * perBaseCoverage));
					}
				}
			}
		}
		//std::cout << klenLog << std::endl;

		firstGraph->setOccurenceCutOff(used_kmerKOcurrenceCutOffs.front());
		klenLog["used_kmerKOcurrenceCutOffs"] = njh::json::toJson(used_kmerKOcurrenceCutOffs);


		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen,
		"- Populating Nodes"));
		klenLog["totalInputBases"] =  njh::json::toJson(totalInputBases);
		klenLog["totalInputSequences"] =  njh::json::toJson(totalInputSequences);
		//kcutLog["numOfPossibleKmers"] =  numOfPossibleKmers;
		firstGraph->populateNodesFromCounts();


		if(extractionPars.debug){
			OutOptions nodeCountsOutOpts(njh::files::make_path(currentKmerDirectory, "nodeCounts.tab.txt"));
			OutputStream nodeCountsOut(nodeCountsOutOpts);
			nodeCountsOut << "node\tcount" << std::endl;
			auto allKs = getVectorOfMapKeys(firstGraph->kCounts_);
			njh::sort(allKs);
			for(const auto & k : allKs){
				nodeCountsOut << k << "\t" << firstGraph->kCounts_[k] << std::endl;
			}
		}
		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen,
		"- Adding Edges With No Tandems"));
		if (bfs::exists(notandemR1Fnp) && bfs::exists(notandemR2Fnp) ) {
			PairedRead pSeq;
			SeqInput noTandemPairReader(
					SeqIOOptions::genPairedIn(notandemR1Fnp, notandemR2Fnp));
			noTandemPairReader.openIn();
			while (noTandemPairReader.readNextRead(pSeq)) {
				firstGraph->threadThroughSequence(pSeq);
			}
		}
		if (bfs::exists(notandemSingleFnp) ) {
			seqInfo seq;
			SeqInput noTandemSingleReader(SeqIOOptions::genFastqIn(notandemSingleFnp));
			noTandemSingleReader.openIn();
			while (noTandemSingleReader.readNextRead(seq)) {
				firstGraph->threadThroughSequence(seq);
			}
		}
		watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
				"klen_", currentKLen,
		"- Adding Edges With Tandems"));
		//pairs
		if(bfs::exists(tandemR1Fnp)){
			SeqInput tandemPairReader(SeqIOOptions::genPairedIn(tandemR1Fnp, tandemR2Fnp ));
			PairedRead pSeq;
			tandemPairReader.openIn();
			while(tandemPairReader.readNextRead(pSeq)){
				++totalInputSequences;
				//first mate
				auto firstExpansionPositions = expansionHelperInitial.getExpansionPositions(pSeq.seqBase_.seq_);
				auto secondExpansionPositions = expansionHelperInitial.getExpansionPositions(pSeq.mateSeqBase_.seq_);
				firstGraph->threadThroughSequenceAdjust(pSeq, firstExpansionPositions, secondExpansionPositions);
//					firstGraph->threadThroughSequenceAdjust(pSeq.seqBase_, firstExpansionPositions, pSeq.seqBase_.name_);
//					firstGraph->threadThroughSequenceAdjust(pSeq.mateSeqBase_,  secondExpansionPositions, pSeq.seqBase_.name_ + "_mate");
			}
		}
		//singles
		if(bfs::exists(tandemSingleFnp)){
			SeqInput singleReader(SeqIOOptions::genFastqIn(tandemSingleFnp ));
			seqInfo seq;
			singleReader.openIn();
			while(singleReader.readNextRead(seq)){
				//single
				{
					++totalInputSequences;
					auto expansionPositions = expansionHelperInitial.getExpansionPositions(seq.seq_);
					firstGraph->threadThroughSequenceAdjust(seq,expansionPositions, seq.name_);
				}
			}
		}
		uint32_t kCutIter = 0;
		std::vector<OptimizationReconResult> allCurrentKCutOffOptResults;
		//get max k occurence cut
		auto maxKCut = vectorMaximum(used_kmerKOcurrenceCutOffs);
		for(const auto kmerOccurenceCutOff : used_kmerKOcurrenceCutOffs){
			auto genPrefixForWatchLapNameKCut = [&watch,&currentKLen,&kmerOccurenceCutOff,
																			 &maxKCut](){
				return njh::pasteAsStr(
						njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
						"klen_", currentKLen,
				"-kcut_", njh::leftPadNumStr(kmerOccurenceCutOff, maxKCut));
			};
			bfs::path currentKCutOffDir = njh::files::makeDir(currentKmerDirectory, njh::files::MkdirPar(sampName + njh::pasteAsStr("_kcut-", leftPadNumStr(kmerOccurenceCutOff, maxKCut)) ));
			Json::Value & kcutOffLog = klenLog[sampName + njh::pasteAsStr("_kcut-", leftPadNumStr(kmerOccurenceCutOff, maxKCut))];
			kcutOffLog["kOccurenceCutOff"] =  kmerOccurenceCutOff;
			std::stringstream kCutOffLogErrorLog;
			uint32_t shortTipIter = 0;
			std::vector<OptimizationReconResult> allCurrentShortTipOptRunResults;
			//remove nodes for current occurrence cut off, ok to do for the first graph since all following graphs will have this cut off or above
			watch.startNewLap(
					njh::pasteAsStr(genPrefixForWatchLapNameKCut(),
			"- Pruning nodes"));
			firstGraph->setOccurenceCutOff(kmerOccurenceCutOff);
			//remove off edges and nodes
			firstGraph->turnOffNodesBelowCutOff();
			firstGraph->removeOffNodes();
			watch.startNewLap(
					njh::pasteAsStr(genPrefixForWatchLapNameKCut(),
			"- Pruning edges"));
			firstGraph->turnOffEdgesBelowCutOff();
			firstGraph->removeOffEdges();
			watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapNameKCut(),
			"- Checking node sanity"));
			firstGraph->edgeSanityCheckThrow();
			for(const auto shortTipNumberIter : extractionPars.shortTipNumbers){
				uint32_t shortTipNumber = shortTipNumberIter;
				if(!extractionPars.trimShortTips_){
					shortTipNumber = 0;
				}


				const auto maxShortTip = vectorMaximum(extractionPars.shortTipNumbers);
				auto genPrefixForWatchLapName = [&watch,&currentKLen,&kmerOccurenceCutOff,
																				 &maxKCut,&shortTipNumber,&maxShortTip](){
					return njh::pasteAsStr(
							njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": ",
							"klen_", currentKLen,
					"-kcut_", njh::leftPadNumStr(kmerOccurenceCutOff, maxKCut),
					"-shorttip_", njh::leftPadNumStr(shortTipNumber, maxShortTip) );
				};
				bfs::path currentShortTipNumberDir = njh::files::makeDir(
						currentKCutOffDir,
						njh::files::MkdirPar(
								sampName
										+ njh::pasteAsStr("_shortTip-",
												leftPadNumStr(shortTipNumber, maxShortTip))));
				Json::Value & shortTipLog = kcutOffLog[sampName + njh::pasteAsStr("_shortTip-", leftPadNumStr(shortTipNumber, maxShortTip))];
				shortTipLog["shortTipNumber"] =  shortTipNumber;

				uint32_t tipCutOff = std::numeric_limits<uint32_t>::max();
				if (!extractionPars.trimAllShortTips_) {

					if (nodeCountForAutoCutOff < 10 && !extractionPars.forceTipsByFreq_) {
						tipCutOff = std::numeric_limits<uint32_t>::max();
					} else {
						tipCutOff = std::max<uint32_t>(
								std::round(nodeCountForAutoCutOff * .30),
								kmerOccurenceCutOff + 2);
					}
				}

				shortTipLog["shortTipNumber-tipCutOff"] =  tipCutOff;

				//headlessTaillessLenCutOff
				auto headlessTaillessLenCutOff = extractionPars.headlessTailessLenCutOff;
				if(0 == headlessTaillessLenCutOff){
					headlessTaillessLenCutOff = medReadLeng;
				}else if(headlessTaillessLenCutOff < currentKLen){
					headlessTaillessLenCutOff = currentKLen + shortTipNumber;
				}
				shortTipLog["headlessTaillessLenCutOff"] =  headlessTaillessLenCutOff;
				std::stringstream shortTipLogErrorLog;
				OptimizationReconResult optRunRes(OptimizationReconResult::Params(currentKLen, kmerOccurenceCutOff, shortTipNumber),
						OptimizationReconResult::Dirs(currentKmerDirectory, currentKCutOffDir, currentShortTipNumberDir),
						kCutIter, shortTipIter);
				//std::cout << "currentKLen: " << currentKLen << ", kmerOccurenceCutOff: " << kmerOccurenceCutOff << std::endl;
				try {
					watch.startNewLap(
							njh::pasteAsStr(genPrefixForWatchLapName(), "- Copy Graph"));
					std::shared_ptr<KmerPathwayGraph> currentGraph;
					//copy graph so the entire setup doesn't have to be done again
					if(kmerOccurenceCutOff == maxKCut &&
							shortTipNumberIter == maxShortTip &&
							!extractionPars.writeOutFinalConnections_){
						//last time graph is going to be used so can just use the input since it won't be needed again
						currentGraph = firstGraph;
					}else{
						//copy the graph
						currentGraph = std::make_shared<KmerPathwayGraph>(firstGraph->copyGraph());
					}
					if (extractionPars.verbose) {
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(), "- Getting Stats"));
						table outNodeStates = currentGraph->getNodeStateCounts();
						outNodeStates.outPutContentOrganized(std::cout);
						table outTailCount = currentGraph->getTailCounts();
						outTailCount.outPutContentOrganized(std::cout);
					}
					watch.startNewLap(
							njh::pasteAsStr(genPrefixForWatchLapName(), "- Removing Orphan Nodes"));
					//remove orphan nodes
					currentGraph->removeHeadlessTaillessNodes();
					if(extractionPars.debug){
						std::ofstream outDot(njh::files::make_path(currentShortTipNumberDir, "initialGraph.dot").string());
						currentGraph->writeDot(outDot);
					}

					KmerGraphDebugWriter graphWriter(currentShortTipNumberDir,*currentGraph, debugEstimatingGraph);
					graphWriter.writeEdgeInfo_ = extractionPars.debugWriteEdgeInfo_;
					if(extractionPars.initialBreakingSinglePaths_){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
						"- Break Single Linked Paths Read Threading"));
						currentGraph->breakSingleLinkedPathsReadThreading();
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs("initialGraph_breakSingleLinkedPathsReadThreading");
						}
					}
					watch.startNewLap(
							njh::pasteAsStr(genPrefixForWatchLapName(), "- Collapse paths"));
					currentGraph->collapseSingleLinkedPaths(true);
//					if(kmerOccurenceCutOff > 2){
//						exit(1);
//					}
					if (extractionPars.debug) {
						graphWriter.writeOutDotsAndSeqs("collapsedInitialGraph");
					}
					if(extractionPars.removeHeadlessTaillessAlongTheWay){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
						"- Remove Headless Tailless Nodes"));
						currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff);
					}
					if(extractionPars.debug){
						graphWriter.writeOutDotsAndSeqs("collapsedInitialGraph-rmHeadTailLess");
					}
					watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(), "- Break Self Pointing Paths After Initial Collapse"));
					//break self pointing paths;
					if(currentGraph->breakSelfPointingPaths()){
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs("collapsedInitialGraph-rmHeadTailLess-breakSelfPointingPaths");
						}
						currentGraph->collapseSingleLinkedPaths();
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("collapsedInitialGraph-rmHeadTailLess-breakSelfPointingPaths-collapsed"));
						}
					}

					if(extractionPars.collapseOneBaseIndelsBeforeDisentanglement_){
						//collapse homopolymer nodes
						if(extractionPars.collapseOneBaseIndelsNodes_){
							watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(), "- Collapsing One Base Indel nodes immediately after initial collapse"));
							//if(currentGraph->collapseOneBaseIndelsNodes()){
							uint32_t collapseOneBaseIndelsNodesCount = 0;
							while(currentGraph->collapseOneBaseIndelsNodes(debugEstimatingGraph)){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseOneBaseIndelsNodes-",collapseOneBaseIndelsNodesCount ));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseOneBaseIndelsNodes-",collapseOneBaseIndelsNodesCount ));
								}
								++collapseOneBaseIndelsNodesCount;
							}
							uint32_t collapseOneBaseIndelsNodesCountComplex = 0;

							if(currentGraph->collapseOneBaseIndelsNodesComplex()){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseOneBaseIndelsNodesComplex-",collapseOneBaseIndelsNodesCountComplex));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseOneBaseIndelsNodesComplex-",collapseOneBaseIndelsNodesCountComplex, "-collapsed"));
								}
								++collapseOneBaseIndelsNodesCountComplex;
							}
						}
					}


//					if(extractionPars.collapseNodesWithAllowableError_){
//						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(), "- Collapsing nodes with allowable error immediately after initial collapse"));
//						//if(currentGraph->collapseOneBaseIndelsNodes()){
//						uint32_t collapseNodesWithErrorsCount = 0;
//						while(currentGraph->collapseBubbleNodesWithError(extractionPars.errorsToAllow_, extractionPars.collapseNodesWithAllowableErrorFreqMultiplier_)){
//							if(extractionPars.debug){
//								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseNodesWithAllowableError-",collapseNodesWithErrorsCount ));
//							}
//							currentGraph->collapseSingleLinkedPaths();
//							if(extractionPars.debug){
//								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseNodesWithAllowableError-",collapseNodesWithErrorsCount ));
//							}
//							++collapseNodesWithErrorsCount;
//						}
//					}
					uint32_t headlessTaillessBeforeDisentaglement = 0;
					uint32_t headlessTaillessBeforeDisentaglementBelowLen = 0;

					for(const auto & n : currentGraph->nodes_){
						if(n->on_){
							if(n->tailless() && n->headless()){
								++headlessTaillessBeforeDisentaglement;
								if(n->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
									++headlessTaillessBeforeDisentaglementBelowLen;
								}
							}
						}
					}

					if(extractionPars.startDisentanglementConservative_){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),"- Disentangle Internal Nodes Conservative"));
						uint32_t disentagleCount = 0;
						KmerPathwayGraph::disentangleInternalNodesPars disPars{true, shortTipNumber, tipCutOff, headlessTaillessLenCutOff, extractionPars};
						while(currentGraph->disentangleInternalNodes(disPars)){
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-", disentagleCount));
							}
							//currentGraph->breakSingleHeadSingleTailNodesLowCoverage();
//							if(extractionPars.trimShortTips_){
//								bool removedTips = currentGraph->removeShortTips_after_disentangleInternalNodes(shortTipNumber, tipCutOff);
//								if(removedTips){
//									if(extractionPars.debug){
//										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-", disentagleCount, "-removeShortTips"));
//									}
//								}
//							}
							currentGraph->collapseSingleLinkedPaths();
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-", disentagleCount, "-collapsed"));
							}
							//break self pointing paths;
							if(currentGraph->breakSelfPointingPaths()){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-", disentagleCount, "-collapsed-breakSelfPointingPaths"));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-", disentagleCount, "-collapsed-breakSelfPointingPaths-collapsed"));
								}
							}

							if(extractionPars.collapseOneBaseIndelsNodes_){
								watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
							"- Collapsing One Base Indel nodes while Split Internals Conservative"));
								//if(currentGraph->collapseOneBaseIndelsNodes()){
								uint32_t collapseOneBaseIndelsNodesCount = 0;
								while(currentGraph->collapseOneBaseIndelsNodes(debugEstimatingGraph)){
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-collapseOneBaseIndelsNodes-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCount ));
									}
									currentGraph->collapseSingleLinkedPaths();
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-collapseOneBaseIndelsNodes-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCount , "-collapsed"));
									}
									++collapseOneBaseIndelsNodesCount;
								}

								uint32_t collapseOneBaseIndelsNodesCountComplex = 0;
								if(currentGraph->collapseOneBaseIndelsNodesComplex()){
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-collapseOneBaseIndelsNodesComplex-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCountComplex));
									}
									currentGraph->collapseSingleLinkedPaths();
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-Conservative-collapseOneBaseIndelsNodesComplex-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCountComplex, "-collapsed"));
									}
									++collapseOneBaseIndelsNodesCountComplex;
								}
							}


							++disentagleCount;
							if(extractionPars.debug){
								std::cout << __FILE__ << " " << __LINE__ << std::endl;
								std::cout << "conservative disentagleCount: " << disentagleCount << std::endl;
							}
						}
					}
					watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),"- Disentangle Internal Nodes"));
					uint32_t disentagleCount = 0;
					KmerPathwayGraph::disentangleInternalNodesPars disPars{false, shortTipNumber, tipCutOff, headlessTaillessLenCutOff, extractionPars};
					while(currentGraph->disentangleInternalNodes(disPars)){
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-", disentagleCount));
						}
						//currentGraph->breakSingleHeadSingleTailNodesLowCoverage();
						if(extractionPars.trimShortTips_){
							bool removedTips = currentGraph->removeShortTips_after_disentangleInternalNodes(shortTipNumber, tipCutOff);
							if(removedTips){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-", disentagleCount, "-removeShortTips"));
								}
							}
						}
						currentGraph->collapseSingleLinkedPaths();
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-", disentagleCount, "-collapsed"));
						}
						//break self pointing paths;
						if(currentGraph->breakSelfPointingPaths()){
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-", disentagleCount, "-collapsed-breakSelfPointingPaths"));
							}
							currentGraph->collapseSingleLinkedPaths();
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-", disentagleCount, "-collapsed-breakSelfPointingPaths-collapsed"));
							}
						}
						if(extractionPars.collapseOneBaseIndelsNodes_){
							watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
						"- Collapsing One Base Indel nodes while Split Internals Conservative"));
							//if(currentGraph->collapseOneBaseIndelsNodes()){
							uint32_t collapseOneBaseIndelsNodesCount = 0;
							while(currentGraph->collapseOneBaseIndelsNodes(debugEstimatingGraph)){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-collapseOneBaseIndelsNodes-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCount ));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-collapseOneBaseIndelsNodes-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCount , "-collapsed"));
								}
								++collapseOneBaseIndelsNodesCount;
							}

							uint32_t collapseOneBaseIndelsNodesCountComplex = 0;
							if(currentGraph->collapseOneBaseIndelsNodesComplex()){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-collapseOneBaseIndelsNodesComplex-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCountComplex));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-collapseOneBaseIndelsNodesComplex-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCountComplex, "-collapsed"));
								}
								++collapseOneBaseIndelsNodesCountComplex;
							}
						}

						++disentagleCount;
						if(extractionPars.debug){
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << "disentagleCount: " << disentagleCount << std::endl;
						}
					}

					if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){
						currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff);
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-", disentagleCount, "-collapsed-rmHeadTailless"));
						}
					}


					uint32_t headlessTaillessAfterDisentaglement = 0;
					uint32_t headlessTaillessAfterDisentaglementBelowLen = 0;

					for(const auto & n : currentGraph->nodes_){
						if(n->on_){
							if(n->tailless() && n->headless()){
								++headlessTaillessAfterDisentaglement;
								if(n->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
									++headlessTaillessAfterDisentaglementBelowLen;
								}
							}
						}
					}
					uint32_t headlessTaillessDifference = 0;
					if(headlessTaillessAfterDisentaglement > headlessTaillessBeforeDisentaglement){
						headlessTaillessDifference = headlessTaillessAfterDisentaglement - headlessTaillessBeforeDisentaglement;
					}

					uint32_t headlessTaillessDifferenceBelowLen = 0;
					if(headlessTaillessAfterDisentaglementBelowLen > headlessTaillessBeforeDisentaglementBelowLen){
						headlessTaillessDifferenceBelowLen = headlessTaillessAfterDisentaglementBelowLen - headlessTaillessBeforeDisentaglementBelowLen;
					}

					if(extractionPars.collapseOneBaseIndelsNodes_){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Collapsing One Base Indel nodes"));
						//if(currentGraph->collapseOneBaseIndelsNodes()){
						uint32_t collapseOneBaseIndelsNodesCount = 0;
						while(currentGraph->collapseOneBaseIndelsNodes(debugEstimatingGraph)){
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-collapseOneBaseIndelsNodes-",collapseOneBaseIndelsNodesCount ));
							}
							currentGraph->collapseSingleLinkedPaths();
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-collapseOneBaseIndelsNodes-",collapseOneBaseIndelsNodesCount ));
							}
							++collapseOneBaseIndelsNodesCount;
						}
						uint32_t collapseOneBaseIndelsNodesCountComplex = 0;

						if(currentGraph->collapseOneBaseIndelsNodesComplex()){
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-collapseOneBaseIndelsNodesComplex-",collapseOneBaseIndelsNodesCountComplex));
							}
							currentGraph->collapseSingleLinkedPaths();
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-collapseOneBaseIndelsNodesComplex-",collapseOneBaseIndelsNodesCountComplex, "-collapsed"));
							}
							++collapseOneBaseIndelsNodesCountComplex;
						}
					}

					if(extractionPars.collapseNodesWithAllowableError_){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(), "- Collapsing nodes with allowable error"));
						//if(currentGraph->collapseOneBaseIndelsNodes()){
						uint32_t collapseNodesWithErrorsCount = 0;
						while(currentGraph->collapseBubbleNodesWithError(extractionPars.errorsToAllow_, extractionPars.collapseNodesWithAllowableErrorFreqMultiplier_)){
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseNodesWithAllowableError-",collapseNodesWithErrorsCount ));
							}
							currentGraph->collapseSingleLinkedPaths();
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseNodesWithAllowableError-",collapseNodesWithErrorsCount ));
							}
							++collapseNodesWithErrorsCount;
						}
					}

					if(extractionPars.trimShortTips_){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Removing Short tips"));
						uint32_t removeShortTipRunNumber = 0;
						while(currentGraph->removeShortTips(shortTipNumber, tipCutOff)){
							if(extractionPars.debug){
								std::cout << "Removing short tips run: " << removeShortTipRunNumber << std::endl;
							}
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-removedShortTips", "-", removeShortTipRunNumber));
							}
							//currentGraph->breakSingleHeadSingleTailNodesLowCoverage();
							currentGraph->collapseSingleLinkedPaths();
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-removedShortTips", "-", removeShortTipRunNumber, "-collapsed"));
							}
							//break self pointing paths;
							if(currentGraph->breakSelfPointingPaths()){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-removedShortTips", "-", removeShortTipRunNumber, "-collapsed-breakSelfPointingPaths"));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-removedShortTips", "-", removeShortTipRunNumber, "-collapsed-breakSelfPointingPaths-collapsed"));
								}
							}
							++removeShortTipRunNumber;
						}

						if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){
							currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff);
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-removedShortTips", "-", removeShortTipRunNumber, "-collapsed-rmHeadTailless"));
							}
						}
						{
							uint32_t headlessTaillessBeforeDisentaglementAfterTip = 0;
							uint32_t headlessTaillessBeforeDisentaglementBelowLenAfterTip = 0;

							for(const auto & n : currentGraph->nodes_){
								if(n->on_){
									if(n->tailless() && n->headless()){
										++headlessTaillessBeforeDisentaglementAfterTip;
										if(n->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
											++headlessTaillessBeforeDisentaglementBelowLenAfterTip;
										}
									}
								}
							}
							if(extractionPars.startDisentanglementConservative_){
								watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),"- Disentangle Internal Nodes Conservative After tip removal"));
								uint32_t disentagleCount = 0;
								KmerPathwayGraph::disentangleInternalNodesPars disPars{true, shortTipNumber, tipCutOff, headlessTaillessLenCutOff, extractionPars};
								while(currentGraph->disentangleInternalNodes(disPars)){
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-", disentagleCount));
									}
									//currentGraph->breakSingleHeadSingleTailNodesLowCoverage();
//									if(extractionPars.trimShortTips_){
//										bool removedTips = currentGraph->removeShortTips_after_disentangleInternalNodes(shortTipNumber, tipCutOff);
//										if(removedTips){
//											if(extractionPars.debug){
//												graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-", disentagleCount, "-removeShortTips"));
//											}
//										}
//									}
									currentGraph->collapseSingleLinkedPaths();
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-", disentagleCount, "-collapsed"));
									}
									//break self pointing paths;
									if(currentGraph->breakSelfPointingPaths()){
										if(extractionPars.debug){
											graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-", disentagleCount, "-collapsed-breakSelfPointingPaths"));
										}
										currentGraph->collapseSingleLinkedPaths();
										if(extractionPars.debug){
											graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-", disentagleCount, "-collapsed-breakSelfPointingPaths-collapsed"));
										}
									}
									if(extractionPars.collapseOneBaseIndelsNodes_){
										watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
									"- Collapsing One Base Indel nodes while Split Internals Conservative"));
										//if(currentGraph->collapseOneBaseIndelsNodes()){
										uint32_t collapseOneBaseIndelsNodesCount = 0;
										while(currentGraph->collapseOneBaseIndelsNodes(debugEstimatingGraph)){
											if(extractionPars.debug){
												graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-collapseOneBaseIndelsNodes-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCount ));
											}
											currentGraph->collapseSingleLinkedPaths();
											if(extractionPars.debug){
												graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-collapseOneBaseIndelsNodes-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCount , "-collapsed"));
											}
											++collapseOneBaseIndelsNodesCount;
										}

										uint32_t collapseOneBaseIndelsNodesCountComplex = 0;
										if(currentGraph->collapseOneBaseIndelsNodesComplex()){
											if(extractionPars.debug){
												graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-collapseOneBaseIndelsNodesComplex-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCountComplex));
											}
											currentGraph->collapseSingleLinkedPaths();
											if(extractionPars.debug){
												graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-Conservative-collapseOneBaseIndelsNodesComplex-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCountComplex, "-collapsed"));
											}
											++collapseOneBaseIndelsNodesCountComplex;
										}
									}
									++disentagleCount;
									if(extractionPars.debug){
										std::cout << __FILE__ << " " << __LINE__ << std::endl;
										std::cout << "after tip removal conservative disentagleCount: " << disentagleCount << std::endl;
									}
								}
							}

							watch.startNewLap(
									njh::pasteAsStr(genPrefixForWatchLapName(),
					    "- Connecting easy adjacent bubbles after tip removal"));
							KmerPathwayGraph::disentangleInternalNodesPars disPars{false, shortTipNumber, tipCutOff, headlessTaillessLenCutOff, extractionPars};
							while(currentGraph->disentangleInternalNodes(disPars)){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-", disentagleCount));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-", disentagleCount,"-collapsed"));
								}
								//break self pointing paths;
								if(currentGraph->breakSelfPointingPaths()){
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-", disentagleCount,"-collapsed-breakSelfPointingPaths"));
									}
									currentGraph->collapseSingleLinkedPaths();
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval", "-", disentagleCount, "-collapsed-breakSelfPointingPaths-collapsed"));
									}
								}
								if(extractionPars.collapseOneBaseIndelsNodes_){
									watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
								"- Collapsing One Base Indel nodes while Split Internals Conservative"));
									//if(currentGraph->collapseOneBaseIndelsNodes()){
									uint32_t collapseOneBaseIndelsNodesCount = 0;
									while(currentGraph->collapseOneBaseIndelsNodes(debugEstimatingGraph)){
										if(extractionPars.debug){
											graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodes-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCount ));
										}
										currentGraph->collapseSingleLinkedPaths();
										if(extractionPars.debug){
											graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodes-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCount , "-collapsed"));
										}
										++collapseOneBaseIndelsNodesCount;
									}

									uint32_t collapseOneBaseIndelsNodesCountComplex = 0;
									if(currentGraph->collapseOneBaseIndelsNodesComplex()){
										if(extractionPars.debug){
											graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodesComplex-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCountComplex));
										}
										currentGraph->collapseSingleLinkedPaths();
										if(extractionPars.debug){
											graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodesComplex-", disentagleCount, "-collapsed-",collapseOneBaseIndelsNodesCountComplex, "-collapsed"));
										}
										++collapseOneBaseIndelsNodesCountComplex;
									}
								}
								++disentagleCount;
								if(extractionPars.debug){
									std::cout << __FILE__ << " " << __LINE__ << std::endl;
									std::cout << "after tip removal disentagleCount: " << disentagleCount << std::endl;
								}
							}

							if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){
								currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff);
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-rmHeadTailless"));
								}
							}
							uint32_t headlessTaillessAfterDisentaglementAfterTip = 0;
							uint32_t headlessTaillessAfterDisentaglementBelowLenAfterTip = 0;

							for(const auto & n : currentGraph->nodes_){
								if(n->on_){
									if(n->tailless() && n->headless()){
										++headlessTaillessAfterDisentaglementAfterTip;
										if(n->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
											++headlessTaillessAfterDisentaglementBelowLenAfterTip;
										}
									}
								}
							}
							if(headlessTaillessAfterDisentaglementBelowLenAfterTip > headlessTaillessBeforeDisentaglementAfterTip){
								headlessTaillessDifference = headlessTaillessAfterDisentaglementBelowLenAfterTip - headlessTaillessBeforeDisentaglementAfterTip;
							}

							if(headlessTaillessAfterDisentaglementBelowLenAfterTip > headlessTaillessBeforeDisentaglementBelowLenAfterTip){
								headlessTaillessDifferenceBelowLen = headlessTaillessAfterDisentaglementBelowLenAfterTip - headlessTaillessBeforeDisentaglementBelowLenAfterTip;
							}
						}

						if(extractionPars.collapseOneBaseIndelsNodes_){
							watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
						"- Collapsing One Base Indel nodes after short tip removal"));
							//if(currentGraph->collapseOneBaseIndelsNodes()){
							uint32_t collapseOneBaseIndelsNodesCount = 0;
							while(currentGraph->collapseOneBaseIndelsNodes(debugEstimatingGraph)){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodes-",collapseOneBaseIndelsNodesCount ));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodes-",collapseOneBaseIndelsNodesCount , "-collapsed"));
								}
								++collapseOneBaseIndelsNodesCount;
							}

							uint32_t collapseOneBaseIndelsNodesCountComplex = 0;
							if(currentGraph->collapseOneBaseIndelsNodesComplex()){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodesComplex-",collapseOneBaseIndelsNodesCountComplex));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodesComplex-",collapseOneBaseIndelsNodesCountComplex, "-collapsed"));
								}
								++collapseOneBaseIndelsNodesCountComplex;
							}
						}
						if(extractionPars.collapseNodesWithAllowableError_){
							watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(), "- Collapsing nodes with allowable error after short tip removal"));
							//if(currentGraph->collapseOneBaseIndelsNodes()){
							uint32_t collapseNodesWithErrorsCount = 0;
							while(currentGraph->collapseBubbleNodesWithError(extractionPars.errorsToAllow_, extractionPars.collapseNodesWithAllowableErrorFreqMultiplier_)){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseNodesWithAllowableError-",collapseNodesWithErrorsCount ));
								}
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-collapsedInitialGraph-collapseNodesWithAllowableError-",collapseNodesWithErrorsCount ));
								}
								++collapseNodesWithErrorsCount;
							}
						}
						if(extractionPars.collapseOneBaseIndelsNodes_ || extractionPars.collapseNodesWithAllowableError_){
							if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){
								if(currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff)){
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-afterCollapsingErrorNodes-rmHeadTailless"));
									}
								}
							}
						}
					}else{
						if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){
							if(currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff)){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-rmHeadTailless"));
								}
							}
						}
					}


					{
						//one final disentanglement
						{
							uint32_t headlessTaillessBeforeDisentaglementFinal = 0;
							uint32_t headlessTaillessBeforeDisentaglementBelowLenFinal = 0;

							for(const auto & n : currentGraph->nodes_){
								if(n->on_){
									if(n->tailless() && n->headless()){
										++headlessTaillessBeforeDisentaglementFinal;
										if(n->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
											++headlessTaillessBeforeDisentaglementBelowLenFinal;
										}
									}
								}
							}
							if(extractionPars.startDisentanglementConservative_){
								watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),"- Disentangle Internal Nodes Conservative After tip removal"));
								uint32_t disentagleCount = 0;
								KmerPathwayGraph::disentangleInternalNodesPars disPars{true, shortTipNumber, tipCutOff, headlessTaillessLenCutOff, extractionPars};
								while(currentGraph->disentangleInternalNodes(disPars)){
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-Conservative-", disentagleCount));
									}
									//currentGraph->breakSingleHeadSingleTailNodesLowCoverage();
//									if(extractionPars.trimShortTips_){
//										bool removedTips = currentGraph->removeShortTips_after_disentangleInternalNodes(shortTipNumber, tipCutOff);
//										if(removedTips){
//											if(extractionPars.debug){
//												graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-Conservative-", disentagleCount, "-removeShortTips"));
//											}
//										}
//									}
									currentGraph->collapseSingleLinkedPaths();
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-Conservative-", disentagleCount, "-collapsed"));
									}
									//break self pointing paths;
									if(currentGraph->breakSelfPointingPaths()){
										if(extractionPars.debug){
											graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-Conservative-", disentagleCount, "-collapsed-breakSelfPointingPaths"));
										}
										currentGraph->collapseSingleLinkedPaths();
										if(extractionPars.debug){
											graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-Conservative-", disentagleCount, "-collapsed-breakSelfPointingPaths-collapsed"));
										}
									}
									++disentagleCount;
									if(extractionPars.debug){
										std::cout << __FILE__ << " " << __LINE__ << std::endl;
										std::cout << "after tip removal conservative disentagleCount: " << disentagleCount << std::endl;
									}
								}
							}

							watch.startNewLap(
									njh::pasteAsStr(genPrefixForWatchLapName(),
					    "- Connecting easy adjacent bubbles after tip removal"));
							KmerPathwayGraph::disentangleInternalNodesPars disPars{false, shortTipNumber, tipCutOff, headlessTaillessLenCutOff, extractionPars};
							while(currentGraph->disentangleInternalNodes(disPars)){
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-", disentagleCount));
								}
								//currentGraph->breakSingleHeadSingleTailNodesLowCoverage();
								currentGraph->collapseSingleLinkedPaths();
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-", disentagleCount,"-collapsed"));
								}
								//break self pointing paths;
								if(currentGraph->breakSelfPointingPaths()){
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-", disentagleCount,"-collapsed-breakSelfPointingPaths"));
									}
									currentGraph->collapseSingleLinkedPaths();
									if(extractionPars.debug){
										graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-", disentagleCount, "-collapsed-breakSelfPointingPaths-collapsed"));
									}
								}
								++disentagleCount;
								if(extractionPars.debug){
									std::cout << __FILE__ << " " << __LINE__ << std::endl;
									std::cout << "after tip removal disentagleCount: " << disentagleCount << std::endl;
								}
							}

							if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){
								currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff);
								if(extractionPars.debug){
									graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-", disentagleCount,"-collapsed-rmHeadTailless"));
								}
							}
							uint32_t headlessTaillessAfterDisentaglementFinal = 0;
							uint32_t headlessTaillessAfterDisentaglementBelowLenFinal = 0;

							for(const auto & n : currentGraph->nodes_){
								if(n->on_){
									if(n->tailless() && n->headless()){
										++headlessTaillessAfterDisentaglementFinal;
										if(n->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
											++headlessTaillessAfterDisentaglementBelowLenFinal;
										}
									}
								}
							}
							if(headlessTaillessAfterDisentaglementBelowLenFinal > headlessTaillessBeforeDisentaglementFinal){
								headlessTaillessDifference = headlessTaillessAfterDisentaglementBelowLenFinal - headlessTaillessBeforeDisentaglementFinal;
							}

							if(headlessTaillessAfterDisentaglementBelowLenFinal > headlessTaillessBeforeDisentaglementBelowLenFinal){
								headlessTaillessDifferenceBelowLen = headlessTaillessAfterDisentaglementBelowLenFinal - headlessTaillessBeforeDisentaglementBelowLenFinal;
							}
						}
					}

					//count tailless and countless before outputting sequence
					//StitchResForKCut kcutRes(kmerOccurenceCutOff, currentKCutDir, kCutIter);
					std::unordered_set<std::string> finalNodeReadCounts;
					//uint64_t totalBasesUsed = 0;
					//determine optimization results;
					for(const auto & node : currentGraph->nodes_){
						if(node->on_){
							for(const auto & nameIdx : node->inReadNamesIdx_){
								std::string rName = currentGraph->readNames_[nameIdx];
								if(std::string::npos != rName.find("-multipleOccurenceKmer-")){
									rName = rName.substr(0, rName.find("-multipleOccurenceKmer-"));
								}
								finalNodeReadCounts.emplace(njh::replaceString(rName, "_mate", ""));
							}
							if(!extractionPars.throwAwayConservedAddNodesDuringDisentaglement){
								bool skipCurrentNode = false;
								for(const auto & otherNode : currentGraph->nodes_){
									if(otherNode->on_){
										if(otherNode->k_.size() > node->k_.size() && node->tailless() && node->headless() && std::string::npos != otherNode->k_.find(node->k_)){
											//we found this node perfectly within another node and this node is headless and tailless
											//this would happen if we kept a conserved piece and didn't throw it away, don't want this to count towards optimization in this case
											skipCurrentNode = true;
											break;
										}
									}
								}
								if(skipCurrentNode){
									continue;
								}
							}
							if(node->tailless()){
								++optRunRes.taillessCount_;
								if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
									++optRunRes.taillessCountBelowLen_;
								}
							}
							if(node->headless()){
								++optRunRes.headlessCount_;
								if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
									++optRunRes.headlessCountBelowLen_;
								}
							}
							if(node->tailless() || node->headless()){
								++optRunRes.optimalCount_;
								if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
									++optRunRes.optimalCountBelowLen_;
								}
							}
							if(!node->tailless() && !node->headless()){
								if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
									++optRunRes.internalNodesCountBelowLen_;
								}
							}
						}
					}
					if(extractionPars.throwAwayConservedAddNodesDuringDisentaglement && extractionPars.adjustForHeadLessAddedAfterDisentanglement){
						if(extractionPars.debug){
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << "headlessTaillessDifference:" << headlessTaillessDifference << std::endl;
							std::cout << "headlessTaillessDifferenceBelowLen:" << headlessTaillessDifferenceBelowLen << std::endl;
						}
						if(headlessTaillessDifference > 0 && optRunRes.optimalCount_ != headlessTaillessDifference){
							if(optRunRes.taillessCount_ > headlessTaillessDifference){
								optRunRes.taillessCount_ -=	headlessTaillessDifference;
							}else{
								optRunRes.taillessCount_ = 0; //i don't think this should be possible
							}
							if(optRunRes.headlessCount_ > headlessTaillessDifference){
								optRunRes.headlessCount_ -=	headlessTaillessDifference;
							}else{
								optRunRes.headlessCount_ = 0; //i don't think this should be possible
							}
							if(optRunRes.optimalCount_ > headlessTaillessDifference){
								optRunRes.optimalCount_ -=	headlessTaillessDifference;
							}else{
								optRunRes.optimalCount_ = 0; //i don't think this should be possible
							}
						}
						if(headlessTaillessDifferenceBelowLen > 0 && optRunRes.optimalCountBelowLen_ != headlessTaillessDifferenceBelowLen){
							if(optRunRes.taillessCountBelowLen_ > headlessTaillessDifferenceBelowLen){
								optRunRes.taillessCountBelowLen_ -=	headlessTaillessDifferenceBelowLen;
							}else{
								optRunRes.taillessCountBelowLen_ = 0; //i don't think this should be possible
							}
							if(optRunRes.headlessCountBelowLen_ > headlessTaillessDifferenceBelowLen){
								optRunRes.headlessCountBelowLen_ -=	headlessTaillessDifferenceBelowLen;
							}else{
								optRunRes.headlessCountBelowLen_ = 0; //i don't think this should be possible
							}
							if(optRunRes.optimalCountBelowLen_ > headlessTaillessDifferenceBelowLen){
								optRunRes.optimalCountBelowLen_ -=	headlessTaillessDifferenceBelowLen;
							}else{
								optRunRes.optimalCountBelowLen_ = 0; //i don't think this should be possible
							}
						}
					}
					//error add taking the min beteween splitting ends once and the previous
					//exit(1);
//					if( 2 == optRunRes.runParams_.kcut_ &&
//						 25 == optRunRes.runParams_.klen_&&
//							2 == optRunRes.runParams_.shortNumber_){
//						//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//						std::cout << optRunRes.toJson() << std::endl;
//					}

					if (extractionPars.splitEndsOnce_) {
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Splitting End Nodes Once"));
						currentGraph->splitEndNodes(extractionPars.maxSplitEndNodeSize_);
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitEndsOnce"));
						}
						currentGraph->collapseSingleLinkedPaths();
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitEndsOnce-collapsed"));
						}
						if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){

							currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff);
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitEndsOnce-collapsed-rmHeadTailless"));
							}
						}


						OptimizationReconResult optRunResAfterOneEndSplit(OptimizationReconResult::Params(currentKLen, kmerOccurenceCutOff, shortTipNumber),
								OptimizationReconResult::Dirs(currentKmerDirectory, currentKCutOffDir, currentShortTipNumberDir),
								kCutIter, shortTipIter);

						//determine optimization results;
						for(const auto & node : currentGraph->nodes_){
							if(node->on_){
								if(!extractionPars.throwAwayConservedAddNodesDuringDisentaglement){
									bool skipCurrentNode = false;
									for(const auto & otherNode : currentGraph->nodes_){
										if(otherNode->on_){
											if(otherNode->k_.size() > node->k_.size() && node->tailless() && node->headless() && std::string::npos != otherNode->k_.find(node->k_)){
												//we found this node perfectly within another node and this node is headless and tailless
												//this would happen if we kept a conserved piece and didn't throw it away, don't want this to count towards optimization in this case
												skipCurrentNode = true;
												break;
											}
										}
									}
									if(skipCurrentNode){
										continue;
									}
								}
								if(node->tailless()){
									++optRunResAfterOneEndSplit.taillessCount_;
									if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
										++optRunResAfterOneEndSplit.taillessCountBelowLen_;
									}
								}
								if(node->headless()){
									++optRunResAfterOneEndSplit.headlessCount_;
									if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
										++optRunResAfterOneEndSplit.headlessCountBelowLen_;
									}
								}
								if(node->tailless() || node->headless()){
									++optRunResAfterOneEndSplit.optimalCount_;
									if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
										++optRunResAfterOneEndSplit.optimalCountBelowLen_;
									}
								}
								if(!node->tailless() && !node->headless()){
									if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
										++optRunResAfterOneEndSplit.internalNodesCountBelowLen_;
									}
								}
							}
						}

						if(extractionPars.throwAwayConservedAddNodesDuringDisentaglement && extractionPars.adjustForHeadLessAddedAfterDisentanglement){
							if(extractionPars.debug){
								std::cout << __FILE__ << " " << __LINE__ << std::endl;
								std::cout << "headlessTaillessDifference:" << headlessTaillessDifference << std::endl;
								std::cout << "headlessTaillessDifferenceBelowLen:" << headlessTaillessDifferenceBelowLen << std::endl;
							}
							if(headlessTaillessDifference > 0 && optRunResAfterOneEndSplit.optimalCount_ != headlessTaillessDifference){
								if(optRunResAfterOneEndSplit.taillessCount_ > headlessTaillessDifference){
									optRunResAfterOneEndSplit.taillessCount_ -=	headlessTaillessDifference;
								}else{
									optRunResAfterOneEndSplit.taillessCount_ = 0; //i don't think this should be possible
								}
								if(optRunResAfterOneEndSplit.headlessCount_ > headlessTaillessDifference){
									optRunResAfterOneEndSplit.headlessCount_ -=	headlessTaillessDifference;
								}else{
									optRunResAfterOneEndSplit.headlessCount_ = 0; //i don't think this should be possible
								}
								if(optRunResAfterOneEndSplit.optimalCount_ > headlessTaillessDifference){
									optRunResAfterOneEndSplit.optimalCount_ -=	headlessTaillessDifference;
								}else{
									optRunResAfterOneEndSplit.optimalCount_ = 0; //i don't think this should be possible
								}
							}
							if(headlessTaillessDifferenceBelowLen > 0 && optRunResAfterOneEndSplit.optimalCountBelowLen_ != headlessTaillessDifferenceBelowLen){
								if(optRunResAfterOneEndSplit.taillessCountBelowLen_ > headlessTaillessDifferenceBelowLen){
									optRunResAfterOneEndSplit.taillessCountBelowLen_ -=	headlessTaillessDifferenceBelowLen;
								}else{
									optRunResAfterOneEndSplit.taillessCountBelowLen_ = 0; //i don't think this should be possible
								}
								if(optRunResAfterOneEndSplit.headlessCountBelowLen_ > headlessTaillessDifferenceBelowLen){
									optRunResAfterOneEndSplit.headlessCountBelowLen_ -=	headlessTaillessDifferenceBelowLen;
								}else{
									optRunResAfterOneEndSplit.headlessCountBelowLen_ = 0; //i don't think this should be possible
								}
								if(optRunResAfterOneEndSplit.optimalCountBelowLen_ > headlessTaillessDifferenceBelowLen){
									optRunResAfterOneEndSplit.optimalCountBelowLen_ -=	headlessTaillessDifferenceBelowLen;
								}else{
									optRunResAfterOneEndSplit.optimalCountBelowLen_ = 0; //i don't think this should be possible
								}
							}
						}
						if(extractionPars.useOptAfterOneEndSplit_ && optRunResAfterOneEndSplit.optimalCount_ < optRunRes.optimalCount_){
							optRunRes = optRunResAfterOneEndSplit;
						}
					}


					if (extractionPars.splitEnds_) {
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Splitting End Nodes"));
						uint32_t splitEndsCount = 0;
						while (currentGraph->splitEndNodes(extractionPars.maxSplitEndNodeSize_)) {
							if (currentGraph->debug_) {
								double failCount = 0;
								std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
								for (const auto & n : currentGraph->nodes_) {
									if (1 == n->tailCount() && 1 == n->headCount()) {
										nodesToProcess.push_back(n);
									}
								}
								for (const auto & n : nodesToProcess) {
									if (n->headCount() == 1 && n->tailCount() == 1) {
										//next should only have a head count of 1 for it to be here
										std::set<uint32_t> readsEntering(
												n->getFirstOnHeadEdge()->inReadNamesIdx_.begin(),
												n->getFirstOnHeadEdge()->inReadNamesIdx_.end());
										//and to reach here tailEdges_ should only be one as well
										std::set<uint32_t> readsLeaving(
												n->getFirstOnTailEdge()->inReadNamesIdx_.begin(),
												n->getFirstOnTailEdge()->inReadNamesIdx_.end());
										std::vector<uint32_t> readsThatEnteredLeaving;
										std::set_intersection(readsEntering.begin(), readsEntering.end(),
												readsLeaving.begin(), readsLeaving.end(),
												std::back_inserter(readsThatEnteredLeaving));
										if (readsThatEnteredLeaving.size()
												!= std::min(readsEntering.size(), readsLeaving.size())
												&& readsThatEnteredLeaving.size() < currentGraph->occurenceCutOff_) {
											std::cout << "For next uid: " << n->uid_ << std::endl;
											std::cout << "\tReadsEnteringNumber    : " << readsEntering.size()
													<< std::endl;
											std::cout << "\tReadsLeavingNumber     : " << readsLeaving.size()
													<< std::endl;
											std::cout << "\tReadsLeavingThatEntered: "
													<< readsThatEnteredLeaving.size() << std::endl;
											std::cout << std::endl;
											++failCount;
										}
									}
								}
								std::cout << failCount << "/" << nodesToProcess.size() << " "
										<< failCount / nodesToProcess.size() << std::endl;
								std::cout << std::endl;
							}
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitEnds", "-", splitEndsCount));
							}
							currentGraph->collapseSingleLinkedPaths();
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitEnds", "-", splitEndsCount, "-collapsed"));
							}
							++splitEndsCount;
						}
					}
					if(extractionPars.debug){
						bfs::path groupedSeqsDir = njh::files::makeDir(currentShortTipNumberDir,njh::files::MkdirPar{ "groupedNodesBySameEdges"});
						currentGraph->writeOutNodesInGroups(groupedSeqsDir, true);
					}
					if(extractionPars.splitTailed_){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Splitting Tailed"));
						currentGraph->splitMultitailedNodes();
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitTailed-collapsed"));
						}
						currentGraph->collapseSingleLinkedPaths();
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitTailed-collapsed"));
						}
						uint32_t splitEndsCountAfterSplitTailed = 0;
						while(currentGraph->splitEndNodes(extractionPars.maxSplitEndNodeSize_)){
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitTailed-collapsed-splitEndNodes","-", splitEndsCountAfterSplitTailed));
							}
							currentGraph->collapseSingleLinkedPaths();
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitTailed-collapsed-splitEndNodes","-", splitEndsCountAfterSplitTailed, "-collapsed"));
							}
							++splitEndsCountAfterSplitTailed;
						}
					}

					std::unordered_set<std::string> kmersInOutSeqs;
					for(const auto & node : currentGraph->nodes_){
						if(node->on_){
							if(node->k_.size() > extractionPars.estimatorKlen){
								for(uint32_t pos = 0; pos < len(node->k_) + 1 - extractionPars.estimatorKlen; ++pos){
									kmersInOutSeqs.emplace(node->k_.substr(pos, extractionPars.estimatorKlen));
								}
							}
						}
					}

					if(extractionPars.trimTipsOfLowEntropyNodes_){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Collapsing One Base Indel nodes"));
						//if(currentGraph->collapseOneBaseIndelsNodes()){
						uint32_t trimTipsOfLowEntropyNodesCount = 0;
						while(currentGraph->trimTipsOfLowEntropyNodes(extractionPars.trimTipsOfLowEntropyNodesCutOff_)){
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-removeLowEntropyTips-", trimTipsOfLowEntropyNodesCount ));
							}
							//currentGraph->collapseSingleLinkedPaths();
//							if(extractionPars.debug){
//								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-collapseOneBaseIndelsNodes-",collapseOneBaseIndelsNodesCount ));
//							}
							++trimTipsOfLowEntropyNodesCount;
						}
						if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){
							currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff);
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-removeLowEntropyTips-rmHeadTailless"));
							}
						}
					}

					//trimming end nodes
					if (extractionPars.trimEdgesNodeTipsWithLowEntropy_) {
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Trimming Edges of End Nodes Once by low entropy"));
						uint32_t trimEndNodesByEntropyCount = 0;
						if(extractionPars.debug){
							std::cout << "trimEndNodesByEntropyCount: " << trimEndNodesByEntropyCount << std::endl;
						}
						while(currentGraph->trimEdgesNodeTipsWithLowEntropy(extractionPars.trimEdgesOfEndNodesPars_)){

							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-trimEndNodesByEntropy-", trimEndNodesByEntropyCount));
							}
							++trimEndNodesByEntropyCount;
							if(extractionPars.debug){
								std::cout << "trimEndNodesByEntropyCount: " << trimEndNodesByEntropyCount << std::endl;
							}
						}
						if(0 == trimEndNodesByEntropyCount){
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-final-trimEndNodesByEntropy-", trimEndNodesByEntropyCount));
							}
						}
						if(extractionPars.removeHeadlessTaillessAlongTheWay && extractionPars.removeHeadlessTaillessAfterDisentaglement){
							currentGraph->removeHeadlessTaillessNodes(headlessTaillessLenCutOff);
							if(extractionPars.debug){
								graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-splitInternals-afterTipTemoval-trimEndNodesByEntropy-rmHeadTailless"));
							}
						}
					}

					if(extractionPars.debug){
						graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("fullCollapse-final"));
					}
					if(extractionPars.writeOutFinalConnections_){
//						std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
//						for(const auto & n : currentGraph->nodes_){
//							if(!n->tailless() && !n->headless() && (n->tailCount() > 1 || n->headCount() > 1)){
//								nodesToProcess.emplace_back(n);
//							}
//						}
//						for(const auto & n : nodesToProcess){
//							auto finalConsFnp = njh::files::make_path(currentShortTipNumberDir, "finalConnectionsInfo.tab.txt");
//							OutputStream finalConsOut(finalConsFnp);
//						}
					}
//					OptimizationReconResult optRunResAfterFurtherSpliting(OptimizationReconResult::Params(currentKLen, kmerOccurenceCutOff, shortTipNumber),
//							OptimizationReconResult::Dirs(currentKmerDirectory, currentKCutOffDir, currentShortTipNumberDir),
//							kCutIter, shortTipIter);
//					if((extractionPars.splitEndsOnce_ ||
//						 extractionPars.splitEnds_ ||
//						 extractionPars.splitTailed_) &&
//						 extractionPars.removePossibleOutliersKSim){
//						//determine optimization results;
//						for(const auto & node : currentGraph->nodes_){
//							if(node->on_){
//								if(node->tailless()){
//									++optRunResAfterFurtherSpliting.taillessCount_;
//									if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
//										++optRunResAfterFurtherSpliting.taillessCountBelowLen_;
//									}
//								}
//								if(node->headless()){
//									++optRunResAfterFurtherSpliting.headlessCount_;
//									if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
//										++optRunResAfterFurtherSpliting.headlessCountBelowLen_;
//									}
//								}
//								if(node->tailless() || node->headless()){
//									++optRunResAfterFurtherSpliting.optimalCount_;
//									if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
//										++optRunResAfterFurtherSpliting.optimalCountBelowLen_;
//									}
//								}
//								if(!node->tailless() && !node->headless()){
//									if(node->k_.size() < extractionPars.optimizeNodeSizeCutOff_){
//										++optRunResAfterFurtherSpliting.internalNodesCountBelowLen_;
//									}
//								}
//							}
//						}
//					}

					if(extractionPars.writeAllPossibleHaps_){
						auto allPossibleHaps = currentGraph->generateAllPossibleHaps();
						auto allPossibleHapsOpts = SeqIOOptions::genFastaOut(njh::files::make_path(currentShortTipNumberDir, "allPossibleHaps.fasta"));
						uint32_t posHapCount = 0;
						for(auto & possHap : allPossibleHaps){
							auto hapName = njh::pasteAsStr("possHap.", leftPadNumStr<uint32_t>(posHapCount, allPossibleHaps.size()));
							possHap.seqBase_.name_ = hapName;
							for(auto & contigReg : possHap.contigLocations_){
								contigReg.chrom_ = hapName;
							}
							++posHapCount;
						}
						SeqOutput::write(allPossibleHaps, allPossibleHapsOpts);
						OutputStream possHapContigRegsOut(njh::files::make_path(currentShortTipNumberDir, "possHapContigRegs.bed"));
						for(const auto & possHap : allPossibleHaps){
							for(const auto & contigReg : possHap.contigLocations_){
								possHapContigRegsOut << contigReg.toDelimStr() << std::endl;
							}
						}
					}
					if(extractionPars.writeEstimatedMajorHaps_){
						auto tempGraph = currentGraph->copyGraph();
						tempGraph.turnOffLowestEdges();
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("estimatedMajorHaps-turnOffLowestEdges"));
						}
						tempGraph.collapseSingleLinkedPaths();
						if(extractionPars.debug){
							graphWriter.writeOutDotsAndSeqs(njh::pasteAsStr("estimatedMajorHaps-turnOffLowestEdges-collapsed"));
						}
						tempGraph.sortNodesBySize();
						uint32_t nodeCount = 0;
						SeqOutput estMajorHapWriter(SeqIOOptions::genFastaOut(njh::files::make_path(currentShortTipNumberDir, "estimatedMajorHaps-final-collapsed.fasta")));
						estMajorHapWriter.openOut();
						for(const auto & node : tempGraph.nodes_){
							seqInfo seq("contig" + leftPadNumStr<uint32_t>(nodeCount, currentGraph->nodes_.size()), node->k_);
							seq.cnt_ = node->inReadNamesIdx_.size();
							seq.updateName();
							MetaDataInName meta;
							meta.addMeta("sample", sampName);
							meta.addMeta("length", len(seq));
							meta.addMeta("headless", node->headless());
							meta.addMeta("tailless", node->tailless());
							if(extractionPars.addGroups_){
								meta.addMeta("group", node->group_);
							}
							meta.resetMetaInName(seq.name_, seq.name_.find("_t"));
							estMajorHapWriter.write(seq);
							++nodeCount;
						}
					}
					watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Creating out seqs"));
					auto outSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(currentShortTipNumberDir, "output.fasta"));
					std::vector<seqInfo> outSeqs;
					currentGraph->sortNodesBySize();
					OutOptions outInfoOpts(njh::files::make_path(currentShortTipNumberDir, "outputInfo.tab.txt"));
					OutputStream outInfoFile(outInfoOpts);
					outInfoFile << "name\treadCount\tfraction" << "\n";
					OutOptions outInfoNamesOpts(njh::files::make_path(currentShortTipNumberDir, "inputNames.tab.txt"));
					std::unique_ptr<OutputStream> outInfoNamesFilePtr;
					if(extractionPars.debug){
						outInfoNamesFilePtr = std::make_unique<OutputStream>(outInfoNamesOpts);
						(*outInfoNamesFilePtr) << "name\treadNames" << "\n";
					}
					uint32_t nodeCount = 0;
					double totalCount = 0;
					for(const auto & node : currentGraph->nodes_){
						totalCount += node->inReadNamesIdx_.size();
					}
					MetaDataInName inputMeta;
					if(nullptr != meta && meta->hasSample(sampName)){
						inputMeta = meta->getMetaForSample(sampName, getVectorOfMapKeys(meta->groupData_));
					}
					std::unordered_map<std::string, std::unordered_set<std::string>> outSeqNamesToReadNames;
					uint32_t trimEndsByLength = extractionPars.trimEndsBy;
					if(0 == extractionPars.trimEndsBy ){
						trimEndsByLength = currentKLen + 5;
					}
					if(extractionPars.addGroups_){
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
						"- Creating out seqs - creating groups"));
						currentGraph->resetGroups();
						if(extractionPars.writeGroupInputNames_){
							watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
							"- Creating out seqs - writing out group names"));
							auto inputNamesDirFnp = njh::files::make_path(currentShortTipNumberDir, "inputNamesPerGroup");
							njh::files::makeDir(njh::files::MkdirPar{inputNamesDirFnp});
							std::unordered_map<uint32_t, std::unordered_set<std::string>> namesPerGroup;
							for(const auto & node : currentGraph->nodes_){
								for(const auto & nameIdx : node->inReadNamesIdx_){
									namesPerGroup[node->group_].emplace(njh::replaceString(currentGraph->readNames_[nameIdx], "_mate", ""));
								}
							}
							auto maxGroup = vectorMaximum(njh::getVecOfMapKeys(namesPerGroup));
							for(const auto & groupNames : namesPerGroup){
								OutOptions nameOutOpts(njh::files::make_path(inputNamesDirFnp, njh::leftPadNumStr(groupNames.first, maxGroup) + ".txt"));
								OutputStream nameOut(nameOutOpts);
								std::set<std::string> outNames(groupNames.second.begin(),groupNames.second.end() );
								nameOut << njh::conToStr(outNames, "\n") << std::endl;
							}
						}
						watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
						"- Creating out seqs - rest after creating groups"));
					}
					for(const auto & node : currentGraph->nodes_){
						//if not keeping cycles, and the group is not set (which only happens with cycles) and we did add gropups
						if(!extractionPars.keepCycles_ && std::numeric_limits<uint32_t>::max() == node->group_ && extractionPars.addGroups_){
							continue;
						}
						seqInfo seq("contig" + leftPadNumStr<uint32_t>(nodeCount, currentGraph->nodes_.size()), node->k_);
						if(extractionPars.addSeqOfSingleHeadAndTailSeqs_){
							if(1 == node->headCount()){
								for(const auto & edge : node->headEdges_){
									if(edge->on_){
										std::string head = edge->head_.lock()->k_;
										seq.prepend(head.substr(0, head.size() - currentGraph->klen_ + 1));
										if(extractionPars.trimEnds && edge->head_.lock()->headless() && len(seq) > (trimEndsByLength)){
											readVecTrimmer::trimOffForwardBases(seq, trimEndsByLength);
										}
										break;
									}
								}
							}
							if(1 == node->tailCount()){
								for(const auto & edge : node->tailEdges_){
									if(edge->on_){
										std::string tail = edge->tail_.lock()->k_;
										seq.append(tail.substr(currentGraph->klen_ - 1));
										if(extractionPars.trimEnds && node->tailless() && len(seq) > (trimEndsByLength)){
											readVecTrimmer::trimOffEndBases(seq, trimEndsByLength);
										}
										break;
									}
								}
							}
						}
						if(extractionPars.trimEnds && node->headless() && len(seq) > (trimEndsByLength)){
							readVecTrimmer::trimOffForwardBases(seq, trimEndsByLength);
						}
						if(extractionPars.trimEnds && node->tailless() && len(seq) > (trimEndsByLength)){
							readVecTrimmer::trimOffEndBases(seq, trimEndsByLength);
						}
						seq.cnt_ = node->inReadNamesIdx_.size();
						seq.updateName();
						MetaDataInName meta;
						meta.addMeta("sample", sampName);
						meta.addMeta("length", len(seq));
						meta.addMeta("headless", node->headless());
						meta.addMeta("tailless", node->tailless());
						if(extractionPars.addGroups_){
							meta.addMeta("group", node->group_);
						}
						meta.addMeta(inputMeta, true);
						meta.resetMetaInName(seq.name_, seq.name_.find("_t"));
						outSeqs.emplace_back(seq);
						if(extractionPars.debug || extractionPars.filterOffOutlierInputSeqs){
							std::unordered_set<std::string> readNamesForSeq;
							for(const auto & nameIdx : node->inReadNamesIdx_){
								readNamesForSeq.emplace(njh::replaceString(currentGraph->readNames_[nameIdx], "_mate", ""));
							}
							outSeqNamesToReadNames[seq.name_] = readNamesForSeq;
						}
						++nodeCount;
					}
					if(extractionPars.addSeqOfSingleHeadAndTailSeqs_){
						njh::sort(outSeqs, [](const seqInfo & seq1, const seqInfo & seq2){
							return len(seq1) > len(seq2);
						});
						uint32_t newNodeCount = 0;
						std::regex contigPat{std::string(R"(^contig[0-9]+)")};
						std::unordered_map<std::string, std::unordered_set<std::string>> newOutSeqNamesToReadNames;
						for(auto & seq : outSeqs){
							std::string oldName = seq.name_;
							std::string newContigNamePrefix = "contig" + leftPadNumStr<uint32_t>(newNodeCount, currentGraph->nodes_.size());
							seq.name_ = std::regex_replace(seq.name_, contigPat, newContigNamePrefix	);
							++newNodeCount;
							if(extractionPars.debug || extractionPars.filterOffOutlierInputSeqs){
								newOutSeqNamesToReadNames[seq.name_] = outSeqNamesToReadNames[oldName];
							}
						}
						if(extractionPars.debug || extractionPars.filterOffOutlierInputSeqs){
							outSeqNamesToReadNames = newOutSeqNamesToReadNames;
						}
					}
					for(const auto & seq : outSeqs){
						outInfoFile << seq.name_ << "\t" << seq.cnt_ << "\t" << seq.cnt_/totalCount << "\n";
						if(extractionPars.debug){
							(*outInfoNamesFilePtr) << seq.name_ << "\t" << njh::conToStr(outSeqNamesToReadNames[seq.name_], ",") << "\n";
						}
					}
					SeqOutput::write(outSeqs, outSeqOpts);
					watch.startNewLap(njh::pasteAsStr(genPrefixForWatchLapName(),
					"- Processing out seqs"));
					uint64_t maxLen  = 0;
					readVec::getMaxLength(outSeqs, maxLen);
					//filtering
					auto writeOutDebugSeqs = [&currentShortTipNumberDir](const std::vector<seqInfo> & keptDebugSeqs, const std::string & suffix){
						OutOptions keptOutInfoOpts(njh::files::make_path(currentShortTipNumberDir, "output" + suffix + ".tab.txt"));
						OutputStream keptOutInfoFile(keptOutInfoOpts);
						keptOutInfoFile << "name\treadCount\tfraction" << "\n";
						double total = 0;
						for(const auto & seq : keptDebugSeqs){
							total += getSeqBase(seq).cnt_;
						}
						for(const auto & seq : keptDebugSeqs){
							keptOutInfoFile << getSeqBase(seq).name_ << "\t" << getSeqBase(seq).cnt_ << "\t" << getSeqBase(seq).cnt_/total << std::endl;;
						}
						auto outSeqOptsKept = SeqIOOptions::genFastaOut(njh::files::make_path(currentShortTipNumberDir, "output" + suffix + ".fasta"));
						SeqOutput::write(keptDebugSeqs, outSeqOptsKept);
					};
					//filter on length
					double finalReadLengthCutOff = extractionPars.lenCutOff;
					double finalHeadlessTaillessCutOff = extractionPars.taillessOrHeaddLesslenCutOff;

					if(0 == finalReadLengthCutOff){
						/**@todo need a better length cut off*/
						//finalReadLengthCutOff = bestResult.readLengthMedian_;
						finalReadLengthCutOff = currentKLen - 1;
					}
					if(0 == finalHeadlessTaillessCutOff){
						finalHeadlessTaillessCutOff = finalReadLengthCutOff;
					}
					if(finalReadLengthCutOff > 0 && !outSeqs.empty()){
						std::vector<seqInfo> seqsAboveLenCutOff;
						auto outSeqBelowLenOpts = SeqIOOptions::genFastaOut(njh::files::make_path(currentShortTipNumberDir, "output_belowCutOff.fasta"));
						SeqOutput writer(outSeqBelowLenOpts);
						for(const auto & seq : outSeqs){
							bool tailless = false;
							bool headless = false;
							if(MetaDataInName::nameHasMetaData(seq.name_)){
								MetaDataInName seqMeta(seq.name_);
								if(seqMeta.containsMeta("headless") && seqMeta.getMeta<bool>("headless")){
									headless = true;
								}
								if(seqMeta.containsMeta("tailless") && seqMeta.getMeta<bool>("tailless")){
									tailless = true;
								}
							}
							if(((headless || tailless) && len(seq) > finalHeadlessTaillessCutOff)  || len(seq)  > finalReadLengthCutOff){
								seqsAboveLenCutOff.emplace_back(seq);
							} else {
								writer.openWrite(seq);
							}
						}
						outSeqs = seqsAboveLenCutOff;
						if(extractionPars.debug){
							writeOutDebugSeqs(outSeqs, "_aboveLenCutOff");
						}
						if(outSeqs.empty()){
							shortTipLogErrorLog << "No outputs above cut off of " << finalReadLengthCutOff << std::endl;
						}
					}
					if(extractionPars.removePossibleOutliersKSim && !outSeqs.empty()){
						auto outSeqOptsOutliers = SeqIOOptions::genFastaOut(njh::files::make_path(currentShortTipNumberDir, "output_possibleOutliersWithKSim.fasta"));
						SeqOutput outliersWriter(outSeqOptsOutliers);
						std::vector<kmerInfo> kinfos;
						std::vector<kmerInfo> kinfosFront;
						std::vector<kmerInfo> kinfosBacks;
						for(const auto & refSeq : extractionPars.trimSeqs.empty() ? extractionPars.inputSeqs : extractionPars.trimSeqs){
							kinfos.emplace_back(refSeq.seq_, extractionPars.outlierKLen, false);
							if(len(refSeq) > 100){
								kinfosFront.emplace_back(refSeq.seq_.substr(0,75), extractionPars.outlierKLen, false);
								kinfosBacks.emplace_back(refSeq.seq_.substr(refSeq.seq_.size()-75), extractionPars.outlierKLen, false);
							}
						}
						std::vector<seqInfo> keepSeqs;
						std::vector<seqInfo> possibleOutlier;
						for(const auto & outSeq : outSeqs){
							kmerInfo outSeqKinfo(outSeq.seq_, extractionPars.outlierKLen, false);
							bool keep = false;
							for(const auto & refInfo : kinfos){
								if(refInfo.compareKmers(outSeqKinfo).second >= extractionPars.outliersKSimCutOff){
									keep = true;
									break;
								}
							}
							//below is an addition for short reconstructions to make sure not to remove ends pieces
							//that extended beyond ref seq area due to expanded region recruitment
							if(!keep && len(outSeq) > 100 && len(outSeq) < 400){
								kmerInfo outSeqKinfoBack(outSeq.seq_.substr(outSeq.seq_.size() -75), extractionPars.outlierKLen, false);
								for(const auto & refInfo : kinfosFront){
									if(refInfo.compareKmers(outSeqKinfoBack).second >= extractionPars.outliersKSimCutOff){
										keep = true;
										break;
									}
								}
							}
							//below is an addition for short reconstructions to make sure not to remove ends pieces
							//that extended beyond ref seq area due to expanded region recruitment
							if(!keep && len(outSeq) > 100 && len(outSeq) < 400){
								kmerInfo outSeqKinfoFront(outSeq.seq_.substr(0,75), extractionPars.outlierKLen, false);
								for(const auto & refInfo : kinfosBacks){
									if(refInfo.compareKmers(outSeqKinfoFront).second >= extractionPars.outliersKSimCutOff){
										keep = true;
										break;
									}
								}
							}
							if(!keep){
								possibleOutlier.emplace_back(outSeq);
							}else{
								keepSeqs.emplace_back(outSeq);
							}
						}
						if(!possibleOutlier.empty()){
							if(extractionPars.keepGroupsForOutlierFilter){
								std::vector<seqInfo> keepGroupSeqs;
								std::set<std::string> keepGroups;
								for(const auto & seq : keepSeqs){
									MetaDataInName seqMeta(seq.name_);
									keepGroups.emplace(seqMeta.getMeta("group"));
								}
								possibleOutlier.clear();
								for(const auto & seq : outSeqs){
									MetaDataInName seqMeta(seq.name_);
									if(njh::in(seqMeta.getMeta("group"), keepGroups)){
										keepGroupSeqs.emplace_back(seq);
									}else{
										possibleOutlier.emplace_back(seq);
									}
								}
								outSeqs = keepGroupSeqs;
							}else{
								outSeqs = keepSeqs;
							}
							if(!possibleOutlier.empty()){
								//update optimization
								if ((  extractionPars.splitEndsOnce_
										|| extractionPars.splitEnds_
										|| extractionPars.splitTailed_)) {
									//below is currently not used as it is not optimal to use the optimization after splitting ends
									//optRunRes = optRunResAfterFurtherSpliting;
								}
								// have to check to see if greater than zero as number can increase after spliting ends
								for(const auto & posOut : possibleOutlier){
									MetaDataInName meta(posOut.name_);
									if(meta.getMeta<bool>("tailless")){
										if(optRunRes.taillessCount_ > 0){
											--optRunRes.taillessCount_;
										}
										if(meta.getMeta<uint32_t>("length") < extractionPars.optimizeNodeSizeCutOff_){
											if(optRunRes.taillessCountBelowLen_ > 0){
												--optRunRes.taillessCountBelowLen_;
											}
										}
									}
									if(meta.getMeta<bool>("headless")){
										if(optRunRes.headlessCount_ > 0){
											--optRunRes.headlessCount_;
										}
										if(meta.getMeta<uint32_t>("length") < extractionPars.optimizeNodeSizeCutOff_){
											if(optRunRes.headlessCountBelowLen_ > 0){
												--optRunRes.headlessCountBelowLen_;
											}
										}
									}
									if(meta.getMeta<bool>("tailless") || meta.getMeta<bool>("headless")){
										if(optRunRes.optimalCount_ > 0){
											--optRunRes.optimalCount_;
										}
										if(meta.getMeta<uint32_t>("length") < extractionPars.optimizeNodeSizeCutOff_){
											if(optRunRes.optimalCountBelowLen_ > 0){
												--optRunRes.optimalCountBelowLen_;
											}
										}
									}
								}
								outliersWriter.openWrite(possibleOutlier);
								if(extractionPars.filterOffOutlierInputSeqs){
									//not updating final read name counts because it would seem like the run did even worse
									std::unordered_set<std::string> keepSeqNames;
									for(const auto & seq : outSeqs){
										for(const auto & name : outSeqNamesToReadNames[seq.name_]){
											keepSeqNames.emplace(njh::replaceString(name, "_mate", ""));
										}
										//keepSeqNames.insert(outSeqNamesToReadNames[seq.name_].begin(), outSeqNamesToReadNames[seq.name_].end());
									}
									OutputStream keepSeqNamesOut(njh::files::make_path(currentShortTipNumberDir, "keepSeqNames_fromOutlierFilter.txt"));
									keepSeqNamesOut << njh::conToStr(keepSeqNames, "\n") << std::endl;
								}
							}
						}
						if(extractionPars.debug){
							writeOutDebugSeqs(outSeqs, "_passKmerComp");
						}
						if(outSeqs.empty()){
							shortTipLogErrorLog << "No outputs passing kmer comparison filter" << std::endl;
						}
					}
					if(extractionPars.removePossibleOutliersWithMuscle && !outSeqs.empty()){
						Muscler mRunner;
						auto split = mRunner.splitSeqsOnOverlapInMultipleAln(outSeqs, extractionPars.trimSeqs.empty() ? extractionPars.inputSeqs : extractionPars.trimSeqs, extractionPars.outliersCutOff, true);
						outSeqs = split.first;
						if(split.second.size() > 0){
							auto outSeqOptsOutliers = SeqIOOptions::genFastaOut(njh::files::make_path(currentShortTipNumberDir, "output_possibleOutliersWithMuscle.fasta"));
							SeqOutput::write(split.second, outSeqOptsOutliers);
						}
						if(extractionPars.debug){
							writeOutDebugSeqs(outSeqs, "_passMuscleFilter");
						}
						if(outSeqs.empty()){
							shortTipLogErrorLog << "No outputs passing muscle alignment comparison filter" << std::endl;
						}
					}
					shortTipLog["finalNodeReadCounts"] = njh::json::toJson(finalNodeReadCounts.size());

					if(extractionPars.calcPercentUsedByKmerUsage && totalInputSequences > 200){

						uint32_t totalReads = 0;
						uint32_t totalReadsAboveCutOff = 0;

						//calculating percentage used
						std::vector<double> percentages;
						std::vector<double> percentagesByBases;

						auto usedSinglesOpts =  SeqIOOptions::genFastqIn(njh::files::make_path(finalCurrentDir, "filteredSingles.fastq"));
						auto usedPairedOpts = SeqIOOptions::genPairedIn(njh::files::make_path(finalCurrentDir, "filteredExtractedPairs_R1.fastq"), njh::files::make_path(finalCurrentDir, "filteredExtractedPairs_R2.fastq"));
						if(usedSinglesOpts.inExists()){
							SeqInput reader(usedSinglesOpts);
							seqInfo seq;
							reader.openIn();
							while(reader.readNextRead(seq)){
								uint32_t kmersFound = 0;
								uint32_t totalKmers = 0;
								for(uint32_t pos = 0; pos < len(seq) + 1 - extractionPars.estimatorKlen; ++pos){
									auto k = seq.seq_.substr(pos, extractionPars.estimatorKlen);
									if(debugEstimatingGraph.kCounts_[k] > 1){
										if(njh::in(k, kmersInOutSeqs)){
											++kmersFound;
										}
										++totalKmers;
									}
								}
								if(totalKmers > 0){
									auto percentOfKmersUsed = static_cast<double>(kmersFound)/totalKmers;
									++totalReads;
									if(percentOfKmersUsed > extractionPars.percentageOfKmersUsedCutOff){
										++totalReadsAboveCutOff;
									}
									percentages.emplace_back(percentOfKmersUsed);
								}
							}
						}
						if(usedPairedOpts.inExists()){
							SeqInput reader(usedPairedOpts);
							PairedRead pSeq;
							reader.openIn();
							while(reader.readNextRead(pSeq)){
								uint32_t kmersFound = 0;
								uint32_t totalKmers = 0;
								if(len(pSeq.seqBase_) > extractionPars.estimatorKlen){

									for(uint32_t pos = 0; pos < len(pSeq.seqBase_) + 1 - extractionPars.estimatorKlen; ++pos){
										auto k = pSeq.seqBase_.seq_.substr(pos, extractionPars.estimatorKlen);
										if(debugEstimatingGraph.kCounts_[k] > 1){
											if(njh::in(k, kmersInOutSeqs)){
												++kmersFound;
											}
											++totalKmers;
										}
									}
								}
								if(len(pSeq.mateSeqBase_) > extractionPars.estimatorKlen){
									for(uint32_t pos = 0; pos < len(pSeq.mateSeqBase_) + 1 - extractionPars.estimatorKlen; ++pos){
										auto k = pSeq.mateSeqBase_.seq_.substr(pos, extractionPars.estimatorKlen);
										if(debugEstimatingGraph.kCounts_[k] > 1){
											if(njh::in(k, kmersInOutSeqs)){
												++kmersFound;
											}
											++totalKmers;
										}
									}
								}
								if(totalKmers > 0){
									auto percentOfKmersUsed = static_cast<double>(kmersFound)/totalKmers;
									++totalReads;
									if(percentOfKmersUsed > extractionPars.percentageOfKmersUsedCutOff){
										++totalReadsAboveCutOff;
									}
									percentages.emplace_back(percentOfKmersUsed);
								}
							}
						}
						if(extractionPars.debug){
							OutputStream percentagesUsedOut(njh::files::make_path(currentShortTipNumberDir, "percentagesUsedByKmer.txt"));
							percentagesUsedOut << njh::conToStr(percentages, "\n") << std::endl;
							OutputStream percentagesUsedInfoOut(njh::files::make_path(currentShortTipNumberDir, "percentagesUsedInfo.txt"));
							percentagesUsedInfoOut << "totalInputSequences\ttotalReads\ttotalReadsAboveCutOff\tpercentageAbove\t"<< std::endl;
							percentagesUsedInfoOut<< totalInputSequences
									<< "\t" << totalReads
									<< "\t" << totalReadsAboveCutOff
									<< "\t" << static_cast<double>(totalReadsAboveCutOff)/totalReads
									<< std::endl;
						}
						optRunRes.percentOfInputUsed_ = static_cast<double>(totalReadsAboveCutOff)/totalReads;
						shortTipLog["percentOfInputReadsUsed"] = njh::json::toJson(optRunRes.percentOfInputUsed_);
						shortTipLog["totalReadsAboveKmersUsedCutOff"] = njh::json::toJson(totalReadsAboveCutOff);
					}else{
						optRunRes.percentOfInputUsed_ = static_cast<double>(finalNodeReadCounts.size())/totalInputSequences;
						shortTipLog["percentOfInputReadsUsed"] = njh::json::toJson(optRunRes.percentOfInputUsed_);
					}


					//kcutRes.foundTrimmed_ = foundTrimmed;
					//kcutRes.numberOfFinalFilteredSeqs_ = finalFilteredSeqs.size();
					optRunRes.numberOfFinalFilteredSeqs_ = outSeqs.size();
					optRunRes.readLengthAverage_= vectorMean(readLengths);
					optRunRes.readLengthMedian_ = medReadLeng;
					allOptRunResults.emplace_back(optRunRes);
					allCurrentKCutOffOptResults.emplace_back(optRunRes);
					allCurrentShortTipOptRunResults.emplace_back(optRunRes);

					shortTipLog["headlessCount"] = njh::json::toJson(optRunRes.headlessCount_);
					shortTipLog["taillessCount"] = njh::json::toJson(optRunRes.taillessCount_);
					shortTipLog["headless_plus_taillessCount"] = njh::json::toJson(optRunRes.headlessCount_ + optRunRes.taillessCount_);
					shortTipLog["optimalCount"] = njh::json::toJson(optRunRes.optimalCount_);
					shortTipLog["headlessCountBelowLen"] = njh::json::toJson(optRunRes.headlessCountBelowLen_);
					shortTipLog["taillessCountBelowLen"] = njh::json::toJson(optRunRes.taillessCountBelowLen_);
					shortTipLog["headlessCountBelowLen_plus_taillessCountBelowLen"] = njh::json::toJson(optRunRes.headlessCountBelowLen_ + optRunRes.taillessCountBelowLen_);
					shortTipLog["optimalCountBelowLen"] = njh::json::toJson(optRunRes.optimalCountBelowLen_);
					//kmerOccurenceCutOff += extractionPars.optimizeKcutStep;
					shortTipLog["errorLog"] = njh::json::toJson(shortTipLogErrorLog.str());
					if(extractionPars.debug){
						auto usedSinglesFnp = singletOuts.getPriamryOutName();
						auto usedPairedR1Fnp = pairedOpts.getPriamryOutName();
						auto usedPairedR2Fnp = pairedOpts.getSecondaryOutName();
						std::vector<double> readPercentagesUsed;
						if(bfs::exists(usedPairedR1Fnp)){
							SeqInput reader(SeqIOOptions::genPairedIn(usedPairedR1Fnp, usedPairedR2Fnp));
							reader.openIn();
							PairedRead seq;
							while(reader.readNextRead(seq)){
								if(seq.seqBase_.seq_.size() > currentGraph->klen_){
									uint32_t kmersUsed = 0;
									for (auto pos : iter::range(seq.seqBase_.seq_.size() - currentGraph->klen_ + 1)) {
										if(currentGraph->kCounts_[seq.seqBase_.seq_.substr(pos, currentGraph->klen_)] > currentGraph->occurenceCutOff_){
											++kmersUsed;
										}
									}
									readPercentagesUsed.push_back(static_cast<double>(kmersUsed)/(seq.seqBase_.seq_.size() - currentGraph->klen_ + 1));
								}
							}
						}
						if(bfs::exists(usedSinglesFnp)){
							SeqInput reader(SeqIOOptions::genFastqIn(usedSinglesFnp));
							reader.openIn();
							seqInfo seq;
							while(reader.readNextRead(seq)){
								if(seq.seq_.size() > currentGraph->klen_){
									uint32_t kmersUsed = 0;
									for (auto pos : iter::range(seq.seq_.size() - currentGraph->klen_ + 1)) {
										if(currentGraph->kCounts_[seq.seq_.substr(pos, currentGraph->klen_)] > currentGraph->occurenceCutOff_){
											++kmersUsed;
										}
									}
									readPercentagesUsed.push_back(static_cast<double>(kmersUsed)/(seq.seq_.size() - currentGraph->klen_ + 1));
								}
							}
						}
						auto & readPercentagesUsedAboveCounts = shortTipLog["readPercentagesUsedAboveCounts"];
						for(const auto cutOff : iter::range(0.20, 1.0, .10)){
							double aboveCutOff = std::count_if(readPercentagesUsed.begin(), readPercentagesUsed.end(), [&cutOff](const double val){
								return val >=cutOff;
							});
							readPercentagesUsedAboveCounts["Above-" + estd::to_string(cutOff)] = njh::json::toJson(aboveCutOff/readPercentagesUsed.size());
						}
						auto readPercentagesUsedStats = getStatsOnVec(readPercentagesUsed);
						shortTipLog["readPercentagesUsedStats"] = njh::json::toJson(readPercentagesUsedStats);
					}
					if(extractionPars.writeOutFinalDot_){
						OutputStream rectDotOut(
								njh::files::make_path(currentShortTipNumberDir,
										"output_final.dot"));
						currentGraph->writeRectangleWithEstimatedCovDot(rectDotOut,
								debugEstimatingGraph);
						if (extractionPars.writeOutDotsByGroups_) {
							currentGraph->writeRectangleWithEstimatedCovDotByGroup(
									njh::files::make_path(currentShortTipNumberDir,
											"output_final"), debugEstimatingGraph,
									extractionPars.writeOutDotsByGroupsSizeCutOff_);
						}
					}
					auto outSeqFinalOpts = SeqIOOptions::genFastaOut(njh::files::make_path(currentShortTipNumberDir, "output_finalFiltered.fasta"));
					SeqOutput::write(outSeqs, outSeqFinalOpts);
	//				std::multimap<uint32_t, double> approxCoverages;
	//				uint32_t longestLength = 0;
	//				double maxApproxCoverage = 1;
	//				for(const auto & node : currentGraph->nodes_){
	//					double nCov = static_cast<double>(node->cnt_)/(node->k_.length() - node->kLen_ + 1);
	//					if(longestLength < node->k_.length()){
	//						longestLength = node->k_.length();
	//						maxApproxCoverage = nCov;
	//					}
	//					approxCoverages.emplace(node->k_.length(), nCov);
	//				}
	//
	//				kcutLog["approxCoverages"] = njh::json::toJson(approxCoverages);
	//				kcutLog["maxApproxCoverage"] = njh::json::toJson(maxApproxCoverage);
	//				kcutLog["longestLength"] = njh::json::toJson(longestLength);
//					if(2 == optRunRes.runParams_.kcut_ &&
//							25 == optRunRes.runParams_.klen_&&
//							2 == optRunRes.runParams_.shortNumber_){
//
//						std::cout << optRunRes.toJson() << std::endl;
//						exit(1);
//					}
					if(!extractionPars.forceFullShortTipOptimize_ && shortTipIter >= extractionPars.shortTipOptimizationAttempts){
						auto bestOptResults = OptimizationReconResult::getBestResults(
								allCurrentShortTipOptRunResults,
								extractionPars.percentOfInputReadsUsedForOptimization,
								extractionPars.percentOfInputReadsUsedForOptimizationStepDown,
								extractionPars.useFullOptimalCount_);
						if(!bestOptResults.empty()){
							//sort while optimizing first for lowest occurrence cut off, and then shortest tip cut off and then lowest kmer size
							njh::sort(bestOptResults,OptimizationReconResult::sortFunc);
							OptimizationReconResult bestRes = bestOptResults.front();
							if(optRunRes.shortTipIterNumber_ - bestRes.shortTipIterNumber_ >= extractionPars.shortTipOptimizationAttempts){
								break; //break out of the kmerOccurenceCutOff optimization loop, haven't seen any improvements in kCutOptimizationAttempts (default 5) iters
							}
						}
					}
					++shortTipIter;
					if(!extractionPars.trimShortTips_){
						break;
					}
				} catch (std::exception & exp) {
					//kmerOccurenceCutOff += extractionPars.optimizeKcutStep;
					if(extractionPars.debug){
						std::cerr << exp.what() << std::endl;
					}
					if(extractionPars.exitOnException_){
						exit(1);
					}
					shortTipLog["exception"] = njh::json::toJson(exp.what());
					if(!extractionPars.trimShortTips_){
						break;
					}
				}
			}
			if(!extractionPars.forceFullKcutOptimize_ && kCutIter >= extractionPars.kCutOptimizationAttempts){
				auto bestOptResults = OptimizationReconResult::getBestResults(
						allCurrentKCutOffOptResults,
						extractionPars.percentOfInputReadsUsedForOptimization,
						extractionPars.percentOfInputReadsUsedForOptimizationStepDown,
						extractionPars.useFullOptimalCount_);
//				std::cout << "extractionPars.useFullOptimalCount_: " << njh::colorBool(extractionPars.useFullOptimalCount_) << std::endl;
//				std::cout << "bestOptResults.size(): " << bestOptResults.size() << std::endl;
//				for(const auto & res : bestOptResults){
//					std::cout << res.runParams_.toJson() << std::endl;
//					std::cout << "res.optimalCountBelowLen_: " << res.optimalCountBelowLen_ << std::endl;
//					std::cout << "res.optimalCount_: " << res.optimalCount_ << std::endl;
//				}
				if(!bestOptResults.empty()){
					//sort while optimizing first for lowest occurrence cut off, and then shortest tip cut off and then lowest kmer size
					njh::sort(bestOptResults, OptimizationReconResult::sortFunc);
					OptimizationReconResult bestRes = bestOptResults.front();
//					std::cout << "kmerOccurenceCutOff:" << kmerOccurenceCutOff << std::endl;
//					std::cout << "kCutIter: " << kCutIter << std::endl;
//					std::cout << "bestOptResults: " << bestOptResults.size() << std::endl;
//					std::cout << bestRes.toJson() << std::endl;
					if(kCutIter - bestRes.kcutIterNumber_ >= extractionPars.kCutOptimizationAttempts){
						break; //break out of the kmerOccurenceCutOff optimization loop, haven't seen any improvements in kCutOptimizationAttempts (default 5) iters
					}
				}
			}
			++kCutIter;
		}

		//currentKLen += extractionPars.optimizeKLenStep;
	}
//	{
//		double percentUsed = extractionPars.percentOfInputReadsUsedForOptimization;
//		while(nullptr == bestResult && percentUsed > 0){
//
//			for(auto & klenRes : resultsForKLens){
//				for(const auto & kcutRes : klenRes.kcutResults){
//					if(kcutRes.percentOfInputUsed_ > percentUsed &&
//							0 != kcutRes.numberOfFinalFilteredSeqs_ &&
//							nullptr == bestResult){
//						klenRes.bestKcutResult = std::make_unique<StitchResForKCut>(kcutRes);
//						bestResult = std::make_unique<StitchResForKLen>(klenRes);
//					}else{
//						if(kcutRes.percentOfInputUsed_ > percentUsed
//								&& 0 != kcutRes.numberOfFinalFilteredSeqs_
//								&& (extractionPars.useFullOptimalCount_ ? kcutRes.optimalCount_ <= bestResult->bestKcutResult->optimalCount_ : kcutRes.optimalCountBelowLen_ <= bestResult->bestKcutResult->optimalCountBelowLen_)){
////								(kcutRes.foundTrimmed_ ||
////								kcutRes.optimalCount_ < bestResult->bestKcutResult->optimalCount_)){
//								//res.headlessCount_ + res.taillessCount_ < bestResult->headlessCount_ + bestResult->taillessCount_)){
//							//optimize for lowest kcut
//							if(extractionPars.useFullOptimalCount_ ? kcutRes.optimalCount_ == bestResult->bestKcutResult->optimalCount_ : kcutRes.optimalCountBelowLen_ == bestResult->bestKcutResult->optimalCountBelowLen_){
//								if(kcutRes.kcut_ >= bestResult->bestKcutResult->kcut_){
//									continue;
//								}
//							}
//							if(klenRes.klen_ == bestResult->klen_){
//								bestResult->bestKcutResult = std::make_unique<StitchResForKCut>(kcutRes);
//							}else{
//								klenRes.bestKcutResult = std::make_unique<StitchResForKCut>(kcutRes);
//								bestResult = std::make_unique<StitchResForKLen>(klenRes);
//							}
//						}
//					}
//				}
//			}
//			//std::cout << "percentUsed: " << percentUsed << std::endl;
//			percentUsed -= extractionPars.percentOfInputReadsUsedForOptimizationStepDown;
//		}
//	}
//	std::cout << "extractionPars.percentOfInputReadsUsedForOptimization: " << extractionPars.percentOfInputReadsUsedForOptimization << std::endl;
//	std::cout << "allOptRunResults.size(): " << allOptRunResults.size() << std::endl;
	auto bestOptResults = OptimizationReconResult::getBestResults(
			allOptRunResults,
			extractionPars.percentOfInputReadsUsedForOptimization,
			extractionPars.percentOfInputReadsUsedForOptimizationStepDown,
			extractionPars.useFullOptimalCount_);

//	std::cout << "bestOptResults.empty(): " << njh::boolToStr(bestOptResults.empty()) << std::endl;

	if(!bestOptResults.empty()){
		//sort while optimizing first for lowest occurrence cut off, and then shortest tip cut off and then lowest kmer size
		njh::sort(bestOptResults,OptimizationReconResult::sortFunc);
		OptimizationReconResult bestResult = bestOptResults.front();
		ret.log_["readLengthAverage"] = njh::json::toJson(bestResult.readLengthAverage_);
		ret.log_["readLengthMedian"] = njh::json::toJson(bestResult.readLengthMedian_);
		ret.log_["bestShortTipNumber"] = bestResult.runParams_.shortNumber_;
		ret.log_["bestKmerOccurenceCutOff"] = bestResult.runParams_.kcut_;
		ret.log_["bestKmerLength"] = bestResult.runParams_.klen_;
		ret.log_["bestKmerParsName"] = njh::json::toJson(njh::pasteAsStr("klen_", bestResult.runParams_.klen_, "-kcut_", bestResult.runParams_.kcut_, "-shortTip_", bestResult.runParams_.shortNumber_));
		ret.log_["bestResult"] = njh::json::toJson(bestResult);
		OutputStream bestOptInfoOut(njh::files::make_path(finalCurrentDir, "optimizationInfoBest.json"));
		bestOptInfoOut << njh::json::toJson(bestResult)<< std::endl;
		OutputStream bestOptAllInfoOut(njh::files::make_path(finalCurrentDir, "optimizationInfoAll.json"));
		bestOptAllInfoOut << njh::json::toJson(allOptRunResults)<< std::endl;
		{
			watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": Getting Estimate Counts") );
			//KmerPathwayGraph estimatingGraph(bestResult.runParams_.klen_);
			KmerPathwayGraph estimatingGraph(std::min(extractionPars.estimatorKlen, *std::min_element(extractionPars.kmerLengths.begin(),extractionPars.kmerLengths.end() )));

			{
				estimatingGraph.setOccurenceCutOff(bestResult.runParams_.kcut_);
				estimatingGraph.debug_ = extractionPars.graphDebug_;
				estimatingGraph.verbose_ = extractionPars.graphVerbose_;
//				auto usedSinglesFnp = njh::files::make_path(bestResult.runDirs_.klenDir_,  "extractedSingles.fastq");
//				auto usedPairedR1Fnp = njh::files::make_path(bestResult.runDirs_.klenDir_, "extractedPairs_R1.fastq");
//				auto usedPairedR2Fnp = njh::files::make_path(bestResult.runDirs_.klenDir_, "extractedPairs_R2.fastq");
				auto usedSinglesFnp =  njh::files::make_path(finalCurrentDir, "filteredSingles.fastq");
				auto usedPairedR1Fnp = njh::files::make_path(finalCurrentDir, "filteredExtractedPairs_R1.fastq");
				auto usedPairedR2Fnp = njh::files::make_path(finalCurrentDir, "filteredExtractedPairs_R2.fastq");
				if(bfs::exists(usedPairedR1Fnp)){
					SeqInput reader(SeqIOOptions::genPairedIn(usedPairedR1Fnp, usedPairedR2Fnp));
					reader.openIn();
					PairedRead seq;
					while(reader.readNextRead(seq)){
						kmerInfo firstInfo(seq.seqBase_.seq_, estimatingGraph.klen_, false);
						estimatingGraph.increaseKCounts(seq.seqBase_.seq_);
						//estimatingGraph.increaseKCounts(seq.mateSeqBase_.seq_);
						if(seq.mateSeqBase_.seq_.size() > estimatingGraph.klen_){
							for (auto pos : iter::range(seq.mateSeqBase_.seq_.size() - estimatingGraph.klen_ + 1)) {
								estimatingGraph.kCounts_[seq.mateSeqBase_.seq_.substr(pos, estimatingGraph.klen_)] += 1;
//								if(!njh::in(seq.mateSeqBase_.seq_.substr(pos, estimatingGraph.klen_), firstInfo.kmers_)){
//									estimatingGraph.kCounts_[seq.mateSeqBase_.seq_.substr(pos, estimatingGraph.klen_)] += 1;
//								}
							}
						}
					}
				}
				if(bfs::exists(usedSinglesFnp)){
					SeqInput reader(SeqIOOptions::genFastqIn(usedSinglesFnp));
					reader.openIn();
					seqInfo seq;
					while(reader.readNextRead(seq)){
						estimatingGraph.increaseKCounts(seq.seq_);
					}
				}
			}
			//process best result
			watch.startNewLap(njh::pasteAsStr(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 10000), ": Processing out seqs - collapsing and trimming"));
			//read in the final filtered seqs
			auto finalSeqOpts = SeqIOOptions::genFastaIn(njh::files::make_path(bestResult.runDirs_.shortTipDir_, "output_finalFiltered.fasta"));
			finalSeqOpts.processed_ = true;
			auto outSeqs = SeqInput::getSeqVec<seqInfo>(finalSeqOpts);
			//sort by length
			njh::sort(outSeqs, [](const seqInfo & seq1, const seqInfo & seq2){
				return len(seq1) == len(seq2) ? seq1.cnt_ > seq2.cnt_ : len(seq1) > len(seq2);
			});
			auto writeOutDebugSeqs = [&bestResult](const std::vector<seqInfo> & keptDebugSeqs, const std::string & suffix){
				OutOptions keptOutInfoOpts(njh::files::make_path(bestResult.runDirs_.shortTipDir_, "output" + suffix + ".tab.txt"));
				OutputStream keptOutInfoFile(keptOutInfoOpts);
				keptOutInfoFile << "name\treadCount\tfraction" << "\n";
				double total = 0;
				for(const auto & seq : keptDebugSeqs){
					total += getSeqBase(seq).cnt_;
				}
				for(const auto & seq : keptDebugSeqs){
					keptOutInfoFile << getSeqBase(seq).name_ << "\t" << getSeqBase(seq).cnt_ << "\t" << getSeqBase(seq).cnt_/total << std::endl;;
				}
				auto outSeqOptsKept = SeqIOOptions::genFastaOut(njh::files::make_path(bestResult.runDirs_.shortTipDir_, "output" + suffix + ".fasta"));
				SeqOutput::write(keptDebugSeqs, outSeqOptsKept);
			};
			if(extractionPars.trimToInputSeqs && !outSeqs.empty()){
				double finalReadLengthCutOff = extractionPars.lenCutOff;
				if(0 == finalReadLengthCutOff){
					/**@todo need a better length cut off*/
					//finalReadLengthCutOff = bestResult.readLengthMedian_;
					finalReadLengthCutOff = bestResult.runParams_.klen_ - 1;
				}

				Muscler mRunner;
				if(extractionPars.trimWithGlobalAln){
					uint64_t maxLen = 0;
					readVec::getMaxLength(extractionPars.trimSeqs.empty() ? extractionPars.inputSeqs : extractionPars.trimSeqs, maxLen);
					readVec::getMaxLength(outSeqs, maxLen);
					aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);
					alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(workingDir, "trimAlnCache").string(), extractionPars.verbose);
					std::vector<kmerInfo> inputSeqsKmerInfos;
					for(const auto & input : extractionPars.trimSeqs.empty() ? extractionPars.inputSeqs : extractionPars.trimSeqs){
						inputSeqsKmerInfos.emplace_back(input.seq_, 7, false);
					}
					readVecTrimmer::trimSeqToRefByGlobalAln(outSeqs,
							extractionPars.trimSeqs.empty() ? extractionPars.inputSeqs : extractionPars.trimSeqs,
									inputSeqsKmerInfos,
									alignerObj	);
					alignerObj.processAlnInfoOutputNoCheck(njh::files::make_path(workingDir, "trimAlnCache").string(), extractionPars.verbose);
				}else{
					mRunner.trimSeqsToMultiAlnRef(outSeqs, extractionPars.trimSeqs.empty() ? extractionPars.inputSeqs : extractionPars.trimSeqs, extractionPars.mtPars);
				}
				std::vector<uint32_t> refReadLens;
				for(const auto & seq : extractionPars.trimSeqs.empty() ? extractionPars.inputSeqs : extractionPars.trimSeqs){
					refReadLens.emplace_back(len(seq));
				}
				auto minRefReadLen = *std::min_element(refReadLens.begin(), refReadLens.end());
				//collapse to unique seqs after trimming
				std::vector<seqInfo> uniqueSeqs;
				//sort by length
				njh::sort(outSeqs, [](const seqInfo & seq1, const seqInfo & seq2){
					return len(seq1) == len(seq2) ? seq1.cnt_ > seq2.cnt_ : len(seq1) > len(seq2);
				});
				for(auto & seq : outSeqs){
					MetaDataInName meta(seq.name_);
					meta.addMeta("length", len(seq), true);
					meta.addMeta("trimStatus", seq.on_, true);
					meta.resetMetaInName(seq.name_);
					seq.on_ = true;
					bool unique = true;
					for( auto & uniSeq : uniqueSeqs){
						if(uniSeq.seq_ == seq.seq_){
							unique = false;
							break;
						}
					}
					if(unique){
						uniqueSeqs.emplace_back(seq);
					}
				}
				outSeqs = uniqueSeqs;
				if(finalReadLengthCutOff > 0 && !outSeqs.empty() && finalReadLengthCutOff < minRefReadLen){
					std::vector<seqInfo> seqsAboveLenCutOff;
					auto outSeqBelowLenOpts = SeqIOOptions::genFastaOut(njh::files::make_path(bestResult.runDirs_.klenDir_, "output_belowCutOffAfterTrimming.fasta"));
					SeqOutput writer(outSeqBelowLenOpts);
					for(const auto & seq : outSeqs){
						if(len(seq)  > finalReadLengthCutOff){
							seqsAboveLenCutOff.emplace_back(seq);
						}else{
							writer.openWrite(seq);
						}
					}
					outSeqs = seqsAboveLenCutOff;
					if(extractionPars.debug){
						writeOutDebugSeqs(outSeqs, "_aboveLenCutOffAfterTrimming");
					}
					if(outSeqs.empty()){
						std::stringstream finalProcessingMessageLog;
						finalProcessingMessageLog<< "No outputs above cut off of " << finalReadLengthCutOff << " after trimming "<< std::endl;;
						ret.log_["final processing message"] =finalProcessingMessageLog.str();
					}
				}
			}
			//trim with coverage
			/**@todo need a better trimmer for coverage,  */
			if (extractionPars.trimEndsByCoverage) {
				for (auto & seq : outSeqs) {
					if (len(seq) > estimatingGraph.klen_ + 1) {
						if (estimatingGraph.kCounts_[seq.seq_.substr(0,
								estimatingGraph.klen_)]
								<= estimatingGraph.occurenceCutOff_ * 2) {
							readVecTrimmer::trimOffForwardBases(seq, estimatingGraph.klen_);
						}
					}
					if (len(seq) > estimatingGraph.klen_ + 1) {
						if (estimatingGraph.kCounts_[seq.seq_.substr(
								len(seq) - estimatingGraph.klen_, estimatingGraph.klen_)]
								<= estimatingGraph.occurenceCutOff_ * 2) {
							readVecTrimmer::trimOffEndBases(seq, estimatingGraph.klen_);
						}
					}
				}
			}

			std::unordered_map<std::string, double> coverageForReadSorting;
			for(auto & seq : outSeqs){
				auto estCov = CoverageEstimator::estimateCov(seq, estimatingGraph);
				MetaDataInName meta(seq.name_);
				double covForSorting = 0;
				if(std::numeric_limits<double>::max() != estCov.minCov_.avgCov_ ){
					meta.addMeta("estimatedPerBaseCoverage", estCov.minCov_.avgCov_, true);
					covForSorting = estCov.minCov_.avgCov_;
					if(extractionPars.debug){
						{
							OutputStream seqCovCountsOut(OutOptions(njh::files::make_path(bestResult.runDirs_.klenDir_, seq.name_ + "_coverage.txt")));
							seqCovCountsOut << "coverage" << "\n";
							seqCovCountsOut << njh::conToStr(estCov.allCounts_, "\n") << std::endl;
						}
						OutputStream seqCovCountsDetailedOut(OutOptions(njh::files::make_path(bestResult.runDirs_.klenDir_, seq.name_ + "_coverageDetailed.txt")));
						seqCovCountsDetailedOut << "pos\tcoverage\tminPos\tsdCov\tavgCov\tstart\tend\tallCounts.size()" << "\n";
						for(const auto & estCovPos : iter::range(estCov.coverageInfos_.size())){
							auto pos = estCov.coverageInfos_[estCovPos].pos_;
							seqCovCountsDetailedOut << pos
									<< "\t" << estCov.allCounts_[pos]
									<< "\t" << estCov.coverageInfos_[estCovPos].minPos_
									<< "\t" << estCov.coverageInfos_[estCovPos].sdCov_
									<< "\t" << estCov.coverageInfos_[estCovPos].avgCov_
									<< "\t" << (pos + 1 > estimatingGraph.klen_ ? pos + 1 - estimatingGraph.klen_: 0)
									<< "\t" << pos + 1
									<< "\t" << estCov.allCounts_.size() << std::endl;
						}
					}
				} else{
					meta.addMeta("estimatedPerBaseCoverage", "NA", true);
				}

				meta.resetMetaInName(seq.name_);
				coverageForReadSorting[seq.name_] = covForSorting;
			}


			uint32_t finalTotalCountAbove = 0;
			std::vector<seqInfo> finalFilteredSeqs;

			if(extractionPars.collapseFinalClusters_){
				//sort by coverage and then tie break with read counts
				//this is mostly being done for collapsing on homopolymer indels and especially in SWGA if there is insertions having reads sorted by length
				//will lead to the longer sequence (which is incorrect due to it's insertions) will be kept, so need to sort by coverage and read count
				njh::sort(outSeqs, [&coverageForReadSorting](const seqInfo & seq1, const seqInfo & seq2){
					if(coverageForReadSorting[seq1.name_] == coverageForReadSorting[seq2.name_]){
						return seq1.cnt_ > seq2.cnt_;
					}else{
						return coverageForReadSorting[seq1.name_] > coverageForReadSorting[seq2.name_];
					}
				});

				uint64_t maxLen = 0;
				readVec::getMaxLength(outSeqs, maxLen);
				//throw away any seqs that differ only in homopolymers and can be completely found in the other sequences
				bool countEndGaps = false;
				aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), countEndGaps);
				alignerObj.processAlnInfoInputNoCheck(njh::files::make_path(workingDir, "trimAlnCache").string(), extractionPars.verbose);
				alignerObj.weighHomopolymers_ = true;
				comparison passableErrors;
				passableErrors.oneBaseIndel_ = 0.75;
				passableErrors.twoBaseIndel_ = 0.75;
				passableErrors.largeBaseIndel_ = 0.75;
				passableErrors.distances_.query_.coverage_ = 1;
				std::vector<kmerInfo> finalFilteredSeqsKInfos;
				std::vector<kmerInfo> outSeqsKInfos;
				for(const auto & outSeq : outSeqs){
					outSeqsKInfos.emplace_back(outSeq.seq_, extractionPars.outlierKLen, false);
				}
				for(const auto & finalSeqPos : iter::range(outSeqs.size())){
					const auto & finalSeq = outSeqs[finalSeqPos];
					bool pass = true;
					for(const auto & otherPos : iter::range(finalFilteredSeqs.size())){
						if(outSeqsKInfos[finalSeqPos].compareKmers(finalFilteredSeqsKInfos[otherPos]).second < .80){
							continue;
						}
						const auto & other = finalFilteredSeqs[otherPos];
						if(other.seq_ == finalSeq.seq_){
							pass = false;
							break;
						}else{
							alignerObj.alignCacheGlobal(other, finalSeq);
							alignerObj.profileAlignment(other, finalSeq, false, true, false);
							if (passableErrors.passErrorProfile(alignerObj.comp_) &&
									alignerObj.comp_.distances_.query_.coverage_ >= passableErrors.distances_.query_.coverage_) {
								pass = false;
								break;
							}
						}
					}
					if(pass){
						finalTotalCountAbove += finalSeq.cnt_;
						finalFilteredSeqs.emplace_back(finalSeq);
						finalFilteredSeqsKInfos.emplace_back(finalSeq.seq_, extractionPars.outlierKLen, false);
					}
				}
				alignerObj.processAlnInfoOutputNoCheck(njh::files::make_path(workingDir, "trimAlnCache").string(), extractionPars.verbose);
			} else {
				if(!extractionPars.doNotCollapseIdenticalContigs_){
					for (const auto & finalSeqPos : iter::range(outSeqs.size())) {
						const auto & finalSeq = outSeqs[finalSeqPos];
						bool pass = true;
						for (const auto & otherPos : iter::range(finalFilteredSeqs.size())) {
							const auto & other = finalFilteredSeqs[otherPos];
							if (other.seq_ == finalSeq.seq_) {
								pass = false;
								break;
							}
						}
						if (pass) {
							finalTotalCountAbove += finalSeq.cnt_;
							finalFilteredSeqs.emplace_back(finalSeq);
						}
					}
				} else {
					for(const auto & finalSeq : outSeqs){
						finalFilteredSeqs.emplace_back(finalSeq);
						finalTotalCountAbove += finalSeq.cnt_;
					}
				}
			}



			auto outSeqOptsAboveFrac = SeqIOOptions::genFastaOut(njh::files::make_path(bestResult.runDirs_.klenDir_, "output_aboveCutOff.fasta"));
			OutOptions outInfoAboveFracOpts(njh::files::make_path(bestResult.runDirs_.klenDir_, "outputInfo_aboveCutOff.tab.txt"));
			std::ofstream outInfoAboveFracFile;
			outInfoAboveFracOpts.openFile(outInfoAboveFracFile);
			outInfoAboveFracFile << "name\treadCount\tfraction" << "\n";
			SeqOutput aboveWriter(outSeqOptsAboveFrac);
			aboveWriter.openOut();
			for(auto & seq : finalFilteredSeqs){
				outInfoAboveFracFile << seq.name_ << "\t" << seq.cnt_ << "\t" << seq.cnt_/finalTotalCountAbove << std::endl;;
				aboveWriter.write(seq);
			}
		}

		auto filesToMoveTop = njh::files::filesInFolder(bestResult.runDirs_.klenDir_);
		for(const auto & fnp : filesToMoveTop){
			if(bfs::is_regular_file(fnp)){
				if(!extractionPars.keepOptimizedSubDirs){
					bfs::rename(
							njh::files::make_path(bestResult.runDirs_.klenDir_, fnp.filename()),
							njh::files::make_path(finalCurrentDir, fnp.filename()));
				}else{
					bfs::copy_file(
							njh::files::make_path(bestResult.runDirs_.klenDir_, fnp.filename()),
							njh::files::make_path(finalCurrentDir, fnp.filename()));
					//set the mod time to the what was copied from
					bfs::last_write_time(njh::files::make_path(finalCurrentDir, fnp.filename()),
					bfs::last_write_time(njh::files::make_path(bestResult.runDirs_.klenDir_, fnp.filename())));
				}
			}
		}
		auto filesToMoveKCut = njh::files::filesInFolder(bestResult.runDirs_.shortTipDir_);
		for(const auto & fnp : filesToMoveKCut){
			if(bfs::is_regular_file(fnp)){
				if(!extractionPars.keepOptimizedSubDirs){
					bfs::rename(
							njh::files::make_path(bestResult.runDirs_.shortTipDir_, fnp.filename()),
							njh::files::make_path(finalCurrentDir, fnp.filename()));
				}else{
					bfs::copy_file(
							njh::files::make_path(bestResult.runDirs_.shortTipDir_, fnp.filename()),
							njh::files::make_path(finalCurrentDir, fnp.filename()));
					//set the mod time to the what was copied from
					bfs::last_write_time(njh::files::make_path(finalCurrentDir, fnp.filename()),
					bfs::last_write_time(njh::files::make_path(bestResult.runDirs_.shortTipDir_, fnp.filename())));
				}
			}
		}
		if(!extractionPars.keepOptimizedSubDirs){
			std::set<std::string> klenDirs;
			for(const auto & res : allOptRunResults){
				klenDirs.emplace(res.runDirs_.klenDir_.string());
			}
			for(const auto & klenDir : klenDirs){
				njh::files::rmDirForce(klenDir);
			}
		}

		if(extractionPars.filterOffOutlierInputSeqs && bfs::exists(njh::files::make_path(finalCurrentDir, "keepSeqNames_fromOutlierFilter.txt"))){
			//read in names to keep
			std::unordered_set<std::string> keepNames;
			InputStream keepNamesIn{njh::files::make_path(finalCurrentDir, "keepSeqNames_fromOutlierFilter.txt")};
			std::string line;
			while(njh::files::crossPlatGetline(keepNamesIn, line)){
				keepNames.emplace(line)	;
			}
			SeqIOOptions outlierFilteredPairedOpts = SeqIOOptions::genPairedOut(
					njh::files::make_path(finalCurrentDir, "temp_filteredExtractedPairs"));
			SeqIOOptions outlierFilteredSingletOuts = SeqIOOptions::genFastqOut(
					njh::files::make_path(finalCurrentDir, "temp_filteredSingles"));
			SeqIOOptions outlierFilteredOff_pairedOpts = SeqIOOptions::genPairedOut(
					njh::files::make_path(finalCurrentDir, "outlierFilteredOff_extractedPairs"));
			SeqIOOptions outlierFilteredOff_singletOuts = SeqIOOptions::genFastqOut(
					njh::files::make_path(finalCurrentDir, "outlierFilteredOff_extractedSingles"));
			SeqInput pairedReader(findTanStitchPars.usedFilteredPairedOpts);
			SeqInput singleReader(findTanStitchPars.usedFilteredSinglesOpts);
			//std::cout << "usedFilteredPairedOpts: " << usedFilteredPairedOpts.firstName_ <<  std::endl;
			//std::cout << "usedFilteredPairedOpts: " << usedFilteredPairedOpts.secondName_ <<  std::endl;
			if(pairedReader.ioOptions_.inExists()){
				pairedReader.openIn();
				PairedRead pseq;
				SeqOutput keepWriter(outlierFilteredPairedOpts);
				SeqOutput filterWriter(outlierFilteredOff_pairedOpts);
				while(pairedReader.readNextRead(pseq)){
					if(njh::in(pseq.seqBase_.name_, keepNames)){
						keepWriter.openWrite(pseq);
					}else{
						filterWriter.openWrite(pseq);
					}
				}
				pairedReader.closeIn();
				bfs::remove(pairedReader.ioOptions_.firstName_);
				bfs::remove(pairedReader.ioOptions_.secondName_);
				if(keepWriter.outOpen()){
					auto primOutFnp = keepWriter.getPrimaryOutFnp();
					auto secdOutFnp = keepWriter.getSecondaryOutFnp();
					keepWriter.closeOut();
					bfs::rename(primOutFnp, pairedReader.ioOptions_.firstName_);
					bfs::rename(secdOutFnp, pairedReader.ioOptions_.secondName_);
				}
			}
			//std::cout << "usedFilteredSinglesOpts: " << usedFilteredSinglesOpts.firstName_ <<  std::endl;
			if(singleReader.ioOptions_.inExists()){
				singleReader.openIn();
				seqInfo seq;
				SeqOutput keepWriter(outlierFilteredSingletOuts);
				SeqOutput filterWriter(outlierFilteredOff_singletOuts);
				while(singleReader.readNextRead(seq)){
					if(njh::in(seq.name_, keepNames)){
						keepWriter.openWrite(seq);
					}else{
						filterWriter.openWrite(seq);
					}
				}
				singleReader.closeIn();
				bfs::remove(singleReader.ioOptions_.firstName_);
				if(keepWriter.outOpen()){
					auto primOutFnp = keepWriter.getPrimaryOutFnp();
					keepWriter.closeOut();
					bfs::rename(primOutFnp, singleReader.ioOptions_.firstName_);
				}
			}
		}
	} else{
		ret.log_["bestKmerOccurenceCutOff"] = njh::json::toJson(extractionPars.kmerKOcurrenceCutOffs.front());
		ret.log_["bestKmerLength"] = njh::json::toJson(extractionPars.kmerLengths.front());
		ret.log_["message"] = njh::json::toJson("no parameters worked");
	}
	ret.log_["time"] = njh::json::toJson(watch);
	return ret;
}
}  // namespace njhseq
