/*
 * WayFindingPreFilteringUtils.cpp
 *
 *  Created on: Oct 17, 2019
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


#include "WayFindingPreFilteringUtils.hpp"


namespace njhseq {

std::string getPossibleSampleNameFromFnp(const bfs::path & fnp){
	std::string bName = bfs::basename(fnp);
	if(std::string::npos != bName.find(".")){
		bName = bName.substr(0, bName.find("."));
	}
	return bName;
}

std::string getPossibleSampleNameFromSeqName(const std::string & seqName, const std::string & defaultSamp){
	if(MetaDataInName::nameHasMetaData(seqName)){
		auto meta = MetaDataInName(seqName);
		if(meta.containsMeta("sample")){
			return meta.getMeta("sample");
		}
	}
	return defaultSamp;
}




preprocessSeqsForWayFindingRes preprocessSeqsForWayFinding(
		const preprocessSeqsForWayFindingPars & outPars,
		const BamExtractor::ExtractedFilesOpts & inOpts,
		const HaploPathFinder::PathFinderCorePars & extractionPars,
		const std::string & sampName){

	preprocessSeqsForWayFindingRes ret;
	SeqOutput writer(outPars.filteredPairedOpts);
	SeqOutput singleWriter(outPars.filteredSingletOuts);
	SeqOutput unused_writer(outPars.filteredOff_pairedOpts);
	SeqOutput unused_singleWriter(outPars.filteredOff_singletOuts);

	SeqOutput dup_writer(outPars.filteredOffDups_pairedOpts);
	SeqOutput dup_singleWriter(outPars.filteredOffDups_singletOuts);

	ReadCheckerQualCheck readQualChecker(extractionPars.qualCheck_, extractionPars.qualCheckCutOff_, false);
	/**@todo take the time to program a more efficient way of doing this*/
	std::unordered_map<std::string, std::unordered_set<std::string>> alreadyAddedSeqBySamplePair;
	std::unordered_map<std::string, std::unordered_set<std::string>> alreadyAddedSeqBySampleSingles;

	if (inOpts.inPairs_.inExists()) {
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << inOpts.inPairs_.firstName_ << std::endl;
//		std::cout << inOpts.inPairs_.secondName_ << std::endl;

		PairedRead pSeq;
		SeqInput pairedReader(inOpts.inPairs_);
		pairedReader.openIn();
		while (pairedReader.readNextRead(pSeq)) {
			readVec::getMaxLength(pSeq.seqBase_, ret.maxInputSeqLen);
			readVec::getMaxLength(pSeq.mateSeqBase_, ret.maxInputSeqLen);
			if(extractionPars.trimOnQual_){
				readVecTrimmer::trimAtLstripQualScore(pSeq.seqBase_, extractionPars.trimOnQual_);
				readVecTrimmer::trimAtRstripQualScore(pSeq.seqBase_, extractionPars.trimOnQual_);
				readVecTrimmer::trimAtLstripQualScore(pSeq.mateSeqBase_, extractionPars.trimOnQual_);
				readVecTrimmer::trimAtRstripQualScore(pSeq.mateSeqBase_, extractionPars.trimOnQual_);
			}

			//moved checking for Ns later
			//bool firstSkipped = std::string::npos != pSeq.seqBase_.seq_.find("N");
			//bool secondSkipped = std::string::npos != pSeq.mateSeqBase_.seq_.find("N");
			bool firstSkipped = false;
			bool secondSkipped = false;
			//quality filter, will turn read off if it doesn't pass checks
			if(!firstSkipped){
				readQualChecker.checkRead(pSeq.seqBase_);
				firstSkipped = !pSeq.seqBase_.on_;
			}

			if(!secondSkipped){
				readQualChecker.checkRead(pSeq.mateSeqBase_);
				secondSkipped = !pSeq.mateSeqBase_.on_;
			}
			if(firstSkipped){
				++ret.filteredR1Qual_;
			}
			if(secondSkipped){
				++ret.filteredR2Qual_;
			}
			bool dup = false;
			if(extractionPars.removeDuplicatedSequences_){
				if(!firstSkipped && !secondSkipped){
					std::string seqSampName = getPossibleSampleNameFromSeqName(pSeq.seqBase_.name_, sampName);
					VecStr seqs{pSeq.seqBase_.seq_, pSeq.mateSeqBase_.seq_};
					njh::sort(seqs);
					auto checkStr = njh::pasteAsStr(seqs.front(), "-", seqs.back());
					if (njh::in(checkStr, alreadyAddedSeqBySamplePair[seqSampName])) {
						dup = true;
						++ret.filteredPairsDups_;
					} else {
						alreadyAddedSeqBySamplePair[seqSampName].emplace(checkStr);
					}
				}
			}
			//check mate
			if(dup){
				dup_writer.openWrite(pSeq);
			}else	if (!firstSkipped && !secondSkipped) {
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
	if (inOpts.inPairsMateUnmapped_.inExists()) {
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << inOpts.inPairsMateUnmapped_.firstName_ << std::endl;
//		std::cout << inOpts.inPairsMateUnmapped_.secondName_ << std::endl;
		PairedRead pSeq;
		SeqInput pairedReader(inOpts.inPairsMateUnmapped_);
		pairedReader.openIn();
		while (pairedReader.readNextRead(pSeq)) {
			readVec::getMaxLength(pSeq.seqBase_, ret.maxInputSeqLen);
			readVec::getMaxLength(pSeq.mateSeqBase_, ret.maxInputSeqLen);
			if(extractionPars.trimOnQual_){
				readVecTrimmer::trimAtLstripQualScore(pSeq.seqBase_, extractionPars.trimOnQual_);
				readVecTrimmer::trimAtRstripQualScore(pSeq.seqBase_, extractionPars.trimOnQual_);
				readVecTrimmer::trimAtLstripQualScore(pSeq.mateSeqBase_, extractionPars.trimOnQual_);
				readVecTrimmer::trimAtRstripQualScore(pSeq.mateSeqBase_, extractionPars.trimOnQual_);
			}
			//moved checking for Ns later
			//bool firstSkipped = std::string::npos != pSeq.seqBase_.seq_.find("N");
			//bool secondSkipped = std::string::npos != pSeq.mateSeqBase_.seq_.find("N");
			bool firstSkipped = false;
			bool secondSkipped = false;
			//quality filter, will turn read off if it doesn't pass checks

			if(!firstSkipped){
				readQualChecker.checkRead(pSeq.seqBase_);
				firstSkipped = !pSeq.seqBase_.on_;
			}
			if(!secondSkipped){
				readQualChecker.checkRead(pSeq.mateSeqBase_);
				secondSkipped = !pSeq.mateSeqBase_.on_;
			}
			if(firstSkipped){
				++ret.filteredR1Qual_;
			}
			if(secondSkipped){
				++ret.filteredR2Qual_;
			}
			bool dup = false;
			if(extractionPars.removeDuplicatedSequences_){
				if(!firstSkipped && !secondSkipped){
					std::string seqSampName = getPossibleSampleNameFromSeqName(pSeq.seqBase_.name_, sampName);
					VecStr seqs{pSeq.seqBase_.seq_, pSeq.mateSeqBase_.seq_};
					njh::sort(seqs);
					auto checkStr = njh::pasteAsStr(seqs.front(), "-", seqs.back());
					if (njh::in(checkStr, alreadyAddedSeqBySamplePair[seqSampName])) {
						dup = true;
						++ret.filteredPairsDups_;
					} else {
						alreadyAddedSeqBySamplePair[seqSampName].emplace(checkStr);
					}
				}
			}
			//check mate
			if(dup){
				dup_writer.openWrite(pSeq);
			}else	if (!firstSkipped && !secondSkipped) {
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
	if (inOpts.inUnpaired_.inExists()) {
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << inOpts.inUnpaired_.firstName_ << std::endl;
		seqInfo seq;
		SeqInput singleReader(inOpts.inUnpaired_);
		singleReader.openIn();
		while (singleReader.readNextRead(seq)) {
			readVec::getMaxLength(seq, ret.maxInputSeqLen);
			if (extractionPars.trimOnQual_) {
				readVecTrimmer::trimAtLstripQualScore(seq,
						extractionPars.trimOnQual_);
				readVecTrimmer::trimAtRstripQualScore(seq,
						extractionPars.trimOnQual_);
			}
			//quality filter, will turn read off if it doesn't pass checks
			readQualChecker.checkRead(seq);
			if (seq.on_) {
			//if (seq.on_ && std::string::npos == seq.seq_.find("N")) { //moved checking for N's later
				if (extractionPars.removeDuplicatedSequences_) {
					std::string seqSampName = getPossibleSampleNameFromSeqName(seq.name_, sampName);
					if (njh::in(seq.seq_, alreadyAddedSeqBySampleSingles[seqSampName])) {
						dup_singleWriter.openWrite(seq);
						++ret.filteredSinglesDups_;
					} else {
						alreadyAddedSeqBySampleSingles[seqSampName].emplace(seq.seq_);
						singleWriter.openWrite(seq);
					}
				}else{
					singleWriter.openWrite(seq);
				}
			}else{
				++ret.filteredSinglesQual_;
				unused_singleWriter.openWrite(seq);
			}
		}
	}
	if (writer.outOpen()) {
		ret.usedFilteredPairedR1Fnp = writer.getPrimaryOutFnp();
		ret.usedFilteredPairedR2Fnp = writer.getSecondaryOutFnp();
	} else {
		writer.openOut();
		ret.usedFilteredPairedR1Fnp = writer.getPrimaryOutFnp();
		ret.usedFilteredPairedR2Fnp = writer.getSecondaryOutFnp();
		writer.closeOut();
		bfs::remove(ret.usedFilteredPairedR1Fnp);
		bfs::remove(ret.usedFilteredPairedR2Fnp);
	}
	if(singleWriter.outOpen()){
		ret.usedFilteredSinglesFnp = singleWriter.getPrimaryOutFnp();
	}else {
		singleWriter.openOut();
		ret.usedFilteredSinglesFnp = singleWriter.getPrimaryOutFnp();
		singleWriter.closeOut();
		bfs::remove(ret.usedFilteredSinglesFnp);
	}
	writer.closeOut();
	singleWriter.closeOut();
	unused_writer.closeOut();
	unused_singleWriter.closeOut();
	return ret;
}


repairSeqsAfterRecruitmentResults repairSeqsAfterRecruitment(const repairSeqsAfterRecruitmentPars & rePairPars){
	repairSeqsAfterRecruitmentResults ret;
	//re-pair reads
	if (bfs::exists(rePairPars.singlesFnp_) ) {
		auto singleReads = SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastqIn(rePairPars.singlesFnp_));
		if(len(singleReads) > 1){
			auto reWriteSeqOut = SeqIOOptions(SeqIOOptions::genFastqOut(rePairPars.singlesFnp_));
			reWriteSeqOut.out_.overWriteFile_ = true;
			std::regex r1Pat{"_R1.fastq$"};
			auto pairedOut = SeqIOOptions::genPairedOut(std::regex_replace(rePairPars.r1Fnp_.string(), r1Pat, ""));
			pairedOut.out_.append_ = true;
			SeqOutput rePairWriter(pairedOut);
			SeqOutput reWritten(reWriteSeqOut);
			readVecSorter::sortByName(singleReads, false);
			uint32_t pos = 1;
			uint32_t singlesWritten = 0;
			while (pos < len(singleReads)) {
				if (singleReads[pos - 1].name_ == singleReads[pos].name_) {
					PairedRead pSeq(singleReads[pos - 1], singleReads[pos]);
					rePairWriter.openWrite(pSeq);
					pos += 2;
					++ret.numberOfPairsRePaired_;
				} else {
					reWritten.openWrite(singleReads[pos - 1]);
					++singlesWritten;
					++pos;
					if (len(singleReads) == pos) {
						reWritten.openWrite(singleReads[pos - 1]);
					}
				}
			}
			if(0 == singlesWritten){
				bfs::remove(rePairPars.singlesFnp_);
			}
		}
	}
	return ret;
}


writeOutTandemsAndOptionallyStitchRes writeOutTandemsAndOptionallyStitch(
		const writeOutTandemsAndOptionallyStitchPars & runPars,
		const uint64_t maxInputSeqsLen,
		const HaploPathFinder::PathFinderCorePars & extractionPars
		){
	writeOutTandemsAndOptionallyStitchRes empty;
	auto currentRes = writeOutTandemsAndOptionallyStitch(runPars, maxInputSeqsLen, extractionPars, empty);
	return currentRes;
}

writeOutTandemsAndOptionallyStitchRes writeOutTandemsAndOptionallyStitch(
		const writeOutTandemsAndOptionallyStitchPars & runPars,
		const uint64_t maxInputSeqsLen,
		const HaploPathFinder::PathFinderCorePars & extractionPars,
		const writeOutTandemsAndOptionallyStitchRes & previouslyDeterminedTandems) {
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	writeOutTandemsAndOptionallyStitchRes ret;
	uint32_t initialKmerLength = extractionPars.kmerLengths.front();
	if(!extractionPars.originalKmerLengths.empty()){
		initialKmerLength = extractionPars.originalKmerLengths.front();
	}
	SeqInput pairedReader(runPars.usedFilteredPairedOpts);
	SeqInput singleReader(runPars.usedFilteredSinglesOpts);

	std::mutex allTandemMut;

	if (pairedReader.ioOptions_.inExists()) {
		auto currentPairParams = extractionPars.pairProcessorParams_;
//		currentPairParams.debug_ = true;
//		currentPairParams.verbose_ = true;
		if(std::numeric_limits<uint32_t>::max() == currentPairParams.minOverlap_){
			currentPairParams.minOverlap_ = initialKmerLength;
		}
		concurrent::AlignerPool alnPool(maxInputSeqsLen,
				gapScoringParameters(10,1,0,0,0,0), substituteMatrix(2,-4), extractionPars.numThreads);
//			alnPool.inAlnDir_ = njh::files::make_path(workingDir, "trimAlnCache");
//			alnPool.outAlnDir_ = njh::files::make_path(workingDir, "trimAlnCache");

		alnPool.initAligners();

		pairedReader.openIn();




		MultiSeqIO multiWriter;
		multiWriter.addReader("stitched", SeqIOOptions::genFastqOut(njh::files::make_path(runPars.workingDir_, "stitched_filteredSingles")));
		multiWriter.addReader("not_stitched", SeqIOOptions::genPairedOut(njh::files::make_path(runPars.workingDir_, "not_stitched_filteredExtractedPairs")));


		PairedReadProcessor::ProcessedResultsCounts pairProcessCountsBeforeReOrientingAll;
		PairedReadProcessor::ProcessedResultsCounts pairProcessCountsAfterReOrientingAll;
		PairedReadProcessor::ProcessedResultsCounts pairProcessCountsCombined;

		std::function<void()> countTandemsPaired =
				[&extractionPars,&pairProcessCountsBeforeReOrientingAll,

				&pairProcessCountsAfterReOrientingAll,
				&allTandemMut,&pairedReader,&alnPool,
				&currentPairParams,&initialKmerLength,&multiWriter,
				&previouslyDeterminedTandems,&ret]() {
			PairedRead pSeq;
			std::vector<TandemRepeatPlusAdjustedSize> allCurrentTandems;
			auto alignerObj = alnPool.popAligner();
			PairedReadProcessor pProcessor(currentPairParams);

			PairedReadProcessor::ProcessedResultsCounts pairProcessCountsBeforeReOrienting;
			PairedReadProcessor::ProcessedResultsCounts pairProcessCountsAfterReOrienting;

			std::unordered_map<std::string, std::vector<TandemRepeatPlusAdjustedSize>> currentProcessedInfo;
			while (pairedReader.readNextReadLock(pSeq)) {
				std::vector<TandemRepeatPlusAdjustedSize> tandemsForPair;
				bool previouslyProcessed = false;
				//check to see if this pair has been previously processed
				if(!njh::in(pSeq.seqBase_.name_, previouslyDeterminedTandems.previouslyDeterminedTandemsForPairs_)){
					auto firstTandems =  aligner::findTandemRepeatsInSequence(pSeq.seqBase_.seq_, 2, -2, -7, initialKmerLength * 2);
					auto secondTandems = aligner::findTandemRepeatsInSequence(pSeq.mateSeqBase_.seq_,  2, -2, -7, initialKmerLength * 2);
					for (const auto & t : firstTandems) {
						uint32_t adjustedSize = t.getSize();
						if (t.startPos_ > 3 && t.startPos_ < t.repeat_.size()) {
							std::string begSeq = pSeq.seqBase_.seq_.substr(0, t.startPos_);
							std::string tPortionSeq = t.repeat_.substr(
									t.repeat_.size() - t.startPos_);
							if (numberOfMismatches(begSeq, tPortionSeq) <= 0) {
								adjustedSize += t.startPos_;
							}
						}
						if (pSeq.seqBase_.seq_.size() - t.stopPos_ > 3
								&& pSeq.seqBase_.seq_.size() - t.stopPos_ < t.repeat_.size()) {
							std::string backSeq = pSeq.seqBase_.seq_.substr(t.stopPos_);
							std::string tPortionSeq = t.repeat_.substr(0,
									pSeq.seqBase_.seq_.size() - t.stopPos_);
							if (numberOfMismatches(backSeq, tPortionSeq) <= 0) {
								adjustedSize += pSeq.seqBase_.seq_.size() - t.stopPos_;
							}
						}
						tandemsForPair.emplace_back(t, adjustedSize);
					}
					for (const auto & t : secondTandems) {
						//make sure not to double count tandems
						bool inFirstTandems = false;
						for(const auto & firstT : firstTandems){
							if(firstT.repeat_ == t.repeat_){
								inFirstTandems = true;
								break;
							}
						}
						if(inFirstTandems){
							continue;
						}
						uint32_t adjustedSize = t.getSize();
						if (t.startPos_ > 3 && t.startPos_ < t.repeat_.size()) {
							std::string begSeq = pSeq.mateSeqBase_.seq_.substr(0, t.startPos_);
							std::string tPortionSeq = t.repeat_.substr(
									t.repeat_.size() - t.startPos_);
							if (numberOfMismatches(begSeq, tPortionSeq) <= 0) {
								adjustedSize += t.startPos_;
							}
						}
						if (pSeq.mateSeqBase_.seq_.size() - t.stopPos_ > 3
								&& pSeq.mateSeqBase_.seq_.size() - t.stopPos_
										< t.repeat_.size()) {
							std::string backSeq = pSeq.mateSeqBase_.seq_.substr(t.stopPos_);
							std::string tPortionSeq = t.repeat_.substr(0,
									pSeq.mateSeqBase_.seq_.size() - t.stopPos_);
							if (numberOfMismatches(backSeq, tPortionSeq) <= 0) {
								adjustedSize += pSeq.mateSeqBase_.seq_.size() - t.stopPos_;
							}
						}
						tandemsForPair.emplace_back(t, adjustedSize);
					}
					currentProcessedInfo[pSeq.seqBase_.name_] = tandemsForPair;
				}else{
					previouslyProcessed = true;
					tandemsForPair = previouslyDeterminedTandems.previouslyDeterminedTandemsForPairs_.at(pSeq.seqBase_.name_);
				}

				bool longerThanKLen = false;
				for(const auto & t : tandemsForPair){
					if(t.adjustedPartialSize_ >= initialKmerLength){
						/**@todo should change this to minimum overlap size */
						longerThanKLen = true;
						break;
					}
				}
				bool addTandems = true;
				if(extractionPars.stitchPairs){
					/**@todo should change this making sure there's no tandem in the overlap*/
					if(!longerThanKLen && !previouslyProcessed){
						++pairProcessCountsBeforeReOrienting.total;
						auto pairRes = pProcessor.processPairedEnd(pSeq, pairProcessCountsBeforeReOrienting, *alignerObj);
						if(PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2 == pairRes.status_ && extractionPars.reOrientPairsForStitching){
							auto pSeqCopy = pSeq;
							pSeqCopy.seqBase_.reverseComplementRead(false, true);
							pSeqCopy.mateSeqBase_.reverseComplementRead(false, true);
							//auto previousLen = len(*pairRes.combinedSeq_);
							++pairProcessCountsAfterReOrienting.total;
							pairRes = pProcessor.processPairedEnd(pSeqCopy, pairProcessCountsAfterReOrienting, *alignerObj);
							if(nullptr != pairRes.combinedSeq_){
								pairRes.combinedSeq_->reverseComplementRead(false, true);
//								if(pairRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::NONE &&
//									 pairRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP &&
//									 pairRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2  ){
//								}
							}
						}
						if(!njh::in(pairRes.status_, extractionPars.acceptableOverlapStatuses) || (!extractionPars.reOrientPairsForStitching && pairRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP)){
//						if( pairRes.status_ == PairedReadProcessor::ReadPairOverLapStatus::NONE ||
//								pairRes.status_ == PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP ||
//								pairRes.status_ == PairedReadProcessor::ReadPairOverLapStatus::R1BEGINSINR2 ){
							multiWriter.openWrite("not_stitched", pSeq);
						}else{
							multiWriter.openWrite("stitched", pairRes.combinedSeq_);
							addTandems = false;
						}
					}else{
						multiWriter.openWrite("not_stitched", pSeq);
					}
				}
				if(addTandems){
					addOtherVec(allCurrentTandems, tandemsForPair);
				}
			}
			{
				std::lock_guard<std::mutex> tandemLock(allTandemMut);
				for(const auto & currentTandem : allCurrentTandems){
					if(std::string::npos == currentTandem.repeat_ .repeat_.find("N")){
						ret.allTandems.emplace_back(currentTandem);
					}
				}
				for(const auto & info : currentProcessedInfo){
					ret.previouslyDeterminedTandemsForPairs_.emplace(info);
				}
				pairProcessCountsBeforeReOrientingAll.addOther(pairProcessCountsBeforeReOrienting);
				pairProcessCountsAfterReOrientingAll.addOther(pairProcessCountsAfterReOrienting);
			}
		};
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		njh::concurrent::runVoidFunctionThreaded(countTandemsPaired, extractionPars.numThreads);
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		multiWriter.closeOutAll();
		if(extractionPars.stitchPairs){
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			OutputStream stitchResultsOut(njh::files::make_path(runPars.workingDir_, "stitchedResults.json"));

			if(extractionPars.reOrientPairsForStitching){
				pairProcessCountsCombined.total = pairProcessCountsBeforeReOrientingAll.total;
				pairProcessCountsCombined.perfectOverlapCombined += pairProcessCountsBeforeReOrientingAll.perfectOverlapCombined;
				pairProcessCountsCombined.r1EndsInR2Combined += pairProcessCountsBeforeReOrientingAll.r1EndsInR2Combined;
				pairProcessCountsCombined.perfectOverlapCombined += pairProcessCountsAfterReOrientingAll.perfectOverlapCombined;
				pairProcessCountsCombined.r1EndsInR2Combined += pairProcessCountsAfterReOrientingAll.r1EndsInR2Combined;
				pairProcessCountsCombined.overlapFail +=
						pairProcessCountsBeforeReOrientingAll.overlapFail +
						pairProcessCountsBeforeReOrientingAll.r1AllInR2Combined +
						pairProcessCountsBeforeReOrientingAll.r2AllInR1Combined;
				pairProcessCountsCombined.overlapFail +=
						pairProcessCountsAfterReOrientingAll.overlapFail +
						pairProcessCountsAfterReOrientingAll.r1AllInR2Combined +
						pairProcessCountsAfterReOrientingAll.r2AllInR1Combined +
						pairProcessCountsAfterReOrientingAll.r1BeginsInR2Combined;
				Json::Value stitchResults;
				stitchResults["beforeReorienting"] = pairProcessCountsBeforeReOrientingAll.toJsonCounts();
				stitchResults["afterReorienting"] = pairProcessCountsAfterReOrientingAll.toJsonCounts();
				stitchResults["combined"] = pairProcessCountsCombined.toJsonCounts();
				stitchResultsOut << stitchResults << std::endl;
			}else{
				stitchResultsOut << pairProcessCountsBeforeReOrientingAll.toJsonCounts() << std::endl;
			}
			auto notStitchedR1Fnp = njh::files::make_path(runPars.workingDir_,   "not_stitched_filteredExtractedPairs_R1.fastq");
			auto notStitchedR2Fnp = njh::files::make_path(runPars.workingDir_,   "not_stitched_filteredExtractedPairs_R2.fastq");
			auto stitchedSinglesFnp = njh::files::make_path(runPars.workingDir_, "stitched_filteredSingles.fastq");
			//if the stitched file exists then stitching happened
			if(bfs::exists(stitchedSinglesFnp)){
				bfs::remove(runPars.usedFilteredPairedOpts.firstName_);
				bfs::remove(runPars.usedFilteredPairedOpts.secondName_);
				if(bfs::exists(notStitchedR1Fnp)){
					bfs::rename(notStitchedR1Fnp, runPars.usedFilteredPairedOpts.firstName_);
					bfs::rename(notStitchedR2Fnp, runPars.usedFilteredPairedOpts.secondName_);
				}
				{
					OutOptions singlesOutOpts(runPars.usedFilteredSinglesOpts.firstName_);
					singlesOutOpts.append_ = true;
					OutputStream singlesOut(singlesOutOpts);
					InputStream stitchedSinglesIn(stitchedSinglesFnp);
					singlesOut << stitchedSinglesIn.rdbuf();
				}
				bfs::remove(stitchedSinglesFnp);
			}else if(bfs::exists(notStitchedR1Fnp)){
				bfs::remove(notStitchedR1Fnp);
				bfs::remove(notStitchedR2Fnp);
			}
//				bfs::remove(prePars.usedFilteredPairedR1Fnp);
//				bfs::remove(prePars.usedFilteredPairedR2Fnp);
//				if(bfs::exists(notStitchedR1Fnp)){
//					bfs::rename(notStitchedR1Fnp, prePars.usedFilteredPairedR1Fnp);
//					bfs::rename(notStitchedR2Fnp, prePars.usedFilteredPairedR2Fnp);
//				}
//				{
//					OutOptions singlesOutOpts(prePars.usedFilteredSinglesFnp);
//					singlesOutOpts.append_ = true;
//					OutputStream singlesOut(singlesOutOpts);
//					InputStream stitchedSinglesIn(stitchedSinglesFnp);
//					singlesOut << stitchedSinglesIn.rdbuf();
//				}
//				bfs::remove(stitchedSinglesFnp);
//			}else if(bfs::exists(notStitchedR1Fnp)){
//				bfs::remove(notStitchedR1Fnp);
//				bfs::remove(notStitchedR2Fnp);
//			}
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if (singleReader.ioOptions_.inExists()) {
		singleReader.openIn();
		std::function<void()> countTandemsSingles = [&ret,&allTandemMut,
																								 &singleReader,&initialKmerLength,
																								 &previouslyDeterminedTandems](){
			seqInfo seq;
			std::vector<TandemRepeatPlusAdjustedSize> allCurrentTandems;
			std::unordered_map<std::string, std::vector<TandemRepeatPlusAdjustedSize>> currentProcessedInfo;
			while (singleReader.readNextReadLock(seq)) {
				std::vector<TandemRepeatPlusAdjustedSize> tandemsForSingle;
				if(!njh::in(seq.name_, previouslyDeterminedTandems.previouslyDeterminedTandemsForSingles_)){
					auto firstTandems = aligner::findTandemRepeatsInSequence(seq.seq_, 2, -2, -7, initialKmerLength * 2);
					for (const auto & t : firstTandems) {
						uint32_t adjustedSize = t.getSize();
						if (t.startPos_ > 3 && t.startPos_ < t.repeat_.size()) {
							std::string begSeq = seq.seq_.substr(0, t.startPos_);
							std::string tPortionSeq = t.repeat_.substr(
									t.repeat_.size() - t.startPos_);
							if (numberOfMismatches(begSeq, tPortionSeq) <= 0) {
								adjustedSize += t.startPos_;
							}
						}
						if (seq.seq_.size() - t.stopPos_ > 3
								&& seq.seq_.size() - t.stopPos_ < t.repeat_.size()) {
							std::string backSeq = seq.seq_.substr(t.stopPos_);
							std::string tPortionSeq = t.repeat_.substr(0,
									seq.seq_.size() - t.stopPos_);
							if (numberOfMismatches(backSeq, tPortionSeq) <= 0) {
								adjustedSize += seq.seq_.size() - t.stopPos_;
							}
						}
						tandemsForSingle.emplace_back(t, adjustedSize);
					}
					currentProcessedInfo[seq.name_] = tandemsForSingle;
				} else {
					tandemsForSingle = previouslyDeterminedTandems.previouslyDeterminedTandemsForSingles_.at(seq.name_);
				}
				addOtherVec(allCurrentTandems, tandemsForSingle);
			}
			{
				std::lock_guard<std::mutex> tandemLock(allTandemMut);
				for(const auto & currentTandem : allCurrentTandems){
					if(std::string::npos == currentTandem.repeat_ .repeat_.find("N")){
						ret.allTandems.emplace_back(currentTandem);
					}
				}
				for(const auto & info : currentProcessedInfo){
					ret.previouslyDeterminedTandemsForSingles_.emplace(info);
				}
			}
		};
		njh::concurrent::runVoidFunctionThreaded(countTandemsSingles,extractionPars.numThreads);
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

	{
		njh::sort(ret.allTandems, [](const TandemRepeatPlusAdjustedSize & tr1, const TandemRepeatPlusAdjustedSize & tr2){
			if(tr1.repeat_.getSize() == tr2.repeat_.getSize()){
				if(tr1.repeat_.repeat_ == tr2.repeat_.repeat_){
					return tr1.adjustedPartialSize_ < tr2.adjustedPartialSize_;
				}else{
					return tr1.repeat_.repeat_ < tr2.repeat_.repeat_;
				}
			}else{
				return tr1.repeat_.getSize() < tr2.repeat_.getSize();
			}
		});
		OutputStream outTandemInfoOut(runPars.tandemInfoOpts_);
		outTandemInfoOut << "repeat\tsize\tadjustedPartialSize" << std::endl;
		for (const auto & t : ret.allTandems) {
			outTandemInfoOut << t.repeat_.repeat_
					<< "\t" << t.repeat_.getSize()
					<< "\t" << t.adjustedPartialSize_ << "\n";
		}
	}
	{
		//add in the previously added tandem information
		for(const auto & info : previouslyDeterminedTandems.previouslyDeterminedTandemsForPairs_){
			ret.previouslyDeterminedTandemsForPairs_.emplace(info);
		}
		for(const auto & info : previouslyDeterminedTandems.previouslyDeterminedTandemsForSingles_){
			ret.previouslyDeterminedTandemsForSingles_.emplace(info);
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	return ret;
}






}  // namespace njhseq