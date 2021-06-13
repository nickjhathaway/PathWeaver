/*
 * RegionInvestigatorInBam.cpp
 *
 *  Created on: Feb 7, 2020
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



#include "RegionInvestigatorInBam.hpp"
#include <njhseq/concurrency/pools/BamReaderPool.hpp>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecTrimmer.hpp>

namespace njhseq {


	BamRegionInvestigator::BamCountSpecficRegionsPars::BamCountSpecficRegionsPars(){
		pairPars.minOverlap_ = 8;
	}

	void BamRegionInvestigator::BamCountSpecficRegionsPars::setDefaults(seqSetUp & setUp){
		setUp.setOption(forcePlusStrand, "--forcePlusStrand", "force Plus Strand orientation, otherwise use the strand orientation in bed file");
		setUp.setOption(baseQuality, "--baseQuality", "Base Quality Cut Off");
		setUp.setOption(mappingQuality, "--mappingQuality", "Mapping Quality Cut Off");
		setUp.setOption(matchIDCutOff, "--matchIDCutOff", "Match ID Cut Off");
		setUp.setOption(numThreads, "--numThreads", "Number Threads");
		setUp.setOption(pairPars.minOverlap_, "--minOverlap", "Min Overlap in overlap when stitching pairs");
		setUp.setOption(pairPars.errorAllowed_, "--errorAllowed", "Error allowed in overlap when stitching pairs");
		setUp.setOption(pairPars.hardMismatchCutOff_, "--hardMismatchCutOff", "Hard Mismatch Cut Off allowed in overlaop when stitching pairs");
		setUp.setOption(pairPars.lqMismatchCutOff, "--lqMismatchCutOff", "low quality Mismatch Cut Off allowed in overlaop when stitching pairs");
		setUp.setOption(pairPars.hqMismatchCutOff, "--hqMismatchCutOff", "hq Mismatch Cut Off allowed in overlaop when stitching pairs");
		setUp.setOption(extendAmount, "--extendAmount", "extend this amount around the region to get better trimming");
		setUp.setOption(countDuplicates, "--countDuplicates", "Skip reads marked as duplicates");
		setUp.setOption(totalCountCutOff, "--totalCountCutOff", "Total Count Cut Off");
		setUp.setOption(perBaseCountCutOff, "--perBaseCountCutOff", "Per Base Count Cut Off");

	}


BamRegionInvestigator::RegionInfo::RegionInfo(
		const GenomicRegion & region) :

		region_(region) {
}





double BamRegionInvestigator::RegionInfo::getPerBaseCov() const {
	return static_cast<double>(coverage_) / region_.getLen();
}

double BamRegionInvestigator::RegionInfo::getProperPairFrac() const {
	return totalPairedReads_ == 0 ? 0.0 : totalProperPairedReads_ / static_cast<double>(totalPairedReads_);
}
double BamRegionInvestigator::RegionInfo::getProperPairMateUnmappedFrac() const {
	return totalPairedReads_ == 0 ? 0.0 : (totalProperPairedReads_ + totalOneMateUnmappedImproperPairs_) / static_cast<double>(totalPairedReads_);
}




BamRegionInvestigator::BamRegionInvestigator(const BamRegionInvestigatorPars & pars) :
		pars_(pars) {

}

std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> BamRegionInvestigator::getCoverageOnFromBam(
		const bfs::path & bamFnp,
		const std::vector<GenomicRegion> & regions) const {
	std::vector<std::shared_ptr<RegionInfo>> pairs;
	std::mutex pairsMut;
	njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);
	std::function<void()> getCov = [&bamFnp,&pairs,&pairsMut,
																	&regionsQueue,this](){
		std::vector<std::shared_ptr<RegionInfo>> currentBamRegionsPairs;
		GenomicRegion region;
		BamTools::BamReader bReader;
		bReader.Open(bamFnp.string());
		loadBamIndexThrow(bReader, __PRETTY_FUNCTION__);
		checkBamOpenThrow(bReader, bamFnp);
		BamTools::BamAlignment bAln;
		auto refData = bReader.GetReferenceData();
		while(regionsQueue.getVal(region)){
			currentBamRegionsPairs.emplace_back(std::make_shared<RegionInfo>(getCoverageForRegion(bReader, region)));
		}
		{
			std::lock_guard<std::mutex> lock(pairsMut);
			addOtherVec(pairs, currentBamRegionsPairs);
		}
	};

	njh::concurrent::runVoidFunctionThreaded(getCov, pars_.numThreads_);
	//sort
	njh::sort(pairs,[](const std::shared_ptr<RegionInfo> & p1, const std::shared_ptr<RegionInfo> & p2){
		if(p1->region_.chrom_ == p2->region_.chrom_) {
			if(p1->region_.start_ == p2->region_.start_) {
				return p1->region_.end_ < p2->region_.end_;
			} else {
				return p1->region_.start_ < p2->region_.start_;
			}
		} else {
			return p1->region_.chrom_ < p2->region_.chrom_;
		}
	});
	return pairs;
}

BamRegionInvestigator::RegionInfo BamRegionInvestigator::getCoverageForRegion(
		BamTools::BamReader & bReader, const GenomicRegion & region) const {
	RegionInfo ret(region);
	setBamFileRegionThrow(bReader, region);
	BamAlnsCache bCache;
	BamTools::BamAlignment bAln;
	auto refData = bReader.GetReferenceData();
	while(bReader.GetNextAlignmentCore(bAln)){
		if(bAln.IsPrimaryAlignment() && bAln.IsPaired()){
				++ret.totalPairedReads_;
			if(bAln.IsProperPair()){
				++ret.totalProperPairedReads_;
			}else if(!bAln.IsProperPair() && ((bAln.IsMapped() && !bAln.IsMateMapped()) || (!bAln.IsMapped() && bAln.IsMateMapped())) ){
				++ret.totalOneMateUnmappedImproperPairs_;
			}
		}
		if(bAln.IsMapped() && bAln.IsPrimaryAlignment()){
			if(bAln.IsDuplicate() && !pars_.countDups_){
				continue;
			}
			if(bAln.MapQuality <  pars_.mapQualityCutOff_){
				continue;
			}
			//only try to find and adjust for the mate's overlap if is mapped and could possibly fall within the region
			if(bAln.IsPaired() && bAln.IsMateMapped() && bAln.MatePosition < ret.region_.end_){
				bAln.BuildCharData();
				if(bCache.has(bAln.Name)){
					//get and adjust current read region
					GenomicRegion alnRegion(bAln, refData);
					alnRegion.start_ = std::max(alnRegion.start_, ret.region_.start_);
					alnRegion.end_ = std::min(alnRegion.end_,     ret.region_.end_);
					//get and adjust the mate's read region
					auto mate = bCache.get(bAln.Name);
					GenomicRegion mateRegion(*mate, refData);
					mateRegion.start_ = std::max(mateRegion.start_, ret.region_.start_);
					mateRegion.end_ = std::min(mateRegion.end_,     ret.region_.end_);
					//now that the two regions have been adjusted to be just what overlaps the current region
					//coverage will be the length of the two regions minus the overlap between the two
					auto cov = alnRegion.getLen() + mateRegion.getLen() - alnRegion.getOverlapLen(mateRegion);
					ret.coverage_ += cov;
					if(bAln.MapQuality < pars_.mapQualityCutOffForMultiMap_ || mate->MapQuality < pars_.mapQualityCutOffForMultiMap_){
						ret.multiMapCoverage_ += cov;
					}
					bCache.remove(bAln.Name);
				}else{
					bCache.add(bAln);
				}
			} else {
				/**@todo this doesn't take into account gaps, so base coverage isn't precise right here, more of an approximation, should improve */
				GenomicRegion alnRegion(bAln, refData);
				ret.coverage_ += ret.region_.getOverlapLen(alnRegion);
				if(bAln.MapQuality < pars_.mapQualityCutOffForMultiMap_ ){
					ret.multiMapCoverage_ += ret.region_.getOverlapLen(alnRegion);
				}
			}
		}
	}
	//save alignments where mate couldn't be found, this could be for several reasons like mate was dup or mate's mapping quality or other filtering reasons
	for(const auto & name : bCache.getNames()){
		auto aln = bCache.get(name);
		GenomicRegion alnRegion(*aln, refData);
		ret.coverage_ += ret.region_.getOverlapLen(alnRegion);
		if(aln->MapQuality < pars_.mapQualityCutOffForMultiMap_ ){
			ret.multiMapCoverage_ += ret.region_.getOverlapLen(alnRegion);
		}
	}
	return ret;
}



BamRegionInvestigator::RegionInfo BamRegionInvestigator::getCoverageAndFullySpanningForRegion(
		BamTools::BamReader &bReader, const GenomicRegion &region,
		const seqInfo &refSeq, const BamCountSpecficRegionsPars &spanningReadsPar,
		aligner &alignerObj) const {
	RegionInfo ret(getCoverageForRegion(bReader, region));
	updateFullSpanningCountForRegion(ret, refSeq, spanningReadsPar, bReader,
			alignerObj);
	return ret;
}



BamRegionInvestigator::ReadCountsRes BamRegionInvestigator::getTotalAndFullySpanningReadCounts(
		const bfs::path & bamFnp, const std::vector<GenomicRegion> & regions,
		const BamCountSpecficRegionsPars & spanningReadsPar,
		const bfs::path & genomeFnpInput) const {

	ReadCountsRes ret;
	njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);
	njhseq::concurrent::BamReaderPool bamPool(bamFnp, spanningReadsPar.numThreads);
	bamPool.openBamFile();


	std::mutex spanRCountsMut;
	auto genomeFnp2bit = genomeFnpInput;
	if(!njh::endsWith(genomeFnp2bit.string(), ".2bit")){
		if(njh::endsWith(genomeFnp2bit.string(), ".fasta") || njh::endsWith(genomeFnp2bit.string(), ".fna")){
			genomeFnp2bit.replace_extension(".2bit");
		}
	}
	TwoBit::TwoBitFile topRReader(genomeFnp2bit);
	std::unordered_map<std::string, seqInfo> regionSeqs;
	uint64_t maxlenForRegions = 0;
	auto chromLengths = topRReader.getSeqLens();
	for(const auto & reg : regions){
		auto regCopy = reg;
		BedUtility::extendLeftRight(regCopy, spanningReadsPar.extendAmount, spanningReadsPar.extendAmount, njh::mapAt(chromLengths, regCopy.chrom_));
		regionSeqs[reg.createUidFromCoordsStrand()] = regCopy.extractSeq(topRReader);
		readVec::getMaxLength(regionSeqs[reg.createUidFromCoordsStrand()], maxlenForRegions);
	}

	std::function<void()> extractReadsForRegion = [&regionsQueue,
																&bamPool,&spanningReadsPar,
																&ret,
																&spanRCountsMut,
																&regionSeqs,&maxlenForRegions,this
																](){

		//it's unlikely there will reads that are longer 1000 for now, so work around from trying to request too much memory for aligners
		aligner alignerObj(std::min<uint64_t>(std::max<uint64_t>(maxlenForRegions, 500), 1000), gapScoringParameters(5,1,0,0,0,0));
		GenomicRegion currentRegion;
		auto currentBReader = bamPool.popReader();
		auto refData = currentBReader->GetReferenceData();
		BamTools::BamAlignment bAln;
		std::unordered_map<std::string, uint32_t> currentSpanningReadCounts;
		std::unordered_map<std::string, uint32_t> currentTotalReadCounts;
		while(regionsQueue.getVal(currentRegion)){
			const auto & refSeq = regionSeqs.at(currentRegion.createUidFromCoordsStrand());
			auto currentCount = getFullSpanningCountForRegion(currentRegion, refSeq, spanningReadsPar, *currentBReader, alignerObj);
			currentSpanningReadCounts[currentRegion.createUidFromCoordsStrand()]+=currentCount.spanningReadCount_;
			currentTotalReadCounts[currentRegion.createUidFromCoordsStrand()]+=currentCount.totalReadCount_;;
		}
		{
			std::lock_guard<std::mutex> lock(spanRCountsMut);
			for(const auto & spanCount : currentSpanningReadCounts){
				ret.spanningReadCounts_[spanCount.first] = spanCount.second;
			}

			for(const auto & totalCount : currentTotalReadCounts){
				ret.totalReadCounts_[totalCount.first] = totalCount.second;
			}
		}
	};
	njh::concurrent::runVoidFunctionThreaded(extractReadsForRegion, spanningReadsPar.numThreads);
	return ret;
}


void BamRegionInvestigator::getTotalAndFullySpanningReadCounts(const bfs::path & bamFnp,
		const std::vector<GenomicRegion> & regions,
		const BamCountSpecficRegionsPars & spanningReadsPar,
		const bfs::path & genomeFnp,
		std::vector<std::shared_ptr<RegionInfo>> & pairs) const{
	std::unordered_map<std::string, std::shared_ptr<RegionInfo>> pairsByUID;
	for(const auto & p : pairs){
		pairsByUID[p->region_.createUidFromCoordsStrand()] = p;
	}
	getTotalAndFullySpanningReadCounts(bamFnp, regions, spanningReadsPar, genomeFnp, pairsByUID);
}


void BamRegionInvestigator::getTotalAndFullySpanningReadCounts(const bfs::path & bamFnp,
		const std::vector<GenomicRegion> & regions,
		const BamCountSpecficRegionsPars & spanningReadsPar,
		const bfs::path & genomeFnp,
		std::unordered_map<std::string, std::shared_ptr<RegionInfo>> & pairsByUid) const{
	auto counts = getTotalAndFullySpanningReadCounts(bamFnp, regions,
			spanningReadsPar, genomeFnp);
	for (const auto & totalCount : counts.totalReadCounts_) {
		njh::mapAt(pairsByUid, totalCount.first)->totalReads_ = totalCount.second;
	}
	for (const auto & totalSpanCount : counts.spanningReadCounts_) {
		njh::mapAt(pairsByUid, totalSpanCount.first)->totalFullySpanningReads_ =
				totalSpanCount.second;
	}
}

void BamRegionInvestigator::updateFullSpanningCountForRegion(
		RegionInfo &currentInfo, const seqInfo &refSeq,
		const BamCountSpecficRegionsPars &spanningReadsPar,
		BamTools::BamReader &currentBReader, aligner &alignerObj) const {

	readVecTrimmer::GlobalAlnTrimPars trimPars;
	trimPars.needJustOneEnd_ = false;
	trimPars.startInclusive_ = spanningReadsPar.extendAmount;
	trimPars.endInclusive_ = len(refSeq) - 1 - spanningReadsPar.extendAmount;

	BamAlnsCache cache;
	BamTools::BamAlignment bAln;
	auto refData = currentBReader.GetReferenceData();
	uint64_t maxlen = 0;
	PairedReadProcessor pProcess(spanningReadsPar.pairPars);
	PairedReadProcessor::ProcessedResultsCounts processCounts;
	setBamFileRegionThrow(currentBReader, currentInfo.region_);
	while(currentBReader.GetNextAlignment(bAln)){
		if(bAln.IsPrimaryAlignment() && bAln.IsMapped()){
			if(bAln.IsPaired() && bAln.IsMateMapped() && bAln.MatePosition < currentInfo.region_.end_){
				if(bAln.Position < bAln.MatePosition){
					cache.add(bAln);
				} else {
					if(cache.has(bAln.Name)) {
						++currentInfo.totalReads_;
						auto mate = cache.get(bAln.Name);
						//see if stitching is even plausible
						if(mate->GetEndPosition() >  bAln.Position){
							if(mate->Position <= currentInfo.region_.start_ &&
								 bAln.GetEndPosition() >= currentInfo.region_.end_ &&
								 mate->IsReverseStrand() != bAln.IsReverseStrand()){
								//stitching plausible and spanning
								//spanning read
								seqInfo firstMate;
								seqInfo secondMate;
								bool reverseCompFirstMate = false;
								if(bAln.IsFirstMate()){
									firstMate = bamAlnToSeqInfo(bAln, false);
									secondMate = bamAlnToSeqInfo(*mate, false);
									secondMate.reverseComplementRead(false, true);
									reverseCompFirstMate = bAln.IsReverseStrand();
								}else{
									firstMate = bamAlnToSeqInfo(*mate, false);
									secondMate = bamAlnToSeqInfo(bAln, false);
									secondMate.reverseComplementRead(false, true);
									reverseCompFirstMate = mate->IsReverseStrand();
								}
								readVec::getMaxLength(firstMate, maxlen);
								readVec::getMaxLength(secondMate, maxlen);
								PairedRead pseq(firstMate, secondMate);
								pseq.mateRComplemented_ = false;
								auto pairRes = pProcess.processPairedEnd(pseq, processCounts, alignerObj);
								if(nullptr != pairRes.combinedSeq_){
									readVec::getMaxLength(*pairRes.combinedSeq_, maxlen);
									alignerObj.parts_.setMaxSize(maxlen);
									if(reverseCompFirstMate != currentInfo.region_.reverseSrand_){
										pairRes.combinedSeq_->reverseComplementRead(false, true);
									}
									readVecTrimmer::trimSeqToRefByGlobalAln(*pairRes.combinedSeq_, refSeq,trimPars, alignerObj);
									if(pairRes.combinedSeq_->on_){
										//spanning read
										++currentInfo.totalFullySpanningReads_;
									}
								}else if(bAln.Position <= currentInfo.region_.start_ && bAln.GetEndPosition() >= currentInfo.region_.end_){
									//spanning read
									seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
									if(bAln.IsReverseStrand() != currentInfo.region_.reverseSrand_){
										querySeq.reverseComplementRead(false, true);
									}
									readVec::getMaxLength(querySeq, maxlen);
									alignerObj.parts_.setMaxSize(maxlen);
									readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
									if(querySeq.on_){
										//spanning read
										++currentInfo.totalFullySpanningReads_;
									}
								}else if(mate->Position <= currentInfo.region_.start_ && mate->GetEndPosition() >= currentInfo.region_.end_){
									//spanning read
									seqInfo querySeq = bamAlnToSeqInfo(*mate, false);
									if(mate->IsReverseStrand() != currentInfo.region_.reverseSrand_){
										querySeq.reverseComplementRead(false, true);
									}
									readVec::getMaxLength(querySeq, maxlen);
									alignerObj.parts_.setMaxSize(maxlen);
									readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
									if(querySeq.on_){
										//spanning read
										++currentInfo.totalFullySpanningReads_;
									}
								}
							}
						}
						cache.remove(bAln.Name);
					} else {
						++currentInfo.totalReads_;
						//doesn't have mate, it's mate probably doesn't map to this region
						if(bAln.Position <= currentInfo.region_.start_ && bAln.GetEndPosition() >= currentInfo.region_.end_){
							//spanning read
							seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
							if(bAln.IsReverseStrand() != currentInfo.region_.reverseSrand_){
								querySeq.reverseComplementRead(false, true);
							}
							readVec::getMaxLength(querySeq, maxlen);
							alignerObj.parts_.setMaxSize(maxlen);
							readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
							if(querySeq.on_){
								//spanning read
								++currentInfo.totalFullySpanningReads_;
							}
						}
					}
				}
			}else{
				++currentInfo.totalReads_;
				if(bAln.Position <= currentInfo.region_.start_ && bAln.GetEndPosition() >= currentInfo.region_.end_){
					//spanning read
					seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
					if(bAln.IsReverseStrand() != currentInfo.region_.reverseSrand_){
						querySeq.reverseComplementRead(false, true);
					}
					readVec::getMaxLength(querySeq, maxlen);
					alignerObj.parts_.setMaxSize(maxlen);
					readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq, trimPars, alignerObj);
					if(querySeq.on_){
						//spanning read
						++currentInfo.totalFullySpanningReads_;
					}
				}
			}
		}
	}
}

BamRegionInvestigator::ReadCountsResSingle BamRegionInvestigator::getFullSpanningCountForRegion(
		const GenomicRegion & currentRegion, const seqInfo & refSeq,
		const BamCountSpecficRegionsPars & spanningReadsPar,
		BamTools::BamReader &currentBReader, aligner & alignerObj) const {

	BamRegionInvestigator::ReadCountsResSingle ret;

	readVecTrimmer::GlobalAlnTrimPars trimPars;
	trimPars.needJustOneEnd_ = false;
	trimPars.startInclusive_ = spanningReadsPar.extendAmount;
	trimPars.endInclusive_ = len(refSeq) - 1 - spanningReadsPar.extendAmount;

	BamAlnsCache cache;
	BamTools::BamAlignment bAln;
	auto refData = currentBReader.GetReferenceData();
	uint64_t maxlen = 0;
	PairedReadProcessor pProcess(spanningReadsPar.pairPars);
	PairedReadProcessor::ProcessedResultsCounts processCounts;
	setBamFileRegionThrow(currentBReader, currentRegion);
	while(currentBReader.GetNextAlignment(bAln)){
		if(bAln.IsPrimaryAlignment() && bAln.IsMapped()){
			if(bAln.IsPaired() && bAln.IsMateMapped() && bAln.MatePosition < currentRegion.end_){
				if(bAln.Position < bAln.MatePosition){
					cache.add(bAln);
				} else {
					if(cache.has(bAln.Name)) {
						++ret.totalReadCount_;
						auto mate = cache.get(bAln.Name);
						//see if stitching is even plausible
						if(mate->GetEndPosition() >  bAln.Position){
							if(mate->Position <= currentRegion.start_ &&
								 bAln.GetEndPosition() >= currentRegion.end_ &&
								 mate->IsReverseStrand() != bAln.IsReverseStrand()){
								//stitching plausible and spanning
								//spanning read
								seqInfo firstMate;
								seqInfo secondMate;
								bool reverseCompFirstMate = false;
								if(bAln.IsFirstMate()){
									firstMate = bamAlnToSeqInfo(bAln, false);
									secondMate = bamAlnToSeqInfo(*mate, false);
									secondMate.reverseComplementRead(false, true);
									reverseCompFirstMate = bAln.IsReverseStrand();
								}else{
									firstMate = bamAlnToSeqInfo(*mate, false);
									secondMate = bamAlnToSeqInfo(bAln, false);
									secondMate.reverseComplementRead(false, true);
									reverseCompFirstMate = mate->IsReverseStrand();
								}
								readVec::getMaxLength(firstMate, maxlen);
								readVec::getMaxLength(secondMate, maxlen);
								PairedRead pseq(firstMate, secondMate);
								pseq.mateRComplemented_ = false;
								auto pairRes = pProcess.processPairedEnd(pseq, processCounts, alignerObj);
								if(nullptr != pairRes.combinedSeq_){
									readVec::getMaxLength(*pairRes.combinedSeq_, maxlen);
									alignerObj.parts_.setMaxSize(maxlen);
									if(reverseCompFirstMate != currentRegion.reverseSrand_){
										pairRes.combinedSeq_->reverseComplementRead(false, true);
									}
									readVecTrimmer::trimSeqToRefByGlobalAln(*pairRes.combinedSeq_, refSeq,trimPars, alignerObj);
									if(pairRes.combinedSeq_->on_){
										//spanning read
										++ret.spanningReadCount_;
									}
								}else if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
									//spanning read
									seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
									if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
										querySeq.reverseComplementRead(false, true);
									}
									readVec::getMaxLength(querySeq, maxlen);
									alignerObj.parts_.setMaxSize(maxlen);
									readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
									if(querySeq.on_){
										//spanning read
										++ret.spanningReadCount_;
									}
								}else if(mate->Position <= currentRegion.start_ && mate->GetEndPosition() >= currentRegion.end_){
									//spanning read
									seqInfo querySeq = bamAlnToSeqInfo(*mate, false);
									if(mate->IsReverseStrand() != currentRegion.reverseSrand_){
										querySeq.reverseComplementRead(false, true);
									}
									readVec::getMaxLength(querySeq, maxlen);
									alignerObj.parts_.setMaxSize(maxlen);
									readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
									if(querySeq.on_){
										//spanning read
										++ret.spanningReadCount_;
									}
								}
							}
						}
						cache.remove(bAln.Name);
					} else {
						++ret.totalReadCount_;
						//doesn't have mate, it's mate probably doesn't map to this region
						if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
							//spanning read
							seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
							if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
								querySeq.reverseComplementRead(false, true);
							}
							readVec::getMaxLength(querySeq, maxlen);
							alignerObj.parts_.setMaxSize(maxlen);
							readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq,trimPars, alignerObj);
							if(querySeq.on_){
								//spanning read
								++ret.spanningReadCount_;
							}
						}
					}
				}
			}else{
				++ret.totalReadCount_;
				if(bAln.Position <= currentRegion.start_ && bAln.GetEndPosition() >= currentRegion.end_){
					//spanning read
					seqInfo querySeq = bamAlnToSeqInfo(bAln, false);
					if(bAln.IsReverseStrand() != currentRegion.reverseSrand_){
						querySeq.reverseComplementRead(false, true);
					}
					readVec::getMaxLength(querySeq, maxlen);
					alignerObj.parts_.setMaxSize(maxlen);
					readVecTrimmer::trimSeqToRefByGlobalAln(querySeq, refSeq, trimPars, alignerObj);
					if(querySeq.on_){
						//spanning read
						++ret.spanningReadCount_;
					}
				}
			}
		}
	}
	return ret;
}






std::vector<std::shared_ptr<BamRegionInvestigator::RegionInfo>> BamRegionInvestigator::getCoverageAndFullSpanningReads(
		const bfs::path & bamFnp,
		const std::vector<GenomicRegion> & regions,
		const BamCountSpecficRegionsPars & spanningReadsPar,
		const bfs::path & genomeFnp) const{

	std::vector<std::shared_ptr<RegionInfo>> pairs;
	std::mutex pairsMut;
	njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);

	auto genomeFnp2bit = genomeFnp;
	if(!njh::endsWith(genomeFnp2bit.string(), ".2bit")){
		if(njh::endsWith(genomeFnp2bit.string(), ".fasta") || njh::endsWith(genomeFnp2bit.string(), ".fna")){
			genomeFnp2bit.replace_extension(".2bit");
		}
	}
	TwoBit::TwoBitFile topRReader(genomeFnp2bit);
	std::unordered_map<std::string, seqInfo> regionSeqs;
	uint64_t maxlenForRegions = 0;
	auto chromLengths = topRReader.getSeqLens();
	for(const auto & reg : regions){
		auto regCopy = reg;
		BedUtility::extendLeftRight(regCopy, spanningReadsPar.extendAmount, spanningReadsPar.extendAmount, njh::mapAt(chromLengths, regCopy.chrom_));
		regionSeqs[reg.createUidFromCoordsStrand()] = regCopy.extractSeq(topRReader);
		readVec::getMaxLength(regionSeqs[reg.createUidFromCoordsStrand()], maxlenForRegions);
	}

	std::function<void()> getCov = [&bamFnp,&pairs,&pairsMut,&maxlenForRegions,
																	&regionSeqs,&spanningReadsPar,
																	&regionsQueue,this](){
		std::vector<std::shared_ptr<RegionInfo>> currentBamRegionsPairs;
		GenomicRegion region;
		BamTools::BamReader bReader;
		bReader.Open(bamFnp.string());
		loadBamIndexThrow(bReader, __PRETTY_FUNCTION__);
		checkBamOpenThrow(bReader, bamFnp);

		aligner alignerObj(std::min<uint64_t>(std::max<uint64_t>(maxlenForRegions, 500), 1000), gapScoringParameters(5,1,0,0,0,0));

		BamTools::BamAlignment bAln;
		auto refData = bReader.GetReferenceData();
		while(regionsQueue.getVal(region)){
			const auto & refSeq = regionSeqs.at(region.createUidFromCoordsStrand());
			currentBamRegionsPairs.emplace_back(std::make_shared<RegionInfo>(getCoverageAndFullySpanningForRegion(bReader, region, refSeq, spanningReadsPar, alignerObj)));
		}
		{
			std::lock_guard<std::mutex> lock(pairsMut);
			addOtherVec(pairs, currentBamRegionsPairs);
		}
	};

	njh::concurrent::runVoidFunctionThreaded(getCov, pars_.numThreads_);
	//sort
	njh::sort(pairs,[](const std::shared_ptr<RegionInfo> & p1, const std::shared_ptr<RegionInfo> & p2){
		if(p1->region_.chrom_ == p2->region_.chrom_) {
			if(p1->region_.start_ == p2->region_.start_) {
				return p1->region_.end_ < p2->region_.end_;
			} else {
				return p1->region_.start_ < p2->region_.start_;
			}
		} else {
			return p1->region_.chrom_ < p2->region_.chrom_;
		}
	});
	return pairs;
}


void BamRegionInvestigator::writeCovInfo(
		const std::vector<std::shared_ptr<RegionInfo>> & pairs,
		const std::string & sampName, const OutOptions & outOpts) const {
	OutputStream regionInfoOut(outOpts);
	regionInfoOut << "#chrom\tstart\tend\tname\tlength\tstrand\tperBaseCoverage";
	if(pars_.mapQualityCutOffForMultiMap_  > pars_.mapQualityCutOff_){
		regionInfoOut << "\tperBaseCoverageFromMapQOf<=" << pars_.mapQualityCutOffForMultiMap_;
	}
	regionInfoOut << "\t#OfPairedReads\tproperPairs\timproperMateUnmapped\tproperPairFrac\tproperPairAndMateUnmappedFrac";

	regionInfoOut << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : pairs){

		auto bedOut = p->region_.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		regionInfoOut << "\textraField"<<t;
	}
	regionInfoOut << "\n";
	for(const auto & p : pairs){
		auto bedOut = p->region_.genBedRecordCore();
		regionInfoOut << bedOut.toDelimStr();
		regionInfoOut << "\t" << p->coverage_/static_cast<double>(p->region_.getLen());
		if(pars_.mapQualityCutOffForMultiMap_  > pars_.mapQualityCutOff_){
			regionInfoOut << "\t" << p->multiMapCoverage_/static_cast<double>(p->region_.getLen());
		}
		regionInfoOut << "\t" << p->totalPairedReads_
				<< '\t' << p->totalProperPairedReads_
				<< '\t' << p->totalOneMateUnmappedImproperPairs_
				<< '\t' << p->getProperPairFrac()
				<< '\t' << p->getProperPairMateUnmappedFrac();
		regionInfoOut << "\t" << sampName;
		for(const auto & extra : bedOut.extraFields_){
			regionInfoOut << "\t" << extra;
		}
		uint32_t leftOverExtra = maxExtraFields - bedOut.extraFields_.size();
		for(uint32_t t = 0; t < leftOverExtra; ++t){
			regionInfoOut << "\t";
		}
		regionInfoOut << std::endl;
	}
}

void BamRegionInvestigator::writeBasicInfo(
		const std::vector<std::shared_ptr<RegionInfo>> & pairs,
		const std::string & sampName, const OutOptions & outOpts) const {
	OutputStream regionInfoOut(outOpts);
	regionInfoOut << "#chrom\tstart\tend\tname\tlength\tstrand\treadTotal\tfullySpanningReads\tperBaseCoverage";
	if(pars_.mapQualityCutOffForMultiMap_  > pars_.mapQualityCutOff_){
		regionInfoOut << "\tperBaseCoverageFromMapQOf<=" << pars_.mapQualityCutOffForMultiMap_;
	}
	regionInfoOut << "\t#OfPairedReads\tproperPairs\timproperMateUnmapped\tproperPairFrac\tproperPairAndMateUnmappedFrac";
	regionInfoOut << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : pairs){

		auto bedOut = p->region_.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		regionInfoOut << "\textraField"<<t;
	}
	regionInfoOut << "\n";
	for(const auto & p : pairs){
		auto bedOut = p->region_.genBedRecordCore();
		regionInfoOut << bedOut.toDelimStr();
		regionInfoOut << "\t" << p->totalReads_;
		regionInfoOut << "\t" << p->totalFullySpanningReads_;
		regionInfoOut << "\t" << p->coverage_/static_cast<double>(p->region_.getLen());
		if(pars_.mapQualityCutOffForMultiMap_  > pars_.mapQualityCutOff_){
			regionInfoOut << "\t" << p->multiMapCoverage_/static_cast<double>(p->region_.getLen());
		}
		regionInfoOut << "\t" << p->totalPairedReads_
				<< '\t' << p->totalProperPairedReads_
				<< '\t' << p->totalOneMateUnmappedImproperPairs_
				<< '\t' << p->getProperPairFrac()
				<< '\t' << p->getProperPairMateUnmappedFrac();
		regionInfoOut << "\t" << sampName;
		for(const auto & extra : bedOut.extraFields_){
			regionInfoOut << "\t" << extra;
		}
		uint32_t leftOverExtra = maxExtraFields - bedOut.extraFields_.size();
		for(uint32_t t = 0; t < leftOverExtra; ++t){
			regionInfoOut << "\t";
		}
		regionInfoOut << std::endl;
	}
}

void BamRegionInvestigator::writeBasicInfoWithHapRes(
		const std::vector<std::shared_ptr<RegionInfo>> & pairs,
		const std::string & sampName, const OutOptions & outOpts) const {
	OutputStream regionInfoOut(outOpts);
	regionInfoOut << "#chrom\tstart\tend\tname\tlength\tstrand\tsuccess\tuniqHaps\treadTotal\treadTotalUsed\tfullySpanningReads\tperBaseCoverage";
	if(pars_.mapQualityCutOffForMultiMap_  > pars_.mapQualityCutOff_){
		regionInfoOut << "\tperBaseCoverageFromMapQOf<=" << pars_.mapQualityCutOffForMultiMap_;
	}
	regionInfoOut << "\t#OfPairedReads\tproperPairs\timproperMateUnmapped\tproperPairFrac\tproperPairAndMateUnmappedFrac";
	regionInfoOut << "\tsample";
	uint32_t maxExtraFields = 0;
	for(const auto & p : pairs){

		auto bedOut = p->region_.genBedRecordCore();
		if(bedOut.extraFields_.size() > maxExtraFields){
			maxExtraFields = bedOut.extraFields_.size();
		}
	}
	for(uint32_t t = 0; t < maxExtraFields; ++t){
		regionInfoOut << "\textraField"<<t;
	}
	regionInfoOut << "\n";
	std::unordered_map<std::string, double> regionLengthSumByUid;
	for(const auto & p : pairs){
		regionLengthSumByUid[p->region_.uid_] += p->region_.getLen();
	}
	for(const auto & p : pairs){
		auto bedOut = p->region_.genBedRecordCore();
		regionInfoOut << bedOut.toDelimStr();
		regionInfoOut << "\t" << njh::boolToStr(p->infoCalled_);
		regionInfoOut << "\t" << p->uniqHaps_;
		regionInfoOut << "\t" << p->totalReads_;
		regionInfoOut << "\t" << p->totalFinalReads_;
		regionInfoOut << "\t" << p->totalFullySpanningReads_;
		regionInfoOut << "\t" << p->coverage_/static_cast<double>(regionLengthSumByUid[p->region_.uid_] );
		if(pars_.mapQualityCutOffForMultiMap_  > pars_.mapQualityCutOff_){
			regionInfoOut << "\t" << p->multiMapCoverage_/static_cast<double>(regionLengthSumByUid[p->region_.uid_] );
		}
		regionInfoOut << "\t" << p->totalPairedReads_
				<< '\t' << p->totalProperPairedReads_
				<< '\t' << p->totalOneMateUnmappedImproperPairs_
				<< '\t' << p->getProperPairFrac()
				<< '\t' << p->getProperPairMateUnmappedFrac();
		regionInfoOut << "\t" << sampName;
		for(const auto & extra : bedOut.extraFields_){
			regionInfoOut << "\t" << extra;
		}
		uint32_t leftOverExtra = maxExtraFields - bedOut.extraFields_.size();
		for(uint32_t t = 0; t < leftOverExtra; ++t){
			regionInfoOut << "\t";
		}
		regionInfoOut << std::endl;
	}
}


}  // namespace njhseq
