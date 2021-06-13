#pragma once

/*
 * RegionInvestigatorInBam.hpp
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


#include <TwoBit.h>

#include <njhseq/BamToolsUtils.h>
#include <njhseq/objects/BioDataObject.h>
#include <njhseq/programUtils/seqSetUp.hpp>
#include <SeekDeep/objects/IlluminaUtils/PairedReadProcessor.hpp>

namespace njhseq {


class BamRegionInvestigator {
public:

	struct BamCountSpecficRegionsPars{
		BamCountSpecficRegionsPars();
		uint32_t extendAmount = 30;
		uint32_t numThreads = 1;
		uint32_t mappingQuality = 20;
		uint32_t baseQuality = 25;
		double matchIDCutOff = 0.70;

		uint32_t totalCountCutOff = 5;
		uint32_t perBaseCountCutOff = 3;
		bool forcePlusStrand = false;
		bool countDuplicates = false;
		uint32_t insertLengthCutOff = 1000;
		PairedReadProcessor::ProcessParams pairPars;

		void setDefaults(seqSetUp & setUp);
	};




	struct RegionInfo {

		RegionInfo(const GenomicRegion & region);

		GenomicRegion region_; //!< region
		uint32_t coverage_ = 0; //!< number of bases that fall in this region
		uint32_t multiMapCoverage_ = 0; //!< number of bases that fall in this region with a mapping quality less than multimapping qual cut off

		uint32_t uniqHaps_{0}; //!< number of unique haplotypes called by PathWeaver for this region
		bool infoCalled_{false}; //!< whether PathWeaver could be called for this region
		uint32_t totalReads_{0}; //!< the total number of reads that fall in this region
		uint32_t totalFinalReads_{0}; //!< the total number of reads that are used in the final analysis
		uint32_t totalFullySpanningReads_{0}; //!< the total number of reads that fully span this region

		uint32_t totalPairedReads_{0}; //!< total number of paired reads that fall in this region, will count both mates of a pair twice
		uint32_t totalProperPairedReads_{0};//!< total number of proper pairs within region
		uint32_t totalOneMateUnmappedImproperPairs_{0}; //!< total number of paired reads that are improper due to one mate being unmapped

		/**@brief calculate the estimated per base read coverage for this region
		 *
		 * @return the per base coverage for this region
		 */
		double getPerBaseCov() const;
		/**@brief the fraction of the total pairs that are proper pairs
		 *
		 * @return fraction of proper pairs range [0-1]
		 */
		double getProperPairFrac() const;
		/**@brief the fraction of the total pairs that are proper pairs or one mate was unmapped
		 *
		 * @return fraction of proper pairs range [0-1]
		 */
		double getProperPairMateUnmappedFrac() const;
	};

	struct BamRegionInvestigatorPars {
		bool countDups_ { false };
		uint32_t numThreads_ { 1 };
		uint32_t mapQualityCutOff_ { 0 }; //!< Don't include alignments below this map quality, non-inclusive
		uint32_t mapQualityCutOffForMultiMap_ { 5 };//!< Include alignments below this map quality in the multimapping count colum
	};

	BamRegionInvestigator(const BamRegionInvestigatorPars & pars);
	BamRegionInvestigatorPars pars_;

	std::vector<std::shared_ptr<RegionInfo>> getCoverageOnFromBam(
			const bfs::path & bamFnp,
			const std::vector<GenomicRegion> & regions) const;

	RegionInfo getCoverageForRegion(
			BamTools::BamReader & bReader, const GenomicRegion & region) const;

	RegionInfo getCoverageAndFullySpanningForRegion(
			BamTools::BamReader &bReader, const GenomicRegion &region,
			const seqInfo &refSeq, const BamCountSpecficRegionsPars &spanningReadsPar,
			aligner &alignerObj) const;


	struct ReadCountsRes{
		std::map<std::string, uint32_t> spanningReadCounts_;
		std::map<std::string, uint32_t> totalReadCounts_;
	};

	struct ReadCountsResSingle {
		uint32_t spanningReadCount_{0};
		uint32_t totalReadCount_{0};
	};

	ReadCountsRes getTotalAndFullySpanningReadCounts(const bfs::path & bamFnp,
			const std::vector<GenomicRegion> & regions,
			const BamCountSpecficRegionsPars & spanningReadsPar,
			const bfs::path & genomeFnp) const;


	void getTotalAndFullySpanningReadCounts(const bfs::path & bamFnp,
			const std::vector<GenomicRegion> & regions,
			const BamCountSpecficRegionsPars & spanningReadsPar,
			const bfs::path & genomeFnp,
			std::vector<std::shared_ptr<RegionInfo>> & pairs) const;

	void getTotalAndFullySpanningReadCounts(const bfs::path & bamFnp,
			const std::vector<GenomicRegion> & regions,
			const BamCountSpecficRegionsPars & spanningReadsPar,
			const bfs::path & genomeFnp,
			std::unordered_map<std::string, std::shared_ptr<RegionInfo>> & pairsByUid) const;

	ReadCountsResSingle getFullSpanningCountForRegion(const GenomicRegion & currentRegion, const seqInfo & refSeq,
			const BamCountSpecficRegionsPars & spanningReadsPar,
			BamTools::BamReader &currentBReader, aligner & alignerObj) const;

	void updateFullSpanningCountForRegion(RegionInfo & currentRegion, const seqInfo & refSeq,
			const BamCountSpecficRegionsPars & spanningReadsPar,
			BamTools::BamReader &currentBReader, aligner & alignerObj) const;


	std::vector<std::shared_ptr<RegionInfo>> getCoverageAndFullSpanningReads(
			const bfs::path & bamFnp,
			const std::vector<GenomicRegion> & regions,
			const BamCountSpecficRegionsPars & spanningReadsPar,
			const bfs::path & genomeFnp) const;



	void writeCovInfo(const std::vector<std::shared_ptr<RegionInfo>> & pairs,
			const std::string & sampName,
			const OutOptions & outOpts) const;

	void writeBasicInfo(const std::vector<std::shared_ptr<RegionInfo>> & pairs,
			const std::string & sampName,
			const OutOptions & outOpts) const;

	void writeBasicInfoWithHapRes(const std::vector<std::shared_ptr<RegionInfo>> & pairs,
			const std::string & sampName,
			const OutOptions & outOpts) const;

};






}  // namespace njhseq
