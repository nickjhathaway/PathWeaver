#pragma once
/*
 * HaplotypeLocator.hpp
 *
 *  Created on: Mar 31, 2017
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

#include <njhseq/utils.h>
#include <njhseq/objects/seqObjects/BaseObjects/seqInfo.hpp>
#include <njhseq/objects/dataContainers/tables/table.hpp>


namespace njhseq {
class HaplotypeLocator {
public:

	class LocationsCounter {
	public:
		LocationsCounter(const std::string & uid);
		std::string uid_;
		std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> locationsTotals_;

		void increaseCount(const std::string & country, const std::string & site,
				uint32_t count);

		void increaseCount(const std::string & country, uint32_t count);

		uint32_t getTotal() const;

		bool hasLoc(const std::string & country, const std::string & site) const;
		uint32_t getCount(const std::string & country,
				const std::string & site) const;
	};

	class SeqCollection {
	public:
		SeqCollection(const std::string & uid);
		std::string uid_;
		std::vector<seqInfo> seqs_;
	};

	HaplotypeLocator(const std::string & name,
			const std::shared_ptr<std::vector<seqInfo>> & inputSeqs_);
	std::string name_;
	std::shared_ptr<std::vector<seqInfo>> inputSeqs_;
	std::shared_ptr<std::vector<seqInfo>> refSeqs_;
	LocationsCounter totalLocCounter_;

	std::shared_ptr<std::vector<seqInfo>> collapsedSeqs_;
	std::unordered_map<std::string, std::string> alternateUids_;
	std::unordered_map<std::string, LocationsCounter> hapLocCounters_;
	std::unordered_map<std::string, SeqCollection> collapsingSeqs_;

	void addRefSeqs(const std::vector<seqInfo> & refSeqs);
	void collapseInputSeqs();

	void updateAlternativeUids(const std::vector<seqInfo> & alts);

	table getCountryTotals() const;
	table getHaplotypeCountsInCountries() const;
	table getHaplotypeCounts() const;

};

}  // namespace njhseq
