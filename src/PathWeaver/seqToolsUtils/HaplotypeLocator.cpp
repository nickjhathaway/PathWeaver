/*
 * HaplotypeLocator.cpp
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

#include "HaplotypeLocator.hpp"

namespace njhseq {

HaplotypeLocator::LocationsCounter::LocationsCounter(const std::string & uid) :
		uid_(uid) {
}

void HaplotypeLocator::LocationsCounter::increaseCount(
		const std::string & country, const std::string & site, uint32_t count) {
	locationsTotals_[country][site] += count;
}

void HaplotypeLocator::LocationsCounter::increaseCount(
		const std::string & country, uint32_t count) {
	increaseCount(country, "na", count);
}

uint32_t HaplotypeLocator::LocationsCounter::getTotal() const {
	uint32_t ret = 0;
	for (const auto & country : locationsTotals_) {
		for (const auto & site : country.second) {
			ret += site.second;
		}
	}
	return ret;
}

bool HaplotypeLocator::LocationsCounter::hasLoc(const std::string & country,
		const std::string & site) const {
	auto search = locationsTotals_.find(country);
	if (search != locationsTotals_.end()) {
		auto siteSearch = search->second.find(site);
		if (siteSearch != search->second.end()) {
			return true;
		}
	}
	return false;
}
uint32_t HaplotypeLocator::LocationsCounter::getCount(
		const std::string & country, const std::string & site) const {
	if (hasLoc(country, site)) {
		return locationsTotals_.at(country).at(site);
	}
	return 0;
}

HaplotypeLocator::SeqCollection::SeqCollection(const std::string & uid) :
		uid_(uid) {
}

HaplotypeLocator::HaplotypeLocator(const std::string & name,
		const std::shared_ptr<std::vector<seqInfo>> & inputSeqs_) :
		name_(name), inputSeqs_(inputSeqs_), totalLocCounter_("all") {

}

void HaplotypeLocator::addRefSeqs(const std::vector<seqInfo> & refSeqs) {
	refSeqs_ = std::make_shared<std::vector<seqInfo>>(refSeqs);
}

void HaplotypeLocator::collapseInputSeqs() {
	if (nullptr == inputSeqs_) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, inputSeqs_ isn't set" << "\n";
		throw std::runtime_error { ss.str() };
	}
	//clear any info
	collapsedSeqs_ = std::make_shared<std::vector<seqInfo>>();
	alternateUids_.clear();
	hapLocCounters_.clear();
	collapsingSeqs_.clear();
	for (const auto seq : *inputSeqs_) {
		auto search = collapsingSeqs_.find(seq.seq_);
		if (collapsingSeqs_.end() != search) {
			collapsingSeqs_.at(seq.seq_).seqs_.emplace_back(seq);
		} else {
			SeqCollection col(estd::to_string(collapsingSeqs_.size()));
			col.seqs_.emplace_back(seq);
			collapsingSeqs_.emplace(seq.seq_, col);
		}
	}
	auto keys = getVectorOfMapKeys(collapsingSeqs_);
	njh::sort(keys,
			[this](const std::string & key1, const std::string& key2) {
				return collapsingSeqs_.at(key1).seqs_.size() > collapsingSeqs_.at(key2).seqs_.size();
			});
	uint32_t count = 0;
	uint32_t total = 0;
	std::unordered_map<std::string, std::string> refSeqToName;
	if(nullptr != refSeqs_){
		for (const auto & seq : *refSeqs_) {
			refSeqToName[seq.seq_] = seq.name_;
		}
	}

	for (const auto & key : keys) {
		std::string uid = "";
		if (njh::in(key, refSeqToName)) {
			uid = refSeqToName[key];
		} else {
			uid = name_ + "."
					+ leftPadNumStr<uint32_t>(count, collapsingSeqs_.size());
			++count;
		}
		collapsedSeqs_->emplace_back(uid, key);
		collapsedSeqs_->back().cnt_ = collapsingSeqs_.at(key).seqs_.size();
		collapsingSeqs_.at(key).uid_ = uid;
		total += collapsingSeqs_.at(key).seqs_.size();
		LocationsCounter counter(uid);
		for (const auto & seq : collapsingSeqs_.at(key).seqs_) {
			MetaDataInName meta(seq.name_);
			totalLocCounter_.increaseCount(meta.getMeta("country"),
					meta.getMeta("site"), 1);
			counter.increaseCount(meta.getMeta("country"), meta.getMeta("site"), 1);
		}
		hapLocCounters_.emplace(uid, counter);
	}
	njh::for_each(*collapsedSeqs_,
			[&total](seqInfo & seq) {seq.setFractionByCount(total);});
	//add in any ref seqs that weren't found
	if (nullptr != refSeqs_) {
		for (const auto & ref : *refSeqs_) {
			if (!njh::in(ref.name_, hapLocCounters_)) {
				collapsedSeqs_->emplace_back(ref.name_, ref.seq_);
			}
		}
	}
	//by default make alternative uid with own name
	for (const auto & seq : *collapsedSeqs_) {
		auto search = alternateUids_.find(seq.name_);
		if (alternateUids_.end() == search) {
			alternateUids_[seq.name_] = seq.name_;
		}
	}
}

void HaplotypeLocator::updateAlternativeUids(
		const std::vector<seqInfo> & alts) {
	for (const auto & seq : alts) {
		auto search = collapsingSeqs_.find(seq.seq_);
		if (collapsingSeqs_.end() != search) {
			alternateUids_[search->second.uid_] = seq.name_;
		}
	}
}

table HaplotypeLocator::getCountryTotals() const {
	table ret(VecStr { "country", "site", "count" });
	auto countryKeys = getVectorOfMapKeys(totalLocCounter_.locationsTotals_);
	njh::sort(countryKeys);
	for (const auto & country : countryKeys) {
		auto siteKeys = getVectorOfMapKeys(
				totalLocCounter_.locationsTotals_.at(country));
		njh::sort(siteKeys);
		for (const auto & site : siteKeys) {
			ret.addRow(country, site,
					totalLocCounter_.locationsTotals_.at(country).at(site));
		}
	}
	return ret;
}

table HaplotypeLocator::getHaplotypeCountsInCountries() const {
	table ret(VecStr { "country", "site", "haplotype", "alternativeName", "count",
			"frac" });
	auto countryKeys = getVectorOfMapKeys(totalLocCounter_.locationsTotals_);
	njh::sort(countryKeys);
	for (const auto & country : countryKeys) {
		auto siteKeys = getVectorOfMapKeys(
				totalLocCounter_.locationsTotals_.at(country));
		njh::sort(siteKeys);
		for (const auto & site : siteKeys) {
			std::unordered_map<std::string, uint32_t> hCounts;
			double total = 0;
			auto hapKeys = getVectorOfMapKeys(hapLocCounters_);
			njh::sort(hapKeys);
			for (const auto & hapCountKey : hapKeys) {
				if (hapLocCounters_.at(hapCountKey).hasLoc(country, site)) {
					hCounts[hapCountKey] = hapLocCounters_.at(hapCountKey).getCount(
							country, site);
					total += hapLocCounters_.at(hapCountKey).getCount(country, site);
				}
			}
			for (const auto & hCount : hCounts) {
				ret.addRow(country, site, hCount.first, alternateUids_.at(hCount.first),
						hCount.second, hCount.second / total);
			}
		}
	}
	return ret;
}

table HaplotypeLocator::getHaplotypeCounts() const {
	table ret(VecStr { "haplotype", "alternativeName", "country", "site", "count",
			"frac" });
	auto hapKeys = getVectorOfMapKeys(hapLocCounters_);
	njh::sort(hapKeys);
	for (const auto & hapCountKey : hapKeys) {
		double total = hapLocCounters_.at(hapCountKey).getTotal();
		auto countryKeys = getVectorOfMapKeys(
				hapLocCounters_.at(hapCountKey).locationsTotals_);
		njh::sort(countryKeys);
		for (const auto & country : countryKeys) {
			auto siteKeys = getVectorOfMapKeys(
					hapLocCounters_.at(hapCountKey).locationsTotals_.at(country));
			njh::sort(siteKeys);
			for (const auto & site : siteKeys) {
				ret.addRow(hapCountKey, alternateUids_.at(hapCountKey), country, site,
						hapLocCounters_.at(hapCountKey).locationsTotals_.at(country).at(
								site),
						hapLocCounters_.at(hapCountKey).locationsTotals_.at(country).at(
								site) / total);
			}
		}
	}
	return ret;
}

}  // namespace njhseq
