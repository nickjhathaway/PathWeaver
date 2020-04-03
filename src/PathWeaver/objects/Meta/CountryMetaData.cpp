/*
 * CountryMetaData.cpp
 *
 *  Created on: Mar 19, 2017
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

#include "CountryMetaData.hpp"

namespace njhseq {



CountryMetaData::CountryMetaData(const bfs::path & metaFnp,
		const VecStr & samples) :
		meta_(metaFnp, std::set<std::string> { samples.begin(), samples.end() }) {

	meta_.checkForFieldsThrow(VecStr{"country","site"});
	//process fields;
	auto transformer = [](const std::string & subField){

		std::string ret;
		ret = njh::replaceString(njh::strToLowerRet(subField), " ","_");
		ret = njh::replaceString(ret, ".", "-");
		return ret;
	};
	meta_.transformSubFields("country", transformer);
	meta_.transformSubFields("site", transformer);

}

std::unordered_map<std::string, CountryMetaData::CountrySite> CountryMetaData::getLocsBySamples() const {
	std::unordered_map<std::string, CountrySite> ret;
	for (const auto & samp : meta_.samples_) {
		auto loc = CountrySite { meta_.groupData_.at("country")->getGroupForSample(
				samp), meta_.groupData_.at("site")->getGroupForSample(samp) };
		ret[samp] = loc;
	}
	return ret;
}

} /* namespace njhseq */
