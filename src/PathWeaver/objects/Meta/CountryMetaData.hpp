#pragma once
/*
 * CountryMetaData.hpp
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

#include <njhseq/utils.h>
#include <njhseq/objects/Meta/MultipleGroupMetaData.hpp>

namespace njhseq {

/**@brief class to hold meta data about samples with country and site locations
 *
 */
class CountryMetaData{
public:
	/**@brief a struct to hold country and site info
	 *
	 */
	struct CountrySite{
		std::string country_;/**< country name*/
		std::string site_;/**< sit name, normally a city*/
	};

	/**@brief Construct with meta data filename and the samples in the dataset
	 *
	 * @param metaFnp a path to a meta data file, should have a column named sample and each additional is another meta field, required fields, country and site
	 * @param samples the samples in the dataset
	 */
	CountryMetaData(const bfs::path & metaFnp, const VecStr & samples);

	MultipleGroupMetaData meta_;/**< meta data*/
	std::unordered_map<std::string, CountrySite> getLocsBySamples() const; /**< get the locations by sample*/



};



} /* namespace njhseq */
