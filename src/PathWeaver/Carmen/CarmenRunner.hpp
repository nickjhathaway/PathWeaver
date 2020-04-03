#pragma once
/*
 * CarmenRunner.hpp
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

#include <seqServer.h>
#include <njhseq/utils.h>

namespace njhseq {

class CarmenRunner: public njhseq::SeqApp {
	typedef njhseq::SeqApp super;

	bfs::path serverResourceDir_;
	bfs::path masterDir_;
	bfs::path sitesGpsFnp_;

public:
	class LocationInfo {
	public:
		std::string country_;
		std::string site_;
		std::string region_;
		double lat_;
		double long_;

		double total_ = 0;

		class HapInfo {
		public:
			HapInfo(const std::string & hapName, const VecStr & altNames,
					uint32_t count, double frac) :
					hapName_(hapName), altNames_(altNames), count_(count), frac_(frac) {
			}
			std::string hapName_;
			VecStr altNames_;
			uint32_t count_;
			double frac_;
		};
		std::unordered_map<std::string, HapInfo> hapCounts_;

		Json::Value exportJson() const;

	};
	std::unordered_map<std::string, LocationInfo> locInfos_;

	CarmenRunner(const Json::Value & config);


	virtual std::vector<std::shared_ptr<restbed::Resource>> getAllResources();

	virtual VecStr requiredOptions() const;


	void getSitesJsonHandler(std::shared_ptr<restbed::Session> session);

	void getWorldJsonHandler(std::shared_ptr<restbed::Session> session);
	void getAllNamesHandler(std::shared_ptr<restbed::Session> session);

	void getCombinedNwkHandler(std::shared_ptr<restbed::Session> session);
	void getCombinedJsonHandler(std::shared_ptr<restbed::Session> session);

	void getCombinedSeqsPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getCombinedSeqsHandler(std::shared_ptr<restbed::Session> session);

	//page handlers
	void mainPageHandler(std::shared_ptr<restbed::Session> session);

	//getting general info
	std::shared_ptr<restbed::Resource> getWorldJson();

	std::shared_ptr<restbed::Resource> getAllNames();

	std::shared_ptr<restbed::Resource> getSitesJson();
	std::shared_ptr<restbed::Resource> getCombinedJson();
	std::shared_ptr<restbed::Resource> getCombinedSeqs();
	std::shared_ptr<restbed::Resource> getCombinedNwk();

	std::shared_ptr<restbed::Resource> mainPage();



	void redirect(std::shared_ptr<restbed::Session> session,
			std::string errorMessage);


	static void CarmenRunnerErrorHandler(const int statusCode,
			const std::exception& exception,
			const std::shared_ptr<restbed::Session>& session);

};

} /* namespace njhseq */
