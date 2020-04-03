/*
 * CarmenRunner.cpp
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

#include "CarmenRunner.hpp"

namespace njhseq {

void checkForColsThrow(const table & tab, const VecStr & cols,
		const std::string & funcName) {
	if (!tab.containsColumns(cols)) {
		std::stringstream ss;
		ss << funcName << ", error, table " << " should have the following columns"
				<< "\n";
		ss << njh::conToStr(cols, ",") << "\n";
		ss << "only has: " << njh::conToStr(tab.columnNames_, ",") << "\n";
		throw std::runtime_error { ss.str() };
	}
}

std::vector<std::string> countryColorsHex = {"#c3d290",
		"#6e45cf",
		"#95da45",
		"#c750cb",
		"#5fd67a",
		"#592e82",
		"#d3c847",
		"#697dcd",
		"#d5903c",
		"#462a4e",
		"#4f863b",
		"#cc4183",
		"#6fcbb1",
		"#d24232",
		"#91bfd9",
		"#793036",
		"#ccae9b",
		"#35382b",
		"#c88bbd",
		"#6b6e3a",
		"#d87370",
		"#577184",
		"#925b32"};


Json::Value CarmenRunner::LocationInfo::exportJson() const {
	Json::Value ret;
	ret["country"] = njh::json::toJson(country_);
	ret["site"] = njh::json::toJson(site_);
	ret["size"] = njh::json::toJson(total_);
	ret["coords"] = njh::json::toJson(std::vector<double> { long_, lat_ });
	ret["region"] = njh::json::toJson(region_);
	return ret;
}

CarmenRunner::CarmenRunner(const Json::Value & config) :
		SeqApp(config) {
	serverResourceDir_ = njh::appendAsNeededRet(config["resources"].asString(),
			"/");

	masterDir_ = config["masterDir"].asString();
	sitesGpsFnp_  = config["sitesGpsFnp"].asString();


	jsFiles_->addFiles(
	 njh::files::gatherFiles(njh::files::make_path(serverResourceDir_, "carmen/js"), ".js"));

	cssFiles_->addFiles(
			njh::files::gatherFiles(
					njh::files::make_path(serverResourceDir_, "carmen/css"), ".css"));

	addScripts(njh::files::make_path(serverResourceDir_, "carmen"));

	table siteTab(sitesGpsFnp_, "\t", true);
	checkForColsThrow(siteTab,VecStr { "site", "country", "longitude",
		"latitude" }, __PRETTY_FUNCTION__ );



	for (const auto & row : siteTab.content_) {
		LocationInfo locInfo;
		locInfo.country_ = row[siteTab.getColPos("country")];
		locInfo.site_ = row[siteTab.getColPos("site")];
		if(siteTab.containsColumn("region")){
			locInfo.region_ = row[siteTab.getColPos("region")];
		}
		locInfo.lat_ = njh::lexical_cast<double>(row[siteTab.getColPos("latitude")]);
		locInfo.long_ = njh::lexical_cast<double>(row[siteTab.getColPos("longitude")]);
		locInfos_.emplace(locInfo.country_ + "-" + locInfo.site_, locInfo);
	}

	auto processInputNames = [](const std::string & input){
		return njh::replaceString(njh::strToLowerRet(input), " ", "_");
	};

	table siteTotalsTab(njh::files::make_path(masterDir_, "locationsTotals.tab.txt"), "\t", true);
	checkForColsThrow(siteTotalsTab, VecStr { "site", "country", "count"}, __PRETTY_FUNCTION__ );
	VecStr missing;
	for (const auto & row : siteTotalsTab.content_) {
		std::string country = row[siteTotalsTab.getColPos("country")];
		std::string site = row[siteTotalsTab.getColPos("site")];
		if("na" == processInputNames(country) &&
				"na" == processInputNames(site)){
			continue;
		}
		std::string uid = processInputNames(country) + "-" + processInputNames(site);
		std::cout << country << std::endl;
		std::cout << processInputNames(country) << std::endl << std::endl;

		if(!njh::in(uid, locInfos_)){
			missing.emplace_back(uid);
		}else{
			locInfos_.at(uid).total_ = njh::lexical_cast<double>(row[siteTotalsTab.getColPos("count")]);
		}
	}

	if(!missing.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ <<  ", error " << njh::conToStr(missing, ",") << " not found in gps locations input" << "\n";
		ss << "Optiions are: " << njh::conToStr(njh::getVecOfMapKeys(locInfos_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}

	table hapInCountriesCountTab(njh::files::make_path(masterDir_, "haplotypeInCountriesCount.tab.txt"), "\t", true);
	auto hapInCountriesCountTabCountrySplit = hapInCountriesCountTab.splitTableOnColumn("country");
	checkForColsThrow(hapInCountriesCountTab, VecStr { "country", "site", "haplotype", "alternativeName", "count", "frac"}, __PRETTY_FUNCTION__ );
	for (const auto & countryTab : hapInCountriesCountTabCountrySplit) {
		auto hapInCountriesCountTabSiteSplit = countryTab.second.splitTableOnColumn("site");
		for(const auto & siteTab : hapInCountriesCountTabSiteSplit){
			std::string country = countryTab.first;
			std::string site = siteTab.first;
			if("na" == processInputNames(country) &&
					"na" == processInputNames(site)){
				continue;
			}
			std::string uid = processInputNames(country) + "-" + processInputNames(site);
			if(!njh::in(uid, locInfos_)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ <<  ", error " << uid << " not found in gps locations input" << "\n";
				ss << "Optiions are: " << njh::conToStr(njh::getVecOfMapKeys(locInfos_), ", ") << "\n";
				throw std::runtime_error{ss.str()};
			}
			for(const auto & row : siteTab.second.content_){
				LocationInfo::HapInfo hInfo(row[siteTab.second.getColPos("haplotype")],
						VecStr{row[siteTab.second.getColPos("alternativeName")]},
						njh::lexical_cast<double>(row[siteTab.second.getColPos("count")]),
						njh::lexical_cast<double>(row[siteTab.second.getColPos("frac")])
						);
				locInfos_.at(uid).hapCounts_.emplace(hInfo.hapName_, hInfo);
			}
		}
	}


	seqs_->updateAddCache("combined",
						SeqIOOptions::genFastaIn(njh::files::make_path(masterDir_, "combined.fasta")));


}



std::vector<std::shared_ptr<restbed::Resource>> CarmenRunner::getAllResources() {
	auto ret = super::getAllResources();
	//main page
	ret.emplace_back(mainPage());
	//basic info
	ret.emplace_back(getWorldJson());
	ret.emplace_back(getSitesJson());

	ret.emplace_back(getCombinedSeqs());
	ret.emplace_back(getCombinedNwk());
	ret.emplace_back(getCombinedJson());

	ret.emplace_back(getAllNames());


	return ret;
}

VecStr CarmenRunner::requiredOptions() const {
	return concatVecs(super::requiredOptions(), VecStr { "numThreads",
			"resources", "sitesGpsFnp", "masterDir"});
}

void CarmenRunner::getAllNamesHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();

	//input names
	//ref names
	//loc uids
		//countries
		//sites
	Json::Value ret;

	std::set<std::string> inputNames;
	std::set<std::string> refNames;

	seqInfo seq;
	auto inputOpts = SeqIOOptions::genFastaIn(njh::files::make_path(masterDir_, "input.fasta"));
	SeqInput inputReader(inputOpts);
	inputReader.openIn();
	while(inputReader.readNextRead(seq)){
		inputNames.insert(seq.name_);
	}
	auto refOpts = SeqIOOptions::genFastaIn(njh::files::make_path(masterDir_, "refSeqs.fasta"));
	SeqInput refReader(refOpts);
	refReader.openIn();
	while(refReader.readNextRead(seq)){
		refNames.insert(seq.name_);
	}
	std::unordered_map<std::string, std::set<std::string>> countries;
	for (const auto & loc : locInfos_) {
		countries[loc.second.country_].insert(loc.second.site_);
	}
	ret["inputNames"] = njh::json::toJson(inputNames);
	ret["refNames"] = njh::json::toJson(refNames);
	//ret["countries"] = njh::json::toJson(countries);
	auto & countriesJson = ret["countries"];
	for(const auto & country :countries){
		Json::Value countryJson;
		countryJson["country"] = njh::json::toJson(country.first);
		countryJson["sites"] = njh::json::toJson(country.second);
		countriesJson.append(countryJson);
	}
	auto countryNames = getVectorOfMapKeys(countries);
	njh::sort(countryNames);
	auto refColors = getColorsForNames(VecStr{refNames.begin(), refNames.end()});
	auto countryColors = getColorsForNames(countryNames);
	std::unordered_map<std::string, std::string> handSelectedCountryColors;
	uint32_t countryCount = 0;
	for(const auto & country : countryNames){
		if(countryCount >= countryColorsHex.size()){
			handSelectedCountryColors[country] = njh::color().getHexStr();
		}else{
			handSelectedCountryColors[country] = countryColorsHex[countryCount];
		}
		++countryCount;
	}
	ret["refColors"] = njh::json::toJson(refColors);
	ret["countryColors"] = njh::json::toJson(handSelectedCountryColors);

	std::string body = njh::json::writeAsOneLine(ret);
	auto header = HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, header);
}



void CarmenRunner::getSitesJsonHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	njh::randomGenerator gen;
	Json::Value ret;
	std::vector<double> sizes;
	for (const auto & loc : locInfos_) {
		sizes.emplace_back(loc.second.total_);
	}
	std::set<std::string> hapNames;
	for (const auto & loc : locInfos_) {
		for(const auto & hap : loc.second.hapCounts_){
			hapNames.insert(hap.first);
		}
	}
	auto hapColors = getColorsForNames(VecStr{hapNames.rbegin(), hapNames.rend()});
	scale<double> sizeScale(std::make_pair(*std::min_element(sizes.begin(), sizes.end()),
			*std::max_element(sizes.begin(), sizes.end()) ),
			std::make_pair(1000,20000));
	for (const auto & loc : locInfos_) {
		if(loc.second.total_ > 0){
			Json::Value siteJson;
			siteJson = loc.second.exportJson();
			siteJson["scaledSize"] = sizeScale.get(siteJson["size"].asDouble());
			Json::Value & proportions = siteJson["proportions"];
			uint32_t hapCount = 0;
			for(const auto & hap : loc.second.hapCounts_){
				proportions[hapCount]["name"] = hap.second.hapName_;
				proportions[hapCount]["altName"] = hap.second.altNames_.front();
				proportions[hapCount]["proportion"] = hap.second.count_;
				proportions[hapCount]["frac"] = hap.second.frac_;
				proportions[hapCount]["color"] = hapColors[hap.first].getHexStr();
				++hapCount;
			}
			ret.append(siteJson);
		}
	}
	std::string body = njh::json::writeAsOneLine(ret);
	auto header = HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, header);
}

void CarmenRunner::getWorldJsonHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto worldJsonFnp = njh::files::make_path(config_["seqServerCore"].asString(),
			"json/d3map/world-50m.json");
	std::string body = njh::files::get_file_contents(worldJsonFnp, messFac_->debug_);
	auto header = HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, header);
}


void CarmenRunner::getCombinedNwkHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto combinedNwk = njh::files::make_path(masterDir_, "combined.nwk");
	std::string body = "";
	if(bfs::exists(combinedNwk)){
		body = njh::files::get_file_contents(combinedNwk, false);
	}
	auto header = HeaderFactory::initiatePlainTxtHeader(body);
	session->close(restbed::OK, body, header);
}

//page handlers
void CarmenRunner::mainPageHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto body = genHtmlDoc(rootName_, pages_.at("mainPage.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}

void CarmenRunner::getCombinedJsonHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto combinedJsonFnp = njh::files::make_path(masterDir_, "combined.json");
	auto hapCountsFnp = njh::files::make_path(masterDir_, "haplotypeCounts.tab.txt");

	Json::Value inputJson = njh::json::parseFile(combinedJsonFnp.string());
	std::set<std::string> hapNames;
	for (const auto & loc : locInfos_) {
		for(const auto & hap : loc.second.hapCounts_){
			hapNames.insert(hap.first);
		}
	}

	table hapCountsTab(hapCountsFnp, "\t", true);
	auto splitOnName = hapCountsTab.splitTableOnColumn("alternativeName");

	std::vector<double> totals;
	for(const auto & name : splitOnName){
		double nameTotal = 0;
		for(const auto & row : name.second.content_){
			auto count = estd::stou(row[name.second.getColPos("count")]);
			nameTotal+= count;
		}
		totals.emplace_back(std::sqrt(nameTotal/3.14));
	}

	std::unordered_map<std::string, std::set<std::string>> countries;
	for (const auto & loc : locInfos_) {
		countries[loc.second.country_].insert(loc.second.site_);
	}

	auto countryNames = getVectorOfMapKeys(countries);
	njh::sort(countryNames);
	auto countryColors = getColorsForNames(countryNames);
	std::unordered_map<std::string, std::string> handSelectedCountryColors;
	uint32_t countryCount = 0;

	for(const auto & country : countryNames){
		if(countryCount >= countryColorsHex.size()){
			handSelectedCountryColors[country] = njh::color().getHexStr();
		}else{
			handSelectedCountryColors[country] = countryColorsHex[countryCount];
		}
		++countryCount;
	}

	scale<double> sizeScale(std::make_pair(*std::min_element(totals.begin(), totals.end()),
			*std::max_element(totals.begin(), totals.end())),
			std::make_pair(10.0, 30.0));

	for(const auto & name : splitOnName){

		std::unordered_map<std::string, uint32_t> countryCounts;
		double nameTotal = 0;
		for(const auto & row : name.second.content_){
			auto count = estd::stou(row[name.second.getColPos("count")]);
			countryCounts[row[name.second.getColPos("country")]] += count;
			nameTotal+= count;
		}
		auto & pieces = inputJson[name.first]["pieces"];
		for(const auto & country : countryCounts){
			Json::Value countryVal;
			//countryVal["color"] = njh::json::toJson(countryColors[country.first].getHexStr());
			countryVal["color"] = njh::json::toJson(handSelectedCountryColors[country.first]);

			countryVal["country"] = njh::json::toJson(country.first);
			countryVal["proportion"] = njh::json::toJson(country.second);
			countryVal["frac"] = njh::json::toJson(country.second/nameTotal);
			pieces.append(countryVal);
		}
		inputJson[name.first]["size"] = njh::json::toJson(sizeScale.get(std::sqrt(nameTotal/3.14)));

		//inputJson[name.first]["size"] = njh::json::toJson(sizeScale.get(nameTotal));

	}


	std::set<std::string> refNames;
	auto refOpts = SeqIOOptions::genFastaIn(
			njh::files::make_path(masterDir_, "refSeqs.fasta"));
	SeqInput refReader(refOpts);
	refReader.openIn();
	seqInfo seq;
	while (refReader.readNextRead(seq)) {
		refNames.insert(seq.name_);
	}

	auto refColors = getColorsForNames(
			VecStr { refNames.begin(), refNames.end() });


	auto hapColors = getColorsForNames(VecStr{hapNames.rbegin(), hapNames.rend()});
	for(auto & val : inputJson){
		if(!njh::in(val["name"].asString(), splitOnName)){
			val["size"] = njh::json::toJson(10.0);
			auto & pieces = val["pieces"];
			Json::Value countryVal;
			countryVal["country"] = njh::json::toJson("none");
			countryVal["proportion"] = njh::json::toJson(1);
			countryVal["frac"] = njh::json::toJson(1);
			countryVal["color"] = njh::json::toJson(njh::color().getHexStr());
			pieces.append(countryVal);
		}
		if(val["refSeq"].asBool()){
			val["refColor"] = njh::json::toJson(refColors[val["refName"].asString()].getHexStr());
		}else{
			val["refColor"] = njh::json::toJson(njh::color().getHexStr());
		}
		val["color"] = njh::json::toJson(hapColors[val["originalName"].asString()].getHexStr());
	}

	std::string body = njh::json::writeAsOneLine(inputJson);
	auto header = HeaderFactory::initiateAppJsonHeader(body);

	session->close(restbed::OK, body, header);
}

std::shared_ptr<restbed::Resource> CarmenRunner::getCombinedJson(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getCombinedJson" } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
		getCombinedJsonHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> CarmenRunner::getCombinedSeqs(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getCombinedSeqs" } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getCombinedSeqsHandler(session);
					}));
	return resource;
}



std::shared_ptr<restbed::Resource> CarmenRunner::getCombinedNwk() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getCombinedNwk" } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
		getCombinedNwkHandler(session);
					}));
	return resource;
}



void CarmenRunner::getCombinedSeqsPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	Json::Value postData;
	if(!body.empty()){
		postData = njh::json::parse(std::string(body.begin(), body.end()));
	}
	Json::Value ret;

	std::string searchTerm = "combined";
	if (njh::in(searchTerm, seqs_->cache_)) {
		uint32_t sesUid = std::numeric_limits<uint32_t>::max();
		//check to see if there is a session already started associated with this seq
		if (!postData.isMember("sessionUID")) {
			sesUid = startSeqCacheSession();
		} else {
			sesUid = postData["sessionUID"].asUInt();
		}
		ret = seqsBySession_[sesUid]->getJson(searchTerm);
		ret["sessionUID"] = njh::json::toJson(sesUid);
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": couldn't find " << searchTerm
				<< " Couldn't find :" << searchTerm << ", options are: " << njh::conToStr(njh::getVecOfMapKeys(seqs_->cache_), ", ")
				<< std::endl;
	}

	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);

}
void CarmenRunner::getCombinedSeqsHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes &)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
		getCombinedSeqsPostHandler(ses, body);
					}));
}

//getting general info
std::shared_ptr<restbed::Resource> CarmenRunner::getWorldJson() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getWorldJson" } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getWorldJsonHandler(session);
					}));
	return resource;
}



std::shared_ptr<restbed::Resource> CarmenRunner::getAllNames(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getAllNames" } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getAllNamesHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> CarmenRunner::getSitesJson() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getSitesJson" } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getSitesJsonHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> CarmenRunner::mainPage() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						mainPageHandler(session);
					}));
	return resource;
}

void CarmenRunner::redirect(std::shared_ptr<restbed::Session> session,
		std::string errorMessage) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	std::cerr << errorMessage << std::endl;
	auto body = genHtmlDoc(rootName_, pages_.at("redirectPage.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}

void CarmenRunner::CarmenRunnerErrorHandler(const int statusCode,
		const std::exception& exception,
		const std::shared_ptr<restbed::Session>& session) {
	std::cerr << "statusCode: " << statusCode << std::endl;
	std::cerr << exception.what() << std::endl;
	if (session->is_open()) {
		session->close(statusCode, exception.what(), { { "Server", "Restbed" } });
	}
}

} /* namespace njhseq */
