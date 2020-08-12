/*
 * CarmenServerRunner.cpp
 *
 *  Created on: May 18, 2015
 *      Author: nickhathaway
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

#include "CarmenServer.hpp"

#include <seqServer.h>
#include "PathWeaver/Carmen.h"

namespace njhseq {
CarmenServerRunner::CarmenServerRunner()
    : njh::progutils::ProgramRunner(
          {
					 addFunc("runCarmen", runCarmen, false),
           },
          "CarmenServer") {}

int CarmenServerRunner::runCarmen(const njh::progutils::CmdArgs & inputCommands){
	bfs::path sitesFnp = "";
	bfs::path masterDir = "";
	uint32_t numThreads = 1;
	SeqAppCorePars seqServerCorePars;
	seqServerCorePars.name_ = "carmen0";
	seqServerCorePars.port_ = 10000;
#if defined(PATHWEAVERDEBUG)
	bfs::path resourceDirName = njh::files::make_path(PathWeaverDebug_INSTALLDIR,
			"etc/serverResources");
#else
	bfs::path resourceDirName = njh::files::make_path(PathWeaver_INSTALLDIR,
			"etc/serverResources");
#endif

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(sitesFnp, "--sitesFnp", "Table with columns, site,country,longitude,latitude", true);
	setUp.setOption(masterDir, "--masterDir", "Master Directory containing input information", true);
	setUp.setOption(numThreads, "--numThreads", "Number of Threads");
	setUp.setOption(resourceDirName, "--resourceDirName",
			"Name of the resource Directory where the js and html is located",
			!bfs::exists(resourceDirName));
	resourceDirName = njh::appendAsNeededRet(resourceDirName.string(), "/");
	seqServerCorePars.setCoreOptions(setUp);
	setUp.finishSetUp(std::cout);

  //check for html/js/css resource files
	//checkExistenceThrow(resourceDirName);
  //checkExistenceThrow(njh::files::make_path(resourceDirName,"map/"));


  Json::Value appConfig;
  seqServerCorePars.addCoreOpts(appConfig);
  appConfig["resources"] = njh::json::toJson(resourceDirName);
  appConfig["sitesGpsFnp"] = njh::json::toJson(sitesFnp);
  appConfig["masterDir"] = njh::json::toJson(masterDir);

	if (setUp.pars_.verbose_) {
		std::cout << seqServerCorePars.getAddress() << std::endl;
	}

  CarmenRunner server(appConfig);
	auto resources = server.getAllResources();

	auto settings = std::make_shared<restbed::Settings>();
	settings->set_port(seqServerCorePars.port_);
	settings->set_bind_address(seqServerCorePars.bindAddress_);
	settings->set_default_header("Connection", "close");
	settings->set_worker_limit(4);
	restbed::Service service;
	service.set_error_handler(CarmenRunner::CarmenRunnerErrorHandler);
	for(const auto & resource : resources){
		service.publish(resource);
	}
	try {
		service.start(settings);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
	}


	return 0;
}


} /* namespace njhseq */
