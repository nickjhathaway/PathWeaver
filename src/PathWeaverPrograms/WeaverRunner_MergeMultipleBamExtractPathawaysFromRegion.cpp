/*
 * WeaverRunner_MergeMultipleBamExtractPathawaysFromRegion.cpp
 *
 *  Created on: Aug 29, 2021
 *      Author: nick
 *
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

#include "WeaverRunner.hpp"

#include <njhseq/BamToolsUtils.h>
#include <njhseq/objects/BioDataObject.h>
#include <njhseq/GenomeUtils.h>

#include <njhcpp/concurrency/LockableJsonLog.hpp>
#include <njhseq/concurrency/pools/BamReaderPool.hpp>

#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/seqContainers.h>
#include <njhseq/objects/seqObjects/Clusters/cluster.hpp>

#include "PathWeaver/objects/bam/RegionInvestigatorInBam.hpp"
#include "PathWeaver/objects/Meta/CountryMetaData.hpp"
#include "PathWeaver/seqToolsUtils/HaplotypeLocator.hpp"
#include "PathWeaver/PathFinding.h"



namespace njhseq {

int WeaverRunner::MergeMultipleBamExtractPathwaysFromRegion(
		const njh::progutils::CmdArgs & inputCommands) {

	std::vector<bfs::path> weavedResultsDirs;
	bfs::path mergingBed;

	seqSetUp setUp(inputCommands);
	setUp.setOption(weavedResultsDirs, "--resultDirs", "A list of directory with PathWeaver results", true);
	setUp.setOption(mergingBed, "--mergingBed", "take only results matching regions in this bed file", true);
	setUp.processDirectoryOutputName(true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	auto regions = getBed3s(mergingBed);


	VecStr warnings;
	std::string basicInfoHeader;
	for(const auto & resDir : weavedResultsDirs){
		bfs::path allFinalFasta = njh::files::make_path(resDir, "final", "allFinal.fasta");
		bfs::path allBasicInfo  = njh::files::make_path(resDir, "final", "basicInfoPerRegion.tab.txt");
		bfs::path coiCounts = njh::files::make_path(resDir, "final", "coiCounts.tab.txt");
		bfs::path allPartialFasta = njh::files::make_path(resDir, "partial", "allPartial.fasta");
		auto addWarningsPRN = [&warnings](const bfs::path & p){
			if(!bfs::exists(p)){
				warnings.emplace_back(njh::pasteAsStr(p, " needs to exist"));
			}
		};
		addWarningsPRN(allFinalFasta);
		addWarningsPRN(allBasicInfo);
		addWarningsPRN(coiCounts);
		addWarningsPRN(allPartialFasta);
		if(bfs::exists(allBasicInfo)){
			auto firstLine = njh::files::getFirstLine(allBasicInfo);
			if(basicInfoHeader.empty()){
				basicInfoHeader = njh::pasteAsStr(firstLine);
			}else{
				if(basicInfoHeader != firstLine){
					warnings.emplace_back(njh::pasteAsStr(allBasicInfo, " header doesn't match other headers"));
				}
			}
		}
	}

	if(!warnings.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << "\n";
		ss << njh::conToStr(warnings, "\n") << "\n";
		throw std::runtime_error{ss.str()};
	}

	//prep
	bfs::path finalDir =    njh::files::make_path(setUp.pars_.directoryName_, "final");
	bfs::path partialDir  = njh::files::make_path(setUp.pars_.directoryName_, "partial");
	njh::files::makeDir(njh::files::MkdirPar(finalDir));
	njh::files::makeDir(njh::files::MkdirPar(partialDir));
	bfs::path allFinalFasta =   njh::files::make_path(setUp.pars_.directoryName_, "final", "allFinal.fasta");
	bfs::path allBasicInfo  =   njh::files::make_path(setUp.pars_.directoryName_, "final", "basicInfoPerRegion.tab.txt");
	bfs::path coiCounts =       njh::files::make_path(setUp.pars_.directoryName_, "final", "coiCounts.tab.txt");
	bfs::path allPartialFasta = njh::files::make_path(setUp.pars_.directoryName_, "partial", "allPartial.fasta");

	SeqOutput allFinalFasta_writer(SeqIOOptions::genFastaOut(allFinalFasta));
	allFinalFasta_writer.openOut();

	SeqOutput allPartialFasta_writer(SeqIOOptions::genFastaOut(allPartialFasta));
	allPartialFasta_writer.openOut();

	OutputStream allBasicInfo_out(allBasicInfo);
	allBasicInfo_out << basicInfoHeader << std::endl;

	OutputStream coiCounts_out(coiCounts);
	coiCounts_out << "coi\tcount" << std::endl;


	std::map<uint32_t, uint32_t> uniqHapCounts;

	std::set<std::string> regionsUIDs;
	for(const auto & reg : regions){
		auto UID = reg->genUIDFromCoords();
		if(njh::in(UID, regionsUIDs)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "already have " << UID << "\n";
			throw std::runtime_error{ss.str()};
		}
		regionsUIDs.emplace(UID);
	}
	std::set<std::string> readInUIDs;
	std::string sample;

	for(const auto & resDir : weavedResultsDirs){
		std::set<std::string> names;
		bfs::path allBasicInfo  = njh::files::make_path(resDir, "final", "basicInfoPerRegion.tab.txt");
		TableReader basicReader(TableIOOpts::genTabFileIn(allBasicInfo, true));
		VecStr row;
		while (basicReader.getNextRow(row)) {
			if("" == sample){
				sample = row[basicReader.header_.getColPos("sample")];
			}else{
				if(sample !=  row[basicReader.header_.getColPos("sample")]){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "previous sample name: " << sample << " doesn't match: " << row[basicReader.header_.getColPos("sample")]<< "\n";
					ss << njh::conToStr(row, "\t") << "\n";
					throw std::runtime_error{ss.str()};
				}
			}
			std::string uid = njh::pasteAsStr(row[basicReader.header_.getColPos("#chrom")], "-",
					row[basicReader.header_.getColPos("start")], "-",
					row[basicReader.header_.getColPos("end")]);
			if(njh::in(uid, regionsUIDs)){
				if(njh::in(uid, readInUIDs)){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << "already read in region: " << uid << "\n";
					ss << njh::conToStr(row, "\t") << "\n";
					throw std::runtime_error{ss.str()};
				}
				allBasicInfo_out << njh::conToStr(row, "\t") << std::endl;
				++uniqHapCounts[njh::StrToNumConverter::stoToNum<uint32_t>(row[basicReader.header_.getColPos("uniqHaps")])];
				readInUIDs.emplace(uid);
				names.emplace(row[basicReader.header_.getColPos("name")]);
			}
		}

		//
		bfs::path allFinalFasta = njh::files::make_path(resDir, "final", "allFinal.fasta");
		if(0 != bfs::file_size(allFinalFasta)){
			SeqInput reader(SeqIOOptions::genFastaIn(allFinalFasta));
			reader.openIn();
			seqInfo seq;
			while(reader.readNextRead(seq)){
				MetaDataInName seqMeta(seq.name_);
				if(njh::in(seqMeta.getMeta("regionUID"), names)){
					allFinalFasta_writer.write(seq);
				}
			}
		}
		bfs::path allPartialFasta = njh::files::make_path(resDir, "partial", "allPartial.fasta");
		if(0 != bfs::file_size(allPartialFasta)){
			SeqInput reader(SeqIOOptions::genFastaIn(allPartialFasta));
			reader.openIn();
			seqInfo seq;
			while(reader.readNextRead(seq)){
				MetaDataInName seqMeta(seq.name_);
				if(njh::in(seqMeta.getMeta("regionUID"), names)){
					allPartialFasta_writer.write(seq);
				}
			}
		}
	}
	for(const auto & coiCount : uniqHapCounts){
		coiCounts_out << coiCount.first << "\t" << coiCount.second << std::endl;
	}

	return 0;

}

} // namespace njhseq
