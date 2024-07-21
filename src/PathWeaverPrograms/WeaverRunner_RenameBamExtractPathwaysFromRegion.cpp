//
// Created by Nicholas Hathaway on 7/19/24.
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

int WeaverRunner::RenameBamExtractPathwaysFromRegion(
		const njh::progutils::CmdArgs & inputCommands) {

	bool subset = false;
	bfs::path weavedResultsDir;

	bfs::path renamingFile;

	seqSetUp setUp(inputCommands);
	setUp.setOption(subset, "--subset", "subset the file with the renaming file, will only take targets that are in the renaming file");

	setUp.setOption(weavedResultsDir, "--resultDir", "The directory with PathWeaver results", true);
	setUp.setOption(renamingFile, "--renamingFile", "Renaming File, needs to have columns 'new' and 'old', old is the name of the target region and new is the new name", true);
	setUp.processDirectoryOutputName(true);

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);


	{
		VecStr warnings;
		bfs::path allFinalFasta = njh::files::make_path(weavedResultsDir, "final", "allFinal.fasta");
		bfs::path allBasicInfo  = njh::files::make_path(weavedResultsDir, "final", "basicInfoPerRegion.tab.txt");
		bfs::path coiCounts = njh::files::make_path(weavedResultsDir, "final", "coiCounts.tab.txt");
		bfs::path allPartialFasta = njh::files::make_path(weavedResultsDir, "partial", "allPartial.fasta");
		auto addWarningsPRN = [&warnings](const bfs::path & p){
			if(!bfs::exists(p)){
				warnings.emplace_back(njh::pasteAsStr(p, " needs to exist"));
			}
		};
		addWarningsPRN(allFinalFasta);
		addWarningsPRN(allBasicInfo);
		addWarningsPRN(coiCounts);
		addWarningsPRN(allPartialFasta);
		

		if(!warnings.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "\n";
			ss << njh::conToStr(warnings, "\n") << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	table renamingTab(renamingFile, "\t", true);
	renamingTab.checkForColumnsThrow(VecStr{"old", "new"}, __PRETTY_FUNCTION__);

	table basic_info_table( njh::files::make_path(weavedResultsDir, "final", "basicInfoPerRegion.tab.txt"), "\t", true);

	{
		VecStr warnings;
		auto renamingOldCol= njh::vecToSet(renamingTab.getColumn("old"));
		auto nameCol= njh::vecToSet(basic_info_table.getColumn("name"));


		VecStr targetsInRenamingNotInBasic;
		VecStr targetsInBasicNotInRenaming;
		VecStr shared;
		njh::decompose_sets(renamingOldCol.begin(), renamingOldCol.end(),
												nameCol.begin(), nameCol.end(),
												std::back_inserter(targetsInRenamingNotInBasic),
												std::back_inserter(targetsInBasicNotInRenaming),
												std::back_inserter(shared));
		if(!targetsInRenamingNotInBasic.empty()) {
			warnings.emplace_back(
				njh::pasteAsStr("The following names were found in the renaming file, ", renamingFile, " but not in the basic info file, ", njh::files::make_path(weavedResultsDir, "final", "basicInfoPerRegion.tab.txt"),
					"\n",
					targetsInRenamingNotInBasic)
			);
		}
		if(!targetsInBasicNotInRenaming.empty() && !subset) {
			warnings.emplace_back(
				njh::pasteAsStr("The following names were not found in the renaming file, ", renamingFile, " but were in the basic info file, ", njh::files::make_path(weavedResultsDir, "final", "basicInfoPerRegion.tab.txt"),
					"\n",
					targetsInBasicNotInRenaming)
			);
		}
		if(!warnings.empty()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << "\n";
			ss << njh::conToStr(warnings, "\n") << "\n";
			throw std::runtime_error{ss.str()};
		}
	}

	//key
	std::unordered_map<std::string, std::string> renamingKey;
	for(const auto & row : renamingTab) {
		if(njh::in(row[renamingTab.getColPos("old")], renamingKey)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " already have a key for name: " <<  row[renamingTab.getColPos("old")] << "\n";
			throw std::runtime_error{ss.str()};
		}
		renamingKey[row[renamingTab.getColPos("old")]] = row[renamingTab.getColPos("new")];
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


	//rename basic info
	OutputStream allBasicInfo_out(allBasicInfo);
	for(auto & row : basic_info_table) {
		if(!subset || njh::in(row[basic_info_table.getColPos("name")], renamingKey)) {
			row[basic_info_table.getColPos("name")] = renamingKey[row[basic_info_table.getColPos("name")]];
		}
	}
	basic_info_table.outPutContents(allBasicInfo_out, "\t");

	//renaming shouldn't affect coi counts so just copy over
	bfs::copy(njh::files::make_path(weavedResultsDir, "final", "coiCounts.tab.txt"), coiCounts);


	SeqOutput allFinalFasta_writer(SeqIOOptions::genFastaOut(allFinalFasta));
	allFinalFasta_writer.openOut();
	//
	bfs::path input_allFinalFasta = njh::files::make_path(weavedResultsDir, "final", "allFinal.fasta");
	if(0 != bfs::file_size(input_allFinalFasta)){
		SeqInput reader(SeqIOOptions::genFastaIn(input_allFinalFasta));
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			MetaDataInName seqMeta(seq.name_);
			auto old_regionUID = seqMeta.getMeta("regionUID");
			if(!subset || njh::in(old_regionUID, renamingKey)) {
				seqMeta.addMeta("regionUID", renamingKey[old_regionUID], true);
				seqMeta.resetMetaInName(seq.name_);
				allFinalFasta_writer.write(seq);
			}
		}
	}

	SeqOutput allPartialFasta_writer(SeqIOOptions::genFastaOut(allPartialFasta));
	allPartialFasta_writer.openOut();
	bfs::path input_allPartialFasta = njh::files::make_path(weavedResultsDir, "partial", "allPartial.fasta");
	if(0 != bfs::file_size(input_allPartialFasta)){
		SeqInput reader(SeqIOOptions::genFastaIn(input_allPartialFasta));
		reader.openIn();
		seqInfo seq;
		while(reader.readNextRead(seq)){
			MetaDataInName seqMeta(seq.name_);
			auto old_regionUID = seqMeta.getMeta("regionUID");
			if(!subset || njh::in(old_regionUID, renamingKey)) {
				seqMeta.addMeta("regionUID", renamingKey[old_regionUID], true);
				seqMeta.resetMetaInName(seq.name_);
				allPartialFasta_writer.write(seq);
			}
		}
	}

	return 0;

}

} // namespace njhseq
