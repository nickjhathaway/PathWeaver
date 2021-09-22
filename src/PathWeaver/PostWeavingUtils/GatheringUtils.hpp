#pragma once

/*
 * GatheringUtils.hpp
 *
 *  Created on: Jan 23, 2021
 *      Author: nick
 */


//#include <cppitertools/sorted.hpp>
#include <SeekDeepPrograms/SeekDeepProgram/SeekDeepSetUp.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/Meta.h>




namespace njhseq {



class SeqGatheringFromPathWeaver{
public:

	struct SeqGatheringFromPathWeaverCorePars{
		std::string sampleField;
		std::string targetField;
		std::string countField;
		bool overWrite = false;
		uint32_t numThreads = 1;
		std::shared_ptr<MultipleGroupMetaData> meta;
	};

	SeqGatheringFromPathWeaver(const SeqGatheringFromPathWeaverCorePars & corePars):corePars_(corePars){

	}
	SeqGatheringFromPathWeaverCorePars corePars_;

	struct gatherSeqsAndSortByTargetRes{


		std::set<std::string> allMetaFields;
		std::unordered_map<std::string, std::pair<uint32_t, uint32_t>> seqsLocations;
		std::set<std::string> missingOutput;
		bfs::path allSeqFnp;
		std::set<std::string> allSamples;
		std::unordered_map<std::string, std::string> targetKey;//a key because "." can't be analyses names

	};

	struct gatherSeqsAndSortByTargetPars{
		bfs::path allSeqFnp;
		std::vector<bfs::path> directories;
		std::set<std::string> targets;
		bool addPartial = false;
		std::unordered_map<std::string, std::vector<seqWithKmerInfo>> trimSeqs;
		uint32_t minInputSeqLen = 25;
	};

	gatherSeqsAndSortByTargetRes gatherSeqsAndSortByTarget(const gatherSeqsAndSortByTargetPars & pars);


	struct processedGatherSeqsMetaRes{
		std::unordered_map<std::string, std::pair<uint32_t, uint32_t>> seqsLocations;
		bfs::path allSeqFnp;

		std::set<std::string> allSamples;
		std::unordered_map<std::string,std::unordered_map<std::string, VecStr>> failedToCombine;
		//std::unordered_map<std::string, std::string> targetKey;//a key because "." can't be analyses names

	};

	struct processedGatherSeqsMetaPars{

		bfs::path allSeqFnp;
		bfs::path outMetaFnp;
		bool keepCommonSeqsWhenFiltering = true;
		bool overWrite = false;
	};

	processedGatherSeqsMetaRes processedGatherSeqsMeta(
			const processedGatherSeqsMetaPars & pars,
			const gatherSeqsAndSortByTargetRes & gatheredRes);

};







}  // namespace njhseq

