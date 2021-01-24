/*
 * GatheringUtils.cpp
 *
 *  Created on: Jan 23, 2021
 *      Author: nick
 */



#include "GatheringUtils.hpp"

namespace njhseq {

SeqGatheringFromPathWeaver::gatherSeqsAndSortByTargetRes SeqGatheringFromPathWeaver::gatherSeqsAndSortByTarget(const gatherSeqsAndSortByTargetPars & pars){
		gatherSeqsAndSortByTargetRes ret;


		//key1 = target, key2= sample, value = the sequences
		std::unordered_map<std::string, std::vector<std::shared_ptr<seqInfo>>> allSeqsByTarget;



		njh::concurrent::LockableVec<bfs::path> inputDirQueue(pars.directories);
		std::mutex allSeqsByTargetMut;

		std::function<void()> collectSeqs = [&inputDirQueue,
																				 &allSeqsByTarget, &allSeqsByTargetMut,
																				 &pars,&ret,this](){
			bfs::path inputDir;
			VecStr currentMissingOutput;
			std::unordered_map<std::string, std::vector<std::shared_ptr<seqInfo>>> currentAllSeqsByTarget;
			while(inputDirQueue.getVal(inputDir)){
				auto finalSeqFnp = njh::files::make_path(inputDir,"final", "allFinal.fasta");
				auto coiPerBedLocationFnp = njh::files::make_path(inputDir,"final", "basicInfoPerRegion.tab.txt");
				if(bfs::exists(finalSeqFnp) && 0 != bfs::file_size(finalSeqFnp) && bfs::exists(coiPerBedLocationFnp)){
					std::unordered_map<std::string, uint32_t> readTotals;
					TableReader readTab(TableIOOpts::genTabFileIn(coiPerBedLocationFnp,true));
					VecStr row;
					uint32_t usedTotal = 0;
					while(readTab.getNextRow(row)){
						if("NA" != row[readTab.header_.getColPos("readTotal")]){
							if(!pars.targets.empty() && !njh::in(row[readTab.header_.getColPos("name")], pars.targets)){
								continue;
							}
							if("true" == row[readTab.header_.getColPos("success")]){
								usedTotal+= njh::StrToNumConverter::stoToNum<uint32_t>(row[readTab.header_.getColPos("readTotalUsed")]);
							}
							readTotals[row[readTab.header_.getColPos("name")]] = njh::StrToNumConverter::stoToNum<uint32_t>(row[readTab.header_.getColPos("readTotal")]);
						}
					}
					if(readTotals.empty() || 0 == usedTotal ){
						currentMissingOutput.emplace_back(inputDir.string());
					}else{
						std::unordered_map<std::string, std::vector<std::shared_ptr<seqInfo>>> seqsByTarget;
						std::set<std::string> samplesInPrcoessFiles;
						{
							auto inputOpts = SeqIOOptions::genFastaIn(finalSeqFnp);
							inputOpts.processed_ = true;
							SeqInput reader(inputOpts);
							reader.openIn();
							seqInfo seq;
							while(reader.readNextRead(seq)){
								MetaDataInName seqMeta(seq.name_);
								auto rawTarName = seqMeta.getMeta(corePars_.targetField);
								if(!pars.targets.empty() && !njh::in(rawTarName, pars.targets)){
									continue;
								}
								auto tarName = njh::replaceString(rawTarName, ".", "-");
								//targetKey[tarName] = rawTarName;
								seqsByTarget[tarName].emplace_back(std::make_shared<seqInfo>(seq));
								samplesInPrcoessFiles.emplace(seqMeta.getMeta(corePars_.sampleField));
							}
						}
						if(pars.addPartial){
							auto partialSeqFnp = njh::files::make_path(inputDir,"partial", "allPartial.fasta");
							if(bfs::exists(partialSeqFnp) && 0 != bfs::file_size(partialSeqFnp)){
								auto inputOpts = SeqIOOptions::genFastaIn(partialSeqFnp);
								inputOpts.processed_ = true;
								SeqInput reader(inputOpts);
								reader.openIn();
								seqInfo seq;

								while(reader.readNextRead(seq)){
									MetaDataInName seqMeta(seq.name_);
									auto rawTarName = seqMeta.getMeta(corePars_.targetField);
									if(seqMeta.containsMeta("trimStatus") && "true" == seqMeta.getMeta("trimStatus")){
										if(!pars.targets.empty() && !njh::in(rawTarName, pars.targets)){
											continue;
										}
										auto tarName = njh::replaceString(rawTarName, ".", "-");
										//targetKey[tarName] = rawTarName;
										seqsByTarget[tarName].emplace_back(std::make_shared<seqInfo>(seq));
										samplesInPrcoessFiles.emplace(seqMeta.getMeta(corePars_.sampleField));
									}
								}
							}
						}
						if(1 != samplesInPrcoessFiles.size()){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << " found more than 1 sample name in " << inputDir << ", found: " << njh::conToStr(samplesInPrcoessFiles, ",")<< "\n";
							throw std::runtime_error{ss.str()};
						}
						for( auto & tar : seqsByTarget){
							double totalOfCountField = 0;
							for(const auto & seq : tar.second){
								MetaDataInName seqMeta(seq->name_);
								if("reads" == corePars_.countField){
									totalOfCountField += seq->cnt_;
								}else{
									totalOfCountField += seqMeta.getMeta<double>(corePars_.countField);
								}
							}
							std::vector<std::shared_ptr<seqInfo>> tarSeqs;
							for(auto & seq : tar.second){
								MetaDataInName seqMeta(seq->name_);
								double count = 0;
								if("reads" == corePars_.countField){
									count = seq->cnt_;
								}else{
									count = seqMeta.getMeta<double>(corePars_.countField);
								}
								if(nullptr != corePars_.meta){
									auto sampleName = seqMeta.getMeta(corePars_.sampleField);
									auto newMeta = corePars_.meta->getMetaForSample(sampleName, njh::getVecOfMapKeys(corePars_.meta->groupData_));
									//newMeta.addMeta("sample", sampleName);
									seqMeta.addMeta(newMeta, true);
									seqMeta.resetMetaInName(seq->name_);
								}
								njh::addVecToSet(getVectorOfMapKeys(seqMeta.meta_), ret.allMetaFields);
								seq->cnt_ = round((count/totalOfCountField) * readTotals[seqMeta.getMeta(corePars_.targetField)]);
								seq->frac_ = count/totalOfCountField;
								seq->updateName();
								tarSeqs.emplace_back(seq);
							}
							addOtherVec(currentAllSeqsByTarget[tar.first], tarSeqs);
						}
					}
				} else {
					currentMissingOutput.emplace_back(inputDir.string());
				}
			}
			{
				std::lock_guard<std::mutex> lock(allSeqsByTargetMut);
				njh::addVecToSet(currentMissingOutput, ret.missingOutput);
				for(auto & tarSeqs : currentAllSeqsByTarget){
					addOtherVec(allSeqsByTarget[tarSeqs.first], tarSeqs.second);
				}
			}
		};

		njh::concurrent::runVoidFunctionThreaded(collectSeqs, corePars_.numThreads);
		{
			auto allSeqOpts = SeqIOOptions::genFastaOut(pars.allSeqFnp);
			allSeqOpts.out_.overWriteFile_ = corePars_.overWrite;
			SeqOutput writer(allSeqOpts);
			writer.openOut();
			uint32_t seqCounts = 0;
			auto targetKeys = getVectorOfMapKeys(allSeqsByTarget);
			njh::sort(targetKeys);

			for(const auto & target : targetKeys){
				ret.seqsLocations[target] = std::make_pair(seqCounts, allSeqsByTarget[target].size());
				writer.write(allSeqsByTarget[target]);
				seqCounts += allSeqsByTarget[target].size();
				allSeqsByTarget[target].clear();
				allSeqsByTarget[target].resize(0);
			}
			ret.allSeqFnp = allSeqOpts.out_.outName();
		}

		SeqInput::buildIndex(SeqIOOptions::genFastaIn(ret.allSeqFnp));


		return ret;
	};


SeqGatheringFromPathWeaver::processedGatherSeqsMetaRes SeqGatheringFromPathWeaver::processedGatherSeqsMeta(
		const processedGatherSeqsMetaPars & pars,
		const gatherSeqsAndSortByTargetRes & gatheredRes){
	processedGatherSeqsMetaRes ret;
	if (!gatheredRes.allMetaFields.empty() &&
			njh::in("ExperimentSample", gatheredRes.allMetaFields)
			&& njh::in("BiologicalSample", gatheredRes.allMetaFields)
			&& njh::in("sample", gatheredRes.allMetaFields)) {
		auto tarKeys = getVectorOfMapKeys(gatheredRes.seqsLocations);
		njh::sort(tarKeys);
		njh::concurrent::LockableQueue<std::string> tarKeysQueue(tarKeys);
		std::mutex allSamplesMut;
		auto processAllSeqFnp = njh::files::prependFileBasename(gatheredRes.allSeqFnp, "processed_");
		ret.allSeqFnp = processAllSeqFnp;
		uint32_t seqCount = 0;
		auto allSeqOpts = SeqIOOptions::genFastaOut(ret.allSeqFnp);
		SeqOutput writer(allSeqOpts);
		writer.openOut();
		auto seqInputOpts = SeqIOOptions::genFastaIn(gatheredRes.allSeqFnp);
		std::function<void()> filterBasedOnMeta = [&tarKeysQueue,
																							 &pars, &gatheredRes,&ret,
																							 &writer,&seqCount,
																							 &allSamplesMut,&seqInputOpts,
																							 this](){
			SeqInput reader(seqInputOpts);

			std::string tarKey="";
			std::set<std::string> currentAllSamples;
			std::unordered_map<std::string,std::unordered_map<std::string, VecStr>> currentFailedToCombine;
			while(tarKeysQueue.getVal(tarKey)){
				auto targetSeqs = reader.getReads<seqInfo>(gatheredRes.seqsLocations.at(tarKey).first, gatheredRes.seqsLocations.at(tarKey).second);

				std::vector<uint32_t> allSeqsPositions(targetSeqs.size());
				njh::iota<uint32_t>(allSeqsPositions, 0);
				for(const auto & metaField : VecStr{"ExperimentSample", "BiologicalSample"}){
					MetaFieldSeqFilterer::MetaFieldSeqFiltererPars filterPars;
					filterPars.tieBreakerPreferenceField_ = "PreferredSample";
					filterPars.groupingMetaField_ = corePars_.sampleField;
					filterPars.filteringMetaField_ = metaField;
					filterPars.keepCommonSeqsWhenFiltering_ = pars.keepCommonSeqsWhenFiltering;
					MetaFieldSeqFilterer metaFilterer(filterPars);
					auto res = metaFilterer.filterSeqs(targetSeqs, allSeqsPositions);
					allSeqsPositions = res.filteredAllSeqs;
					if(!res.failedGroups.empty()){
						currentFailedToCombine[tarKey][metaField] = res.failedGroups;
					}
				}
				for(const auto & filtPos : allSeqsPositions){
					MetaDataInName meta(targetSeqs[filtPos].name_);
					meta.removeMeta("ExperimentSample");
					meta.removeMeta(corePars_.sampleField);
					auto sampleName = meta.getMeta("BiologicalSample");
					meta.addMeta("sample", sampleName);
					currentAllSamples.emplace(sampleName);
					meta.removeMeta("BiologicalSample");
					meta.resetMetaInName(targetSeqs[filtPos].name_);
				}
				{
					std::lock_guard<std::mutex> lock(allSamplesMut);
					ret.seqsLocations[tarKey] = std::make_pair(seqCount, allSeqsPositions.size());
					for(const auto & filtPos : allSeqsPositions){
						writer.write(targetSeqs[filtPos]);
					}
					ret.targetKey[tarKey] = MetaDataInName(targetSeqs.front().name_).getMeta(corePars_.targetField);
					seqCount += allSeqsPositions.size();
				}
			}
			{
				std::lock_guard<std::mutex> lock(allSamplesMut);
				ret.allSamples.insert(currentAllSamples.begin(), currentAllSamples.end());
				for(const auto & tar : currentFailedToCombine){
					ret.failedToCombine[tar.first] = tar.second;
				}
			}
		};

		njh::concurrent::runVoidFunctionThreaded(filterBasedOnMeta, corePars_.numThreads);

		writer.closeOut();
		SeqInput::buildIndex(SeqIOOptions::genFastaIn(ret.allSeqFnp));

	} else {
		ret.seqsLocations = gatheredRes.seqsLocations;
		ret.allSeqFnp = gatheredRes.allSeqFnp;
	}
	if (!gatheredRes.allMetaFields.empty() &&
			njh::in("ExperimentSample", gatheredRes.allMetaFields)
			&& njh::in("BiologicalSample", gatheredRes.allMetaFields)
			&& njh::in("sample", gatheredRes.allMetaFields)) {
		//redo the meta
		table metaTab(corePars_.meta->groupingsFile_, "\t", true);
		const auto inputMetaColumns = metaTab.columnNames_;
		metaTab.deleteColumn("ExperimentSample");
		metaTab.deleteColumn("sample");
		auto bioSampColumnPos = getPositionsOfTarget(metaTab.columnNames_, "BiologicalSample");
		metaTab.columnNames_[bioSampColumnPos.front()] = "sample";
		metaTab.setColNamePositions();

		if(njh::in(std::string("PreferredSample"), metaTab.columnNames_)){
			auto nonPreferredSamples = metaTab.extractByComp("PreferredSample", [](const auto & val){return "FALSE" == val;});
			metaTab = metaTab.extractByComp("PreferredSample", [](const auto & val){return "TRUE" == val;});
			metaTab.deleteColumn("PreferredSample");
			metaTab.setColNamePositions();
			nonPreferredSamples.deleteColumn("PreferredSample");
			nonPreferredSamples = nonPreferredSamples.getUniqueRows();
			nonPreferredSamples.setColNamePositions();

			auto allInputBiologicalSamples = metaTab.getColumnLevels("sample");
			std::unordered_map<std::string, uint32_t> countsOfSampleMissingFromMeta;
			for(const auto & row : nonPreferredSamples){
				if(!njh::in(row[nonPreferredSamples.getColPos("sample")], allInputBiologicalSamples)){
					++countsOfSampleMissingFromMeta[row[nonPreferredSamples.getColPos("sample")]];
				}
			}
			VecStr samplesWithMultipleInputs;
			for(const auto & count : countsOfSampleMissingFromMeta){
				if(count.second > 1){
					samplesWithMultipleInputs.emplace_back(count.first);
				} //createSharedPathwaysFromRefSeqs
			}
			if(!samplesWithMultipleInputs.empty()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error " << " the following samples have multiple entries: " << njh::conToStr(samplesWithMultipleInputs, ",")<< "\n";
				throw std::runtime_error{ss.str()};
			}
			for(const auto & row : nonPreferredSamples){
				if(!njh::in(row[nonPreferredSamples.getColPos("sample")], allInputBiologicalSamples)){
					metaTab.addRow(row);
				}
			}
		}
		metaTab = metaTab.getUniqueRows();
		//add in missing samples that had no preferred sample
		OutputStream metaOut(pars.outMetaFnp);
		metaTab.outPutContents(metaOut, "\t");
	}else{
		if(nullptr != corePars_.meta && "" != corePars_.meta->groupingsFile_){
			bfs::copy_file(corePars_.meta->groupingsFile_, pars.outMetaFnp);
		}
	}
	return ret;
};

}  // namespace njhseq

