/*
 * WeaverRunner_pathFinders_BamExtractPathawaysFromRegion.cpp
 *
 *  Created on: May 20, 2018
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

#include "WeaverRunner.hpp"
#include <TwoBit.h>

#include <njhseq/BamToolsUtils.h>
#include <njhseq/objects/BioDataObject.h>
#include <njhseq/GenomeUtils.h>

#include <njhcpp/concurrency/LockableJsonLog.hpp>
#include <njhseq/objects/seqContainers.h>

#include "PathWeaver/objects/bam/RegionInvestigatorInBam.hpp"
#include "PathWeaver/objects/Meta/CountryMetaData.hpp"
#include "PathWeaver/seqToolsUtils/HaplotypeLocator.hpp"
#include "PathWeaver/PathFinding.h"

namespace njhseq {

int WeaverRunner::BamExtractPathawaysFromRegion(
		const njh::progutils::CmdArgs & inputCommands) {
	BamRegionInvestigator::BamRegionInvestigatorPars brInvestPars;

	BamRegionInvestigator::BamCountSpecficRegionsPars spanningReadsPar;
	BioCmdsUtils::LastZPars lzPars;
	lzPars.identity = 80;

	HaploPathFinder::ExtractParams masterPars;
	bfs::path bedFile = "";
	bfs::path outputDir;
	bfs::path refDir = "";

	bool keepTemporaryFiles = false;
	bool writeOutBaseCoveragePerSeq = false;
	bool reDetermineReadCounts = false;

	bool filterOnBaseCoverage = false;
	double baseCoverageCutOff = 10;


	bool ignoreFailedTrim = false;

	bool expandEndRegions = false;
	uint32_t expandRegionBy = 25;


	bool noTrim = false;

	uint32_t subNumThreads = 1;

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	masterPars.pFinderPars_.autoDetermineKCuts = true;
	setUp.setOption(noTrim, "--noTrim", "No trimming");
	setUp.setOption(subNumThreads, "--subNumThreads",
			"Sub Num Threads for each region, total number of threads will be --numThreads * --subNumThreads, e.g. total =100 for --numThreads 10 --subNumThreads 10");

	//optimization opts
	masterPars.pFinderPars_.setOptimizationParameters(setUp);

	//pair stitching parameters
	masterPars.pFinderPars_.reOrientPairsForStitching = true;
	masterPars.pFinderPars_.setStitchParameters(setUp);

	//group info
	masterPars.pFinderPars_.setAddingGroupInfoOpts(setUp);

	//possible and major haplotypes
	masterPars.pFinderPars_.setPossibleHapsOpts(setUp);

	//pre-process
	//tandem repeat
	masterPars.pFinderPars_.setTandemRepeatHandlingOptions(setUp);
	//duplicate seqs handling
	masterPars.pFinderPars_.removeDuplicatedSequences_ = false;
	masterPars.pFinderPars_.setDuplicatedSeqHandlingOpts(setUp);

	//quality trimming and filtering
	masterPars.pFinderPars_.setQualityTrimAndFiltOpts(setUp);

	//graph options
	//splitting
	masterPars.pFinderPars_.setFurtherSplittingOpts(setUp);
	//processing nodes, onebase indel nodes, remove headless and tailless nodes
	masterPars.pFinderPars_.setNodeProcessingOpts(setUp);

	//genome options
	bool needGenome = true;
	masterPars.setGenomeOpts(setUp, needGenome);
	//run modes
	masterPars.pFinderPars_.setRunningOpts(setUp);
	//outliers
	masterPars.pFinderPars_.setOutlierRemovalOptions(setUp);

	//post processing of collapsing contigs;
	masterPars.pFinderPars_.setPostProcessContigHandlingOpts(setUp);

	//recuritment
	masterPars.bamExtractPars_.keepMarkedDuplicate_ = false;
	masterPars.bamExtractPars_.filterOffLowEntropyOrphansRecruits_ = false;
	//masterPars.bamExtractPars_.softClipPercentageCutOff_ = 0.30;
	masterPars.bamExtractPars_.softClipPercentageCutOff_ = 1.00;
	masterPars.bamExtractPars_.removeImproperPairs_ = true;
	masterPars.bamExtractPars_.keepImproperMateUnmapped_ = true;
	masterPars.setBamExtractOpts(setUp);

	//adding meta to the final seqs
	setUp.setOption(masterPars.metaDataFnp, "--metaDataFnp", "Name of the meta data fnp");

	bool doNotExapndRegions = false;
	setUp.setOption(doNotExapndRegions, "--doNotExapndRegions", "Do Not Expand Ends");
	expandEndRegions = !doNotExapndRegions;
	setUp.setOption(expandRegionBy,   "--expandRegionBy", "Expand Region By this length for the inital pull down");

	//// Specific for here
	bool dontAddSeqOfSingleHeadAndTailSeqs = false;
	setUp.setOption(dontAddSeqOfSingleHeadAndTailSeqs, "--doNotAddHeadAndTailSeqs", "Don't Add Head And Tail sequences when creating output nodes");
	masterPars.pFinderPars_.addSeqOfSingleHeadAndTailSeqs_ = !dontAddSeqOfSingleHeadAndTailSeqs;

	setUp.setOption(ignoreFailedTrim, "--ignoreFailedTrim", "Ignore failed trim status sequences");

	setUp.setOption(masterPars.pFinderPars_.trimToCircularGenome, "--trimToCircularGenome", "Trim contigs to the reference assuming ref is circular");

	//getting base coverage
	setUp.setOption(filterOnBaseCoverage, "--filterOnBaseCoverage", "Filter On Base Coverage");
	setUp.setOption(baseCoverageCutOff, "--baseCoverageCutOff", "Base Coverage Cut Off");
	setUp.setOption(reDetermineReadCounts, "--reDetermineReadCounts", "Re-determine Read Counts by mapping the recruited sequences back to final sequences");
	setUp.setOption(writeOutBaseCoveragePerSeq, "--writeOutBaseCoveragePerSeq", "Write Out Base Coverage Per Seq");
	setUp.setOption(keepTemporaryFiles, "--keepTemporaryFiles", "Keep Temporary Files, by default only the final calls are kept");


	setUp.setOption(brInvestPars.mapQualityCutOffForMultiMap_, "--mapQualityCutOffForMultiMapForCovInfo", "map Quality Cut Off For Coverage");
	setUp.setOption(brInvestPars.mapQualityCutOff_, "--mapQualityCutOffForCov", "map Quality Cut Off For Coverage");

	setUp.setOption(lzPars.identity, "--lastzIdentity", "Identity used for lastz when getting reference sequences");
	setUp.setOption(lzPars.coverage, "--lastzCoverage", "Coverage used for lastz when getting reference sequences");
	setUp.setOption(bedFile, "--bed", "Bed file of multiple regions to process", true);
	//setUp.setOption(refDir, "--refDir", "Reference directory created by PathWeaver extractRefSeqsFromGenomes", true);
	setUp.setOption(refDir, "--refDir", "Reference directory created by elucidator extractRefSeqsFromGenomes", true);

	setUp.processDefaultReader( { "--bam" }, true);
	setUp.processDirectoryOutputName(true);
	bool writeOutLogs = false;
	setUp.setOption(writeOutLogs, "--writeOutLogs", "Write Out Logs");

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);




	spanningReadsPar.numThreads = masterPars.pFinderPars_.numThreads;
	spanningReadsPar.countDuplicates = masterPars.bamExtractPars_.keepMarkedDuplicate_;


	masterPars.pFinderPars_.verbose = setUp.pars_.verbose_;
	masterPars.pFinderPars_.debug = setUp.pars_.debug_;

	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);
	masterPars.region_ = regions.front();
	sortGRegionsByStart(regions);


	std::unique_ptr<MultipleGroupMetaData> meta;
	auto sampName = getPossibleSampleNameFromFnp(
					setUp.pars_.ioOptions_.firstName_);
	if ("" != masterPars.metaDataFnp) {

		meta = std::make_unique<MultipleGroupMetaData>(masterPars.metaDataFnp,
				std::set<std::string> { sampName });
	} else {
		meta = std::make_unique<MultipleGroupMetaData>(table { VecStr { "sample",
				"holder" } }, std::set<std::string> { "empty" });
	}

	auto genomeFnp2bit = masterPars.genomeFnp;
	genomeFnp2bit.replace_extension(".2bit");

	njh::concurrent::LockableQueue<GenomicRegion> regionsQueue(regions);

	njh::LockableJsonLog jLog(njh::files::make_path(setUp.pars_.directoryName_, "extractionLog.json"));

	masterPars.outputDir = setUp.pars_.directoryName_;
	bool regionFailed = false;
	std::stringstream ss;
	ss << __PRETTY_FUNCTION__ << ", error, couldn't find ref seqs for the following directories: " << "\n";
	for(const auto & region : regions){
		auto refSeqs = njh::files::make_path(refDir, region.createUidFromCoordsStrand(), "allRefs.fasta");
		if(!bfs::exists(refSeqs)){
			ss << 	region.createUidFromCoordsStrand() << "\n";
		}
	}
	if(regionFailed){
		throw std::runtime_error{ss.str()};
	}
	std::string coverageField = "estimatedPerBaseCoverage";
	if(reDetermineReadCounts){
		coverageField = "meanBaseCoverage";
	}
	//
	bfs::path finalDirectory = njh::files::makeDir(masterPars.outputDir, njh::files::MkdirPar("final"));
	bfs::path partialDirectory = njh::files::makeDir(masterPars.outputDir, njh::files::MkdirPar("partial"));

	//get coverage and spanning read info
	brInvestPars.numThreads_ = masterPars.pFinderPars_.numThreads;
	brInvestPars.countDups_ = masterPars.bamExtractPars_.keepMarkedDuplicate_;
	BamRegionInvestigator brInvestor(brInvestPars);
	auto regInfos = brInvestor.getCoverageAndFullSpanningReads(
			setUp.pars_.ioOptions_.firstName_, regions, spanningReadsPar, genomeFnp2bit);

	brInvestor.writeBasicInfo(regInfos, sampName, OutOptions(njh::files::make_path(finalDirectory, "perBaseCoveragePerRegion.bed")));

	std::unordered_map<std::string,
			std::shared_ptr<BamRegionInvestigator::RegionInfo>> regInfosByUID;
	for (const auto & p : regInfos) {
		regInfosByUID[p->region_.createUidFromCoordsStrand()] = p;
	}

	std::function<void(HaploPathFinder::ExtractParams)> extractPathway =
			[&meta,&regionsQueue,&jLog,&refDir,
			 &finalDirectory, &writeOutBaseCoveragePerSeq,&baseCoverageCutOff,
			 &filterOnBaseCoverage,&reDetermineReadCounts,
			 &setUp,&coverageField,
			 &ignoreFailedTrim, &partialDirectory,&genomeFnp2bit,
			 &writeOutLogs,&expandEndRegions,&expandRegionBy,
			 &regInfosByUID,&noTrim,&subNumThreads](HaploPathFinder::ExtractParams inputPars) {
				GenomicRegion region;
				BamExtractor bExtractor(setUp.pars_.verbose_);
				bExtractor.debug_ = setUp.pars_.debug_;
				const bfs::path bamFnp = setUp.pars_.ioOptions_.firstName_;
				TwoBit::TwoBitFile tReader(genomeFnp2bit);
				while(regionsQueue.getVal(region)) {
					HaploPathFinder::ExtractParams parsForRegion = inputPars;
					parsForRegion.pFinderPars_.numThreads = subNumThreads;
					//region.setUidWtihCoords();
					auto regionDir = njh::files::makeDir(parsForRegion.outputDir,
							njh::files::MkdirPar(region.createUidFromCoordsStrand()));
//					auto refAlignsDir = njh::files::makeDir(regionDir,
//							njh::files::MkdirPar("refAlignments"));
					uint64_t maxLen = 350;

					auto refSeqsFnp = njh::files::make_path(refDir, region.createUidFromCoordsStrand(), "allRefs.fasta");
					if(!noTrim){
						if(bfs::exists(refSeqsFnp)){
							auto refSeqsOpts = SeqIOOptions::genFastaIn(refSeqsFnp);
							parsForRegion.pFinderPars_.inputSeqs = SeqInput::getSeqVec<seqInfo>(refSeqsOpts, maxLen);
							parsForRegion.pFinderPars_.trimToInputSeqs = true;
							SeqOutput::write(parsForRegion.pFinderPars_.inputSeqs,
									SeqIOOptions::genFastaOut(
											njh::files::make_path(regionDir, "allRefs.fasta")));
						} else {
							parsForRegion.pFinderPars_.inputSeqs = {region.extractSeq(tReader)};
							parsForRegion.pFinderPars_.trimToInputSeqs = true;
							SeqOutput::write(parsForRegion.pFinderPars_.inputSeqs,
									SeqIOOptions::genFastaOut(
											njh::files::make_path(regionDir, "allRefs.fasta")));
							readVec::getMaxLength(parsForRegion.pFinderPars_.inputSeqs, maxLen);
						}
					} else {
						parsForRegion.pFinderPars_.trimToInputSeqs = false;
					}

					Json::Value logValue;
					std::string sampName = getPossibleSampleNameFromFnp(bamFnp);
					auto extractionDir = njh::files::makeDir(regionDir, njh::files::MkdirPar(sampName + "_extraction"));
					try {
						auto extractedOpts = setUp.pars_.ioOptions_;
						extractedOpts.out_.outFilename_ = njh::files::make_path(extractionDir, "extracted");
						auto regionUsed = region;
						if(expandEndRegions){
							regionUsed.start_ = regionUsed.start_ > expandRegionBy ? regionUsed.start_ - expandRegionBy : 0;
							regionUsed.end_ +=  expandRegionBy;
						}
						auto regionExtracted = bExtractor.extractReadsWtihCrossRegionMapping(
								extractedOpts, {regionUsed}, parsForRegion.bamExtractPars_);
						{
							OutOptions extractedReadsStatsFileOpts(njh::files::make_path(extractionDir, "extractedReadsStats.tab.txt"));
							OutputStream extractedReadsStatsFileOut(extractedReadsStatsFileOpts);
							regionExtracted.log(extractedReadsStatsFileOut, setUp.pars_.ioOptions_.firstName_);
						}

						if(parsForRegion.pFinderPars_.autoDetermineKCuts){
							auto perBaseCoverage = regInfosByUID.at(region.createUidFromCoordsStrand())->getPerBaseCov();
							if(perBaseCoverage > 100 || parsForRegion.pFinderPars_.forceDetermineKCuts){
								parsForRegion.pFinderPars_.kmerKOcurrenceCutOffs.clear();
								for(const auto & perc : parsForRegion.pFinderPars_.autoDetermineKCutsPercentages){
									parsForRegion.pFinderPars_.kmerKOcurrenceCutOffs.emplace_back(round(perBaseCoverage * perc));
								}
							}
						}

						auto PWRes = PathFinderFromSeqs(regionExtracted, regionDir, sampName, parsForRegion.pFinderPars_, meta);
						logValue = PWRes.log_;
						auto outputOpts = SeqIOOptions::genFastaIn(njh::files::make_path(regionDir, sampName, "output_aboveCutOff.fasta"), true);

						if(outputOpts.inExists()){
							if(logValue.isMember("bestTotalInputReads")){
								regInfosByUID.at(region.createUidFromCoordsStrand())->totalFinalReads_ = logValue["bestTotalInputReads"].asUInt();
							}
							SeqInput reader(outputOpts);
							if(!reader.isFirstEmpty()){
								std::vector<std::shared_ptr<seqInfo>> outputSeqs = reader.readAllReadsPtrs<seqInfo>();
								double medianAvgRedLen = logValue["readLengthMedian"].asDouble();
								//std::cout << "outputSeqs.size() : " << outputSeqs.size() << std::endl;
								//get trim status
								bool failFilters = outputSeqs.empty();
								std::vector<uint32_t> removeSeqPositions;
								if(noTrim){
									failFilters = true;
								}else{
									for(const auto  seqPos : iter::range(outputSeqs.size())){
										const auto & seq = outputSeqs[seqPos];
										MetaDataInName seqMeta(seq->name_);
										//std::cout << seq->name_ << " " << njh::colorBool(seqMeta.getMeta<bool>("trimStatus"))<< std::endl;
										if(!seqMeta.getMeta<bool>("trimStatus")){
											if(ignoreFailedTrim){
												removeSeqPositions.push_back(seqPos);
											}else{
												failFilters = true;
												break;
											}
										}
									}
								}

								if(ignoreFailedTrim){
									for(const auto pos : iter::reversed(removeSeqPositions)){
										outputSeqs.erase(outputSeqs.begin() + pos);
									}
									failFilters = outputSeqs.empty();
								}


								//add region meta info
								for(auto & seq : outputSeqs){
									MetaDataInName seqMeta(seq->name_);
									seqMeta.addMeta("regionUID", region.uid_);
									seqMeta.addMeta("regionCoords", region.createUidFromCoordsStrand());
									seqMeta.resetMetaInName(seq->name_);
								}

								if(!failFilters){

									auto finalSeqsOnClusters = cluster::convertVectorToClusterVector<cluster>(outputSeqs);
									auto finalKeptSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(regionDir, sampName, "finalKeptTrimedSeqs.fasta"));
									SeqOutput::write(finalSeqsOnClusters, finalKeptSeqOpts);
									auto finalInfoRegionDir = njh::files::makeDir(finalDirectory,
																njh::files::MkdirPar(region.createUidFromCoordsStrand()));
									if(reDetermineReadCounts){
										readVec::getMaxLength(finalSeqsOnClusters, maxLen);

										auto r1Fnp =      njh::files::make_path(regionDir, sampName, "filteredExtractedPairs_R1.fastq");
										auto r2Fnp =      njh::files::make_path(regionDir, sampName, "filteredExtractedPairs_R2.fastq");
										auto singlesFnp = njh::files::make_path(regionDir, sampName, "filteredSingles.fastq");
										auto getMaxLenFromFile = [&maxLen](const bfs::path & seqFnp){
											if(bfs::exists(seqFnp)){
												SeqInput reader(SeqIOOptions::genFastqIn(seqFnp));
												reader.openIn();
												seqInfo seq;
												while(reader.readNextRead(seq)){
													readVec::getMaxLength(seq, maxLen);
												}
											}
										};
										getMaxLenFromFile(r1Fnp);
										getMaxLenFromFile(r2Fnp);
										getMaxLenFromFile(singlesFnp);
										bool countEndGaps = true;
										aligner alignerObj(maxLen, gapScoringParameters(5,1,5,1,5,1), substituteMatrix(1,-1), countEndGaps);
										alignerObj.weighHomopolymers_ = true;
										comparison allowableErrors;
										allowableErrors.oneBaseIndel_ = 0.2;
										allowableErrors.twoBaseIndel_ = 0.2;
										allowableErrors.largeBaseIndel_ = 0.2;



										std::vector<refVariants> refVariationInfo;
										for (const auto & refPos : iter::range(finalSeqsOnClusters.size())) {
											refVariationInfo.emplace_back(finalSeqsOnClusters[refPos].seqBase_);
										}
										for (const auto & refPos : iter::range(finalSeqsOnClusters.size())) {
											for (const auto & refSubPos : iter::range(finalSeqsOnClusters.size())) {
												if (refPos == refSubPos) {
													continue;
												}
												refVariationInfo[refPos].addVariant(finalSeqsOnClusters[refSubPos].seqBase_,
														alignerObj, false);
											}
										}
										alignerObj.countEndGaps_ = false;
										alignerObj.setGapScoring(gapScoringParameters(5,1,0,0,0,0));
										allowableErrors.lqMismatches_ = 2;
										allowableErrors.hqMismatches_ = 1;
										allowableErrors.distances_.query_.coverage_ = parsForRegion.bamExtractPars_.percInRegion_;
										uint64_t minLen = 0;
										readVec::getMinLength(finalSeqsOnClusters, minLen);
										if(minLen < medianAvgRedLen){
											allowableErrors.distances_.query_.coverage_ = minLen/medianAvgRedLen;
										}


										std::vector<kmerInfo> refInfos;
										double matchingKmerCutOff = 0.50;
										uint32_t matchingKmerLen = 5;


										for (const auto& clus : finalSeqsOnClusters){
											refInfos.emplace_back(clus.seqBase_.seq_, matchingKmerLen, false);
										}

										uint32_t unmappable    = 0;
										uint32_t indeterminate = 0;
										uint32_t multiMapping  = 0;
										std::unordered_map<std::string, std::set<std::string>> uniqueSeqs;
										std::unordered_map<std::string, std::vector<seqInfo>> multiMappingSeqs;
										std::set<std::string> allReadInSamples;
										std::unordered_map<std::string, std::map<uint32_t, std::unordered_map<std::string, double>>> baseCoverrage;

										auto increaseMapCounts = [&unmappable,&indeterminate,&uniqueSeqs,&multiMappingSeqs, &refInfos,&finalSeqsOnClusters,&alignerObj,
																						 &allowableErrors,&refVariationInfo,&matchingKmerLen,&allReadInSamples,&multiMapping,
																						 &baseCoverrage](const seqInfo & reMapSeq, const double matchingKmerCutOff){
											allReadInSamples.emplace(reMapSeq.name_);
											//best score to be used to get all the references that map within a certain distance to this latter
											double bestScore = 0;
											kmerInfo reMappingSeqInfo(reMapSeq.seq_, matchingKmerLen, false);
											std::vector<uint64_t> bestRefs;
											uint32_t subSize = len(reMapSeq) * allowableErrors.distances_.query_.coverage_;
											for (const auto & refPos : iter::range(finalSeqsOnClusters.size())) {
	//											if(reMappingSeqInfo.compareSubKmersToFull(refInfos[refPos],
	//															0,
	//															subSize).second < matchingKmerCutOff &&
	//
	//													reMappingSeqInfo.compareSubKmersToFull(refInfos[refPos],
	//															len(reMapSeq) - subSize,
	//															subSize).second < matchingKmerCutOff) {
	//
	//												continue;
	//											}
												if(reMappingSeqInfo.compareKmers(refInfos[refPos]).second < matchingKmerCutOff &&
													 reMappingSeqInfo.compareSubKmersToFull(refInfos[refPos], 0, subSize).second < matchingKmerCutOff &&
													 reMappingSeqInfo.compareSubKmersToFull(refInfos[refPos], len(reMapSeq) - subSize, subSize).second < matchingKmerCutOff) {
													continue;
												}
												//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
												const auto & ref = finalSeqsOnClusters[refPos];
												alignerObj.alignCacheGlobal(ref.seqBase_,reMapSeq);
												alignerObj.profilePrimerAlignment(ref.seqBase_,reMapSeq);
												if (alignerObj.comp_.distances_.query_.coverage_ >= allowableErrors.distances_.query_.coverage_ &&
														allowableErrors.passErrorProfile(alignerObj.comp_)) {
													//if the current score is a certain amount within the best score and doesn't have large gaps then put it
													//into the vector for possible match
													if (alignerObj.parts_.score_> bestScore) {
														bestScore = alignerObj.parts_.score_;
														bestRefs.clear();
														bestRefs.emplace_back(refPos);
													}else if(alignerObj.parts_.score_ == bestScore){
														bestRefs.emplace_back(refPos);
													}
												}
											}

											if (bestRefs.size() == 1) {
												//if only matching one, then place with that one
												uniqueSeqs[finalSeqsOnClusters[bestRefs.front()].seqBase_.name_].emplace(reMapSeq.name_);
												alignerObj.alignCacheGlobal(finalSeqsOnClusters[bestRefs.front()].seqBase_,reMapSeq);
												for(const auto seqAlignPos : iter::range(alignerObj.alignObjectA_.seqBase_.seq_.size())){
													if(alignerObj.alignObjectA_.seqBase_.seq_[seqAlignPos] != '-' &&
														 alignerObj.alignObjectB_.seqBase_.seq_[seqAlignPos] != '-'){
														baseCoverrage[finalSeqsOnClusters[bestRefs.front()].seqBase_.name_][alignerObj.getSeqPosForAlnAPos(seqAlignPos)][reMapSeq.name_] = 1;
													}
												}
											} else if (bestRefs.size() == 0) {
												//if no mappable candidates found, put into un-mappable
												++unmappable;
											} else {
												VecStr matchingRefNames;
												for(const auto & refPos : bestRefs) {
													matchingRefNames.emplace_back(finalSeqsOnClusters[refPos].seqBase_.name_);
												}
												std::vector<uint64_t> matchingRefs;
												for (const auto & refPos : bestRefs) {
													const auto & ref = finalSeqsOnClusters[refPos];
													//get the current snp info for the current ref to the other matching refs
													auto seqSnpPosBases = refVariationInfo[refPos].getVariantSnpLociMap(matchingRefNames, 0);
													alignerObj.alignCacheGlobal(ref.seqBase_, reMapSeq);
													alignerObj.profilePrimerAlignment(ref.seqBase_, reMapSeq);
													//determine the snps for this current read to this ref
													std::unordered_map<uint32_t, char> currentSnps;
													for (const auto & m : alignerObj.comp_.distances_.mismatches_) {
														currentSnps[m.second.refBasePos] = m.second.seqBase;
													}
	//												printOutMapContents(currentSnps, "\t", std::cout);
	//												std::cout << std::endl;
	//												for(const auto & seqSnps : seqSnpPosBases){
	//													std::cout << seqSnps.first << std::endl;
	//													std::cout << "\t" << njh::conToStr(seqSnps.second) << std::endl;
	//												}
													bool pass = true;
													uint32_t queryStart = alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_first_not_of("-"));
													uint32_t queryStop = alignerObj.getSeqPosForAlnAPos(alignerObj.alignObjectB_.seqBase_.seq_.find_last_not_of("-"));
	//												std::cout << "queryStart: " << queryStart << std::endl;
	//												std::cout << "queryStop:  "  << queryStop << std::endl;
													//iterate over segregating snps locations
													for(const auto & loc : seqSnpPosBases) {
														if(loc.first < queryStart || loc.first > queryStop){
															continue;
														}
														//determine if current read has a snp at location
														auto search = currentSnps.find(loc.first);
														if(search != currentSnps.end()) {
															//if it does have a snp here, check to see if it is a known variant, if it is the read doesn't pass
															if(njh::in(search->second, loc.second)) {
																pass = false;
																break;
															}
														}
													}
													if(pass) {
														auto seqSnpPos = refVariationInfo[refPos].getVariantSnpLoci(matchingRefNames,0);
	//													std::cout << "snp loci" << std::endl;
	//													std::cout << "\t" << njh::conToStr(seqSnpPos) << std::endl;

														for(const auto & loc : seqSnpPos) {
															if(loc < queryStart || loc > queryStop){
																continue;
															}
															//determine if current read has a snp at location
															if(currentSnps.find(loc) == currentSnps.end()) {
																//if there is a gap here, don't believe it either
																if(alignerObj.alignObjectB_.seqBase_.seq_[alignerObj.getAlignPosForSeqAPos(loc)] == '-') {
																	pass = false;
																	break;
																}
															}
														}
													}
	//												std::cout << njh::colorBool(pass) << std::endl;
	//												if(!pass){
	//													std::stringstream ss;
	//													std::cerr << __PRETTY_FUNCTION__ << "  " << __LINE__ << " error, no pass: name: " << reMapSeq.name_ <<  " in "  << region.uid_  << std::endl;
	//													//throw std::runtime_error{ss.str()};
	//													exit(1);
	//												}
													//if all the segregating check out, put this ref in the matchingRefs
													if(pass) {
														matchingRefs.emplace_back(refPos);
													}
												}
												if(matchingRefs.size() == 1) {
													//if only one matching ref, add to it
													uniqueSeqs[finalSeqsOnClusters[matchingRefs.front()].seqBase_.name_].emplace(reMapSeq.name_);
													alignerObj.alignCacheGlobal(finalSeqsOnClusters[matchingRefs.front()].seqBase_,reMapSeq);
													for(const auto seqAlignPos : iter::range(alignerObj.alignObjectA_.seqBase_.seq_.size())){
														if(alignerObj.alignObjectA_.seqBase_.seq_[seqAlignPos] != '-' &&
															 alignerObj.alignObjectB_.seqBase_.seq_[seqAlignPos] != '-'){
															baseCoverrage[finalSeqsOnClusters[matchingRefs.front()].seqBase_.name_][alignerObj.getSeqPosForAlnAPos(seqAlignPos)][reMapSeq.name_] = 1;
														}
													}
												} else if(matchingRefs.size() == 0) {
													//if none of the ref work out, put it in indeterminate
													++indeterminate;
												} else {
													//if there is more than one matching ref, put it in ties
													njh::sort(matchingRefs);
													std::string refMixName = njh::conToStr(matchingRefs, "_");
													multiMappingSeqs[refMixName].emplace_back(reMapSeq);
													++multiMapping;
												}
											}

										};

										auto checkFileIncreaseCount = [&increaseMapCounts,&matchingKmerCutOff](const bfs::path & seqFnp){
											if(bfs::exists(seqFnp)){
												seqInfo reMappingSeq;
												SeqInput remapReader(SeqIOOptions::genFastqIn(seqFnp));
												remapReader.openIn();
												while(remapReader.readNextRead(reMappingSeq)){
													increaseMapCounts(reMappingSeq, matchingKmerCutOff);
												}
											}
										};
										//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
										checkFileIncreaseCount(r1Fnp);
										checkFileIncreaseCount(r2Fnp);
										checkFileIncreaseCount(singlesFnp);
										//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
										std::unordered_map<std::string, double> uniqueMapFrac;
										std::unordered_map<std::string, double> finalTotals;
										double totalUniques = 0;
										for(const auto & uni : uniqueSeqs){
											totalUniques += uni.second.size();
											finalTotals[uni.first] = uni.second.size();
										}
										for(const auto & clus : finalSeqsOnClusters){
											if(njh::in(clus.seqBase_.name_, uniqueSeqs)){
												uniqueMapFrac[clus.seqBase_.name_] = uniqueSeqs[clus.seqBase_.name_].size()/totalUniques;
											}else{
												/**@todo might want to use a different amount or just get rid of the read, worried about low coverage*/
												uniqueMapFrac[clus.seqBase_.name_] = 1/totalUniques;
											}
										}
										//filter in case mate was already found
										std::unordered_map<std::string, std::vector<seqInfo>> filteredMultiMappingSeqs;
										for(const auto & multi : multiMappingSeqs){
											for(const auto & seqName : multi.second){
												bool found = false;
												for(const auto & uni : uniqueSeqs){
													if(njh::in(seqName.name_, uni.second)){
														found = true;
														break;
													}
												}
												if(!found){
													filteredMultiMappingSeqs[multi.first].emplace_back(seqName);
												}
											}
										}
										for(const auto & multi : filteredMultiMappingSeqs){
											auto refPoss = vecStrToVecNum<uint32_t>(tokenizeString(multi.first, "_"));
											double totalFrac = 0;
											for(const auto & pos : refPoss){
												totalFrac += uniqueMapFrac[finalSeqsOnClusters[pos].seqBase_.name_];
											}
											auto names = readVec::getNames(multi.second);
											std::set<std::string> namesSet(names.begin(), names.end());

											for(const auto & pos : refPoss){
							//											std::cout << "pos: " << pos << std::endl;
							//											std::cout << "multi.second.size(): " << multi.second.size() << std::endl;
							//											std::cout << "multi.second.size() * (uniqueMapFrac[finalSeqsOnClusters[pos].seqBase_.name_]/totalFrac): " <<  multi.second.size() * (uniqueMapFrac[finalSeqsOnClusters[pos].seqBase_.name_]/totalFrac) << std::endl;
												finalTotals[finalSeqsOnClusters[pos].seqBase_.name_] += namesSet.size() * (uniqueMapFrac[finalSeqsOnClusters[pos].seqBase_.name_]/totalFrac);
												for(const auto & reMapSeq : multi.second){
													alignerObj.alignCacheGlobal(finalSeqsOnClusters[pos].seqBase_, reMapSeq);
													for(const auto seqAlignPos : iter::range(alignerObj.alignObjectA_.seqBase_.seq_.size())){
														if(alignerObj.alignObjectA_.seqBase_.seq_[seqAlignPos] != '-' &&
															 alignerObj.alignObjectB_.seqBase_.seq_[seqAlignPos] != '-'){
															baseCoverrage[finalSeqsOnClusters[pos].seqBase_.name_][alignerObj.getSeqPosForAlnAPos(seqAlignPos)].emplace(reMapSeq.name_, uniqueMapFrac[finalSeqsOnClusters[pos].seqBase_.name_]/totalFrac);
														}
													}
												}
											}
										}

										for(auto & clus : finalSeqsOnClusters){
											clus.seqBase_.cnt_ = std::round(finalTotals[clus.seqBase_.name_]);
										}

										readVec::allSetFractionByTotalCount(finalSeqsOnClusters);
										/**@todo consider putting a fraction cut off to include final
										 *
										 */
										readVecSorter::sort(finalSeqsOnClusters);



										std::unordered_map<std::string, uint32_t> nameKey;
										uint32_t namePos = 0;
										for(const auto & clus : finalSeqsOnClusters){
											nameKey[clus.seqBase_.name_] = namePos;
											++namePos;
										}
										std::unordered_map<uint32_t, std::vector<double>> baseCoveragesFinal;
										for(const auto & b : baseCoverrage){
											for(const auto & pos : b.second){
												double baseTotal = 0;
												for(const auto & readWeight : pos.second){
													baseTotal += readWeight.second;
												}
												baseCoveragesFinal[nameKey[b.first]].emplace_back(baseTotal);
											}
										}


										for(const auto & bCov : baseCoveragesFinal){
											double meanCoverage = 0;
											if(bCov.second.size() > 60){
												//if(finalSeqsOnClusters[bCov.first].seqBase_.seq_.size() > 60){
												meanCoverage = vectorMean(getSubVector(bCov.second, 25, finalSeqsOnClusters[bCov.first].seqBase_.seq_.size() - 50));
											}else{
												meanCoverage = vectorMean(bCov.second);
											}
	//										if(std::string::npos != finalSeqsOnClusters[bCov.first].seqBase_.name_.find(region.uid_)){
	//											uint32_t pos = 0;
	//											for(const auto & b : bCov.second){
	//												std::cout  << pos << "\t"<< b << std::endl;
	//												++pos;
	//											}
	//											std::cout << std::endl;
	//											std::cout << vectorMean(getSubVector(bCov.second, 25, finalSeqsOnClusters[bCov.first].seqBase_.seq_.size() - 50)) << std::endl;
	//											std::cout << vectorMean(bCov.second) << std::endl;
	//											uint32_t below = 0;
	//											for(auto pos : iter::range<uint32_t>(25, finalSeqsOnClusters[bCov.first].seqBase_.seq_.size() - 50)){
	//												if(bCov.second[pos] < baseCoverageCutOff){
	//													++below;
	//												}
	//											}
	//											std::cout << getPercentageString(below, finalSeqsOnClusters[bCov.first].seqBase_.seq_.size() - 50 - 25) << std::endl;
	//											std::cout << std::endl;
	//										}
											MetaDataInName nameMeta(finalSeqsOnClusters[bCov.first].seqBase_.name_);
											nameMeta.addMeta("meanBaseCoverage", meanCoverage, true);
											nameMeta.resetMetaInName(finalSeqsOnClusters[bCov.first].seqBase_.name_);
										}
										Json::Value mapCounts;
										mapCounts["unmappable"]    = njh::json::toJson(unmappable);
										mapCounts["indeterminate"] = njh::json::toJson(indeterminate);
										mapCounts["totalUniques"]  = njh::json::toJson(totalUniques);
										mapCounts["multiMapping"] = njh::json::toJson(multiMapping);
										mapCounts["totalReadIn"]  =  njh::json::toJson(allReadInSamples.size());
										OutputStream mapCountsOut(OutOptions(njh::files::make_path(finalInfoRegionDir,"mapCounts.json")));
										mapCountsOut << mapCounts << std::endl;
										if(writeOutBaseCoveragePerSeq){
											OutOptions baseCoverageOutOpts(njh::files::make_path(finalInfoRegionDir,"baseCoverage.tab.txt"));
											OutputStream baseCoverageOut(baseCoverageOutOpts);
											baseCoverageOut << "name\tposition\tcoverage" << std::endl;
											for(const auto & b : baseCoverrage){
												for(const auto & pos : b.second){
													double baseTotal = 0;
													for(const auto & readWeight : pos.second){
														baseTotal += readWeight.second;
													}
													baseCoverageOut << finalSeqsOnClusters[nameKey[b.first]].seqBase_.name_
															<< "\t" << pos.first
															<< "\t" << baseTotal << std::endl;
												}
											}
										}
									}


									if(filterOnBaseCoverage){
										uint32_t offCount = 0;
										for( auto & seq : finalSeqsOnClusters){
											MetaDataInName nameMeta(seq.seqBase_.name_);
											if(filterOnBaseCoverage && nameMeta.getMeta<double>(coverageField) < baseCoverageCutOff){
												seq.seqBase_.on_ = false;
											}
											if(!seq.seqBase_.on_){
												++offCount ;
											}
										}
										if(finalSeqsOnClusters.size() == offCount){
											std::stringstream ss;
											ss << __PRETTY_FUNCTION__ << " error, no final seqs above mean base coverage cut off of " << baseCoverageCutOff << "\n";
											throw std::runtime_error{ss.str()};
										}else{
											auto finalSeqBellowMeanBaseCovOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalInfoRegionDir,"seqBelowCoverageCutOff.fasta"));
											SeqOutput cutOffWriter(finalSeqBellowMeanBaseCovOpts);
											std::vector<uint32_t> toBeRemoved;
											for(const auto & seqPos : iter::range(finalSeqsOnClusters.size())){
												const auto & seq = finalSeqsOnClusters[seqPos];
												if(!seq.seqBase_.on_){
													cutOffWriter.openWrite(seq);
													toBeRemoved.emplace_back(seqPos);
												}
											}
											for(auto seqPos : iter::reversed(toBeRemoved)){
												finalSeqsOnClusters.erase(finalSeqsOnClusters.begin() + seqPos);
											}
										}
									}
									double totalMeanCoverage = 0;
									double totalReadCnt = 0;
									for(auto & clus : finalSeqsOnClusters){
										MetaDataInName metaName(clus.seqBase_.name_);
										totalMeanCoverage += metaName.getMeta<double>(coverageField);
										totalReadCnt+=clus.seqBase_.cnt_;
									}
									OutOptions outFinalInfoOpts(njh::files::make_path(finalInfoRegionDir,"finalOutputInfo.tab.txt"));
									OutputStream outFinalInfo(outFinalInfoOpts);
									outFinalInfo << "ReadName\treadCount\treadFraction\tlength\t"<< coverageField << "\tfractionByCoverage";
									outFinalInfo << "\n";
									uint32_t count = 0;
									for(auto & clus : finalSeqsOnClusters){
										MetaDataInName metaName(clus.seqBase_.name_);
										metaName.addMeta("length", len(clus), true);
										clus.seqBase_.name_ = region.uid_ + "." + leftPadNumStr<uint32_t>(count, finalSeqsOnClusters.size()) + metaName.createMetaName();
										clus.updateName();
										outFinalInfo << clus.seqBase_.name_
												<< "\t" << clus.seqBase_.cnt_
												<< "\t" << clus.seqBase_.cnt_/totalReadCnt
												<< "\t" << len(clus)
												<< "\t" << metaName.getMeta<double>(coverageField)
												<< "\t" << metaName.getMeta<double>(coverageField)/totalMeanCoverage
												<< "\n";
										++count;
									}
									auto finalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalInfoRegionDir,"finalOutput.fasta"));
									SeqOutput::write(finalSeqsOnClusters, finalSeqOpts);
									//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
								}else{
									auto partialInfoRegionDir = njh::files::makeDir(partialDirectory,
											njh::files::MkdirPar(region.createUidFromCoordsStrand()));
									SeqInput originalSeqReader(outputOpts);
									SeqOutput originalWriter(SeqIOOptions::genFastaOut(njh::files::make_path(partialInfoRegionDir, "output.fasta")));
									seqInfo seq;
									originalSeqReader.openIn();
									originalWriter.openOut();

									while(originalSeqReader.readNextRead(seq)){
										MetaDataInName seqMeta(seq.name_);
										seqMeta.addMeta("regionUID", region.uid_);
										seqMeta.addMeta("regionCoords", region.createUidFromCoordsStrand());
										seqMeta.resetMetaInName(seq.name_);
										originalWriter.write(seq);
									}
								}
							}
						}
					} catch (std::exception & e) {
						try {
							auto finalInfoRegionDir = njh::files::make_path(finalDirectory,
														region.createUidFromCoordsStrand());
							if(bfs::exists(finalInfoRegionDir)){
								njh::files::rmDirForce(finalInfoRegionDir);
							}
						} catch (std::exception & e2) {
							logValue["exception-" + region.createUidFromCoordsStrand()] = e.what() + std::string("\n")+ e2.what();
						}
						logValue["exception-" + region.createUidFromCoordsStrand()] = e.what();
					}

					if(writeOutLogs){
						jLog.addToLog(region.createUidFromCoordsStrand(), logValue);
					}
				}
			};


	if (masterPars.pFinderPars_.numThreads <= 1) {
		extractPathway(masterPars);
	} else {
		std::vector<std::thread> threads;
		for (uint32_t t = 0; t < masterPars.pFinderPars_.numThreads; ++t) {
			threads.emplace_back(std::thread(extractPathway, masterPars));
		}
		njh::concurrent::joinAllThreads(threads);
	}

	jLog.writeLog();


	//make a table of coi counts
	//make a table of bed6 plus the coi at that location and total read count

	//make a table that is basically a concatenation of the finalOutputInfo files
	//mkae a table of call rates
	std::unordered_map<uint32_t, std::uint32_t> coiCounts;
	OutputStream finalTableOut(OutOptions(njh::files::make_path(finalDirectory, "allFinalTable.tab.txt")));
	OutputStream coiTableOut(OutOptions(njh::files::make_path(finalDirectory, "coiCounts.tab.txt")));
	auto allFinalSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(finalDirectory, "allFinal.fasta"));
	auto allPartialSeqOpts = SeqIOOptions::genFastaOut(njh::files::make_path(partialDirectory, "allPartial.fasta"));

	SeqOutput allFinalWriter(allFinalSeqOpts);
	SeqOutput allPartialWriter(allPartialSeqOpts);

	seqInfo seq;
	allFinalWriter.openOut();
	allPartialWriter.openOut();
	//collect meta fields
	std::set<std::string> metaFields;
	std::vector<uint32_t> readTotals;
	for(const auto & reg : regions){
		auto finalOutputSeqFnp = njh::files::make_path(finalDirectory, reg.createUidFromCoordsStrand(), "finalOutput.fasta");
		if(bfs::exists(finalOutputSeqFnp)){
			SeqInput finalReader(SeqIOOptions::genFastaIn(finalOutputSeqFnp));
			finalReader.openIn();
			while(finalReader.readNextRead(seq)){
				allFinalWriter.write(seq);
			}
		}
		auto partialOutputSeqFnp = njh::files::make_path(partialDirectory, reg.createUidFromCoordsStrand(), "output.fasta");
		if(bfs::exists(partialOutputSeqFnp)){
			SeqInput partialReader(SeqIOOptions::genFastaIn(partialOutputSeqFnp));
			partialReader.openIn();
			while(partialReader.readNextRead(seq)){
				allPartialWriter.write(seq);
			}
		}
		auto finalOutputInfoFnp = njh::files::make_path(finalDirectory, reg.createUidFromCoordsStrand(), "finalOutputInfo.tab.txt");
		if(bfs::exists(finalOutputInfoFnp)){
			table finalInfoTab(finalOutputInfoFnp, "\t", true);
			uint32_t readTotal = 0;
			for(const auto & row : finalInfoTab.content_){
				readTotal += njh::lexical_cast<uint32_t>(row[finalInfoTab.getColPos("readCount")]);
				MetaDataInName metaInReadName(row[finalInfoTab.getColPos("ReadName")]);
				for(const auto & m : metaInReadName.meta_){
					if("length" != m.first && coverageField != m.first){
						metaFields.insert(m.first);
					}
				}
			}
			readTotals.emplace_back(readTotal);
		}
	}
	allFinalWriter.closeOut();
	allPartialWriter.closeOut();

	//auto readTotalMedian = vectorMedianRef(readTotals);


	finalTableOut << "ReadName\treadCount\treadFraction\tlength\t" << coverageField << "\tfractionByCoverage\t" << njh::conToStr(metaFields, "\t") << "\n";

	for(const auto & reg : regions){
		auto regUID = reg.createUidFromCoordsStrand();
		auto finalOutputInfoFnp = njh::files::make_path(finalDirectory, regUID, "finalOutputInfo.tab.txt");
		if(bfs::exists(finalOutputInfoFnp)){
			table finalInfoTab(finalOutputInfoFnp, "\t", true);
			//add to coi count
			++coiCounts[finalInfoTab.nRow()];
			//add to final table
			uint32_t readTotal = 0;
			for(const auto & row : finalInfoTab.content_){
				finalTableOut << njh::conToStr(row, "\t");
				readTotal += njh::lexical_cast<uint32_t>(row[finalInfoTab.getColPos("readCount")]);
				//add meta
				MetaDataInName metaInReadName(row[finalInfoTab.getColPos("ReadName")]);
				for(const auto & field : metaFields){
					finalTableOut << "\t" << (metaInReadName.containsMeta(field) ? metaInReadName.getMeta(field) : "NA");
				}
				finalTableOut << std::endl;
			}
			//output to called calledTableOut
			regInfosByUID[regUID]->uniqHaps_ = finalInfoTab.nRow();
			regInfosByUID[regUID]->infoCalled_ = true;
		}else{
			//output to bed coi file
			//output to called calledTableOut
			regInfosByUID[regUID]->uniqHaps_ = 0;
			regInfosByUID[regUID]->infoCalled_ = false;
		}
	}

	brInvestor.writeBasicInfoWithHapRes(regInfos, sampName, OutOptions(njh::files::make_path(finalDirectory, "basicInfoPerRegion.tab.txt")));


	table coiCountsTab(coiCounts, VecStr{"coi", "count"});
	coiCountsTab.sortTable("coi", false);
	coiCountsTab.outPutContents(coiTableOut, "\t");

	if(!keepTemporaryFiles){
		auto perBaseCovFnp = njh::files::make_path(finalDirectory, "perBaseCoveragePerRegion.bed");
		if(bfs::exists(perBaseCovFnp)){
			bfs::remove(perBaseCovFnp);
		}
		for(const auto & region : regions){
			auto regionDir = njh::files::make_path(finalDirectory,region.createUidFromCoordsStrand());
			if(bfs::exists(regionDir)){
				njh::files::rmDirForce(regionDir);
			}
		}
		for(const auto & region : regions){
			auto regionDir = njh::files::make_path(partialDirectory,region.createUidFromCoordsStrand());
			if(bfs::exists(regionDir)){
				njh::files::rmDirForce(regionDir);
			}
		}
		for(const auto & region : regions){
			auto regionDir = njh::files::make_path(masterPars.outputDir,region.createUidFromCoordsStrand());
			if(bfs::exists(regionDir)){
				njh::files::rmDirForce(regionDir);
			}
		}
	}
	return 0;
}



}  // namespace njhseq
