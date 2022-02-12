/*
 * WeaverRunner_PathFinderAcrossRegions.cpp
 *
 *  Created on: Jun 1, 2017
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

#include <njhcpp/concurrency/LockableJsonLog.hpp>

#include <TwoBit.h>

#include <njhseq/BamToolsUtils.h>
#include <njhseq/objects/BioDataObject.h>
#include <njhseq/GenomeUtils.h>
#include <njhseq/objects/seqContainers.h>

#include "PathWeaver/objects/Meta/CountryMetaData.hpp"
#include "PathWeaver/seqToolsUtils/HaplotypeLocator.hpp"
#include "PathWeaver/PathFinding.h"
#include "PathWeaver/objects/bam/RegionInvestigatorInBam.hpp"


namespace njhseq {

void setExtractPathWaysDefault(HaploPathFinder::ExtractParams & pars, seqSetUp & setUp, bool needGenome = true ){
	//optimization opts
	pars.pFinderPars_.setOptimizationParameters(setUp);

	//pair stitching parameters
	pars.pFinderPars_.setStitchParameters(setUp);

	//group info
	pars.pFinderPars_.setAddingGroupInfoOpts(setUp);

	//possible and major haplotypes
	pars.pFinderPars_.setPossibleHapsOpts(setUp);

	//pre-process
	//tandem repeat
	pars.pFinderPars_.setTandemRepeatHandlingOptions(setUp);
	//duplicate seqs handling
	pars.pFinderPars_.setDuplicatedSeqHandlingOpts(setUp);

	//quality trimming and filtering
	pars.pFinderPars_.setQualityTrimAndFiltOpts(setUp);

	//graph options
	//splitting
	pars.pFinderPars_.setFurtherSplittingOpts(setUp);
	//processing nodes, onebase indel nodes, remove headless and tailless nodes
	pars.pFinderPars_.setNodeProcessingOpts(setUp);

	//genome options
	pars.setGenomeOpts(setUp, needGenome);
	//run modes
	pars.pFinderPars_.setRunningOpts(setUp);
	//outliers
	pars.pFinderPars_.setOutlierRemovalOptions(setUp);

	//post processing of collapsing contigs;
	pars.pFinderPars_.setPostProcessContigHandlingOpts(setUp);

	//adding meta to the final seqs
	setUp.setOption(pars.metaDataFnp, "--metaDataFnp", "Name of the meta data fnp");

	//// specific for here
	//length cut offs
	pars.pFinderPars_.setFinalLengthCutOffs(setUp);
	//trimming
	pars.pFinderPars_.setPostProcessTrimmingOpts(setUp);
	//adding on of headless or tailess nodes
	setUp.setOption(pars.pFinderPars_.addSeqOfSingleHeadAndTailSeqs_, "--addHeadAndTailSeqs", "Add Head And Tail sequences when creating output nodes");

	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapInfo_ = gapScoringParameters(5,1,0,0,0,0);

}


int WeaverRunner::SeqsExtractPathaways(const njh::progutils::CmdArgs & inputCommands){
	HaploPathFinder::ExtractParams pars;
	std::string sampName = "";
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	double fracCutOff = 0;

	setUp.setOption(fracCutOff, "--fracCutOff", "Fraction Cut Off");
	setUp.setOption(sampName,   "--sampName", "Samp Name", true);
	pars.pFinderPars_.reOrientPairsForStitching = false;
	setExtractPathWaysDefault(pars, setUp, false);
	setUp.processRefFilename(pars.pFinderPars_.trimToInputSeqs);
	bfs::path fastq1Fnp = "";
	bfs::path fastq2Fnp = "";
	bfs::path fastqFnp = "";
	bool setFastq1 = setUp.setOption(fastq1Fnp, "--fastq1", "Fastq first mate File");
	setUp.setOption(fastq2Fnp, "--fastq2", "Fastq second mate File", setFastq1);
	bool revCompMate = false;
	setUp.setOption(revCompMate, "--revCompMate", "Reverse Complement Sequences in mate file");
	setUp.setOption(fastqFnp, "--fastq", "Fastq File", !setFastq1);

	if (setFastq1) {
		setUp.pars_.ioOptions_.firstName_ = fastq1Fnp;
		setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQPAIRED;
		setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQPAIRED;
	} else {
		setUp.pars_.ioOptions_.firstName_ = fastqFnp;
		setUp.pars_.ioOptions_.inFormat_ = SeqIOOptions::inFormats::FASTQ;
		setUp.pars_.ioOptions_.outFormat_ = SeqIOOptions::outFormats::FASTQ;
	}
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	pars.pFinderPars_.verbose = setUp.pars_.verbose_;
	//pars.pFinderPars_.debug = setUp.pars_.debug_;
	std::unique_ptr<MultipleGroupMetaData> meta;
	if ("" != pars.metaDataFnp) {
		meta = std::make_unique<MultipleGroupMetaData>(pars.metaDataFnp,
				std::set<std::string> { sampName });
	} else {
		meta = std::make_unique<MultipleGroupMetaData>(table { VecStr { "sample",
				"holder" } }, std::set<std::string> { "empty" });
	}



	auto pairedIn = SeqIOOptions::genPairedIn(fastq1Fnp, fastq2Fnp);
	pairedIn.revComplMate_ = revCompMate;
	auto singleIn = SeqIOOptions::genFastqIn(fastqFnp);
	BamExtractor::ExtractedFilesOpts inputOpts;
	inputOpts.inPairs_ = pairedIn;
	inputOpts.inUnpaired_ = singleIn;
	if("" != setUp.pars_.refIoOptions_.firstName_.string()){
		pars.pFinderPars_.inputSeqs = SeqInput::getSeqVec<seqInfo>(setUp.pars_.refIoOptions_);
	}
	Json::Value currentLog;
	auto PWRes = PathFinderFromSeqsDev(inputOpts, setUp.pars_.directoryName_, sampName, pars.pFinderPars_, meta);
	currentLog[sampName] = PWRes.log_;
	OutOptions logOutOpts(njh::files::make_path(setUp.pars_.directoryName_, "extractionLog.json"));
	auto logFilePtr = logOutOpts.openFile();
	(*logFilePtr) << currentLog << std::endl;

	return 0;
}

int WeaverRunner::ExtractPathWaysReadsFallingInMultipleRegions(const njh::progutils::CmdArgs & inputCommands){
	bool development = false;

	HaploPathFinder::ExtractParams pars;

	std::string sampName = "";
	bool optimizeAgainAndRecruitBeforeFinal = false;
	uint32_t maxIteration = 20;
	uint32_t finalLenCutOff = 0;
	uint32_t remappRecruitmentContigLenCutOff = 80;
	bool expandEndRegions = false;
	uint32_t expandRegionBy = 25;

	bfs::path bedFile = "";
	bool useBowtie2 = false;

	bool addRefsForSingleRegions = false;
	bfs::path refSeqExtractsDir = "";
	BioCmdsUtils::LastZPars lzPars;
	lzPars.identity = 80;
	lzPars.coverage = 90;
	bool keepBestOnly = false;
	bool keepIntermediatefiles = false;

	BamExtractor::extractReadsFromBamToSameOrientationContigsPars bamContigsExtractPars;

	bool redetermineReadCounts = false;
 	comparison mapAllowableErrors;
 	mapAllowableErrors.distances_.query_.coverage_ = 0.25;
 	mapAllowableErrors.lqMismatches_ = 2;
 	mapAllowableErrors.hqMismatches_ = 1;
 	mapAllowableErrors.oneBaseIndel_ = 0.1;
 	mapAllowableErrors.twoBaseIndel_ = 0.1;
 	mapAllowableErrors.largeBaseIndel_ = 0;


 	bool trimFinalEnds = false;
 	uint32_t trimFinalEndsBy = 0;


 	uint32_t finalHeadlessTaillessCutOff = 0;

	BamRegionInvestigator::BamRegionInvestigatorPars brInvestPars;
	BamRegionInvestigator::BamCountSpecficRegionsPars spanningReadsPar;

	uint32_t bwaRemappingSeedLength = 19;

	seqSetUp setUp(inputCommands);





	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(bedFile, "--bed", "Bed file, first entry is used", true);

	if(bfs::exists(bedFile)){
		auto beds = getBeds(bedFile);
		if(beds.size() == 1){
			maxIteration = 0;
			pars.pFinderPars_.autoDetermineKCuts = true;
			pars.pFinderPars_.trimToInputSeqs = true;
			if(beds.front()->length() >=100){
				pars.pFinderPars_.headlessTailessLenCutOff = 100;
			}
		} else {
			bool allAbove100 = true;
			for(const auto & b : beds){
				if(b->length() > 100){
					allAbove100  = false;
					break;
				}
			}
			if(allAbove100){
				pars.pFinderPars_.headlessTailessLenCutOff = 100;
			}
			pars.pFinderPars_.autoDetermineKCutsOnTotalBaseCount = true;
		}

	}else{
		setUp.failed_ = true;
		setUp.addWarning(njh::pasteAsStr(bedFile, " doesn't exist"));
	}
	setUp.setOption(maxIteration, "--maxIteration", "Maximum Iterations to be done");

	{
		//when only one region and no iterations, the default should be to remove markedDups and to keep duplicated sequence
		if(bfs::exists(bedFile)){
			auto beds = getBeds(bedFile);
			if(beds.size() == 1 && 0 == maxIteration){
				pars.bamExtractPars_.keepMarkedDuplicate_ = false;
				pars.pFinderPars_.removeDuplicatedSequences_ = false;
				pars.bamExtractPars_.removeImproperPairs_ = true;
				pars.bamExtractPars_.keepImproperMateUnmapped_ = true;
			}else{
				pars.bamExtractPars_.keepMarkedDuplicate_ = true;
				pars.pFinderPars_.removeDuplicatedSequences_ = true;
			}
		}
	}
	pars.pFinderPars_.reOrientPairsForStitching = true;
	setExtractPathWaysDefault(pars, setUp);
	finalHeadlessTaillessCutOff = pars.pFinderPars_.headlessTailessLenCutOff;
	setUp.setOption(finalHeadlessTaillessCutOff, "--finalHeadlessTaillessCutOff", "Final Cut for headless and tailless nodes to use for the -finalPass");
	setUp.setOption(bwaRemappingSeedLength, "--bwaRemappingSeedLength", "bwa mem Re-mapping Seed Length");



	pars.setBamExtractOpts(setUp);

	//pars.unmappedPars_.maxAlnSize_ = pars.bamExtractPars_.minAlnMapSize_;
	pars.unmappedPars_.maxAlnSize_ = 0;
	setUp.setOption(pars.unmappedPars_.maxAlnSize_, "--maxAlnSizeForAddRecruit", "max Aln Size For additional recruitment");
	setUp.setOption(pars.unmappedPars_.minQuerySize_, "--minQuerySize", "minimum query size for additional recruitment");
	setUp.setOption(pars.unmappedPars_.filterOffLowEntropyShortAlnsRecruitsCutOff_, "--filterOffLowEntropyShortAlnsRecruitsCutOff", "Filter Off Low Entropy Short Alns Recruits Cut Off");

	uint32_t expandRegionsToAvoidBy = 500;
	setUp.setOption(expandRegionsToAvoidBy, "--expandRegionsToAvoidby", "Expand Regions To Avoid By");


	setUp.setOption(brInvestPars.mapQualityCutOffForMultiMap_, "--mapQualityCutOffForMultiMapForCovInfo", "map Quality Cut Off For Coverage");
	setUp.setOption(brInvestPars.mapQualityCutOff_,            "--mapQualityCutOffForCov", "map Quality Cut Off For Coverage");



	setUp.setOption(development, "--development", "development mode, saves a lot more files and information");
	setUp.setOption(remappRecruitmentContigLenCutOff, "--remappRecruitmentContigLenCutOff", "Remapping Recruitment Contig Len Cut Off, this is to prevent using super small contigs for mapping recruitment");
	setUp.setOption(optimizeAgainAndRecruitBeforeFinal,   "--optimizeAgainAndRecruitBeforeFinal",   "Optimize Again And Recruit Before Final");
	setUp.setOption(pars.pFinderPars_.optAfterFirstRecruit, "--optAfterFirstRecruit", "Optimize parameters again after initial recruitment");
	setUp.setOption(pars.pFinderPars_.optimizeOnIteration,  "--optimizeOnIteration", "Optimize parameters on each iteration ");

	setUp.setOption(trimFinalEnds,   "--trimFinalEnds",   "Trim Final Ends");
	setUp.setOption(trimFinalEndsBy, "--trimFinalEndsBy", "Trim Final Ends By this much ");

	setUp.setOption(useBowtie2, "--useBowtie2", "Use bowtie2 for the remapping");

	bamContigsExtractPars.throwAwayUnmmpaedMates = false;
	setUp.setOption(bamContigsExtractPars.throwAwayUnmmpaedMates, "--throwAwayRemappedUnMappedMate", "Throw Away Re-mapped Unmapped Mate");
	setUp.setOption(bamContigsExtractPars.centerClipCutOff, "--centerClipCutOffForMappedToContigs", "Center Clip Cut Off For Mapped To Contigs to prevent un-specific mapping to contigs");
	setUp.setOption(bamContigsExtractPars.forSoftClipFilterDistanceToEdges, "--forSoftClipFilterDistanceToEdgesForMappedToContigs", "For Soft Clip Filter Distance To Edges For Mapped To Contigs");

	setUp.setOption(redetermineReadCounts, "--redetermineReadCounts", "Re-determine Read Counts by remapping final ");

	setUp.setOption(lzPars.identity, "--identity", "Identity to use for when searching with lastz");
	setUp.setOption(lzPars.coverage, "--coverage", "Coverage to use for when searching with lastz");
	setUp.setOption(keepBestOnly, "--keepBestOnly", "Keep best hits only");

	setUp.setOption(keepIntermediatefiles, "--keepIntermediatefiles", "Keep the intermediate files from ther iterative remapping process");

	setUp.setOption(addRefsForSingleRegions, "--addRefsForSingleRegions", "Add Refs seqs For Single Regions");
	setUp.setOption(refSeqExtractsDir, "--refSeqExtractsDir", "A directory where reference sequences have already been extracted to");



	setUp.setOption(pars.pFinderPars_.removePossibleOutliersWithMuscle, "--removePossibleOutliersWithMuscle", "Remove Possible Outliers based off of the alignment to the input region seqs");
	setUp.setOption(pars.pFinderPars_.outliersCutOff, "--outliersCutOff", "outliers Cut Off value for when removing possible outliers is turned on");
	setUp.setOption(pars.pFinderPars_.mtPars.spanningCutOff, "--spanningCutOff", "Spanning Cut Off");
	setUp.setOption(pars.pFinderPars_.mtPars.baseCutOff, "--gapCutOff", "gap Cut Off");
	setUp.setOption(pars.pFinderPars_.mtPars.hardGapCutOff, "--hardGapCutOff", "Hard Gap Cut Off");
	setUp.setOption(pars.pFinderPars_.mtPars.streakLenCutOff, "--streakLenCutOff", "streak Len Cut Off");

	setUp.setOption(pars.trimSeqBedFnp, "--trimSeqBedFnp", "Trim to these seqs instead");


	bool doNotExapndRegions = false;
	setUp.setOption(doNotExapndRegions, "--doNotExapndRegions", "Do Not Expand Ends");
	expandEndRegions = !doNotExapndRegions;
	setUp.setOption(expandRegionBy,   "--expandRegionBy", "Expand Region By this length for the inital pull down");


	bool thrownAwayDiscordant = false;
	setUp.setOption(thrownAwayDiscordant, "--thrownAwayDiscordant", "Throw away discordant Pairs when mapping unmapped sequences");
	bool keepDiscordant = !thrownAwayDiscordant;


	setUp.setOption(finalLenCutOff, "--finalLenCutOff", "Final Length Cut Off to be included in the muscle trimming");
	if(0 != pars.pFinderPars_.lenCutOff && 0 == finalLenCutOff){
		finalLenCutOff = pars.pFinderPars_.lenCutOff;
	}
	setUp.processReadInNames( { "--bam" }, true);
	sampName = getPossibleSampleNameFromFnp(setUp.pars_.ioOptions_.firstName_);
	setUp.setOption(sampName, "--sampName", "Sample Name");
	setUp.processDirectoryOutputName(true);

	seqInfo realignToSeq("realign");
	setUp.processSeq(realignToSeq, "--realignToSeq", "realign raw extract To input Seq");
	std::string extraBwaArgsForReAlign = "";

	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	pars.pFinderPars_.verbose = setUp.pars_.verbose_;
	//pars.pFinderPars_.debug = setUp.pars_.debug_;
	spanningReadsPar.numThreads = pars.pFinderPars_.numThreads;
	spanningReadsPar.countDuplicates = pars.bamExtractPars_.keepMarkedDuplicate_;


	if(setUp.pars_.verbose_){
		std::cout << "Kmer Lengths      : " << njh::conToStr(pars.pFinderPars_.kmerLengths) << std::endl;
		std::cout << "Occurence Cut Offs: " << njh::conToStr(pars.pFinderPars_.kmerKOcurrenceCutOffs) << std::endl;
		std::cout << "Short Tip Number  : " << njh::conToStr(pars.pFinderPars_.shortTipNumbers) << std::endl;
	}


	if("" != realignToSeq.seq_){
		BamExtractor bExtractor(setUp.pars_.verbose_);
		auto realignmentDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"realignment"});
		OutOptions extractRegion(njh::files::make_path(setUp.pars_.directoryName_));

		auto realignmentDirGenome = njh::files::makeDir(realignmentDir, njh::files::MkdirPar{"genome"});
		auto realignmentBamDir = njh::files::makeDir(realignmentDir, njh::files::MkdirPar{"bam"});

		auto genomeFnp = njh::files::make_path(realignmentDirGenome, "realign.fasta");
		SeqOutput::write(std::vector<seqInfo>{realignToSeq}, SeqIOOptions::genFastaOut(genomeFnp));
		BioCmdsUtils bRunner(setUp.pars_.verbose_);
		auto indexRunLog = bRunner.RunBwaIndex(genomeFnp);
		BioCmdsUtils::checkRunOutThrow(indexRunLog, __PRETTY_FUNCTION__);
		auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);

		for(const auto & region : inputRegions){
			OutOptions seqOutOpts(njh::files::make_path(realignmentDir, "rawExtract"));
			seqOutOpts.append_ = true;
			bExtractor.writeExtractReadsFromBamRegion(setUp.pars_.ioOptions_.firstName_, region, pars.bamExtractPars_.percInRegion_, seqOutOpts);
		}
		bfs::path inputSingles = njh::files::make_path(njh::files::make_path(realignmentDir, "rawExtract.fastq"));
		bfs::path inputPairedFirstMates = njh::files::make_path(njh::files::make_path(realignmentDir, "rawExtract_R1.fastq"));;
		bfs::path inputPairedSecondMates = njh::files::make_path(njh::files::make_path(realignmentDir, "rawExtract_R2.fastq"));;
		njh::files::checkExistenceThrow(genomeFnp,__PRETTY_FUNCTION__);
		bfs::path outputFnp = njh::files::make_path(realignmentBamDir, "realigned.sorted.bam");

		bfs::path outputFnpBai = outputFnp.string() + ".bai";
		bfs::path singlesSortedBam = njh::files::make_path(realignmentBamDir, inputSingles.filename().string() + ".sorted.bam");
		bfs::path pairedSortedBam = njh::files::make_path(realignmentBamDir,  inputPairedFirstMates.filename().string() + ".sorted.bam");

		bool useSambamba = false;
		std::stringstream singlesCmd;
		auto singlesBwaLogFnp = bfs::path(singlesSortedBam.string() + ".bwa.log");
		singlesCmd << "bwa mem  -M -t " << pars.pFinderPars_.numThreads
				<< " -R " << R"("@RG\tID:)" << bfs::basename(inputSingles) << "" << R"(\tSM:)"
				<< sampName << R"(")"
				<< " "   << extraBwaArgsForReAlign
				<< " "   << genomeFnp
				<< " "   << inputSingles
				<< " 2> " << singlesBwaLogFnp;
		if(useSambamba){
			singlesCmd << " | sambamba view -S /dev/stdin -o /dev/stdout -f bam | sambamba sort -t " << pars.pFinderPars_.numThreads << " -o " << singlesSortedBam << " /dev/stdin";
		}else{
			singlesCmd << " | samtools sort -@ " << pars.pFinderPars_.numThreads << " -o " << singlesSortedBam;
		}

		std::stringstream pairedCmd;
		auto pairedBwaLogFnp = bfs::path(pairedSortedBam.string() + ".bwa.log");
		pairedCmd << "bwa mem  -M -t " << pars.pFinderPars_.numThreads
				<< " -R " << R"("@RG\tID:)" << bfs::basename(inputPairedFirstMates) << "" << R"(\tSM:)"
				<< sampName << R"(")"
				<< " " << extraBwaArgsForReAlign
				<< " " << genomeFnp
				<< " " << inputPairedFirstMates
				<< " " << inputPairedSecondMates
				<< " 2> " << pairedBwaLogFnp;
		if (useSambamba) {
			pairedCmd << " | sambamba view -S /dev/stdin -o /dev/stdout -f bam | sambamba sort -t " << pars.pFinderPars_.numThreads << " -o " << pairedSortedBam << " /dev/stdin";
		} else {
			pairedCmd << " | samtools sort -@ " << pars.pFinderPars_.numThreads << " -o " << pairedSortedBam;
		}


		std::stringstream bamtoolsMergeAndIndexCmd;
		if (useSambamba) {
			bamtoolsMergeAndIndexCmd << "sambamba merge " << outputFnp << " " << pairedSortedBam;
			if(bfs::exists(inputSingles)){
				bamtoolsMergeAndIndexCmd <<  " " << singlesSortedBam;
			}
		}else{
			bamtoolsMergeAndIndexCmd << "bamtools merge " << " -in " << pairedSortedBam;

			if(bfs::exists(inputSingles)){
				bamtoolsMergeAndIndexCmd << " -in " <<  singlesSortedBam;
			}
			bamtoolsMergeAndIndexCmd << " -out " << outputFnp
					<< " && samtools index " << outputFnp;
		}


		if(true){
			bfs::path logFnp = njh::files::make_path(realignmentBamDir, "alignTrimoOutputs_" + sampName + "_" + njh::getCurrentDate() + "_log.json");
			logFnp = njh::files::findNonexitantFile(logFnp);
			OutOptions logOpts(logFnp);
			std::ofstream logFile;
			logOpts.openFile(logFile);
			std::unordered_map<std::string, njh::sys::RunOutput> runOutputs;
			if(bfs::exists(inputSingles)){
				auto singlesRunOutput = njh::sys::run({singlesCmd.str()});
				BioCmdsUtils::checkRunOutThrow(singlesRunOutput, __PRETTY_FUNCTION__);
				runOutputs["bwa-singles"] = singlesRunOutput;
			}
			if(bfs::exists(inputPairedFirstMates)){
				auto pairedRunOutput = njh::sys::run({pairedCmd.str()});
				BioCmdsUtils::checkRunOutThrow(pairedRunOutput, __PRETTY_FUNCTION__);
				runOutputs["bwa-paired"] = pairedRunOutput;
			}
			if(bfs::exists(singlesSortedBam) && bfs::exists(pairedSortedBam)){
				auto bamtoolsMergeAndIndexRunOutput = njh::sys::run({bamtoolsMergeAndIndexCmd.str()});
				BioCmdsUtils::checkRunOutThrow(bamtoolsMergeAndIndexRunOutput, __PRETTY_FUNCTION__);
				runOutputs["bamtools-merge-index"] = bamtoolsMergeAndIndexRunOutput;
			} else {
				if(bfs::exists(pairedSortedBam)){
					bfs::rename(pairedSortedBam, outputFnp);
					if(useSambamba){
						bfs::rename(pairedSortedBam.string() + ".bai", outputFnp.string() + ".bai");
					} else {
						std::stringstream ss;
						ss << "samtools index " << outputFnp;
						auto indexRunOutput = njh::sys::run({ss.str()});
						BioCmdsUtils::checkRunOutThrow(indexRunOutput, __PRETTY_FUNCTION__);
						runOutputs["index"] = indexRunOutput;
					}
				}else if(bfs::exists(singlesSortedBam) ){
					bfs::rename(singlesSortedBam, outputFnp);
					if(useSambamba){
						bfs::rename(singlesSortedBam.string() + ".bai", outputFnp.string() + ".bai");
					} else {
						std::stringstream ss;
						ss << "samtools index " << outputFnp;
						auto indexRunOutput = njh::sys::run({ss.str()});
						BioCmdsUtils::checkRunOutThrow(indexRunOutput, __PRETTY_FUNCTION__);
						runOutputs["index"] = indexRunOutput;
					}
				}
			}
			logFile << njh::json::toJson(runOutputs) << std::endl;
			if(true){
				if(bfs::exists(pairedSortedBam)){
					bfs::remove(pairedSortedBam);
				}
				if(bfs::exists(singlesSortedBam)){
					bfs::remove(singlesSortedBam);
				}
			}
		}
		setUp.pars_.ioOptions_.firstName_ = outputFnp;
		pars.genomeFnp = genomeFnp;
		pars.genomeDir = realignmentDirGenome;
		pars.primaryGenome = "realign";
		{
			OutputStream bedOut(njh::files::make_path(realignmentDir, "region.bed"));
			trimAtFirstWhitespace(realignToSeq.name_);
			bedOut << realignToSeq.name_
					<< "\t" << 0
					<< "\t" << len(realignToSeq)
					<< "\t" << realignToSeq.name_
					<< "\t" << len(realignToSeq)
					<< "\t" << '+' << std::endl;
		}
		bedFile = njh::files::make_path(realignmentDir, "region.bed");
	}

	std::unique_ptr<MultipleGroupMetaData> meta;
	if ("" != pars.metaDataFnp) {
		meta = std::make_unique<MultipleGroupMetaData>(pars.metaDataFnp,
				std::set<std::string> { sampName });
	} else {
		meta = std::make_unique<MultipleGroupMetaData>(table { VecStr { "sample",
				"holder" } }, std::set<std::string> { "empty" });
	}

	bfs::path logDir = njh::files::make_path(setUp.pars_.directoryName_, "logs");
	bfs::path extractionStatsDir = njh::files::make_path(setUp.pars_.directoryName_, "extractionStats");

	njh::files::makeDir(njh::files::MkdirPar{logDir});
	njh::files::makeDir(njh::files::MkdirPar{extractionStatsDir});

	Muscler muscleRunner;

	BamExtractor bExtractor(setUp.pars_.verbose_);
	bExtractor.debug_ = setUp.pars_.debug_;

	auto inputRegions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);

	OutOptions originalInputBedOpts(njh::files::make_path(setUp.pars_.directoryName_, "inputRegions.bed"));
	auto originalInputBed = originalInputBedOpts.openFile();
	for(const auto & region : inputRegions){
		(*originalInputBed ) << region.genBedRecordCore().toDelimStr() << std::endl;
	}

	HaploPathFinder hFinder(pars);
	auto regions = inputRegions;

	//if(pars.pFinderPars_.trimEnds || expandEndRegions){
	if(expandEndRegions){
		/**@todo reasses the need for expanding region */
		uint32_t expandRegionLength = pars.pFinderPars_.kmerLengths.front() + 5;
		if (0 != expandRegionBy) {
			expandRegionLength = expandRegionBy;
		}
		for(auto & region : regions){
			region.start_ = region.start_ > expandRegionLength ? region.start_ - expandRegionLength : 0;
			region.end_ +=  expandRegionLength;
		}
		OutputStream usedInputBedOut(njh::files::make_path(setUp.pars_.directoryName_, "usedExpandedRegions.bed"));
		for(const auto & region : regions){
			usedInputBedOut << region.genBedRecordCore().toDelimStr() << std::endl;
		}
	}
	uint32_t minInputLen = std::numeric_limits<uint32_t>::max();

	if (inputRegions.size() == 1 && addRefsForSingleRegions) {
		auto inputRefSeqs  = njh::files::make_path(refSeqExtractsDir, inputRegions.front().createUidFromCoordsStrand(), "allRefs.fasta");
		if(bfs::exists(inputRefSeqs)){
			auto refSeqsOpts = SeqIOOptions::genFastaIn(inputRefSeqs);
			if(!pars.pFinderPars_.rawInputSeqs){
				refSeqsOpts.lowerCaseBases_ = "upper";
			}
			SeqInput reader(refSeqsOpts);
			pars.pFinderPars_.inputSeqs = reader.readAllReads<seqInfo>();
		}else{
			auto refAlignsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("refAlignmentsToDtermineRegions"));
			std::unique_ptr<MultiGenomeMapper> gMapper;
			//set up genome mapper;
			gMapper = std::make_unique<MultiGenomeMapper>(pars.genomeDir, pars.primaryGenome);
			//set threads;
			if(pars.pFinderPars_.numThreads >= 4){
				gMapper->pars_.numThreads_ = 2;
				pars.pFinderPars_.numThreads = pars.pFinderPars_.numThreads/2;
			}
			//set up selected genomes
			gMapper->setSelectedGenomes(pars.selectedGenomes);
			gMapper->init();
			MultiGenomeMapper::getRefSeqsWithPrimaryGenomePars extractPars;
			extractPars.lzPars = lzPars;
			extractPars.keepBestOnly = keepBestOnly;
			uint64_t maxLen = inputRegions.front().getLen();
			aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0));
			pars.pFinderPars_.inputSeqs = gMapper->getRefSeqsWithPrimaryGenome(inputRegions.front(), refAlignsDir, extractPars, alignerObj);
			if(!pars.pFinderPars_.rawInputSeqs){
				readVec::handelLowerCaseBases(pars.pFinderPars_.inputSeqs, "upper");
			}
		}
	} else {
		for (const auto & region : inputRegions) {
			TwoBit::TwoBitFile tReader(
					hFinder.gMapper_->genomes_.at(hFinder.pars_.primaryGenome)->fnpTwoBit_);
			pars.pFinderPars_.inputSeqs.emplace_back(region.extractSeq(tReader));
		}
		if(!pars.pFinderPars_.rawInputSeqs){
			readVec::handelLowerCaseBases(pars.pFinderPars_.inputSeqs, "upper");
		}
	}

	if("" != pars.trimSeqBedFnp ){
		auto trimInputRegions = gatherRegions(pars.trimSeqBedFnp.string(), "", setUp.pars_.verbose_);
		for (const auto & region : trimInputRegions) {
			TwoBit::TwoBitFile tReader(
					hFinder.gMapper_->genomes_.at(hFinder.pars_.primaryGenome)->fnpTwoBit_);
			pars.pFinderPars_.trimSeqs.emplace_back(region.extractSeq(tReader));
		}
		if(!pars.pFinderPars_.rawInputSeqs){
			readVec::handelLowerCaseBases(pars.pFinderPars_.trimSeqs, "upper");
		}
		auto trimInputRegionsSeqsOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "trimInputRegions.fasta"));
		SeqOutput::write(pars.pFinderPars_.trimSeqs, trimInputRegionsSeqsOpts);
		OutOptions trimInputBedOpts(njh::files::make_path(setUp.pars_.directoryName_, "trimInputRegions.bed"));
		auto trimInputBed = trimInputBedOpts.openFile();
		for(const auto & region : trimInputRegions){
			(*trimInputBed ) << region.genBedRecordCore().toDelimStr() << std::endl;
		}
	}
	for(const auto & seq : pars.pFinderPars_.trimSeqs.empty() ? pars.pFinderPars_.inputSeqs : pars.pFinderPars_.trimSeqs){
		if (len(seq) < minInputLen) {
			minInputLen = len(seq);
		}
	}
	if(finalLenCutOff == 0){
		finalLenCutOff = pars.pFinderPars_.kmerLengths.front();
		if(minInputLen > finalLenCutOff){
			finalLenCutOff = minInputLen - pars.pFinderPars_.kmerLengths.front();
		}
	}

	auto inputRegionsSeqsOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, "inputRegions.fasta"));
	SeqOutput::write(pars.pFinderPars_.inputSeqs, inputRegionsSeqsOpts);

	auto extractedOpts = setUp.pars_.ioOptions_;
	extractedOpts.out_.outFilename_ = njh::files::make_path(
			setUp.pars_.directoryName_, "extracted");

	auto regionExtracted = bExtractor.extractReadsWtihCrossRegionMapping(extractedOpts, regions, pars.bamExtractPars_);
	{
		OutOptions extractedReadsStatsFileOpts(njh::files::make_path(extractionStatsDir, "extractedReadsStats.tab.txt"));
		OutputStream extractedReadsStatsFileOut(extractedReadsStatsFileOpts);
		regionExtracted.log(extractedReadsStatsFileOut, setUp.pars_.ioOptions_.firstName_);
	}
	bool iterate = true;
	bool noSplittingTailed = false;
	bool noSplittingEnd = false;
	if(pars.pFinderPars_.splitToRecruit_){
		noSplittingTailed = !pars.pFinderPars_.splitTailed_;
		if(noSplittingTailed){
			pars.pFinderPars_.splitTailed_ = true;
		}
		noSplittingEnd = !pars.pFinderPars_.splitEnds_;
		if(noSplittingEnd){
			pars.pFinderPars_.splitEnds_ = true;
		}
	}
	bool noAddHeadTail = false;
	if(pars.pFinderPars_.addHeadTailSeqsToRecruit){
		noAddHeadTail = !pars.pFinderPars_.addSeqOfSingleHeadAndTailSeqs_;
		if(noAddHeadTail){
			pars.pFinderPars_.addSeqOfSingleHeadAndTailSeqs_ = true;
		}
	}

	//get coverage and spanning read info
	brInvestPars.numThreads_ = pars.pFinderPars_.numThreads;
	brInvestPars.countDups_ = pars.bamExtractPars_.keepMarkedDuplicate_;

	BamRegionInvestigator brInvestor(brInvestPars);
	auto regInfos = brInvestor.getCoverageAndFullSpanningReads(
			setUp.pars_.ioOptions_.firstName_, inputRegions, spanningReadsPar, pars.genomeFnp);

	brInvestor.writeBasicInfo(regInfos, sampName, OutOptions(njh::files::make_path(extractionStatsDir, "perBaseCoveragePerRegion.bed")));



	Json::Value fullLog;
	std::string initialSampleName = sampName + "-" + (0 == maxIteration ? "0" : leftPadNumStr<uint32_t>(0, maxIteration) ) ;


	std::vector<uint32_t> occurrencCutOffAfterInitialIter;


	if((pars.pFinderPars_.kmerKOcurrenceCutOffs.size()> 1 || pars.pFinderPars_.forceDetermineKCuts)  && (pars.pFinderPars_.autoDetermineKCuts || pars.pFinderPars_.forceDetermineKCuts)){
		std::vector<double> perBaseCoverages;
		for(const auto & regInfo : regInfos){
			perBaseCoverages.emplace_back(regInfo->getPerBaseCov());
		}
		auto perBaseCoverage = vectorMean(perBaseCoverages);
		if(perBaseCoverage > 100 || pars.pFinderPars_.forceDetermineKCuts){
			pars.pFinderPars_.kmerKOcurrenceCutOffs.clear();
			for(const auto & perc : pars.pFinderPars_.autoDetermineKCutsPercentages){
				occurrencCutOffAfterInitialIter.emplace_back(round(perBaseCoverage * perc));
			}
			if(maxIteration > 0 && !pars.pFinderPars_.initialAutoDetermineKCutsPercentages.empty()){
				for(const auto & perc : pars.pFinderPars_.initialAutoDetermineKCutsPercentages){
					pars.pFinderPars_.kmerKOcurrenceCutOffs.emplace_back(round(perBaseCoverage * perc));
				}
			}else{
				pars.pFinderPars_.kmerKOcurrenceCutOffs = occurrencCutOffAfterInitialIter;
			}
		}
	}

	if(maxIteration > 0 && !pars.pFinderPars_.initialKmerKOcurrenceCutOffs.empty()){
		occurrencCutOffAfterInitialIter = pars.pFinderPars_.kmerKOcurrenceCutOffs;
		pars.pFinderPars_.kmerKOcurrenceCutOffs = pars.pFinderPars_.initialKmerKOcurrenceCutOffs;
	}

	auto autoDetermineKCutsPercentagesAfterInitialIter = pars.pFinderPars_.autoDetermineKCutsPercentages;
	if(maxIteration > 0 && !pars.pFinderPars_.initialAutoDetermineKCutsPercentages.empty()){
		pars.pFinderPars_.autoDetermineKCutsPercentages = pars.pFinderPars_.initialAutoDetermineKCutsPercentages;
	}

	std::vector<uint32_t> kmerLengthsAfterInitialIter;
	if(maxIteration > 0 && !pars.pFinderPars_.initialKmerLengths.empty()){
		kmerLengthsAfterInitialIter = pars.pFinderPars_.kmerLengths;
		pars.pFinderPars_.kmerLengths = pars.pFinderPars_.initialKmerLengths;
	}

	std::vector<uint32_t> shortTipNumbersInitialIter;
	if(maxIteration > 0 && !pars.pFinderPars_.initialShortTipNumbers.empty()){
		shortTipNumbersInitialIter = pars.pFinderPars_.shortTipNumbers;
		pars.pFinderPars_.shortTipNumbers = pars.pFinderPars_.initialShortTipNumbers;
	}


	auto initialPwRes = PathFinderFromSeqsDev(regionExtracted, setUp.pars_.directoryName_, initialSampleName, pars.pFinderPars_, meta);

	std::vector<bfs::path> iterationDirectories;
	iterationDirectories.emplace_back(njh::files::make_path(setUp.pars_.directoryName_, initialSampleName));

	if(!occurrencCutOffAfterInitialIter.empty()){
		pars.pFinderPars_.kmerKOcurrenceCutOffs = occurrencCutOffAfterInitialIter;
	}

	if(!kmerLengthsAfterInitialIter.empty()){
		pars.pFinderPars_.kmerLengths = kmerLengthsAfterInitialIter;
	}

	if(!shortTipNumbersInitialIter.empty()){
		pars.pFinderPars_.shortTipNumbers = shortTipNumbersInitialIter;
	}

	if(maxIteration > 0 && !pars.pFinderPars_.initialAutoDetermineKCutsPercentages.empty()){
		pars.pFinderPars_.autoDetermineKCutsPercentages = autoDetermineKCutsPercentagesAfterInitialIter;
	}

	fullLog["initial"] = initialPwRes.log_;

	uint32_t bestKmerCut = initialPwRes.log_["bestKmerOccurenceCutOff"].asInt();
	uint32_t bestKmerLength = initialPwRes.log_["bestKmerLength"].asInt();
	uint32_t bestShortTip = initialPwRes.log_["bestShortTipNumber"].asInt();



	if(0 == maxIteration){
		iterate = false;
	}


	uint32_t iterNumber = 1;
	std::string lastIterationSampName = initialSampleName;
	auto lastIterationPwRes = initialPwRes;
	if(pars.pFinderPars_.removePossibleOutliersKSim && pars.pFinderPars_.filterOffOutlierInputSeqs){
		//reoptimize if seqs for possible outliers were removed
		auto lastIteractionPaired = SeqIOOptions::genPairedIn(
		njh::files::make_path(
				setUp.pars_.directoryName_, lastIterationSampName, "filteredExtractedPairs_R1.fastq"),
		njh::files::make_path(
				setUp.pars_.directoryName_, lastIterationSampName, "filteredExtractedPairs_R2.fastq"));
		auto lastIteractionSingles = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "filteredSingles.fastq"));

		std::string currentIterationSampleName = njh::pasteAsStr(lastIterationSampName, "-reoptimized");
		BamExtractor::ExtractedFilesOpts iterationOpts;
		iterationOpts.inPairs_ = lastIteractionPaired;
		iterationOpts.inUnpaired_ = lastIteractionSingles;
		auto currentPwRes = PathFinderFromSeqsDev(iterationOpts, setUp.pars_.directoryName_, currentIterationSampleName, pars.pFinderPars_, meta, lastIterationPwRes);
		iterationDirectories.emplace_back(njh::files::make_path(setUp.pars_.directoryName_, currentIterationSampleName));
		bestKmerCut = currentPwRes.log_["bestKmerOccurenceCutOff"].asInt();
		bestKmerLength = currentPwRes.log_["bestKmerLength"].asInt();
		bestShortTip = currentPwRes.log_["bestShortTipNumber"].asInt();
		if(!keepIntermediatefiles){
			//clean up
			iterationOpts.removeAllInFiles();

			std::vector<bfs::path> oldInputFiles{
				"extractedPairsWithNoTandems_R1.fastq",
				"extractedPairsWithNoTandems_R2.fastq",
				"extractedPairsWithTandems_R1.fastq",
				"extractedPairsWithTandems_R2.fastq",
				"extractedPairs_R1.fastq",
				"extractedPairs_R2.fastq",
				"extractedSingles.fastq",
				"extractedSinglesWithNoTandems.fastq",
				"extractedSinglesWithTandems.fastq",
				"filteredExtractedPairs_R1.fastq",
				"filteredExtractedPairs_R2.fastq",
				"filteredSingles.fastq",
			  "filteredOff_extractedSingles.fastq",
				"filteredOff_extractedPairs_R1.fastq",
				"filteredOff_extractedPairs_R2.fastq",
			  "outlierFilteredOff_extractedPairs",
			  "outlierFilteredOff_extractedSingles"};

			for(const auto & fn : oldInputFiles){
				auto fnp = njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, fn);
				if(bfs::exists(fnp)){
					bfs::remove(fnp);
				}
			}
		}
		fullLog[currentIterationSampleName] = currentPwRes.log_;
		lastIterationSampName = currentIterationSampleName;
		lastIterationPwRes = currentPwRes;
	}

	//get unmapped reads
	njhseq::BamExtractor::ExtractedFilesOpts unmappedFiles;
	auto unmappedOpts = setUp.pars_.ioOptions_;
	unmappedOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, "unmapped-" + leftPadNumStr<uint32_t>(iterNumber, maxIteration));
	unmappedOpts.out_.outExtention_ = ".bam";
	if(iterate){
		ReadCheckerQualCheck qualChecker(pars.pFinderPars_.qualCheck_, pars.pFinderPars_.qualCheckCutOff_, false);

		auto regionsToAvoid = inputRegions;
		njh::for_each(regionsToAvoid, [&expandRegionsToAvoidBy](GenomicRegion & region){
			BedUtility::extendLeftRight(region, expandRegionsToAvoidBy, expandRegionsToAvoidBy);
		});
		unmappedFiles = bExtractor.writeUnMappedSeqsAndSmallAlnsWithFilters(unmappedOpts, qualChecker, pars.unmappedPars_, regionsToAvoid);
		OutOptions unmappedInitialPassOutOpts(njh::files::make_path(extractionStatsDir, "initialUnmappedPass.tab.txt"));
		OutputStream unmappedInitialPassOut(unmappedInitialPassOutOpts);
		unmappedFiles.log(unmappedInitialPassOut,setUp.pars_.ioOptions_.firstName_);
	}



	SeqIOOptions reextractedSeqsOpts;
	SeqIOOptions reextractedSeqsOptsSingles;
	BamExtractor::ExtractedFilesOpts rextractedSeqs;
	BamExtractor::ExtractedFilesOpts rextractedSeqsSingles;

//	uint32_t originalKmerCutOff = pars.pFinderPars_.kOccurenceCutOff;
//	uint32_t originalKmerLength = pars.pFinderPars_.klen;
//	uint32_t originalKmerCutOffOptStop = pars.pFinderPars_.optimizeKcutStop;
//	uint32_t originalKmerLenOptStop = pars.pFinderPars_.optimizeKLenStop;
	auto originalKmerCutOffs = pars.pFinderPars_.kmerKOcurrenceCutOffs;
	auto originalKmerLengths = pars.pFinderPars_.kmerLengths;
	pars.pFinderPars_.originalKmerLengths = pars.pFinderPars_.kmerLengths;
	auto originalShorttipNumbers = pars.pFinderPars_.shortTipNumbers;

	bool alreadyTriedToOptimizedAgain = false;

	uint32_t iterCountAfterReOptimization = 0;
	while(iterate){
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "iterate: " << njh::colorBool(iterate) << std::endl;
//		std::cout << "keepDiscordant: " << njh::colorBool(keepDiscordant) << std::endl;
//		std::cout << "iterNumber: " << iterNumber << std::endl;
		//need to add re-pairing
		pars.pFinderPars_.needToRePair_ = true;

		if (iterNumber > 1) {
			//collect the previously unmapped reads
			unmappedFiles = njhseq::BamExtractor::ExtractedFilesOpts();
			unmappedFiles.inPairsUnMapped_ = SeqIOOptions::genPairedIn(
					njh::files::make_path(setUp.pars_.directoryName_, "unmapped-" + leftPadNumStr<uint32_t>(iterNumber, maxIteration) + "_R1.fastq"),
					njh::files::make_path(setUp.pars_.directoryName_, "unmapped-" + leftPadNumStr<uint32_t>(iterNumber, maxIteration) + "_R2.fastq"));
			unmappedFiles.inUnpairedUnMapped_ = SeqIOOptions::genFastqIn(
					njh::files::make_path(setUp.pars_.directoryName_, "unmapped-" + leftPadNumStr<uint32_t>(iterNumber, maxIteration) + ".fastq")
					);
			//group the un-mapping
			std::vector<bfs::path> firstMates;
			std::vector<bfs::path> secondMates;
			std::vector<bfs::path> singlesFnps;
			if(rextractedSeqs.inPairsUnMapped_.inExists()){
				firstMates.emplace_back(rextractedSeqs.inPairsUnMapped_.firstName_);
				secondMates.emplace_back(rextractedSeqs.inPairsUnMapped_.secondName_);
			}

			if(rextractedSeqs.inDiscordant_.inExists() && !keepDiscordant){
				firstMates.emplace_back(rextractedSeqs.inDiscordant_.firstName_);
				secondMates.emplace_back(rextractedSeqs.inDiscordant_.secondName_);
			}

			if(rextractedSeqs.inUnpairedUnMapped_.inExists()){
				singlesFnps.emplace_back(rextractedSeqs.inUnpairedUnMapped_.firstName_);
			}

			if(rextractedSeqsSingles.inUnpairedUnMapped_.inExists()){
				singlesFnps.emplace_back(rextractedSeqsSingles.inUnpairedUnMapped_.firstName_);
			}
//			std::cout << iterNumber << std::endl;
//			std::cout << "firstMates: " << njh::conToStr(firstMates, ", ") << std::endl;
//			std::cout << "singlesFnps: " << njh::conToStr(singlesFnps, ", ") << std::endl;
//			std::cout << std::endl;
			if(firstMates.empty() && singlesFnps.empty()){
				//no reads unmapped, break out
				iterate = false;
				break;
			} else {
				if(!firstMates.empty()){
					OutOptions firstMateOuts(unmappedFiles.inPairsUnMapped_.firstName_);
					firstMateOuts.overWriteFile_ = true;
					concatenateFiles(firstMates, firstMateOuts);
					if(secondMates.empty()){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << " error on iteration " << iterNumber << " first mates not empty but second mates is " << "\n";
						throw std::runtime_error{ss.str()};
					}
					OutOptions secondMateOuts(unmappedFiles.inPairsUnMapped_.secondName_);
					secondMateOuts.overWriteFile_ = true;
					concatenateFiles(secondMates, secondMateOuts);
				}
				if(!singlesFnps.empty()){
					OutOptions singlesOuts(unmappedFiles.inUnpairedUnMapped_.firstName_);
					singlesOuts.overWriteFile_ = true;
					concatenateFiles(singlesFnps, singlesOuts);
				}
			}
		}


		// Index the last output files
		auto rawLastIterationOuputFnp = njh::files::make_path(
				setUp.pars_.directoryName_, lastIterationSampName,
				"output_aboveCutOff.fasta");
		if(!bfs::exists(rawLastIterationOuputFnp)){
			std::stringstream ss;
			ss << rawLastIterationOuputFnp << " doesn't exist, likely there were no reads above cut off" << std::endl;
			throw std::runtime_error{ss.str()};
		}
		auto renamedLastIterationOuputFnp = njh::files::make_path(
				setUp.pars_.directoryName_, lastIterationSampName,
				"renamed_output_aboveCutOff.fasta");
		//need to rename in order to get proper mapping
		auto lastIterationSeqOpts =	SeqIOOptions::genFastaInOut(rawLastIterationOuputFnp, renamedLastIterationOuputFnp);
		std::unordered_map<std::string, std::string> lastIterationNameKey;
		{
			seqInfo seq;
			SeqIO renamingSeqIO(lastIterationSeqOpts);
			renamingSeqIO.openIn();
			renamingSeqIO.openOut();
			uint32_t seqCount = 0;
			while(renamingSeqIO.readNextRead(seq)){
				//don't use small contigs for mapping recruitment
				if(len(seq) < remappRecruitmentContigLenCutOff){
					continue;
				}
				lastIterationNameKey[seq.name_] = estd::to_string(seqCount);
				seq.name_ = estd::to_string(seqCount);
				renamingSeqIO.write(seq);
				++seqCount;
			}
		}




		Json::Value bwaLog;
		Json::Value bowtie2Log;

		if (!useBowtie2) {
			std::stringstream bwaIndexCmd;
			bwaIndexCmd << "bwa index " << renamedLastIterationOuputFnp;
			auto bwaIndexOutput = njh::sys::run({bwaIndexCmd.str()});
			BioCmdsUtils::checkRunOutThrow(bwaIndexOutput, __PRETTY_FUNCTION__);
			bwaLog["bwa-index"] = bwaIndexOutput.toJson();
		} else {
			BioCmdsUtils bRunner;
			auto bowtie2Output = bRunner.RunBowtie2Index(renamedLastIterationOuputFnp);
			bowtie2Log["bowtie2-build"] = bowtie2Output.toJson();
		}

		//Gather singles to remap
		std::vector<bfs::path> singlesReMappInput;
		if (rextractedSeqs.inThrownAwayUnmappedMate_.inExists()) {
			singlesReMappInput.emplace_back(rextractedSeqs.inThrownAwayUnmappedMate_.firstName_);
		}
		if (unmappedFiles.inUnpairedUnMapped_.inExists()) {
			singlesReMappInput.emplace_back(unmappedFiles.inUnpairedUnMapped_.firstName_);
		}

		if(iterNumber == 1){
			//if the first iteration grab the original thrown away mates
			if (bfs::exists(regionExtracted.inThrownAwayUnmappedMate_.firstName_)) {
				singlesReMappInput.emplace_back(regionExtracted.inThrownAwayUnmappedMate_.firstName_);
			}
			if (bfs::exists(regionExtracted.inFilteredSingles_.firstName_)) {
				singlesReMappInput.emplace_back(regionExtracted.inFilteredSingles_.firstName_);
			}
		}
//		std::cout << iterNumber << std::endl;
//		std::cout << "singlesReMappInput: " << njh::conToStr(singlesReMappInput, ", ") << std::endl;
//		std::cout << std::endl;

		bfs::path sortedReMappedBamFnp = njh::files::make_path(setUp.pars_.directoryName_,
				"reMapped-paired-"   + leftPadNumStr<uint32_t>(iterNumber, maxIteration) + ".bam");
		bfs::path sortedReMappedBamFnpUnpaired = njh::files::make_path(setUp.pars_.directoryName_,
				"reMapped-singles-" + leftPadNumStr<uint32_t>(iterNumber, maxIteration) + ".bam");
		// re map the single end sequences
		if(!singlesReMappInput.empty()){
			bfs::path sinlgesToRemapFnp = njh::files::make_path(setUp.pars_.directoryName_, "remappingSingles-" + leftPadNumStr<uint32_t>(iterNumber, maxIteration) + ".fastq");
			OutOptions sinlgesToRemapOpts(sinlgesToRemapFnp);
			concatenateFiles(singlesReMappInput, sinlgesToRemapOpts);

			if(useBowtie2){
				std::stringstream bowtie2MapCmds;
				auto bowtie2Prefix = bfs::path(renamedLastIterationOuputFnp).replace_extension("");
				bowtie2MapCmds << "bowtie2 -p " << pars.pFinderPars_.numThreads << " "
						<< "-x " << bowtie2Prefix  << " "
						<< "-U " << sinlgesToRemapFnp
						<< " 2> " << njh::files::make_path(logDir, njh::pasteAsStr("remappingSingles_" , leftPadNumStr(iterNumber, maxIteration) , "_bowtie2.log.txt"))
						<< " | samtools view - -o " << sortedReMappedBamFnpUnpaired;
				auto bowtie2MappCmdOutput = njh::sys::run({bowtie2MapCmds.str()});
				BioCmdsUtils::checkRunOutThrow(bowtie2MappCmdOutput, __PRETTY_FUNCTION__);
				bowtie2Log["bowtie2-map-sinlges"] = bowtie2MappCmdOutput.toJson();
			}else{
				std::stringstream bwaMappCmd;
				//bwaMappCmd << "bwa mem -M  -t " << pars.pFinderPars_.numThreads << " "
				bwaMappCmd << "bwa mem -M -k " << bwaRemappingSeedLength << " -t " << pars.pFinderPars_.numThreads << " "
						<< renamedLastIterationOuputFnp  << " "
						<< sinlgesToRemapFnp
						<< " 2> " << njh::files::make_path(logDir, njh::pasteAsStr("remappingSingles_" , leftPadNumStr(iterNumber, maxIteration) , "_bwa.log.txt"))
						<< " | samtools view - -o " << sortedReMappedBamFnpUnpaired;
				auto bwaMappCmdOutput = njh::sys::run({bwaMappCmd.str()});
				BioCmdsUtils::checkRunOutThrow(bwaMappCmdOutput, __PRETTY_FUNCTION__);
				bwaLog["bwa-map-sinlges"] = bwaMappCmdOutput.toJson();
			}

			if(!keepIntermediatefiles){
				//clean up
				for(const auto & fnp : singlesReMappInput){
					bfs::remove(fnp);
				}
				bfs::remove(sinlgesToRemapFnp);
			}
		}

		// re map the unmapped paired sequences
		std::vector<SeqIOOptions> remappingPairs;
		if(unmappedFiles.inPairsUnMapped_.inExists()){
			remappingPairs.emplace_back(unmappedFiles.inPairsUnMapped_);
		}
		if(iterNumber == 1){
			//if the first iteration grab the original filtered off pairs
			if (bfs::exists(regionExtracted.inFilteredPairs_.firstName_)) {
				remappingPairs.emplace_back(regionExtracted.inFilteredPairs_);
			}
		}
		if(setUp.pars_.debug_){
			std::cout << "regionExtracted.inFilteredPairs_" << regionExtracted.inFilteredPairs_.firstName_ << std::endl;
			std::cout << "iteration: " << iterNumber << std::endl;
			std::cout << "remappingPairs.size(): " << remappingPairs.size() << std::endl;
			for (const auto & remapPair : remappingPairs) {
				std::cout << "\t" << remapPair.firstName_ << std::endl;
				std::cout << "\t" << remapPair.secondName_ << std::endl;
			}
		}

		if(!remappingPairs.empty()){
			SeqIOOptions currentRemappingPairs;
			if (remappingPairs.size() > 1) {
				bfs::path firstMatesToRemapFnp = njh::files::make_path(setUp.pars_.directoryName_, "remappingPairs-" + leftPadNumStr<uint32_t>(iterNumber, maxIteration) + "_R1.fastq");
				bfs::path secondMatesToRemapFnp = njh::files::make_path(setUp.pars_.directoryName_, "remappingPairs-" + leftPadNumStr<uint32_t>(iterNumber, maxIteration) + "_R2.fastq");
				std::vector<bfs::path> firstMates;
				std::vector<bfs::path> secondMates;
				for (const auto & remapPair : remappingPairs) {
					firstMates.emplace_back(remapPair.firstName_);
					secondMates.emplace_back(remapPair.secondName_);
				}
				OutOptions firstMatesToRemapOpts(firstMatesToRemapFnp);
				concatenateFiles(firstMates, firstMatesToRemapOpts);

				OutOptions secondMatesToRemapOpts(secondMatesToRemapFnp);
				concatenateFiles(secondMates, secondMatesToRemapOpts);
				currentRemappingPairs = SeqIOOptions::genPairedIn(firstMatesToRemapFnp, secondMatesToRemapFnp);
			} else if (1 == remappingPairs.size()) {
				currentRemappingPairs = remappingPairs.front();
			}
			if(useBowtie2){
				std::stringstream bowtie2MapCmds;
				auto bowtie2Prefix = bfs::path(renamedLastIterationOuputFnp).replace_extension("");
				bowtie2MapCmds << "bowtie2 -p " << pars.pFinderPars_.numThreads << " "
						<< "-x " << bowtie2Prefix  << " "
						<< "-1 " << currentRemappingPairs.firstName_ << " "
						<< "-2 " << currentRemappingPairs.secondName_
						<< " 2> " << njh::files::make_path(logDir, njh::pasteAsStr("remappingPairs_" , leftPadNumStr(iterNumber, maxIteration) , "_bowtie2.log.txt"))
						<< " | samtools view - -o " << sortedReMappedBamFnp;
				auto bowtie2MappCmdOutput = njh::sys::run({bowtie2MapCmds.str()});
				BioCmdsUtils::checkRunOutThrow(bowtie2MappCmdOutput, __PRETTY_FUNCTION__);
				bowtie2Log["bowtie2-map-paired"] = bowtie2MappCmdOutput.toJson();
			}else{
				std::stringstream bwaMappCmd;
				//bwaMappCmd << "bwa mem -M  -t " << pars.pFinderPars_.numThreads << " "
				bwaMappCmd << "bwa mem -M -k " << bwaRemappingSeedLength << " -t " << pars.pFinderPars_.numThreads << " "
						<< renamedLastIterationOuputFnp  << " "
						<< currentRemappingPairs.firstName_ << " "
						<< currentRemappingPairs.secondName_
						<< " 2> " << njh::files::make_path(logDir, njh::pasteAsStr( "remappingPairs_", leftPadNumStr(iterNumber, maxIteration) , "_bwa.log.txt"))
						<< " | samtools view - -o " << sortedReMappedBamFnp;
				auto bwaMappCmdOutput = njh::sys::run({bwaMappCmd.str()});
				BioCmdsUtils::checkRunOutThrow(bwaMappCmdOutput, __PRETTY_FUNCTION__);
				bwaLog["bwa-map-paired"] = bwaMappCmdOutput.toJson();
			}
			if(!keepIntermediatefiles){
				//clean up
				if(bfs::exists(unmappedFiles.inPairsUnMapped_.firstName_)){
					bfs::remove(unmappedFiles.inPairsUnMapped_.firstName_);
					bfs::remove(unmappedFiles.inPairsUnMapped_.secondName_);
				}
				if(currentRemappingPairs.inExists()){
					bfs::remove(currentRemappingPairs.firstName_);
					bfs::remove(currentRemappingPairs.secondName_);
				}
			}
		}
		if(!keepIntermediatefiles){
			//clean up
			if(bfs::exists(unmappedOpts.out_.outName())){
				bfs::remove(unmappedOpts.out_.outName());
			}
		}
		if(useBowtie2){
			OutOptions bowtie2LogFileOpts(njh::files::make_path(logDir, "bowtie2-log-" + leftPadNumStr(iterNumber, maxIteration)  + ".json"));
			OutputStream bowtie2LogFile(bowtie2LogFileOpts);
			bowtie2LogFile << bowtie2Log << std::endl;
		} else	{
			OutOptions bwaLogFileOpts(njh::files::make_path(logDir, "bwa-log-" + leftPadNumStr(iterNumber, maxIteration)  + ".json"));
			OutputStream bwaLogFile(bwaLogFileOpts);
			bwaLogFile << bwaLog << std::endl;
		}
		std::string currentIterationSampleName = sampName + "-"
				+ leftPadNumStr(iterNumber, maxIteration) ;
		//extract all seqs



		//merge with the last iteration seqs
		//singles
		auto lastIteractionSingles = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "filteredSingles.fastq"));
		auto lastIteractionSinglesMarkedDups = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "filteredOffDups_extractedSingles.fastq"));

		auto nextInteractionSinglesInput = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, "singles_for_" + currentIterationSampleName + ".fastq"));

		std::vector<bfs::path> singlesInputs;

		if(!keepIntermediatefiles){
			//clean up
			rextractedSeqs.removeAllInFiles();
			rextractedSeqsSingles.removeAllInFiles();
		}
		//singles from mapping unmapped singles
		if (bfs::exists(sortedReMappedBamFnpUnpaired)) {
			reextractedSeqsOptsSingles.firstName_ = sortedReMappedBamFnpUnpaired;
			reextractedSeqsOptsSingles.out_.outFilename_ = njh::files::make_path(
					setUp.pars_.directoryName_, "reExtracted_fromSingles_" + leftPadNumStr(iterNumber, maxIteration));
			rextractedSeqsSingles = bExtractor.extractReadsFromBamToSameOrientationContigs(reextractedSeqsOptsSingles, bamContigsExtractPars);
			if (rextractedSeqsSingles.inUnpaired_.inExists()) {
				singlesInputs.emplace_back(rextractedSeqsSingles.inUnpaired_.firstName_);
			}
			OutOptions outSinglesReMappedFileOpts(njh::files::make_path(extractionStatsDir, leftPadNumStr(iterNumber, maxIteration)  + "_unmappedSinglesExtractedReadStats.tab.txt"));
			OutputStream outSinglesReMappedFile(outSinglesReMappedFileOpts);
			rextractedSeqsSingles.log(outSinglesReMappedFile,sortedReMappedBamFnpUnpaired);
			if(!keepIntermediatefiles){
				//clean up
				bfs::remove(sortedReMappedBamFnpUnpaired);
				auto bamWithMappedSingles = njh::files::make_path(setUp.pars_.directoryName_, "mapped_reExtracted_fromSingles_" + leftPadNumStr(iterNumber, maxIteration) + ".bam");
				bfs::remove(bamWithMappedSingles);
			}
		}



		//paired
		std::vector<bfs::path> firstMatesInput;
		std::vector<bfs::path> secondMatesInput;
		auto lastIteractionPaired = SeqIOOptions::genPairedIn(
				njh::files::make_path(
				    setUp.pars_.directoryName_, lastIterationSampName, "filteredExtractedPairs_R1.fastq"),
				njh::files::make_path(
						setUp.pars_.directoryName_, lastIterationSampName, "filteredExtractedPairs_R2.fastq"));

		auto nextInteractionPairedInput = SeqIOOptions::genPairedIn(njh::files::make_path(setUp.pars_.directoryName_, "pairs_for_" + currentIterationSampleName + "_R1.fastq"),
				njh::files::make_path(setUp.pars_.directoryName_, "pairs_for_" + currentIterationSampleName + "_R2.fastq"));

		if(bfs::exists(sortedReMappedBamFnp)){
			reextractedSeqsOpts.firstName_ = sortedReMappedBamFnp;
			reextractedSeqsOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, "reExtracted_fromPairs_" + leftPadNumStr(iterNumber, maxIteration));
			rextractedSeqs = bExtractor.extractReadsFromBamToSameOrientationContigs(reextractedSeqsOpts, bamContigsExtractPars);

			if (rextractedSeqs.inPairs_.inExists()) {
				firstMatesInput.emplace_back(rextractedSeqs.inPairs_.firstName_);
				secondMatesInput.emplace_back(rextractedSeqs.inPairs_.secondName_);
			}

			if(rextractedSeqs.inDiscordant_.inExists() && keepDiscordant){
				/**@todo make it so it's just chrom discordants we keep*/
				firstMatesInput.emplace_back(rextractedSeqs.inDiscordant_.firstName_);
				secondMatesInput.emplace_back(rextractedSeqs.inDiscordant_.secondName_);
			}
			if (rextractedSeqs.inPairsMateUnmapped_.inExists()) {
				firstMatesInput.emplace_back(rextractedSeqs.inPairsMateUnmapped_.firstName_);
				secondMatesInput.emplace_back(rextractedSeqs.inPairsMateUnmapped_.secondName_);
			}

			if (rextractedSeqs.inUnpaired_.inExists()) {
				singlesInputs.emplace_back(rextractedSeqs.inUnpaired_.firstName_);
			}
			OutOptions outPairsReMappedFileOpts(njh::files::make_path(extractionStatsDir, leftPadNumStr(iterNumber, maxIteration)  + "_unmappedPairsExtractedReadStats.tab.txt"));
			OutputStream outPairsReMappedFile(outPairsReMappedFileOpts);
			rextractedSeqs.log(outPairsReMappedFile, sortedReMappedBamFnp);
			if(!keepIntermediatefiles){
				//clean up
				bfs::remove(sortedReMappedBamFnp);
				auto bamWithMappedPairs = njh::files::make_path(setUp.pars_.directoryName_, "mapped_reExtracted_fromPairs_" + leftPadNumStr(iterNumber, maxIteration) + ".bam");
				bfs::remove(bamWithMappedPairs);
			}
		}
//		std::cout << iterNumber << std::endl;
//		std::cout << "firstMatesInput: " << njh::conToStr(firstMatesInput, ", ") << std::endl;
//		std::cout << "singlesInputs: " << njh::conToStr(singlesInputs, ", ") << std::endl;
//		std::cout << std::endl;
		if(pars.pFinderPars_.filterOffOutlierInputSeqs && pars.pFinderPars_.removePossibleOutliersKSim && iterNumber <= 1){

		}else{
			if(!pars.pFinderPars_.optimizeOnIteration ){
				if(pars.pFinderPars_.optAfterFirstRecruit && iterNumber <=1){

				}else{
					pars.pFinderPars_.kmerLengths = {bestKmerLength};
					pars.pFinderPars_.kmerKOcurrenceCutOffs = {bestKmerCut};
					pars.pFinderPars_.shortTipNumbers = {bestShortTip};
				}
			}
//			pars.pFinderPars_.kOccurenceCutOff = bestKmerCut;
//			pars.pFinderPars_.klen = bestKmerLength;
//	//		if(!pars.optimizeOnIteration && iterNumber > 1){
//			if (!pars.pFinderPars_.optimizeOnIteration) {
//				pars.pFinderPars_.optimizeKcutStop = bestKmerCut;
//				pars.pFinderPars_.optimizeKLenStop = bestKmerLength;
//			}
		}
//		bool print = false;
		if (firstMatesInput.empty() && singlesInputs.empty()) {
			//no reads extracted
			if(optimizeAgainAndRecruitBeforeFinal && !alreadyTriedToOptimizedAgain){
				//if(pars.filterOffOutlierInputSeqs && pars.removePossibleOutliersKSim && iterNumber <= 2){
				if(iterNumber <= 1){
					iterate = false;
					break;
				}else{
					alreadyTriedToOptimizedAgain = true;
					pars.pFinderPars_.kmerLengths = originalKmerLengths;
					pars.pFinderPars_.kmerKOcurrenceCutOffs = originalKmerCutOffs;
					pars.pFinderPars_.shortTipNumbers = originalShorttipNumbers;
					//with the re-optimization, add back in the filtered off seqs if they were
					for(const auto & iterDir : iterationDirectories){
						auto outlierFilteredOff_extractedSinglesFnp = njh::files::make_path(iterDir, "outlierFilteredOff_extractedSingles.fastq");
						auto outlierFilteredOff_extractedPairs_R1Fnp = njh::files::make_path(iterDir, "outlierFilteredOff_extractedPairs_R1.fastq");
						auto outlierFilteredOff_extractedPairs_R2Fnp = njh::files::make_path(iterDir, "outlierFilteredOff_extractedPairs_R2.fastq");
						//pairs filtered off from previous
						if(bfs::exists(outlierFilteredOff_extractedPairs_R1Fnp)){
							firstMatesInput.insert(firstMatesInput.begin(), outlierFilteredOff_extractedPairs_R1Fnp);
							secondMatesInput.insert(secondMatesInput.begin(), outlierFilteredOff_extractedPairs_R2Fnp);
						}
						//singles filtered off from previous
						if (bfs::exists(outlierFilteredOff_extractedSinglesFnp)){
							singlesInputs.insert(singlesInputs.begin(), outlierFilteredOff_extractedSinglesFnp);
						}
					}
				}
			} else{
				iterate = false;
				break;
			}
		}



		//pairs from previous
		if(lastIteractionPaired.inExists()){
			firstMatesInput.insert(firstMatesInput.begin(), lastIteractionPaired.firstName_);
			secondMatesInput.insert(secondMatesInput.begin(), lastIteractionPaired.secondName_);
		}
		//singles from previous
		if (lastIteractionSingles.inExists()){
			singlesInputs.insert(singlesInputs.begin(), lastIteractionSingles.firstName_);
		}
		if (lastIteractionSinglesMarkedDups.inExists()){
			singlesInputs.insert(singlesInputs.begin(), lastIteractionSinglesMarkedDups.firstName_);
		}


//		if(print){
//			std::cout << "first mates: " << std::endl;
//			for(const auto & firstMate : firstMatesInput){
//				std::cout << firstMate << std::endl;
//			}
//			std::cout << "singles: " << std::endl;
//			for(const auto & single : singlesInputs){
//				std::cout << single << std::endl;
//			}
//			exit(1);
//		}


		if(!firstMatesInput.empty()){
			concatenateFiles(firstMatesInput, OutOptions(nextInteractionPairedInput.firstName_));
			concatenateFiles(secondMatesInput, OutOptions(nextInteractionPairedInput.secondName_));
			if(!keepIntermediatefiles){
				//clean up
				for(const auto & fnp : firstMatesInput){
					if(fnp != lastIteractionPaired.firstName_){
						bfs::remove(fnp);
					}
				}
				for(const auto & fnp : secondMatesInput){
					if(fnp != lastIteractionPaired.secondName_){
						bfs::remove(fnp);
					}
				}
			}
		}
		if(!singlesInputs.empty()){
			concatenateFiles(singlesInputs, OutOptions(nextInteractionSinglesInput.firstName_));
			if(!keepIntermediatefiles){
				//clean up
				for(const auto & fnp : singlesInputs){
					if (fnp != lastIteractionSingles.firstName_
							&& fnp
									!= njh::files::make_path(setUp.pars_.directoryName_,
											"thrownAwayMate_extracted.fastq")
							&& fnp
									!= njh::files::make_path(setUp.pars_.directoryName_,
											"filteredSingles_extracted.fastq")) {
						bfs::remove(fnp);
					}
				}
			}
		}
		BamExtractor::ExtractedFilesOpts iterationOpts;
		iterationOpts.inPairs_ = nextInteractionPairedInput;
		iterationOpts.inUnpaired_ = nextInteractionSinglesInput;




		auto currentPwRes = PathFinderFromSeqsDev(iterationOpts, setUp.pars_.directoryName_, currentIterationSampleName, pars.pFinderPars_, meta, lastIterationPwRes);
		iterationDirectories.emplace_back(njh::files::make_path(setUp.pars_.directoryName_, currentIterationSampleName));
		bestKmerCut = currentPwRes.log_["bestKmerOccurenceCutOff"].asInt();
		bestKmerLength = currentPwRes.log_["bestKmerLength"].asInt();
		bestShortTip = currentPwRes.log_["bestShortTipNumber"].asInt();

		if(!keepIntermediatefiles){
			//clean up
			iterationOpts.removeAllInFiles();

			std::vector<bfs::path> oldInputFiles{
				"extractedPairsWithNoTandems_R1.fastq",
				"extractedPairsWithNoTandems_R2.fastq",
				"extractedPairsWithTandems_R1.fastq",
				"extractedPairsWithTandems_R2.fastq",
				"extractedPairs_R1.fastq",
				"extractedPairs_R2.fastq",
				"extractedSingles.fastq",
				"extractedSinglesWithNoTandems.fastq",
				"extractedSinglesWithTandems.fastq",
				"filteredExtractedPairs_R1.fastq",
				"filteredExtractedPairs_R2.fastq",
				"filteredSingles.fastq",
			  "filteredOff_extractedSingles.fastq",
				"filteredOff_extractedPairs_R1.fastq",
				"filteredOff_extractedPairs_R2.fastq",
			  "filteredOffDups_extractedSingles.fastq",
			  "filteredOffDups_extractedPairs_R1.fastq",
				"filteredOffDups_extractedPairs_R2.fastq",
			  "outlierFilteredOff_extractedPairs",
			  "outlierFilteredOff_extractedSingles"};

			for(const auto & fn : oldInputFiles){
				auto fnp = njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, fn);
				if(bfs::exists(fnp)){
					bfs::remove(fnp);
				}
			}
		}
		fullLog[currentIterationSampleName] = currentPwRes.log_;
		lastIterationSampName = currentIterationSampleName;
		lastIterationPwRes = currentPwRes;
		if(alreadyTriedToOptimizedAgain){
			++iterCountAfterReOptimization;
		}

		if(alreadyTriedToOptimizedAgain && 1 == iterCountAfterReOptimization ){
			if(pars.pFinderPars_.removePossibleOutliersKSim && pars.pFinderPars_.filterOffOutlierInputSeqs){
				//reoptimize if seqs for possible outliers were removed
				auto lastIteractionPaired = SeqIOOptions::genPairedIn(
				njh::files::make_path(
						setUp.pars_.directoryName_, lastIterationSampName, "filteredExtractedPairs_R1.fastq"),
				njh::files::make_path(
						setUp.pars_.directoryName_, lastIterationSampName, "filteredExtractedPairs_R2.fastq"));
				auto lastIteractionSingles = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "filteredSingles.fastq"));

				std::string currentIterationSampleName = njh::pasteAsStr(lastIterationSampName, "-reoptimized");
				BamExtractor::ExtractedFilesOpts iterationOpts;
				iterationOpts.inPairs_ = lastIteractionPaired;
				iterationOpts.inUnpaired_ = lastIteractionSingles;
				auto currentPwRes = PathFinderFromSeqsDev(iterationOpts, setUp.pars_.directoryName_, currentIterationSampleName, pars.pFinderPars_, meta, lastIterationPwRes);
				iterationDirectories.emplace_back(njh::files::make_path(setUp.pars_.directoryName_, currentIterationSampleName));
				bestKmerCut = currentPwRes.log_["bestKmerOccurenceCutOff"].asInt();
				bestKmerLength = currentPwRes.log_["bestKmerLength"].asInt();
				bestShortTip = currentPwRes.log_["bestShortTipNumber"].asInt();
				if(!keepIntermediatefiles){
					//clean up
					iterationOpts.removeAllInFiles();

					std::vector<bfs::path> oldInputFiles{
						"extractedPairsWithNoTandems_R1.fastq",
						"extractedPairsWithNoTandems_R2.fastq",
						"extractedPairsWithTandems_R1.fastq",
						"extractedPairsWithTandems_R2.fastq",
						"extractedPairs_R1.fastq",
						"extractedPairs_R2.fastq",
						"extractedSingles.fastq",
						"extractedSinglesWithNoTandems.fastq",
						"extractedSinglesWithTandems.fastq",
						"filteredExtractedPairs_R1.fastq",
						"filteredExtractedPairs_R2.fastq",
						"filteredSingles.fastq",
					  "filteredOff_extractedSingles.fastq",
						"filteredOff_extractedPairs_R1.fastq",
						"filteredOff_extractedPairs_R2.fastq",
					  "filteredOffDups_extractedSingles.fastq",
					  "filteredOffDups_extractedPairs_R1.fastq",
						"filteredOffDups_extractedPairs_R2.fastq",
					  "outlierFilteredOff_extractedPairs",
					  "outlierFilteredOff_extractedSingles"};

					for(const auto & fn : oldInputFiles){
						auto fnp = njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, fn);
						if(bfs::exists(fnp)){
							bfs::remove(fnp);
						}
					}
				}
				fullLog[currentIterationSampleName] = currentPwRes.log_;
				lastIterationSampName = currentIterationSampleName;
				lastIterationPwRes = currentPwRes;
			}
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "iterate: " << njh::colorBool(iterate) << std::endl;
//		std::cout << "iterNumber: " << iterNumber << std::endl;
		if(iterNumber >= maxIteration){
			break;
		}
		++iterNumber;
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "iterate: " << njh::colorBool(iterate) << std::endl;
//		std::cout << "iterNumber: " << iterNumber << std::endl;
	}


	if(!keepIntermediatefiles){
		//clean up
		rextractedSeqs.removeAllInFiles();
		rextractedSeqsSingles.removeAllInFiles();
	}


	if(noSplittingTailed){
		pars.pFinderPars_.splitTailed_ = false;
	}
	if(noSplittingEnd){
		pars.pFinderPars_.splitEnds_ = false;
	}
	if(noAddHeadTail){
		pars.pFinderPars_.addSeqOfSingleHeadAndTailSeqs_ = false;
	}

	//reset the original cut off and lengths to re-optimize once more
	if (0 == maxIteration) {
		pars.pFinderPars_.kmerLengths = {bestKmerLength};
		pars.pFinderPars_.kmerKOcurrenceCutOffs = {bestKmerCut};
		pars.pFinderPars_.shortTipNumbers = {bestShortTip};
	} else {
		pars.pFinderPars_.kmerLengths = originalKmerLengths;
		pars.pFinderPars_.kmerKOcurrenceCutOffs = originalKmerCutOffs;
		pars.pFinderPars_.shortTipNumbers = originalShorttipNumbers;
	}



	std::string currentIterationSampleName = sampName + "-finalPass";
	BamExtractor::ExtractedFilesOpts iterationOpts;
//	iterationOpts.inPairs_ =   SeqIOOptions::genPairedIn(njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "extractedPairs_R1.fastq"),
//																									     njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "extractedPairs_R2.fastq"));
//	iterationOpts.inUnpaired_ = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "extractedSingles.fastq"));
	iterationOpts.inPairs_ =   SeqIOOptions::genPairedIn(njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "filteredExtractedPairs_R1.fastq"),
																									     njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "filteredExtractedPairs_R2.fastq"));
	iterationOpts.inUnpaired_ = SeqIOOptions::genFastqIn(njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, "filteredSingles.fastq"));

	//auto currentLog = PathFinderFromSeqsDev(iterationOpts, setUp.pars_.directoryName_, currentIterationSampleName, pars, meta);
	if(trimFinalEnds){
		pars.pFinderPars_.trimEnds = true;
		pars.pFinderPars_.trimEndsBy = trimFinalEndsBy;
	}
	pars.pFinderPars_.headlessTailessLenCutOff = finalHeadlessTaillessCutOff;
	auto currentPwRes = PathFinderFromSeqsDev(iterationOpts, setUp.pars_.directoryName_, currentIterationSampleName, pars.pFinderPars_, meta, lastIterationPwRes);
	auto finalPassDir = njh::files::make_path(setUp.pars_.directoryName_, currentIterationSampleName);
	fullLog[currentIterationSampleName] = currentPwRes.log_;


	fullLog["final-bestKmerOccurenceCutOff"] = currentPwRes.log_["bestKmerOccurenceCutOff"].asInt();
	fullLog["final-bestKmerLength"] = currentPwRes.log_["bestKmerLength"].asInt();
	fullLog["final-bestShortTipNumber"] = currentPwRes.log_["bestShortTipNumber"].asInt();


	if(!keepIntermediatefiles){
		std::vector<bfs::path> oldInputFiles{
			"extractedPairsWithNoTandems_R1.fastq",
			"extractedPairsWithNoTandems_R2.fastq",
			"extractedPairsWithTandems_R1.fastq",
			"extractedPairsWithTandems_R2.fastq",
			"extractedPairs_R1.fastq",
			"extractedPairs_R2.fastq",
			"extractedSingles.fastq",
			"extractedSinglesWithNoTandems.fastq",
			"extractedSinglesWithTandems.fastq",
			"filteredExtractedPairs_R1.fastq",
			"filteredExtractedPairs_R2.fastq",
			"filteredSingles.fastq",
		  "filteredOff_extractedSingles.fastq",
			"filteredOff_extractedPairs_R1.fastq",
			"filteredOff_extractedPairs_R2.fastq",
		  "filteredOffDups_extractedSingles.fastq",
		  "filteredOffDups_extractedPairs_R1.fastq",
			"filteredOffDups_extractedPairs_R2.fastq",
		  "outlierFilteredOff_extractedPairs",
		  "outlierFilteredOff_extractedSingles"};


		for(const auto & fn : oldInputFiles){
			auto fnp = njh::files::make_path(setUp.pars_.directoryName_, lastIterationSampName, fn);
			if(bfs::exists(fnp)){
				bfs::remove(fnp);
			}
		}
	}
	lastIterationSampName = currentIterationSampleName;
	lastIterationPwRes = currentPwRes;


	fullLog["lastIterationSampName"] = lastIterationSampName;
	//read in final iteration seqs
	auto finalSeqInOpts = SeqIOOptions::genFastaIn(njh::files::make_path(
			setUp.pars_.directoryName_, lastIterationSampName,
			"output_aboveCutOff.fasta"));
	finalSeqInOpts.processed_ = true;
	uint64_t maxLen = 0;
	std::vector<readObject> finalSeqs;
	int returnStatus = 0;
	if (finalSeqInOpts.inExists()
			&& 0 != bfs::file_size(finalSeqInOpts.firstName_)) {
		auto finalSeqsRaw = SeqInput::getReferenceSeq(finalSeqInOpts, maxLen);
		//get final seqs above length cut off
		for (const auto & finalSeq : finalSeqsRaw) {
			if(MetaDataInName::nameHasMetaData(finalSeq.seqBase_.name_)){
				MetaDataInName finalSeqMeta(finalSeq.seqBase_.name_);
				if(finalSeqMeta.containsMeta("trimStatus") && !finalSeqMeta.getMeta<bool>("trimStatus")){
					continue;
				}
			}
			if (len(finalSeq) >= finalLenCutOff) {
				finalSeqs.emplace_back(finalSeq);
			}
		}
	} else {
		returnStatus = 1;
		std::cerr << "No sequences reconstructed for given parameters" << std::endl;
	}
	if(!finalSeqs.empty()){
		readVec::allSetFractionByTotalCount(finalSeqs);
		//update length meta
		for (auto & clus : finalSeqs) {
			MetaDataInName clusMeta(clus.seqBase_.name_);
			clusMeta.addMeta("length", len(clus), true);
			clusMeta.addMeta("sample", sampName, true);
			clusMeta.resetMetaInName(clus.seqBase_.name_);
		}
		//sort and put longest contig on top
		njh::sort(finalSeqs,
				[](const readObject & clus1, const readObject clus2) {
					return len(clus1) > len(clus2);
				});
		if(redetermineReadCounts){
			//pairs
			bfs::path lastIterationDir = njh::files::make_path(
					setUp.pars_.directoryName_, lastIterationSampName);
			auto r1Fnp = njh::files::make_path(lastIterationDir,       "filteredExtractedPairs_R1.fastq");
			auto r2Fnp = njh::files::make_path(lastIterationDir,       "filteredExtractedPairs_R2.fastq");
			auto singlesFnp = njh::files::make_path(lastIterationDir,  "filteredSingles.fastq");
			uint64_t maxLen = 0;
			readVec::getMaxLength(finalSeqs, maxLen);

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

			aligner alignerObj(maxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2), false);

			std::vector<refVariants> refVariationInfo;
			for (const auto refPos : iter::range(finalSeqs.size())) {
				refVariationInfo.emplace_back(finalSeqs[refPos].seqBase_);
			}
			for (const auto refPos : iter::range(finalSeqs.size())) {
				for (const auto refSubPos : iter::range(finalSeqs.size())) {
					if (refPos == refSubPos) {
						continue;
					}
					refVariationInfo[refPos].addVariant(finalSeqs[refSubPos].seqBase_,
							alignerObj, false);
				}
			}
			alignerObj.weighHomopolymers_ = true;
			alignerObj.countEndGaps_ = false;

			//
			std::vector<kmerInfo> refInfos;
			double matchingKmerCutOff = 0.90;
			uint32_t matchingKmerLen = 5;
			for (const auto& clus : finalSeqs){
				refInfos.emplace_back(clus.seqBase_.seq_, matchingKmerLen, false);
			}

			uint32_t unmappable = 0;
			uint32_t indeterminate = 0;
			uint32_t multiMapping = 0;
			std::unordered_map<std::string, std::set<std::string>> uniqueSeqs;
			//std::unordered_map<std::string, std::set<std::string>> multiMappingSeqs;
			std::unordered_map<std::string, std::vector<seqInfo>> multiMappingSeqs;
			std::set<std::string> allReadInSamples;
			std::unordered_map<std::string, std::map<uint32_t, std::unordered_map<std::string, double>>> baseCoverrage;
			bfs::path finalDirectoryOut = setUp.pars_.directoryName_;
			auto increaseMapCounts = [&unmappable,&indeterminate,&uniqueSeqs,&multiMappingSeqs, & refInfos,&finalSeqs,&alignerObj,
															 &mapAllowableErrors,&refVariationInfo,&matchingKmerLen,&matchingKmerCutOff,&allReadInSamples,&multiMapping,
															 &baseCoverrage](const seqInfo & reMapSeq){
				allReadInSamples.emplace(reMapSeq.name_);
				//best score to be used to get all the references that map within a certain distance to this latter
				double bestScore = 0;
				kmerInfo reMappingSeqInfo(reMapSeq.seq_, matchingKmerLen, false);
				std::vector<uint64_t> bestRefs;
				uint32_t subSize = len(reMapSeq) * mapAllowableErrors.distances_.query_.coverage_;


				for (const auto refPos : iter::range(finalSeqs.size())) {
					if(reMappingSeqInfo.compareKmers(refInfos[refPos]).second < matchingKmerCutOff &&
						 reMappingSeqInfo.compareSubKmersToFull(refInfos[refPos], 0, subSize).second < matchingKmerCutOff &&
						 reMappingSeqInfo.compareSubKmersToFull(refInfos[refPos], len(reMapSeq) - subSize, subSize).second < matchingKmerCutOff) {
						continue;
					}
//					if(reMappingSeqInfo.compareSubKmersToFull(refInfos[refPos],
//									0,
//									subSize).second < matchingKmerCutOff &&
//
//							reMappingSeqInfo.compareSubKmersToFull(refInfos[refPos],
//									len(reMapSeq) - subSize,
//									subSize).second < matchingKmerCutOff) {
//
//						continue;
//					}
//					std::cout << "subSize: " << subSize << std::endl;
//					std::cout << "refSeqBackPos: " << (subSize<outSeqs[refPos].seq_.size() ? outSeqs[refPos].seq_.size() - subSize: 0) << "/" << outSeqs[refPos].seq_.size() << std::endl;
//					std::cout << "mapSeqBackPos: " << len(reMapSeq) - subSize << "/" << len(reMapSeq) << std::endl;
//					std::cout << "full: " << reMappingSeqInfo.compareKmers(refInfos[refPos]).second  << std::endl;
//					std::cout << "front: " << reMappingSeqInfo.compareKmers(refInfos[refPos], 0, subSize<outSeqs[refPos].seq_.size() ? outSeqs[refPos].seq_.size() - subSize: 0, subSize).second << std::endl;
//					std::cout << "back:  " << reMappingSeqInfo.compareKmers(refInfos[refPos], len(reMapSeq) - subSize, 0, subSize).second << " " << njh::colorBool(reMappingSeqInfo.compareKmers(refInfos[refPos], len(reMapSeq) - subSize, 0, subSize).second >+ matchingKmerCutOff) << std::endl;

					const auto & ref = finalSeqs[refPos].seqBase_;
					alignerObj.alignCacheGlobal(ref,reMapSeq);
					alignerObj.profileAlignment(ref,reMapSeq, false, true, false);
					if (alignerObj.comp_.distances_.query_.coverage_ >= mapAllowableErrors.distances_.query_.coverage_ &&
							mapAllowableErrors.passErrorProfile(alignerObj.comp_)) {
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
					uniqueSeqs[finalSeqs[bestRefs.front()].seqBase_.name_].emplace(reMapSeq.name_);
					alignerObj.alignCacheGlobal(finalSeqs[bestRefs.front()].seqBase_,reMapSeq);
					for(const auto seqAlignPos : iter::range(alignerObj.alignObjectA_.seqBase_.seq_.size())){
						if(alignerObj.alignObjectA_.seqBase_.seq_[seqAlignPos] != '-' &&
							 alignerObj.alignObjectB_.seqBase_.seq_[seqAlignPos] != '-'){
							baseCoverrage[finalSeqs[bestRefs.front()].seqBase_.name_][alignerObj.getSeqPosForAlnAPos(seqAlignPos)][reMapSeq.name_] = 1;
						}
					}
					//tempWriter.openWrite(reMapSeq);
				} else if (bestRefs.size() == 0) {
					//if no mappable candidates found, put into un-mappable
					++unmappable;
				} else {
					VecStr matchingRefNames;
					for(const auto & refPos : bestRefs) {
						matchingRefNames.emplace_back(finalSeqs[refPos].seqBase_.name_);
					}
					std::vector<uint64_t> matchingRefs;
					for (const auto & refPos : bestRefs) {
						const auto & ref = finalSeqs[refPos].seqBase_;
						//get the current snp info for the current ref to the other matching refs
						auto seqSnpPosBases = refVariationInfo[refPos].getVariantSnpLociMap(matchingRefNames, 0);
						alignerObj.alignCacheGlobal(ref, reMapSeq);
						alignerObj.profileAlignment(ref,reMapSeq, false, false, false);
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
						uniqueSeqs[finalSeqs[matchingRefs.front()].seqBase_.name_].emplace(reMapSeq.name_);
						alignerObj.alignCacheGlobal(finalSeqs[matchingRefs.front()].seqBase_,reMapSeq);
						for(const auto seqAlignPos : iter::range(alignerObj.alignObjectA_.seqBase_.seq_.size())){
							if(alignerObj.alignObjectA_.seqBase_.seq_[seqAlignPos] != '-' &&
								 alignerObj.alignObjectB_.seqBase_.seq_[seqAlignPos] != '-'){
								baseCoverrage[finalSeqs[matchingRefs.front()].seqBase_.name_][alignerObj.getSeqPosForAlnAPos(seqAlignPos)][reMapSeq.name_] = 1;
							}
						}
						//tempWriter.openWrite(reMapSeq);
					} else if(matchingRefs.size() == 0) {
						//if none of the ref work out, put it in indeterminate
						++indeterminate;
					} else {
						//tempWriter.openWrite(reMapSeq);
						//if there is more than one matching ref, put it in ties
						njh::sort(matchingRefs);
						std::string refMixName = njh::conToStr(matchingRefs, "_");
						multiMappingSeqs[refMixName].emplace_back(reMapSeq);
						++multiMapping;
					}
				}
			};
			seqInfo reMappingSeq;
			if(bfs::exists(r1Fnp)){
				SeqInput remapReader(SeqIOOptions::genFastqIn(r1Fnp));
				remapReader.openIn();
				while(remapReader.readNextRead(reMappingSeq)){
					increaseMapCounts(reMappingSeq);
				}
			}
			if(bfs::exists(r2Fnp)){
				SeqInput remapReader(SeqIOOptions::genFastqIn(r2Fnp));
				remapReader.openIn();
				while(remapReader.readNextRead(reMappingSeq)){
					increaseMapCounts(reMappingSeq);
				}
			}
			if(bfs::exists(singlesFnp)) {
				SeqInput remapReader(SeqIOOptions::genFastqIn(singlesFnp));
				remapReader.openIn();
				while(remapReader.readNextRead(reMappingSeq)) {
					increaseMapCounts(reMappingSeq);
				}
			}

			std::unordered_map<std::string, double> uniqueMapFrac;
			std::unordered_map<std::string, double> finalTotals;
			double totalUniques = 0;
			for(const auto & uni : uniqueSeqs){
				totalUniques += uni.second.size();
				finalTotals[uni.first] = uni.second.size();
			}
			for(const auto & outSeq : finalSeqs){
				if(njh::in(outSeq.seqBase_.name_, uniqueSeqs)){
					uniqueMapFrac[outSeq.seqBase_.name_] = uniqueSeqs[outSeq.seqBase_.name_].size()/totalUniques;
				}else{
					/**@todo might want to use a different amount or just get rid of the read, worried about low coverage*/
					uniqueMapFrac[outSeq.seqBase_.name_] = 1/totalUniques;
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
					totalFrac += uniqueMapFrac[finalSeqs[pos].seqBase_.name_];
				}
				auto names = readVec::getNames(multi.second);
				std::set<std::string> namesSet(names.begin(), names.end());

				for(const auto & pos : refPoss){
//											std::cout << "pos: " << pos << std::endl;
//											std::cout << "multi.second.size(): " << multi.second.size() << std::endl;
//											std::cout << "multi.second.size() * (uniqueMapFrac[finalSeqsClusters[pos].seqBase_.name_]/totalFrac): " <<  multi.second.size() * (uniqueMapFrac[finalSeqsClusters[pos].seqBase_.name_]/totalFrac) << std::endl;
					finalTotals[finalSeqs[pos].seqBase_.name_] += namesSet.size() * (uniqueMapFrac[finalSeqs[pos].seqBase_.name_]/totalFrac);
					for(const auto & reMapSeq : multi.second){
						alignerObj.alignCacheGlobal(finalSeqs[pos].seqBase_, reMapSeq);
						for(const auto seqAlignPos : iter::range(alignerObj.alignObjectA_.seqBase_.seq_.size())){
							if(alignerObj.alignObjectA_.seqBase_.seq_[seqAlignPos] != '-' &&
								 alignerObj.alignObjectB_.seqBase_.seq_[seqAlignPos] != '-'){
								baseCoverrage[finalSeqs[pos].seqBase_.name_][alignerObj.getSeqPosForAlnAPos(seqAlignPos)].emplace(reMapSeq.name_, uniqueMapFrac[finalSeqs[pos].seqBase_.name_]/totalFrac);
							}
						}
					}
				}
			}
			uint32_t totalReadCount = 0;
			for(auto & outSeq : finalSeqs){
				outSeq.seqBase_.cnt_ = std::round(finalTotals[outSeq.seqBase_.name_]);
				totalReadCount+= outSeq.seqBase_.cnt_;
			}
			std::unordered_map<std::string, uint32_t> nameKey;
			uint32_t pos = 0;
			for(auto & outSeq : finalSeqs){
				std::string oldName = outSeq.seqBase_.name_;
				outSeq.seqBase_.frac_ = outSeq.seqBase_.cnt_/totalReadCount;
				outSeq.updateName();
				nameKey[oldName] = pos;
				++pos;
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
				//if(outSeqs[bCov.first].seq_.size() > 60){
					meanCoverage = vectorMean(getSubVector(bCov.second, 25, finalSeqs[bCov.first].seqBase_.seq_.size() - 50));
				}else{
					meanCoverage = vectorMean(bCov.second);
				}
				MetaDataInName nameMeta(finalSeqs[bCov.first].seqBase_.name_);
				nameMeta.addMeta("meanBaseCoverage", meanCoverage, true);
				nameMeta.resetMetaInName(finalSeqs[bCov.first].seqBase_.name_);
			}


			OutOptions baseCoverageOutOpts(njh::files::make_path(finalDirectoryOut,"baseCoverage.tab.txt"));
			OutputStream baseCoverageOut(baseCoverageOutOpts);
			baseCoverageOut << "name\tposition\tcoverage" << std::endl;
			for(const auto & b : baseCoverrage){
				for(const auto & pos : b.second){

					double baseTotal = 0;
					for(const auto & readWeight : pos.second){
						baseTotal += readWeight.second;
					}
					baseCoverageOut << finalSeqs[nameKey[b.first]].seqBase_.name_
							<< "\t" << pos.first
							<< "\t" << baseTotal << std::endl;
				}
			}
			Json::Value mapCounts;
			mapCounts["unmappable"]    = njh::json::toJson(unmappable);
			mapCounts["indeterminate"] = njh::json::toJson(indeterminate);
			mapCounts["totalUniques"]  = njh::json::toJson(totalUniques);
			mapCounts["multiMapping"] = njh::json::toJson(multiMapping);
			mapCounts["totalReadIn"]  =  njh::json::toJson(allReadInSamples.size());
			OutputStream mapCountsOut(OutOptions(njh::files::make_path(finalDirectoryOut,"redeteringCountsMapCounts.json")));
			mapCountsOut << mapCounts << std::endl;

		}


		//write out final seqs; printTandems
		OutOptions finalInfoOpts(
				njh::files::make_path(setUp.pars_.directoryName_,
						"finalOutputInfo.tab.txt"));
		std::ofstream finalInfoFile;
		finalInfoOpts.openFile(finalInfoFile);
		finalInfoFile << "name\tlength\treadCount\tfraction" << "\n";
		auto finalSeqOutOpts = SeqIOOptions::genFastaOut(
				njh::files::make_path(setUp.pars_.directoryName_, "finalOutput.fasta"));
		SeqOutput finalWriter(finalSeqOutOpts);
		finalWriter.openOut();

		for (const auto & clus : finalSeqs) {
			finalWriter.write(clus);
			finalInfoFile << clus.seqBase_.name_ << "\t" << len(clus) << "\t"
					<< clus.seqBase_.cnt_ << "\t" << clus.seqBase_.frac_ << std::endl;
		}
	} else {
		if(0 == returnStatus ){
			returnStatus = 1;
			std::cerr << "No reads above the final read length cut off" << std::endl;
		}
	}

	//compress original files to save space
	{

		std::vector<bfs::path> originalExtractionFiles{
			"extracted.fastq",
			"extracted_R1.fastq",
			"extracted_R2.fastq",
			"inverse_extracted_R1.fastq",
			"inverse_extracted_R2.fastq",
			"mateUnmapped_extracted_R1.fastq",
			"mateUnmapped_extracted_R2.fastq",
		  "thrownAwayMate_extracted.fastq",
			"filteredSingles_extracted.fastq",
		  "filteredPairs_extracted_R1.fastq",
		  "filteredPairs_extracted_R2.fastq",
			"filteredSoftClipSingles_extracted.fastq"};

		njh::concurrent::LockableQueue<bfs::path> fileQueque(originalExtractionFiles);
		if(development){

			bfs::path originalExtractionFilesDir = njh::files::make_path(setUp.pars_.directoryName_, "originalExtractionFiles");
			njh::files::makeDir(njh::files::MkdirPar{originalExtractionFilesDir});
			auto compressFile = [&fileQueque,&setUp,&originalExtractionFilesDir](){
				bfs::path fnp;
				while(fileQueque.getVal(fnp)){
					auto fnpInDir = njh::files::make_path(setUp.pars_.directoryName_, fnp);
					if(bfs::exists(fnpInDir)){
						auto fnpInDirGz = njh::files::make_path(originalExtractionFilesDir, fnp.string()  + ".gz");
						OutputStream outGz(fnpInDirGz);
						InputStream in(fnpInDir);
						outGz << in.rdbuf();
					}
					if(bfs::exists(fnpInDir)){
						bfs::remove(fnpInDir);
					}
				}
			};
			std::vector<std::thread> threads;
			for(uint32_t t = 0; t < pars.pFinderPars_.numThreads; ++t){
				threads.emplace_back(std::thread(compressFile));
			}
			njh::concurrent::joinAllThreads(threads);
			auto bamFnp = njh::files::make_path(setUp.pars_.directoryName_, "extracted.bam");
			auto mvBamFnp = njh::files::make_path(originalExtractionFilesDir, "extracted.bam");
			bfs::rename(bamFnp, mvBamFnp);
		}else{
			bfs::path fnp;
			while(fileQueque.getVal(fnp)){
				auto fnpInDir = njh::files::make_path(setUp.pars_.directoryName_, fnp);
				if(bfs::exists(fnpInDir)){
					bfs::remove(fnpInDir);
				}
			}
			auto bamFnp = njh::files::make_path(setUp.pars_.directoryName_, "extracted.bam");
			if(bfs::exists(bamFnp)){
				bfs::remove(bamFnp);
			}
		}
	}

	if(development){
		OutOptions logOutOpts(njh::files::make_path(logDir, "extractionLog.json"));
		auto logFilePtr = logOutOpts.openFile();
		(*logFilePtr) << fullLog << std::endl;
	}else{

		//clean up iteration directories
		for(const auto & iterDir : iterationDirectories){
			njh::files::rmDirForce(iterDir);
		}
		//clean up final pass directory
		{
			auto filesInFinalPass = njh::files::filesInFolder(finalPassDir);
			VecStr filesToKeep = { "output_aboveCutOff.fasta",
					"outputInfo_aboveCutOff.tab.txt", "outputInfo.tab.txt",
					"output.fasta", "optimizationInfoBest.json",
					"optimizationInfoAll.json", "output_final.dot" };
			for(const auto & f : filesInFinalPass){
				if(!njh::in(f.filename().string(), filesToKeep)){
					njh::files::rmDirForce(f);
				}
			}
		}
		std::vector<bfs::path> filesToCleanUp {
			  "trimAlnCache", "logs",
				"inputRegions.bed", "inputRegions.fasta",
				"usedExpandedRegions.bed", "realignment"};
//		std::vector<bfs::path> filesToCleanUp { "trimAlnCache", "logs",
//				"extractionStats", "inputRegions.bed", "inputRegions.fasta",
//				"usedExpandedRegions.bed" };
		for(const auto & fnp : filesToCleanUp){
			auto fullFnp = njh::files::make_path(setUp.pars_.directoryName_, fnp);
			if(bfs::exists(fullFnp)){
				njh::files::rmDirForce(fullFnp);
			}
		}
	}



	return returnStatus;
}





} // namespace njhseq
