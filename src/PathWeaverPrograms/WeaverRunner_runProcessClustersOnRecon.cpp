/*
 * WeaverRunner_runProcessClustersOnRecon.cpp
 *
 *  Created on: Nov 19, 2020
 *      Author: nick
 */




#include "WeaverRunner.hpp"


#include <cppitertools/sorted.hpp>
#include <SeekDeepPrograms/SeekDeepProgram/SeekDeepSetUp.hpp>
#include <njhseq/objects/dataContainers/tables/TableReader.hpp>
#include <njhseq/objects/Meta.h>

#include "PathWeaver/PostWeavingUtils/GatheringUtils.hpp"

#include <chrono>

namespace njhseq {

int WeaverRunner::rawGatherSeqs(const njh::progutils::CmdArgs & inputCommands) {
	SeqGatheringFromPathWeaver::gatherSeqsAndSortByTargetPars rawGatherPars;
	SeqGatheringFromPathWeaver::processedGatherSeqsMetaPars processedPars;
	SeqGatheringFromPathWeaver::SeqGatheringFromPathWeaverCorePars corePars;


	std::set<std::string> samples;
	std::set<std::string> targets;
	bfs::path inputDirectory = "./";
	std::string pat = "";
	bfs::path groupingsFile;

	corePars.countField = "estimatedPerBaseCoverage";
	corePars.numThreads = 1;
	rawGatherPars.addPartial = false;


	seqSetUp setUp(inputCommands);
	// parameters
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(pat, "--pat","The results directory pattern to process, directories must end with this, the prefix to this pattern will be treated as the sample name", true);
	setUp.setOption(inputDirectory, "--inputDirectory", "Input Directory to search");
	setUp.setOption(samples, "--samples", "Process input from only these samples");
	setUp.setOption(targets, "--targets", "Process input for only these targets");
	setUp.setOption(rawGatherPars.addPartial, "--addPartial", "Add Partial");
	setUp.setOption(corePars.numThreads, "--numThreads", "number threads");
	setUp.setOption(corePars.countField, "--countField", "countField");
	setUp.setOption(processedPars.keepCommonSeqsWhenFiltering, "--keepCommonSeqsWhenFiltering", "Keep Common Seqs When Filtering");

	setUp.setOption(groupingsFile, "--groupingsFile,--meta", "Groupings File");

	setUp.processDirectoryOutputName(pat+ "_gatheredSeqs", true);
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	corePars.sampleField = "sample";
	corePars.targetField = "regionUID";

	std::vector<bfs::path> directories;
	if (samples.empty()) {
		auto allFiles = njh::files::listAllFiles(inputDirectory, false, {std::regex { ".*" + pat + "$" } });
		for (const auto &f : allFiles) {
			if (f.second) {
				directories.emplace_back(f.first);
			}
		}
	} else {
		for (const auto &samp : samples) {
			directories.emplace_back(njh::files::make_path(inputDirectory, samp + pat));
		}
	}

	if(directories.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no directories found in " << inputDirectory << " ending with " << pat << "\n";
		throw std::runtime_error{ss.str()};
	}


	std::shared_ptr<MultipleGroupMetaData> meta;
	if("" != groupingsFile){
		meta = std::make_shared<MultipleGroupMetaData>(groupingsFile);
	}

	processedPars.outMetaFnp = njh::files::make_path(setUp.pars_.directoryName_, "metaFnp.tab.txt");
	auto rawAllSeqsFnp = njh::files::make_path(setUp.pars_.directoryName_, "rawAllSeqsFile.fasta");
	auto targetsLocsFnp = njh::files::make_path(setUp.pars_.directoryName_, "targetLocations.tab.txt");
	//
	rawGatherPars.allSeqFnp = rawAllSeqsFnp;
	rawGatherPars.directories = directories;
	corePars.meta = meta;
	rawGatherPars.targets = targets;

	SeqGatheringFromPathWeaver gatherer(corePars);

	auto rawGatherRes = gatherer.gatherSeqsAndSortByTarget(rawGatherPars);
	{
		OutputStream targetsLocsOut(targetsLocsFnp);
		VecStr tarKeys = getVectorOfMapKeys(rawGatherRes.seqsLocations);
		njh::sort(tarKeys);
		for(const auto & tar : tarKeys){
			targetsLocsOut << tar
					<< "\t" << rawGatherRes.seqsLocations[tar].first
					<< "\t" << rawGatherRes.seqsLocations[tar].second
					<< std::endl;
		}
	}

	auto processedGatherRes = gatherer.processedGatherSeqsMeta(processedPars, rawGatherRes);
	if(processedGatherRes.allSeqFnp != rawGatherRes.allSeqFnp){
		VecStr tarKeys = getVectorOfMapKeys(processedGatherRes.seqsLocations);
		njh::sort(tarKeys);
		OutputStream processedTargetsLocsOut(njh::files::make_path(setUp.pars_.directoryName_, "processedTargetLocations.tab.txt"));

		for(const auto & tar : tarKeys){
			processedTargetsLocsOut << tar
					<< "\t" << rawGatherRes.seqsLocations[tar].first
					<< "\t" << rawGatherRes.seqsLocations[tar].second
					<< std::endl;
		}
	}

	if(!processedGatherRes.failedToCombine.empty()){
		OutputStream failedToCombineOut(njh::files::make_path(setUp.pars_.directoryName_, "failedToCombine.txt"));
		failedToCombineOut << "target\tmetaField\tsamples" << std::endl;
		for(const auto & target : iter::sorted(njh::getVecOfMapKeys(processedGatherRes.failedToCombine))){
			for(const auto & field : processedGatherRes.failedToCombine.at(target)){
				failedToCombineOut << target
						<< "\t" << field.first
						<< "\t" << njh::conToStr(field.second, ",") << std::endl;
			}
		}
	}



	return 0;
}


int WeaverRunner::runProcessClustersOnRecon(const njh::progutils::CmdArgs & inputCommands) {
	std::set<std::string> samples;
	std::set<std::string> targets;
	bfs::path inputDirectory = "./";
	std::string pat = "";
	std::string countField = "estimatedPerBaseCoverage";

	bfs::path trimBedFnp = "";
	bfs::path genomeFnp = "";

	bool addPartial = false;
	bool keepCommonSeqsWhenFiltering = false;

	SeekDeepSetUp setUp(inputCommands);
	// parameters
	processClustersPars masterPopClusPars;


	bool swgaSampleClusErrorSet = false;
	bool swgaPopClusErrorSet = false;
	bool skipRBind = false;

	setUp.processDebug();
	setUp.processVerbose();
	bool trimBedSet = setUp.setOption(trimBedFnp, "--trimBedFnp", "Bed File of trim locations, 4th column must match the name of the input targets");
	setUp.setOption(genomeFnp, "--genome2bit", "genome 2bit file for when supplying trim locations", trimBedSet);
	if(trimBedSet){
		if(!bfs::exists(trimBedFnp)){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr(trimBedFnp, " doesn't exist"));
		}
		if(!bfs::exists(genomeFnp)){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr(genomeFnp, " doesn't exist"));
		}
	}
	setUp.setOption(skipRBind, "--skipRBind", "Skip doing large rBind of all basic info files");

	setUp.setOption(keepCommonSeqsWhenFiltering, "--keepCommonSeqsWhenFiltering", "Keep Common Seqs When Filtering");
	setUp.setOption(pat, "--pat","The results directory pattern to process, directories must end with this, the prefix to this pattern will be treated as the sample name", true);
	setUp.setOption(inputDirectory, "--inputDirectory", "Input Directory to search");
	setUp.setOption(samples, "--samples", "Process input from only these samples");
	setUp.setOption(targets, "--targets", "Process input for only these targets");
	setUp.setOption(addPartial, "--addPartial", "Add Partial");

	setUp.setOption(countField, "--countField", "counts field (set equal to \"reads\" to use the total reads for a hap");
	setUp.processDirectoryOutputName(pat+ "_populationClustering", true);


	masterPopClusPars.transPars.setOptions(setUp);

	setUp.setOption(masterPopClusPars.chiCutOff, "--chiCutOff",
			"The Fraction of a cluster to determine if it chimeric", false, "Chimeras");
	setUp.setOption(masterPopClusPars.customCutOffs, "--custumCutOffs",
			"Two Column Table, first column is sample name, second is a custom frac cut off, if sample not found will default to --fracCutOff", false, "Filtering");
	setUp.setOption(masterPopClusPars.previousPopFilename, "--previousPop", "previousPopFilename", false, "Population");
	setUp.processComparison(masterPopClusPars.previousPopErrors, "previousPop");
	setUp.setOption(masterPopClusPars.groupingsFile, "--groupingsFile,--meta",
			"A file to sort samples into different groups", false, "Meta");

	setUp.pars_.colOpts_.verboseOpts_.verbose_ = setUp.pars_.verbose_;
	setUp.pars_.colOpts_.verboseOpts_.debug_ = setUp.pars_.debug_;

	setUp.setOption(masterPopClusPars.noErrorsSet, "--noErrors", "Collapse parameters with no errors", false, "Clustering");


	setUp.setOption(masterPopClusPars.strictErrorsSetHq1, "--strictErrors-hq1", "Collapse parameters with a several low quality mismatches and 1 high quality mismatch", false, "Clustering");
	if(masterPopClusPars.strictErrorsSetHq1){
		masterPopClusPars.strictErrorsSet = true;
		masterPopClusPars.hqMismatches = 1;
	}


	setUp.setOption(swgaSampleClusErrorSet, "--swgaSampleClusErrorSet", "Collapse parameters for each sample with SWGA data (collapsing on homopolymer indel)s", false, "Clustering");

	setUp.setOption(masterPopClusPars.strictErrorsSet, "--strictErrors", "Collapse parameters with a several low quality mismatches", false, "Clustering");
	setUp.setOption(masterPopClusPars.hqMismatches, "--hq", "Number of high quality mismatches to allow", false, "Clustering");
	setUp.setOption(masterPopClusPars.stopAfter, "--stopAfter", "Number of top haplotypes to check", false, "Clustering");

	setUp.setOption(masterPopClusPars.parameters, "--par", "ParametersFileName", !swgaSampleClusErrorSet && !masterPopClusPars.noErrorsSet && !masterPopClusPars.strictErrorsSet && !masterPopClusPars.strictErrorsSetHq1, "Clustering");

	setUp.setOption(masterPopClusPars.binParameters, "--binPar", "bin Parameters Filename", false, "Clustering");


	masterPopClusPars.preFiltCutOffs.sampleMinReadCount = 250;
	bool sampMinSet = setUp.setOption(masterPopClusPars.preFiltCutOffs.sampleMinReadCount, "--sampleMinTotalReadCutOff",
			"Sample Minimum Total Read Cut Off, if the total read count for the sample is below this it will be thrown out", false, "Filtering");

	masterPopClusPars.preFiltCutOffs.replicateMinReadCount = masterPopClusPars.preFiltCutOffs.sampleMinReadCount;
	bool repMinSet = setUp.setOption(masterPopClusPars.preFiltCutOffs.replicateMinReadCount, "--replicateMinTotalReadCutOff",
				"Replicate Minimum Total Read Cut Off, if the total read count for the replicate is below this it will be thrown out", false, "Filtering");
	if(repMinSet && !sampMinSet){
		masterPopClusPars.preFiltCutOffs.sampleMinReadCount = masterPopClusPars.preFiltCutOffs.replicateMinReadCount;
	}
	setUp.setOption(masterPopClusPars.runsRequired, "--runsRequired", "Number of PCR runs Required for a haplotype to be kept", false, "Filtering");
	masterPopClusPars.preFiltCutOffs.clusterSizeCutOff = 10;
	setUp.setOption(masterPopClusPars.preFiltCutOffs.clusterSizeCutOff, "--clusterCutOff", "Input Cluster Size Cut Off", false, "Filtering");
	setUp.setOption(masterPopClusPars.excludeSamples, "--excludeSamples", "Samples to Exclude from analysis", false, "Filtering");
	setUp.setOption(masterPopClusPars.lowLevelPopFiltPars_.removeCommonlyLowFreqHaplotypes_, "--excludeCommonlyLowFreqHaplotypes", "Remove Commonly Low Freq Haplotypes", false, "Filtering");
	setUp.setOption(masterPopClusPars.lowLevelPopFiltPars_.lowFreqHaplotypeFracCutOff_, "--lowFreqHaplotypeFracCutOff", "Low Freq Haplotype Frac Cut Off", false, "Filtering");
	setUp.setOption(masterPopClusPars.experimentNames.controlSamples_, "--controlSamples", "Samples that shouldn't be included in frequency filtering calcs", false, "Filtering");
	setUp.setOption(masterPopClusPars.collapseLowFreqOneOffs, "--excludeLowFreqOneOffs",
			"Collapse any haplotypes that are low frequency compared to another haplotype (determined by lowFreqMultiplier) and only differs by 1 base", false, "Filtering");
	setUp.setOption(masterPopClusPars.lowFreqMultiplier, "--oneOffLowFreqMultiplier",
			"Low Freq Multiplier used for --excludeLowFreqOneOffs, considered low frequency if haplotype frac is less than its fraction times this number than the other haplotype", false, "Filtering");

	setUp.setOption(masterPopClusPars.lowLevelPopFiltPars_.removeOneSampOnlyOneOffHaps_, "--removeOneSampOnlyOneOffHaps",
			"Remove haplotypes that are below --oneSampOnlyOneOffHapsFrac fraction(default 0.20) that only appear in one sample that is one off of another haplotype within sample", false, "Filtering");
	setUp.setOption(masterPopClusPars.lowLevelPopFiltPars_.oneSampOnlyOneOffHapsFrac_, "--oneSampOnlyOneOffHapsFrac",
			"Fraction for --removeOneSampOnlyOneOffHaps", false, "Filtering");

	setUp.setOption(masterPopClusPars.lowLevelPopFiltPars_.removeOneSampOnlyHaps_, "--removeOneSampOnlyHaps",
			"Remove haplotypes that are below --OneSampOnlyHapsFrac fraction(default 0.20) that only appear in one sample that is one off of another haplotype within sample", false, "Filtering");
	setUp.setOption(masterPopClusPars.lowLevelPopFiltPars_.oneSampOnlyHapsFrac_, "--oneSampOnlyHapsFrac",
			"Fraction for --removeOneSampOnlyHaps", false, "Filtering");

	setUp.setOption(masterPopClusPars.rescuePars_.majorHaplotypeFracForRescue_, "--majorHaplotypeFracForRescue", "In order to be considered a major haplotype in a sample for comparing during rescue");
	setUp.setOption(masterPopClusPars.rescuePars_.rescueExcludedChimericHaplotypes, "--rescueExcludedChimericHaplotypes", "Rescue Excluded chimeric Haplotypes if they appear as a major haplotype in another sample");
	setUp.setOption(masterPopClusPars.rescuePars_.rescueExcludedOneOffLowFreqHaplotypes, "--rescueExcludedOneOffLowFreqHaplotypes", "Rescue Excluded one off low freq Haplotypes if they appear as a major haplotype in another sample");
	setUp.setOption(masterPopClusPars.rescuePars_.rescueExcludedLowFreqHaplotypes, "--rescueExcludedLowFreqHaplotypes", "Rescue Excluded Haplotypes for falling below frequency cut off if they appear as a major haplotype in another sample");
	setUp.setOption(masterPopClusPars.rescueMatchingExpected, "--rescueMatchingExpected", "Rescue Haplotypes that match expected sequences if they have been read using --ref");
	setUp.setOption(masterPopClusPars.fracCutoff, "--fracCutOff", "Final cluster Fraction Cut off", false, "Filtering");
	if(masterPopClusPars.withinReplicateFracCutOff > masterPopClusPars.fracCutoff){
		masterPopClusPars.withinReplicateFracCutOff = masterPopClusPars.fracCutoff;
	}
	setUp.setOption(masterPopClusPars.withinReplicateFracCutOff, "--withinReplicateFracCutOff", "Within Replicate Frac Cut Off, this is done before filtering sequences for appearing in all replicates", false, "Filtering");

	masterPopClusPars.differentPar = setUp.setOption(masterPopClusPars.parametersPopulation, "--popPar",
			"Parameters For Population Collapse", false, "Population");
	struct ClusteringParametersPars {
		bool noErrorsSet = false;
		bool strictErrorsSet = false;
		uint32_t stopAfter = 100;
		uint32_t hqMismatches = 0;
	};
	ClusteringParametersPars popClusParsPars;
	setUp.setOption(popClusParsPars.noErrorsSet, "--pop-noErrors", "Collapse parameters with no errors in population clustering", false, "Population");
	setUp.setOption(popClusParsPars.strictErrorsSet, "--pop-strictErrors", "Collapse parameters with a several low quality mismatches in population clustering", false, "Population");
	setUp.setOption(popClusParsPars.hqMismatches, "--pop-hq", "Number of high quality mismatches to allow in population clustering", false, "Population");
	setUp.setOption(popClusParsPars.stopAfter, "--pop-stopAfter", "Number of top haplotypes to check in population clustering", false, "Population");
	setUp.setOption(swgaPopClusErrorSet, "--pop-swgaErrorSet", "Collapse parameters for population clustering with SWGA data (collapsing on homopolymer indel)s", false, "Clustering");

	setUp.setOption(setUp.pars_.chiOpts_.checkChimeras_, "--recheckChimeras", "Check Input sequences for possible Chimeras", false, "Chimeras");
	setUp.setOption(masterPopClusPars.keepChimeras, "--keepChimeras", "KeepChimeras", false, "Chimeras");
	setUp.setOption(setUp.pars_.chiOpts_.parentFreqs_, "--parFreqs", "Chimeric Parent Frequency multiplier cutoff", false, "Chimeras");

	setUp.setOption(masterPopClusPars.numThreads, "--numThreads", "Number of threads to use");
	uint32_t numThreadsForSample = 1;
	setUp.setOption(numThreadsForSample, "--numThreadsForSample", "Number of threads to use for each population clustering");

	setUp.setOption(masterPopClusPars.writeOutAllInfoFile, "--writeOutAllInfoFile", "Write Out All Info File that contains information on all clusters including excluded ones");
	setUp.setOption(setUp.pars_.colOpts_.clusOpts_.converge_, "--converge", "Keep clustering at each iteration until there is no more collapsing, could increase run time significantly", false, "Clustering");
	//setUp.setOption(pars.plotRepAgreement, "--plotRepAgreement", "Plot Rep Agreement");
	if(swgaPopClusErrorSet || swgaSampleClusErrorSet){
		setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_ = true;
	}
	setUp.processAlignerDefualts();


	if (!setUp.failed_) {
		if ("" != masterPopClusPars.parameters) {
			// read in the parameters from the parameters file
			if (masterPopClusPars.onPerId) {
				masterPopClusPars.iteratorMap = setUp.processIteratorMapOnPerId(masterPopClusPars.parameters);
			} else {
				masterPopClusPars.iteratorMap = setUp.processIteratorMap(masterPopClusPars.parameters);
			}
		} else if(swgaSampleClusErrorSet){
			masterPopClusPars.iteratorMap = CollapseIterations::genIlluminaDefaultParsWithHqsCollapseHomopolymers(100,0);
		} else if (masterPopClusPars.noErrorsSet) {
			if (masterPopClusPars.hqMismatches > 0) {
				masterPopClusPars.iteratorMap =
						CollapseIterations::genStrictNoErrorsDefaultParsWithHqs(
								masterPopClusPars.stopAfter, masterPopClusPars.hqMismatches);
			} else {
				masterPopClusPars.iteratorMap = CollapseIterations::genStrictNoErrorsDefaultPars(
						masterPopClusPars.stopAfter);
			}
		} else if (masterPopClusPars.strictErrorsSet) {
			masterPopClusPars.iteratorMap = CollapseIterations::genStrictDefaultParsWithHqs(
					masterPopClusPars.stopAfter, masterPopClusPars.hqMismatches, masterPopClusPars.illumina);
		}

		if (masterPopClusPars.binParameters != "") {
			if (masterPopClusPars.onPerId) {
				masterPopClusPars.binIteratorMap = setUp.processIteratorMapOnPerId(masterPopClusPars.binParameters);
			} else {
				masterPopClusPars.binIteratorMap = setUp.processIteratorMap(masterPopClusPars.binParameters);
			}
		}
		if (masterPopClusPars.differentPar || popClusParsPars.noErrorsSet || popClusParsPars.strictErrorsSet || swgaPopClusErrorSet) {
			if (popClusParsPars.noErrorsSet) {
				if (popClusParsPars.hqMismatches > 0) {
					masterPopClusPars.popIteratorMap =
							CollapseIterations::genStrictNoErrorsDefaultParsWithHqs(
									popClusParsPars.stopAfter, popClusParsPars.hqMismatches);
				} else {
					masterPopClusPars.popIteratorMap = CollapseIterations::genStrictNoErrorsDefaultPars(
							popClusParsPars.stopAfter);
				}
			} else if(swgaPopClusErrorSet){
				masterPopClusPars.popIteratorMap = CollapseIterations::genIlluminaDefaultParsWithHqsCollapseHomopolymers(100,0);
			} else if (popClusParsPars.strictErrorsSet) {
				masterPopClusPars.popIteratorMap = CollapseIterations::genStrictDefaultParsWithHqs(
						popClusParsPars.stopAfter, popClusParsPars.hqMismatches, masterPopClusPars.illumina);
			} else {
				if (masterPopClusPars.onPerId) {
					masterPopClusPars.popIteratorMap = setUp.processIteratorMapOnPerId(
							masterPopClusPars.parametersPopulation);
				} else {
					masterPopClusPars.popIteratorMap = setUp.processIteratorMap(masterPopClusPars.parametersPopulation);
				}
			}
		} else {
			masterPopClusPars.popIteratorMap = masterPopClusPars.iteratorMap;
		}
	}
	setUp.finishSetUp(std::cout);

	setUp.startARunLog(setUp.pars_.directoryName_);
	std::vector<bfs::path> directories;
	if (samples.empty()) {
		auto allFiles = njh::files::listAllFiles(inputDirectory, false, {std::regex { ".*" + pat + "$" } });
		for (const auto &f : allFiles) {
			if (f.second) {
				directories.emplace_back(f.first);
			}
		}
	} else {
		for (const auto &samp : samples) {
			directories.emplace_back(njh::files::make_path(inputDirectory, samp + pat));
		}
	}

	if(directories.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error no directories found in " << inputDirectory << " ending with " << pat << "\n";
		throw std::runtime_error{ss.str()};
	}
	auto reportsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"reports"});
	auto infoDir =    njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar{"info"});
	std::unordered_set<std::string> allTargets;
	//gather all basic info
	{
		OutputStream allBasicInfo(njh::files::make_path(reportsDir, "allBasicInfo.tab.txt.gz"));
		std::shared_ptr<TableReader> firstTable;

		for(const auto & dir : directories){
			auto basicFnp = njh::files::make_path(dir, "final", "basicInfoPerRegion.tab.txt");
			if(bfs::exists(basicFnp)){
				TableIOOpts currentTabOpts(InOptions(basicFnp), "\t", true);
				firstTable = std::make_shared<TableReader>(currentTabOpts);
				if(!skipRBind){
					firstTable->header_.outPutContents(allBasicInfo, "\t");
				}
				VecStr row;
				while(firstTable->getNextRow(row)){
//					if(!skipRBind){
//						allBasicInfo << njh::conToStr(row, "\t") << "\n";
//					}
					allTargets.emplace(row[firstTable->header_.getColPos("name")]);
				}
				break;
			}
		}
		std::mutex allBasicInfoMut;
		if(!skipRBind){
			njh::concurrent::LockableVec<bfs::path> dirQueue(directories);

			std::function<void()> gatherAllBasicInfoFiles= [&allBasicInfo,&allBasicInfoMut,&dirQueue, &firstTable](){
				bfs::path dir;
				while(dirQueue.getVal(dir)){
					auto basicFnp = njh::files::make_path(dir, "final", "basicInfoPerRegion.tab.txt");
					std::stringstream tabOut;

					if(bfs::exists(basicFnp)){
						TableIOOpts currentTabOpts(InOptions(basicFnp), "\t", true);
						TableReader currentTable(currentTabOpts);
						if(currentTable.header_.columnNames_ != firstTable->header_.columnNames_){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "table " << basicFnp << " header's doesn't match other's header"<< "\n";
							ss << "Expected Header: " << njh::conToStr(firstTable->header_.columnNames_, ",") << "\n";
							ss << "Current  Header: " << njh::conToStr(currentTable.header_.columnNames_, ",") << "\n";
							throw std::runtime_error{ss.str()};
						}
						VecStr row;
						while(currentTable.getNextRow(row)){
							tabOut << njh::conToStr(row, "\t") << "\n";
						}
					}
					{
						std::lock_guard<std::mutex> lock(allBasicInfoMut);
						allBasicInfo << tabOut.str();
					}
				}
			};
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			njh::concurrent::runVoidFunctionThreaded(gatherAllBasicInfoFiles, masterPopClusPars.numThreads);
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(skipRBind){
		bfs::remove(njh::files::make_path(reportsDir, "allBasicInfo.tab.txt.gz"));
	}


	std::shared_ptr<MultipleGroupMetaData> meta;
	if("" != masterPopClusPars.groupingsFile){
		meta = std::make_shared<MultipleGroupMetaData>(masterPopClusPars.groupingsFile);
	}




	std::string sampleField = "sample";
	std::string targetField = "regionUID";

	auto sampleMetaDataFnp = njh::files::make_path(infoDir, "sampleMetaData.tab.txt");
	njh::stopWatch watch;
	watch.setLapName("Gathering Raw Seqs");
	auto rawAllSeqsFnp = njh::files::make_path(infoDir, "rawAllSeqsFile.fasta");

	//
	SeqGatheringFromPathWeaver::SeqGatheringFromPathWeaverCorePars gatherCorePars;
	SeqGatheringFromPathWeaver::gatherSeqsAndSortByTargetPars rawGatherPars;
	SeqGatheringFromPathWeaver::processedGatherSeqsMetaPars processGatherPars;

	gatherCorePars.sampleField = sampleField;
	gatherCorePars.targetField = targetField;
	gatherCorePars.countField = countField;
	gatherCorePars.meta = meta;
	gatherCorePars.numThreads = masterPopClusPars.numThreads;

	rawGatherPars.addPartial = addPartial;
	rawGatherPars.allSeqFnp = rawAllSeqsFnp;
	rawGatherPars.directories = directories;
	rawGatherPars.targets = targets;

	if(trimBedSet){
		TwoBit::TwoBitFile tReader(genomeFnp);
		auto locations = bedPtrsToGenomicRegs(getBeds(trimBedFnp));
		for(const auto & loc : locations){
			rawGatherPars.trimSeqs[loc.uid_].emplace_back(seqWithKmerInfo(loc.extractSeq(tReader),7, false));
		}
	}


	SeqGatheringFromPathWeaver seqGatherer(gatherCorePars);
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	auto rawGatherRes = seqGatherer.gatherSeqsAndSortByTarget(rawGatherPars);
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(!rawGatherRes.missingOutput.empty()){
		OutputStream missingOut(njh::files::make_path(reportsDir, "missingDataForSamples.txt"));
		missingOut << njh::conToStr(rawGatherRes.missingOutput, "\n") << std::endl;
	}

	processGatherPars.keepCommonSeqsWhenFiltering = keepCommonSeqsWhenFiltering;
	processGatherPars.allSeqFnp = rawGatherRes.allSeqFnp;
	processGatherPars.outMetaFnp = sampleMetaDataFnp;

	watch.startNewLap("Filtering Seqs On Meta");


	auto processedGatherRes = seqGatherer.processedGatherSeqsMeta(processGatherPars, rawGatherRes);

	{
		OutputStream targetsGatheredReport(njh::files::make_path(reportsDir, "seqsPerTargetGathered.txt"));
		auto targetKeys = getVectorOfMapKeys(processedGatherRes.seqsLocations);
		njh::sort(targetKeys);
		targetsGatheredReport << "target\tinputSeqs" << std::endl;
		for(const auto & key : targetKeys){
			targetsGatheredReport
			<< key << "\t" << processedGatherRes.seqsLocations[key].second
			<< std::endl;
		}
		VecStr noDataForTargets;
		for(const auto & tar : allTargets){
			auto rawTargetName = rawGatherRes.targetKey[tar];
			if(!njh::in(rawTargetName, targetKeys) && (targets.empty() || njh::in(rawTargetName, targets))){
				noDataForTargets.emplace_back(tar);
			}
		}
		njh::sort(noDataForTargets);
		for(const auto & tar : noDataForTargets){
			targetsGatheredReport << tar << "\t" << 0 << std::endl;
		}
	}
	if(bfs::exists(processGatherPars.outMetaFnp)){
		masterPopClusPars.groupingsFile = processGatherPars.outMetaFnp.string();
	}

	auto allSeqsFnp = processedGatherRes.allSeqFnp;
	auto allSeqsFnpIn = SeqIOOptions::genFastaIn(allSeqsFnp);
	SeqInput::buildIndex(allSeqsFnpIn);
	allSeqsFnpIn.processed_ = true;

	{
		OutputStream outKey(njh::files::make_path(infoDir, "targetToPNameKey.tab.txt"));
		outKey << "p_name\t" << targetField << std::endl;
		auto keys = getVectorOfMapKeys(rawGatherRes.targetKey);
		njh::sort(keys);
		for(const auto & key : keys){
			outKey << key << '\t' << rawGatherRes.targetKey[key] << std::endl;
		}
	}
	{
		auto targetsLocsFnp = njh::files::make_path(infoDir, "targetLocations.tab.txt");

		OutputStream targetsLocsOut(targetsLocsFnp);
		VecStr tarKeys = getVectorOfMapKeys(rawGatherRes.seqsLocations);
		njh::sort(tarKeys);
		for(const auto & tar : tarKeys){
			targetsLocsOut << tar
					<< "\t" << rawGatherRes.seqsLocations[tar].first
					<< "\t" << rawGatherRes.seqsLocations[tar].second
					<< std::endl;
		}
	}
	if(processedGatherRes.allSeqFnp != rawGatherRes.allSeqFnp){
		VecStr tarKeys = getVectorOfMapKeys(processedGatherRes.seqsLocations);
		njh::sort(tarKeys);
		OutputStream processedTargetsLocsOut(njh::files::make_path(infoDir, "processedTargetLocations.tab.txt"));

		for(const auto & tar : tarKeys){
			processedTargetsLocsOut << tar
					<< "\t" << rawGatherRes.seqsLocations[tar].first
					<< "\t" << rawGatherRes.seqsLocations[tar].second
					<< std::endl;
		}
	}

	if(!processedGatherRes.failedToCombine.empty()){
		OutputStream failedToCombineOut(njh::files::make_path(infoDir, "failedToCombine.txt"));
		failedToCombineOut << "target\tmetaField\tsamples" << std::endl;
		for(const auto & target : iter::sorted(njh::getVecOfMapKeys(processedGatherRes.failedToCombine))){
			for(const auto & field : processedGatherRes.failedToCombine.at(target)){
				failedToCombineOut << target
						<< "\t" << field.first
						<< "\t" << njh::conToStr(field.second, ",") << std::endl;
			}
		}
	}

	if(setUp.pars_.verbose_){
		watch.logLapTimes(std::cout, true, 6, true);
		watch.setLapName("clustering");
	}


	SeqIOOptions seqOpts;
	seqOpts.inFormat_ = SeqIOOptions::inFormats::FASTA;
	seqOpts.outFormat_ = SeqIOOptions::outFormats::FASTAGZ;

	auto tarNames = getVectorOfMapKeys(processedGatherRes.seqsLocations);
	njh::sort(tarNames);
	njh::concurrent::LockableQueue<std::string> tarNamesQueue(tarNames);
	const auto masterPopClusParsCopyConst = masterPopClusPars;
	std::function<void()> runPopClusOnTar = [&tarNamesQueue,&masterPopClusParsCopyConst,
																					 &setUp,&processedGatherRes,
																					 &allSeqsFnpIn,
																					 &seqOpts,&numThreadsForSample](){
		std::string tar = "";

	while(tarNamesQueue.getVal(tar)){
		uint64_t maxLen = 0;

		//read in seqs for target
		std::unordered_map<std::string, std::unordered_map<std::string, njhseq::collapse::SampleCollapseCollection::RepSeqs>> seqsForSample;
		{
			SeqInput targetSeqReader(allSeqsFnpIn);
			auto seqs = targetSeqReader.getReads<seqInfo>(processedGatherRes.seqsLocations.at(tar).first, processedGatherRes.seqsLocations.at(tar).second);
//			std::cout << "target: " << tar << std::endl;
//			std::cout << '\t' << seqs.size() << std::endl;

			for(const auto & seq : seqs){
//				std::cout << "\t" << seq.name_ << std::endl;
				readVec::getMaxLength(seq, maxLen);
				MetaDataInName meta(seq.name_);
				auto sample = meta.getMeta("sample");
				if(njh::in(sample, seqsForSample)){
					seqsForSample[sample].at(sample).repSeqs_.emplace_back(seq);
				}else{
					seqsForSample[sample].emplace(sample, njhseq::collapse::SampleCollapseCollection::RepSeqs(sample, std::vector<seqInfo>{seq}));
				}
			}
		}

		//translation helpers
		std::unique_ptr<TranslatorByAlignment> translator;
		if("" != masterPopClusParsCopyConst.transPars.gffFnp_){
			translator = std::make_unique<TranslatorByAlignment>(masterPopClusParsCopyConst.transPars);
		}

		auto currentPars = masterPopClusParsCopyConst;
		currentPars.experimentNames.populationName_ = tar;
		//make directory
		auto currentMasterDir = njh::files::make_path(setUp.pars_.directoryName_, tar);
		njh::files::makeDir(njh::files::MkdirPar{currentMasterDir});
		currentPars.masterDir = currentMasterDir;
		//population seqs;
//		std::vector<seqInfo> globalPopSeqs;
//		if("" != pars.popSeqsFnp){
//			globalPopSeqs = SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastaIn(pars.popSeqsFnp));
//		}
//
		// reading expected sequences to compare to
		//bool checkingExpected = setUp.pars_.refIoOptions_.firstName_ != "";
//		std::vector<readObject> expectedSeqs;
//		if (checkingExpected) {
//			expectedSeqs = SeqInput::getReferenceSeq(setUp.pars_.refIoOptions_, maxLen);
//		}

//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		// create aligner class object
		aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
				KmerMaps(setUp.pars_.colOpts_.kmerOpts_.kLength_),
				setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
				setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
		alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
		njhseq::concurrent::AlignerPool alnPool(alignerObj,numThreadsForSample );
		alnPool.initAligners();
		alnPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		// create collapserObj used for clustering
		collapser collapserObj(setUp.pars_.colOpts_);
		collapserObj.opts_.kmerOpts_.checkKmers_ = false;
		// output info about the read In reads
		collapse::SampleCollapseCollection sampColl(seqOpts, currentPars.masterDir,
				currentMasterDir,
				PopNamesInfo(currentPars.experimentNames.populationName_, processedGatherRes.allSamples, masterPopClusParsCopyConst.experimentNames.controlSamples_),
				currentPars.preFiltCutOffs);

		sampColl.keepSampleInfoInMemory_ = true;

		if("" != currentPars.groupingsFile){
			sampColl.addGroupMetaData(currentPars.groupingsFile);
		}
		//process custom cut offs
		std::unordered_map<std::string, double> customCutOffsMap = collapse::SampleCollapseCollection::processCustomCutOffs(currentPars.customCutOffs, VecStr{processedGatherRes.allSamples.begin(), processedGatherRes.allSamples.end()}, currentPars.fracCutoff);
		std::unordered_map<std::string, double> customCutOffsMapPerRep = collapse::SampleCollapseCollection::processCustomCutOffs(currentPars.customCutOffs, VecStr{processedGatherRes.allSamples.begin(), processedGatherRes.allSamples.end()}, currentPars.withinReplicateFracCutOff);
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		{
			njh::concurrent::LockableQueue<std::string> sampleQueue(processedGatherRes.allSamples);
//			std::cout << njh::conToStr(processedGatherRes.allSamples, ",") << std::endl;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			std::function<void()> setupClusterSamples = [
																									 &sampleQueue,&alnPool,&collapserObj,&currentPars,&setUp,
																	&sampColl,&customCutOffsMap,
																	&customCutOffsMapPerRep, &seqsForSample](){
				std::string samp = "";
				auto currentAligner = alnPool.popAligner();
				while(sampleQueue.getVal(samp)){
					if(setUp.pars_.verbose_){
						std::cout << "Starting: " << samp << std::endl;
					}
					if(!njh::in(samp, seqsForSample)){
						sampColl.lowRepCntSamples_.emplace_back(samp);
						//doesn't contain this sample, continue onwards
						continue;
					}

					sampColl.setUpSample(samp,seqsForSample.at(samp), *currentAligner, collapserObj, setUp.pars_.chiOpts_);

					sampColl.clusterSample(samp, *currentAligner, collapserObj, currentPars.iteratorMap);

					sampColl.sampleCollapses_.at(samp)->markChimeras(currentPars.chiCutOff);

					//exclude clusters that don't have the necessary replicate number
					//defaults to the number of input replicates if none supplied

					if (0 != currentPars.runsRequired) {
						sampColl.sampleCollapses_.at(samp)->excludeBySampNum(currentPars.runsRequired, true);

					} else {
						sampColl.sampleCollapses_.at(samp)->excludeBySampNum(sampColl.sampleCollapses_.at(samp)->input_.info_.infos_.size(), true);
					}

					if(currentPars.collapseLowFreqOneOffs){
						bool skipChimeras = true;
						sampColl.sampleCollapses_.at(samp)->excludeLowFreqOneOffs(true, currentPars.lowFreqMultiplier, *currentAligner,
								skipChimeras, customCutOffsMap.at(samp));
					}

					sampColl.sampleCollapses_.at(samp)->excludeFractionAnyRep(customCutOffsMapPerRep.at(samp), true);
					sampColl.sampleCollapses_.at(samp)->excludeFraction(customCutOffsMap.at(samp), true);

					if (!currentPars.keepChimeras) {
						//now exclude all marked chimeras
						sampColl.sampleCollapses_.at(samp)->excludeChimerasNoReMark(true);
					}


					std::string sortBy = "fraction";
					sampColl.sampleCollapses_.at(samp)->renameClusters(sortBy);


//					if (!expectedSeqs.empty()) {
//						sampColl.sampleCollapses_.at(samp)->excluded_.checkAgainstExpected(
//								expectedSeqs, *currentAligner, false);
//						sampColl.sampleCollapses_.at(samp)->collapsed_.checkAgainstExpected(
//								expectedSeqs, *currentAligner, false);
//						if(setUp.pars_.debug_){
//							std::cout << "sample: " << samp << std::endl;
//						}
//						for(const auto & clus : sampColl.sampleCollapses_.at(samp)->collapsed_.clusters_){
//							if(setUp.pars_.debug_){
//								std::cout << clus.seqBase_.name_ << " : " << clus.expectsString << std::endl;
//							}
//							if("" ==  clus.expectsString ){
//								std::stringstream ss;
//								ss << __PRETTY_FUNCTION__ << ": Error, expects string is blank" << std::endl;
//								ss << clus.seqBase_.name_ << std::endl;
//								throw std::runtime_error{ss.str()};
//							}
//						}
//					}

					if(!sampColl.keepSampleInfoInMemory_){
						sampColl.dumpSample(samp);
					}


					if(setUp.pars_.verbose_){
						std::cout << "Ending: " << samp << std::endl;
					}
				}
			};
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			njh::concurrent::runVoidFunctionThreaded(setupClusterSamples,numThreadsForSample );
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
		if(setUp.pars_.verbose_){
			std::cout << njh::bashCT::boldGreen("Pop Clustering") << std::endl;
		}


//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//first population clustering createSharedPathwaysFromReads
//		auto popInput = sampColl.createPopInput();
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "popInput.size() : "<< popInput.size() << std::endl;
		{

			auto populationInput = sampColl.createPopInput();
			//std::cout << populationInput.size() << std::endl;
//			std::cout <<njh::bashCT::boldRed("Sleeping......") << std::endl;;
//			using namespace std::chrono_literals;
//			std::this_thread::sleep_for(100000s);
			sampColl.doPopulationClustering(populationInput, alignerObj, collapserObj, currentPars.popIteratorMap);

		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if (setUp.pars_.debug_) {
			std::cout << "sampColl.popCollapse_->input_.info_.numberOfClusters_: " << sampColl.popCollapse_->input_.info_.numberOfClusters_ << std::endl;
			std::cout << "sampColl.popCollapse_->input_.info_.totalReadCount_: " << sampColl.popCollapse_->input_.info_.totalReadCount_ << std::endl;
		}
		if(setUp.pars_.verbose_){
			std::cout << njh::bashCT::boldGreen("Done Initial Population Clustering") << std::endl;
		}
		if(setUp.pars_.verbose_){

		}


		if(currentPars.rescuePars_.performResuce()){
			sampColl.conductResuceOperations(currentPars.rescuePars_, alignerObj, collapserObj, currentPars.popIteratorMap);
		}
		sampColl.performLowLevelFilters(currentPars.lowLevelPopFiltPars_, alignerObj, collapserObj, currentPars.popIteratorMap);

//		if(currentPars.rescueMatchingExpected && !expectedSeqs.empty()){
//			sampColl.rescueMatchingSeqs(expectedSeqs, alignerObj, collapserObj, currentPars.popIteratorMap);
//		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(setUp.pars_.verbose_){
			std::cout << njh::bashCT::boldRed("Done Pop Clustering") << std::endl;
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//if ("" != currentPars.previousPopFilename && !currentPars.noPopulation) {
		if ("" != currentPars.previousPopFilename) {
			auto previousPopSeqsRaw = getSeqs<readObject>(currentPars.previousPopFilename);
			//collapse indentical seqs
			std::vector<readObject> previousPopSeqs;
			std::vector<std::set<std::string>> allNamesForPreviousPops;
			for(const auto & seq : previousPopSeqsRaw){
				bool matchPrevious = false;
				for(const auto  other : iter::enumerate(previousPopSeqs)){
					if(other.element.seqBase_.seq_ == seq.seqBase_.seq_){
						matchPrevious = true;
						allNamesForPreviousPops[other.first].emplace(seq.seqBase_.name_);
						break;
					}
				}
				if(!matchPrevious){
					previousPopSeqs.emplace_back(seq);
					allNamesForPreviousPops.emplace_back(std::set<std::string>{seq.seqBase_.name_});
				}
			}
			//rename
			for(const auto  other : iter::enumerate(previousPopSeqs)){
				if(allNamesForPreviousPops[other.index].size() > 1){
					other.element.seqBase_.name_ = njh::conToStr(allNamesForPreviousPops[other.index], ":");
				}
			}
			sampColl.renamePopWithSeqs(previousPopSeqs, currentPars.previousPopErrors);
		}

//		if (!expectedSeqs.empty()) {
//			sampColl.comparePopToRefSeqs(expectedSeqs, alignerObj);
//		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		std::vector<seqInfo> popSeqsPerSamp;
		std::vector<seqInfo> outPopSeqsPerSamp;

		std::unordered_map<std::string, uint32_t> sampCountsForPopHaps;
		std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;

		uint32_t totalPopCount = 0;
		std::set<std::string> samplesCount;
//		std::cout << njh::conToStr(processedGatherRes.allSamples) << std::endl;
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		popSeqsPerSamp = sampColl.genOutPopSeqsPerSample();
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		outPopSeqsPerSamp = popSeqsPerSamp;
		for(const auto & popClus : sampColl.popCollapse_->collapsed_.clusters_){
			sampCountsForPopHaps[popClus.seqBase_.name_] = popClus.sampleClusters().size();
			auto samples = getVectorOfMapKeys(popClus.sampleClusters());
			sampNamesForPopHaps[popClus.seqBase_.name_].insert(samples.begin(), samples.end());
			totalPopCount += popClus.sampleClusters().size();
		}
		for(const auto & seq : popSeqsPerSamp){
			MetaDataInName seqMeta(seq.name_);
			samplesCount.emplace(seqMeta.getMeta("sample"));
		}

//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		std::map<std::string, std::map<std::string, MetaDataInName>> knownAAMeta;
		//       seqName               transcript   amino acid positions and amino acid
		std::map<std::string, std::map<std::string, std::string>> knownAATyped;
		std::map<std::string, std::map<std::string, std::string>> fullAATyped;
		if("" != currentPars.transPars.gffFnp_){

			auto variantInfoDir =  njh::files::make_path(sampColl.masterOutputDir_, "variantInfo");
			njh::files::makeDir(njh::files::MkdirPar{variantInfoDir});
			translator->pars_.keepTemporaryFiles_ = true;
			translator->pars_.workingDirtory_ = variantInfoDir;
			auto popSeqsOpts = seqOpts;
			SeqIOOptions popOutOpts(
					njh::files::make_path(variantInfoDir,
							"PopSeqs" + seqOpts.getOutExtension()),
							seqOpts.outFormat_);
			SeqOutput::write(sampColl.popCollapse_->collapsed_.clusters_, popOutOpts);

			popSeqsOpts.firstName_ = popOutOpts.out_.outName();
			auto translatedRes = translator->run(popSeqsOpts, sampNamesForPopHaps, currentPars.variantCallerRunPars);
			SeqOutput transwriter(SeqIOOptions::genFastaOut(njh::files::make_path(variantInfoDir, "translatedInput.fasta")));
			for(const auto & seqName : translatedRes.translations_){
				for(const auto & transcript : seqName.second){
					transwriter.openWrite(transcript.second.translation_);
				}
			}
			SeqInput popReader(popSeqsOpts);
			auto popSeqs = popReader.readAllReads<seqInfo>();
			std::unordered_map<std::string, uint32_t> popSeqsPosition;
			for(const auto popPos : iter::range(popSeqs.size())){
				popSeqsPosition[popSeqs[popPos].name_] = popPos;
			}
			OutputStream popBedLocs(njh::files::make_path(variantInfoDir, "PopSeqs.bed"));
			for(const auto & seqLocs : translatedRes.seqAlns_){
				for(const auto & loc : seqLocs.second){
					popBedLocs << loc.gRegion_.genBedRecordCore().toDelimStrWithExtra() << std::endl;
				}
			}

			for(const auto & pop : popSeqs){
				if(!njh::in(pop.name_, translatedRes.seqAlns_)){
					totalPopCount -= sampCountsForPopHaps[pop.name_];
					popBedLocs << "*"
							<< "\t" << "*"
							<< "\t" << "*"
							<< "\t" << pop.name_
							<< "\t" << "*"
							<< "\t" << "*" << std::endl;
				}
			}

			{
				//protein
				for(auto & varPerTrans : translatedRes.proteinVariants_){

					varPerTrans.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidVariable.tab.txt")), varPerTrans.first, true	);
					varPerTrans.second.writeOutSNPsAllInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidsAll.tab.txt")), varPerTrans.first, true	);

					std::set<uint32_t> knownMutationsLocations;
					for(const auto & snpPos : varPerTrans.second.allBases){
						if(njh::in(snpPos.first + 1, translator->knownAminoAcidPositions_[varPerTrans.first])){
							knownMutationsLocations.emplace(snpPos.first);
						}
					}
					if(!varPerTrans.second.variablePositons_.empty()){
						GenomicRegion variableRegion = varPerTrans.second.getVariableRegion();
						variableRegion.start_ = variableRegion.start_ +1;
						OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_variableRegion.bed"))));
						bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
					}

					std::set<uint32_t> allLocations(knownMutationsLocations.begin(), knownMutationsLocations.end());
					for(const auto & variablePos : varPerTrans.second.snpsFinal){
						allLocations.emplace(variablePos.first);
					}
					std::map<std::string, MetaDataInName> aaMeta;


					for(auto & seqName : translatedRes.translations_){
						if(njh::in(varPerTrans.first, seqName.second)){
							for(const auto & variablePos : allLocations){
								auto aa = seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_, variablePos)];
								aaMeta[seqName.first].addMeta(estd::to_string(variablePos), aa, false);
							}
						}
						for(const auto & knownLoc : knownMutationsLocations){
							auto aa = seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_, knownLoc)];
							knownAAMeta[seqName.first.substr(0, seqName.first.rfind("_f"))][varPerTrans.first].addMeta(estd::to_string(knownLoc), aa, false);
						}
					}

					for (auto & seqName : translatedRes.translations_) {
						if (njh::in(varPerTrans.first, seqName.second)) {
							VecStr allAAPosCoded;
							std::string popName = seqName.first.substr(0, seqName.first.rfind("_f"));
							std::string transcript = varPerTrans.first;
							for (const auto & loc : allLocations) {
								//std::cout << __FILE__ << " " << __LINE__ << std::endl;
								auto aa =seqName.second[varPerTrans.first].queryAlnTranslation_.seq_[getAlnPosForRealPos(seqName.second[varPerTrans.first].refAlnTranslation_.seq_,loc)];
								allAAPosCoded.emplace_back(njh::pasteAsStr(loc + 1, "-", aa));
								//std::cout << __FILE__ << " " << __LINE__ << std::endl;
							}
							if(!allAAPosCoded.empty()){
								fullAATyped[popName][transcript] = njh::conToStr(allAAPosCoded, ":");
							}else{
								fullAATyped[popName][transcript] = "NONE";
							}
						}
					}
					if(!knownMutationsLocations.empty()){
						varPerTrans.second.writeOutSNPsInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-protein_aminoAcidKnownMutations.tab.txt")), varPerTrans.first,knownMutationsLocations, true	);

						for(const auto & pop : knownAAMeta){
							std::string popName = pop.first;
							std::string transcript = varPerTrans.first;
							std::vector<std::string> aaPos;
							for(const auto & variablePos : knownMutationsLocations){
								aaPos.emplace_back(estd::to_string(variablePos + 1) + "-" + pop.second.at(transcript).getMeta(estd::to_string(variablePos)));
							}
							knownAATyped[popName][transcript] = njh::conToStr(aaPos, ":");
						}
					}
					OutputStream outPopHapAminos(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerTrans.first +  "-popHapToAmino.tab.txt")));
					outPopHapAminos << "h_PopUID" ;
					VecStr aminoPositionsHeader;
					for(const auto & variablePos : allLocations){
						outPopHapAminos << "\t" << varPerTrans.first << "-aa" << variablePos + 1;
						aminoPositionsHeader.emplace_back(njh::pasteAsStr(varPerTrans.first, "-aa",  variablePos + 1));
					}
					outPopHapAminos << std::endl;
					std::unordered_map<std::string, std::string> popHapAminoTyped;
					std::unordered_map<std::string, VecStr> popHapAminoTypedRow;

					for(const auto & pop : aaMeta){
						outPopHapAminos << pop.first.substr(0, pop.first.rfind("_f")) ;
						std::string typed = varPerTrans.first + "=";
						std::vector<std::string> aaPos;
						std::vector<std::string> aa;

						for(const auto & variablePos : allLocations){
							outPopHapAminos << "\t" << pop.second.getMeta(estd::to_string(variablePos));
							aaPos.emplace_back(estd::to_string(variablePos + 1) + "-" + pop.second.getMeta(estd::to_string(variablePos)));
							aa.emplace_back(pop.second.getMeta(estd::to_string(variablePos)));

						}
						if (!aaPos.empty()) {
							typed += njh::conToStr(aaPos, ":");
						} else {
							typed += "NONE";
						}
						//typed +=";";

						popHapAminoTyped[pop.first.substr(0, pop.first.rfind("_f"))] = typed;
						popHapAminoTypedRow[pop.first.substr(0, pop.first.rfind("_f"))] = aa;
						outPopHapAminos << std::endl;
					}
					auto popMetaTable = seqsToMetaTable(popSeqsPerSamp);
					popMetaTable.deleteColumn("seq");
					popMetaTable.deleteColumn("count");
					addOtherVec(popMetaTable.columnNames_, aminoPositionsHeader);
					for(auto & row : popMetaTable){
						MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
						auto popName = row[popMetaTable.getColPos("PopUID")];
						VecStr aminos;
						if(njh::in(popName, popHapAminoTypedRow)){
							aminos = popHapAminoTypedRow[popName];
						}else{
							//pop uid was untypable, didn't align
							aminos = VecStr{aminoPositionsHeader.size(), std::string("NA")};
						}
						addOtherVec(row, aminos);
					}

					popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(variantInfoDir, varPerTrans.first + "-popSeqsWithMetaAndVariableAAInfoTable.tab.txt")));

	//				for(auto & popHapSamp : outPopSeqsPerSamp){
	//					MetaDataInName meta(popHapSamp.name_);
	//					auto typed = popHapAminoTyped[meta.getMeta("PopUID")].substr(popHapAminoTyped[meta.getMeta("PopUID")].find("=") + 1);
	//					if("" ==typed){
	//						typed = "NONE";
	//					}
	//					meta.addMeta(varPerTrans.first, typed);
	//					meta.resetMetaInName(popHapSamp.name_);
	//				}
				}
			}

			{
				//snps
				OutputStream outSnpDepthPerSample(njh::files::make_path(variantInfoDir, njh::pasteAsStr("snpDepthPerSample.tab.txt")));
				outSnpDepthPerSample << "AnalysisName\tsample\tchromosome\tposition\trefBase\tbase\treadDepth" ;
				VecStr metaLevels;
				if(nullptr != sampColl.groupMetaData_){
					metaLevels = getVectorOfMapKeys(sampColl.groupMetaData_->groupData_);
					for(const auto & meta : metaLevels){
						outSnpDepthPerSample << "\t" << meta;
					}
				}
				outSnpDepthPerSample << std::endl;

				for( auto & varPerChrom : translatedRes.seqVariants_){
					varPerChrom.second.writeOutSNPsFinalInfo(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-SNPs.tab.txt")), varPerChrom.first, false);
					varPerChrom.second.writeOutSNPsAllInfo(  njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-allBases.tab.txt")), varPerChrom.first, false);
					if(!varPerChrom.second.variablePositons_.empty()){
						GenomicRegion variableRegion = varPerChrom.second.getVariableRegion();
						OutputStream bedVariableRegionOut(OutOptions(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-chromosome_variableRegion.bed"))));
						bedVariableRegionOut << variableRegion.genBedRecordCore().toDelimStrWithExtra() << std::endl;
					}
					std::map<std::string, MetaDataInName> snpMeta;
					for(auto & seqName : translatedRes.seqAlns_){
						for(const auto & variablePos : varPerChrom.second.snpsFinal){

							if(variablePos.first < seqName.second.front().gRegion_.start_ || variablePos.first >= seqName.second.front().gRegion_.end_){
								snpMeta[seqName.first].addMeta(estd::to_string(variablePos.first), "X", false);
							} else {
								auto aa = seqName.second.front().alnQuerySeq_.seq_[getAlnPosForRealPos(seqName.second.front().alnRefSeq_.seq_, variablePos.first - seqName.second.front().gRegion_.start_)];
								snpMeta[seqName.first].addMeta(estd::to_string(variablePos.first), aa, false);
							}
						}
					}
	//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
					OutputStream outPopHapAminos(njh::files::make_path(variantInfoDir, njh::pasteAsStr(varPerChrom.first +  "-popHapToSNPs.tab.txt")));
					outPopHapAminos << "h_PopUID" ;
					VecStr aminoPositionsHeader;
					for(const auto & variablePos : varPerChrom.second.snpsFinal){
						outPopHapAminos << "\t" << varPerChrom.first << "-" << variablePos.first;
						aminoPositionsHeader.emplace_back(njh::pasteAsStr(varPerChrom.first, "-",  variablePos.first));
					}
					outPopHapAminos << std::endl;
					std::unordered_map<std::string, std::string> popHapAminoTyped;
					std::unordered_map<std::string, VecStr> popHapAminoTypedRow;

					for(const auto & pop : snpMeta){
						outPopHapAminos << pop.first.substr(0, pop.first.rfind("_f")) ;
						std::string typed = varPerChrom.first + "=";
						std::vector<std::string> aaPos;
						std::vector<std::string> aa;

						for(const auto & variablePos : varPerChrom.second.snpsFinal){
							outPopHapAminos << "\t" << pop.second.getMeta(estd::to_string(variablePos.first));
							aaPos.emplace_back(estd::to_string(variablePos.first) + "-" + pop.second.getMeta(estd::to_string(variablePos.first)));
							aa.emplace_back(pop.second.getMeta(estd::to_string(variablePos.first)));
						}
						if (!aaPos.empty()) {
							typed += njh::conToStr(aaPos, ":");
						} else {
							typed += "NONE";
						}
						//typed +=";";
						popHapAminoTyped[pop.first.substr(0, pop.first.rfind("_f"))] = typed;
						popHapAminoTypedRow[pop.first.substr(0, pop.first.rfind("_f"))] = aa;
						outPopHapAminos << std::endl;
					}


	//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
					auto popMetaTable = seqsToMetaTable(popSeqsPerSamp);
					popMetaTable.deleteColumn("seq");
					popMetaTable.deleteColumn("count");
					addOtherVec(popMetaTable.columnNames_, aminoPositionsHeader);

					std::map<std::string, std::map<std::string, std::map<uint32_t, std::map<std::string, double>>>> snpsPerSampleDepth;

					for(auto & row : popMetaTable){
						MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
						auto popName = row[popMetaTable.getColPos("PopUID")];
						VecStr aminos;
						if(njh::in(popName, popHapAminoTypedRow)){
							aminos = popHapAminoTypedRow[popName];
							uint32_t snpPosition = 0;
							for(const auto & snp : aminos){
								if("X" == snp){
									continue;
								}
								if("" == snp){
									continue;
								}
								if(" " == snp){
									continue;
								}
								auto sample = row[popMetaTable.getColPos("sample")];
								auto chrom = aminoPositionsHeader[snpPosition].substr(0, aminoPositionsHeader[snpPosition].rfind("-"));
	//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//							std::cout << "snpPosition: " << snpPosition << std::endl;
	//							std::cout << "aminoPositionsHeader[snpPosition]: " <<  aminoPositionsHeader[snpPosition] << std::endl;
	//							std::cout << "aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind(\"-\") + 1): " << aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind("-") + 1)<< std::endl;
	//							std::cout << "row[popMetaTable.getColPos(\"readCount\")]: " << row[popMetaTable.getColPos("readCount")] << std::endl;
								//std::cout << aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind("-") + 1) << std::endl;
								auto position = njh::StrToNumConverter::stoToNum<uint32_t>(aminoPositionsHeader[snpPosition].substr(aminoPositionsHeader[snpPosition].rfind("-") + 1));
								double depth = njh::StrToNumConverter::stoToNum<double>(row[popMetaTable.getColPos("readCount")]);
	//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
								snpsPerSampleDepth[sample][chrom][position][snp] += depth;
								++snpPosition;
							}
						} else {
							//pop uid was untypable, didn't align
							aminos = VecStr{aminoPositionsHeader.size(), std::string("NA")};
						}
						addOtherVec(row, aminos);
					}
					//std::cout << __FILE__ << " " << __LINE__ << std::endl;
					popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(variantInfoDir, varPerChrom.first + "-popSeqsWithMetaAndVariableSNPInfoTable.tab.txt")));


					for(const auto & sample : snpsPerSampleDepth){
						for(const auto & chrom : sample.second){
							for(const auto & position : chrom.second){
								for(const auto & snp : position.second){
									outSnpDepthPerSample
											<< currentPars.experimentNames.populationName_
											<< "\t" << sample.first
											<< "\t" << chrom.first
											<< "\t" << position.first
											<< "\t" << varPerChrom.second.getBaseForGenomicRegion(position.first) //translatedRes.baseForPosition_[varPerChrom.first][position.first]
											<< "\t" << snp.first
											<< "\t" << snp.second ;

									if (nullptr != sampColl.groupMetaData_) {
										auto metaForSample = sampColl.groupMetaData_->getMetaForSample(sample.first, metaLevels);
										for (const auto & meta : metaLevels) {
											outSnpDepthPerSample << "\t" << metaForSample.getMeta(meta);
										}
									}
									outSnpDepthPerSample << std::endl;
								}
							}
						}
					}
	//				for(auto & popHapSamp : outPopSeqsPerSamp){
	////					std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//					MetaDataInName meta(popHapSamp.name_);
	//					auto typed = popHapAminoTyped[meta.getMeta("PopUID")].substr(popHapAminoTyped[meta.getMeta("PopUID")].find("=") + 1);
	//					if("" ==typed){
	//						typed = "NONE";
	//					}
	//					meta.addMeta(varPerChrom.first, typed);
	//					meta.resetMetaInName(popHapSamp.name_);
	//				}
				}
			}
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for(auto & clus : sampColl.popCollapse_->collapsed_.clusters_){
			auto popName = clus.seqBase_.name_.substr(0, clus.seqBase_.name_.rfind("_f"));
			std::string typed = "";
	//		std::cout << "clus.seqBase_.name_: " << clus.seqBase_.name_ << std::endl;
	//		std::cout << "clus.meta_.toJson(): " << clus.meta_.toJson() << std::endl;
	//		std::cout << "popName: " << popName << std::endl;
			for(const auto & tran : fullAATyped[popName]){
				if("" != typed){
					typed += ";";
				}
				typed += tran.first + "=" + tran.second;
	//			std::cout << "\t" << tran.first << ":"<< tran.second << std::endl;
	//			std::cout << "\t" << "typed:" << typed << std::endl;
			}
			clus.meta_.addMeta("h_AATyped", typed);
		}
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		sampColl.printSampleCollapseInfo(
				njh::files::make_path(sampColl.masterOutputDir_,
						"selectedClustersInfo.tab.txt.gz"));
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;

//		table hapIdTab = sampColl.genHapIdTable();
//		hapIdTab.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(sampColl.masterOutputDir_,
//				"hapIdTable.tab.txt.gz"), true));
		sampColl.dumpPopulation();
//		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		SeqOutput::write(outPopSeqsPerSamp,SeqIOOptions(njh::files::make_path(sampColl.masterOutputDir_, "population", "popSeqsWithMetaWtihSampleName"), seqOpts.outFormat_));

		auto popMetaTable = seqsToMetaTable(outPopSeqsPerSamp);
		for(auto & row : popMetaTable){
			MetaDataInName::removeMetaDataInName(row[popMetaTable.getColPos("name")]);
		}

		popMetaTable.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(sampColl.masterOutputDir_, "population", "popSeqsWithMetaWtihSampleNameTable.tab.txt.gz")));
		sampColl.createCoreJsonFile();

	}

	};

	//njh::concurrent::runVoidFunctionThreaded(runPopClusOnTar, 1);
	njh::concurrent::runVoidFunctionThreaded(runPopClusOnTar, masterPopClusPars.numThreads);


	//gather all selected cluster info
	{

		OutputStream allSelectedInfo(njh::files::make_path(reportsDir, "allSelectedClustersInfo.tab.txt.gz"));
		std::shared_ptr<TableReader> firstTable;
		for(const auto & tar : tarNames){

			auto resultsFnp = njh::files::make_path(setUp.pars_.directoryName_, tar, "selectedClustersInfo.tab.txt.gz");
			if(bfs::exists(resultsFnp)){
				TableIOOpts currentTabOpts(InOptions(resultsFnp), "\t", true);
				if(nullptr == firstTable){
					firstTable = std::make_shared<TableReader>(currentTabOpts);
					firstTable->header_.outPutContents(allSelectedInfo, "\t");
					VecStr row;
					while(firstTable->getNextRow(row)){
						allSelectedInfo << njh::conToStr(row, "\t") << "\n";
					}
				}else{
					TableReader currentTable(currentTabOpts);
					if(currentTable.header_.columnNames_ != firstTable->header_.columnNames_){
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ", error " << "table " << resultsFnp << " header's doesn't match other's header"<< "\n";
						ss << "Expected Header: " << njh::conToStr(firstTable->header_.columnNames_, ",") << "\n";
						ss << "Current  Header: " << njh::conToStr(currentTable.header_.columnNames_, ",") << "\n";
						throw std::runtime_error{ss.str()};
					}
					VecStr row;
					while(currentTable.getNextRow(row)){
						allSelectedInfo << njh::conToStr(row, "\t") << "\n";
					}
				}
			}
		}
	}

	//zip all seqs file
	if(rawGatherRes.allSeqFnp != processedGatherRes.allSeqFnp){
		//if filtered, there will be two files
		bfs::path outAllSeqsFnp = rawGatherRes.allSeqFnp.string() + ".gz";
		njh::gzZipFile(IoOptions(InOptions(rawGatherRes.allSeqFnp), OutOptions(outAllSeqsFnp)));
		bfs::remove(rawGatherRes.allSeqFnp);
	}
	bfs::path outAllSeqsFnp = allSeqsFnp.string() + ".gz";
	njh::gzZipFile(IoOptions(InOptions(allSeqsFnp), OutOptions(outAllSeqsFnp)));
	bfs::remove(allSeqsFnp);


	return 0;
}





}  // namespace njhseq



