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
	setUp.setOption(rawGatherPars.minInputSeqLen, "--minInputSeqLen", "Min Input Seq Len");


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

	setUp.setOption(masterPopClusPars.collapseVarCallPars.variantCallerRunPars.occurrenceCutOff, "--occurrenceCutOff", "Occurrence Cut Off, don't report variants below this count");
	setUp.setOption(masterPopClusPars.collapseVarCallPars.variantCallerRunPars.lowVariantCutOff, "--lowVariantCutOff", "Low Variant Cut Off, don't report variants below this fraction");
  masterPopClusPars.collapseVarCallPars.calcPopMeasuresPars.lowVarFreq = masterPopClusPars.collapseVarCallPars.variantCallerRunPars.lowVariantCutOff;
  masterPopClusPars.collapseVarCallPars.transPars.setOptions(setUp);
  setUp.setOption(masterPopClusPars.collapseVarCallPars.calcPopMeasuresPars.getPairwiseComps, "--getPairwiseComps", "get Pairwise comparison metrics");
  bool noDiagAlnPairwiseComps = false;
  setUp.setOption(masterPopClusPars.collapseVarCallPars.noDiagAlnPairwiseComps, "--noDiagAlnPairwiseComps", "Use diagonal Alignment for Pairwise Comparisons");
  masterPopClusPars.collapseVarCallPars.calcPopMeasuresPars.diagAlnPairwiseComps = !noDiagAlnPairwiseComps;

  masterPopClusPars.collapseVarCallPars.ignoreSubFields = std::set<std::string>{"site:LabCross", "site:LabControl","site:LabContaminated"};
  setUp.setOption(masterPopClusPars.collapseVarCallPars.ignoreSubFields, "--ignoreSubFields", "Meta Sub Field values to ignore when calculating variants, e.g. --ignoreSubFields \"isFieldSample:TRUE,PreferredSample:FALSE\"");

  masterPopClusPars.collapseVarCallPars.calcPopMeasuresPars.numThreads = masterPopClusPars.numThreads;
 // masterPopClusPars.collapseVarCallPars.alnCacheDir = pars_.alnInfoDirName_;



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
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

	{
		OutputStream allBasicInfo(njh::files::make_path(reportsDir, "allBasicInfo.tab.txt.gz"));
		std::shared_ptr<TableReader> firstTable;
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;

		for(const auto & dir : directories){
			auto basicFnp = njh::files::make_path(dir, "final", "basicInfoPerRegion.tab.txt");
			if(bfs::exists(basicFnp)){
				TableIOOpts currentTabOpts(InOptions(basicFnp), "\t", true);
				firstTable = std::make_shared<TableReader>(currentTabOpts);
				firstTable->doNotCheckRowSizes = true;
				if(!skipRBind){
					firstTable->header_.outPutContents(allBasicInfo, "\t");
				}
				VecStr row;
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				while(firstTable->getNextRow(row)){
//					if(!skipRBind){
//						allBasicInfo << njh::conToStr(row, "\t") << "\n";
//					}
					allTargets.emplace(row[firstTable->header_.getColPos("name")]);
				}
				break;
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;

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
						currentTable.doNotCheckRowSizes = true;
						if(currentTable.header_.columnNames_ != firstTable->header_.columnNames_){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error " << "table " << basicFnp << " header's doesn't match other's header"<< "\n";
							ss << "Expected Header: " << njh::conToStr(firstTable->header_.columnNames_, ",") << "\n";
							ss << "Current  Header: " << njh::conToStr(currentTable.header_.columnNames_, ",") << "\n";
							throw std::runtime_error{ss.str()};
						}
						VecStr row;
//						std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			njh::concurrent::runVoidFunctionThreaded(gatherAllBasicInfoFiles, masterPopClusPars.numThreads);
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(skipRBind){
		bfs::remove(njh::files::make_path(reportsDir, "allBasicInfo.tab.txt.gz"));
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;


	std::shared_ptr<MultipleGroupMetaData> meta;
	if("" != masterPopClusPars.groupingsFile){
		meta = std::make_shared<MultipleGroupMetaData>(masterPopClusPars.groupingsFile);
	}



	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

	std::string sampleField = "sample";
	std::string targetField = "regionUID";

	auto sampleMetaDataFnp = njh::files::make_path(infoDir, "sampleMetaData.tab.txt");
	njh::stopWatch watch;
	watch.setLapName("Gathering Raw Seqs");
	auto rawAllSeqsFnp = njh::files::make_path(infoDir, "rawAllSeqsFile.fasta");
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

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
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	auto rawGatherRes = seqGatherer.gatherSeqsAndSortByTarget(rawGatherPars);
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(!rawGatherRes.missingOutput.empty()){
		OutputStream missingOut(njh::files::make_path(reportsDir, "missingDataForSamples.txt"));
		missingOut << njh::conToStr(rawGatherRes.missingOutput, "\n") << std::endl;
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

	processGatherPars.keepCommonSeqsWhenFiltering = keepCommonSeqsWhenFiltering;
	processGatherPars.allSeqFnp = rawGatherRes.allSeqFnp;
	processGatherPars.outMetaFnp = sampleMetaDataFnp;

	watch.startNewLap("Filtering Seqs On Meta");

//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	auto processedGatherRes = seqGatherer.processedGatherSeqsMeta(processGatherPars, rawGatherRes);
//	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

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
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

	auto allSeqsFnp = processedGatherRes.allSeqFnp;
	auto allSeqsFnpIn = SeqIOOptions::genFastaIn(allSeqsFnp);
	sleep(3); //have to sleep before building index
	SeqInput::buildIndex(allSeqsFnpIn);
	allSeqsFnpIn.processed_ = true;
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;

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

//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		// create aligner class object
		aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
				KmerMaps(setUp.pars_.colOpts_.kmerOpts_.kLength_),
				setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
				setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
		alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);
		njhseq::concurrent::AlignerPool alnPool(alignerObj,numThreadsForSample );
		alnPool.initAligners();
		alnPool.outAlnDir_ = setUp.pars_.outAlnInfoDirName_;
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		{
			njh::concurrent::LockableQueue<std::string> sampleQueue(processedGatherRes.allSamples);
//			std::cout << njh::conToStr(processedGatherRes.allSamples, ",") << std::endl;
//			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
//			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			njh::concurrent::runVoidFunctionThreaded(setupClusterSamples,numThreadsForSample );
//			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		}
		if(setUp.pars_.verbose_){
			std::cout << njh::bashCT::boldGreen("Pop Clustering") << std::endl;
		}


//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//first population clustering createSharedPathwaysFromReads
//		auto popInput = sampColl.createPopInput();
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//		std::cout << "popInput.size() : "<< popInput.size() << std::endl;
		{

			auto populationInput = sampColl.createPopInput();
			//std::cout << populationInput.size() << std::endl;
//			std::cout <<njh::bashCT::boldRed("Sleeping......") << std::endl;;
//			using namespace std::chrono_literals;
//			std::this_thread::sleep_for(100000s);
			sampColl.doPopulationClustering(populationInput, alignerObj, collapserObj, currentPars.popIteratorMap);

		}
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(setUp.pars_.verbose_){
			std::cout << njh::bashCT::boldRed("Done Pop Clustering") << std::endl;
		}
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		std::vector<seqInfo> popSeqsPerSamp;
		std::vector<seqInfo> outPopSeqsPerSamp;

		std::unordered_map<std::string, uint32_t> sampCountsForPopHaps;
		std::unordered_map<std::string, std::unordered_set<std::string>> sampNamesForPopHaps;

		uint32_t totalPopCount = 0;
		std::set<std::string> samplesCount;
//		std::cout << njh::conToStr(processedGatherRes.allSamples) << std::endl;
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		popSeqsPerSamp = sampColl.genOutPopSeqsPerSample();
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
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

	  auto varCallDirPath = njh::files::make_path(sampColl.masterOutputDir_,  "variantCalling");
	  if(!masterPopClusParsCopyConst.collapseVarCallPars.transPars.lzPars_.genomeFnp.empty()){
	  	auto collapseVarCallParsCopy = masterPopClusParsCopyConst.collapseVarCallPars;

	    //if genome set call variants against genome
	    collapseVarCallParsCopy.identifier = tar;
	    collapseVarCallParsCopy.outputDirectory = varCallDirPath;
	    collapseAndCallVariants(collapseVarCallParsCopy, outPopSeqsPerSamp);
	  }

	  std::map<std::string, std::string> fullAATyped;
	  //if typing file exists, read it in and set in map
	  auto seqTypingFnp = njh::files::make_path(varCallDirPath, "variantCalls/seqsAATyped.tab.txt.gz");
	  if(bfs::exists(seqTypingFnp)){
	    std::unordered_map<std::string, std::string> nameKey;
	    auto seqNamesFnp = njh::files::make_path(varCallDirPath, "uniqueSeqs_meta.tab.txt.gz");
	    {
	      TableReader seqNamesTab(TableIOOpts::genTabFileIn(seqNamesFnp, true));
	      VecStr row;
	      while(seqNamesTab.getNextRow(row)){
	        nameKey[row[seqNamesTab.header_.getColPos("CollapsedName")]] = row[seqNamesTab.header_.getColPos("PopUID")];
	      }
	    }
	    {
	      VecStr row;
	      TableReader seqTypingTab(TableIOOpts::genTabFileIn(seqTypingFnp, true));
	      while(seqTypingTab.getNextRow(row)){
	      	auto collapsedName = row[seqTypingTab.header_.getColPos("name")];
	        auto typed = row[seqTypingTab.header_.getColPos("fullTyped")];
	        auto popName = nameKey[collapsedName];
	        fullAATyped[popName] = typed;
	      }
	    }
	  }

//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for(auto & clus : sampColl.popCollapse_->collapsed_.clusters_){
			auto popName = clus.seqBase_.name_.substr(0, clus.seqBase_.name_.rfind("_f"));
			std::string typed = fullAATyped[popName];
			clus.meta_.addMeta("h_AATyped", typed);
		}
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		sampColl.printSampleCollapseInfo(
				njh::files::make_path(sampColl.masterOutputDir_,
						"selectedClustersInfo.tab.txt.gz"));
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;

//		table hapIdTab = sampColl.genHapIdTable();
//		hapIdTab.outPutContents(TableIOOpts::genTabFileOut(njh::files::make_path(sampColl.masterOutputDir_,
//				"hapIdTable.tab.txt.gz"), true));
		sampColl.dumpPopulation();
//		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
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



