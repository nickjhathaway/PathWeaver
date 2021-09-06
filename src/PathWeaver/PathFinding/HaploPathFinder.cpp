/*
 * HaploPathFinder.cpp
 *
 *  Created on: May 22, 2017
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

#include "HaploPathFinder.hpp"

namespace njhseq {

uint32_t countEndHomopolymerLength(const std::string & seq){
	uint32_t ret = 1;
	if(seq.size() > 1){
		uint32_t pos = seq.size() - 1;
		while(pos > 0 && seq[pos - 1] == seq.back()){
			++ret;
			--pos;
		}
	}
	return ret;
}



Json::Value HaploPathFinder::toJson() const{
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["pars_"] = njh::json::toJson(pars_);
	ret["gMapper_"] = njh::json::toJson(gMapper_);


	return ret;
}


HaploPathFinder::ExtractParams::ExtractParams(){
	unmappedPars_.minQuerySize_ = 50;
}

Json::Value HaploPathFinder::ExtractParams::toJson() const{
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["inputDir"] = njh::json::toJson(inputDir);
	ret["extractionDir"] = njh::json::toJson(extractionDir);
	ret["regionDir"] = njh::json::toJson(regionDir);
	ret["overWriteDir"] = njh::json::toJson(overWriteDir);
	ret["metaDataFnp"] = njh::json::toJson(metaDataFnp);
	ret["genomeDir"] = njh::json::toJson(genomeDir);
	ret["primaryGenome"] = njh::json::toJson(primaryGenome);
	ret["selectedGenomes"] = njh::json::toJson(selectedGenomes);
	ret["genomeFnp"] = njh::json::toJson(genomeFnp);
	ret["outputDir"] = njh::json::toJson(outputDir);
	ret["trimSeqBedFnp"] = njh::json::toJson(trimSeqBedFnp);
	ret["bamExtractPars_"] = njh::json::toJson(bamExtractPars_);
	ret["unmappedPars_"] = njh::json::toJson(unmappedPars_);

	ret["region_"] = njh::json::toJson(region_);
	ret["pFinderPars_"] = njh::json::toJson(pFinderPars_);

	ret["maxPerBaseCoverage"] = njh::json::toJson(maxPerBaseCoverage);


	return ret;
}


void HaploPathFinder::ExtractParams::setGenomeOpts(seqSetUp & setUp, bool genomeRequired){
	setUp.setOption(genomeDir, "--genomeDir", "Name of the genome file fnp", genomeRequired);
	setUp.setOption(primaryGenome, "--primaryGenome", "Name of the primary genome, should be in --genomeDir", genomeRequired);
	if(genomeRequired){
		setGenomeFnp();
	}
	setUp.setOption(selectedGenomes, "--selectedGenomes", "Name of the other genomes in --genomeDir to be read in, leave blank to just do all fastas");


}



HaploPathFinder::PathFinderCorePars::PathFinderCorePars(){
	errorsToAllow_.hqMismatches_ = 1;
	trimEdgesOfEndNodesPars_.entropyCutOff = 1.50;
	pairProcessorParams_.minOverlap_ = 10;


	seqEdgeEntropyTrimPars_.entropyCutOff = 0.75;
	seqEdgeEntropyTrimPars_.kLen = 1;
	seqEdgeEntropyTrimPars_.windowSize = 20;
	seqEdgeEntropyTrimPars_.windowStep = 1;


	circularTrimPars_.mark_ = true;
	circularTrimPars_.preferHeader_ = true;

}

void HaploPathFinder::PathFinderCorePars::setTandemRepeatHandlingOptions(seqSetUp & setUp){
	bool dontCollapseFinalDeterminedTandems = false;
	setUp.setOption(dontCollapseFinalDeterminedTandems, "--dontCollapseFinalDeterminedTandems", "Don't Collapse Final Determined Tandems");
  collapseFinalDeterminedTandems_ = !dontCollapseFinalDeterminedTandems;

  finalTandemOccurenceCutOff_ = kmerKOcurrenceCutOffs.front(); /**@todo add up tandem counts when collapsing*/

  setUp.setOption(finalTandemOccurenceCutOff_, "--finalTandemOccurenceCutOff", "Don't Collapse Final Determined Tandems");
	setUp.setOption(useSmallerThanKlenExpandedPositions_, "--useSmallerThanKlenExpandedPositions", "Use Smaller Than Klen Expanded Positions");

}

void HaploPathFinder::PathFinderCorePars::setOptimizationParameters(seqSetUp & setUp){


	//for estimating coverage

	setUp.setOption(estimatorKlen, "--estimatorKlen", "K-mer length to estimate coverage");

	//kmer length and optimization
	optKLenStart = 31;
	optKLenStop = 41;
	optKLenStep = 10;
	kmerLengths = {31, 41};
	bool optKLenStartSet = setUp.setOption(optKLenStart,"--optKLenStart","Optimize k-mer Length Start");
	bool optKLenStopSet =  setUp.setOption(optKLenStop, "--optKLenStop", "Optimize k-mer Length Stop");
	bool optKLenStepSet =  setUp.setOption(optKLenStep, "--optKLenStep", "Optimize k-mer Length Step");
	if(optKLenStartSet || optKLenStopSet || optKLenStepSet){
		if(optKLenStop < optKLenStart){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr("optKLenStop: ",optKLenStop, " shouldn't be less than optKLenStart: ", optKLenStart));
		}else{
			kmerLengths.clear();
			for(const auto klen : iter::range(optKLenStart, optKLenStop + 1, optKLenStep)){
				kmerLengths.emplace_back(klen);
			}
		}
	}
	setUp.setOption(kmerLengths, "--kmerLengths", "K-mer Lengths to use in reconstructions");
	njh::sort(kmerLengths);

	//kmer occurrence cut off and optimization
	optKcutStart = 2;
	//optKcutStop = 5;
	optKcutStop = 5;
	optKcutStep = 1;
	//kmerKOcurrenceCutOffs = {2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40};
	//kmerKOcurrenceCutOffs = {2,3,4,5};
	kmerKOcurrenceCutOffs = {2,3,4, 5};
	bool optKcutStartSet = setUp.setOption(optKcutStart,"--optKcutStart", "Optimize k-mer occurrence cut off start");
	bool optKcutStopSet =  setUp.setOption(optKcutStop, "--optKcutStop",  "Optimize k-mer occurrence cut off stop");
	bool optKcutStepSet =  setUp.setOption(optKcutStep, "--optKcutStep",  "Optimize k-mer occurrence cut off step");
	if(optKcutStartSet || optKcutStopSet || optKcutStepSet){
		if(optKcutStop < optKcutStart){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr("optKcutStop: ",optKcutStop, " shouldn't be less than optKcutStart: ", optKcutStart));
		}else{
			kmerKOcurrenceCutOffs.clear();
			for(const auto kcut : iter::range(optKcutStart, optKcutStop + 1, optKcutStep)){
				kmerKOcurrenceCutOffs.emplace_back(kcut);
			}
		}
	}
	setUp.setOption(kmerKOcurrenceCutOffs, "--kmerKOcurrenceCutOffs", "K-mer occurrence cut offs to use in reconstructions");
	njh::sort(kmerKOcurrenceCutOffs);

	setUp.setOption(initialKmerKOcurrenceCutOffs, "--initialKmerKOcurrenceCutOffs", "K-mer occurrence cut offs to use in the first reconstructions when doing multiple iterations");
	njh::sort(initialKmerKOcurrenceCutOffs);

	setUp.setOption(initialShortTipNumbers, "--initialShortTipNumbers", "initial Short Tip Numbers");
	njh::sort(initialShortTipNumbers);

	setUp.setOption(initialKmerLengths, "--initialKmerLengths", "Initial Kmer Lengths");
	njh::sort(initialKmerLengths);


	setUp.setOption(forceDetermineKCuts, "--forceDetermineKCuts", "Force determine occurrence cut offs based of percentage of per base coverage for a region, percentages set by --autoDetermineKCutsPercentages");
	if (autoDetermineKCuts) {
		bool noAutoDKCuts = false;
		setUp.setOption(noAutoDKCuts, "--noAutoDetermineKCuts",
				"Auto determine occurrence cut offs based of percentage of per base coverage for a region for regions with at least 100x coverage, percentages set by --autoDetermineKCutsPercentages");
		autoDetermineKCuts = !noAutoDKCuts;
	} else {
		if(setUp.setOption(autoDetermineKCuts, "--autoDetermineKCuts",
				"Auto determine occurrence cut offs based of percentage of per base coverage for a region, percentages set by --autoDetermineKCutsPercentages")){
			if(autoDetermineKCutsOnNodeCount){
				//if setting auto determine counts based on per base coverage and the auto determine node count is set by default to be true, now set it to be false.
				autoDetermineKCutsOnNodeCount = false;
			}
			if(autoDetermineKCutsOnTotalBaseCount){
				autoDetermineKCutsOnTotalBaseCount = false;
			}
		}
	}



	if (autoDetermineKCutsOnNodeCount) {
		bool noAutoDKCuts = false;
		setUp.setOption(noAutoDKCuts, "--noAutoDetermineKCutsOnNodeCount",
				"Auto determine occurrence cut offs based of percentage k-mer node counts with at least 100x coverage, percentages set by --autoDetermineKCutsPercentages");
		autoDetermineKCutsOnNodeCount = !noAutoDKCuts;
	} else {
		setUp.setOption(autoDetermineKCutsOnNodeCount, "--autoDetermineKCutsOnNodeCount",
				"Auto determine occurrence cut offs based of percentage k-mer node counts with at least 100x coverage, percentages set by --autoDetermineKCutsPercentages");
	}
	if (autoDetermineKCutsOnTotalBaseCount) {
		bool noAutoDKCuts = false;
		setUp.setOption(noAutoDKCuts, "--noAutoDetermineKCutsOnTotalBaseCount",
				"Auto determine occurrence cut offs based of percentage k-mer total bases divided by the input region with at least 100x coverage, percentages set by --autoDetermineKCutsPercentages");
		autoDetermineKCutsOnTotalBaseCount = !noAutoDKCuts;
	} else {
		setUp.setOption(autoDetermineKCutsOnTotalBaseCount, "--autoDetermineKCutsOnTotalBaseCount",
				"Auto determine occurrence cut offs based of percentage k-mer node counts with at least 100x coverage, percentages set by --autoDetermineKCutsPercentages");
	}



	setUp.setOption(autoDetermineKCutsPercentages, "--autoDetermineKCutsPercentages", "Percentages to use when auto determining occurrence cut offs");
	setUp.setOption(initialAutoDetermineKCutsPercentages, "--initialAutoDetermineKCutsPercentages", "Percentages to use when auto determining occurrence cut offs for the first step");


	//other optimization parameters
	setUp.setOption(keepOptimizedSubDirs, "--keepOptimizedSubDirs", "keepOptimizedSubDirs");
	bool useMinLenOptimalCount = false;
	setUp.setOption(useMinLenOptimalCount, "--useMinLenOptimalCount", "Use Min Length Optimal Count");
	useFullOptimalCount_ = !useMinLenOptimalCount;
	setUp.setOption(optimizeNodeSizeCutOff_, "--optimizeNodeSizeCutOff", "Optimize Node Size Cut Off");

	//optimization attempts
	setUp.setOption(forceFullKcutOptimize_, "--forceFullKcutOptimize", "Force Full Kcut Optimize");
	setUp.setOption(kCutOptimizationAttempts, "--kCutOptimizationAttempts", "Number of k occurrence cut off optimization attempts");
	setUp.setOption(forceFullShortTipOptimize_, "--forceFullShortTipOptimize", "Force Full short tip Optimize");
	setUp.setOption(shortTipOptimizationAttempts, "--shortTipOptimizationAttempts", "Number of short tip optimization attempts");


	//optimization input used
	setUp.setOption(percentOfInputReadsUsedForOptimization, "--percentOfInputReadsUsedForOptimization", "Percent Of Input Reads Used For Optimization");
	setUp.setOption(percentOfInputReadsUsedForOptimizationStepDown, "--percentOfInputReadsUsedForOptimizationStepDown", "Percent Of Input Reads Used For Optimization Step Down");
	bool doNotCalcPercentUsedByKmerUsage = false;
	setUp.setOption(doNotCalcPercentUsedByKmerUsage, "--doNotCalcPercentUsedByKmerUsage", "Don't Calc Percent Used By Kmer Usage");
	calcPercentUsedByKmerUsage = !doNotCalcPercentUsedByKmerUsage;
//	setUp.setOption(calcPercentUsedByKmerUsage, "--calcPercentUsedByKmerUsage", "Calculate Percent of reads Used By Kmer Usage");
	setUp.setOption(percentageOfKmersUsedCutOff, "--percentageOfKmersUsedCutOff", "Percentage Of Kmers Used Cut Off for --calcPercentUsedByKmerUsage");


	//trimming of short tips
	bool optShortTipNumberStartSet = setUp.setOption(optShortTipNumberStart,"--optShortTipNumberStart", "Optimize Short Tip Pathway length to trim off start");
	bool optShortTipNumberStopSet  = setUp.setOption(optShortTipNumberStop, "--optShortTipNumberStop",  "Optimize Short Tip Pathway length to trim off stop");
	bool optShortTipNumberStepSet  = setUp.setOption(optShortTipNumberStep, "--optShortTipNumberStep",  "Optimize Short Tip Pathway length to trim off step");
	if(optShortTipNumberStartSet || optShortTipNumberStopSet || optShortTipNumberStepSet){
		if(optShortTipNumberStop < optShortTipNumberStart){
			setUp.failed_ = true;
			setUp.addWarning(njh::pasteAsStr("optShortTipNumberStop: ",optShortTipNumberStop, " shouldn't be less than optShortTipNumberStart: ", optShortTipNumberStart));
		}else{
			shortTipNumbers.clear();
			for(const auto shortTip : iter::range(optShortTipNumberStart, optShortTipNumberStop + 1, optShortTipNumberStep)){
				shortTipNumbers.emplace_back(shortTip);
			}
		}
	}
	setUp.setOption(shortTipNumbers, "--shortTipNumbers", "Short Tip Pathway length to trim off");
	njh::sort(shortTipNumbers);

	bool noTrimShortTips = false;
	setUp.setOption(noTrimShortTips, "--noTrimShortTips", "Don't Trim Short Tips that are only a few kmers");
	trimShortTips_ = !noTrimShortTips;
	if (!trimAllShortTips_) {
		setUp.setOption(trimAllShortTips_, "--trimAllShortTips", "Trim Short Tips all short tips");
	} else {
		bool trimShortTipsBasedOnFreq = false;
		setUp.setOption(trimShortTipsBasedOnFreq, "--trimShortTipsBasedOnFreq", "Trim Short Tips only if they occur at ~30% of the median kmer frequency");
		trimAllShortTips_ = !trimShortTipsBasedOnFreq;
	}
}

void HaploPathFinder::PathFinderCorePars::setStitchParameters(seqSetUp & setUp){
	bool noStitchPairs = false;
	setUp.setOption(noStitchPairs, "--noStitchPairs", "Don't attempt to stitch pairs");
	stitchPairs = !noStitchPairs;

	if(reOrientPairsForStitching){
		bool noReOrientPairsForStitching = false;
		setUp.setOption(reOrientPairsForStitching, "--noReOrientPairsForStitching", "Don't reorient pairs when stitching, this should only be turned on if adaptors have been trimmed off");
		reOrientPairsForStitching = !noReOrientPairsForStitching;
	}else{
		setUp.setOption(reOrientPairsForStitching, "--reOrientPairsForStitching", "Reorient pairs when stitching, this should only be turned on if adaptors have been trimmed off");
	}


	setUp.setOption(pairProcessorParams_.minOverlap_, "--stitchPairsMinOverlap", "stitch pairs min overlap");
	setUp.setOption(pairProcessorParams_.errorAllowed_, "--stitchErrorAllowed", "Stitch Error Allowed");
	setUp.setOption(pairProcessorParams_.hardMismatchCutOff_, "--stitchHardMismatchCutOff", "Stitch Hard Mismatch Cut Off");
	setUp.setOption(pairProcessorParams_.lqMismatchCutOff, "--stitchLqMismatchCutOff", "Stitch Lq Mismatch Cut Off");
	setUp.setOption(pairProcessorParams_.hqMismatchCutOff, "--stitchHqMismatchCutOff", "Stitch Hq Mismatch Cut Off");
}

void HaploPathFinder::PathFinderCorePars::setOutlierRemovalOptions(seqSetUp & setUp){
	setUp.setOption(removePossibleOutliersKSim, "--removePossibleOutliersKSim", "Remove Possible Outliers based off of a kmer similarity score");
	setUp.setOption(outlierKLen, "--outlierKLen", "Outlier K Len to use");
	setUp.setOption(outliersKSimCutOff, "--outliersKSimCutOff", "Outliers Kmer Similarity Cut Off, a similarity score below this will remove that contig");
	setUp.setOption(filterOffOutlierInputSeqs, "--filterOffOutlierInputSeqs", "Filter Off Outlier Input Seqs");

	//setUp.setOption(keepGroupsForOutlierFilter, "--keepGroupsForOutlierFilter", "Keep Groups For Outlier Filter");
	bool disregardGroupsForOutliersRemoval = false;
	setUp.setOption(disregardGroupsForOutliersRemoval, "--disredgardGroupsForOutliersRemoval", "Disregard Groups For Outliers Removal");
	keepGroupsForOutlierFilter =!disregardGroupsForOutliersRemoval;
}

void HaploPathFinder::PathFinderCorePars::setQualityTrimAndFiltOpts(seqSetUp & setUp){
	setUp.setOption(qualTrim_, "--qualTrim", "Qual to trim at if trimming on quals, will trim at first occurence");
	setUp.setOption(trimOnQual_, "--trimOnQual", "Whether to trim input on quality, which quality score is trimmed at is controlled by qualTrim flag");
	setUp.setOption(qualCheckCutOff_, "--qualCheckCutOff", "Fraction of a read has to be above this of the given quality by --qualCheck");
	setUp.setOption(qualCheck_, "--qualCheck", "Per base quality check to for filtering based off of quality");
	setUp.setOption(trimSeqsWithNs_, "--trimSeqsWithNs", "Trim Seqs With Ns to the largest sub seq without Ns");

	if(preFilterReadsOnEntropy_){
		bool noPreFilterReadsOnEntropy = false;
		setUp.setOption(noPreFilterReadsOnEntropy, "--noPreFilterReadsOnEntropy", "no pre-filtering reads on entropy");
		preFilterReadsOnEntropy_ = !noPreFilterReadsOnEntropy;
	}else{
		setUp.setOption(preFilterReadsOnEntropy_, "--preFilterReadsOnEntropy", "pre-filter reads on entropy");
	}
	setUp.setOption(preFilterReadsOnEntropyCutOff_, "--preFilterReadsOnEntropyCutOff", "Pre-filtering reads on entropy cut off");
	if(markPreFilterInfo_){
		bool noMarkPreFilterInfo = false;
		setUp.setOption(noMarkPreFilterInfo, "--noMarkPreFilterInfo", "no mark pre-filtered reads with failed check info");
		markPreFilterInfo_ = !noMarkPreFilterInfo;
	}else{
		setUp.setOption(markPreFilterInfo_, "--markPreFilterInfo", "mark pre-filtered reads with failed check info");
	}

	if(!trimSeqsEdgesForLowEntropy_){
		setUp.setOption(trimSeqsEdgesForLowEntropy_, "--trimSeqsEdgesForLowEntropy", "Trim Seqs Edges For Low Entropy");
	}else{
		bool noTrimSeqsEdgesForLowEntropy = false;
		setUp.setOption(noTrimSeqsEdgesForLowEntropy, "--noTrimSeqsEdgesForLowEntropy", "No Trim Seqs Edges For Low Entropy");
		trimSeqsEdgesForLowEntropy_ = !noTrimSeqsEdgesForLowEntropy;
	}
	setUp.setOption(seqEdgeEntropyTrimPars_.entropyCutOff, "--trimSeqsEdgesForLowEntropy-entropyCutOff", "Entropy Cut Off for when --trimSeqsEdgesForLowEntropy is true");
	setUp.setOption(seqEdgeEntropyTrimPars_.windowSize, "--trimSeqsEdgesForLowEntropy-windowSize", "Window Size for when --trimSeqsEdgesForLowEntropy is true");
	setUp.setOption(seqEdgeEntropyTrimPars_.windowStep, "--trimSeqsEdgesForLowEntropy-windowStep", "Window Step for when --trimSeqsEdgesForLowEntropy is true");


	setUp.setOption(kmerCommonLocKmerLength, "--kmerCommonLocKmerLength", "Kmer Common Loc Kmer Length");
	setUp.setOption(kmerCommonLocOccurenceCutOff, "--kmerCommonLocOccurenceCutOff", "kmer Common Loc Occurence Cut Off");
	setUp.setOption(kmerCommonLocStdCutOff, "--kmerCommonLocStdCutOff", "kmer Common Loc Std Cut Off");
	setUp.setOption(kmerCommonLocWithin, "--kmerCommonLocWithin", "kmer Common Loc Within");

	if(!filterbyKmerCommonLoc_){
		setUp.setOption(filterbyKmerCommonLoc_, "--filterByKmerCommonLoc", "filter by Kmer Common Loc");
	}else{
		bool nofilterbyKmerCommonLoc = false;
		setUp.setOption(nofilterbyKmerCommonLoc, "--noFilterByKmerCommonLoc", "Don't filter by Kmer Common Loc");
		filterbyKmerCommonLoc_ = !nofilterbyKmerCommonLoc;
	}

//
}


void HaploPathFinder::PathFinderCorePars::setFurtherSplittingOpts(seqSetUp & setUp){
	//splitting

	if(splitEndsOnce_){
		bool doNotSplitEndsOnce = false;
		setUp.setOption(doNotSplitEndsOnce, "--doNotSplitEndsOnce", "Do Not Split Ends Nodes Once");
		splitEndsOnce_ = !doNotSplitEndsOnce;
	}else{
		setUp.setOption(splitEndsOnce_, "--splitEndsOnce", "Split Ends Nodes Once");
	}



	setUp.setOption(splitTailed_, "--splitTailed", "Split Tailed Nodes");
	setUp.setOption(splitEnds_, "--splitEnds", "Split End Nodes");
	setUp.setOption(maxSplitEndNodeSize_, "--maxSplitEndNodeSize", "Max Split End Node Size");
	setUp.setOption(keepCycles_, "--keepCycles", "Keep Cycles Pathways seqs, this have a higher change of being false");
	setUp.setOption(splitToRecruit_, "--splitToRecruit", "Split To Recruit");
	setUp.setOption(addHeadTailSeqsToRecruit, "--addHeadTailSeqsToRecruit", "Add Head Tail Seqs To Recruit, this allows longer pieces to map to");
	setUp.setOption(initialBreakingSinglePaths_, "--initialBreakingSinglePaths_", "Initial Breaking Single Paths");
}

void HaploPathFinder::PathFinderCorePars::setHeadlessTaillessNodesOpts(seqSetUp & setUp){
	//headlessTailessLenCutOff = kmerLengths.empty() ? 0 : kmerLengths.front() + 1;
	setUp.setOption(headlessTailessLenCutOff, "--headlessTailessLenCutOff", "headless Tailess Len Cut Off (inclusive)");
	bool keepHeadlessTaillessAlongTheWay = false;
	setUp.setOption(keepHeadlessTaillessAlongTheWay, "--keepHeadlessTaillessAlongTheWay", "Keep small Headless Tailless Along The Way, these represent reads that might have been recruited accidentally");
	removeHeadlessTaillessAlongTheWay = !keepHeadlessTaillessAlongTheWay;

	bool keepHeadlessTaillessAfterDisentaglement = false;
	setUp.setOption(keepHeadlessTaillessAfterDisentaglement, "--keepHeadlessTaillessAfterDisentaglement", "If removing headless and tailless along the way, do so again after disentanglement");
	removeHeadlessTaillessAfterDisentaglement = !keepHeadlessTaillessAfterDisentaglement;

}

void HaploPathFinder::PathFinderCorePars::setNodeProcessingOpts(seqSetUp & setUp){
	bool noCollapseOneBaseIndelsNodes = false;
	setUp.setOption(noCollapseOneBaseIndelsNodes, "--noCollapseOneBaseIndelsNodes", "Don't Combine nodes that share the same heads and tails that only differ by a 1 base indel in a long homopolymer");
	collapseOneBaseIndelsNodes_ = ! noCollapseOneBaseIndelsNodes;
	setUp.setOption(collapseOneBaseIndelsBeforeDisentanglement_, "--collapseOneBaseIndelsBeforeDisentanglement", "Collapse One Base Indels Before Disentanglement");
	setUp.setOption(collapsePossibleSimpleLoops_, "--collapsePossibleSimpleLoops", "Collapse Possible Simple Loops");


	setUp.setOption(oneBaseIndelError_, "--oneBaseIndelError", "The amount of one base indel error to allow when collapsing nodes");
	setUp.setOption(twoBaseIndelError_, "--twoBaseIndelError", "The amount of two base indel error to allow when collapsing nodes");
	setUp.setOption(largeBaseIndelError_, "--largeBaseIndelError", "The amount of >=3 base indel error to allow when collapsing nodes");
	setUp.setOption(homopolymerIndelCollapseFreqMultiplier_, "--homopolymerIndelCollapseFreqMultiplier", "Frequency Multiplier to use for collapsing");
	setUp.setOption(nodeBridgingReadPercCutOff_, "--nodeBridgingReadPercCutOff", "The percent of reads entering and leaving a tailed and multi head node to consider connecting the heads and tails together");

	bool doNotStartDisentanglementConservative = false;
	setUp.setOption(doNotStartDisentanglementConservative, "--doNotStartDisentanglementConservative", "Do not start disentanglement conservatively");
	startDisentanglementConservative_ = !doNotStartDisentanglementConservative;

	setUp.setOption(adjustForHeadLessAddedAfterDisentanglement,    "--adjustForHeadLessAddedAfterDisentanglement", "Adjust For Head Less Added After Disentanglement");
	setUp.setOption(throwAwayConservedAddNodesDuringDisentaglement,"--throwAwayConservedAddNodesDuringDisentaglement", "Throw Away Conserved Added Nodes During Disentaglement, this would decrease the amount of sequenced output that is shared between mix strains to be given to the minor strains");

	setUp.setOption(trimTipsOfLowEntropyNodes_,       "--trimTipsOfLowEntropyNodes",       "Trim Tips Of Low Entropy Nodes");
	setUp.setOption(trimTipsOfLowEntropyNodesCutOff_, "--trimTipsOfLowEntropyNodesCutOff", "Trim Tips Of Low Entropy Nodes Cut Off to use");

	setUp.setOption(trimEdgesNodeTipsWithLowEntropy_,       "--trimEdgesNodeTipsWithLowEntropy",       "Trim the edges of the Tips Of Low Entropy Nodes");
	setUp.setOption(trimEdgesOfEndNodesPars_.entropyCutOff, "--trimEdgesNodeTipsWithLowEntropyCutOff", "trim Edges Node Tips With Low Entropy Cut Off");
	setUp.setOption(trimEdgesOfEndNodesPars_.windowSize,    "--trimEdgesNodeTipsWithLowEntropyWindowSize", "trim Edges Node Tips With Low Entropy Window Size");
	setUp.setOption(trimEdgesOfEndNodesPars_.windowStep,    "--trimEdgesNodeTipsWithLowEntropyWindowStep", "trim Edges Node Tips With Low Entropy Window Step");
	setUp.setOption(trimEdgesOfEndNodesPars_.kLen,          "--trimEdgesNodeTipsWithLowEntropykLen", "trim Edges Node Tips With Low Entropy kmer Len");


	setUp.setOption(collapseNodesWithAllowableError_, "--collapseNodesWithAllowableError", "Collapse bubbles (nodes that have the same heads and tails) that differ by a certain amount of error and frequency");
	setUp.setOption(collapseNodesWithAllowableErrorFreqMultiplier_, "--collapseNodesWithAllowableErrorFreqMultiplier", "Collapse Nodes With Allowable Error Frequency Multiplier");
	setUp.processComparison(errorsToAllow_, "--nodeAllowableError");



	setHeadlessTaillessNodesOpts(setUp);
}

void HaploPathFinder::PathFinderCorePars::setRunningOpts(seqSetUp & setUp){
	setUp.setOption(graphDebug_, "--graphDebug", "Set graph to debug mode");
	setUp.setOption(debugWriteEdgeInfo_, "--debugWriteEdgeInfo", "Write edge info when graph is set to debug mode");

	setUp.setOption(graphVerbose_, "--graphVerbose", "Set graph to verbose mode");
	setUp.setOption(writeOutDotsByGroups_, "--writeOutDotsByGroups", "Write Out Final Dot file for the final contigs with a dot file for each group");
	setUp.setOption(writeOutDotsByGroupsSizeCutOff_, "--writeOutDotsByGroupsSizeCutOff", "Write Out Final Dot file for the final contigs with a dot file for each group size of group cut off(groups with this count and above will be written)");

	setUp.setOption(writeOutFinalDot_, "--writeOutFinalDot", "Write Out Final Dot file for the final contigs");
	setUp.setOption(exitOnException_, "--exitOnException", "Exit On Exception");
	setUp.setOption(numThreads, "--numThreads", "Number of CPUs to utilize");
	setUp.setOption(writeNodeCounts_, "--writeNodeCounts", "write Node Counts");

	setUp.setOption(rawInputSeqs, "--rawInputSeqs", "By Default the input seqs are automatically made upper case, use this flag to keep the input as it was");

}

void HaploPathFinder::PathFinderCorePars::setPossibleHapsOpts(seqSetUp & setUp){
	setUp.setOption(writeAllPossibleHaps_, "--writeAllPossibleHaps", "Write All Possible Haps");
	setUp.setOption(writeEstimatedMajorHaps_, "--writeEstimatedMajorHaps", "Write Estimated Major Haps");
	setUp.setOption(writeOutFinalConnections_, "--writeOutFinalConnections", "Write out possible connections for multitailed and multiheaded nodes utilizing firstgraph");
}


void HaploPathFinder::PathFinderCorePars::setAddingGroupInfoOpts(seqSetUp & setUp){
	bool noGroupInfo = false;
	setUp.setOption(noGroupInfo, "--noAddGroupInfo", "Don't add group information");
	addGroups_ = !noGroupInfo;
	setUp.setOption(writeGroupInputNames_, "--writeInputNamesPerGroup", "Write Input Read Names Per Group");
}


void HaploPathFinder::PathFinderCorePars::setPostProcessTrimmingOpts(seqSetUp & setUp){
	if (trimToInputSeqs) {
		bool noTrimToInputSeqs = false;
		setUp.setOption(noTrimToInputSeqs, "--noTrimEachIteration", "Do not trim to the input references after each iteration");
		trimToInputSeqs = !noTrimToInputSeqs;

	} else {
		setUp.setOption(trimToInputSeqs, "--trimEachIteration", "Trim to the input references after each iteration");
	}

	bool noTrimWithGlobalAln = false;
	setUp.setOption(noTrimWithGlobalAln, "--noTrimWithGlobalAln", "Trim to the input references with global alignment to input seqs rather than with muscle");
	trimWithGlobalAln = !noTrimWithGlobalAln;
	setUp.setOption(trimToCircularGenome, "--trimToCircularGenome", "Trim contigs to the reference assuming ref is circular");
	setUp.setOption(circularTrimPars_.extend_, "--trimToCircularGenomeExtendAmount", "When trimming to circular reference, extend by this much which can help with trimming");
	setUp.setOption(trimEnds, "--trimEnds", "Trim the ends of called haplotypes");
	setUp.setOption(trimEndsByCoverage, "--trimEndsByCoverage", "Trim Ends By low coverage");
	setUp.setOption(trimEndsBy, "--trimEndsBy", "Trim Ends By this length when trimming, rather than the kmer length");
}

void HaploPathFinder::PathFinderCorePars::setPostProcessContigHandlingOpts(seqSetUp & setUp){
	setUp.setOption(collapseFinalClusters_, "--collapseFinalClusters", "Collapse final clusters for homopolymer errors");
	setUp.setOption(doNotCollapseIdenticalContigs_, "--doNotCollapseIdenticalContigs", "Do not collapse final contigs for exact matches which might happen after post processing");
}




void HaploPathFinder::PathFinderCorePars::setDuplicatedSeqHandlingOpts(seqSetUp & setUp){
	if(removeDuplicatedSequences_){
		bool keepDuplicatedSeqs = false;
		setUp.setOption(keepDuplicatedSeqs, "--keepDuplicatedSeqs", "Keep Duplicated Seqs");
		removeDuplicatedSequences_ = !keepDuplicatedSeqs;
	}else{
		setUp.setOption(removeDuplicatedSequences_, "--removeDuplicatedSequences", "Remove identical sequences");
	}
//	setUp.setOption(trimRemoveDuplicatedSequences_, "--trimRemoveDuplicatedSequences", "Use Only unique input sequence use 1 less than");
}



void HaploPathFinder::ExtractParams::setBamExtractOpts(seqSetUp & setUp){

	setUp.setOption(maxPerBaseCoverage, "--maxPerBaseCoverage", "Sub sample reads to this read depth, helps to cut down on memory usage and compute time", false, "BamExtracting");

	if(setUp.setOption(bamExtractPars_.percentSubSample_, "--percentSubSample", "A number between 0 and 1 for sub sampling", false, "BamExtracting", njh::progutils::ProgramSetUp::flagCheckFrom0To1<double>("percentSubSample"))){
		maxPerBaseCoverage = std::numeric_limits<double>::max();
	}

	if (bamExtractPars_.keepMarkedDuplicate_) {
		bool removeMakredDups = false;
		setUp.setOption(removeMakredDups, "--removeMakredDuplicate",
				"Remove alignments marked as duplicate", false, "BamExtracting");
		bamExtractPars_.keepMarkedDuplicate_ = !removeMakredDups;
	} else {
		setUp.setOption(bamExtractPars_.keepMarkedDuplicate_,
				"--keepMarkedDuplicate", "Keep alignments marked as duplicate", false, "BamExtracting");
	}
	if (bamExtractPars_.filterOffLowEntropyOrphansRecruits_) {
		bool noFilterOffLowEntropyOrphansRecruits = false;
		setUp.setOption(noFilterOffLowEntropyOrphansRecruits,
				"--noFilterOffLowEntropyOrphansRecruits",
				"No Filter Off Low Entropy Orphans Recruits", false, "BamExtracting");
		bamExtractPars_.filterOffLowEntropyOrphansRecruits_ =
				!noFilterOffLowEntropyOrphansRecruits;
	} else {
		setUp.setOption(bamExtractPars_.filterOffLowEntropyOrphansRecruits_,
				"--filterOffLowEntropyOrphansRecruits",
				"Filter Off Low Entropy Orphans Recruits", false, "BamExtracting");
	}

	if (bamExtractPars_.removeImproperPairs_) {
		bool keepImproperMates = false;
		setUp.setOption(keepImproperMates, "--keepImproperMates",
				"Keep Pairs marked as improper (mate doesn't map, inverse mapping etc)", false, "BamExtracting");
		bamExtractPars_.removeImproperPairs_ = !keepImproperMates;
	} else {
		setUp.setOption(bamExtractPars_.removeImproperPairs_,
				"--removeImproperPairs",
				"Filter Off Pairs marked as improper (mate doesn't map, inverse mapping etc)", false, "BamExtracting");
	}

	if (bamExtractPars_.removeImproperPairs_) {
		bool removeImproperMateUnmapped = false;
		setUp.setOption(removeImproperMateUnmapped, "--removeImproperMateUnmapped",
				"When Filtering Off Pairs marked as improper, also remove the improper pairs where it's due to one mate not mapping)", false, "BamExtracting");
		bamExtractPars_.keepImproperMateUnmapped_ = !removeImproperMateUnmapped;
	} else {
		setUp.setOption(bamExtractPars_.keepImproperMateUnmapped_,
				"--keepImproperMateUnmapped",
				"When Filtering Off Pairs marked as improper,  keep the improper pairs where it's due to one mate not mapping", false, "BamExtracting");
	}



	setUp.setOption(bamExtractPars_.filterOffLowEntropyOrphansRecruitsCutOff_,  "--filterOffLowEntropyOrphansRecruitsCutOff",  "Filter Off Low Entropy Orphans Recruits Cut Off", false, "BamExtracting");
	setUp.setOption(bamExtractPars_.entropyKlen_,  "--entropyKlen",  "Entropy Klen", false, "BamExtracting");


	setUp.setOption(bamExtractPars_.tryToFindOrphansMate_,  "--tryToFindOrphanMates",  "Try To Find Orphan Mates", false, "BamExtracting");
	setUp.setOption(bamExtractPars_.throwAwayUnmappedMate_, "--throwAwayUnmappedMate", "Throw Away Unmapped Mate", false, "BamExtracting");
	setUp.setOption(bamExtractPars_.percInRegion_,          "--percInRegion",          "Percent of bases in Region to be included", false, "BamExtracting");
	setUp.setOption(bamExtractPars_.minAlnMapSize_, "--minAlnMapSize", "min Aln Map Size for initial recruitment", false, "BamExtracting");
	setUp.setOption(bamExtractPars_.softClipPercentageCutOff_, "--softClipPercentageCutOff", "The minimum percentage of the bases that can be soft clipped to include alignment", false, "BamExtracting");

}


void HaploPathFinder::PathFinderCorePars::setFinalLengthCutOffs(seqSetUp & setUp){
	lenCutOff = 0;
	setUp.setOption(lenCutOff, "--lenCutOff", "Length Cut Off");
	taillessOrHeaddLesslenCutOff = lenCutOff;
	setUp.setOption(taillessOrHeaddLesslenCutOff, "--taillessOrHeaddLesslenCutOff", "Length Cut Off for headless or tailless nodes");
}





Json::Value HaploPathFinder::PathFinderCorePars::toJson() const{
	/**@todo this needs updating with new members 2020/02/19;
	 *
	 */
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	//ret["connectorCutOff"] = njh::json::toJson(connectorCutOff);
	ret["trimEnds"] = njh::json::toJson(trimEnds);
	ret["trimEndsByCoverage"] = njh::json::toJson(trimEndsByCoverage);
	ret["verbose"] = njh::json::toJson(verbose);
	//ret["debug"] = njh::json::toJson(debug);
	ret["graphVerbose_"] = njh::json::toJson(graphVerbose_);
	ret["graphDebug_"] = njh::json::toJson(graphDebug_);
	ret["removeDuplicatedSequences_"] = njh::json::toJson(removeDuplicatedSequences_);
	ret["trimRemoveDuplicatedSequences_"] = njh::json::toJson(trimRemoveDuplicatedSequences_);
	ret["qualTrim_"] = njh::json::toJson(qualTrim_);
	ret["trimOnQual_"] = njh::json::toJson(trimOnQual_);
	ret["trimShortTips_"] = njh::json::toJson(trimShortTips_);
	ret["collapseOneBaseIndelsNodes_"] = njh::json::toJson(collapseOneBaseIndelsNodes_);
	ret["initialBreakingSinglePaths_"] = njh::json::toJson(initialBreakingSinglePaths_);
	ret["splitTailed_"] = njh::json::toJson(splitTailed_);
	ret["splitEnds_"] = njh::json::toJson(splitEnds_);
	ret["splitEndsOnce_"] = njh::json::toJson(splitEndsOnce_);
	ret["keepCycles_"] = njh::json::toJson(keepCycles_);
	ret["addSeqOfSingleHeadAndTailSeqs_"] = njh::json::toJson(addSeqOfSingleHeadAndTailSeqs_);

	ret["collapseFinalDeterminedTandems_"] = njh::json::toJson(collapseFinalDeterminedTandems_);
	ret["finalTandemOccurenceCutOff_"] = njh::json::toJson(finalTandemOccurenceCutOff_);
	ret["addGroups_"] = njh::json::toJson(addGroups_);
	ret["writeGroupInputNames_"] = njh::json::toJson(writeGroupInputNames_);
	ret["removeHeadlessTaillessAlongTheWay"] = njh::json::toJson(removeHeadlessTaillessAlongTheWay);
	ret["headlessTailessLenCutOff"] = njh::json::toJson(headlessTailessLenCutOff);
	ret["collapseFinalClusters_"] = njh::json::toJson(collapseFinalClusters_);
	ret["writeAllPossibleHaps_"] = njh::json::toJson(writeAllPossibleHaps_);
	ret["writeEstimatedMajorHaps_"] = njh::json::toJson(writeEstimatedMajorHaps_);
	ret["numThreads"] = njh::json::toJson(numThreads);
	ret["lenCutOff"] = njh::json::toJson(lenCutOff);
	ret["taillessOrHeaddLesslenCutOff"] = njh::json::toJson(taillessOrHeaddLesslenCutOff);

	ret["trimEndsBy"] = njh::json::toJson(trimEndsBy);
	ret["removePossibleOutliersWithMuscle"] = njh::json::toJson(removePossibleOutliersWithMuscle);
	ret["outliersCutOff"] = njh::json::toJson(outliersCutOff);
	ret["removePossibleOutliersKSim"] = njh::json::toJson(removePossibleOutliersKSim);
	ret["outlierKLen"] = njh::json::toJson(outlierKLen);
	ret["outliersKSimCutOff"] = njh::json::toJson(outliersKSimCutOff);
	ret["keepGroupsForOutlierFilter"] = njh::json::toJson(keepGroupsForOutlierFilter);
	ret["filterOffOutlierInputSeqs"] = njh::json::toJson(filterOffOutlierInputSeqs);
	ret["stitchPairs"] = njh::json::toJson(stitchPairs);
	ret["reOrientPairsForStitching"] = njh::json::toJson(reOrientPairsForStitching);
	//ret["pairProcessorParams_"] = njh::json::toJson(pairProcessorParams_);
	ret["trimToInputSeqs"] = njh::json::toJson(trimToInputSeqs);
	ret["trimWithGlobalAln"] = njh::json::toJson(trimWithGlobalAln);
	ret["splitToRecruit_"] = njh::json::toJson(splitToRecruit_);
	ret["addHeadTailSeqsToRecruit"] = njh::json::toJson(addHeadTailSeqsToRecruit);
	//ret["mtPars"] = njh::json::toJson(mtPars);
	ret["optKLenStart"] = njh::json::toJson(optKLenStart);
	ret["optKLenStop"] = njh::json::toJson(optKLenStop);
	ret["optKLenStep"] = njh::json::toJson(optKLenStep);
	ret["kmerLengths"] = njh::json::toJson(kmerLengths);
	ret["optKcutStart"] = njh::json::toJson(optKcutStart);
	ret["optKcutStop"] = njh::json::toJson(optKcutStop);
	ret["optKcutStep"] = njh::json::toJson(optKcutStep);
	ret["kmerKOcurrenceCutOffs"] = njh::json::toJson(kmerKOcurrenceCutOffs);
	ret["forceFullKcutOptimize_"] = njh::json::toJson(forceFullKcutOptimize_);
	ret["kCutOptimizationAttempts"] = njh::json::toJson(kCutOptimizationAttempts);
	ret["optimizeNodeSizeCutOff_"] = njh::json::toJson(optimizeNodeSizeCutOff_);
	ret["useFullOptimalCount_"] = njh::json::toJson(useFullOptimalCount_);
	ret["optimizeOnIteration"] = njh::json::toJson(optimizeOnIteration);
	ret["keepOptimizedSubDirs"] = njh::json::toJson(keepOptimizedSubDirs);
	ret["percentOfInputReadsUsedForOptimization"] = njh::json::toJson(percentOfInputReadsUsedForOptimization);
	ret["percentOfInputReadsUsedForOptimizationStepDown"] = njh::json::toJson(percentOfInputReadsUsedForOptimizationStepDown);


	ret["optShortTipNumberStart"] = njh::json::toJson(optShortTipNumberStart);
	ret["optShortTipNumberStop"] = njh::json::toJson(optShortTipNumberStop);
	ret["optShortTipNumberStep"] = njh::json::toJson(optShortTipNumberStep);
	ret["shortTipNumbers"] = njh::json::toJson(shortTipNumbers);
	ret["forceFullShortTipOptimize_"] = njh::json::toJson(forceFullShortTipOptimize_);
	ret["shortTipOptimizationAttempts"] = njh::json::toJson(shortTipOptimizationAttempts);

	return ret;
}



void HaploPathFinder::ExtractParams::setGenomeFnp() {
	genomeFnp = njh::files::make_path(genomeDir,
			njh::appendAsNeededRet(primaryGenome, ".fasta"));
}

HaploPathFinder::HaploPathFinder(ExtractParams inPars) :
		pars_(inPars) {
	/**@todo make the input just the genome paramters
	 *
	 */
	//make sure the genome fnp is set;
	pars_.setGenomeFnp();

	//set up genome mapper;
	gMapper_ = std::make_unique<MultiGenomeMapper>(pars_.genomeDir, pars_.primaryGenome);
	//set up selected genomes
	gMapper_->setSelectedGenomes(pars_.selectedGenomes);
	gMapper_->pars_.numThreads_ = pars_.pFinderPars_.numThreads;
	gMapper_->init();


}


}  // namespace njhseq


