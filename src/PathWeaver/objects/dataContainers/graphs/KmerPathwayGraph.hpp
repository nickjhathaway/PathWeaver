#pragma once
/*
 * KmerPathwayGraph.hpp
 *
 *  Created on: Feb 9, 2017
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

#include <njhseq/utils.h>
#include "PathWeaver/PathFinding/HaploPathFinder.hpp"



namespace njhseq {


uint32_t countEndHomopolymerLength(const std::string & seq);



class KmerPathwayGraph {
public:
	KmerPathwayGraph(uint32_t klen);
	KmerPathwayGraph(uint32_t klen, uint32_t occurenceCutOff);
	uint32_t klen_;
	uint32_t occurenceCutOff_;

	uint32_t numThreads_{1};

	double homopolymerIndelCollapseFreqMultiplier_{3};
	comparison allowableErrorForHPIndexCollapse_;

	//bool debug_ = false;
	//bool verbose_ = false;

	double bridgingCutOff_ = .00;



	std::vector<std::string> readNames_;

	void setOccurenceCutOff(uint32_t cutOff);

	class edge;
	class node {
	public:

		node(const std::string & k, uint32_t cnt, uint32_t klen);
		std::string k_;
		std::string uid_;
		uint32_t cnt_;
		uint32_t kLen_;

		uint32_t group_ = std::numeric_limits<uint32_t>::max();

		std::vector<std::shared_ptr<edge>> headEdges_;
		std::vector<std::shared_ptr<edge>> tailEdges_;

		uint32_t visitCount_ = 0;
		bool on_ = true;

		std::unordered_set<uint32_t> inReadNamesIdx_;

		void resetVisitCount();

		void addHead(const std::shared_ptr<edge> & e);

		void addTail(const std::shared_ptr<edge> & e);

		bool headless() const;
		uint32_t headCount() const;
		std::shared_ptr<edge> getFirstOnHeadEdge() const;
		std::shared_ptr<edge> getLastOnHeadEdge() const;

		bool tailless() const;
		uint32_t tailCount() const;
		std::shared_ptr<edge> getFirstOnTailEdge() const;
		std::shared_ptr<edge> getLastOnTailEdge() const;

		std::shared_ptr<node> getHeadNode(const std::string & k) const;

		std::shared_ptr<node> getTailNode(const std::string & k) const;

		void turnOffAllEdges();
	};
	class edge {
	public:
		edge(const std::shared_ptr<node> & head, const std::shared_ptr<node> & tail,
				uint32_t cnt, const std::unordered_set<uint32_t> & readNames);

		std::weak_ptr<node> head_;
		std::weak_ptr<node> tail_;
		uint32_t cnt_;
		std::unordered_set<uint32_t> inReadNamesIdx_;
		bool on_ = true;

		std::string createUid() const;
	};

	std::vector<std::shared_ptr<node>> nodes_;
	std::unordered_map<std::string, uint32_t> nodePositions_;
	std::unordered_map<std::string, uint32_t> kCounts_;

	void populateNodesFromCounts();

	bool hasNode(const std::string & nodeName) const;

	void addNode(const std::string & k, uint32_t cnt, uint32_t kLen);

	void addEdge(const std::string & head, const std::string & tailName,
			double cnt, const uint32_t & readNameIdx);


	void removeOffEdges();
	void turnOffEdgesBelowCutOff(uint32_t countCutOff);
	void turnOffEdgesBelowCutOff();

	void turnOffLowestEdges();

	void turnOffNodesBelowCutOff(uint32_t countCutOff);
	void turnOffNodesBelowCutOff();


	void edgeSanityCheckThrow() const;

	bool hasCycle()const;
	void collapseSingleLinkedPaths(bool initialCollapse = false);
	void collapseSingleLinkedPathsForPossibleLoops();

	void breakSingleLinkedPathsReadThreading();

	bool removeShortTips(uint32_t shortNumber, uint32_t cntCutOff);
	bool removeShortTips_after_disentangleInternalNodes(uint32_t shortNumber, uint32_t cntCutOff);



	void removeOffNodes();
	void removeNullNodes();
	bool removeHeadlessTaillessNodes();
	bool removeHeadlessTaillessNodes(uint32_t lenCutOff);

	bool breakSingleHeadSingleTailNodesLowCoverage();

	struct disentangleInternalNodesPars {
		disentangleInternalNodesPars(bool conservative,
				uint32_t shortTipNumber,
				uint32_t shortTipCutOff,
				uint32_t headlessTaillessCutOff,
				const HaploPathFinder::PathFinderCorePars & extractionPars) :
				conservative_(conservative), shortTipNumber_(shortTipNumber), shortTipCutOff_(
						shortTipCutOff), headlessTaillessCutOff_(headlessTaillessCutOff) {
			trimShortTips_ = extractionPars.trimShortTips_;
			throwAwayConservedAddNodesDuringDisentaglement_ = extractionPars.throwAwayConservedAddNodesDuringDisentaglement;
			removeHeadlessTaillessAfterDisentaglement_ = extractionPars.removeHeadlessTaillessAfterDisentaglement;
		}
		bool conservative_ = false;
		uint32_t shortTipNumber_ { 2 };
		uint32_t shortTipCutOff_ { std::numeric_limits < uint32_t > ::max() };
		uint32_t headlessTaillessCutOff_ {0};
		bool throwAwayConservedAddNodesDuringDisentaglement_{false};
		bool trimShortTips_{false};

		bool removeHeadlessTaillessAfterDisentaglement_{false};
	};

	bool disentangleInternalNodes(const disentangleInternalNodesPars & pars);

	bool breakSelfPointingPathsKeepOtherEdges();
  bool breakSelfPointingPaths();
  bool hasSelfPointingPaths();

	bool splitEndNodes();
	bool splitEndNodes(uint32_t maxLen);

	bool trimEdgesNodeTipsWithLowEntropy(const readVecTrimmer::TrimEdgesByLowEntropyPars & pars);



	bool splitMultitailedNodes();

	void resetNodePositions();
	void resetNodeUids();

	table getNodeStateCounts() const;
	table getTailCounts() const;

	Json::Value genSimpleGraphJson() const;

	void writeDot(std::ostream & out) const;
	void writeRectangleDot(std::ostream & out, bool noLabels = false) const;
	void writeRectangleWithEstimatedCovDot(std::ostream & out, KmerPathwayGraph & estimatingCovGraph, bool noLabels = false) const;
	void writeRectangleWithEstimatedCovDotByGroup(const OutOptions & prefixPath, KmerPathwayGraph & estimatingCovGraph, uint32_t groupSizeCutOff = 2, bool noLabels = false) const;


	static VecStr writeEdgesHeader();
	void writeEdges(std::ostream & out) const;
	static VecStr writeEdgesWithNamesHeader();
	void writeEdgesWithNames(std::ostream & out) const;
	static VecStr writeNodesHeader();
	void writeNodes(std::ostream & out) const;


	void writeOutNodesInGroups(const bfs::path & outDirname, bool overWrite) const;

	void sortNodesBySize();

	//void sortNodes(std::function<bool(const node&,const node &)> func);
	void sortNodes(std::function<bool(const std::shared_ptr<node>&,const std::shared_ptr<node> &)> func);


	void increaseKCounts(const std::string & seq);
	struct StartEndPos{
		StartEndPos(uint32_t start, uint32_t end);
		uint32_t start_;
		uint32_t end_;

		uint32_t len() const;


		static std::vector<StartEndPos> merge(std::vector<StartEndPos> positions, uint32_t expand = 1);
	};



	void increaseKCountsAdjust(const std::string & seq, const std::vector<StartEndPos> & adjustPositions);

	std::unordered_map<std::string, uint32_t> threadThroughSequence(const seqInfo & seq);
	std::unordered_map<std::string, uint32_t> threadThroughSequence(const seqInfo & seq, const std::string & seqName);

	std::unordered_map<std::string, uint32_t> threadThroughSequenceMate(const PairedRead & pSeq,
			std::unordered_map<std::string, uint32_t> & firstMateCounts,
			const std::string & threadingSeqNameFirstMate,
			uint32_t threadingSeqNameFirstMateIdx);


	void threadThroughSequence(const PairedRead & pSeq);

	std::unordered_map<std::string, uint32_t> threadThroughSequenceAdjust(const seqInfo & seq,
			const std::vector<StartEndPos> & adjustPositions,
			const std::string & seqName);

	std::unordered_map<std::string, uint32_t> threadThroughSequenceAdjustMate(const PairedRead & pseq,
			const std::vector<StartEndPos> & adjustPositions,
			std::unordered_map<std::string, uint32_t> & firstMateCounts,
						const std::string & threadingSeqNameFirstMate,
						uint32_t threadingSeqNameFirstMateIdx);

	void threadThroughSequenceAdjust(const PairedRead & pseq,
			const std::vector<StartEndPos> & adjustPositionsFirstMate,
			const std::vector<StartEndPos> & adjustPositionsSecondMate);

	void threadThroughSequenceHelper(const std::string & firstKmer,
			const std::shared_ptr<node> & headNode,
			const std::string & nextKmer,
			const std::shared_ptr<node> & tailNode,
			std::unordered_map<std::string, uint32_t> & internalCount,
			const std::string & threadingSeqName,
			const uint32_t threadingSeqNameIdx);

	void threadThroughSequenceHelperMate(const std::string & firstKmer,
			const std::shared_ptr<node> & headNode,
			const std::string & nextKmer,
			const std::shared_ptr<node> & tailNode,
			std::unordered_map<std::string, uint32_t> & internalCount,
			const std::string & threadingSeqName,
			const uint32_t threadingSeqNameIdx,
			std::unordered_map<std::string, uint32_t> & internalCountMate,
			const std::string & threadingSeqNameMate,
			const uint32_t threadingSeqNameMateIdx);


	bool collapseBubbleNodesWithError(const comparison & errorAllowed, double freqMultiCutOff);


	bool collapseOneBaseIndelsNodes(KmerPathwayGraph & estimatingCovGraph);
	//bool collapseOneBaseIndelsNodes3Nodes(KmerPathwayGraph & estimatingCovGraph);
	bool collapseOneBaseIndelsNodesComplex();

	bool trimTipsOfLowEntropyNodes(double entropyCutOff);


	static bool skipInputSeqForKCountOld(const std::string & seq, uint32_t kLen);

	static bool skipInputSeqForKCount(const std::string & seq, uint32_t kLen);

	std::vector<seqInfo> nodesToSeqs(bool addSeqOfSingleHeadAndTailSeqs = false) const;

	void resetGroups() const;
	void resetGroupsLoopAware() const;





	KmerPathwayGraph copyGraph() const;


	struct PossibleHap {
		PossibleHap();
		PossibleHap(const node & firstNode);
		PossibleHap(const seqInfo & seqBase, const Bed6RecordCore & firstContigLoc);
		seqInfo seqBase_;
		std::vector<Bed6RecordCore> contigLocations_;
		void initWithFirstContig(const node & nextNode);
		void addContig(const node & nextNode, const std::string & edgeUID);
		std::unordered_set<std::string> edgesUids_;
	};

	std::vector<PossibleHap> generateAllPossibleHaps() const;


};



}  // namespace njhseq
