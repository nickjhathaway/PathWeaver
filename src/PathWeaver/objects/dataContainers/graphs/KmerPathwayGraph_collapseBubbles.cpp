/*
 * KmerPathwayGraph_collapseBubbles.cpp
 *
 *  Created on: Feb 27, 2020
 *      Author: nicholashathaway
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




#include "KmerPathwayGraph.hpp"
#include "PathWeaver/PathFinding/CoverageEstimator.hpp"


namespace njhseq {


bool KmerPathwayGraph::collapseBubbleNodesWithError(const comparison & errorAllowed, double freqMultiCutOff){

	//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//frequency of the second node needs to be this multipler times it's own coverage less than the higher node
	//collapse single base errors in homopolymer runs

	std::map<std::string, std::vector<std::shared_ptr<node>>> groupedNodes;
	//gather together the nodes that have the same heads and tails and are also singlely linked forward and backwards
	for(const auto & n : nodes_){
		VecStr tailConnectorNames;
		VecStr headConnectorNames;
		if(( 1 == n->headCount() && 1== n->tailCount() ) ||
				(1 == n->headCount() && n->tailless() )      ||
				(n->headless()       && 1== n->tailCount() ) ){
			for(const auto & h : n->headEdges_){
				if(h->on_){
					headConnectorNames.emplace_back(h->head_.lock()->uid_);
				}
			}
			for(const auto & t : n->tailEdges_){
				if(t->on_){
					tailConnectorNames.emplace_back(t->tail_.lock()->uid_);
				}
			}
		}
		if(!headConnectorNames.empty() || !tailConnectorNames.empty() ){
			njh::sort(headConnectorNames);
			njh::sort(tailConnectorNames);
			groupedNodes[njh::pasteAsStr("head:",njh::conToStr(headConnectorNames, ","), ";tail:", njh::conToStr(tailConnectorNames, ","))].emplace_back(n);
		}else{
			//leave out the headless and tailless nodes
			//groupedNodes[n->uid_].emplace_back(n);
		}
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	std::cout << "groupedNodes.size(): " << groupedNodes.size() << std::endl;
//	struct CovInfoPerPos {
//		CovInfoPerPos(){
//		}
//		CovInfoPerPos(uint32_t pos, uint32_t minPos, double sdCov, double avgCov) :
//				pos_(pos), minPos_(minPos), sdCov_(sdCov), avgCov_(avgCov) {
//		}
//		uint32_t pos_ {std::numeric_limits<uint32_t>::max()};
//		uint32_t minPos_ {std::numeric_limits<uint32_t>::max()};
//		double sdCov_ {std::numeric_limits<double>::max()};
//		double avgCov_ {std::numeric_limits<double>::max()};
//	};
	bool modifiedNodes = false;
//	uint32_t roundingSize = 5;
	//iterate over the groups
	// uint32_t groupCount = 0;
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(auto & group : groupedNodes){

//		std::cout << "\tgroup.second.size(): " << group.second.size() << std::endl;
		// if there are only two nodes in the group, process them
		if(2 == group.second.size()){
			// ++groupCount;
			if(uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length()) > 20){
				continue;
			}
//			std::cout << "groupCount: " << groupCount << std::endl;
//			for(const auto & n : group.second){
//				std::cout << "\t" << n->k_.size() << std::endl;
//			}

			//sort the two nodes so that the front node is the higher read count than the 2nd one
			njh::sort(group.second, [](const std::shared_ptr<node> & node1, const std::shared_ptr<node> & node2){
				if(node1->cnt_ == node2->cnt_){
					return node1->k_.length() < node2->k_.length();
				}else{
					return node1->cnt_ > node2->cnt_;
				}
			});

			/**@todo optimize this coverage comparison, this starts to become non-optimal way of
			 *  estimating coverage when there are one base indel which is bad since that's the whole point of this function lol*/
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			// get the approximate coverage of the two nodes
			double approxPerBaseCoverageNode1 = static_cast<double>(group.second.front()->cnt_)/(group.second.front()->k_.length() - klen_ + 1);
			double approxPerBaseCoverageNode2 = static_cast<double>(group.second.back()->cnt_)/(group.second.back()->k_.length() - klen_ + 1);
//			std::cout << "\tNode1 length: " << group.second.front()->k_.length() << std::endl;
//			std::cout << "\tNode2 length: " << group.second.back()->k_.length() << std::endl;
//			std::cout << "\tNode1 approxPerBaseCoverageNode1: " << approxPerBaseCoverageNode1 << std::endl;
//			std::cout << "\tNode2 approxPerBaseCoverageNode2: " << approxPerBaseCoverageNode2 << std::endl;
//			std::cout << "\tNode1 group.second.front()->cnt_: " << group.second.front()->cnt_ << std::endl;
//			std::cout << "\tNode2 group.second.back()->cnt_ : " << group.second.back()->cnt_ << std::endl;
//			auto estCovNode1 = CoverageEstimator::estimateCov(group.second.front()->k_, MetaDataInName(), estimatingCovGraph);
//			auto estCovNode2 = CoverageEstimator::estimateCov(group.second.back()->k_, MetaDataInName(), estimatingCovGraph);
//			auto estCovNode1 = CoverageEstimator::estimateCov(group.second.front()->k_, MetaDataInName(), *this);
//			auto estCovNode2 = CoverageEstimator::estimateCov(group.second.back()->k_, MetaDataInName(), *this);
//			std::cout << "\tNode1 estCovNode1.minCov_.avgCov_: " << estCovNode1.minCov_.avgCov_ << std::endl;
//			std::cout << "\tNode2 estCovNode2.minCov_.avgCov_: " << estCovNode2.minCov_.avgCov_ << std::endl;
//			if(group.second.front()->k_.length() == 80 && group.second.back()->k_.length() == 81){
//				std::cout << "pos\tcoverage" << std::endl;
//				for(const auto & pos : iter::range(estCovNode1.allCounts_.size())){
//					std::cout << pos << "\t" << estCovNode1.allCounts_[pos] << std::endl;
//				}
//				std::cout << "pos\tcoverage\tminPos\tsdCov\tavgCov\tstart\tend\tallCounts.size()" << "\n";
//				for(const auto & estCovPos : iter::range(estCovNode1.coverageInfos_.size())){
//					auto pos = estCovNode1.coverageInfos_[estCovPos].pos_;
//					std::cout << pos
//							<< "\t" << estCovNode1.allCounts_[pos]
//							<< "\t" << estCovNode1.coverageInfos_[estCovPos].minPos_
//							<< "\t" << estCovNode1.coverageInfos_[estCovPos].sdCov_
//							<< "\t" << estCovNode1.coverageInfos_[estCovPos].avgCov_
//							<< "\t" << (pos + 1 > estimatingCovGraph.klen_ ? pos + 1 - estimatingCovGraph.klen_: 0)
//							<< "\t" << pos + 1
//							<< "\t" << estCovNode1.allCounts_.size() << std::endl;
//				}
//				std::cout << "pos\tcoverage" << std::endl;
//				for(const auto & pos : iter::range(estCovNode2.allCounts_.size())){
//					std::cout << pos << "\t" << estCovNode2.allCounts_[pos] << std::endl;
//				}
//				std::cout << "pos\tcoverage\tminPos\tsdCov\tavgCov\tstart\tend\tallCounts.size()" << "\n";
//				for(const auto & estCovPos : iter::range(estCovNode2.coverageInfos_.size())){
//					auto pos = estCovNode2.coverageInfos_[estCovPos].pos_;
//					std::cout << pos
//							<< "\t" << estCovNode2.allCounts_[pos]
//							<< "\t" << estCovNode2.coverageInfos_[estCovPos].minPos_
//							<< "\t" << estCovNode2.coverageInfos_[estCovPos].sdCov_
//							<< "\t" << estCovNode2.coverageInfos_[estCovPos].avgCov_
//							<< "\t" << (pos + 1 > estimatingCovGraph.klen_ ? pos + 1 - estimatingCovGraph.klen_: 0)
//							<< "\t" << pos + 1
//							<< "\t" << estCovNode2.allCounts_.size() << std::endl;
//				}
//			}
//
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << "group.second.front()->k_: " << group.second.front()->k_ << std::endl;
//			std::cout << "group.second.back()->k_: " << group.second.back()->k_ << std::endl;
//
//			std::cout << njh::bashCT::cyan;
//			std::cout << "\tapproxPerBaseCoverageNode1: " << approxPerBaseCoverageNode1<< std::endl;
//			std::cout << "\tapproxPerBaseCoverageNode2: " << approxPerBaseCoverageNode2<< std::endl;
//			std::cout << "\tapproxPerBaseCoverageNode1/approxPerBaseCoverageNode2: " << approxPerBaseCoverageNode1/approxPerBaseCoverageNode2<< std::endl;
//			std::cout << njh::bashCT::reset;
			if(group.second.front()->k_.length() > 2 * klen_ && group.second.back()->k_.length() > 2 * klen_){
				//node1
				auto estCovNode1 = CoverageEstimator::estimateCov(group.second.front()->k_, MetaDataInName(), *this);
				approxPerBaseCoverageNode1 = estCovNode1.minCov_.avgCov_;
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
			}
			if(group.second.front()->k_.length() > 2 * klen_ && group.second.back()->k_.length() > 2 * klen_){
				//node2
				auto estCovNode2 = CoverageEstimator::estimateCov(group.second.back()->k_, MetaDataInName(), *this);
				approxPerBaseCoverageNode2 = estCovNode2.minCov_.avgCov_;
//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
			}
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			std::cout << njh::bashCT::red;
//			std::cout << "\tapproxPerBaseCoverageNode1: " << approxPerBaseCoverageNode1<< std::endl;
//			std::cout << "\tapproxPerBaseCoverageNode2: " << approxPerBaseCoverageNode2<< std::endl;
//			std::cout << "\tapproxPerBaseCoverageNode1/approxPerBaseCoverageNode2: " << approxPerBaseCoverageNode1/approxPerBaseCoverageNode2<< std::endl;
//			std::cout << njh::bashCT::reset;
//			std::cout << std::endl;
			//if the coverage of the second node is approximately homopolymerIndelCollapseFreqMultiplier_ times less than node one, compare them
			if(approxPerBaseCoverageNode1/approxPerBaseCoverageNode2 >= freqMultiCutOff){
//				std::cout << "\t" << "uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length()) : " << uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length())  << std::endl;
				if(uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length()) <= 20){
	//					std::cout << "std::max(group.second.front()->k_.length(), group.second.back()->k_.length()): " << std::max(group.second.front()->k_.length(), group.second.back()->k_.length()) << std::endl;
					//compare the two nodes
					aligner alignerObj(std::max(group.second.front()->k_.length(), group.second.back()->k_.length()), gapScoringParameters(5,1,5,1,5,1));
					//aligner alignerObj(std::max(group.second.front()->k_.length(), group.second.back()->k_.length()), gapScoringParameters(5,1,0,0,0,0));

					alignerObj.countEndGaps_ = true;
					alignerObj.weighHomopolymers_ = true;
					alignerObj.alignRegGlobal(seqInfo("1", group.second.front()->k_),
							seqInfo("2", group.second.back()->k_));
					alignerObj.profileAlignment(
							seqInfo("1", group.second.front()->k_),
							seqInfo("2", group.second.back()->k_), false, false, false);
	//					//uncomment for debugging purposes
	//					alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
	//					alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
	//					std::cout << "allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_): " << njh::colorBool(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)) << std::endl;
					if(errorAllowed.passErrorProfile(alignerObj.comp_)){
						modifiedNodes = true;
						//add in the other node's read names
						group.second.front()->inReadNamesIdx_.insert(
								group.second.back()->inReadNamesIdx_.begin(),
								group.second.back()->inReadNamesIdx_.end());
						if(!group.second.front()->headless()){
	//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//							std::cout << "\tgroup.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.size(): " << group.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.size() << std::endl;
							group.second.front()->getFirstOnHeadEdge()->cnt_ += group.second.back()->getFirstOnHeadEdge()->cnt_;
							group.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.insert(
									group.second.back()->getFirstOnHeadEdge()->inReadNamesIdx_.begin(),
									group.second.back()->getFirstOnHeadEdge()->inReadNamesIdx_.end());
	//							std::cout << "\tgroup.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.size(): " << group.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.size() << std::endl;
						}
						if(!group.second.front()->tailless()){
	//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//							std::cout << "\tgroup.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.size(): " << group.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.size() << std::endl;
							group.second.front()->getFirstOnTailEdge()->cnt_ += group.second.back()->getFirstOnTailEdge()->cnt_;
							group.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.insert(
									group.second.back()->getFirstOnTailEdge()->inReadNamesIdx_.begin(),
									group.second.back()->getFirstOnTailEdge()->inReadNamesIdx_.end());
	//							std::cout << "\tgroup.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.size(): " << group.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.size() << std::endl;

						}
						// consider adding the count of the other node
						// now turn off the node and it's edges
						group.second.back()->on_ = false;
						if(!group.second.back()->headless()){
							group.second.back()->getFirstOnHeadEdge()->on_ = false;
						}
						if(!group.second.back()->tailless()){
							group.second.back()->getFirstOnTailEdge()->on_ = false;
						}
					}
				}
			}
		}
	}
//	//
	if (modifiedNodes) {
		//remove off nodes;
		removeOffNodes();
		removeOffEdges();
		resetNodePositions();
		return true;
	}
	return false;
}



bool KmerPathwayGraph::collapseOneBaseIndelsNodes(KmerPathwayGraph & estimatingCovGraph){
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	//frequency of the second node needs to be this multipler times it's own coverage less than the higher node
	//collapse single base errors in homopolymer runs

	std::map<std::string, std::vector<std::shared_ptr<node>>> groupedNodes;
  //gather together the nodes that have the same heads and tails and are also singlely linked forward and backwards
	for(const auto & n : nodes_){
		VecStr tailConnectorNames;
		VecStr headConnectorNames;
		if(( 1 == n->headCount() && 1== n->tailCount() ) ||
				(1 == n->headCount() && n->tailless() )      ||
				(n->headless()       && 1== n->tailCount() ) ){
			for(const auto & h : n->headEdges_){
				if(h->on_){
					headConnectorNames.emplace_back(h->head_.lock()->uid_);
				}
			}
			for(const auto & t : n->tailEdges_){
				if(t->on_){
					tailConnectorNames.emplace_back(t->tail_.lock()->uid_);
				}
			}
		}
		if(!headConnectorNames.empty() || !tailConnectorNames.empty() ){
			njh::sort(headConnectorNames);
			njh::sort(tailConnectorNames);
			groupedNodes[njh::pasteAsStr("head:",njh::conToStr(headConnectorNames, ","), ";tail:", njh::conToStr(tailConnectorNames, ","))].emplace_back(n);
		}else{
			//leave out the headless and tailless nodes
			//groupedNodes[n->uid_].emplace_back(n);
		}
	}
//	std::cout << "groupedNodes.size(): " << groupedNodes.size() << std::endl;
//	struct CovInfoPerPos {
//		CovInfoPerPos(){
//		}
//		CovInfoPerPos(uint32_t pos, uint32_t minPos, double sdCov, double avgCov) :
//				pos_(pos), minPos_(minPos), sdCov_(sdCov), avgCov_(avgCov) {
//		}
//		uint32_t pos_ {std::numeric_limits<uint32_t>::max()};
//		uint32_t minPos_ {std::numeric_limits<uint32_t>::max()};
//		double sdCov_ {std::numeric_limits<double>::max()};
//		double avgCov_ {std::numeric_limits<double>::max()};
//	};
	bool modifiedNodes = false;
//	uint32_t roundingSize = 5;
	//iterate over the groups
	// uint32_t groupCount = 0;
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for(auto & group : groupedNodes){

//		std::cout << "\tgroup.second.size(): " << group.second.size() << std::endl;
		// if there are only two nodes in the group, process them
		if(3 == group.second.size() || 2 == group.second.size()){
			// ++groupCount;

//			std::cout << "groupCount: " << groupCount << std::endl;
//			for(const auto & n : group.second){
//				std::cout << "\t" << n->k_.size() << std::endl;
//			}

			//sort the two nodes so that the front node is the higher read count than the 2nd one
			njh::sort(group.second, [](const std::shared_ptr<node> & node1, const std::shared_ptr<node> & node2){
				if(node1->cnt_ == node2->cnt_){
					return node1->k_.length() < node2->k_.length();
				}else{
					return node1->cnt_ > node2->cnt_;
				}
			});
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			for(const auto & node : group.second){
//				std::cout << node->k_ << std::endl;
//				std::cout << node->k_.length() << std::endl;
//				std::cout << node->cnt_ << std::endl;
//				std::cout << std::endl;
//			}

			/**@todo optimize this coverage comparison, this starts to become non-optimal way of
			 *  estimating coverage when there are one base indel which is bad since that's the whole point of this function lol*/
//			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			// get the approximate coverage of the two nodes
			for(uint32_t nodePos = 1; nodePos < group.second.size(); ++nodePos){
				if(uAbsdiff(group.second[0]->k_.length(), group.second[1]->k_.length()) > 10){
					continue;
				}
				double approxPerBaseCoverageNode1 = static_cast<double>(group.second[0]->cnt_)/(group.second[0]->k_.length() - klen_ + 1);
				double approxPerBaseCoverageNode2 = static_cast<double>(group.second[nodePos]->cnt_)/(group.second[nodePos]->k_.length() - klen_ + 1);

	//			std::cout << "\tNode1 length: " << group.second[0]->k_.length() << std::endl;
	//			std::cout << "\tNode2 length: " << group.second[nodePos]->k_.length() << std::endl;
	//			std::cout << "\tNode1 approxPerBaseCoverageNode1: " << approxPerBaseCoverageNode1 << std::endl;
	//			std::cout << "\tNode2 approxPerBaseCoverageNode2: " << approxPerBaseCoverageNode2 << std::endl;
	//			std::cout << "\tNode1 group.second[0]->cnt_: " << group.second[0]->cnt_ << std::endl;
	//			std::cout << "\tNode2 group.second[nodePos]->cnt_ : " << group.second[nodePos]->cnt_ << std::endl;
	//			auto estCovNode1 = CoverageEstimator::estimateCov(group.second[0]->k_, MetaDataInName(), estimatingCovGraph);
	//			auto estCovNode2 = CoverageEstimator::estimateCov(group.second[nodePos]->k_, MetaDataInName(), estimatingCovGraph);
	//			auto estCovNode1 = CoverageEstimator::estimateCov(group.second[0]->k_, MetaDataInName(), *this);
	//			auto estCovNode2 = CoverageEstimator::estimateCov(group.second[nodePos]->k_, MetaDataInName(), *this);
	//			std::cout << "\tNode1 estCovNode1.minCov_.avgCov_: " << estCovNode1.minCov_.avgCov_ << std::endl;
	//			std::cout << "\tNode2 estCovNode2.minCov_.avgCov_: " << estCovNode2.minCov_.avgCov_ << std::endl;
	//			if(group.second[0]->k_.length() == 80 && group.second[nodePos]->k_.length() == 81){
	//				std::cout << "pos\tcoverage" << std::endl;
	//				for(const auto & pos : iter::range(estCovNode1.allCounts_.size())){
	//					std::cout << pos << "\t" << estCovNode1.allCounts_[pos] << std::endl;
	//				}
	//				std::cout << "pos\tcoverage\tminPos\tsdCov\tavgCov\tstart\tend\tallCounts.size()" << "\n";
	//				for(const auto & estCovPos : iter::range(estCovNode1.coverageInfos_.size())){
	//					auto pos = estCovNode1.coverageInfos_[estCovPos].pos_;
	//					std::cout << pos
	//							<< "\t" << estCovNode1.allCounts_[pos]
	//							<< "\t" << estCovNode1.coverageInfos_[estCovPos].minPos_
	//							<< "\t" << estCovNode1.coverageInfos_[estCovPos].sdCov_
	//							<< "\t" << estCovNode1.coverageInfos_[estCovPos].avgCov_
	//							<< "\t" << (pos + 1 > estimatingCovGraph.klen_ ? pos + 1 - estimatingCovGraph.klen_: 0)
	//							<< "\t" << pos + 1
	//							<< "\t" << estCovNode1.allCounts_.size() << std::endl;
	//				}
	//				std::cout << "pos\tcoverage" << std::endl;
	//				for(const auto & pos : iter::range(estCovNode2.allCounts_.size())){
	//					std::cout << pos << "\t" << estCovNode2.allCounts_[pos] << std::endl;
	//				}
	//				std::cout << "pos\tcoverage\tminPos\tsdCov\tavgCov\tstart\tend\tallCounts.size()" << "\n";
	//				for(const auto & estCovPos : iter::range(estCovNode2.coverageInfos_.size())){
	//					auto pos = estCovNode2.coverageInfos_[estCovPos].pos_;
	//					std::cout << pos
	//							<< "\t" << estCovNode2.allCounts_[pos]
	//							<< "\t" << estCovNode2.coverageInfos_[estCovPos].minPos_
	//							<< "\t" << estCovNode2.coverageInfos_[estCovPos].sdCov_
	//							<< "\t" << estCovNode2.coverageInfos_[estCovPos].avgCov_
	//							<< "\t" << (pos + 1 > estimatingCovGraph.klen_ ? pos + 1 - estimatingCovGraph.klen_: 0)
	//							<< "\t" << pos + 1
	//							<< "\t" << estCovNode2.allCounts_.size() << std::endl;
	//				}
	//			}
	//
	//			std::cout << __FILE__ << " " << __LINE__ << std::endl;

				if(group.second[0]->k_.length() > 2 * klen_ && group.second[nodePos]->k_.length() > 2 * klen_){
					//node1
					auto estCovNode1 = CoverageEstimator::estimateCov(group.second[0]->k_, MetaDataInName(), *this);
					approxPerBaseCoverageNode1 = estCovNode1.minCov_.avgCov_;
	//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				}
				if(group.second[0]->k_.length() > 2 * klen_ && group.second[nodePos]->k_.length() > 2 * klen_){
					//node2
					auto estCovNode2 = CoverageEstimator::estimateCov(group.second[nodePos]->k_, MetaDataInName(), *this);
					approxPerBaseCoverageNode2 = estCovNode2.minCov_.avgCov_;
	//				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				}

//
//				std::cout << "\tapproxPerBaseCoverageNode1: " << approxPerBaseCoverageNode1<< std::endl;
//				std::cout << "\tapproxPerBaseCoverageNode2: " << approxPerBaseCoverageNode2<< std::endl;
//				std::cout << "\tapproxPerBaseCoverageNode1/approxPerBaseCoverageNode2: " << approxPerBaseCoverageNode1/approxPerBaseCoverageNode2<< std::endl;
				//if the coverage of the second node is approximately homopolymerIndelCollapseFreqMultiplier_ times less than node one, compare them
				if(approxPerBaseCoverageNode1/approxPerBaseCoverageNode2 >= homopolymerIndelCollapseFreqMultiplier_){
					if(uAbsdiff(group.second[0]->k_.length(), group.second[nodePos]->k_.length()) <= 10){

//						std::cout << "\t" << "uAbsdiff(group.second[0]->k_.length(), group.second[nodePos]->k_.length()) : " << uAbsdiff(group.second[0]->k_.length(), group.second[nodePos]->k_.length())  << std::endl;
						//if(uAbsdiff(group.second[0]->k_.length(), group.second[nodePos]->k_.length()) == 1){}

//											std::cout << "\tstd::max(group.second[0]->k_.length(), group.second[nodePos]->k_.length()): " << std::max(group.second[0]->k_.length(), group.second[nodePos]->k_.length()) << std::endl;
						//compare the two nodes
						aligner alignerObj(std::max(group.second[0]->k_.length(), group.second[nodePos]->k_.length()), gapScoringParameters(3,1,3,1,3,1));
						//aligner alignerObj(std::max(group.second[0]->k_.length(), group.second[nodePos]->k_.length()), gapScoringParameters(5,1,0,0,0,0));

						alignerObj.countEndGaps_ = true;
						alignerObj.weighHomopolymers_ = true;
						alignerObj.alignRegGlobal(seqInfo("1", group.second[0]->k_),
								seqInfo("2", group.second[nodePos]->k_));
						alignerObj.profileAlignment(
								seqInfo("1", group.second[0]->k_),
								seqInfo("2", group.second[nodePos]->k_), false, false, false);
		//					//uncomment for debugging purposes
//							alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//							alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
//							std::cout << "allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_): " << njh::colorBool(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)) << std::endl;
						if(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)){
							modifiedNodes = true;
							//add in the other node's read names
							group.second[0]->inReadNamesIdx_.insert(
									group.second[nodePos]->inReadNamesIdx_.begin(),
									group.second[nodePos]->inReadNamesIdx_.end());
							if(!group.second[0]->headless()){
		//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//							std::cout << "\tgroup.second[0]->getFirstOnHeadEdge()->inReadNamesIdx_.size(): " << group.second[0]->getFirstOnHeadEdge()->inReadNamesIdx_.size() << std::endl;
								group.second[0]->getFirstOnHeadEdge()->cnt_ += group.second[nodePos]->getFirstOnHeadEdge()->cnt_;
								group.second[0]->getFirstOnHeadEdge()->inReadNamesIdx_.insert(
										group.second[nodePos]->getFirstOnHeadEdge()->inReadNamesIdx_.begin(),
										group.second[nodePos]->getFirstOnHeadEdge()->inReadNamesIdx_.end());
		//							std::cout << "\tgroup.second[0]->getFirstOnHeadEdge()->inReadNamesIdx_.size(): " << group.second[0]->getFirstOnHeadEdge()->inReadNamesIdx_.size() << std::endl;
							}
							if(!group.second[0]->tailless()){
		//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
		//							std::cout << "\tgroup.second[0]->getFirstOnTailEdge()->inReadNamesIdx_.size(): " << group.second[0]->getFirstOnTailEdge()->inReadNamesIdx_.size() << std::endl;
								group.second[0]->getFirstOnTailEdge()->cnt_ += group.second[nodePos]->getFirstOnTailEdge()->cnt_;
								group.second[0]->getFirstOnTailEdge()->inReadNamesIdx_.insert(
										group.second[nodePos]->getFirstOnTailEdge()->inReadNamesIdx_.begin(),
										group.second[nodePos]->getFirstOnTailEdge()->inReadNamesIdx_.end());
		//							std::cout << "\tgroup.second[0]->getFirstOnTailEdge()->inReadNamesIdx_.size(): " << group.second[0]->getFirstOnTailEdge()->inReadNamesIdx_.size() << std::endl;

							}
							// consider adding the count of the other node
							// now turn off the node and it's edges
							group.second[nodePos]->on_ = false;
							if(!group.second[nodePos]->headless()){
								group.second[nodePos]->getFirstOnHeadEdge()->on_ = false;
							}
							if(!group.second[nodePos]->tailless()){
								group.second[nodePos]->getFirstOnTailEdge()->on_ = false;
							}
						}
					}
				}
			}
		}
	}
//	//
	if (modifiedNodes) {
		//remove off nodes;
		removeOffNodes();
		removeOffEdges();
		resetNodePositions();
		return true;
	}
	return false;
}

//bool KmerPathwayGraph::collapseOneBaseIndelsNodes(KmerPathwayGraph & estimatingCovGraph){
////	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	//frequency of the second node needs to be this multipler times it's own coverage less than the higher node
//	//collapse single base errors in homopolymer runs
//
//	std::map<std::string, std::vector<std::shared_ptr<node>>> groupedNodes;
//  //gather together the nodes that have the same heads and tails and are also singlely linked forward and backwards
//	for(const auto & n : nodes_){
//		VecStr tailConnectorNames;
//		VecStr headConnectorNames;
//		if(( 1 == n->headCount() && 1== n->tailCount() ) ||
//				(1 == n->headCount() && n->tailless() )      ||
//				(n->headless()       && 1== n->tailCount() ) ){
//			for(const auto & h : n->headEdges_){
//				if(h->on_){
//					headConnectorNames.emplace_back(h->head_.lock()->uid_);
//				}
//			}
//			for(const auto & t : n->tailEdges_){
//				if(t->on_){
//					tailConnectorNames.emplace_back(t->tail_.lock()->uid_);
//				}
//			}
//		}
//		if(!headConnectorNames.empty() || !tailConnectorNames.empty() ){
//			njh::sort(headConnectorNames);
//			njh::sort(tailConnectorNames);
//			groupedNodes[njh::pasteAsStr("head:",njh::conToStr(headConnectorNames, ","), ";tail:", njh::conToStr(tailConnectorNames, ","))].emplace_back(n);
//		}else{
//			//leave out the headless and tailless nodes
//			//groupedNodes[n->uid_].emplace_back(n);
//		}
//	}
////	std::cout << "groupedNodes.size(): " << groupedNodes.size() << std::endl;
////	struct CovInfoPerPos {
////		CovInfoPerPos(){
////		}
////		CovInfoPerPos(uint32_t pos, uint32_t minPos, double sdCov, double avgCov) :
////				pos_(pos), minPos_(minPos), sdCov_(sdCov), avgCov_(avgCov) {
////		}
////		uint32_t pos_ {std::numeric_limits<uint32_t>::max()};
////		uint32_t minPos_ {std::numeric_limits<uint32_t>::max()};
////		double sdCov_ {std::numeric_limits<double>::max()};
////		double avgCov_ {std::numeric_limits<double>::max()};
////	};
//	bool modifiedNodes = false;
////	uint32_t roundingSize = 5;
//	//iterate over the groups
//	uint32_t groupCount = 0;
////	std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	for(auto & group : groupedNodes){
//
////		std::cout << "\tgroup.second.size(): " << group.second.size() << std::endl;
//		// if there are only two nodes in the group, process them
//		if(2 == group.second.size()){
//			++groupCount;
//			if(uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length()) > 10){
//				continue;
//			}
//
////			std::cout << "groupCount: " << groupCount << std::endl;
////			for(const auto & n : group.second){
////				std::cout << "\t" << n->k_.size() << std::endl;
////			}
//
//			//sort the two nodes so that the front node is the higher read count than the 2nd one
//			njh::sort(group.second, [](const std::shared_ptr<node> & node1, const std::shared_ptr<node> & node2){
//				if(node1->cnt_ == node2->cnt_){
//					return node1->k_.length() < node2->k_.length();
//				}else{
//					return node1->cnt_ > node2->cnt_;
//				}
//			});
//
//			/**@todo optimize this coverage comparison, this starts to become non-optimal way of
//			 *  estimating coverage when there are one base indel which is bad since that's the whole point of this function lol*/
////			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			// get the approximate coverage of the two nodes
//			double approxPerBaseCoverageNode1 = static_cast<double>(group.second.front()->cnt_)/(group.second.front()->k_.length() - klen_ + 1);
//			double approxPerBaseCoverageNode2 = static_cast<double>(group.second.back()->cnt_)/(group.second.back()->k_.length() - klen_ + 1);
////			std::cout << "\tNode1 length: " << group.second.front()->k_.length() << std::endl;
////			std::cout << "\tNode2 length: " << group.second.back()->k_.length() << std::endl;
////			std::cout << "\tNode1 approxPerBaseCoverageNode1: " << approxPerBaseCoverageNode1 << std::endl;
////			std::cout << "\tNode2 approxPerBaseCoverageNode2: " << approxPerBaseCoverageNode2 << std::endl;
////			std::cout << "\tNode1 group.second.front()->cnt_: " << group.second.front()->cnt_ << std::endl;
////			std::cout << "\tNode2 group.second.back()->cnt_ : " << group.second.back()->cnt_ << std::endl;
////			auto estCovNode1 = CoverageEstimator::estimateCov(group.second.front()->k_, MetaDataInName(), estimatingCovGraph);
////			auto estCovNode2 = CoverageEstimator::estimateCov(group.second.back()->k_, MetaDataInName(), estimatingCovGraph);
////			auto estCovNode1 = CoverageEstimator::estimateCov(group.second.front()->k_, MetaDataInName(), *this);
////			auto estCovNode2 = CoverageEstimator::estimateCov(group.second.back()->k_, MetaDataInName(), *this);
////			std::cout << "\tNode1 estCovNode1.minCov_.avgCov_: " << estCovNode1.minCov_.avgCov_ << std::endl;
////			std::cout << "\tNode2 estCovNode2.minCov_.avgCov_: " << estCovNode2.minCov_.avgCov_ << std::endl;
////			if(group.second.front()->k_.length() == 80 && group.second.back()->k_.length() == 81){
////				std::cout << "pos\tcoverage" << std::endl;
////				for(const auto & pos : iter::range(estCovNode1.allCounts_.size())){
////					std::cout << pos << "\t" << estCovNode1.allCounts_[pos] << std::endl;
////				}
////				std::cout << "pos\tcoverage\tminPos\tsdCov\tavgCov\tstart\tend\tallCounts.size()" << "\n";
////				for(const auto & estCovPos : iter::range(estCovNode1.coverageInfos_.size())){
////					auto pos = estCovNode1.coverageInfos_[estCovPos].pos_;
////					std::cout << pos
////							<< "\t" << estCovNode1.allCounts_[pos]
////							<< "\t" << estCovNode1.coverageInfos_[estCovPos].minPos_
////							<< "\t" << estCovNode1.coverageInfos_[estCovPos].sdCov_
////							<< "\t" << estCovNode1.coverageInfos_[estCovPos].avgCov_
////							<< "\t" << (pos + 1 > estimatingCovGraph.klen_ ? pos + 1 - estimatingCovGraph.klen_: 0)
////							<< "\t" << pos + 1
////							<< "\t" << estCovNode1.allCounts_.size() << std::endl;
////				}
////				std::cout << "pos\tcoverage" << std::endl;
////				for(const auto & pos : iter::range(estCovNode2.allCounts_.size())){
////					std::cout << pos << "\t" << estCovNode2.allCounts_[pos] << std::endl;
////				}
////				std::cout << "pos\tcoverage\tminPos\tsdCov\tavgCov\tstart\tend\tallCounts.size()" << "\n";
////				for(const auto & estCovPos : iter::range(estCovNode2.coverageInfos_.size())){
////					auto pos = estCovNode2.coverageInfos_[estCovPos].pos_;
////					std::cout << pos
////							<< "\t" << estCovNode2.allCounts_[pos]
////							<< "\t" << estCovNode2.coverageInfos_[estCovPos].minPos_
////							<< "\t" << estCovNode2.coverageInfos_[estCovPos].sdCov_
////							<< "\t" << estCovNode2.coverageInfos_[estCovPos].avgCov_
////							<< "\t" << (pos + 1 > estimatingCovGraph.klen_ ? pos + 1 - estimatingCovGraph.klen_: 0)
////							<< "\t" << pos + 1
////							<< "\t" << estCovNode2.allCounts_.size() << std::endl;
////				}
////			}
////
////			std::cout << __FILE__ << " " << __LINE__ << std::endl;
//
//			if(group.second.front()->k_.length() > 2 * klen_ && group.second.back()->k_.length() > 2 * klen_){
//				//node1
//				auto estCovNode1 = CoverageEstimator::estimateCov(group.second.front()->k_, MetaDataInName(), *this);
//				approxPerBaseCoverageNode1 = estCovNode1.minCov_.avgCov_;
////				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			}
//			if(group.second.front()->k_.length() > 2 * klen_ && group.second.back()->k_.length() > 2 * klen_){
//				//node2
//				auto estCovNode2 = CoverageEstimator::estimateCov(group.second.back()->k_, MetaDataInName(), *this);
//				approxPerBaseCoverageNode2 = estCovNode2.minCov_.avgCov_;
////				std::cout << __FILE__ << " " << __LINE__ << std::endl;
//			}
//
////			std::cout << "\tapproxPerBaseCoverageNode1: " << approxPerBaseCoverageNode1<< std::endl;
////			std::cout << "\tapproxPerBaseCoverageNode2: " << approxPerBaseCoverageNode2<< std::endl;
////			std::cout << "\tapproxPerBaseCoverageNode1/approxPerBaseCoverageNode2: " << approxPerBaseCoverageNode1/approxPerBaseCoverageNode2<< std::endl;
//			//if the coverage of the second node is approximately homopolymerIndelCollapseFreqMultiplier_ times less than node one, compare them
//			if(approxPerBaseCoverageNode1/approxPerBaseCoverageNode2 >= homopolymerIndelCollapseFreqMultiplier_){
//				if(uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length()) <= 10){
//
//	//				std::cout << "\t" << "uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length()) : " << uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length())  << std::endl;
//					//if(uAbsdiff(group.second.front()->k_.length(), group.second.back()->k_.length()) == 1){}
//
//					//					std::cout << "std::max(group.second.front()->k_.length(), group.second.back()->k_.length()): " << std::max(group.second.front()->k_.length(), group.second.back()->k_.length()) << std::endl;
//					//compare the two nodes
//					aligner alignerObj(std::max(group.second.front()->k_.length(), group.second.back()->k_.length()), gapScoringParameters(5,1,5,1,5,1));
//					//aligner alignerObj(std::max(group.second.front()->k_.length(), group.second.back()->k_.length()), gapScoringParameters(5,1,0,0,0,0));
//
//					alignerObj.countEndGaps_ = true;
//					alignerObj.weighHomopolymers_ = true;
//					alignerObj.alignRegGlobal(seqInfo("1", group.second.front()->k_),
//							seqInfo("2", group.second.back()->k_));
//					alignerObj.profileAlignment(
//							seqInfo("1", group.second.front()->k_),
//							seqInfo("2", group.second.back()->k_), false, false, false);
//	//					//uncomment for debugging purposes
//	//					alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
//	//					alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
//	//					std::cout << "allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_): " << njh::colorBool(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)) << std::endl;
//					if(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)){
//						modifiedNodes = true;
//						//add in the other node's read names
//						group.second.front()->inReadNamesIdx_.insert(
//								group.second.back()->inReadNamesIdx_.begin(),
//								group.second.back()->inReadNamesIdx_.end());
//						if(!group.second.front()->headless()){
//	//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	//							std::cout << "\tgroup.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.size(): " << group.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.size() << std::endl;
//							group.second.front()->getFirstOnHeadEdge()->cnt_ += group.second.back()->getFirstOnHeadEdge()->cnt_;
//							group.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.insert(
//									group.second.back()->getFirstOnHeadEdge()->inReadNamesIdx_.begin(),
//									group.second.back()->getFirstOnHeadEdge()->inReadNamesIdx_.end());
//	//							std::cout << "\tgroup.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.size(): " << group.second.front()->getFirstOnHeadEdge()->inReadNamesIdx_.size() << std::endl;
//						}
//						if(!group.second.front()->tailless()){
//	//							std::cout << __FILE__ << " " << __LINE__ << std::endl;
//	//							std::cout << "\tgroup.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.size(): " << group.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.size() << std::endl;
//							group.second.front()->getFirstOnTailEdge()->cnt_ += group.second.back()->getFirstOnTailEdge()->cnt_;
//							group.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.insert(
//									group.second.back()->getFirstOnTailEdge()->inReadNamesIdx_.begin(),
//									group.second.back()->getFirstOnTailEdge()->inReadNamesIdx_.end());
//	//							std::cout << "\tgroup.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.size(): " << group.second.front()->getFirstOnTailEdge()->inReadNamesIdx_.size() << std::endl;
//
//						}
//						// consider adding the count of the other node
//						// now turn off the node and it's edges
//						group.second.back()->on_ = false;
//						if(!group.second.back()->headless()){
//							group.second.back()->getFirstOnHeadEdge()->on_ = false;
//						}
//						if(!group.second.back()->tailless()){
//							group.second.back()->getFirstOnTailEdge()->on_ = false;
//						}
//					}
//				}
//			}
//		}
//	}
////	//
//	if (modifiedNodes) {
//		//remove off nodes;
//		removeOffNodes();
//		removeOffEdges();
//		resetNodePositions();
//		return true;
//	}
//	return false;
//}
//


bool KmerPathwayGraph::collapseOneBaseIndelsNodesComplex(){
		removeOffNodes();
		removeOffEdges();
		resetNodePositions();
	bool modifiedNodes = false;
	std::vector<std::shared_ptr<node>> nodesToProcess;
	for(const auto & n : nodes_){
		if(2 == n->tailCount() && countEndHomopolymerLength(n->k_) > 7){
			auto firstTailEdge = n->getFirstOnTailEdge();
			auto lastTailEdge = n->getLastOnTailEdge();
			if(2 == firstTailEdge->tail_.lock()->tailCount() &&
				 2 == lastTailEdge->tail_.lock()->tailCount() &&
				 1 == firstTailEdge->tail_.lock()->headCount() &&
				 1 == lastTailEdge->tail_.lock()->headCount() ){
				if( 1 == firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->tailCount() &&
						1 == firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->tailCount() &&
						1 == lastTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->tailCount() &&
						1 == lastTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->tailCount() &&
						1 == firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->headCount() &&
						1 == firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->headCount() &&
						1 == lastTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->headCount() &&
						1 == lastTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->headCount() ){
					std::vector<std::string> tailsFront;
					tailsFront.emplace_back(firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
					tailsFront.emplace_back(firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
					std::vector<std::string> tailsBack;
					tailsBack.emplace_back(lastTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
					tailsBack.emplace_back(lastTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
					njh::sort(tailsFront);
					njh::sort(tailsBack);
					if(std::equal(tailsFront.begin(), tailsFront.end(), tailsBack.begin())){
						nodesToProcess.emplace_back(n);
					}
				}
			}
		}
	}
	if(!nodesToProcess.empty()){
		std::set<std::string> uidsModified;
		for(const auto & n : nodesToProcess){
			auto firstTailEdge = n->getFirstOnTailEdge();
			auto lastTailEdge = n->getLastOnTailEdge();
			if(!njh::in(n->uid_, uidsModified) &&
					!njh::in(firstTailEdge->tail_.lock()->uid_, uidsModified) &&
					!njh::in(lastTailEdge->tail_.lock()->uid_, uidsModified) &&
					!njh::in(firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_, uidsModified) &&
					!njh::in(firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->uid_, uidsModified) &&
					!njh::in(lastTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_, uidsModified) &&
					!njh::in(lastTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->uid_, uidsModified) &&
					!njh::in(firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_, uidsModified) &&
					!njh::in(firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_, uidsModified)  ){

				std::shared_ptr<node> node1;
				std::shared_ptr<node> node2;
				if(firstTailEdge->tail_.lock()->cnt_ == lastTailEdge->tail_.lock()->cnt_){
					if(firstTailEdge->tail_.lock()->k_.size() < lastTailEdge->tail_.lock()->k_.size()){
						node1 = firstTailEdge->tail_.lock();
						node2 = lastTailEdge->tail_.lock();
					}else{
						node2 = firstTailEdge->tail_.lock();
						node1 = lastTailEdge->tail_.lock();
					}
				}else{
					if(firstTailEdge->tail_.lock()->cnt_ > lastTailEdge->tail_.lock()->cnt_){
						node1 = firstTailEdge->tail_.lock();
						node2 = lastTailEdge->tail_.lock();
					}else{
						node2 = firstTailEdge->tail_.lock();
						node1 = lastTailEdge->tail_.lock();
					}
				}

				double approxPerBaseCoverageNode1 = static_cast<double>(node1->cnt_)/(node1->k_.length() - klen_ + 1);
				double approxPerBaseCoverageNode2 = static_cast<double>(node2->cnt_)/(node2->k_.length() - klen_ + 1);
				if(approxPerBaseCoverageNode1/approxPerBaseCoverageNode2 >= homopolymerIndelCollapseFreqMultiplier_){
					if(uAbsdiff(firstTailEdge->tail_.lock()->k_.length(), lastTailEdge->tail_.lock()->k_.length()) <= 10){
	//						std::cout << "firstTailEdge->tail_.lock()->k_.length(): " << firstTailEdge->tail_.lock()->k_.length() << std::endl;
	//						std::cout << "lastTailEdge->tail_.lock()->k_.length(): " << lastTailEdge->tail_.lock()->k_.length() << std::endl;
	//						std::cout << "std::max(firstTailEdge->tail_.lock()->k_.length(), lastTailEdge->tail_.lock()->k_.length()): " << std::max(firstTailEdge->tail_.lock()->k_.length(), lastTailEdge->tail_.lock()->k_.length()) << std::endl;
						uint64_t maxlen = std::max(firstTailEdge->tail_.lock()->k_.length(), lastTailEdge->tail_.lock()->k_.length());
						//gapScoringParameters gapPars(5,1,0,0,0,0);
						gapScoringParameters gapPars(5,1,5,1,5,1);
						aligner alignerObj(maxlen, gapPars);
						alignerObj.countEndGaps_ = true;
						alignerObj.weighHomopolymers_ = true;
						alignerObj.alignRegGlobal(
								seqInfo("1", node1->k_),
								seqInfo("2", node2->k_));
						alignerObj.profileAlignment(
								seqInfo("1", node1->k_),
								seqInfo("2", node2->k_), false, false, false);
	//					alignerObj.alignObjectA_.seqBase_.outPutSeqAnsi(std::cout);
	//					alignerObj.alignObjectB_.seqBase_.outPutSeqAnsi(std::cout);
						if(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)){
							std::vector<std::string> tailsFront;
							tailsFront.emplace_back(firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
							tailsFront.emplace_back(firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
							std::unordered_map<std::string, std::vector<std::shared_ptr<node>>> groupedNodes;
							groupedNodes[firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_].emplace_back(firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock());
							groupedNodes[firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_].emplace_back(firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock());
							groupedNodes[lastTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_].emplace_back(lastTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock());
							groupedNodes[lastTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_].emplace_back(lastTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock());
							if(2 == groupedNodes.size()){
								if(2 == groupedNodes[tailsFront.front()].size() &&
									 2 == groupedNodes[tailsFront.back()].size()){
									if( (groupedNodes[tailsFront.front()].front()->k_.size() == groupedNodes[tailsFront.front()].back()->k_.size()) &&
											(groupedNodes[tailsFront.back() ].front()->k_.size() == groupedNodes[tailsFront.back() ].back()->k_.size()) ){
										bool beginPass = false;
										bool endPass = false;
										{
											njh::sort(groupedNodes[tailsFront.front()], [](const std::shared_ptr<node> & node1, const std::shared_ptr<node> & node2){
												if(node1->cnt_ == node2->cnt_){
													return node1->k_.length() < node2->k_.length();
												}else{
													return node1->cnt_ > node2->cnt_;
												}
											});
											double approxPerBaseCoverageNode1 = static_cast<double>(groupedNodes[tailsFront.front()].front()->cnt_)/(groupedNodes[tailsFront.front()].front()->k_.length() - klen_ + 1);
											double approxPerBaseCoverageNode2 = static_cast<double>(groupedNodes[tailsFront.front()].back()->cnt_)/(groupedNodes[tailsFront.front()].back()->k_.length() - klen_ + 1);
											if(approxPerBaseCoverageNode1/approxPerBaseCoverageNode2 >= homopolymerIndelCollapseFreqMultiplier_){
												aligner alignerObj(std::max(
														groupedNodes[tailsFront.front()].front()->k_.length(),
														groupedNodes[tailsFront.front()].back()->k_.length()), gapScoringParameters(3,1,3,1,3,1));
												alignerObj.countEndGaps_ = false;
												alignerObj.weighHomopolymers_ = true;
												alignerObj.alignRegGlobal(
														seqInfo("1", groupedNodes[tailsFront.front()].front()->k_),
														seqInfo("2", groupedNodes[tailsFront.front()].back()->k_));
												alignerObj.profileAlignment(
														seqInfo("1", groupedNodes[tailsFront.front()].front()->k_),
														seqInfo("2", groupedNodes[tailsFront.front()].back()->k_), false, false, false);
												if(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)){
													beginPass = true;
												}
											}
										}
										{
											njh::sort(groupedNodes[tailsFront.back()], [](const std::shared_ptr<node> & node1, const std::shared_ptr<node> & node2){
												if(node1->cnt_ == node2->cnt_){
													return node1->k_.length() < node2->k_.length();
												}else{
													return node1->cnt_ > node2->cnt_;
												}
											});
											double approxPerBaseCoverageNode1 = static_cast<double>(groupedNodes[tailsFront.back()].front()->cnt_)/(groupedNodes[tailsFront.back()].front()->k_.length() - klen_ + 1);
											double approxPerBaseCoverageNode2 = static_cast<double>(groupedNodes[tailsFront.back()].back()->cnt_)/(groupedNodes[tailsFront.back()].back()->k_.length() - klen_ + 1);
											if(approxPerBaseCoverageNode1/approxPerBaseCoverageNode2 >= homopolymerIndelCollapseFreqMultiplier_){
												aligner alignerObj(std::max(
														groupedNodes[tailsFront.back()].front()->k_.length(),
														groupedNodes[tailsFront.back()].back()->k_.length()), gapScoringParameters(3,1,3,1,3,1));
												alignerObj.countEndGaps_ = false;
												alignerObj.weighHomopolymers_ = true;
												alignerObj.alignRegGlobal(
														seqInfo("1", groupedNodes[tailsFront.back()].front()->k_),
														seqInfo("2", groupedNodes[tailsFront.back()].back()->k_));
												alignerObj.profileAlignment(
														seqInfo("1", groupedNodes[tailsFront.back()].front()->k_),
														seqInfo("2", groupedNodes[tailsFront.back()].back()->k_), false, false, false);
												if(allowableErrorForHPIndexCollapse_.passErrorProfile(alignerObj.comp_)){
													endPass = true;
												}
											}
										}
										if(beginPass && endPass){
											uidsModified.emplace(n->uid_);
											uidsModified.emplace(firstTailEdge->tail_.lock()->uid_);
											uidsModified.emplace(lastTailEdge->tail_.lock()->uid_);
											uidsModified.emplace(firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
											uidsModified.emplace(firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->uid_);
											uidsModified.emplace(lastTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
											uidsModified.emplace(lastTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->uid_);
											uidsModified.emplace(firstTailEdge->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
											uidsModified.emplace(firstTailEdge->tail_.lock()->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_);
											modifiedNodes = true;
											//add in the other node's read names
											node1->inReadNamesIdx_.insert(
													node2->inReadNamesIdx_.begin(),
													node2->inReadNamesIdx_.end());
											// add in head read edge names too
											node1->getFirstOnHeadEdge()->inReadNamesIdx_.insert(
												node2->getFirstOnHeadEdge()->inReadNamesIdx_.begin(),
												node2->getFirstOnHeadEdge()->inReadNamesIdx_.end()
											);
											std::unordered_map<std::string, std::shared_ptr<node>> node1Tails;
											node1Tails[node1->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_] = node1->getFirstOnTailEdge()->tail_.lock();
											node1Tails[node1->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_] = node1->getLastOnTailEdge()->tail_.lock();
											{
												//add first
												//add node names
												node1Tails[node2->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_]->inReadNamesIdx_.insert(node2->getFirstOnTailEdge()->tail_.lock()->inReadNamesIdx_.begin(), node2->getFirstOnTailEdge()->tail_.lock()->inReadNamesIdx_.end());
												//add head edge names
												node1Tails[node2->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_]->getFirstOnHeadEdge()->inReadNamesIdx_.insert(node2->getFirstOnTailEdge()->tail_.lock()->getFirstOnHeadEdge()->inReadNamesIdx_.begin(), node2->getFirstOnTailEdge()->tail_.lock()->getFirstOnHeadEdge()->inReadNamesIdx_.end());
												//add tail edge names
												node1Tails[node2->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_]->getFirstOnTailEdge()->inReadNamesIdx_.insert(node2->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->inReadNamesIdx_.begin(), node2->getFirstOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->inReadNamesIdx_.end());
												//turn of edge and head and self
											}
											{
												//add second
												//add node names
												node1Tails[node2->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_]->inReadNamesIdx_.insert(node2->getLastOnTailEdge()->tail_.lock()->inReadNamesIdx_.begin(), node2->getLastOnTailEdge()->tail_.lock()->inReadNamesIdx_.end());
												//add head edge names
												node1Tails[node2->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_]->getFirstOnHeadEdge()->inReadNamesIdx_.insert(node2->getLastOnTailEdge()->tail_.lock()->getFirstOnHeadEdge()->inReadNamesIdx_.begin(), node2->getLastOnTailEdge()->tail_.lock()->getFirstOnHeadEdge()->inReadNamesIdx_.end());
												//add tail edge names
												node1Tails[node2->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->tail_.lock()->uid_]->getFirstOnTailEdge()->inReadNamesIdx_.insert(node2->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->inReadNamesIdx_.begin(), node2->getLastOnTailEdge()->tail_.lock()->getFirstOnTailEdge()->inReadNamesIdx_.end());
												//turn of edge and head and self
											}
											// consider adding the count of the other node
											// now turn off the node and it's edges
											auto firstOnTailEdge = node2->getFirstOnTailEdge();
											auto lastOnTailEdge = node2->getLastOnTailEdge();
											//nodes' tails
											firstOnTailEdge->tail_.lock()->getFirstOnTailEdge()->on_ = false;
											lastOnTailEdge->tail_.lock()->getFirstOnTailEdge()->on_ = false;
											//nodes
											firstOnTailEdge->tail_.lock()->on_ = false;
											lastOnTailEdge->tail_.lock()->on_ = false;
											//edge itself
											firstOnTailEdge->on_ = false;
											lastOnTailEdge->on_  = false;
											//turn node's head and itself off
											node2->getFirstOnHeadEdge()->on_ = false;
											node2->on_ = false;
										}
									}
								}
							}
						}
					}
				} //passed covered if statement
			}
		}
	}
	if(modifiedNodes){
		//remove off nodes;
		removeOffNodes();
		removeOffEdges();
		resetNodePositions();
		return true;
	}
	return false;
}




}// namespace njhseq

