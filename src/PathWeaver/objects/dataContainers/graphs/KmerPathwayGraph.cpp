/*
 * KmerPathwayGraph.cpp
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

#include "KmerPathwayGraph.hpp"
#include "PathWeaver/PathFinding/CoverageEstimator.hpp"
namespace njhseq {
KmerPathwayGraph::KmerPathwayGraph(uint32_t klen) :
		klen_(klen), occurenceCutOff_(1) {
	allowableErrorForHPIndexCollapse_.oneBaseIndel_ = .20;
}
KmerPathwayGraph::KmerPathwayGraph(uint32_t klen,uint32_t occurenceCutOff) :
		klen_(klen), occurenceCutOff_(occurenceCutOff) {
	allowableErrorForHPIndexCollapse_.oneBaseIndel_ = .20;
}
void KmerPathwayGraph::setOccurenceCutOff(uint32_t cutOff){
	occurenceCutOff_ = cutOff;
}



bool KmerPathwayGraph::hasNode(const std::string & nodeName) const {
	for(const auto & node : nodes_){
		if(nodeName == node->k_){
			return true;
		}
	}
	return false;
}
void KmerPathwayGraph::populateNodesFromCounts() {
	nodes_.clear();
	nodePositions_.clear();
	for (const auto & kCount : kCounts_) {
		if (kCount.second > occurenceCutOff_) {
			//addNode(kCount.first, kCount.second, kCount.first.size());
			addNode(kCount.first, kCount.second, klen_);
		}
	}
	resetNodePositions();
}


void KmerPathwayGraph::addNode(const std::string & k, uint32_t cnt, uint32_t kLen){
	nodes_.emplace_back(std::make_shared<node>(k, cnt, kLen));
}
void KmerPathwayGraph::addEdge(const std::string & head,
		const std::string & tailName,
		double cnt,
		const uint32_t & readName){
	//get head node
	//see if tail node already exists and if it does, just add count to the edge
	auto headNode = nodes_[nodePositions_.at(head)];
	bool foundEdge = false;
	for(const auto & tail : headNode->tailEdges_){
		auto tailNode = tail->tail_.lock();
		if(tailName == tailNode->k_){
			foundEdge = true;
			tail->cnt_ += cnt;
			tail->inReadNamesIdx_.emplace(readName);
			break;
		}
	}
	if(!foundEdge){
		auto tailNode = nodes_[nodePositions_.at(tailName)];
		std::shared_ptr<edge> e = std::make_shared<edge>(headNode, tailNode, cnt, std::unordered_set<uint32_t>{readName});
		headNode->addTail(e);
		tailNode->addHead(e);
	}
}


bool KmerPathwayGraph::hasCycle() const {
	/*
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	uint32_t headlessCount = 0;
	for (const auto & n : nodes_) {
		if (n->headless()) {
			++headlessCount;
			if(n->tailCount() > 0){
				nodesToProcess.push_back(n);
			}
		}
	}
	if(0 == headlessCount && nodes_.size() > 0){
		return true;
	}
	auto travelLookingForCycle = [](const node & n, VecStr uids){
	};
	for(const auto & n : nodesToProcess){
		for(const auto & t : n->tailEdges_){
			std::vector<std::string> uids;
			uids.emplace_back(t->tail_.lock()->uid_);
			auto next = n;
		}
	}*/
	return false;
}
bool KmerPathwayGraph::breakSingleHeadSingleTailNodesLowCoverage(){
	resetNodePositions();
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	for (const auto & n : nodes_) {
		if (1 == n->tailCount() && 1 == n->headCount()) {
			nodesToProcess.push_back(n);
		}
	}
	if(nodesToProcess.empty()){
		return false;
	}
	for(auto & n : nodesToProcess){
		//check again in case a previously edited node now makes thi untrue
		if (1 == n->tailCount() && 1 == n->headCount()) {
			auto firstHeadEdge = n->getFirstOnHeadEdge();
			auto firstTailEdge = n->getFirstOnTailEdge();
			//next should only have a head count of 1 for it to be here
			std::set<uint32_t> readsEntering(firstHeadEdge->inReadNamesIdx_.begin(), firstHeadEdge->inReadNamesIdx_.end());
			//and to reach here tailEdges_ should only be one as well
			std::set<uint32_t> readsLeaving (firstTailEdge->inReadNamesIdx_.begin(), firstTailEdge->inReadNamesIdx_.end());
			std::vector<uint32_t> readsThatEnteredLeaving;
			std::set_intersection(
					readsEntering.begin(), readsEntering.end(),
					readsLeaving.begin(), readsLeaving.end(),
					std::back_inserter(readsThatEnteredLeaving));
			if(readsThatEnteredLeaving.size() != std::min(firstHeadEdge->inReadNamesIdx_.size(), firstTailEdge->inReadNamesIdx_.size()) &&
//					readsThatEnteredLeaving.size() < 2){
				readsThatEnteredLeaving.size() <= occurenceCutOff_){
#if defined(PATHWEAVERSUPERDEBUG)
				{
					std::cout << "For next uid: " << n->uid_ << std::endl;
					std::cout << "\tReadsEnteringNumber    : " << firstHeadEdge->inReadNamesIdx_.size() << std::endl;
					std::cout << "\tReadsLeavingNumber     : " << firstTailEdge->inReadNamesIdx_.size() << std::endl;
					std::cout << "\tReadsLeavingThatEntered: " << readsThatEnteredLeaving.size() << std::endl;
					std::cout << std::endl;
				}
#endif
				//create new node, this will be the node with the tails, old node will have the heads
				std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
											KmerPathwayGraph::node>(n->k_, n->cnt_, n->kLen_);
				newNode->inReadNamesIdx_ = n->inReadNamesIdx_;
				//reset tail edges
				newNode->tailEdges_ = n->tailEdges_;
				n->tailEdges_.clear();
				for(auto & tailEdge : newNode->tailEdges_){
					tailEdge->head_ = newNode;
				}
				nodes_.emplace_back(newNode);
			}
		}
	}
	resetNodePositions();
	return true;
}
void KmerPathwayGraph::breakSingleLinkedPathsReadThreading(){
	resetNodePositions();
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	for (const auto & n : nodes_) {
		if (1 != n->headCount() && 1 == n->tailCount()) {
			nodesToProcess.push_back(n);
		} else if (n->tailCount() > 1) {
			for (const auto & tail : n->tailEdges_) {
				if(!tail->on_){
					continue;
				}
				//the head check is to avoid to adding duplicates from the addition above that addes
				//nodes with head counts not equal to 1
				if (1 == tail->tail_.lock()->tailCount() && 1 == tail->tail_.lock()->headCount()) {
					nodesToProcess.push_back(tail->tail_.lock());
				}
			}
		}
	}
//	for(auto & n : nodesToProcess){
//		std::cout << "Need to process " << n->k_ << ", tailCount: " << n->tailCount() << ", headCount: " << n->headCount() << std::endl;
//	}
	for(auto & n : nodesToProcess){
		std::vector<std::shared_ptr<njhseq::KmerPathwayGraph::node>> nexts;
//		std::cout << "On " << n->k_ << std::endl;
		auto next = n->getFirstOnTailEdge()->tail_.lock();
		//a tail count of one indicates that the next node should be added in
		//also head count should be less than 2
		while(next != nullptr && next->headCount() == 1){
			nexts.emplace_back(next);
			if (0 == next->tailCount()) {
				//head a tailless node, end of the line
				next = nullptr;
			} else if (next->tailCount() > 1) {
				//head a multi tailed node, section collapse is done
				next = nullptr;
			} else {
				//can still collapse more set next to the tail edge
				//will only hit here if tailEdges_.size() is 1
				next = next->getFirstOnTailEdge()->tail_.lock();
			}
		}
		std::unordered_map<std::string, std::vector<uint32_t>> readsNamesInNodePositions;
		uint32_t nodePosition = 0;
		for(const auto & readNameIdx : n->inReadNamesIdx_){
			readsNamesInNodePositions[readNames_[readNameIdx]].emplace_back(nodePosition);
		}
//		if("ACCCCAATGCAAATCCTAATGCAAATCCTAATGCC" == n->k_){
//			std::cout << njh::bashCT::bold << njh::bashCT::green << std::endl;
//		}
//		std::cout << n->k_ << std::endl;
		for(const auto & nextNode : nexts){
			++nodePosition;
			for(const auto & readNameIdx : nextNode->inReadNamesIdx_){
				readsNamesInNodePositions[readNames_[readNameIdx]].emplace_back(nodePosition);
			}
//			std::cout << "\t" << nextNode->k_ << std::endl;
			//std::cout << "\t\t" << njh::conToStr(nextNode->inReadNamesIdx_, ",") << std::endl;
		}
//		if("ACCCCAATGCAAATCCTAATGCAAATCCTAATGCC" == n->k_){
//			std::cout << njh::bashCT::reset << std::endl;
//		}
//		OutOptions outOptsColor(bfs::path(n->uid_), ".tab.txt");
//		outOptsColor.overWriteFile_ = true;
//		OutputStream outColors(outOptsColor);
		std::vector<uint32_t> allNodePositions(nexts.size() + 1);
		njh::iota<uint32_t>(allNodePositions, 0);
//		outColors << "readName\t" << njh::conToStr(allNodePositions, "\t") << std::endl;
		std::vector<std::string> readNames = getVectorOfMapKeys(readsNamesInNodePositions);
		njh::sort(readNames, [&readsNamesInNodePositions](const std::string & r1, const std::string & r2){
			if(readsNamesInNodePositions[r1].front() == readsNamesInNodePositions[r2].front()){
				return readsNamesInNodePositions[r1].back() < readsNamesInNodePositions[r2].back() ;
			}else{
				return readsNamesInNodePositions[r1].front() < readsNamesInNodePositions[r2].front() ;
			}
		});
		struct NodeStartStop {
			NodeStartStop(uint32_t firstNode, uint32_t lastNode) :
					firstNode_(firstNode), lastNode_(lastNode) {
			}
			uint32_t firstNode_;
			uint32_t lastNode_;
		};
		std::vector<NodeStartStop> paths;
		for(const auto & readName: readNames){
			njh::sort(readsNamesInNodePositions[readName]);
			const auto & readPositions = readsNamesInNodePositions[readName];
			uint32_t lastNode =  readsNamesInNodePositions[readName].front();
			uint32_t startNode = readsNamesInNodePositions[readName].front();
			if(1 == readsNamesInNodePositions[readName].size()){
				paths.emplace_back(NodeStartStop{startNode, lastNode});
			} else {
				if(readsNamesInNodePositions[readName].size() > 2){
					for(const auto pos : iter::range<uint32_t>(1, readsNamesInNodePositions[readName].size() - 1)){
						if(lastNode +1 != readsNamesInNodePositions[readName][pos]){
							paths.emplace_back(NodeStartStop{startNode, lastNode});
							startNode = readsNamesInNodePositions[readName][pos];
						}
						lastNode = readsNamesInNodePositions[readName][pos];
					}
				}
				if (lastNode + 1 != readsNamesInNodePositions[readName].back()) {
					paths.emplace_back(NodeStartStop { startNode, lastNode });
					startNode = readsNamesInNodePositions[readName].back();
				}
				lastNode = readsNamesInNodePositions[readName].back();
				paths.emplace_back(NodeStartStop{startNode, lastNode});
			}
//			outColors << readName << "\t";
			std::vector<uint32_t> allNodePositionsForRead(nexts.size() + 1, 0);
			for(const auto & pos : readPositions){
				allNodePositionsForRead[pos] = 1;
			}
//			outColors << njh::conToStr(allNodePositionsForRead, "\t") << std::endl;
		}
//		OutOptions outOptsPathsLength(bfs::path(n->uid_ + "_pathLengths"), ".tab.txt");
//		outOptsPathsLength.overWriteFile_ = true;
//		OutputStream outPathsLength(outOptsPathsLength);
//		outPathsLength << "readName\tmeanLen\tmedianLen" << std::endl;
//
//		OutOptions outOptsPaths(bfs::path(n->uid_ + "_consectivePaths"), ".tab.txt");
//		outOptsPaths.overWriteFile_ = true;
//		OutputStream outPaths(outOptsPaths);
//		outPaths << "readName\t" << njh::conToStr(allNodePositions, "\t") << std::endl;
		njh::sort(paths, [](const NodeStartStop & obj1, const NodeStartStop & obj2){
			if(obj1.firstNode_ == obj2.firstNode_){
				return obj1.lastNode_ < obj2.lastNode_;
			}else{
				return obj1.firstNode_ < obj2.firstNode_;
			}
		});
		uint32_t pathNumber = 0;
		std::vector<uint32_t> pathLength;
		for(const auto & path : paths){
			pathLength.emplace_back(path.lastNode_ + 1 - path.firstNode_);
//			outPaths << pathNumber << "\t";
			std::vector<uint32_t> allNodePositionsForRead(nexts.size() + 1, 0);
			for(const auto & pos : iter::range(path.firstNode_, path.lastNode_ + 1)){
				allNodePositionsForRead[pos] = 1;
			}
//			outPaths << njh::conToStr(allNodePositionsForRead, "\t") << std::endl;
			++pathNumber;
		}
		auto meanPathLen = vectorMean(pathLength);
//		auto medPathLen = vectorMedianRef(pathLength);
//		outPathsLength << pathNumber
//				<< "\t" << meanPathLen
//				<< "\t" << medPathLen << std::endl;
//
//		OutOptions outOptsPathsReach(bfs::path(n->uid_ + "_pathReach"), ".tab.txt");
//		outOptsPathsReach.overWriteFile_ = true;
//		OutputStream outPathsReach(outOptsPathsReach);
//		outPathsReach << "Position\totherPosition\tcount" << std::endl;
		struct PosReachInfo {
			PosReachInfo(uint32_t pos, uint32_t otherPos, uint32_t readCount) :
					pos_(pos), otherPos_(otherPos), readCount_(readCount) {
			}
			uint32_t pos_;
			uint32_t otherPos_;
			uint32_t readCount_;
		};
		uint32_t reachCutOff = 1;
		std::vector<PosReachInfo> reachInfos;
		for(const auto pos : iter::range(nexts.size() + 1 ) ){
			//uint32_t otherPostion = round(pos + std::min<double>(meanPathLen, 50));
			uint32_t otherPostion = round(pos + meanPathLen);
			if(otherPostion >= nexts.size() + 1){
				break;
			}
			uint32_t pathNumbers = 0;
			for(const auto & path : paths){
				if(path.firstNode_ > pos){
					continue;
				}
				if(path.lastNode_ >= otherPostion){
					++pathNumbers;
				}
			}
			reachInfos.emplace_back(pos, otherPostion, pathNumbers);
//			outPathsReach << pos
//					<< "\t" << otherPostion
//					<< "\t" << pathNumbers << std::endl;
		}
		if(reachInfos.size() > 2){
			//pair 1) start 2) end (inclusive)
			std::vector<std::pair<uint32_t, uint32_t>> lowCutStreaks;
			uint32_t start = std::numeric_limits<uint32_t>::max();
			uint32_t len = 0;
			for(const auto pos : iter::range(reachInfos.size())) {
				if(reachInfos[pos].readCount_ <= reachCutOff){
					if(std::numeric_limits<uint32_t>::max() == start){
						start = pos;
					}else{
						++len;
					}
				}else{
					if(len > 0){
						lowCutStreaks.emplace_back(std::make_pair(start, start + len) );
					}
					start = std::numeric_limits<uint32_t>::max();
					len = 0;
				}
			}
			if(len > 0){
				lowCutStreaks.emplace_back(std::make_pair(start, start + len) );
			}
//			OutOptions outOptsPathsLowReach(bfs::path(n->uid_ + "_lowReachStreaks"), ".tab.txt");
//			outOptsPathsLowReach.overWriteFile_ = true;
//			OutputStream outPathsLowReach(outOptsPathsLowReach);
//			outPathsLowReach << "streak\tpositions\tfromBeg\ttoEnd" << std::endl;
			uint32_t streakNumber = 0;
			std::vector<std::pair<uint32_t, uint32_t>> lowCutStreaksKept;
			for(const auto & streak : lowCutStreaks){
//				bool containsNonKlenNode = false;
//				for(const auto & pos : iter::range(reachInfos[streak.first].pos_ + 1, reachInfos[streak.second].otherPos_)){
//					if(0 == pos){
//						if(n->k_.size() != klen_){
//							containsNonKlenNode = true;
//							break;
//						}
//					}else{
//						if(nexts[pos - 1]->k_.size() != klen_){
//							containsNonKlenNode = true;
//							break;
//						}
//					}
//				}
//				if(containsNonKlenNode){
//					continue;
//				}
				bool fromBeg = false;
				if(streak.first == 0){
					fromBeg = true;
				}
				bool toEnd = false;
				if(reachInfos.size() ==  streak.second + 1){
					toEnd = true;
				}
				if(!toEnd && !fromBeg){
					//only count streaks that occur in areas that aren't in high coverage
					if(reachInfos[streak.first - 1].readCount_ > 5 && reachInfos[streak.second + 1].readCount_ > 5){
//						outPathsLowReach << streakNumber
//								<< "\t" << streak.first << ":" << streak.second
//								<< "\t" << estd::to_string(fromBeg)
//								<< "\t" << estd::to_string(toEnd) << std::endl;
						lowCutStreaksKept.emplace_back(streak);
						++streakNumber;
					}
				}else{
//					outPathsLowReach << streakNumber
//							<< "\t" << streak.first << ":" << streak.second
//							<< "\t" << estd::to_string(fromBeg)
//							<< "\t" << estd::to_string(toEnd) << std::endl;
					lowCutStreaksKept.emplace_back(streak);
					++streakNumber;
				}
			}
			for(const auto & streak : lowCutStreaksKept){
				bool fromBeg = false;
				if(streak.first == 0){
					fromBeg = true;
				}
				bool toEnd = false;
				if(reachInfos.size() ==  streak.second + 1){
					toEnd = true;
				}
				if(fromBeg){
					uint32_t nextNodePosition = streak.second;
					{
						//add begining, starts with n->
						n->tailEdges_.clear();
						std::vector<std::shared_ptr<KmerPathwayGraph::node>> newNodes;
						std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
								KmerPathwayGraph::node>(
										nexts[0]->k_,
										nexts[0]->cnt_,
										nexts[0]->kLen_);
						for(const auto & name : nexts[0]->inReadNamesIdx_){
							if(njh::in(name, n->inReadNamesIdx_)){
								newNode->inReadNamesIdx_.emplace(name);
							}
						}
						std::shared_ptr<edge> e = std::make_shared<edge>(
								n,
								newNode,
								newNode->inReadNamesIdx_.size(),
								newNode->inReadNamesIdx_);
						n->addTail(e);
						newNode->addHead(e);
						newNodes.emplace_back(newNode);
						for(uint32_t addingNode = 1; addingNode < reachInfos[streak.second].otherPos_; ++ addingNode){
							std::unordered_set<uint32_t> readNames;
							for(const auto & name : nexts[addingNode]->inReadNamesIdx_){
								if(njh::in(name, newNodes.back()->inReadNamesIdx_)){
									readNames.emplace(name);
								}
							}
							if(readNames.size() > reachCutOff){
								std::shared_ptr<KmerPathwayGraph::node> currentNewNode = std::make_shared<
										KmerPathwayGraph::node>(
												nexts[addingNode]->k_,
												nexts[addingNode]->cnt_,
												nexts[addingNode]->kLen_);
								currentNewNode->inReadNamesIdx_ = readNames;
								std::shared_ptr<edge> e = std::make_shared<edge>(
										newNodes.back(),
										currentNewNode,
										currentNewNode->inReadNamesIdx_.size(),
										currentNewNode->inReadNamesIdx_);
								newNodes.back()->addTail(e);
								currentNewNode->addHead(e);
								newNodes.emplace_back(currentNewNode);
							}else{
								break;
							}
						}
						for(const auto & addNode : newNodes){
							nodes_.emplace_back(addNode);
						}
						//turn off beginning nodes
						//n->on_ = false;
						for(const auto pos : iter::range(streak.second)){
							nexts[pos]->on_ = false;
						}
					}
					{
						//break the connection to the rest of the path
						nexts[nextNodePosition]->headEdges_.clear();
						nexts[nextNodePosition - 1]->tailEdges_.clear();
					}
					//right here
				}else if(toEnd){
					uint32_t nextNodePos = reachInfos[streak.first].pos_ - 1;
					{
						std::shared_ptr<KmerPathwayGraph::node> previousNode;
						if(0 == nextNodePos){
							previousNode = n;
						}else{
							previousNode = nexts[nextNodePos - 1];
						}
						//add the connected end portions
						//clear the tail edge of the previous node
						previousNode->tailEdges_.clear();
						nexts[nextNodePos]->headEdges_.clear();
						std::vector<std::shared_ptr<KmerPathwayGraph::node>> newNodes;
						std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
								KmerPathwayGraph::node>(
										nexts[nextNodePos]->k_,
										nexts[nextNodePos]->cnt_,
										nexts[nextNodePos]->kLen_);
						for(const auto & name : nexts[nextNodePos]->inReadNamesIdx_){
							if(njh::in(name, previousNode->inReadNamesIdx_)){
								newNode->inReadNamesIdx_.emplace(name);
							}
						}
						std::shared_ptr<edge> e = std::make_shared<edge>(
								previousNode,
								newNode,
								newNode->inReadNamesIdx_.size(),
								newNode->inReadNamesIdx_);
						previousNode->addTail(e);
						newNode->addHead(e);
						newNodes.emplace_back(newNode);
						for(uint32_t addingNodePos = nextNodePos + 1; addingNodePos < nexts.size() - 1; ++addingNodePos){
							std::unordered_set<uint32_t> readNames;
							for(const auto & readName : nexts[addingNodePos]->inReadNamesIdx_){
								if(njh::in(readName, newNodes.back()->inReadNamesIdx_)){
									readNames.emplace(readName);
								}
							}
							if(readNames.size() > reachCutOff){
								std::shared_ptr<KmerPathwayGraph::node> nextNewNode = std::make_shared<
										KmerPathwayGraph::node>(nexts[addingNodePos]->k_,
												nexts[addingNodePos]->cnt_,
												nexts[addingNodePos]->kLen_);
								nextNewNode->inReadNamesIdx_ = readNames;
								std::shared_ptr<edge> e = std::make_shared<edge>(
										newNodes.back(),
										nextNewNode,
										nextNewNode->inReadNamesIdx_.size(),
										nextNewNode->inReadNamesIdx_);
								newNodes.back()->addTail(e);
								nextNewNode->addHead(e);
								newNodes.emplace_back(nextNewNode);
							}else{
								break;
							}
						}
						for(const auto & addNode : newNodes){
							nodes_.emplace_back(addNode);
						}
					}
					{
						//complete the end
						uint32_t endNodePos = nexts.size() - 1;
						nexts[endNodePos]->headEdges_.clear();
						nexts[endNodePos - 1]->tailEdges_.clear();
						std::vector<std::shared_ptr<KmerPathwayGraph::node>> newNodes;
						std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
								KmerPathwayGraph::node>(
										nexts[endNodePos - 1]->k_,
										nexts[endNodePos - 1]->cnt_,
										nexts[endNodePos - 1]->kLen_);
						for(const auto & name : nexts[endNodePos]->inReadNamesIdx_){
							if(njh::in(name, nexts[endNodePos - 1]->inReadNamesIdx_)){
								newNode->inReadNamesIdx_.emplace(name);
							}
						}
						std::shared_ptr<edge> e = std::make_shared<edge>(
								newNode,
								nexts[endNodePos],
								newNode->inReadNamesIdx_.size(),
								newNode->inReadNamesIdx_);
						nexts[endNodePos]->addHead(e);
						newNode->addTail(e);
						newNodes.emplace_back(newNode);
						if(endNodePos > 1){
							for (uint32_t addingNodePos = endNodePos - 1 - 1; addingNodePos > nextNodePos; --addingNodePos) {
								std::unordered_set<uint32_t> readNames;
								for(const auto & readName : nexts[addingNodePos]->inReadNamesIdx_){
									if(njh::in(readName, newNodes.back()->inReadNamesIdx_)){
										readNames.emplace(readName);
									}
								}
								if(readNames.size() > reachCutOff){
									std::shared_ptr<KmerPathwayGraph::node> nextNewNode = std::make_shared<
											KmerPathwayGraph::node>(
													nexts[addingNodePos]->k_,
													nexts[addingNodePos]->cnt_,
													nexts[addingNodePos]->kLen_);
									nextNewNode->inReadNamesIdx_ = readNames;
									std::shared_ptr<edge> e = std::make_shared<edge>(
											nextNewNode,
											newNodes.back(),
											nextNewNode->inReadNamesIdx_.size(),
											nextNewNode->inReadNamesIdx_);
									newNodes.back()->addHead(e);
									nextNewNode->addTail(e);
									newNodes.emplace_back(nextNewNode);
								}else{
									break;
								}
							}
							for(const auto & addNode : newNodes){
								nodes_.emplace_back(addNode);
							}
						}
					}
					for(const auto & nodePos : iter::range<uint32_t>(nextNodePos, nexts.size() - 1)){
						nexts[nodePos]->on_ = false;
					}
					//
				}else{
					//first node
					uint32_t forwardGoingNodePos = streak.second;
					//last node
					uint32_t goingBackWardsNodePos = reachInfos[streak.first-1].otherPos_ - 1;
					//breakage
					if(0 == forwardGoingNodePos){
						n->tailEdges_.clear();
					}else{
						nexts[forwardGoingNodePos - 1]->tailEdges_.clear();
						nexts[forwardGoingNodePos]->headEdges_.clear();
					}
					nexts[goingBackWardsNodePos]->tailEdges_.clear();
					nexts[goingBackWardsNodePos + 1]->headEdges_.clear();
					for(const auto pos : iter::range(forwardGoingNodePos, goingBackWardsNodePos + 1)){
						nexts[pos]->on_ = false;
					}
					{
						//front portion
						std::vector<std::shared_ptr<KmerPathwayGraph::node>> newNodes;
						std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
								KmerPathwayGraph::node>(
										nexts[forwardGoingNodePos]->k_,
										nexts[forwardGoingNodePos]->cnt_,
										nexts[forwardGoingNodePos]->kLen_);
						std::shared_ptr<KmerPathwayGraph::node> previousNode;
						if(0 == forwardGoingNodePos){
							previousNode = n;
						}else{
							previousNode = nexts[forwardGoingNodePos - 1];
						}
						for(const auto & name : nexts[forwardGoingNodePos]->inReadNamesIdx_){
							if(njh::in(name, previousNode->inReadNamesIdx_)){
								newNode->inReadNamesIdx_.emplace(name);
							}
						}
						std::shared_ptr<edge> e = std::make_shared<edge>(
								previousNode,
								newNode,
								newNode->inReadNamesIdx_.size(),
								newNode->inReadNamesIdx_);
						previousNode->addTail(e);
						newNode->addHead(e);
						newNodes.emplace_back(newNode);
						for(uint32_t addingNodePos = forwardGoingNodePos + 1; forwardGoingNodePos < goingBackWardsNodePos + 1; ++ addingNodePos){
							std::unordered_set<uint32_t> readNames;
							for(const auto & name : nexts[addingNodePos]->inReadNamesIdx_){
								if(njh::in(name, newNodes.back()->inReadNamesIdx_)){
									readNames.emplace(name);
								}
							}
							if(readNames.size() > reachCutOff){
								std::shared_ptr<KmerPathwayGraph::node> currentNewNode = std::make_shared<
										KmerPathwayGraph::node>(
												nexts[addingNodePos]->k_,
												nexts[addingNodePos]->cnt_,
												nexts[addingNodePos]->kLen_);
								currentNewNode->inReadNamesIdx_ = readNames;
								std::shared_ptr<edge> e = std::make_shared<edge>(
										newNodes.back(),
										currentNewNode,
										currentNewNode->inReadNamesIdx_.size(),
										currentNewNode->inReadNamesIdx_);
								newNodes.back()->addTail(e);
								currentNewNode->addHead(e);
								newNodes.emplace_back(currentNewNode);
							}else{
								break;
							}
						}
						for(const auto & addNode : newNodes){
							nodes_.emplace_back(addNode);
						}
					} //end front portion
					{
						//back portion
						std::vector<std::shared_ptr<KmerPathwayGraph::node>> newNodes;
						std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
								KmerPathwayGraph::node>(
										nexts[goingBackWardsNodePos]->k_,
										nexts[goingBackWardsNodePos]->cnt_,
										nexts[goingBackWardsNodePos]->kLen_);
						for(const auto & name : nexts[forwardGoingNodePos]->inReadNamesIdx_){
							if(njh::in(name, nexts[goingBackWardsNodePos + 1]->inReadNamesIdx_)){
								newNode->inReadNamesIdx_.emplace(name);
							}
						}
						std::shared_ptr<edge> e = std::make_shared<edge>(
								newNode,
								nexts[goingBackWardsNodePos + 1],
								newNode->inReadNamesIdx_.size(),
								newNode->inReadNamesIdx_);
						nexts[goingBackWardsNodePos + 1]->addHead(e);
						newNode->addTail(e);
						newNodes.emplace_back(newNode);
						for(uint32_t addingNodePos = goingBackWardsNodePos - 1; addingNodePos >= forwardGoingNodePos; --addingNodePos){
							std::unordered_set<uint32_t> readNames;
							for(const auto & name : nexts[addingNodePos]->inReadNamesIdx_){
								if(njh::in(name, newNodes.back()->inReadNamesIdx_)){
									readNames.emplace(name);
								}
							}
							if(readNames.size() > reachCutOff){
								std::shared_ptr<KmerPathwayGraph::node> currentNewNode = std::make_shared<
										KmerPathwayGraph::node>(
												nexts[addingNodePos]->k_,
												nexts[addingNodePos]->cnt_,
												nexts[addingNodePos]->kLen_);
								currentNewNode->inReadNamesIdx_ = readNames;
								std::shared_ptr<edge> e = std::make_shared<edge>(
										currentNewNode,
										newNodes.back(),
										currentNewNode->inReadNamesIdx_.size(),
										currentNewNode->inReadNamesIdx_);
								newNodes.back()->addHead(e);
								currentNewNode->addTail(e);
								newNodes.emplace_back(currentNewNode);
							}else{
								break;
							}
						}
						for(const auto & addNode : newNodes){
							nodes_.emplace_back(addNode);
						}
					}// end of back portion
				} // end middle breakage
			}
		}
//		std::cout << std::endl;
	}
	removeOffNodes();
//removeNullNodes();
	resetNodePositions();
	//exit(1);
}


void KmerPathwayGraph::collapseSingleLinkedPathsForPossibleLoops(){


	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesWithOneTailedTailsToProcess;
	{
		//grab the all nodes with one tail connected to a node with just one tail
		std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesWithOneTailedTails;
		for(const auto & n : nodes_){
			if(1 == n->tailCount()
					&& 1 == n->getFirstOnTailEdge()->tail_.lock()->tailCount() && 1 == n->getFirstOnTailEdge()->tail_.lock()->headCount()
					&& n->uid_ != n->getFirstOnTailEdge()->tail_.lock()->uid_){
				nodesWithOneTailedTails.push_back(n);
			}
		}
		njh::sort(nodesWithOneTailedTails,[](const std::shared_ptr<KmerPathwayGraph::node> & n1, const std::shared_ptr<KmerPathwayGraph::node> & n2){
			uint32_t n1HeadCntCount = 0;
			uint32_t n2HeadCntCount = 0;
			for(const auto & h : n1->headEdges_){
				if(h->on_){
					n1HeadCntCount += h->cnt_;
				}
			}
			for(const auto & h : n2->headEdges_){
				if(h->on_){
					n2HeadCntCount += h->cnt_;
				}
			}
			if(n1HeadCntCount == n2HeadCntCount){
				return n1->k_.size() < n2->k_.size();
			}else{
				return n1HeadCntCount < n2HeadCntCount;
			}
		});
		if(!nodesWithOneTailedTails.empty()){
			nodesWithOneTailedTailsToProcess.push_back(nodesWithOneTailedTails.front());
		}
	}
	while(!nodesWithOneTailedTailsToProcess.empty()){
		std::vector<uint32_t> nodesToErase;
		for(auto & n : nodesWithOneTailedTailsToProcess){
			//auto firstTailEdge = n->getFirstOnTailEdge();
			//std::cout << "On " << n->k_ << std::endl;
			auto next = n->getFirstOnTailEdge()->tail_.lock();
			//a tail count of one indicates that the next node should be added in
			//also head count should be less than 2
			while(next != nullptr && next->headCount() == 1 && next->uid_ != n->uid_){
				////pushing backing only the back of the next k only works if there has been no processing before
				//std::cout << next->k_ << std::endl;
				//std::cout << next->k_.size() << std::endl;
				n->k_.append(next->k_.substr(klen_ - 1 ));
				n->cnt_+= next->cnt_;
				//n->inReadNamesIdx_.reserve(n->inReadNamesIdx_.size() + next->inReadNamesIdx_.size());
				n->inReadNamesIdx_.insert(next->inReadNamesIdx_.begin(), next->inReadNamesIdx_.end());
				auto toErase = next->uid_;
				if (0 == next->tailCount()) {
					//head a tailless node, end of the line
					//need to erase tail edges
					n->tailEdges_.erase(n->tailEdges_.begin());
					next = nullptr;
				} else if (next->tailCount() > 1) {
					//head a multi tailed node, section collapse is done
					//need to add all the multiple tail edges
					n->tailEdges_.clear();
					for (const auto & tail : next->tailEdges_) {
						if(!tail->on_){
							continue;
						}
						n->tailEdges_.push_back(tail);
						n->tailEdges_.back()->head_ = n;
					}
					next = nullptr;
				} else {
					//can still collapse more set next to the tail edge
					//will only hit here if tailEdges_.size() is 1
					n->tailEdges_.clear();
					for (const auto & tail : next->tailEdges_) {
						if (!tail->on_) {
							continue;
						}
						n->tailEdges_.push_back(tail);
						n->tailEdges_.back()->head_ = n;
					}
					next = next->getFirstOnTailEdge()->tail_.lock();
				}
				nodesToErase.emplace_back(nodePositions_[toErase]);
				nodes_[nodePositions_[toErase]] = nullptr;
			}
		}
		if(!nodesToErase.empty()){
			std::sort(nodesToErase.rbegin(), nodesToErase.rend());
			for(const auto & remove : nodesToErase){
				nodes_.erase(nodes_.begin() + remove);
			}
			resetNodePositions();
		}
		nodesWithOneTailedTailsToProcess.clear();
		{
			//grab the all nodes with one tail connected to a node with just one tail
			std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesWithOneTailedTails;
			for(const auto & n : nodes_){
				if(1 == n->tailCount()
						&& 1 == n->getFirstOnTailEdge()->tail_.lock()->tailCount() && 1 == n->getFirstOnTailEdge()->tail_.lock()->headCount()
						&& n->uid_ != n->getFirstOnTailEdge()->tail_.lock()->uid_){
					nodesWithOneTailedTails.push_back(n);
				}
			}
			njh::sort(nodesWithOneTailedTails,[](const std::shared_ptr<KmerPathwayGraph::node> & n1, const std::shared_ptr<KmerPathwayGraph::node> & n2){
				uint32_t n1HeadCntCount = 0;
				uint32_t n2HeadCntCount = 0;
				for(const auto & h : n1->headEdges_){
					if(h->on_){
						n1HeadCntCount += h->cnt_;
					}
				}
				for(const auto & h : n2->headEdges_){
					if(h->on_){
						n2HeadCntCount += h->cnt_;
					}
				}
				if(n1HeadCntCount == n2HeadCntCount){
					return n1->k_.size() < n2->k_.size();
				}else{
					return n1HeadCntCount < n2HeadCntCount;
				}
			});
			if(!nodesWithOneTailedTails.empty()){
				nodesWithOneTailedTailsToProcess.push_back(nodesWithOneTailedTails.front());
			}
		}
	}
}


void KmerPathwayGraph::collapseSingleLinkedPaths(bool initialCollapse){
	//njh::stopWatch watch;
	//watch.setLapName(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps(), 100) + " - preprocess - remove off nodes");
	removeOffNodes();
//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps(), 100) + " - preprocess - remove off edges");
	removeOffEdges();
////watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - preprocess - reseting node positions");
//no longer need this resting of node positions because it's in the remove off nodes
//resetNodePositions();

//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - finding nodes to process");
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	for (const auto & n : nodes_) {
		if (1 != n->headCount() && 1 == n->tailCount()) {
			nodesToProcess.push_back(n);
		} else if (n->tailCount() > 1) {
			for (const auto & tail : n->tailEdges_) {
				if(!tail->on_){
					continue;
				}
				//the head check is to avoid to adding duplicates from the addition above that addes
				//nodes with head counts not equal to 1
				if (1 == tail->tail_.lock()->tailCount() && 1 == tail->tail_.lock()->headCount()) {
					nodesToProcess.push_back(tail->tail_.lock());
				}
			}
		}
	}
#if defined(PATHWEAVERSUPERDEBUG)
	{
		for(auto & n : nodesToProcess){
			std::cout << "Need to process " << n->k_ << ", tailCount: " << n->tailCount() << ", headCount: " << n->headCount() << std::endl;
		}
	}

	if(!initialCollapse){
		for(const auto & n : nodesToProcess){
			if(n->headCount() == 1 && n->tailCount() == 1){
				//next should only have a head count of 1 for it to be here
				auto firstHeadEdge = n->getFirstOnHeadEdge();
				auto firstTailEdge = n->getFirstOnTailEdge();
				std::set<uint32_t> readsEntering(firstHeadEdge->inReadNamesIdx_.begin(), firstHeadEdge->inReadNamesIdx_.end());
				//and to reach here tailEdges_ should only be one as well
				std::set<uint32_t> readsLeaving (firstTailEdge->inReadNamesIdx_.begin(), firstTailEdge->inReadNamesIdx_.end());
				std::vector<uint32_t> readsThatEnteredLeaving;
				std::set_intersection(
						readsEntering.begin(), readsEntering.end(),
						readsLeaving.begin(), readsLeaving.end(),
						std::back_inserter(readsThatEnteredLeaving));
				if(readsThatEnteredLeaving.size() != std::min(firstHeadEdge->inReadNamesIdx_.size(), firstTailEdge->inReadNamesIdx_.size()) &&
						readsThatEnteredLeaving.size() <= occurenceCutOff_){
					std::cout << "For next uid: " << n->uid_ << std::endl;
					std::cout << "\tReadsEnteringNumber    : " << firstHeadEdge->inReadNamesIdx_.size() << std::endl;
					std::cout << "\tReadsLeavingNumber     : " << firstTailEdge->inReadNamesIdx_.size() << std::endl;
					std::cout << "\tReadsLeavingThatEntered: " << readsThatEnteredLeaving.size() << std::endl;
					std::cout << std::endl;
				}
			}
		}
	}
#endif



	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - processing nodes");
	if(initialCollapse && nodesToProcess.empty() && nodes_.size() > 1){
		collapseSingleLinkedPathsForPossibleLoops();
	} else {
		std::vector<uint32_t> nodesToErase;
		if (nodesToProcess.size() > 1 && numThreads_ > 1) {
			std::vector<uint32_t> nodesToProcessPositions(nodesToProcess.size());
			njh::iota<uint32_t>(nodesToProcessPositions, 0);
			njh::concurrent::LockableQueue<uint32_t> posQueue(nodesToProcessPositions);
			std::mutex nodesToEraseMut;
			std::function<void()> collapsePaths = [this,&nodesToErase,&nodesToEraseMut,&posQueue,&nodesToProcess](){
				uint32_t position = 0;
				std::vector<uint32_t> nodesToEraseThread;
				while(posQueue.getVal(position)){
					auto & n = nodesToProcess[position];
					//auto firstTailEdge = n->getFirstOnTailEdge();
					//std::cout << "On " << n->k_ << std::endl;
					auto next = n->getFirstOnTailEdge()->tail_.lock();
					//a tail count of one indicates that the next node should be added in
					//also head count should be less than 2
					while(next != nullptr && next->headCount() == 1 && next->uid_ != n->uid_){
						////pushing backing only the back of the next k only works if there has been no processing before
						//std::cout << next->k_ << std::endl;
						//std::cout << next->k_.size() << std::endl;
						n->k_.append(next->k_.substr(klen_ - 1 ));
						n->cnt_+= next->cnt_;
						//n->inReadNamesIdx_.reserve(n->inReadNamesIdx_.size() + next->inReadNamesIdx_.size());
						n->inReadNamesIdx_.insert(next->inReadNamesIdx_.begin(), next->inReadNamesIdx_.end());
						auto toErase = next->uid_;
						if (0 == next->tailCount()) {
							//head a tailless node, end of the line
							//need to erase tail edges
							n->tailEdges_.erase(n->tailEdges_.begin());
							next = nullptr;
						} else if (next->tailCount() > 1) {
							//head a multi tailed node, section collapse is done
							//need to add all the multiple tail edges
							n->tailEdges_.clear();
							for (const auto & tail : next->tailEdges_) {
								if(!tail->on_){
									continue;
								}
								n->tailEdges_.push_back(tail);
								n->tailEdges_.back()->head_ = n;
							}
							next = nullptr;
						} else {
							//can still collapse more set next to the tail edge
							//will only hit here if tailEdges_.size() is 1
							n->tailEdges_.clear();
							for (const auto & tail : next->tailEdges_) {
								if (!tail->on_) {
									continue;
								}
								n->tailEdges_.push_back(tail);
								n->tailEdges_.back()->head_ = n;
							}
							next = next->getFirstOnTailEdge()->tail_.lock();
						}
						nodesToEraseThread.emplace_back(nodePositions_[toErase]);
						nodes_[nodePositions_[toErase]] = nullptr;
					}
				}
				{
					std::lock_guard<std::mutex> lock(nodesToEraseMut);
					addOtherVec(nodesToErase, nodesToEraseThread);
				}
			};
			njh::concurrent::runVoidFunctionThreaded(collapsePaths, numThreads_);
		} else {
			for(auto & n : nodesToProcess){
				//auto firstTailEdge = n->getFirstOnTailEdge();
				//std::cout << "On " << n->k_ << std::endl;
				auto next = n->getFirstOnTailEdge()->tail_.lock();
				//a tail count of one indicates that the next node should be added in
				//also head count should be less than 2
				while(next != nullptr && next->headCount() == 1 && next->uid_ != n->uid_){
					////pushing backing only the back of the next k only works if there has been no processing before
					//std::cout << next->k_ << std::endl;
					//std::cout << next->k_.size() << std::endl;
					n->k_.append(next->k_.substr(klen_ - 1 ));
					n->cnt_+= next->cnt_;
					//n->inReadNamesIdx_.reserve(n->inReadNamesIdx_.size() + next->inReadNamesIdx_.size());
					n->inReadNamesIdx_.insert(next->inReadNamesIdx_.begin(), next->inReadNamesIdx_.end());
					auto toErase = next->uid_;
					if (0 == next->tailCount()) {
						//head a tailless node, end of the line
						//need to erase tail edges
						n->tailEdges_.erase(n->tailEdges_.begin());
						next = nullptr;
					} else if (next->tailCount() > 1) {
						//head a multi tailed node, section collapse is done
						//need to add all the multiple tail edges
						n->tailEdges_.clear();
						for (const auto & tail : next->tailEdges_) {
							if(!tail->on_){
								continue;
							}
							n->tailEdges_.push_back(tail);
							n->tailEdges_.back()->head_ = n;
						}
						next = nullptr;
					} else {
						//can still collapse more set next to the tail edge
						//will only hit here if tailEdges_.size() is 1
						n->tailEdges_.clear();
						for (const auto & tail : next->tailEdges_) {
							if (!tail->on_) {
								continue;
							}
							n->tailEdges_.push_back(tail);
							n->tailEdges_.back()->head_ = n;
						}
						next = next->getFirstOnTailEdge()->tail_.lock();
					}
					nodesToErase.emplace_back(nodePositions_[toErase]);
					nodes_[nodePositions_[toErase]] = nullptr;
				}
			}
		}
		if(!nodesToErase.empty()){
		//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - post process - sorting nodes to process");
			std::sort(nodesToErase.rbegin(), nodesToErase.rend());
		//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - post process - set nodes to null");
	//		for(const auto & remove : nodesToErase){
	//			nodes_[remove] = nullptr;
	//		}
		//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - post process - remove old nodes");
			for(const auto & remove : nodesToErase){
				nodes_.erase(nodes_.begin() + remove);
			}
		//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - post process - reset node positions");
			resetNodePositions();
		}
	}



	//exit(1);
//	watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - post process - remove null nodes");
//	removeNullNodes();
//	watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps() + 1, 100) + " - post process - reset node positions");
//	resetNodePositions();

	//watch.startNewLap("end");
	//std::cout << watch.toJson() << std::endl;
	//exit(1);
}


void KmerPathwayGraph::removeOffEdges(){
	//remove off edges
	for(const auto &  n : nodes_){
		{
			//remove tails
			std::vector<uint32_t> toErase;
			for( const auto & tailPos : iter::range(n->tailEdges_.size())){
				if(!n->tailEdges_[tailPos]->on_){
					toErase.emplace_back(tailPos);
				}
			}
			std::sort(toErase.rbegin(), toErase.rend());
			for(const auto pos : toErase){
				n->tailEdges_.erase(n->tailEdges_.begin() + pos);
			}
		}
		{
			//remove heads
			std::vector<uint32_t> toErase;
			for( const auto & headPos : iter::range(n->headEdges_.size())){
				if(!n->headEdges_[headPos]->on_){
					toErase.emplace_back(headPos);
				}
			}
			std::sort(toErase.rbegin(), toErase.rend());
			for(const auto pos : toErase){
				n->headEdges_.erase(n->headEdges_.begin() + pos);
			}
		}
	}
}
void KmerPathwayGraph::turnOffNodesBelowCutOff(uint32_t countCutOff){
	for(const auto &  n : nodes_){
		if(n->cnt_ <= countCutOff){
			n->on_ = false;
		}
	}
}
void KmerPathwayGraph::turnOffNodesBelowCutOff(){
	turnOffNodesBelowCutOff(occurenceCutOff_);
}
void KmerPathwayGraph::turnOffEdgesBelowCutOff(uint32_t countCutOff){
	for(const auto &  n : nodes_){
		for( auto & tail : n->tailEdges_){
			if(tail->cnt_ <= countCutOff){
				tail->on_ = false;
			}
		}
	}
}
void KmerPathwayGraph::turnOffEdgesBelowCutOff(){
	turnOffEdgesBelowCutOff(occurenceCutOff_);
}
void KmerPathwayGraph::edgeSanityCheckThrow() const {
	std::stringstream ss;
	bool failed = false;
	for (const auto & n : nodes_) {
		for (const auto pos : iter::range(n->tailEdges_.size())) {
//
//			std::cout << "pos: " << pos << std::endl;
//			std::cout << n->tailEdges_[pos]->head_.lock()->k_ << std::endl;
			if (n->tailEdges_[pos]->head_.lock()->k_ != n->k_) {
				failed = true;
				ss << "Edge doesn't make sense, "
						<< ", tail edge doesn't point actually head node " << std::endl;
				ss << "actual k " << n->k_ << ", tail edge head's k "
						<< n->tailEdges_[pos]->head_.lock()->k_ << std::endl;
			}
//
//			std::cout << "n->tailEdges_.size(): " << n->tailEdges_.size() << std::endl;
			for (const auto subPos : iter::range(pos, n->tailEdges_.size())) {
				if (pos != subPos) {
//					std::cout << "pos: " << pos << " on_ " << njh::colorBool(n->tailEdges_[pos]->on_) << std::endl;
//					std::cout << "pos: " << pos << " cnt_ " << n->tailEdges_[pos]->cnt_ << std::endl;
//					std::cout << "subPos: " << subPos << " on_" << njh::colorBool(n->tailEdges_[subPos]->on_) << std::endl;
//					std::cout << "subPos: " << subPos << " cnt_ " << n->tailEdges_[subPos]->cnt_ << std::endl;
					if (n->tailEdges_[pos]->tail_.lock()->k_
							== n->tailEdges_[subPos]->tail_.lock()->k_) {
						failed = true;
						ss << "Found duplicate edge" << std::endl;
						ss << n->uid_ << " \t " << n->tailEdges_[pos]->tail_.lock()->k_
								<< std::endl;
					}
				}
			}
		}
	}
	if(failed){
		throw std::runtime_error{ss.str()};
	}
}
void KmerPathwayGraph::removeNullNodes(){
	std::vector<uint32_t> toRemove;
	for(const auto & nodePos : iter::range(nodes_.size())){
		if(nullptr == nodes_[nodePos]){
			toRemove.emplace_back(nodePos);
		}
	}
	if(!toRemove.empty()){
		std::sort(toRemove.rbegin(), toRemove.rend());
		for(const auto & remove : toRemove){
			nodes_.erase(nodes_.begin() + remove);
		}
		resetNodePositions();
	}
}

bool KmerPathwayGraph::removeShortTips_after_disentangleInternalNodes(uint32_t shortNumber, uint32_t cntCutOff){
	bool removedTips = false;
	for (const auto & n : nodes_) {
		if(!n->on_){
			continue;
		}
		if (n->headless() && !n->tailless()) {
			auto tailNode = n->getFirstOnTailEdge()->tail_.lock();
			if(1 == n->tailCount()
					&& tailNode->tailless()
					&& 1 == tailNode->headCount()){
				if(n->k_.size() - klen_ < shortNumber && n->cnt_ < cntCutOff){
					n->getFirstOnTailEdge()->on_ = false;
					removedTips = true;
					for(const auto & otherNode : nodes_){
						if(nullptr != otherNode && otherNode->on_){
							if(tailNode->k_.size() == otherNode->k_.size()
									&& tailNode->uid_ != otherNode->uid_
									&& tailNode->k_ == otherNode->k_){
								tailNode->on_ = false;
								break;
							}
						}
					}
				}
			}
		} else if (!n->headless() && n->tailless()) {
			auto headNode = n->getFirstOnHeadEdge()->head_.lock();
			if(1 == n->headCount() && headNode->headless() && 1 == headNode->tailCount()){
				if(n->k_.size() - klen_ < shortNumber && n->cnt_ < cntCutOff){
					n->getFirstOnHeadEdge()->on_ = false;
					removedTips = true;
					for(const auto & otherNode : nodes_){
						if(nullptr != otherNode && otherNode->on_){
							if(headNode->k_.size() == otherNode->k_.size()
									&& headNode->uid_ != otherNode->uid_
									&& headNode->k_ == otherNode->k_){
								headNode->on_ = false;
								break;
							}
						}
					}
				}
			}
		}
	}
	if(removedTips){
		removeOffEdges();
		removeOffNodes();
	}
	return removedTips;
}

bool KmerPathwayGraph::removeShortTips(uint32_t shortNumber, uint32_t cntCutOff){
	bool removedTips = false;
	for (const auto & n : nodes_) {
		if (n->headless() && !n->tailless()) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << njh::bashCT::red;
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::cout << "end tip:              " << n->uid_ << std::endl;
				std::cout << "tip size:             " << n->k_.size() << std::endl;
				std::cout << "n->cnt_:              " << n->cnt_ << std::endl;
				std::cout << "n->k_.size() - klen_: " << n->k_.size() - klen_ << std::endl;
				std::cout << "readNameCount:        " << n->inReadNamesIdx_.size() << std::endl;
				std::cout << njh::bashCT::reset;
				std::cout << std::endl;
			}
#endif
			//if(n->k_.size() - klen_ < shortNumber && n->cnt_ < cntCutOff){
			if(n->k_.size() - klen_ < shortNumber && n->inReadNamesIdx_.size() < cntCutOff){
				for(auto & e : n->tailEdges_){
					e->on_ = false;
					removedTips = true;
				}
#if defined(PATHWEAVERSUPERDEBUG)
				{
					std::cout << njh::bashCT::cyan;
					std::cout << __FILE__ << " " << __LINE__ << std::endl;
					std::cout << "end tip:              " << n->uid_ << std::endl;
					std::cout << "tip size:             " << n->k_.size() << std::endl;
					std::cout << "n->cnt_:              " << n->cnt_ << std::endl;
					std::cout << "n->k_.size() - klen_: " << n->k_.size() - klen_ << std::endl;
					std::cout << "readNameCount:        " << n->inReadNamesIdx_.size() << std::endl;
					std::cout << njh::bashCT::reset;
					std::cout << std::endl;
				}
#endif
			}
		} else if (!n->headless() && n->tailless()) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << njh::bashCT::red;
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::cout << "end tip:              " << n->uid_ << std::endl;
				std::cout << "tip size:             " << n->k_.size() << std::endl;
				std::cout << "n->cnt_:              " << n->cnt_ << std::endl;
				std::cout << "n->k_.size() - klen_: " << n->k_.size() - klen_ << std::endl;
				std::cout << "readNameCount:        " << n->inReadNamesIdx_.size() << std::endl;
				std::cout << njh::bashCT::reset;
				std::cout << std::endl;
			}
#endif
			//if(n->k_.size() - klen_ < shortNumber && n->cnt_ < cntCutOff){
			if(n->k_.size() - klen_ < shortNumber && n->inReadNamesIdx_.size() < cntCutOff){
				for(auto & e : n->headEdges_){
					e->on_ = false;
					removedTips = true;
				}
#if defined(PATHWEAVERSUPERDEBUG)
				{
					std::cout << njh::bashCT::cyan;
					std::cout << __FILE__ << " " << __LINE__ << std::endl;
					std::cout << "end tip:              " << n->uid_ << std::endl;
					std::cout << "tip size:             " << n->k_.size() << std::endl;
					std::cout << "n->cnt_:              " << n->cnt_ << std::endl;
					std::cout << "n->k_.size() - klen_: " << n->k_.size() - klen_ << std::endl;
					std::cout << "readNameCount:        " << n->inReadNamesIdx_.size() << std::endl;
					std::cout << njh::bashCT::reset;
					std::cout << std::endl;
				}
#endif
			}
		}
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	if(removedTips){
		removeOffEdges();
	}
	return removedTips;
}
void KmerPathwayGraph::removeOffNodes(){
	std::vector<uint32_t> toRemove;
	//njh::stopWatch watch;
	//watch.setLapName(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps(), 100) + " - find off nodes");
	for(const auto & nodePos : iter::range(nodes_.size())){
		const auto & n = nodes_[nodePos];
		if(!n->on_){
			toRemove.emplace_back(nodePos);
			//turn of any edges that might be pointing to this node
			for(auto & head : n->headEdges_){
				head->on_ = false;
			}
			for(auto & tail : n->tailEdges_){
				tail->on_ = false;
			}
		}
	}

	if(!toRemove.empty()){
		//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps(), 100) + " - erase nodes");
		std::sort(toRemove.rbegin(), toRemove.rend());
		for(const auto & remove : toRemove){
			nodes_.erase(nodes_.begin() + remove);
		}
		//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps(), 100) + " - remove off edges");
		removeOffEdges();
		//watch.startNewLap(njh::leftPadNumStr<uint32_t>(watch.getNumberOfLaps(), 100) + " - reset node positions");
		resetNodePositions();
	}
	//watch.startNewLap("end");

	//std::cout << watch.toJson() << std::endl;
}

bool KmerPathwayGraph::removeHeadlessTaillessNodes(){
	return removeHeadlessTaillessNodes(klen_);
}

bool KmerPathwayGraph::removeHeadlessTaillessNodes(uint32_t lenCutOff){
	std::vector<uint32_t> toRemove;
	for(const auto & nodePos : iter::range(nodes_.size())){
		const auto & n = nodes_[nodePos];
		if(n->headless() && n->tailless() && n->k_.length() <=lenCutOff){
			toRemove.emplace_back(nodePos);
		}
	}
	if(!toRemove.empty()){
		std::sort(toRemove.rbegin(), toRemove.rend());
		for(const auto & remove : toRemove){
			nodes_.erase(nodes_.begin() + remove);
		}
		resetNodePositions();
	}
	return !toRemove.empty();
}

void KmerPathwayGraph::resetNodeUids(){
	std::unordered_map<std::string, uint32_t> previouslyFound;
	for( auto & node : nodes_){
		if(njh::in(node->k_, previouslyFound)){
			++previouslyFound[node->k_];
			node->uid_ = node->k_ + estd::to_string(previouslyFound[node->k_]);
		}else{
			previouslyFound[node->k_] = 0;
			node->uid_ = node->k_;
		}
	}
	nodePositions_.clear();
}
void KmerPathwayGraph::resetNodePositions(){
	resetNodeUids();
	nodePositions_.clear();
	for(const auto pos : iter::range(nodes_.size())){
		nodePositions_[nodes_[pos]->uid_] = pos;
	}
}
bool KmerPathwayGraph::splitMultitailedNodes() {
	resetNodePositions();
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	for (const auto & n : nodes_) {
		if (n->tailCount() > 1) {
			nodesToProcess.push_back(n);
		}
	}
	if (nodesToProcess.empty()) {
		return false;
	}
#if defined(PATHWEAVERSUPERDEBUG)
	{
		std::cout << "There are " << nodesToProcess.size() << " nodes to process" << std::endl;
	}
#endif
	bool nodesSplit = false;
	for (const auto & n : nodesToProcess) {
		if (n->tailCount() > 1) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << std::endl;
				std::cout << n->uid_ << std::endl;
			}
#endif
			{
				//split tails;
				std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t> > outgoing;
				for (const auto & tail : n->tailEdges_) {
					if(!tail->on_){
						continue;
					}
					auto tailNode = tail->tail_.lock();
					for (const auto & readName : tail->inReadNamesIdx_) {
						++outgoing[tailNode->k_][readName];
					}
				}
#if defined(PATHWEAVERSUPERDEBUG)
				{
					table outgoingTab(outgoing, VecStr { "head", "readNames", "count" });
					outgoingTab.addColumn(VecStr { "outgoing" }, "direction");
					std::ofstream outOut("outgoing.tab.txt");
					outgoingTab.outPutContents(outOut, "\t");
				}
#endif
				if (!outgoing.empty()) {
					nodesSplit = true;
					auto tails = getVectorOfMapKeys(outgoing);
					for (const auto & tail : tails) {
						std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
								KmerPathwayGraph::node>(n->k_, n->cnt_, n->kLen_);
						std::unordered_set<uint32_t> names;
						for (const auto & name : n->inReadNamesIdx_) {
							bool appearsInTails = false;
							bool appearsInOtherTails = false;
							if (njh::in(name, outgoing[tail])) {
								appearsInTails = true;
							}
							if (!appearsInTails) {
								continue;
							}
							for (const auto & other : outgoing) {
								if (other.first != tail && njh::in(name, other.second)) {
									appearsInOtherTails = true;
#if defined(PATHWEAVERSUPERDEBUG)
									{
										if (appearsInTails) {
											std::cout << name << " appears in other tails"
													<< std::endl;
											std::cout << "other.first != tail "
													<< njh::colorBool(other.first != tail) << std::endl;
											std::cout << "njh::in(name, other.second) "
													<< njh::colorBool(njh::in(name, other.second))
													<< std::endl;
											std::cout << "other: " << other.first << std::endl;
											std::cout << "currentTail: " << tail << std::endl;
											std::cout << "appearsInTails: "
													<< njh::colorBool(appearsInTails) << std::endl;
											printOutMapContents(other.second, "\t", std::cout);
											throw std::runtime_error { "stopping" };
										}
									}
#endif
								}
							}
							if (!appearsInOtherTails && appearsInTails) {
								names.emplace(name);
							}
						}
						//add read names
						newNode->inReadNamesIdx_ = names;
						//add edges and modifying head and tail nodes
						for (auto & tailEdge : n->tailEdges_) {
							if(!tailEdge->on_){
								continue;
							}
							if (tail == tailEdge->tail_.lock()->k_) {
								std::unordered_set<uint32_t> edgeNames;
								for (const auto & edgeName : tailEdge->inReadNamesIdx_) {
									if (njh::in(edgeName, newNode->inReadNamesIdx_)) {
										edgeNames.emplace(edgeName);
									}
								}
								auto addingEdge = std::make_shared<KmerPathwayGraph::edge>(
										newNode, tailEdge->tail_.lock(), edgeNames.size(),
										edgeNames);
								newNode->tailEdges_.push_back(addingEdge);
								tailEdge->tail_.lock()->headEdges_.push_back(addingEdge);
							}
						}
						//add node to graph.nodes_
						nodes_.push_back(newNode);
						//printVector(n->readNames_);
						//printVector(newNode->readNames_	);
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << __PRETTY_FUNCTION__ << std::endl;
							std::cout << njh::bashCT::cyan << "Adding new node "
									<< newNode->k_ << " " << newNode->k_ << njh::bashCT::reset
									<< std::endl;
							std::cout << njh::bashCT::purple << "With head edges:  ";
							for (const auto & head : newNode->headEdges_) {
								std::cout << head->head_.lock()->k_ << " ";
							}
							std::cout << njh::bashCT::reset << std::endl;
							std::cout << njh::bashCT::green << "With tail edges:  ";
							for (const auto & tail : newNode->tailEdges_) {
								std::cout << tail->tail_.lock()->k_ << " ";
							}
							std::cout << njh::bashCT::reset << std::endl;
							if (newNode->headCount() == 1 && newNode->tailCount() == 1) {
								std::cout << "singly linked " << newNode->k_ << " "
										<< newNode->k_ << std::endl;
							}
						}
#endif
					}
					//throw std::runtime_error{"stopping"};
					//turn off node and edges
					for (auto & tailEdge : n->tailEdges_) {
						tailEdge->on_ = false;
					}
					n->on_ = false;
				}
			}
			{
				//split heads
				std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t> > incoming;
				for (const auto & head : n->headEdges_) {
					if(!head->on_){
						continue;
					}
					auto headNode = head->head_.lock();
					for (const auto & readName : head->inReadNamesIdx_) {
						++incoming[headNode->k_][readName];
					}
				}
#if defined(PATHWEAVERSUPERDEBUG)
				{
					table incommingTab(incoming, VecStr { "head", "readNames", "count" });
					incommingTab.addColumn(VecStr { "incomming" }, "direction");
					std::ofstream outOut("incomming.tab.txt");
					incommingTab.outPutContents(outOut, "\t");
				}
#endif
				if (!incoming.empty()) {
					nodesSplit = true;
					auto heads = getVectorOfMapKeys(incoming);
					for (const auto & head : heads) {
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << "head: " << head << std::endl;
						}
#endif
						std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
								KmerPathwayGraph::node>(n->k_, n->cnt_, n->kLen_);
						std::unordered_set<uint32_t> names;
						for (const auto & name : n->inReadNamesIdx_) {
							bool appearsInHeads = false;
							bool appearsInOtherHeads = false;
							if (njh::in(name, incoming[head])) {
								appearsInHeads = true;
							}
							if (!appearsInHeads) {
								continue;
							}
							for (const auto & other : incoming) {
								if (other.first != head && njh::in(name, other.second)) {
									appearsInOtherHeads = true;
#if defined(PATHWEAVERSUPERDEBUG)
									{
										if (appearsInHeads) {
											std::cout << name << " appears in other heads" << std::endl;
											std::cout << "other.first != head "
													<< njh::colorBool(other.first != head) << std::endl;
											std::cout << "njh::in(name, other.second) "
													<< njh::colorBool(njh::in(name, other.second))
													<< std::endl;
											std::cout << "other: " << other.first << std::endl;
											std::cout << "currentHead: " << head << std::endl;
											std::cout << "appearsInHeads: "
													<< njh::colorBool(appearsInHeads) << std::endl;
											printOutMapContents(other.second, "\t", std::cout);
											throw std::runtime_error { "stopping" };
										}
									}
#endif
								}
							}
							if (!appearsInOtherHeads && appearsInHeads) {
								names.emplace(name);
							}
						}
						//throw std::runtime_error{"stopping"};
						//add read names
						newNode->inReadNamesIdx_ = names;
						//add edges and modifying head and tail nodes
						for (auto & headEdge : n->headEdges_) {
							if(!headEdge->on_){
								continue;
							}
							if (head == headEdge->head_.lock()->k_) {
								std::unordered_set<uint32_t> edgeNames;
								for (const auto & edgeName : headEdge->inReadNamesIdx_) {
									if (njh::in(edgeName, newNode->inReadNamesIdx_)) {
										edgeNames.emplace(edgeName);
									}
								}
								auto addingEdge = std::make_shared<KmerPathwayGraph::edge>(
										headEdge->head_.lock(), newNode, edgeNames.size(), edgeNames);
								newNode->headEdges_.push_back(addingEdge);
								headEdge->head_.lock()->tailEdges_.push_back(addingEdge);
							}
						}
						//add node to graph.nodes_
						nodes_.push_back(newNode);
						/*
						 *
						 printVector(n->readNames_);
						 printVector(newNode->readNames_	);
						 */
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << __PRETTY_FUNCTION__ << std::endl;
							std::cout << njh::bashCT::cyan << "Adding new node "
									<< newNode->k_ << " " << newNode->k_ << njh::bashCT::reset
									<< std::endl;
							std::cout << njh::bashCT::purple << "With head edges:  ";
							for (const auto & head : newNode->headEdges_) {
								std::cout << head->head_.lock()->k_ << " ";
							}
							std::cout << njh::bashCT::reset << std::endl;
							std::cout << njh::bashCT::green << "With tail edges:  ";
							for (const auto & tail : newNode->tailEdges_) {
								std::cout << tail->tail_.lock()->k_ << " ";
							}
							std::cout << njh::bashCT::reset << std::endl;
							if (newNode->headCount() == 1 && newNode->tailCount() == 1) {
								std::cout << "singly linked " << newNode->k_ << " "
										<< newNode->k_ << std::endl;
							}
						}
#endif
					}
					//throw std::runtime_error{"stopping"};
					//turn off node and edges
					for (auto & headEdge : n->headEdges_) {
						headEdge->on_ = false;
					}
					n->on_ = false;
				}
			}
		}
		//remove off nodes;
		removeOffNodes();
		removeOffEdges();
	}
	resetNodePositions();
	return nodesSplit;
}

bool KmerPathwayGraph::splitEndNodes() {
	resetNodePositions();
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	for (const auto & n : nodes_) {
		if (n->headless() && n->tailCount() > 1) {
			nodesToProcess.push_back(n);
		} else if (n->tailless() && n->headCount() > 1) {
			nodesToProcess.push_back(n);
		}
	}
	if (nodesToProcess.empty()) {
		return false;
	}
#if defined(PATHWEAVERSUPERDEBUG)
	{
		std::cout << "There are " << nodesToProcess.size() << " nodes to process" << std::endl;
	}
#endif
	bool nodesSplit = false;
	for (const auto & n : nodesToProcess) {
		if (n->headless() && n->tailCount() > 1) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << std::endl;
				std::cout << n->uid_ << std::endl;
			}
#endif
			std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t> > outgoing;
			for (const auto & tail : n->tailEdges_) {
				if(!tail->on_){
					continue;
				}
				auto tailNode = tail->tail_.lock();
				for (const auto & readName : tail->inReadNamesIdx_) {
					++outgoing[tailNode->uid_][readName];
				}
			}
#if defined(PATHWEAVERSUPERDEBUG)
			{
				table outgoingTab(outgoing, VecStr { "head", "readNames", "count" });
				outgoingTab.addColumn(VecStr { "outgoing" }, "direction");
				std::ofstream outOut("outgoing.tab.txt");
				outgoingTab.outPutContents(outOut, "\t");
			}
#endif
			if (!outgoing.empty()) {
				nodesSplit = true;
				auto tails = getVectorOfMapKeys(outgoing);
				for (const auto & tail : tails) {
					std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
							KmerPathwayGraph::node>(n->k_, n->cnt_, n->kLen_);
					std::unordered_set<uint32_t> names;
					for (const auto & name : n->inReadNamesIdx_) {
						bool appearsInTails = false;
						bool appearsInOtherTails = false;
						if (njh::in(name, outgoing[tail])) {
							appearsInTails = true;
						}
						if (!appearsInTails) {
							continue;
						}
						for (const auto & other : outgoing) {
							if (other.first != tail && njh::in(name, other.second)) {
								appearsInOtherTails = true;
#if defined(PATHWEAVERSUPERDEBUG)
								{
									if (appearsInTails) {
										std::cout << name << " appears in other tails" << std::endl;
										std::cout << "other.first != tail "
												<< njh::colorBool(other.first != tail) << std::endl;
										std::cout << "njh::in(name, other.second) "
												<< njh::colorBool(njh::in(name, other.second))
												<< std::endl;
										std::cout << "other: " << other.first << std::endl;
										std::cout << "currentTail: " << tail << std::endl;
										std::cout << "appearsInTails: "
												<< njh::colorBool(appearsInTails) << std::endl;
										printOutMapContents(other.second, "\t", std::cout);
										throw std::runtime_error { "stopping" };
									}
								}
#endif
							}
						}
						if (!appearsInOtherTails && appearsInTails) {
							names.emplace(name);
						}
					}
						//add read names
						newNode->inReadNamesIdx_ = names;
						//add edges and modifying head and tail nodes
						for (auto & tailEdge : n->tailEdges_) {
							if(!tailEdge->on_){
								continue;
							}
							if (tail == tailEdge->tail_.lock()->uid_) {
								std::unordered_set<uint32_t> edgeNames;
								for (const auto & edgeName : tailEdge->inReadNamesIdx_) {
									if (njh::in(edgeName, newNode->inReadNamesIdx_)) {
										edgeNames.emplace(edgeName);
									}
								}
								auto addingEdge = std::make_shared<KmerPathwayGraph::edge>(
										newNode, tailEdge->tail_.lock(), edgeNames.size(),
										edgeNames);
								newNode->tailEdges_.push_back(addingEdge);
								tailEdge->tail_.lock()->headEdges_.push_back(addingEdge);
							}
						}
						//add node to graph.nodes_
						nodes_.push_back(newNode);
						//printVector(n->readNames_);
						//printVector(newNode->readNames_	);
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << __PRETTY_FUNCTION__ << std::endl;
							std::cout << njh::bashCT::cyan << "Adding new node " << newNode->k_
									<< " " << newNode->k_ << njh::bashCT::reset << std::endl;
							std::cout << njh::bashCT::purple << "With head edges:  ";
							for (const auto & head : newNode->headEdges_) {
								std::cout << head->head_.lock()->k_ << " ";
							}
							std::cout << njh::bashCT::reset << std::endl;
							std::cout << njh::bashCT::green << "With tail edges:  ";
							for (const auto & tail : newNode->tailEdges_) {
								std::cout << tail->tail_.lock()->k_ << " ";
							}
							std::cout << njh::bashCT::reset << std::endl;
							if (newNode->headCount() == 1 && newNode->tailCount() == 1) {
								std::cout << "singly linked " << newNode->k_ << " " << newNode->k_
										<< std::endl;
							}
						}
#endif
				}
				//throw std::runtime_error{"stopping"};
				//turn off node and edges
				for (auto & headEdge : n->headEdges_) {
					headEdge->on_ = false;
				}
				for (auto & tailEdge : n->tailEdges_) {
					tailEdge->on_ = false;
				}
				n->on_ = false;
			}
		} else if (n->tailless() && n->headCount() > 1) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << std::endl;
				std::cout << n->uid_ << std::endl;
			}
#endif
			std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t> > incoming;
			for (const auto & head : n->headEdges_) {
				if(!head->on_){
					continue;
				}
				auto headNode = head->head_.lock();
				for (const auto & readName : head->inReadNamesIdx_) {
					++incoming[headNode->uid_][readName];
				}
			}
#if defined(PATHWEAVERSUPERDEBUG)
			{
				table incommingTab(incoming, VecStr { "head", "readNames", "count" });
				incommingTab.addColumn(VecStr { "incomming" }, "direction");
				std::ofstream outOut("incomming.tab.txt");
				incommingTab.outPutContents(outOut, "\t");
			}
#endif
			if (!incoming.empty()) {
				nodesSplit = true;
				auto heads = getVectorOfMapKeys(incoming);
				for (const auto & head : heads) {
#if defined(PATHWEAVERSUPERDEBUG)
					{
						std::cout << "head: " << head << std::endl;
					}
#endif
					std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
							KmerPathwayGraph::node>(n->k_, n->cnt_, n->kLen_);
					std::unordered_set<uint32_t> names;
					for (const auto & name : n->inReadNamesIdx_) {
						bool appearsInHeads = false;
						bool appearsInOtherHeads = false;
						if (njh::in(name, incoming[head])) {
							appearsInHeads = true;
						}
						if (!appearsInHeads) {
							continue;
						}
						for (const auto & other : incoming) {
							if (other.first != head && njh::in(name, other.second)) {
								appearsInOtherHeads = true;
#if defined(PATHWEAVERSUPERDEBUG)
								{
									if (appearsInHeads) {
										std::cout << name << " appears in other heads" << std::endl;
										std::cout << "other.first != head "
												<< njh::colorBool(other.first != head) << std::endl;
										std::cout << "njh::in(name, other.second) "
												<< njh::colorBool(njh::in(name, other.second))
												<< std::endl;
										std::cout << "other: " << other.first << std::endl;
										std::cout << "currentHead: " << head << std::endl;
										std::cout << "appearsInHeads: "
												<< njh::colorBool(appearsInHeads) << std::endl;
										printOutMapContents(other.second, "\t", std::cout);
										throw std::runtime_error { "stopping" };
									}
								}
#endif
							}
						}
						if (!appearsInOtherHeads && appearsInHeads) {
							names.emplace(name);
						}
					}
					//throw std::runtime_error{"stopping"};
					//add read names
					newNode->inReadNamesIdx_ = names;
					//add edges and modifying head and tail nodes
					for (auto & headEdge : n->headEdges_) {
						if(!headEdge->on_){
							continue;
						}
						if (head == headEdge->head_.lock()->uid_) {
							std::unordered_set<uint32_t> edgeNames;
							for (const auto & edgeName : headEdge->inReadNamesIdx_) {
								if (njh::in(edgeName, newNode->inReadNamesIdx_)) {
									edgeNames.emplace(edgeName);
								}
							}
							auto addingEdge = std::make_shared<KmerPathwayGraph::edge>(
									headEdge->head_.lock(), newNode, edgeNames.size(), edgeNames);
							newNode->headEdges_.push_back(addingEdge);
							headEdge->head_.lock()->tailEdges_.push_back(addingEdge);
						}
					}
					//add node to graph.nodes_
					nodes_.push_back(newNode);
					/*
					 *
					 printVector(n->readNames_);
					 printVector(newNode->readNames_	);
					 */
#if defined(PATHWEAVERSUPERDEBUG)
					{
						std::cout << __FILE__ << " " << __LINE__ << std::endl;
						std::cout << __PRETTY_FUNCTION__ << std::endl;
						std::cout << njh::bashCT::cyan << "Adding new node " << newNode->k_
								<< " " << newNode->k_ << njh::bashCT::reset << std::endl;
						std::cout << njh::bashCT::purple << "With head edges:  ";
						for (const auto & head : newNode->headEdges_) {
							std::cout << head->head_.lock()->k_ << " ";
						}
						std::cout << njh::bashCT::reset << std::endl;
						std::cout << njh::bashCT::green << "With tail edges:  ";
						for (const auto & tail : newNode->tailEdges_) {
							std::cout << tail->tail_.lock()->k_ << " ";
						}
						std::cout << njh::bashCT::reset << std::endl;
						if (newNode->headCount() == 1 && newNode->tailCount() == 1) {
							std::cout << "singly linked " << newNode->k_ << " " << newNode->k_
									<< std::endl;
						}
					}
#endif
				}
				//throw std::runtime_error{"stopping"};
				//turn off node and edges
				for (auto & headEdge : n->headEdges_) {
					headEdge->on_ = false;
				}
				for (auto & tailEdge : n->tailEdges_) {
					tailEdge->on_ = false;
				}
				n->on_ = false;
			}
		}
		//remove off nodes;
		removeOffNodes();
		removeOffEdges();
	}
	resetNodePositions();
	/*{
		std::ofstream outTestDot("outTest.dot");
		writeDot(outTestDot);
		throw std::runtime_error{"stopping"};
	}*/
	return nodesSplit;
}
bool KmerPathwayGraph::splitEndNodes(uint32_t maxLen) {
	resetNodePositions();
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	for (const auto & n : nodes_) {
		if(n->k_.size() <= maxLen){
			if (n->headless() && n->tailCount() > 1) {
				nodesToProcess.push_back(n);
			} else if (n->tailless() && n->headCount() > 1) {
				nodesToProcess.push_back(n);
			}
		}
	}
	if (nodesToProcess.empty()) {
		return false;
	}
#if defined(PATHWEAVERSUPERDEBUG)
	{
		std::cout << "There are " << nodesToProcess.size() << " nodes to process" << std::endl;
	}
#endif
	bool nodesSplit = false;
	for (const auto & n : nodesToProcess) {
		if (n->headless() && n->tailCount() > 1) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << std::endl;
				std::cout << n->uid_ << std::endl;
			}
#endif
			std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t> > outgoing;
			for (const auto & tail : n->tailEdges_) {
				if(!tail->on_){
					continue;
				}
				auto tailNode = tail->tail_.lock();
				for (const auto & readName : tail->inReadNamesIdx_) {
					++outgoing[tailNode->uid_][readName];
				}
			}
#if defined(PATHWEAVERSUPERDEBUG)
			{
				table outgoingTab(outgoing, VecStr { "head", "readNames", "count" });
				outgoingTab.addColumn(VecStr { "outgoing" }, "direction");
				std::ofstream outOut("outgoing.tab.txt");
				outgoingTab.outPutContents(outOut, "\t");
			}
#endif
			if (!outgoing.empty()) {
				nodesSplit = true;
				auto tails = getVectorOfMapKeys(outgoing);
				for (const auto & tail : tails) {
					std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
							KmerPathwayGraph::node>(n->k_, n->cnt_, n->kLen_);
					std::unordered_set<uint32_t> names;
					for (const auto & name : n->inReadNamesIdx_) {
						bool appearsInTails = false;
						bool appearsInOtherTails = false;
						if (njh::in(name, outgoing[tail])) {
							appearsInTails = true;
						}
						if (!appearsInTails) {
							continue;
						}
						for (const auto & other : outgoing) {
							if (other.first != tail && njh::in(name, other.second)) {
								appearsInOtherTails = true;
#if defined(PATHWEAVERSUPERDEBUG)
								{
									if (appearsInTails) {
										std::cout << name << " appears in other tails" << std::endl;
										std::cout << "other.first != tail "
												<< njh::colorBool(other.first != tail) << std::endl;
										std::cout << "njh::in(name, other.second) "
												<< njh::colorBool(njh::in(name, other.second))
												<< std::endl;
										std::cout << "other: " << other.first << std::endl;
										std::cout << "currentTail: " << tail << std::endl;
										std::cout << "appearsInTails: "
												<< njh::colorBool(appearsInTails) << std::endl;
										printOutMapContents(other.second, "\t", std::cout);
										throw std::runtime_error { "stopping" };
									}
								}
#endif
							}
						}
						if (!appearsInOtherTails && appearsInTails) {
							names.emplace(name);
						}
					}
						//add read names
						newNode->inReadNamesIdx_ = names;
						//add edges and modifying head and tail nodes
						for (auto & tailEdge : n->tailEdges_) {
							if(!tailEdge->on_){
								continue;
							}
							if (tail == tailEdge->tail_.lock()->uid_) {
								std::unordered_set<uint32_t> edgeNames;
								for (const auto & edgeName : tailEdge->inReadNamesIdx_) {
									if (njh::in(edgeName, newNode->inReadNamesIdx_)) {
										edgeNames.emplace(edgeName);
									}
								}
								auto addingEdge = std::make_shared<KmerPathwayGraph::edge>(
										newNode, tailEdge->tail_.lock(), edgeNames.size(),
										edgeNames);
								newNode->tailEdges_.push_back(addingEdge);
								tailEdge->tail_.lock()->headEdges_.push_back(addingEdge);
							}
						}
						//add node to graph.nodes_
						nodes_.push_back(newNode);
						//printVector(n->readNames_);
						//printVector(newNode->readNames_	);
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << __PRETTY_FUNCTION__ << std::endl;
							std::cout << njh::bashCT::cyan << "Adding new node " << newNode->k_
									<< " " << newNode->k_ << njh::bashCT::reset << std::endl;
							std::cout << njh::bashCT::purple << "With head edges:  ";
							for (const auto & head : newNode->headEdges_) {
								std::cout << head->head_.lock()->k_ << " ";
							}
							std::cout << njh::bashCT::reset << std::endl;
							std::cout << njh::bashCT::green << "With tail edges:  ";
							for (const auto & tail : newNode->tailEdges_) {
								std::cout << tail->tail_.lock()->k_ << " ";
							}
							std::cout << njh::bashCT::reset << std::endl;
							if (newNode->headCount() == 1 && newNode->tailCount() == 1) {
								std::cout << "singly linked " << newNode->k_ << " " << newNode->k_
										<< std::endl;
							}
						}
#endif
				}
				//throw std::runtime_error{"stopping"};
				//turn off node and edges
				for (auto & headEdge : n->headEdges_) {
					headEdge->on_ = false;
				}
				for (auto & tailEdge : n->tailEdges_) {
					tailEdge->on_ = false;
				}
				n->on_ = false;
			}
		} else if (n->tailless() && n->headCount() > 1) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << std::endl;
				std::cout << n->uid_ << std::endl;
			}
#endif
			std::unordered_map<std::string, std::unordered_map<uint32_t, uint32_t> > incoming;
			for (const auto & head : n->headEdges_) {
				if(!head->on_){
					continue;
				}
				auto headNode = head->head_.lock();
				for (const auto & readName : head->inReadNamesIdx_) {
					++incoming[headNode->uid_][readName];
				}
			}
#if defined(PATHWEAVERSUPERDEBUG)
			{
				table incommingTab(incoming, VecStr { "head", "readNames", "count" });
				incommingTab.addColumn(VecStr { "incomming" }, "direction");
				std::ofstream outOut("incomming.tab.txt");
				incommingTab.outPutContents(outOut, "\t");
			}
#endif
			if (!incoming.empty()) {
				nodesSplit = true;
				auto heads = getVectorOfMapKeys(incoming);
				for (const auto & head : heads) {
#if defined(PATHWEAVERSUPERDEBUG)
					{
						std::cout << "head: " << head << std::endl;
					}
#endif
					std::shared_ptr<KmerPathwayGraph::node> newNode = std::make_shared<
							KmerPathwayGraph::node>(n->k_, n->cnt_, n->kLen_);
					std::unordered_set<uint32_t> names;
					for (const auto & name : n->inReadNamesIdx_) {
						bool appearsInHeads = false;
						bool appearsInOtherHeads = false;
						if (njh::in(name, incoming[head])) {
							appearsInHeads = true;
						}
						if (!appearsInHeads) {
							continue;
						}
						for (const auto & other : incoming) {
							if (other.first != head && njh::in(name, other.second)) {
								appearsInOtherHeads = true;
#if defined(PATHWEAVERSUPERDEBUG)
								{
									if (appearsInHeads) {
										std::cout << name << " appears in other heads" << std::endl;
										std::cout << "other.first != head "
												<< njh::colorBool(other.first != head) << std::endl;
										std::cout << "njh::in(name, other.second) "
												<< njh::colorBool(njh::in(name, other.second))
												<< std::endl;
										std::cout << "other: " << other.first << std::endl;
										std::cout << "currentHead: " << head << std::endl;
										std::cout << "appearsInHeads: "
												<< njh::colorBool(appearsInHeads) << std::endl;
										printOutMapContents(other.second, "\t", std::cout);
										throw std::runtime_error { "stopping" };
									}
								}
#endif
							}
						}
						if (!appearsInOtherHeads && appearsInHeads) {
							names.emplace(name);
						}
					}
					//throw std::runtime_error{"stopping"};
					//add read names
					newNode->inReadNamesIdx_ = names;
					//add edges and modifying head and tail nodes
					for (auto & headEdge : n->headEdges_) {
						if(!headEdge->on_){
							continue;
						}
						if (head == headEdge->head_.lock()->uid_) {
							std::unordered_set<uint32_t> edgeNames;
							for (const auto & edgeName : headEdge->inReadNamesIdx_) {
								if (njh::in(edgeName, newNode->inReadNamesIdx_)) {
									edgeNames.emplace(edgeName);
								}
							}
							auto addingEdge = std::make_shared<KmerPathwayGraph::edge>(
									headEdge->head_.lock(), newNode, edgeNames.size(), edgeNames);
							newNode->headEdges_.push_back(addingEdge);
							headEdge->head_.lock()->tailEdges_.push_back(addingEdge);
						}
					}
					//add node to graph.nodes_
					nodes_.push_back(newNode);
					/*
					 *
					 printVector(n->readNames_);
					 printVector(newNode->readNames_	);
					 */
#if defined(PATHWEAVERSUPERDEBUG)
					{
						std::cout << __FILE__ << " " << __LINE__ << std::endl;
						std::cout << __PRETTY_FUNCTION__ << std::endl;
						std::cout << njh::bashCT::cyan << "Adding new node " << newNode->k_
								<< " " << newNode->k_ << njh::bashCT::reset << std::endl;
						std::cout << njh::bashCT::purple << "With head edges:  ";
						for (const auto & head : newNode->headEdges_) {
							std::cout << head->head_.lock()->k_ << " ";
						}
						std::cout << njh::bashCT::reset << std::endl;
						std::cout << njh::bashCT::green << "With tail edges:  ";
						for (const auto & tail : newNode->tailEdges_) {
							std::cout << tail->tail_.lock()->k_ << " ";
						}
						std::cout << njh::bashCT::reset << std::endl;
						if (newNode->headCount() == 1 && newNode->tailCount() == 1) {
							std::cout << "singly linked " << newNode->k_ << " " << newNode->k_
									<< std::endl;
						}
					}
#endif
				}
				//throw std::runtime_error{"stopping"};
				//turn off node and edges
				for (auto & headEdge : n->headEdges_) {
					headEdge->on_ = false;
				}
				for (auto & tailEdge : n->tailEdges_) {
					tailEdge->on_ = false;
				}
				n->on_ = false;
			}
		}
		//remove off nodes;
		removeOffNodes();
		removeOffEdges();
	}
	resetNodePositions();
	/*{
		std::ofstream outTestDot("outTest.dot");
		writeDot(outTestDot);
		throw std::runtime_error{"stopping"};
	}*/
	return nodesSplit;
}


bool KmerPathwayGraph::hasSelfPointingPaths(){
	resetNodePositions();
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	for (const auto & n : nodes_) {
		//check for nodes with 1 head node that points to it's self
		if(n->headCount() == 1 ){
			for( auto & head : n->headEdges_){
				if(head->on_){
					if(head->head_.lock()->uid_ == n->uid_){
						//this node has one head and it's head node is this node
						return true;
					}
				}
			}
		}
		//check tail nodes that point to it's self
		if(n->tailCount() == 1 ){
			for( auto & tail : n->tailEdges_){
				if(tail->on_){
					if(tail->tail_.lock()->uid_ == n->uid_){
						//this node has one tail and it's tail node is this node
						return true;
					}
				}
			}
		}
	}

	return false;
}



bool KmerPathwayGraph::breakSelfPointingPathsKeepOtherEdges() {
	resetNodePositions();
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	bool brokeConnections = false;
	for (const auto & n : nodes_) {
		//check for nodes with 1 head node that points to it's self
		if(n->headCount() == 1 ){
			for( auto & head : n->headEdges_){
				if(head->on_){
					if(head->head_.lock()->uid_ == n->uid_){
						//this node has one head nad it's head node is this node turn of this connection
						head->on_ = false;
						brokeConnections = true;
					}
					break;
				}
			}
		}
		//check tail nodes that point to it's self
		if(n->tailCount() == 1 ){
			for( auto & tail : n->tailEdges_){
				if(tail->on_){
					if(tail->tail_.lock()->uid_ == n->uid_){
						//this node has one head nad it's head node is this node turn of this connection
						tail->on_ = false;
						brokeConnections = true;
					}
					break;
				}
			}
		}
	}
	if(brokeConnections){
		removeOffEdges();
		return true;
	}
	return false;
}

bool KmerPathwayGraph::breakSelfPointingPaths() {
	resetNodePositions();
	std::vector<std::shared_ptr<KmerPathwayGraph::node>> nodesToProcess;
	bool brokeConnections = false;
	for (const auto & n : nodes_) {
		bool nodeHasSelfPointing = false;
		//check for nodes with 1 head node that points to it's self
		if(n->headCount() == 1 ){
			for( auto & head : n->headEdges_){
				if(head->on_){
					if(head->head_.lock()->uid_ == n->uid_){
						//this node has one head nad it's head node is this node turn of this connection
						head->on_ = false;
						brokeConnections = true;
						nodeHasSelfPointing = true;
					}
					break;
				}
			}
		}
		//check tail nodes that point to it's self
		if(n->tailCount() == 1 ){
			for( auto & tail : n->tailEdges_){
				if(tail->on_){
					if(tail->tail_.lock()->uid_ == n->uid_){
						//this node has one head nad it's head node is this node turn of this connection
						tail->on_ = false;
						brokeConnections = true;
						nodeHasSelfPointing = true;
					}
					break;
				}
			}
		}
		//if node had any self pointing paths, break all other edges
		if(nodeHasSelfPointing){
			for(auto & head : n->headEdges_){
				head->on_ = false;
			}
			for(auto & tail : n->tailEdges_){
				tail->on_ = false;
			}
		}
	}
	if(brokeConnections){
		removeOffEdges();
		return true;
	}
	return false;
}

table KmerPathwayGraph::getNodeStateCounts() const {
	uint32_t tailless = 0;
	uint32_t headless = 0;
	uint32_t headlessAndTailLess = 0;
	for (const auto & n : nodes_) {
		if (n->headless()) {
			++headless;
		}
		if (n->tailless()) {
			++tailless;
		}
		if (n->headless() && n->tailless()) {
			++headlessAndTailLess;
		}
	}
	table out(VecStr { "state", "count" });
	out.addRow("tailess", tailless);
	out.addRow("headless", headless);
	out.addRow("headlessAndTailLess", headlessAndTailLess);
	return out;
}
table KmerPathwayGraph::getTailCounts() const{
	std::unordered_map<uint32_t, uint32_t> tailCounts;
	for (const auto & n : nodes_) {
		++tailCounts[n->tailCount()];
	}
	table outTailCount(tailCounts, VecStr{"tails", "counts"});
	outTailCount.sortTable("counts", true);
	return outTailCount;
}
Json::Value KmerPathwayGraph::genSimpleGraphJson() const {
	Json::Value jOut;
	auto & nodes = jOut["nodes"];
	auto & links = jOut["links"];
	for (const auto & node : nodes_) {
		Json::Value jnode;
		jnode["k"] = node->k_;
		jnode["uid"] = node->uid_;
		jnode["size"] = node->cnt_;
		if(node->headless() && node->tailless()){
			jnode["group"] = 2;
		}else if(node->headless()){
			jnode["group"] = 1;
		}else if(node->tailless()){
			jnode["group"] = 3;
		}else{
			jnode["group"] = njh::json::toJson(3 + node->tailCount());
		}
		nodes.append(jnode);
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				Json::Value jlink;
				jlink["source"] = tail->head_.lock()->uid_;
				jlink["target"] = tail->tail_.lock()->uid_;
				jlink["value"] = tail->cnt_;
				links.append(jlink);
			}
		}
	}
	return jOut;
}
VecStr KmerPathwayGraph::writeNodesHeader(){
	return VecStr{"uid", "nodekSize", "kLen", "count", "readCount", "headEdgesCount", "tailEdgesCount"};
}
VecStr KmerPathwayGraph::writeEdgesHeader(){
	return VecStr{"uid", "headNode", "tailNode", "connectorCount", "readCount"};
}
void KmerPathwayGraph::writeEdges(std::ostream & out) const{
	std::set<std::string> edgesUids;
	for(const auto & n : nodes_){
		for(const auto & head : n->headEdges_){
			if(!head->on_){
				continue;
			}
			std::string uid = head->head_.lock()->uid_ + "_" + head->tail_.lock()->uid_;
			if(!njh::in(uid, edgesUids)){
				edgesUids.insert(uid);
				out << uid
						<< "\t" << head->head_.lock()->uid_
						<< "\t" << head->tail_.lock()->uid_
						<< "\t" << head->cnt_
						<< "\t" << head->inReadNamesIdx_.size()
						<< std::endl;
			}
		}
		for(const auto & tail : n->tailEdges_){
			if(!tail->on_){
				continue;
			}
			std::string uid = tail->head_.lock()->uid_ + "_" + tail->tail_.lock()->uid_;
			if(!njh::in(uid, edgesUids)){
				edgesUids.insert(uid);
				out << uid
						<< "\t" << tail->head_.lock()->uid_
						<< "\t" << tail->tail_.lock()->uid_
						<< "\t" << tail->cnt_
						<< "\t" << tail->inReadNamesIdx_.size()
						<< std::endl;
			}
		}
	}
}

VecStr KmerPathwayGraph::writeEdgesWithNamesHeader(){
	return VecStr{"uid", "headNode", "tailNode", "connectorCount", "readCount", "readsNames"};
}

void KmerPathwayGraph::writeEdgesWithNames(std::ostream & out) const{
	std::set<std::string> edgesUids;
	for(const auto & n : nodes_){
		for(const auto & head : n->headEdges_){
			if(!head->on_){
				continue;
			}
			std::string uid = head->head_.lock()->uid_ + "_" + head->tail_.lock()->uid_;
			if(!njh::in(uid, edgesUids)){
				edgesUids.insert(uid);
				std::set<uint32_t> readNamesIdx(head->inReadNamesIdx_.begin(), head->inReadNamesIdx_.end());
				VecStr readNames;
				for(const auto pos : readNamesIdx){
					readNames.emplace_back(readNames_[pos]);
				}
				out << uid
						<< "\t" << head->head_.lock()->uid_
						<< "\t" << head->tail_.lock()->uid_
						<< "\t" << head->cnt_
						<< "\t" << head->inReadNamesIdx_.size()
						<< "\t" << njh::conToStr(readNamesIdx, ",")
						<< "\t" << njh::conToStr(readNames, ",")
						<< std::endl;
			}
		}
		for(const auto & tail : n->tailEdges_){
			if(!tail->on_){
				continue;
			}
			std::string uid = tail->head_.lock()->uid_ + "_" + tail->tail_.lock()->uid_;
			if(!njh::in(uid, edgesUids)){
				edgesUids.insert(uid);
				std::set<uint32_t> readNamesIdx(tail->inReadNamesIdx_.begin(), tail->inReadNamesIdx_.end());
				VecStr readNames;
				for(const auto pos : readNamesIdx){
					readNames.emplace_back(readNames_[pos]);
				}
				out << uid
						<< "\t" << tail->head_.lock()->uid_
						<< "\t" << tail->tail_.lock()->uid_
						<< "\t" << tail->cnt_
						<< "\t" << tail->inReadNamesIdx_.size()
						<< "\t" << njh::conToStr(readNamesIdx, ",")
						<< "\t" << njh::conToStr(readNames, ",")
						<< std::endl;
			}
		}
	}
}



void KmerPathwayGraph::writeNodes(std::ostream & out) const{
	for(const auto & n : nodes_){
		out << n->uid_
				<< "\t" << n->k_.size()
				<< "\t" << klen_
				<< "\t" << n->cnt_
				<< "\t" << n->inReadNamesIdx_.size()
				<< "\t" << n->headCount()
				<< "\t" << n->tailCount()
				<< std::endl;
	}
}
void KmerPathwayGraph::writeOutNodesInGroups(const bfs::path & outDirname, bool overWrite) const{
	std::map<std::string, std::vector<std::shared_ptr<node>>> groupedNodes;
	for(const auto & n : nodes_){
		VecStr connectorNames;
		for(const auto & h : n->headEdges_){
			if(!h->on_){
				continue;
			}
			connectorNames.emplace_back(h->head_.lock()->uid_);
		}
		for(const auto & t : n->tailEdges_){
			if(!t->on_){
				continue;
			}
			connectorNames.emplace_back(t->tail_.lock()->uid_);
		}
		njh::sort(connectorNames);
		if(!connectorNames.empty()){
			groupedNodes[njh::conToStr(connectorNames, ",")].emplace_back(n);
		}else{
			//leave out the headless and tailless nodes
			//groupedNodes[n->uid_].emplace_back(n);
		}
	}
	uint32_t groupId = 0;
	for(const auto & group : groupedNodes){
		OutOptions groupOpts(njh::files::make_path(outDirname, njh::leftPadNumStr<uint32_t>(groupId, groupedNodes.size()) + ".fasta"));
		groupOpts.overWriteFile_ = overWrite;
		OutputStream groupOut(groupOpts);
		for(const auto & n : group.second){
			seqInfo nSeq(n->uid_, n->k_);
			nSeq.outPutSeq(groupOut);
		}
		++groupId;
	}
}

void KmerPathwayGraph::writeRectangleWithEstimatedCovDotByGroup(const OutOptions & prefixPath,
		KmerPathwayGraph & estimatingCovGraph, uint32_t groupSizeCutOff,
		bool noLabels) const {


	double heightNormalizer = 250; // 250 bases will equate 1 inch
	double penWidthNormalizer = 50; // 50 will equate 1 inch

	std::unordered_map<uint32_t, uint32_t> groupCounts;
	for(const auto & n : nodes_){
		if(n->on_){
			++groupCounts[n->group_];
		}
	}

	std::vector<uint32_t> groupsToWrite;
	for(const auto & groupCount : groupCounts){
		if(groupCount.second >= groupSizeCutOff){
			groupsToWrite.emplace_back(groupCount.first);
		}
	}
	auto maxGroup  = vectorMaximum(groupsToWrite);
	for(const auto & group : groupsToWrite){
		VecStr colors = {"#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"};
		VecStr moreColors = { "#ff358f", "#07c652", "#a2009b", "#467f00", "#cb73fd",
				"#f4be45", "#0157d8", "#ff8d3b", "#0056a1", "#dd0d3a", "#01d5e4",
				"#b10049", "#7cda97", "#ff76e0", "#018a5a", "#ff87b8", "#4a5b00",
				"#664092", "#8f7400", "#02aee5", "#9e3500", "#8bd5b8", "#8a306b",
				"#e4b7f0", "#ff97a2" };

		OutOptions groupOpts(bfs::path(njh::pasteAsStr(prefixPath.outFilename_.string(), "_", njh::leftPadNumStr(group, maxGroup), ".dot")));
		groupOpts.transferOverwriteOpts(prefixPath);
		OutputStream out(groupOpts);
		out << "digraph graphname {" << std::endl;
		out << "\t" << "node [fixedsize=true, shape=rect]" << std::endl;
		std::vector<double> estBaseCoverages;
		std::unordered_map<std::string, double> estimatedCoverages;
		for(const auto & node : nodes_){
			if(node->on_ && group == node->group_ ){
				MetaDataInName meta;
				meta.addMeta("headless", node->headless());
				meta.addMeta("tailless", node->tailless());
				auto estCov = CoverageEstimator::estimateCov(node->k_, meta, estimatingCovGraph);
				estBaseCoverages.emplace_back(estCov.minCov_.avgCov_);
				estimatedCoverages[node->uid_] = estCov.minCov_.avgCov_;
			}
		}
		auto maxCov = *std::max_element(estBaseCoverages.begin(), estBaseCoverages.end());


		uint32_t maxColorGroup = 0;
		for (const auto &node : nodes_) {
			if (!node->on_ || group != node->group_) {
				continue;
			}
			uint32_t colorGroup = 0;
			if (node->headless() && node->tailless()) {
				colorGroup = 1;
			} else if (node->headless()) {
				colorGroup = 0;
			} else if (node->tailless()) {
				colorGroup = 2;
			} else {
				colorGroup = 2 + node->tailCount();
			}
			if(colorGroup > maxColorGroup){
				maxColorGroup = colorGroup;
			}
		}
		if(maxColorGroup >= colors.size() + moreColors.size()  ){
			auto hColors = njh::heatmapColors(maxColorGroup + 1);
			njh::reverse(hColors);
			moreColors.clear();
			for(const auto & hColor : hColors){
				moreColors.emplace_back(hColor.getHexStr());
			}
			out << "\t"
					<< "graph [ bgcolor=black, resolution=128, fontname=Arial, fontcolor=white,  fontsize=12 ]; "
					<< std::endl;
			out << "\t" << "node [ fontname=Arial, fontcolor=white, fontsize=11];"
					<< std::endl;
			out << "\t" << "edge [ fontname=Helvetica, fontcolor=white, fontsize=10 ];"
					<< std::endl;
		}
		for(const auto & node : nodes_){
			if(!node->on_ || group != node->group_){
				continue;
			}
			uint32_t colorGroup = 0;
			if(node->headless() && node->tailless()){
				colorGroup = 1;
			}else if(node->headless()){
				colorGroup = 0;
			}else if(node->tailless()){
				colorGroup = 2;
			}else{
				colorGroup = 2 + node->tailCount();
			}
			std::string nodeColor = "";
			if(colorGroup >= colors.size()){
				if(colorGroup >= (colors.size() + moreColors.size())){
					std::stringstream ss;
					ss << "Error in " << __PRETTY_FUNCTION__ << "\n";
					ss << "Oh no, not enough colors!" << "\n";
					ss << "Need " << colorGroup << " but only have " << colors.size() + moreColors.size() << "\n";
					throw std::runtime_error{ss.str()};
				}else{
					nodeColor = moreColors.at(colorGroup - colors.size());
				}
			}else{
				nodeColor = colors.at(colorGroup);
			}
			double nheight = node->k_.length()/heightNormalizer;
			//double approxPerBaseCoverage = static_cast<double>(node->cnt_)/(node->k_.length() - klen_ + 1);
			double estCov = estimatedCoverages[node->uid_];
			double nwidth = (estCov / maxCov) * 5;
			out << "\t" << node->uid_ << "[fixedsize=true,shape=rect,width=" << nwidth
					<< ",height=" << nheight << ",style=filled,fillcolor=\"" << nodeColor
					<< "\", label=\""
					<< (noLabels ?
							"" :
							njh::pasteAsStr("Len=", node->k_.size(), ";\nCov=", roundDecPlaces(estCov, 2), "\n"))
	//				<< (noLabels ?
	//						"" :
	//						njh::pasteAsStr("Len=", node->k_.size(), ";\nCnt=", node->cnt_,
	//								";\nCov=", roundDecPlaces(approxPerBaseCoverage, 2), "\n"))
					<< "\"]" << std::endl;
			;
		}
		for (const auto & node : nodes_){
			if(!node->on_ || group != node->group_){
				continue;
			}
			for(const auto & tail : node->tailEdges_){
				if(tail->on_){
					out << "\t" << tail->head_.lock()->uid_ << " -> " << tail->tail_.lock()->uid_ << "[penwidth=" << tail->cnt_/penWidthNormalizer << ", label=\"" << (noLabels ? "" : estd::to_string(tail->cnt_)) << "\""<< "]"   << std::endl;
				}
			}
		}
		out << "}" << std::endl;
	}
}


void KmerPathwayGraph::writeRectangleWithEstimatedCovDot(std::ostream & out,
		KmerPathwayGraph & estimatingCovGraph,
		bool noLabels) const {
	VecStr colors = {"#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"};
	VecStr moreColors = { "#ff358f", "#07c652", "#a2009b", "#467f00", "#cb73fd",
			"#f4be45", "#0157d8", "#ff8d3b", "#0056a1", "#dd0d3a", "#01d5e4",
			"#b10049", "#7cda97", "#ff76e0", "#018a5a", "#ff87b8", "#4a5b00",
			"#664092", "#8f7400", "#02aee5", "#9e3500", "#8bd5b8", "#8a306b",
			"#e4b7f0", "#ff97a2" };
	out << "digraph graphname {" << std::endl;
	out << "\t" << "node [fixedsize=true, shape=rect]" << std::endl;
	double heightNormalizer = 250; // 250 bases will equate 1 inch
//	double widthNormalizer = 25; // per base coverage of 25 will equate 1 inch
	//double penWidthNormalizer = 50; // 50 will equate 1 inch
	double maxPenWidth = 12.5; // 50 will equate 1 inch
	std::unordered_map<std::string, std::string> nodeLabels;
	uint64_t nodePos = 0;
	for(const auto & node : nodes_){
		if(node->on_){
			if(node->uid_.size() > 16000){
				nodeLabels[node->uid_] = njh::pasteAsStr("node", njh::leftPadNumStr(nodePos, static_cast<uint64_t>(nodes_.size())));
			} else {
				nodeLabels[node->uid_] = node->uid_;
			}
		}
		++nodePos;
	}

//	std::vector<double> avgBaseCoverages;
//	for(const auto & node : nodes_){
//		std::vector<uint32_t> allCounts;
//		for(const auto & kpos : iter::range<uint32_t>(0,node->k_.length() - klen_ + 1)){
//			allCounts.emplace_back(estimatedKCounts[node->k_.substr(kpos, klen_)]);
//		}
//		auto minCov = vectorMinimum(allCounts);
//		avgBaseCoverages.emplace_back(0 == minCov ? 1 : minCov);
//	}
//	auto maxCov = *std::max_element(avgBaseCoverages.begin(), avgBaseCoverages.end());
	std::vector<double> estBaseCoverages;
	std::unordered_map<std::string, double> estimatedCoverages;
	for(const auto & node : nodes_){
		if(node->on_){
			MetaDataInName meta;
			meta.addMeta("headless", node->headless());
			meta.addMeta("tailless", node->tailless());
			auto estCov = CoverageEstimator::estimateCov(node->k_, meta, estimatingCovGraph);
			estBaseCoverages.emplace_back(estCov.minCov_.avgCov_);
			estimatedCoverages[node->uid_] = estCov.minCov_.avgCov_;
		}
	}
	double maxCov = 0.0;
	if(!estBaseCoverages.empty()){
		maxCov = *std::max_element(estBaseCoverages.begin(), estBaseCoverages.end());
	}
	uint32_t maxColorGroup = 0;
	for (const auto &node : nodes_) {
		uint32_t colorGroup = 0;
		if (node->headless() && node->tailless()) {
			colorGroup = 1;
		} else if (node->headless()) {
			colorGroup = 0;
		} else if (node->tailless()) {
			colorGroup = 2;
		} else {
			colorGroup = 2 + node->tailCount();
		}
		if(colorGroup > maxColorGroup){
			maxColorGroup = colorGroup;
		}
	}
	if(maxColorGroup >= colors.size() + moreColors.size() ){
		auto hColors = njh::heatmapColors(maxColorGroup + 1);
		njh::reverse(hColors);
		moreColors.clear();
		for(const auto & hColor : hColors){
			moreColors.emplace_back(hColor.getHexStr());
		}
		out << "\t"
				<< "graph [ bgcolor=black, resolution=128, fontname=Arial, fontcolor=white,  fontsize=12 ]; "
				<< std::endl;
		out << "\t" << "node [ fontname=Arial, fontcolor=white, fontsize=11];"
				<< std::endl;
		out << "\t" << "edge [ fontname=Helvetica, fontcolor=white, fontsize=10 ];"
				<< std::endl;
	}
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		uint32_t colorGroup = 0;
		if(node->headless() && node->tailless()){
			colorGroup = 1;
		}else if(node->headless()){
			colorGroup = 0;
		}else if(node->tailless()){
			colorGroup = 2;
		}else{
			colorGroup = 2 + node->tailCount();
		}
		std::string nodeColor = "";
		if(colorGroup >= colors.size()){
			if(colorGroup >= (colors.size() + moreColors.size())){
				std::stringstream ss;
				ss << "Error in " << __PRETTY_FUNCTION__ << "\n";
				ss << "Oh no, not enough colors!" << "\n";
				ss << "Need " << colorGroup << " but only have " << colors.size() + moreColors.size() << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				nodeColor = moreColors.at(colorGroup - colors.size());
			}
		}else{
			nodeColor = colors.at(colorGroup);
		}
		double nheight = node->k_.length()/heightNormalizer;
		//double approxPerBaseCoverage = static_cast<double>(node->cnt_)/(node->k_.length() - klen_ + 1);
		double estCov = estimatedCoverages[node->uid_];
		double nwidth = (estCov / maxCov) * 5;
		std::unordered_set<std::string> finalNodeReadCounts;
		for(const auto & nameIdx : node->inReadNamesIdx_){
			std::string rName = readNames_[nameIdx];
			if(std::string::npos != rName.find("-multipleOccurenceKmer-")){
				rName = rName.substr(0, rName.find("-multipleOccurenceKmer-"));
			}
			finalNodeReadCounts.emplace(njh::replaceString(rName, "_mate", ""));
		}

		uint32_t readCnt = finalNodeReadCounts.size();
		out << "\t" << nodeLabels[node->uid_] << "[fixedsize=true,shape=rect,width=" << nwidth
				<< ",height=" << nheight << ",style=filled,fillcolor=\"" << nodeColor
				<< "\", label=\""
				<< (noLabels ?
						"" :
						njh::pasteAsStr("Len=", node->k_.size(), ";\nCov=", roundDecPlaces(estCov, 2), "\n", "RCnt=", readCnt))
//				<< (noLabels ?
//						"" :
//						njh::pasteAsStr("Len=", node->k_.size(), ";\nCnt=", node->cnt_,
//								";\nCov=", roundDecPlaces(approxPerBaseCoverage, 2), "\n"))
				<< "\"]" << std::endl;
		;
	}

	uint32_t maxTailCnt = 0;
	for (const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				if(tail->cnt_ > maxTailCnt){
					maxTailCnt = tail->cnt_;
				}
			}
		}
	}

	scale<double> edgePenWidthScale(std::make_pair(0,maxTailCnt), std::make_pair(0,maxPenWidth));
	for (const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				out << "\t" << nodeLabels[tail->head_.lock()->uid_] << " -> " << nodeLabels[tail->tail_.lock()->uid_] << "[penwidth=" <<  edgePenWidthScale.get(tail->cnt_) << ", label=\"" << (noLabels ? "" : estd::to_string(tail->cnt_)) << "\""<< "]"   << std::endl;
			}
		}
	}
	out << "}" << std::endl;
	for(const auto & nodeLabel : nodeLabels){
		if(nodeLabel.first.size() > 16000){
			out << "#" << nodeLabel.second << "=" << nodeLabel.first << std::endl;
		}
	}
}


void KmerPathwayGraph::writeRectangleDot(std::ostream & out, bool noLabels) const{
	VecStr colors = {"#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"};
	VecStr moreColors = { "#ff358f", "#07c652", "#a2009b", "#467f00", "#cb73fd",
			"#f4be45", "#0157d8", "#ff8d3b", "#0056a1", "#dd0d3a", "#01d5e4",
			"#b10049", "#7cda97", "#ff76e0", "#018a5a", "#ff87b8", "#4a5b00",
			"#664092", "#8f7400", "#02aee5", "#9e3500", "#8bd5b8", "#8a306b",
			"#e4b7f0", "#ff97a2" };
	out << "digraph graphname {" << std::endl;
	out << "\t" << "node [fixedsize=true, shape=rect]" << std::endl;
	double heightNormalizer = 250; // 250 bases will equate 1 inch
//	double widthNormalizer = 25; // per base coverage of 25 will equate 1 inch
	double penWidthNormalizer = 50; // 50 will equate 1 inch
	std::vector<double> avgBaseCoverages;
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		double approxPerBaseCoverage = node->cnt_/(node->k_.length() - klen_ + 1);
		avgBaseCoverages.emplace_back(approxPerBaseCoverage);
	}
	auto maxCov = *std::max_element(avgBaseCoverages.begin(), avgBaseCoverages.end());
	uint32_t maxColorGroup = 0;
	for (const auto &node : nodes_) {
		uint32_t colorGroup = 0;
		if (node->headless() && node->tailless()) {
			colorGroup = 1;
		} else if (node->headless()) {
			colorGroup = 0;
		} else if (node->tailless()) {
			colorGroup = 2;
		} else {
			colorGroup = 2 + node->tailCount();
		}
		if(colorGroup > maxColorGroup){
			maxColorGroup = colorGroup;
		}
	}
	if(maxColorGroup >= colors.size() + moreColors.size() ){
		auto hColors = njh::heatmapColors(maxColorGroup + 1);
		njh::reverse(hColors);
		moreColors.clear();
		for(const auto & hColor : hColors){
			moreColors.emplace_back(hColor.getHexStr());
		}
		out << "\t"
				<< "graph [ bgcolor=black, resolution=128, fontname=Arial, fontcolor=white,  fontsize=12 ]; "
				<< std::endl;
		out << "\t" << "node [ fontname=Arial, fontcolor=white, fontsize=11];"
				<< std::endl;
		out << "\t" << "edge [ fontname=Helvetica, fontcolor=white, fontsize=10 ];"
				<< std::endl;
	}
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		uint32_t colorGroup = 0;
		if(node->headless() && node->tailless()){
			colorGroup = 1;
		}else if(node->headless()){
			colorGroup = 0;
		}else if(node->tailless()){
			colorGroup = 2;
		}else{
			colorGroup = 2 + node->tailCount();
		}
		std::string nodeColor = "";
		if(colorGroup >= colors.size()){
			if(colorGroup >= (colors.size() + moreColors.size())){
				std::stringstream ss;
				ss << "Error in " << __PRETTY_FUNCTION__ << "\n";
				ss << "Oh no, not enough colors!" << "\n";
				ss << "Need " << colorGroup << " but only have " << colors.size() + moreColors.size() << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				nodeColor = moreColors.at(colorGroup - colors.size());
			}
		}else{
			nodeColor = colors.at(colorGroup);
		}
		double nheight = node->k_.length()/heightNormalizer;
		double approxPerBaseCoverage = static_cast<double>(node->cnt_)/(node->k_.length() - klen_ + 1);
		double nwidth = (approxPerBaseCoverage / maxCov) * 5;
		out << "\t" << node->uid_ << "[fixedsize=true,shape=rect,width=" << nwidth
				<< ",height=" << nheight << ",style=filled,fillcolor=\"" << nodeColor
				<< "\", label=\""
				<< (noLabels ?
						"" :
						njh::pasteAsStr("Len=", node->k_.size(), ";\nCov=", roundDecPlaces(approxPerBaseCoverage, 2), "\n"))
//				<< (noLabels ?
//						"" :
//						njh::pasteAsStr("Len=", node->k_.size(), ";\nCnt=", node->cnt_,
//								";\nCov=", roundDecPlaces(approxPerBaseCoverage, 2), "\n"))
				<< "\"]" << std::endl;
		;
	}
	for (const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				out << "\t" << tail->head_.lock()->uid_ << " -> " << tail->tail_.lock()->uid_ << "[penwidth=" << tail->cnt_/penWidthNormalizer << ", label=\"" << (noLabels ? "" : estd::to_string(tail->cnt_)) << "\""<< "]"   << std::endl;
			}
		}
	}
	out << "}" << std::endl;
}

void KmerPathwayGraph::writeDot(std::ostream & out) const{
	VecStr colors = {"#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"};
	VecStr moreColors = { "#ff358f", "#07c652", "#a2009b", "#467f00", "#cb73fd",
			"#f4be45", "#0157d8", "#ff8d3b", "#0056a1", "#dd0d3a", "#01d5e4",
			"#b10049", "#7cda97", "#ff76e0", "#018a5a", "#ff87b8", "#4a5b00",
			"#664092", "#8f7400", "#02aee5", "#9e3500", "#8bd5b8", "#8a306b",
			"#e4b7f0", "#ff97a2" };
	out << "digraph graphname {" << std::endl;
	out << "\t" << "node [fixedsize=true regular=true shape=ellipse]" << std::endl;
	uint32_t maxColorGroup = 0;
	for (const auto &node : nodes_) {
		uint32_t colorGroup = 0;
		if (node->headless() && node->tailless()) {
			colorGroup = 1;
		} else if (node->headless()) {
			colorGroup = 0;
		} else if (node->tailless()) {
			colorGroup = 2;
		} else {
			colorGroup = 2 + node->tailCount();
		}
		if(colorGroup > maxColorGroup){
			maxColorGroup = colorGroup;
		}
	}
	if(maxColorGroup >= colors.size() + moreColors.size() ){
		auto hColors = njh::heatmapColors(maxColorGroup + 1);
		njh::reverse(hColors);
		moreColors.clear();
		for(const auto & hColor : hColors){
			moreColors.emplace_back(hColor.getHexStr());
		}
		out << "\t"
				<< "graph [ bgcolor=black, resolution=128, fontname=Arial, fontcolor=white,  fontsize=12 ]; "
				<< std::endl;
		out << "\t" << "node [ fontname=Arial, fontcolor=white, fontsize=11];"
				<< std::endl;
		out << "\t" << "edge [ fontname=Helvetica, fontcolor=white, fontsize=10 ];"
				<< std::endl;
	}
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		uint32_t colorGroup = 0;
		if(node->headless() && node->tailless()){
			colorGroup = 1;
		}else if(node->headless()){
			colorGroup = 0;
		}else if(node->tailless()){
			colorGroup = 2;
		}else{
			colorGroup = 2 + node->tailCount();
		}
		//double diameter = 2 * std::sqrt(node->inReadNamesIdx_.size()/M_PI);
		double diameter = 2 * std::sqrt(node->cnt_/M_PI);
		std::string nodeColor = "";
		if(colorGroup >= colors.size()){
			if(colorGroup >= (colors.size() + moreColors.size())){
				std::stringstream ss;
				ss << "Error in " << __PRETTY_FUNCTION__ << "\n";
				ss << "Oh no, not enough colors!" << "\n";
				ss << "Need " << colorGroup << " but only have " << colors.size() + moreColors.size() << "\n";
				throw std::runtime_error{ss.str()};
			}else{
				nodeColor = moreColors.at(colorGroup - colors.size());
			}
		}else{
			nodeColor = colors.at(colorGroup);
		}
		//out << "\t" << node->uid_  <<"[fixedsize=true,shape=circle,width=" << diameter << ",style=filled,fillcolor=\"" << nodeColor << "\", label=\"" << node->uid_ << "-" << node->inReadNamesIdx_.size()<< "\"]"<< std::endl;
		out << "\t" << node->uid_  <<"[fixedsize=true,shape=circle,width=" << diameter << ",style=filled,fillcolor=\"" << nodeColor << "\", label=\"" << node->uid_ << "-" << node->cnt_ << "\"]"<< std::endl;

	}
	for(const auto & node : nodes_){
		if(!node->on_){
			continue;
		}
		for(const auto & tail : node->tailEdges_){
			if(tail->on_){
				out << "\t" << tail->head_.lock()->uid_ << " -> " << tail->tail_.lock()->uid_ << "[penwidth=" << tail->cnt_ << ", label=\"" << tail->cnt_ << "\""<< "]"   << std::endl;
			}
		}
	}
	out << "}" << std::endl;
}
void KmerPathwayGraph::sortNodesBySize() {
	sortNodes([](const std::shared_ptr<KmerPathwayGraph::node> & node1,
			         const std::shared_ptr<KmerPathwayGraph::node> & node2)->bool {
		return node1->k_.size() > node2->k_.size();
	});
//	sortNodes([](const KmerPathwayGraph::node & node1,
//			         const KmerPathwayGraph::node & node2)->bool {
//		return node1.k_.size() > node2.k_.size();
//	});
}
//void KmerPathwayGraph::sortNodes(std::function<bool(const node&,const node &)> func){
//	std::sort(nodes_.begin(), nodes_.end(), func);
//
//	resetNodePositions();
//}
void KmerPathwayGraph::sortNodes(
		std::function<
				bool(const std::shared_ptr<KmerPathwayGraph::node>&,
						 const std::shared_ptr<KmerPathwayGraph::node>&)> func) {
	std::sort(nodes_.begin(), nodes_.end(), func);
	resetNodePositions();
}

void KmerPathwayGraph::increaseKCounts(const std::string & seq) {
	if(seq.size() > klen_){
		for (auto pos : iter::range(seq.size() - klen_ + 1)) {
			kCounts_[seq.substr(pos, klen_)] += 1;
		}
	}
}


void KmerPathwayGraph::threadThroughSequenceHelperMate(const std::string & firstKmer,
		const std::shared_ptr<node> & headNode,
		const std::string & nextKmer,
		const std::shared_ptr<node> & tailNode,
		std::unordered_map<std::string, uint32_t> & internalCountFirstMate,
		const std::string & threadingSeqNameFirstMate,
		const uint32_t threadingSeqNameFirstMateIdx,
		std::unordered_map<std::string, uint32_t> & internalCountSecondMate,
		const std::string & threadingSeqNameSecondMate,
		const uint32_t threadingSeqNameMateSecondIdx){
	uint32_t firstKmerCountInFirstMate = internalCountFirstMate[firstKmer];
	uint32_t nextKmerCountInFirstMate = internalCountFirstMate[nextKmer];
	std::string threadingSeqNameToUse = threadingSeqNameFirstMate;
	uint32_t threadingSeqNameToUseIdx = threadingSeqNameFirstMateIdx;

	if(0 == firstKmerCountInFirstMate && 0 == nextKmerCountInFirstMate){
	} else {
		threadingSeqNameToUse = threadingSeqNameSecondMate;
		threadingSeqNameToUseIdx = threadingSeqNameMateSecondIdx;
	}
	uint32_t firstKmerCountInSecondMate = internalCountSecondMate[firstKmer];
	uint32_t nextKmerCountInSecondMate = internalCountSecondMate[nextKmer];

	if(0 == firstKmerCountInSecondMate && 0 == nextKmerCountInSecondMate){
		headNode->inReadNamesIdx_.emplace(threadingSeqNameToUseIdx);
		tailNode->inReadNamesIdx_.emplace(threadingSeqNameToUseIdx);
		bool foundEdge = false;
		for(const auto & tail : headNode->tailEdges_){
			auto tailNode = tail->tail_.lock();
			if(nextKmer == tailNode->k_){
				foundEdge = true;
				tail->cnt_ += 1;
				tail->inReadNamesIdx_.emplace(threadingSeqNameToUseIdx);
				break;
			}
		}
		if(!foundEdge){
			std::shared_ptr<edge> e = std::make_shared<edge>(headNode, tailNode, 1, std::unordered_set<uint32_t>{threadingSeqNameToUseIdx});
			headNode->addTail(e);
			tailNode->addHead(e);
		}
	} else {
		std::string maxCountStr = njh::pasteAsStr(threadingSeqNameToUse, "-multipleOccurenceKmer-", std::max(firstKmerCountInSecondMate, nextKmerCountInSecondMate));
		uint32_t maxCountStrIdx = readNames_.size();
		readNames_.emplace_back(maxCountStr);
		headNode->inReadNamesIdx_.emplace(maxCountStrIdx);
		tailNode->inReadNamesIdx_.emplace(maxCountStrIdx);
		bool foundEdge = false;
		for(const auto & tail : headNode->tailEdges_){
			auto tailNode = tail->tail_.lock();
			if(nextKmer == tailNode->k_){
				foundEdge = true;
				tail->cnt_ += 1;
				tail->inReadNamesIdx_.emplace(maxCountStrIdx);
				break;
			}
		}
		if(!foundEdge){
			std::shared_ptr<edge> e = std::make_shared<edge>(headNode, tailNode, 1, std::unordered_set<uint32_t>{maxCountStrIdx});
			headNode->addTail(e);
			tailNode->addHead(e);
		}
	}
	++internalCountSecondMate[firstKmer];
}


void KmerPathwayGraph::threadThroughSequenceHelper(
		const std::string & firstKmer,
		const std::shared_ptr<node> & headNode,
		const std::string & nextKmer,
		const std::shared_ptr<node> & tailNode,
		std::unordered_map<std::string, uint32_t> & internalCount,
		const std::string & threadingSeqName,
		const uint32_t threadingSeqNameIdx
		){

	uint32_t firstKmerCount = internalCount[firstKmer];
	uint32_t nextKmerCount = internalCount[nextKmer];
	if(0 == firstKmerCount && 0 == nextKmerCount){
		headNode->inReadNamesIdx_.emplace(threadingSeqNameIdx);
		tailNode->inReadNamesIdx_.emplace(threadingSeqNameIdx);
		bool foundEdge = false;
		for(const auto & tail : headNode->tailEdges_){
			auto tailNode = tail->tail_.lock();
			if(nextKmer == tailNode->k_){
				foundEdge = true;
				tail->cnt_ += 1;
				tail->inReadNamesIdx_.emplace(threadingSeqNameIdx);
				break;
			}
		}
		if(!foundEdge){
			std::shared_ptr<edge> e = std::make_shared<edge>(headNode, tailNode, 1, std::unordered_set<uint32_t>{threadingSeqNameIdx});
			headNode->addTail(e);
			tailNode->addHead(e);
		}
	}else{
		std::string maxCountStr = njh::pasteAsStr(threadingSeqName, "-multipleOccurenceKmer-", std::max(firstKmerCount, nextKmerCount));
		uint32_t maxCountStrIdx = readNames_.size();
		readNames_.emplace_back(maxCountStr);
		headNode->inReadNamesIdx_.emplace(maxCountStrIdx);
		tailNode->inReadNamesIdx_.emplace(maxCountStrIdx);
		bool foundEdge = false;
		for(const auto & tail : headNode->tailEdges_){
			auto tailNode = tail->tail_.lock();
			if(nextKmer == tailNode->k_){
				foundEdge = true;
				tail->cnt_ += 1;
				tail->inReadNamesIdx_.emplace(maxCountStrIdx);
				break;
			}
		}
		if(!foundEdge){
			std::shared_ptr<edge> e = std::make_shared<edge>(headNode, tailNode, 1, std::unordered_set<uint32_t>{maxCountStrIdx});
			headNode->addTail(e);
			tailNode->addHead(e);
		}
	}
	++internalCount[firstKmer];
}

std::unordered_map<std::string, uint32_t> KmerPathwayGraph::threadThroughSequence(
		const seqInfo & seq, const std::string & threadingSeqName) {
	std::unordered_map<std::string, uint32_t> internalCount;
	if(seq.seq_.size() > klen_){
		uint32_t threadingSeqNameIdx = readNames_.size();
		readNames_.emplace_back(threadingSeqName);
		std::string startKmer = seq.seq_.substr(0, klen_);
		auto lastPosIter = nodePositions_.find(startKmer);
		for (const auto pos : iter::range(seq.seq_.size() - klen_)) {
			std::string firstKmer = seq.seq_.substr(pos, klen_);
			std::string nextKmer = seq.seq_.substr(pos + 1, klen_);
			auto firstNodePosition = lastPosIter;
			auto nextNodePosition =  nodePositions_.find(nextKmer);
			//set the last position iterator as the next node iterator
			lastPosIter = nextNodePosition;
			//check if both nodes pass occurrence check
			bool passOccurenceCheck = nodePositions_.end() != firstNodePosition && nodePositions_.end() != nextNodePosition;
			if(passOccurenceCheck){
				threadThroughSequenceHelper(firstKmer,
						nodes_[firstNodePosition->second],
						nextKmer,
						nodes_[nextNodePosition->second],
						internalCount,
						threadingSeqName,
						threadingSeqNameIdx);
			}
		}
	}
	return internalCount;
}

std::unordered_map<std::string, uint32_t> KmerPathwayGraph::threadThroughSequence(const seqInfo & seq){
	return threadThroughSequence(seq, seq.name_);
}




std::unordered_map<std::string, uint32_t> KmerPathwayGraph::threadThroughSequenceMate(const PairedRead & pairedSeq,
		std::unordered_map<std::string, uint32_t> & firstMateCounts,
		const std::string & threadingSeqNameFirstMate,
		uint32_t threadingSeqNameFirstMateIdx){
	std::unordered_map<std::string, uint32_t> secondMateInternalCount;


	if(pairedSeq.mateSeqBase_.seq_.size() > klen_){
		uint32_t threadingSeqNameSecondMateIdx = readNames_.size();
		std::string threadingSeqNameSecondMate = threadingSeqNameFirstMate + "_mate";
		readNames_.emplace_back(threadingSeqNameSecondMate);
		std::string startKmer = pairedSeq.mateSeqBase_.seq_.substr(0, klen_);
		auto lastPosIter = nodePositions_.find(startKmer);
		for (const auto pos : iter::range(pairedSeq.mateSeqBase_.seq_.size() - klen_)) {
			std::string firstKmer = pairedSeq.mateSeqBase_.seq_.substr(pos, klen_);
			std::string nextKmer = pairedSeq.mateSeqBase_.seq_.substr(pos + 1, klen_);
			auto firstNodePosition = lastPosIter;
			auto nextNodePosition =  nodePositions_.find(nextKmer);
			//set the last position iterator as the next node iterator
			lastPosIter = nextNodePosition;
			//check if both nodes pass occurrence check
			bool passOccurenceCheck = nodePositions_.end() != firstNodePosition && nodePositions_.end() != nextNodePosition;
			if(passOccurenceCheck){
				threadThroughSequenceHelperMate(
						firstKmer,
						nodes_[firstNodePosition->second],
						nextKmer,
						nodes_[nextNodePosition->second],
						firstMateCounts,
						threadingSeqNameFirstMate,
						threadingSeqNameFirstMateIdx,
						secondMateInternalCount,
						threadingSeqNameSecondMate,
						threadingSeqNameSecondMateIdx
						);
			}
		}
	}
	return secondMateInternalCount;
}

void KmerPathwayGraph::threadThroughSequence(const PairedRead & pSeq) {
	if(pSeq.seqBase_.seq_.size() <= klen_){
		threadThroughSequence(pSeq.mateSeqBase_);
	}else{
		uint32_t threadingSeqNameFirstMateIdx = readNames_.size();
		std::string threadingSeqNameFirstMate = pSeq.seqBase_.name_;
		auto firstMateCounts = threadThroughSequence(pSeq.seqBase_, pSeq.seqBase_.name_);
		threadThroughSequenceMate(pSeq, firstMateCounts, threadingSeqNameFirstMate, threadingSeqNameFirstMateIdx);
	}
}

std::vector<seqInfo> KmerPathwayGraph::nodesToSeqs(bool addSeqOfSingleHeadAndTailSeqs) const{
	std::vector<seqInfo> ret;
	for(const auto & node : nodes_){
		std::string seq = node->k_;
		if(addSeqOfSingleHeadAndTailSeqs){
			if(1 == node->headCount()){
				std::string head = "";
				for(const auto & edge : node->headEdges_){
					if(edge->on_){
						head = edge->head_.lock()->k_;
						break;
					}
				}
				seq = head.substr(0, head.size() - klen_ + 1) + seq;
			}
			if(1 == node->tailCount()){
				std::string tail = "";
				for(const auto & edge : node->tailEdges_){
					if(edge->on_){
						tail = edge->tail_.lock()->k_;
						break;
					}
				}
				seq.append(tail.substr(klen_ - 1));
			}
		}
		ret.emplace_back(node->uid_, seq);
	}
	return ret;
}
bool KmerPathwayGraph::skipInputSeqForKCount(const std::string & seq, uint32_t kLen) {
	if(seq.size() <= kLen){
		return true;
	}
	/**@todo like the tandem repeat problem add in expanding around large homopolymers*/
//	auto isHomopolymer =
//			[](const std::string & k) {
//				return std::all_of(k.begin(), k.end(),[&k](const char c) {return k.front() == c;});
//			};
	if (std::string::npos != seq.find('N')) {
		return true;
	} else {
//		for (auto pos : iter::range(seq.size() - kLen + 1)) {
//			if (isHomopolymer(seq.substr(pos, kLen))) {
//				return true;
//			}
//		}
	}
	return false;
}

bool KmerPathwayGraph::skipInputSeqForKCountOld(const std::string & seq, uint32_t kLen) {
	if(seq.size() <= kLen){
		return true;
	}
	/**@todo like the tandem repeat problem add in expanding around large homopolymers*/
	auto isHomopolymer =
			[](const std::string & k) {
				return std::all_of(k.begin(), k.end(),[&k](const char c) {return k.front() == c;});
			};
	if (std::string::npos != seq.find('N')) {
		return true;
	} else {
		for (auto pos : iter::range(seq.size() - kLen + 1)) {
			if (isHomopolymer(seq.substr(pos, kLen))) {
				return true;
			}
		}
	}
	return false;
}

void KmerPathwayGraph::resetGroups() const {
	//really should do a non-recursive group set

	//first reset all nodes back to max
	for(auto & n : nodes_){
		n->group_ = std::numeric_limits<uint32_t>::max();
	}
	std::vector<std::shared_ptr<node>> nodesToProcess;
	for(const auto & n : nodes_){
		if(n->headless()){
			nodesToProcess.emplace_back(n);
		}
	}
	std::function<void(const std::shared_ptr<node> &, uint32_t)> spreadGroup = [&spreadGroup](const std::shared_ptr<node> & spnode, uint32_t groupUID){
		//set self
		spnode->group_ = groupUID;
		//spread to heads
		for(const auto & head : spnode->headEdges_){
			if(head->on_){
				if(std::numeric_limits<uint32_t>::max() == head->head_.lock()->group_){
					spreadGroup(head->head_.lock(), groupUID);
				}
			}
		}
		//spread to tails
		for(const auto & tail : spnode->tailEdges_){
			if(tail->on_){
				if(std::numeric_limits<uint32_t>::max() == tail->tail_.lock()->group_){
					spreadGroup(tail->tail_.lock(), groupUID);
				}
			}
		}
	};
	uint32_t groupId = 0;
	for(const auto & n : nodesToProcess){
		//group hasn't been set yet
		if(std::numeric_limits<uint32_t>::max() == n->group_){
			spreadGroup(n, groupId);
			++groupId;
		}
	}
}

void KmerPathwayGraph::resetGroupsLoopAware() const {
	//really should do a non-recursive group set

	//first reset all nodes back to max
	for(auto & n : nodes_){
		n->group_ = std::numeric_limits<uint32_t>::max();
	}
	std::function<void(const std::shared_ptr<node> &, uint32_t)> spreadGroup = [&spreadGroup](const std::shared_ptr<node> & spnode, uint32_t groupUID){
		//set self
		spnode->group_ = groupUID;
		//spread to heads
		for(const auto & head : spnode->headEdges_){
			if(head->on_){
				if(std::numeric_limits<uint32_t>::max() == head->head_.lock()->group_){
					spreadGroup(head->head_.lock(), groupUID);
				}
			}
		}
		//spread to tails
		for(const auto & tail : spnode->tailEdges_){
			if(tail->on_){
				if(std::numeric_limits<uint32_t>::max() == tail->tail_.lock()->group_){
					spreadGroup(tail->tail_.lock(), groupUID);
				}
			}
		}
	};
	uint32_t groupId = 0;
	for(const auto & n : nodes_){
		//group hasn't been set yet
		if(std::numeric_limits<uint32_t>::max() == n->group_){
			spreadGroup(n, groupId);
			++groupId;
		}
	}
}




//KmerPathwayGraph::Pathway::Pathway() {
//}
//
//size_t KmerPathwayGraph::Pathway::size() const{
//	return edges_.size();
//}
//
//KmerPathwayGraph::Pathway::Pathway(const std::shared_ptr<edge> & edge) :
//		edges_(edge), edgeUids_{edge->createUid()} {
//}
//
//bool KmerPathwayGraph::Pathway::hasEdge(const std::shared_ptr<edge> & edge) const{
//	return hasEdge(edge->createUid());
//
//}
//
//bool KmerPathwayGraph::Pathway::hasEdge(const std::string & edgeUid) const{
//	if (njh::in(edgeUid, edgeUids_)) {
//		return true;
//	}
//	return false;
//}
//
//
//bool KmerPathwayGraph::pathwaysPossible() const{
//	//initiate with the headless nodes
//	for(const auto & n : nodes_){
//		if(n->headless() && !n->tailless()){
//			return true;
//		}
//	}
//	return false;
//}
//
//void KmerPathwayGraph::extendPathways(std::shared_ptr<Pathway> & path, std::vector<std::shared_ptr<Pathway>> & pathways) const{
//	auto nextNode = path->edges_.back()->tail_.lock();
//	if(nextNode->tailless()){
//		pathways.emplace_back(path);
//	}else{
//		if(1 == nextNode->tailCount()){
//			std::string nextEdgeUid = nextNode->getFirstOnTailEdge()->createUid();
//			if(path->hasEdge(nextEdgeUid)){
//				std::stringstream ss;
//				ss << __PRETTY_FUNCTION__ << ", error already have edge: " << nextEdgeUid << "\n";
//				throw std::runtime_error{ss.str()};
//			}
//			path->edgeUids_.insert(nextEdgeUid);
//			path->edges_.emplace_back(nextNode->getFirstOnTailEdge());
//			extendPathways(path, pathways);
//		}else{
//			bool extended = false;
//			for(const auto & tail : nextNode->tailEdges_){
//				std::string nextEdgeUid = tail->createUid();
//
//			}
//		}
//	}
//}
//
//std::vector<std::shared_ptr<KmerPathwayGraph::Pathway>> KmerPathwayGraph::getAllPaths() const {
//	//one way forward, multiple ways forward
//
//	std::vector<std::shared_ptr<KmerPathwayGraph::Pathway>> starters;
//	//initiate with the headless nodes
//	for(const auto & n : nodes_){
//		if(n->headless() && !n->tailless()){
//			for(const auto & tail : n->tailEdges_){
//				starters.emplace_back(std::make_shared<Pathway>(tail));
//			}
//		}
//	}
//
//	if(starters.empty()){
//		std::stringstream ss;
//		ss << __PRETTY_FUNCTION__ << ", no starting points for pathways" << "\n";
//		throw std::runtime_error{ss.str()};
//	}
//
//	std::vector<std::shared_ptr<KmerPathwayGraph::Pathway>> ret;
//
//	for(const auto & path : starters){
//		auto nextNode = path->edges_.back()->tail_.lock();
//
//	}
//
//	return ret;
//
//}
}  // namespace njhseq
