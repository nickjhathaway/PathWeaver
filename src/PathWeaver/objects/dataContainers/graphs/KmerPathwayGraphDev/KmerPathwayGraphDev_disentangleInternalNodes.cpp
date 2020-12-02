/*
 * KmerPathwayGraphDev_disentangleInternalNodes.cpp
 *
 *  Created on: Feb 20, 2020
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




#include "KmerPathwayGraphDev.hpp"

namespace njhseq {

template<typename T>
std::vector<T> getIntersection(const std::set<T> & set1, const std::set<T> & set2){
	std::vector<T> ret;
	std::set_intersection(
			set1.begin(), set1.end(),
			set2.begin(), set2.end(), std::back_inserter(ret));
	return ret;
}


bool KmerPathwayGraphDev::disentangleInternalNodes(const disentangleInternalNodesPars & pars) {
	resetNodePositions();
	//reset visit count, visit counts are used to make sure to wait until the next call to dis-entagnle a node if it's been affected during this call
	njh::for_each(nodes_,[](std::shared_ptr<node> & node){
		node->resetVisitCount();
	});

	//
	uint32_t connectorCutOff = occurenceCutOff_;
	if (occurenceCutOff_ >= 6) {
		//connectorCutOff = round(connectorCutOff * 0.50);
		//connectorCutOff = connectorCutOff * 0.50;
		connectorCutOff = std::max<uint32_t>(3, connectorCutOff * 0.34);
		//connectorCutOff = connectorCutOff * 0.5;
		//connectorCutOff = std::max<uint32_t>(2, connectorCutOff * 0.2);
	} else {
		connectorCutOff = std::min<uint32_t>(3, connectorCutOff);
	}
//	if(!conservative){
//		connectorCutOff = 0;
//	}
#if defined(PATHWEAVERSUPERDEBUG)
	{
		OutOptions nodeCheckOpts(bfs::path("nodeCheck.txt"));
		nodeCheckOpts.overWriteFile_ = true;
		OutputStream nodeCheckOut(nodeCheckOpts);
		nodeCheckOut << "Node check: " << std::endl;
		uint32_t nodeCount = 0;
		for(const auto & n : nodes_){
			++nodeCount;
			if("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCT" == n->nameUid_){
				nodeCheckOut << "nodeCount: " << nodeCount << std::endl;
				nodeCheckOut << "n->k_   :"<< n->k_ << std::endl;
				nodeCheckOut << "n->cnt_ :"<< n->cnt_ << std::endl;
			}
		}
		bool foundDupEdges = false;
		for(const auto & n : nodes_){
			if(n->on_){
				//heads
				std::unordered_map<uint32_t, uint32_t> headUidsCounts;
				for(const auto & head : n->headEdges_){
					if(head->on_){
						++headUidsCounts[head->head_.lock()->nodeUid_];
					}
				}
				for(const auto & huid : headUidsCounts){
					if(huid.second > 1){
						foundDupEdges = true;
						std::cout << "Found Dup heads for " << n->nodeUid_ << " " << n->nameUid_ << std::endl;
						for(const auto & huidPrint : headUidsCounts){
							std::cout << "\t" << huidPrint.first << ": " << huidPrint.second << std::endl;
						}
						break;
					}
				}
				//tails
				std::unordered_map<uint32_t, uint32_t> tailUidsCounts;
				for(const auto & tail : n->tailEdges_){
					if(tail->on_){
						++tailUidsCounts[tail->tail_.lock()->nodeUid_];
					}
				}
				for(const auto & tuid : tailUidsCounts){
					if(tuid.second > 1){
						foundDupEdges = true;
						std::cout << "Found Dup Tails for " << n->nodeUid_ << " " << n->nameUid_ << std::endl;
						for(const auto & tuidPrint : tailUidsCounts){
							std::cout << "\t" << tuidPrint.first << ": " << tuidPrint.second << std::endl;
						}
						break;
					}
				}
			}
		}
		if(foundDupEdges){
			exit(1);
		}
	}
#endif

	std::vector<std::shared_ptr<KmerPathwayGraphDev::node>> nodesToProcess;
	//first add in multi-tailed and multi-headed nodes
	for (const auto & n : nodes_) {
		if (n->headCount() > 1 && n->tailCount() > 1) {
			nodesToProcess.push_back(n);
		}
	}
	//then add in multiheaded or multitailed with either one tail or head
	for (const auto & n : nodes_) {
		if ((n->headCount() > 1 && n->tailCount() == 1) || (n->headCount() == 1 && n->tailCount() > 1)) {
			nodesToProcess.push_back(n);
		}
	}

	//sort by the most number bridging so those are processed first
	std::unordered_map<uint32_t, uint32_t> nodeUIDToBridgeCount;
	for(const auto & n : nodesToProcess){
		std::set<uint32_t> headReads;
		std::set<uint32_t> tailReads;
		for(const auto & head : n->headEdges_){
			if(head->on_){
				headReads.insert(head->inReadNamesIdx_.begin(), head->inReadNamesIdx_.end());
			}
		}
		for(const auto & tail : n->tailEdges_){
			if(tail->on_){
				tailReads.insert(tail->inReadNamesIdx_.begin(), tail->inReadNamesIdx_.end());
			}
		}
		auto bridging = getIntersection(headReads, tailReads);
		nodeUIDToBridgeCount[n->nodeUid_] = bridging.size();
	}
	njh::sort(nodesToProcess,[&nodeUIDToBridgeCount](const auto & node1, const auto & node2){
		return nodeUIDToBridgeCount[node1->nodeUid_] > nodeUIDToBridgeCount[node2->nodeUid_];
	});

	/*
	std::set<std::string> connected;
	for (const auto & n : nodes_) {
		if (n->headCount() > 1 && n->tailCount() > 1) {
			if (!njh::in(n->uid_, connected)) {
				nodesToProcess.push_back(n);
				for (const auto & head : n->headEdges_) {
					connected.emplace(head->head_.lock()->uid_);
				}
				for (const auto & tail : n->tailEdges_) {
					connected.emplace(tail->head_.lock()->uid_);
				}
			}
		}
	}*/
#if defined(PATHWEAVERSUPERDEBUG)
	{
		std::cout << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << std::endl;
		std::cout << "There are " << nodesToProcess.size() << " nodes to process" << std::endl;
	}
#endif
	if (nodesToProcess.empty()) {
		return false;
	}



	bool nodesSplit = false;
	uint32_t nodeProcessCount = 0;
	for (const auto & n : nodesToProcess) {
		//bool printInfo = false;
//		if(n->k_.length() == 136){
//			printInfo = true;
//		}
//		if(n->k_ == "TCCAATCAACGAACATAGGGAACATCCAAAAGAATAC"){
//			printInfo = true;
//		}

//		if (n->k_ == "AAAACAGAAAAATTTTGTCATAACTGTGGGTATGGGTTAGGAAGT") {
//			printInfo = true;
//		}
//
//		if (n->k_ == "AAAACAGAAAAATTTTGTCATAACTGTGGGTATGGGTTAGGAAGTGTTGCACCGAATATTGGATTATTGGGGGGACCGGGAATCTATGTATGGAAAATTGCCGCGTTGGCA") {
//			printInfo = true;
//		}

//		if (n->k_ == "GAAAAAAGAAAATGAAGAAAAAACTACAATATATAAAATTATTAAA") {
//			printInfo = true;
//		}
//		if (n->k_ == "AGTAGTTATAGTTTTGGTTATGGTAATTATTTATTTAATTTTACATTATCGACGAAAAAAGAAAATGAAGAAAAAACTACAATATATAAAATTATTA") {
//			printInfo = true;
//		}
//		if (n->k_ == "TATGTGATTTTTTTAGTTATAGTTTTTATCATAAAATAATATACGTATCACAATAAAAAATGAAAA") {
//			printInfo = true;
//		}
//		if (n->k_ == "TTTTATGTGATTTTTTTAGTTATAGTTTTTATCATAAAATAATATAC") {
//			printInfo = true;
//		}
//
//		if (n->k_ == "GATATAAAATAATATACTCCAACTAATTATAACATACATATATAT") {
//			printInfo = true;
//		}
//
//		if (n->k_ == "ATATAAAATAATATACTCCAACTAATTATAACATACATATATATATATTTTATAGGCACATAATAAA") {
//			printInfo = true;
//		}
//
//		if (n->k_ == "AAATTATTAAATCAATAAATATATGGTTTCTTGATATTAAATTCAATTTAAT") {
//			printInfo = true;
//		}
//		if (n->k_ == "ACACAAAAAATACCAACCACAAGGTCATTAAGCGAATGTGAATTATTTTCACCACAAAATTACGACAACGATCCTGAAATGAAA") {
//			printInfo = true;
//		}
//		if (n->k_ == "ATTACAACACGCATATCCAATAGACCACGAAGGTGCCGAACCCGCACCACAAGAACAAAATTTATTTTCAAGCATTGAAATAGTAGAAAGAAGTAATTATATGGGTAATCCATGGACGGAATATATGGCAAAATATGATATTGAAGAAGTTCATGGTTCAGGTATAAGAGTAGATTTAGGAGAAGATGCTGAAGTAGCTGGAACTCAATATAGACTTCCATCAGGGAAATGTCCAGTATTTGGTAAAGGTATAATTATTGAGAATTCAAATACTACTTTTTTAACACCGGTAGCTACGGAAAATCAAGATTTAAAAGATGGAGGTTTTGCTTTTCCTCCAACAAA") {
//			printInfo = true;
//		}


		//
#if defined(PATHWEAVERSUPERDEBUG)
		{
			std::cout << njh::bashCT::red;
			std::cout << "n->k_: " << n->k_ << std::endl;
			std::cout << "n->uid_: " << n->nodeUid_ << " " << n->nameUid_ << std::endl;
			std::cout << "n->visitCount_: " << n->visitCount_ << std::endl;
			std::cout << njh::bashCT::reset;
		}
#endif

		if(0 != n->visitCount_ ){
			//node has been affect already, wait until next call
			continue;
		}
#if defined(PATHWEAVERSUPERDEBUG)
		{
			bool foundDupEdges = false;
			for(const auto & n : nodes_){
				if(n->on_){
					//heads
					std::unordered_map<uint32_t, uint32_t> headUidsCounts;
					for(const auto & head : n->headEdges_){
						if(head->on_){
							++headUidsCounts[head->head_.lock()->nodeUid_];
						}
					}
					for(const auto & huid : headUidsCounts){
						if(huid.second > 1){
							foundDupEdges = true;
							std::cout << "Found Dup heads for " << n->nodeUid_ << " " << n->nameUid_ << std::endl;
							for(const auto & huidPrint : headUidsCounts){
								std::cout << "\t" << huidPrint.first << ": " << huidPrint.second << std::endl;
							}
							break;
						}
					}
					//tails
					std::unordered_map<uint32_t, uint32_t> tailUidsCounts;
					for(const auto & tail : n->tailEdges_){
						if(tail->on_){
							++tailUidsCounts[tail->tail_.lock()->nodeUid_];
						}
					}
					for(const auto & tuid : tailUidsCounts){
						if(tuid.second > 1){
							foundDupEdges = true;
							std::cout << "Found Dup Tails for " << n->nodeUid_ << " " << n->nameUid_ << std::endl;;
							for(const auto & tuidPrint : tailUidsCounts){
								std::cout << "\t" << tuidPrint.first << ": " << tuidPrint.second << std::endl;
							}
							break;
						}
					}
				}
			}

			if(foundDupEdges){
				std::cout << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << std::endl;
				std::cout << "nodeProcessCount: " << nodeProcessCount << std::endl;
				exit(1);
			}
		}
#endif
		++nodeProcessCount;
//		if (n->headCount() > 1 && n->tailCount() > 1) {
		if (!n->headless() && !n->tailless() && (n->headCount() > 1 || n->tailCount() > 1)) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << std::endl;
				std::cout << n->nodeUid_ << " " << n->nameUid_ << std::endl;;
			}
#endif
			//head, read, that read count
			std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t> > incoming;
			std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t> > outgoing;
			std::unordered_set<uint32_t> bridgingReads;
			std::unordered_set<uint32_t> totalEnteringLeaving;
			std::unordered_set<uint32_t> headReads;
			std::unordered_set<uint32_t> tailReads;
			for (const auto & head : n->headEdges_) {
				if(!head->on_){
					continue;
				}
				auto headNode = head->head_.lock();
				headReads.insert(head->inReadNamesIdx_.begin(), head->inReadNamesIdx_.end());
				for (const auto & readName : head->inReadNamesIdx_) {
					++incoming[headNode->nodeUid_][readName];
					//headReads.emplace(readName);
				}
			}
			//check head integrity of head and tail edges
			for (const auto & tail : n->tailEdges_) {
				if(!tail->on_){
					continue;
				}
				auto tailNode = tail->tail_.lock();
				tailReads.insert(tail->inReadNamesIdx_.begin(), tail->inReadNamesIdx_.end());
				for (const auto & readName : tail->inReadNamesIdx_) {
					++outgoing[tailNode->nodeUid_][readName];
					//tailReads.emplace(readName);
					if(njh::in(readName, headReads)){
						bridgingReads.emplace(readName);
					}
				}
			}
			totalEnteringLeaving.insert(headReads.begin(), headReads.end());
			totalEnteringLeaving.insert(tailReads.begin(), tailReads.end());
#if defined(PATHWEAVERSUPERDEBUG)
			{
				table incommingTab(incoming, VecStr { "head", "readNames", "count" });
				incommingTab.addColumn(VecStr { "incomming" }, "direction");
				table outgoingTab(outgoing, VecStr { "head", "readNames", "count" });
				outgoingTab.addColumn(VecStr { "outgoing" }, "direction");
				outgoingTab.rbind(incommingTab, false);
				std::ofstream outOut(n->nameUid_ + "_goings.tab.txt");
				outgoingTab.outPutContents(outOut, "\t");
			}
#endif

			//if (debug_ || true ) {
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << std::endl;
				std::cout << njh::bashCT::cyan;
				std::cout << "\t" << "internalReads       : " << n->inReadNamesIdx_.size() << std::endl;
				std::cout << "\t" << "bridgingReads       : " << bridgingReads.size() << " " << getPercentageString(bridgingReads.size(), n->inReadNamesIdx_.size())<< std::endl;
				std::cout << "\t" << "headReads           : " << headReads.size() << " " << getPercentageString(headReads.size(), n->inReadNamesIdx_.size())<< std::endl;
				std::cout << "\t" << "tailReads           : " << tailReads.size() << " " << getPercentageString(tailReads.size(), n->inReadNamesIdx_.size())<< std::endl;
				std::cout << "\t" << "n->inReadNamesIdx_.size() * bridgingCutOff: " << n->inReadNamesIdx_.size() * bridgingCutOff_ << std::endl;
				std::cout << "\t" << "(n->inReadNamesIdx_.size() * bridgingCutOff) > bridgingReads.size(): " <<  njh::colorBool((n->inReadNamesIdx_.size() * bridgingCutOff_) > bridgingReads.size()) << std::endl;
				//std::cout << "\t" << "(n->inReadNamesIdx_.size() * bridgingCutOff) < bridgingReads.size(): " <<  njh::colorBool((n->inReadNamesIdx_.size() * bridgingCutOff) < bridgingReads.size()) << std::endl;
				std::cout << njh::bashCT::blue;
				std::cout << "\t" << "totalEnteringLeaving: " << totalEnteringLeaving.size() << std::endl;
				std::cout << "\t" << "bridgingReads       : " << bridgingReads.size() << " " << getPercentageString(bridgingReads.size(), totalEnteringLeaving.size())<< std::endl;
				std::cout << "\t" << "headReads           : " << headReads.size() << " " << getPercentageString(headReads.size(), totalEnteringLeaving.size())<< std::endl;
				std::cout << "\t" << "tailReads           : " << tailReads.size() << " " << getPercentageString(tailReads.size(), totalEnteringLeaving.size())<< std::endl;
				std::cout << "\t" << "totalEnteringLeaving.size() * bridgingCutOff: " << totalEnteringLeaving.size() * bridgingCutOff_ << std::endl;
				std::cout << "\t" << "(totalEnteringLeaving.size() * bridgingCutOff) > bridgingReads.size(): " <<  njh::colorBool((totalEnteringLeaving.size() * bridgingCutOff_) > bridgingReads.size()) << std::endl;
				//std::cout << "\t" << "(totalEnteringLeaving.size() * bridgingCutOff) < bridgingReads.size(): " <<  njh::colorBool((totalEnteringLeaving.size() * bridgingCutOff) < bridgingReads.size()) << std::endl;
				std::cout << njh::bashCT::purple;
				std::cout << "\t" << "occurenceCutOff_: " << occurenceCutOff_ << std::endl;
				std::cout << "\t" << "connectorCutOff: " << connectorCutOff << std::endl;
				std::cout << "\t" << "bridgingReads.size() <= occurenceCutOff_: " << njh::colorBool(bridgingReads.size() <= occurenceCutOff_) << std::endl;
				std::cout << "\t" << "bridgingReads.size() <= connectorCutOff: " << njh::colorBool(bridgingReads.size() <= connectorCutOff) << std::endl;

				std::cout << njh::bashCT::red;
				std::cout << "\t" << "(totalEnteringLeaving.size() * bridgingCutOff_) > bridgingReads.size() || bridgingReads.size() <= occurenceCutOff_: " << njh::colorBool((totalEnteringLeaving.size() * bridgingCutOff_) > bridgingReads.size() || bridgingReads.size() <= occurenceCutOff_) << std::endl;
				std::cout << "\t" << "(totalEnteringLeaving.size() * bridgingCutOff_) > bridgingReads.size() || bridgingReads.size() <= connectorCutOff: " << njh::colorBool((totalEnteringLeaving.size() * bridgingCutOff_) > bridgingReads.size() || bridgingReads.size() <= connectorCutOff) << std::endl;

				std::cout << njh::bashCT::reset;

//				if(!((n->inReadNamesIdx_.size() * bridgingCutOff) < bridgingReads.size()) &&
//						((totalEnteringLeaving.size() * bridgingCutOff) < bridgingReads.size()) ){
//					std::stringstream ss;
//					ss << __PRETTY_FUNCTION__ << ", error " << "hold the phone"<< "\n";
//					throw std::runtime_error{ss.str()};
//				}
//				if(!((n->inReadNamesIdx_.size() * bridgingCutOff) < bridgingReads.size()) &&
//						!((totalEnteringLeaving.size() * bridgingCutOff) < bridgingReads.size()) &&
//						bridgingReads.size() > 10){
//					std::stringstream ss;
//					ss << __PRETTY_FUNCTION__ << ", error " << "hold the phone"<< "\n";
//					throw std::runtime_error{ss.str()};
//				}
			}
#endif
			//if less than 10% of reads entering and leaving a node are bridging then skip
//			if((n->inReadNamesIdx_.size() * bridgingCutOff) > bridgingReads.size()){
//				continue;
//			}
			//if less than 10% of reads entering and leaving a node are bridging then skip

			//bool skip = false;
			if((totalEnteringLeaving.size() * bridgingCutOff_) > bridgingReads.size() || bridgingReads.size() <= connectorCutOff){
			//if((totalEnteringLeaving.size() * bridgingCutOff_) > bridgingReads.size() || bridgingReads.size() <= occurenceCutOff_){
				continue;
			}

			std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> headsToTailsCounts;
#if defined(PATHWEAVERSUPERDEBUG)
			{// || true ) {
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::cout << __PRETTY_FUNCTION__ << std::endl;
				std::cout << njh::bashCT::cyan << "Processing node \n\tuid:" << n->nodeUid_ << " " << n->nameUid_
						                                           << "\n\tk_ :" << n->k_ <<
						njh::bashCT::reset << std::endl;
				std::cout << "Heads: " << njh::bashCT::red << std::endl;
				for(const auto & head : n->headEdges_){
					if(head->on_){
						std::cout << "\tuid: " << head->head_.lock()->nodeUid_ << " " <<  head->head_.lock()->nameUid_ << std::endl;
						std::cout << "\tk_ : " << head->head_.lock()->k_ << std::endl;
					}
				}
				std::cout << njh::bashCT::reset;
				std::cout << "Tails: " << njh::bashCT::blue << std::endl;
				for(const auto & tail : n->tailEdges_){
					if(tail->on_){
						std::cout << "\tuid: " << tail->tail_.lock()->nodeUid_ << " " <<  tail->head_.lock()->nameUid_<< std::endl;
						std::cout << "\tk_ : " << tail->tail_.lock()->k_ << std::endl;
					}
				}
				std::cout << njh::bashCT::reset;
			}
#endif
			for (const auto & head : incoming) {
#if defined(PATHWEAVERSUPERDEBUG)
				{// || true ) {
					std::cout << "\t" << "head: " << head.first << " readsEntering: " << head.second.size() << std::endl;
				}
#endif
				for (const auto & tail : outgoing) {
					uint32_t headReadsCount = 0;
					uint32_t headTotalCount = 0;
					uint32_t tailReadsCount = 0;
					uint32_t tailTotalCount = 0;
					for (const auto & tailCounts : tail.second) {
						if (njh::in(tailCounts.first, head.second)) {
							headTotalCount += head.second.at(tailCounts.first);
							tailTotalCount += tailCounts.second;
							++tailReadsCount;
							++headReadsCount;
						}
					}
					if(tailReadsCount > 0){
						headsToTailsCounts[head.first][tail.first] = tailReadsCount;
					}
					/*
					if (tailReadsCount > occurenceCutOff_) {
						headsToTails[head.first].emplace_back(tail.first);
					}*/
#if defined(PATHWEAVERSUPERDEBUG)
					{
						std::cout << "\t\t" << "tail: " << tail.first << std::endl;
						std::cout << "\t\t\t" << "tailTotal: " << tailReadsCount << " " << tailTotalCount << std::endl;
						std::cout << "\t\t\t" << "headTotal: " << headReadsCount << " " << headTotalCount << std::endl;
					}
#endif
				}
			}

			std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<uint32_t>>> headTailReadsConnecting;
			for (const auto & head : headsToTailsCounts) {
				for (const auto & tail : head.second) {
					for (const auto & name : n->inReadNamesIdx_) {
						bool appearsInHeads = false;
						bool appearsInTails = false;
						if (njh::in(name, incoming[head.first])) {
							appearsInHeads = true;
						}
						if (njh::in(name, outgoing[tail.first])) {
							appearsInTails = true;
						}
						if (appearsInTails && appearsInHeads) {
							headTailReadsConnecting[head.first][tail.first].emplace_back(name);
						}
					}
				}
			}

			std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> uniqueHeadsToTailsCounts;
			for(const auto & headReadCon : headTailReadsConnecting){
				for(const auto & tailReadCon : headReadCon.second){

					for(const auto & readName : tailReadCon.second){
						bool appearsInOtherConnections = false;
						for(const auto & otherHeadReadCon : headTailReadsConnecting){
							for(const auto & otherTailReadCon : otherHeadReadCon.second){
								if(otherHeadReadCon.first == headReadCon.first && otherTailReadCon.first == tailReadCon.first){
									continue;
								}
								if(njh::in(readName, otherTailReadCon.second)){
									appearsInOtherConnections = true;
								}
							}
						}
						if(!appearsInOtherConnections){
							++uniqueHeadsToTailsCounts[headReadCon.first][tailReadCon.first];
						}
					}
				}
			}
#if defined(PATHWEAVERSUPERDEBUG)
			{// || true ) {
				std::cout << "Unique heads to tail connections:  " << std::endl;
				for(const auto & uniHead : uniqueHeadsToTailsCounts){
					std::cout << "head: " <<uniHead.first << std::endl;
					for(const auto & uniTail : uniHead.second){
						std::cout << "\t" << "tail: " << uniTail.first << std::endl;
						std::cout << "\t\t" << "uniqConnectionsCounts: " << uniTail.second << std::endl;
					}
				}
			}
#endif
#if defined(PATHWEAVERSUPERDEBUG)
			{// || true ) {
				std::cout << njh::bashCT::flashing << njh::bashCT::red
						<< "headsToTailsCounts: " << headsToTailsCounts.size()
						<< njh::bashCT::reset << std::endl;
				if (!headsToTailsCounts.empty()) {
					uint32_t headCount = 0;
					for(const auto & head : headsToTailsCounts){
						++headCount;
						if(headCount % 2 == 0){
							std::cout << njh::bashCT::red << std::endl;
						}else{
							std::cout << njh::bashCT::blue << std::endl;
						}
						uint32_t total = 0;
						for(const auto & tail : head.second){
							total+= tail.second;
						}
						std::cout << head.first << " totalBridging: " << total << " totalEntering: " << incoming[head.first].size() << std::endl;
						for(const auto & tail : head.second){
							std::cout << "\t" << tail.first << ": " << tail.second << "/" << total << "; totalEntering: " << incoming[head.first].size() << std::endl;
						}
						std::cout << njh::bashCT::reset << std::endl;
					}
				}


				std::cout << njh::bashCT::flashing << njh::bashCT::red
						<< "uniqueHeadsToTailsCounts: " << uniqueHeadsToTailsCounts.size()
						<< njh::bashCT::reset << std::endl;
				if (!uniqueHeadsToTailsCounts.empty()) {
					uint32_t headCount = 0;
					for(const auto & head : uniqueHeadsToTailsCounts){
						++headCount;
						if(headCount % 2 == 0){
							std::cout << njh::bashCT::red << std::endl;
						}else{
							std::cout << njh::bashCT::blue << std::endl;
						}
						uint32_t total = 0;
						for(const auto & tail : head.second){
							total+= tail.second;
						}
						std::cout << head.first << " totalBridging: " << total << " totalEntering: " << incoming[head.first].size() << std::endl;
						for(const auto & tail : head.second){
							std::cout << "\t" << tail.first << ": " << tail.second << "/" << total << "; totalEntering: " << incoming[head.first].size() << std::endl;
						}
						std::cout << njh::bashCT::reset << std::endl;
					}
				}
			}
#endif
			uint32_t amountBelowCutOff = 0;
			uint32_t amountAboveCutOff = 0;


			std::unordered_map<uint32_t, std::vector<uint32_t>> headsToTails;
			std::set<uint32_t> headsGoingToAdd;
			std::set<uint32_t> tailsGoingToAdd;
			if (!uniqueHeadsToTailsCounts.empty()) {
				for(const auto & head : uniqueHeadsToTailsCounts){
					uint32_t total = 0;
					for(const auto & tail : head.second){
						total+= tail.second;
					}
					for(const auto & tail : head.second){
//						if(tail.second > occurenceCutOff_){
						if(tail.second > connectorCutOff){
							amountAboveCutOff += tail.second;
						}else if(tail.second > 1){
							amountBelowCutOff += tail.second;
						}
					}
					for (const auto & tail : head.second) {
						//if (tail.second > occurenceCutOff_) {
						if (tail.second > connectorCutOff) {
							headsGoingToAdd.emplace(head.first);
							tailsGoingToAdd.emplace(tail.first);
							headsToTails[head.first].emplace_back(tail.first);
						}
					}
					/*
					if(aboveCutOff == 1){
						for(const auto & tail : head.second){
							if(tail.second > occurenceCutOff_){
								headsToTails[head.first].emplace_back(tail.first);
							}
						}
					}else if(total + occurenceCutOff_ >= internalReads.size()){
						for(const auto & tail : head.second){
							if(tail.second > occurenceCutOff_){
								headsToTails[head.first].emplace_back(tail.first);
							}
						}
					}*/
				}
			}
//			std::unordered_map<std::string, std::vector<std::string>> headsToTails;
//
//			if (!headsToTailsCounts.empty()) {
//				for(const auto & head : headsToTailsCounts){
//					uint32_t total = 0;
//					for(const auto & tail : head.second){
//						total+= tail.second;
//					}
//					for(const auto & tail : head.second){
////						if(tail.second > occurenceCutOff_){
//						if(tail.second > connectorCutOff){
//							amountAboveCutOff += tail.second;
//						}else if(tail.second > 1){
//							amountBelowCutOff += tail.second;
//						}
//					}
//					for (const auto & tail : head.second) {
//						//if (tail.second > occurenceCutOff_) {
//						if (tail.second > connectorCutOff) {
//							headsToTails[head.first].emplace_back(tail.first);
//						}
//					}
//					/*
//					if(aboveCutOff == 1){
//						for(const auto & tail : head.second){
//							if(tail.second > occurenceCutOff_){
//								headsToTails[head.first].emplace_back(tail.first);
//							}
//						}
//					}else if(total + occurenceCutOff_ >= internalReads.size()){
//						for(const auto & tail : head.second){
//							if(tail.second > occurenceCutOff_){
//								headsToTails[head.first].emplace_back(tail.first);
//							}
//						}
//					}*/
//				}
//			}
			bool skipBasedOnCloseCutOff = false;
			//if the amount below is equal to 1/3 or more of the one connection then don't break
			if(1 == headsToTails.size() && amountBelowCutOff/static_cast<double>(amountAboveCutOff) >= 0.33){
				skipBasedOnCloseCutOff = true;
			}
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::cout << __PRETTY_FUNCTION__ << std::endl;
				std::cout << "skipBasedOnCloseCutOff: " << njh::colorBool(skipBasedOnCloseCutOff) << std::endl;
				std::cout << "amountBelowCutOff: " << amountBelowCutOff << std::endl;
				std::cout << "amountAboveCutOff: " << amountAboveCutOff << std::endl;
				std::cout << "amountBelowCutOff/static_cast<double>(amountAboveCutOff) >= 0.33: " << njh::colorBool(amountBelowCutOff/static_cast<double>(amountAboveCutOff) >= 0.33) << std::endl;
			}
#endif
			std::set<uint32_t> headsGoingToMiss;
			std::set<uint32_t> tailsGoingToMiss;
			for(const auto & head : n->headEdges_){
				auto headUID = head->head_.lock()->nodeUid_;
				if(!njh::in(headUID, headsGoingToAdd)){
					headsGoingToMiss.emplace(head->head_.lock()->nodeUid_);
				}
			}
			for(const auto & tail : n->tailEdges_){
				auto tailUID = tail->tail_.lock()->nodeUid_;
				if(!njh::in(tailUID, tailsGoingToAdd)){
					tailsGoingToMiss.emplace(tailUID);
				}
			}
			bool skipBasedOnConservative = false;
			if(pars.conservative_){
				skipBasedOnConservative = !headsGoingToMiss.empty() || !tailsGoingToMiss.empty();
			}

			if (!headsToTails.empty() && !skipBasedOnCloseCutOff && !skipBasedOnConservative) {
/*
 * 				std::ofstream outTestFile("test_tail", std::ios::app);
				for (auto & tailEdge : n->tailEdges_) {
					if (tailEdge->tail_.lock()->k_ == "TTCATGAAGGAAAAAATTTAAAAACTTCCAATAAAAAAAAAAATGATGACAATAATTCAAAATTATGCAAAGCTTTAAAATACAGTTTTGCTGATTATGGAGATTTA") {
						outTestFile << std::endl;
						outTestFile << "start" << std::endl;
						printVector( tailEdge->readNames_, ", ", outTestFile);
						outTestFile << "end" << std::endl;
					}
				}
 */
				std::unordered_map<uint32_t, std::unordered_map<uint32_t, std::vector<uint32_t>>> headTailConnections;
				std::unordered_map<uint32_t, uint32_t> headbeingAddedCounts;
				std::unordered_map<uint32_t, uint32_t> tailBeingAddedCounts;
				for (const auto & head : headsToTails) {
					for (const auto & tail : head.second) {
						for (const auto & name : n->inReadNamesIdx_) {
							bool appearsInHeads = false;
							bool appearsInTails = false;
							if (njh::in(name, incoming[head.first])) {
								appearsInHeads = true;
							}
							if (njh::in(name, outgoing[tail])) {
								appearsInTails = true;
							}
							if (appearsInTails && appearsInHeads) {
								headTailConnections[head.first][tail].emplace_back(name);
							}
							++headbeingAddedCounts[head.first];
							++tailBeingAddedCounts[tail];
						}
					}
				}
				std::set<uint32_t> addedHeads;
				std::set<uint32_t> addedTails;
				std::vector<std::shared_ptr<node>> newNodes;
				for (const auto & head : headsToTails) {
					for (const auto & tail : head.second) {
						std::shared_ptr<KmerPathwayGraphDev::node> newNode = std::make_shared<KmerPathwayGraphDev::node>(n->k_, n->cnt_, n->kLen_);
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << __PRETTY_FUNCTION__ << std::endl;
							std::cout << njh::bashCT::cyan << "Creating new node for the following head and tail \n\t" << newNode->nodeUid_  << " "<< newNode->nameUid_ << "\n\t" << newNode->k_ <<
									njh::bashCT::reset << std::endl;
							std::cout << "\t" << njh::bashCT::red <<  "head: " << head.first << njh::bashCT::reset << std::endl;
							std::cout << "\t" << njh::bashCT::blue << "tail: " << tail << njh::bashCT::reset << std::endl;

						}
#endif
						std::unordered_set<uint32_t> namesForInternal;
						std::unordered_set<uint32_t> namesForEdges;
						uint32_t appearsInOtherHeadsCounts = 0;
						uint32_t appearsInOtherTailsCounts = 0;
						for (const auto & name : n->inReadNamesIdx_) {
							bool appearsInHeads = false;
							bool appearsInTails = false;
							bool appearsInOtherTails = false;
							bool appearsInOtherHeads = false;
							if (njh::in(name, incoming[head.first])) {
								appearsInHeads = true;
							}
							if (njh::in(name, outgoing[tail])) {
								appearsInTails = true;
							}


							bool appearsInOtherHeadToTailConnections = false;
							for(const auto & headCon : headTailConnections){
								for(const auto & tailCon : headCon.second){
									if(head.first == headCon.first && tail == tailCon.first){
										continue;
									}
									if(njh::in(name, tailCon.second)){
										appearsInOtherHeadToTailConnections = true;
									}
								}
							}
							if (appearsInOtherHeadToTailConnections) {
								continue;
							} else {
								namesForInternal.emplace(name);
							}
							//if this tail is being added multiple times then the reads in the tail have to appear in the head
							if(tailBeingAddedCounts[tail] > 1 && !appearsInHeads){
								continue;
							}
							//if this head is being added multiple times then the reads in the head have to appear in the tail
							if(headbeingAddedCounts[head.first] > 1 && !appearsInTails){
								continue;
							}
							//getting rid of this, with the check that the reads don't appear in other heads or tails i think it should be find
							//this will then get reads only going through this head and this tail - njh 2020-03-05
//							if (!appearsInTails || !appearsInHeads) {
//								continue;
//							}


							for (const auto & other : outgoing) {
								if (other.first != tail && njh::in(name, other.second)) {
									appearsInOtherTails = true;
									++appearsInOtherTailsCounts;
#if defined(PATHWEAVERSUPERDEBUG)
									{
										if (appearsInTails) {
											std::stringstream errorStream;
											errorStream << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << std::endl;
											errorStream << name << " appears in other tails" << std::endl;
											errorStream << "other.first != tail "
													<< njh::colorBool(other.first != tail) << std::endl;
											errorStream << "njh::in(name, other.second) "
													<< njh::colorBool(njh::in(name, other.second))
													<< std::endl;
											errorStream << "other: " << other.first << std::endl;
											errorStream << "currentTail: " << tail << std::endl;
											errorStream << "appearsInTails: "
													<< njh::colorBool(appearsInTails) << std::endl;
											errorStream << "appearsInHeads: "
													<< njh::colorBool(appearsInHeads) << std::endl;
											printOutMapContents(other.second, "\t", errorStream);
											//throw std::runtime_error {errorStream.str() };
										}
									}
#endif
								}
							}
							for(const auto & other : incoming){
								if(other.first != head.first &&
										njh::in(name, other.second)){
									appearsInOtherHeads = true;
									++appearsInOtherHeadsCounts;
#if defined(PATHWEAVERSUPERDEBUG)
									{
										if(appearsInHeads){
											std::stringstream errorStream;
											errorStream << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << std::endl;
											errorStream << name << " appears in other heads" << std::endl;
											errorStream << "other.first != head " << njh::colorBool(other.first != head.first) << std::endl;
											errorStream << "njh::in(name, other.second) " << njh::colorBool(njh::in(name, other.second)) << std::endl;
											errorStream << "other: " << other.first << std::endl;
											errorStream << "currentHead: " << head.first << std::endl;
											errorStream << "appearsInTails: " << njh::colorBool(appearsInTails) << std::endl;
											errorStream << "appearsInHeads: " << njh::colorBool(appearsInHeads) << std::endl;
											printOutMapContents(other.second, "\t", errorStream);
											//throw std::runtime_error { errorStream.str() };
										}
									}
#endif
								}
							}
							//if(!appearsInOtherTails && !appearsInOtherHeads && appearsInTails && appearsInHeads){
						//	if(!appearsInOtherTails && !appearsInOtherHeads && (appearsInTails || appearsInHeads)){
							if((appearsInTails || appearsInHeads)){
								namesForEdges.emplace(name);
							}
						}
						//add read names
						if(namesForEdges.empty()){
							std::stringstream errorStream;
							errorStream << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << std::endl;
							errorStream << "Error, names is empty: " << newNode->k_ << std::endl;
							errorStream << "Was expecting: " << headsToTailsCounts[head.first][tail] << std::endl;
							errorStream << "head: " << head.first << std::endl;
							errorStream << "tail: " << tail << std::endl;
							errorStream << "n->inReadNamesIdx_.size(): " << n->inReadNamesIdx_.size() << "\n";
							//errorStream << "n->inReadNamesIdx_: " << njh::conToStr(n->inReadNamesIdx_, ", ") << "\n";
							errorStream << "appearsInOtherHeadsCounts: " << appearsInOtherHeadsCounts << "\n";
							errorStream << "appearsInOtherTailsCounts: " << appearsInOtherTailsCounts << "\n";
							throw std::runtime_error{errorStream.str()};
						}
						newNode->inReadNamesIdx_ = namesForInternal;
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << "\tNew node read size: " << newNode->inReadNamesIdx_.size() << std::endl;
						}
#endif
						//add edges and modifying head and tail nodes
						for (auto & headEdge : n->headEdges_) {
							if(!headEdge->on_){
								continue;
							}
							if (head.first == headEdge->head_.lock()->nodeUid_) {
								addedHeads.emplace(head.first);
								std::unordered_set<uint32_t> edgeNames;
								for (const auto & edgeName : headEdge->inReadNamesIdx_) {
									if (njh::in(edgeName, namesForEdges)) {
										edgeNames.emplace(edgeName);
									}
								}
								if(edgeNames.empty()){
									std::stringstream errorStream;
									errorStream << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << std::endl;
									errorStream << "Error, edgeNames is empty: " << newNode->k_ << std::endl;
									errorStream << "for head: " << head.first << std::endl;
									errorStream << "Was expecting: " << incoming[head.first].size() << std::endl;
									errorStream << "InternalNames: " << njh::conToStr(namesForInternal, ",") << std::endl;
									errorStream << "EdgeNames:     " << njh::conToStr(namesForEdges, ",")    << std::endl;
									errorStream << "Available Edge Names: " << njh::conToStr(headEdge->inReadNamesIdx_, ",") << std::endl;
									throw std::runtime_error{errorStream.str()};
								}
								auto addingEdge = std::make_shared<KmerPathwayGraphDev::edge>(
										headEdge->head_.lock(), newNode, edgeNames.size(), edgeNames);
								newNode->headEdges_.push_back(addingEdge);
#if defined(PATHWEAVERSUPERDEBUG)
								{
									std::cout << "\tNew node head " << head.first << " read size: " << addingEdge->inReadNamesIdx_.size() << std::endl;
								}
#endif
								//doing this in one area so it's cleaerer what's being marked, see below
								//headEdge->head_.lock()->visitCount_ += 1;
								headEdge->head_.lock()->tailEdges_.push_back(addingEdge);
								break;
							}
						}
						for (auto & tailEdge : n->tailEdges_) {
							if(!tailEdge->on_){
								continue;
							}
							if (tail == tailEdge->tail_.lock()->nodeUid_) {
								addedTails.emplace(tail);
								std::unordered_set<uint32_t> edgeNames;
								for (const auto & edgeName : tailEdge->inReadNamesIdx_) {
									if (njh::in(edgeName, namesForEdges)) {
										edgeNames.emplace(edgeName);
									}
								}
								if(edgeNames.empty()){
									std::stringstream errorStream;
									errorStream << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << std::endl;
									errorStream << "Error, edgeNames is empty: " << newNode->k_ << std::endl;
									errorStream << "for tail: " << tail << std::endl;
									errorStream << "Was expecting: " << outgoing[tail].size() << std::endl;
									errorStream << "InternalNames: " << njh::conToStr(namesForInternal, ",") << std::endl;
									errorStream << "EdgeNames:     " << njh::conToStr(namesForEdges, ",")    << std::endl;
									errorStream << "Available Edge Names: " << njh::conToStr(tailEdge->inReadNamesIdx_, ",") << std::endl;

									errorStream << "Node check: " << std::endl;
									uint32_t nodeCount = 0;
									for(const auto & n : nodes_){
										++nodeCount;
										if("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCT" == n->nameUid_){
											errorStream << "nodeCount: " << nodeCount << std::endl;
											errorStream << "n->k_   :"<< n->k_ << std::endl;
											errorStream << "n->cnt_ :"<< n->cnt_ << std::endl;
										}
									}
									throw std::runtime_error{errorStream.str()};
								}
								auto addingEdge = std::make_shared<KmerPathwayGraphDev::edge>(newNode,
										tailEdge->tail_.lock(), edgeNames.size(), edgeNames);
								newNode->tailEdges_.push_back(addingEdge);
#if defined(PATHWEAVERSUPERDEBUG)
								{
									std::cout << "\tNew node tail " << tail << " read size: " << addingEdge->inReadNamesIdx_.size() << std::endl;
								}
#endif
								//doing this in one area so it's cleaerer what's being marked, see below
								//tailEdge->tail_.lock()->visitCount_ += 1;
								tailEdge->tail_.lock()->headEdges_.push_back(addingEdge);
								break;
							}
						}
						//add node to graph.nodes_
						nodesSplit = true;
//						nodes_.push_back(newNode);
						//newNode->uid_ = njh::pasteAsStr(newNode->uid_, newNodes.size());
						newNode->nameUid_ = njh::pasteAsStr(newNode->nameUid_, newNodes.size());
						newNode->nodeUid_ = nodes_.size() + newNodes.size();
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << __PRETTY_FUNCTION__ << std::endl;
							std::cout << "New node before adding to new nodes: " << std::endl;
							std::cout << "newNode: " << newNode->nodeUid_ << " " << newNode->nameUid_ << std::endl;
							std::cout << "\thead: " << newNode->getFirstOnHeadEdge()->head_.lock()->nodeUid_  << " " << newNode->getFirstOnHeadEdge()->head_.lock()->nameUid_ << std::endl;
							std::cout << "\ttail: " << newNode->getFirstOnTailEdge()->tail_.lock()->nodeUid_  << " " << newNode->getFirstOnTailEdge()->tail_.lock()->nameUid_  << std::endl;

						}
#endif
						newNodes.emplace_back(newNode);
						/*
						 *
						printVector(n->readNames_);
						printVector(newNode->readNames_	);
						*/
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << __PRETTY_FUNCTION__ << std::endl;
							std::cout << "Current new nodes size: " << newNodes.size() << std::endl;
							uint32_t possibleNewNodePos = 0;
							for (const auto & possilbeNewNode : newNodes){
								if(possibleNewNodePos % 2 == 0){
									std::cout << njh::bashCT::blue << std::endl;
								}else{
									std::cout << njh::bashCT::red << std::endl;
								}
								std::cout << "newNode: " << possilbeNewNode->nodeUid_ << " " << possilbeNewNode->nameUid_ << " " << possibleNewNodePos << std::endl;
								std::cout << "\thead: " << possilbeNewNode->getFirstOnHeadEdge()->head_.lock()->nodeUid_ << " " << possilbeNewNode->getFirstOnHeadEdge()->head_.lock()->nameUid_ << std::endl;
								std::cout << "\ttail: " << possilbeNewNode->getFirstOnTailEdge()->tail_.lock()->nodeUid_ << " " << possilbeNewNode->getFirstOnTailEdge()->tail_.lock()->nameUid_ << std::endl;
								std::cout << njh::bashCT::reset;
								++possibleNewNodePos;
							}
						}
#endif

#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << __PRETTY_FUNCTION__ << std::endl;
							std::cout << njh::bashCT::cyan << "Adding new node " << newNode->nodeUid_ << " " << newNode->nameUid_ << " " << newNode->k_ <<
									njh::bashCT::reset << std::endl;
							std::cout << njh::bashCT::purple << "With head edges:  ";
							for(const auto & head : newNode->headEdges_){
								std::cout << head->head_.lock()->nodeUid_ << " "  << head->head_.lock()->nameUid_ << " ";
							}
							std::cout <<njh::bashCT::reset << std::endl;
							std::cout << njh::bashCT::green << "With tail edges:  ";
							for(const auto & tail : newNode->tailEdges_){
								std::cout << tail->tail_.lock()->nodeUid_ << " " << tail->tail_.lock()->nameUid_<< " ";
							}
							std::cout <<njh::bashCT::reset << std::endl;
						}
#endif
					}
				}
//				for (const auto & newNode : newNodes) {
//					nodes_.emplace_back(newNode);
//				}
				//throw std::runtime_error{"stopping"};
				//turn off node and edges
				bool didntAdd = false;

				for(auto & headEdge : n->headEdges_){
					if(!headEdge->on_){
						continue;
					}
					if(!njh::in(headEdge->head_.lock()->nodeUid_, addedHeads)){
						didntAdd = true;
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << njh::bashCT::flashing;
							std::cout << "Didn't add head: " << headEdge->head_.lock()->k_ << std::endl;
							std::cout << "read names size: " << headEdge->inReadNamesIdx_.size()<< std::endl;
							std::cout << njh::bashCT::reset;
						}
#endif
					}else{
						headEdge->on_ = false;
					}
				}

				for(auto & tailEdge : n->tailEdges_){
					if(!tailEdge->on_){
						continue;
					}
					if(!njh::in(tailEdge->tail_.lock()->nodeUid_, addedTails)){
						didntAdd = true;
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << njh::bashCT::flashing;
							std::cout << "Didn't add tail: " << tailEdge->tail_.lock()->k_ << std::endl;
							std::cout << "read names size: " << tailEdge->inReadNamesIdx_.size()  << std::endl;
							std::cout << njh::bashCT::reset;
						}
#endif
					}else{
						tailEdge->on_ = false;
					}
				}
#if defined(PATHWEAVERSUPERDEBUG)
				{
					std::cout << "Added heads: " << std::endl;
					for(const auto & addedHead : addedHeads){
						std::cout << "\thead:" << addedHead << std::endl;
						std::cout << "\twhich has head count of " << nodes_[addedHead]->headCount() << std::endl;
						std::cout << "\twhich has tail count of " << nodes_[addedHead]->tailCount() << std::endl;
						std::cout << "\theads:" << std::endl;
						for(const auto & head : nodes_[addedHead]->headEdges_){
							if(head->on_){
								std::cout << "\t\t" <<head->head_.lock()->nodeUid_ << " " << head->head_.lock()->nameUid_ << " with count: " << head->cnt_ << std::endl;
							}
						}
						std::cout << "\ttails:" << std::endl;
						for(const auto & tail : nodes_[addedHead]->tailEdges_){
							if(tail->on_){
								std::cout << "\t\t" <<tail->tail_.lock()->nodeUid_ << " " << tail->tail_.lock()->nameUid_ << " with count: " << tail->cnt_ << std::endl;
							}
						}
						std::cout << std::endl;
					}

					std::cout << "Added tails: " << std::endl;
					for(const auto & addedTail : addedTails){
						std::cout << "\ttail:" << addedTail << std::endl;
						std::cout << "\twhich has head count of " << nodes_[addedTail]->headCount() << std::endl;
						std::cout << "\twhich has tail count of " << nodes_[addedTail]->tailCount() << std::endl;
						std::cout << "\theads:" << std::endl;
						for(const auto & head : nodes_[addedTail]->headEdges_){
							if(head->on_){
								std::cout << "\t\t" <<head->head_.lock()->nodeUid_  << " " << head->head_.lock()->nameUid_ << " with count: " << head->cnt_ << std::endl;
							}
						}
						std::cout << "\ttails:" << std::endl;
						for(const auto & tail : nodes_[addedTail]->tailEdges_){
							if(tail->on_){
								std::cout << "\t\t" <<tail->tail_.lock()->nodeUid_  << " " << tail->tail_.lock()->nameUid_ << " with count: " << tail->cnt_ << std::endl;
							}
						}
						std::cout << std::endl;
					}
				}
#endif
				//handling a very specific situation here brought on by tandem repeats normally
				//what should be done is checking how many times through the loop or if reads span the whole thing
				if (newNodes.size() == 2) {
#if defined(PATHWEAVERSUPERDEBUG)
					{
						std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
						std::cout << "\tnewNodes.front()->tailCount(): " << newNodes.front()->tailCount() << std::endl;
						std::cout << "\tnewNodes.front()->headCount(): " << newNodes.front()->headCount() << std::endl;
						std::cout << "\tnewNodes.back()->tailCount(): " << newNodes.back()->tailCount() << std::endl;
						std::cout << "\tnewNodes.back()->headCount(): " << newNodes.back()->headCount() << std::endl;
					}
#endif
					if((1 == newNodes.front()->tailCount() && 1 == newNodes.front()->headCount()) &&
						 (1 == newNodes.back()->tailCount()  && 1 == newNodes.back()->headCount())){
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
							std::cout << "Checking for possible loop nodes" << std::endl;
							uint32_t possibleNewNodePos = 0;
							for (const auto & possilbeNewNode : newNodes){
								if(possibleNewNodePos % 2 == 0){
									std::cout << njh::bashCT::blue << std::endl;
								}else{
									std::cout << njh::bashCT::red << std::endl;
								}
								std::cout << "newNode: " << possilbeNewNode->nodeUid_  << " " << possilbeNewNode->nameUid_ << std::endl;
								std::cout << "\thead: " <<  possilbeNewNode->getFirstOnHeadEdge()->head_.lock()->nodeUid_ << " " << possilbeNewNode->getFirstOnHeadEdge()->head_.lock()->nameUid_ << std::endl;
								std::cout << "\ttail: " <<  possilbeNewNode->getFirstOnTailEdge()->tail_.lock()->nodeUid_ << " " << possilbeNewNode->getFirstOnTailEdge()->tail_.lock()->nameUid_ << std::endl;
								std::cout << njh::bashCT::reset;
								++possibleNewNodePos;
							}
						}
#endif
						if(newNodes.front()->getFirstOnTailEdge()->tail_.lock()->nodeUid_ == newNodes.back()->getFirstOnHeadEdge()->head_.lock()->nodeUid_ &&
							 newNodes.front()->getFirstOnHeadEdge()->head_.lock()->nodeUid_ == newNodes.back()->getFirstOnTailEdge()->tail_.lock()->nodeUid_){
							{
								//tail break
#if defined(PATHWEAVERSUPERDEBUG)
								{
									std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
								}
#endif
								auto loopNode = newNodes.front()->getFirstOnTailEdge()->tail_.lock();
								if(1 == loopNode->tailCount()  && 1 == loopNode->headCount()){
#if defined(PATHWEAVERSUPERDEBUG)
									{
										std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
										std::cout << "Found loop node: " << loopNode->nodeUid_ << " " << loopNode->nameUid_ << std::endl;
									}
#endif
									//old node will keep head, new node will have tail
									std::shared_ptr<KmerPathwayGraphDev::node> newLoopNode = std::make_shared<KmerPathwayGraphDev::node>(loopNode->k_, loopNode->cnt_, loopNode->kLen_);
									newLoopNode->inReadNamesIdx_ = loopNode->inReadNamesIdx_;
									//connect the old tail to the new node (should only be one so loop until you find on edge)
									std::shared_ptr<edge> headEdge;
									for(const auto & head : newNodes.back()->headEdges_ ){
										if(head->on_){
											headEdge = head;
											break;
										}
									}
									headEdge->head_ = newLoopNode;
									newLoopNode->tailEdges_.emplace_back(headEdge);
									//clear the tails of the old node to break the connection
									loopNode->tailEdges_.clear();
									//add new node
									newLoopNode->nameUid_ = njh::pasteAsStr(newLoopNode->nameUid_, newNodes.size());
									newLoopNode->nodeUid_ = nodes_.size() + newNodes.size();
									newNodes.emplace_back(newLoopNode);
								}
							}
							{
								//head break
#if defined(PATHWEAVERSUPERDEBUG)
								{
									std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
								}
#endif
								auto loopNode = newNodes.front()->getFirstOnHeadEdge()->head_.lock();
#if defined(PATHWEAVERSUPERDEBUG)
								{
									std::cout << "\tloopNode->tailCount(): " << loopNode->tailCount() << std::endl;
									std::cout << "\tloopNode->headCount(): " << loopNode->headCount() << std::endl;
								}
#endif
								if(1 == loopNode->tailCount()  && 1 == loopNode->headCount()){
#if defined(PATHWEAVERSUPERDEBUG)
									{
										std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
										std::cout << "Found loop node: " << loopNode->nodeUid_ << " " << loopNode->nameUid_ << std::endl;
									}
#endif
									//old node will keep head, new node will have tail
									std::shared_ptr<KmerPathwayGraphDev::node> newLoopNode = std::make_shared<KmerPathwayGraphDev::node>(loopNode->k_, loopNode->cnt_, loopNode->kLen_);
									newLoopNode->inReadNamesIdx_ = loopNode->inReadNamesIdx_;
									//connect the old tail to the new node (should only be one so loop until you find on edge)
									std::shared_ptr<edge> headEdge;
									for(const auto & head : newNodes.front()->headEdges_ ){
										if(head->on_){
											headEdge = head;
											break;
										}
									}
									headEdge->head_ = newLoopNode;
									newLoopNode->tailEdges_.emplace_back(headEdge);
									//clear the tails of the old node to break the connection
									loopNode->tailEdges_.clear();
									//add new node
									newLoopNode->nameUid_ = njh::pasteAsStr(newLoopNode->nameUid_, newNodes.size());
									newLoopNode->nodeUid_ = nodes_.size() + newNodes.size();
									newNodes.emplace_back(newLoopNode);
								}
							}
						}
//						if(newNodes.front()->getFirstOnTailEdge()->tail_.lock()->uid_ == newNodes.back()->getFirstOnHeadEdge()->head_.lock()->uid_){
//							if(debug_){
//								std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//							}
//							auto loopNode = newNodes.front()->getFirstOnTailEdge()->tail_.lock();
//							if(1 == loopNode->tailCount()  && 1 == loopNode->headCount()){
//								if(debug_){
//									std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//									std::cout << "Found loop node: " << loopNode->uid_ << std::endl;
//								}
//								//old node will keep head, new node will have tail
//								std::shared_ptr<KmerPathwayGraphDev::node> newLoopNode = std::make_shared<KmerPathwayGraphDev::node>(loopNode->k_, loopNode->cnt_, loopNode->kLen_);
//								newLoopNode->inReadNamesIdx_ = loopNode->inReadNamesIdx_;
//								//connect the old tail to the new node (should only be one so loop until you find on edge)
//								std::shared_ptr<edge> headEdge;
//								for(const auto & head : newNodes.back()->headEdges_ ){
//									if(head->on_){
//										headEdge = head;
//										break;
//									}
//								}
//								headEdge->head_ = newLoopNode;
//								newLoopNode->tailEdges_.emplace_back(headEdge);
//								//clear the tails of the old node to break the connection
//								loopNode->tailEdges_.clear();
//								//add new node
//								newNodes.emplace_back(newLoopNode);
//							}
//						} else if (newNodes.front()->getFirstOnHeadEdge()->head_.lock()->uid_ == newNodes.back()->getFirstOnTailEdge()->tail_.lock()->uid_){
//							if(debug_){
//								std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//							}
//							auto loopNode = newNodes.front()->getFirstOnHeadEdge()->head_.lock();
//							if(debug_ ){
//								std::cout << "\tloopNode->tailCount(): " << loopNode->tailCount() << std::endl;
//								std::cout << "\tloopNode->headCount(): " << loopNode->headCount() << std::endl;
//							}
//							if(1 == loopNode->tailCount()  && 1 == loopNode->headCount()){
//								if(debug_){
//									std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
//									std::cout << "Found loop node: " << loopNode->uid_ << std::endl;
//								}
//								//old node will keep head, new node will have tail
//								std::shared_ptr<KmerPathwayGraphDev::node> newLoopNode = std::make_shared<KmerPathwayGraphDev::node>(loopNode->k_, loopNode->cnt_, loopNode->kLen_);
//								newLoopNode->inReadNamesIdx_ = loopNode->inReadNamesIdx_;
//								//connect the old tail to the new node (should only be one so loop until you find on edge)
//								std::shared_ptr<edge> headEdge;
//								for(const auto & head : newNodes.front()->headEdges_ ){
//									if(head->on_){
//										headEdge = head;
//										break;
//									}
//								}
//								headEdge->head_ = newLoopNode;
//								newLoopNode->tailEdges_.emplace_back(headEdge);
//								//clear the tails of the old node to break the connection
//								loopNode->tailEdges_.clear();
//								//add new node
//								newNodes.emplace_back(newLoopNode);
//							}
//						}
					}
				}
				std::unordered_set<uint32_t> readNamesInEdgesOfNewNodes;
				for (auto & newNode : newNodes) {
#if defined(PATHWEAVERSUPERDEBUG)
					{
						if(newNode->headCount() == 1 && newNode->tailCount() == 1){
							std::cout<< "New node has one head and one tail uid: " << newNode->nodeUid_ << " " << newNode->nameUid_  << " k_:" << newNode->k_ << std::endl;
							auto headNodeTailCount = newNode->getFirstOnHeadEdge()->head_.lock()->tailCount();
							auto tailNodeHeadCount = newNode->getFirstOnTailEdge()->tail_.lock()->headCount();
							std::cout << "\tit's single head has a tail count of " << headNodeTailCount << std::endl;
							std::cout << "\tit's single tail has a head count of " << tailNodeHeadCount << std::endl;
						}
					}
#endif
					//mark the new nodes's tails and heads as visited
					newNode->visitCount_ += 1;
					for(const auto & tail : newNode->tailEdges_){
						if(tail->on_){
							readNamesInEdgesOfNewNodes.insert(tail->inReadNamesIdx_.begin(), tail->inReadNamesIdx_.end());
							tail->tail_.lock()->visitCount_ += 1;
						}
					}
					for(const auto & head : newNode->headEdges_){
						if(head->on_){
							readNamesInEdgesOfNewNodes.insert(head->inReadNamesIdx_.begin(), head->inReadNamesIdx_.end());
							head->head_.lock()->visitCount_ += 1;
						}
					}
					newNode->nameUid_ = njh::pasteAsStr(newNode->nameUid_, nodes_.size());
					newNode->nodeUid_ = nodes_.size();
					nodePositions_[newNode->nameUid_] = nodes_.size();
					nodes_.emplace_back(newNode);
				}
				if(didntAdd){
					std::unordered_set<uint32_t> readNamesKeepInKeptHeadsTails;
					std::unordered_set<uint32_t> readNamesKeepNotInOtherHeadToTailCons;
					std::vector<std::shared_ptr<edge>> keptTailsEdges;
					std::vector<std::shared_ptr<edge>> keptHeadsEdges;
					std::vector<uint32_t> keptTails;
					std::vector<uint32_t> keptHeads;
					n->visitCount_ += 1;

					for (const auto & e : n->tailEdges_) {
						if (e->on_) {
							//mark any still on head or tail nodes as visited as to get to here this had to have been node was processed
							e->tail_.lock()->visitCount_ += 1;
							keptTailsEdges.push_back(e);
							keptTails.emplace_back(e->tail_.lock()->nodeUid_);
						}
					}
					for (const auto & e : n->headEdges_) {
						if (e->on_) {
							//mark any still on head or tail nodes as visited as to get to here this had to have been node was processed
							e->head_.lock()->visitCount_ += 1;
							keptHeadsEdges.push_back(e);
							keptHeads.emplace_back(e->head_.lock()->nodeUid_);
						}
					}
					for (const auto & name : n->inReadNamesIdx_) {
						for (const auto & e : keptHeadsEdges) {
							if (njh::in(name, e->inReadNamesIdx_)) {
								readNamesKeepInKeptHeadsTails.emplace(name);
								break;
							}
						}
						for (const auto & e : keptTailsEdges) {
							if (njh::in(name, e->inReadNamesIdx_)) {
								readNamesKeepInKeptHeadsTails.emplace(name);
								break;
							}
						}
					}

					for (const auto & name : n->inReadNamesIdx_) {
						if(!njh::in(name, readNamesInEdgesOfNewNodes)){
							readNamesKeepNotInOtherHeadToTailCons.emplace(name);
						}
					}

					if(readNamesKeepNotInOtherHeadToTailCons.empty()){
#if defined(PATHWEAVERSUPERDEBUG)
						{
							for(const auto & name : n->inReadNamesIdx_){
								for(const auto & kh : keptHeadsEdges){
									if(njh::in(name,kh->inReadNamesIdx_ )){
										for(const auto & otherHead : n->headEdges_){
											if(!otherHead->on_){
												continue;
											}
											if(otherHead->head_.lock()->nodeUid_ != kh->head_.lock()->nodeUid_){
												if(njh::in(name, otherHead->inReadNamesIdx_)){
													std::cout << name << " found in kept head: " << kh->head_.lock()->nodeUid_ << " " << kh->head_.lock()->nameUid_ << " and in " << otherHead->head_.lock()->nodeUid_ << " " << otherHead->head_.lock()->nameUid_ << std::endl;
												}
											}
										}
									}
								}
							}
						}
#endif
						std::stringstream errorStream;
						errorStream << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << std::endl;
						errorStream << "Error, readNamesKeep is empty: " << n->k_ << std::endl;
						errorStream << "Cnt: " << n->cnt_ << std::endl;
						errorStream << "read cnt: " << n->inReadNamesIdx_.size() << std::endl;
						errorStream << "for tails: " << njh::conToStr(keptTails, ", ") << " " << keptTails.size() << std::endl;
						errorStream << "for heads: " << njh::conToStr(keptHeads, ", ") << " " << keptHeads.size() << std::endl;
						throw std::runtime_error{errorStream.str()};
					}
					n->inReadNamesIdx_ = readNamesKeepNotInOtherHeadToTailCons;
					std::vector<std::shared_ptr<edge>> keptHead;
					for(const auto & head : n->headEdges_){
						if(head->on_){
							keptHead.emplace_back(head);
						}
					}
					std::vector<std::shared_ptr<edge>> keptTail;
					for(const auto & tail : n->tailEdges_){
						if(tail->on_){
							keptTail.emplace_back(tail);
						}
					}
#if defined(PATHWEAVERSUPERDEBUG)
					{//if (debug_) { //if (debug_) {// || true ) {
						std::cout << "node read names after size: " << n->inReadNamesIdx_.size()  << std::endl;
						//throw std::runtime_error{"stopping"};
					}
#endif

#if defined(PATHWEAVERSUPERDEBUG)
					{ //if (debug_) {// || true ) {
						std::cout << "\tkeptTail.size(): " << keptTail.size() << std::endl;
						std::cout << "\tn->tailCount(): " << n->tailCount() << ", n->tailEdges_.size(): " << n->tailEdges_.size() << std::endl;
						std::cout << "\tkeptHead.size(): " << keptHead.size() << std::endl;
						std::cout << "\tn->headCount(): " << n->headCount() << ", n->headEdges_.size(): " << n->headEdges_.size() << std::endl;
					}
#endif

					//check to see if an artificial connection was created due to losing some edges due to low coverage
					if (1 == n->tailCount() && 1 == n->headCount()) {
						//next should only have a head count of 1 for it to be here
						std::set<uint32_t> readsEntering(
								keptHead.front()->inReadNamesIdx_.begin(),
								keptHead.front()->inReadNamesIdx_.end());
						//and to reach here tailEdges_ should only be one as well
						std::set<uint32_t> readsLeaving(
								keptTail.front()->inReadNamesIdx_.begin(),
								keptTail.front()->inReadNamesIdx_.end());
						std::vector<uint32_t> readsThatEnteredLeaving;
						std::set_intersection(readsEntering.begin(), readsEntering.end(),
								readsLeaving.begin(), readsLeaving.end(),
								std::back_inserter(readsThatEnteredLeaving));
#if defined(PATHWEAVERSUPERDEBUG)
						{//if (debug_) { //if (debug_) {// || true ) {
							std::cout << "Checking For uid: " << n->nodeUid_ << " " << n->nameUid_ << std::endl;;
							std::cout << "\tReadsEnteringNumber    : "
									<< readsEntering.size() << std::endl;
							std::cout << "\tReadsLeavingNumber     : "
									<< readsLeaving.size() << std::endl;
							std::cout << "\tReadsLeavingThatEntered: "
									<< readsThatEnteredLeaving.size() << std::endl;
							std::cout << std::endl;
						}
#endif
						if (readsThatEnteredLeaving.size()
								!= std::min(readsEntering.size(), readsLeaving.size()) &&
						//					readsThatEnteredLeaving.size() < 2){
								readsThatEnteredLeaving.size() <= connectorCutOff) {
							//readsThatEnteredLeaving.size() <= occurenceCutOff_) {
#if defined(PATHWEAVERSUPERDEBUG)
							{//if (debug_) { //if (debug_) {// || true ) {
								std::cout << njh::bashCT::flashing;
								std::cout << "Breaking For uid: " << n->nodeUid_ << " " << n->nameUid_ << std::endl;;
								std::cout << "\t\tReadsEnteringNumber    : "
										<< readsEntering.size() << std::endl;
								std::cout << "\t\tReadsLeavingNumber     : "
										<< readsLeaving.size() << std::endl;
								std::cout << "\t\tReadsLeavingThatEntered: "
										<< readsThatEnteredLeaving.size() << std::endl;
								std::cout << std::endl;
								std::cout << njh::bashCT::reset;
							}
#endif
//							//create new node, this will be the node with the tails, old node will have the heads
//							std::shared_ptr<KmerPathwayGraphDev::node> newNode =
//									std::make_shared<KmerPathwayGraphDev::node>(n->k_, n->cnt_, n->kLen_);
//							newNode->inReadNamesIdx_ = n->inReadNamesIdx_;
//							//reset tail edges
//							newNode->tailEdges_ = n->tailEdges_;
//							n->tailEdges_.clear();
//							for (auto & tailEdge : newNode->tailEdges_) {
//								tailEdge->head_ = newNode;
//							}
//							newNode->uid_ = njh::pasteAsStr(newNode->uid_, nodes_.size());
//							nodePositions_[newNode->uid_] = nodes_.size();
//							nodes_.emplace_back(newNode);
							//turn off head and tail node to break connection
							//there should only be one head and one tail node
							auto tail = n->getFirstOnTailEdge();
							auto head = n->getFirstOnHeadEdge();
							tail->on_ = false;
							head->on_ = false;
							//it's only here if we added this node already to something else,
							//if we are turning off it's remaining heads and tails then we might as well turn of the node as well
							if(pars.throwAwayConservedAddNodesDuringDisentaglement_){
								n->on_ = false;
							} else {
								//if the head and/or tail nodes are shorter than the short tip cut offs, then throw this away anyways
								auto tailNode = tail->tail_.lock();
								auto headNode = head->head_.lock();
								bool headShort = false;
								bool tailShort = false;
//								std::cout << __FILE__ << " " << __LINE__ << std::endl;
//								std::cout << n->nodeUid_ << " " << n->nameUid_ << std::endl;;
//								std::cout << n->k_.length() << std::endl;
//								std::cout << "pars.trimShortTips_ && tailNode->tailless() && tailNode->k_.size() - klen_ < pars.shortTipNumber_ && tailNode->inReadNamesIdx_.size() < pars.shortTipCutOff_: "
//										<< njh::colorBool(pars.trimShortTips_ && tailNode->tailless() && tailNode->k_.size() - klen_ < pars.shortTipNumber_ && tailNode->inReadNamesIdx_.size() < pars.shortTipCutOff_) << std::endl;
//								std::cout << "pars.trimShortTips_ &&  headNode->headless() && headNode->k_.size() - klen_ < pars.shortTipNumber_ && headNode->inReadNamesIdx_.size() < pars.shortTipCutOff_: "
//										<< njh::colorBool(pars.trimShortTips_ &&  headNode->headless() && headNode->k_.size() - klen_ < pars.shortTipNumber_ && headNode->inReadNamesIdx_.size() < pars.shortTipCutOff_) << std::endl;
//								std::cout << "pars.removeHeadlessTaillessAfterDisentaglement_ && tailNode->tailless() && tailNode->k_.length() <=pars.headlessTaillessCutOff_: "
//										<< njh::colorBool(pars.removeHeadlessTaillessAfterDisentaglement_ && tailNode->tailless() && tailNode->k_.length() <=pars.headlessTaillessCutOff_) << std::endl;
//								std::cout << "pars.removeHeadlessTaillessAfterDisentaglement_ && headNode->headless() && headNode->k_.length() <=pars.headlessTaillessCutOff_: "
//										<< njh::colorBool(pars.removeHeadlessTaillessAfterDisentaglement_ && headNode->headless() && headNode->k_.length() <=pars.headlessTaillessCutOff_) << std::endl;
//								std::cout << std::endl;
//								if(pars.trimShortTips_ && tailNode->tailless() && tailNode->k_.size() - klen_ < pars.shortTipNumber_ && tailNode->inReadNamesIdx_.size() < pars.shortTipCutOff_){
//									tailShort = true;
//									tailNode->on_ = false;
//								}
//								if(pars.trimShortTips_ &&  headNode->headless() && headNode->k_.size() - klen_ < pars.shortTipNumber_ && headNode->inReadNamesIdx_.size() < pars.shortTipCutOff_){
//									headShort= true;
//									headNode->on_ = false;
//								}
//								if(pars.removeHeadlessTaillessAfterDisentaglement_ && tailNode->tailless() && tailNode->k_.length() <=pars.headlessTaillessCutOff_){
//									tailShort = true;
//									//i'll leave the headless and tailless removal to it's own function rather than doing it here
//								}
//								if(pars.removeHeadlessTaillessAfterDisentaglement_ && headNode->headless() && headNode->k_.length() <=pars.headlessTaillessCutOff_){
//									headShort = true;
//									//i'll leave the headless and tailless removal to it's own function rather than doing it here
//								}
								//the above doesn't take into account the fact that some of the nodes that are going to be become headless or tailless are already disconected from other nodes
								//should be more precise about this but can't take the time right now njh 2020-03-26
								if(pars.trimShortTips_ &&  tailNode->k_.size() - klen_ < pars.shortTipNumber_ && tailNode->inReadNamesIdx_.size() < pars.shortTipCutOff_){
									tailShort = true;
									tailNode->on_ = false;
								}
								if(pars.trimShortTips_ &&  headNode->k_.size() - klen_ < pars.shortTipNumber_ && headNode->inReadNamesIdx_.size() < pars.shortTipCutOff_){
									headShort= true;
									headNode->on_ = false;
								}
								if(pars.removeHeadlessTaillessAfterDisentaglement_  && tailNode->k_.length() <=pars.headlessTaillessCutOff_){
									tailShort = true;
									//i'll leave the headless and tailless removal to it's own function rather than doing it here
								}
								if(pars.removeHeadlessTaillessAfterDisentaglement_  && headNode->k_.length() <=pars.headlessTaillessCutOff_){
									headShort = true;
									//i'll leave the headless and tailless removal to it's own function rather than doing it here
								}
								if (tailShort && headShort) {
									n->on_ = false;
								}
							}
						}
					} else if(1 == n->tailCount() && 1 != n->headCount()){
						//one tail
						std::set<uint32_t> readsEntering(
								n->inReadNamesIdx_.begin(),
								n->inReadNamesIdx_.end());
						//and to reach here tailEdges_ should only be one as well
						std::set<uint32_t> readsLeaving(
								keptTail.front()->inReadNamesIdx_.begin(),
								keptTail.front()->inReadNamesIdx_.end());
						std::vector<uint32_t> readsThatEnteredLeaving;
						std::set_intersection(
								readsEntering.begin(), readsEntering.end(),
								readsLeaving.begin(),  readsLeaving.end(),
								std::back_inserter(readsThatEnteredLeaving));
#if defined(PATHWEAVERSUPERDEBUG)
						{ //if (debug_) {// || true ) {
							std::cout << "Checking For tail count 1 uid: " << n->nodeUid_ << " " << n->nameUid_ << std::endl;;
							std::cout << "\tReadsEnteringNumber    : "
									<< readsEntering.size() << std::endl;
							std::cout << "\tReadsLeavingNumber     : "
									<< readsLeaving.size() << std::endl;
							std::cout << "\tReadsLeavingThatEntered: "
									<< readsThatEnteredLeaving.size() << std::endl;
							std::cout << njh::bashCT::red << "Breaking, only 1 tailed and no headed" << njh::bashCT::reset << std::endl;
							std::cout << std::endl;
						}
#endif
						//breaking for one tailed
						auto tail = n->getFirstOnTailEdge();
						tail->on_ = false;
						if(n->headless() && pars.throwAwayConservedAddNodesDuringDisentaglement_){
							//it's only here if we added this node already to something else,
							//if we are turning off it's only remaining heads and tails then we might as well turn of the node as well
							n->on_ = false;
						}else{
							//if the head and/or tail nodes are shorter than the short tip cut offs, then throw this away anyways
							auto tailNode = tail->tail_.lock();
							bool tailShort = false;
//							if(pars.trimShortTips_ && tailNode->tailless() && tailNode->k_.size() - klen_ < pars.shortTipNumber_ && tailNode->inReadNamesIdx_.size() < pars.shortTipCutOff_){
//								tailShort = true;
//								tailNode->on_ = false;
//							}
//							if(pars.removeHeadlessTaillessAfterDisentaglement_ && tailNode->tailless() && tailNode->k_.length() <=pars.headlessTaillessCutOff_){
//								tailShort = true;
//								//i'll leave the headless and tailless removal to it's own function rather than doing it here
//							}
							//the above doesn't take into account the fact that some of the nodes that are going to be become headless or tailless are already disconected from other nodes
							//should be more precise about this but can't take the time right now njh 2020-03-26
							if(pars.trimShortTips_  && tailNode->k_.size() - klen_ < pars.shortTipNumber_ && tailNode->inReadNamesIdx_.size() < pars.shortTipCutOff_){
								tailShort = true;
								tailNode->on_ = false;
							}
							if(pars.removeHeadlessTaillessAfterDisentaglement_ && tailNode->k_.length() <=pars.headlessTaillessCutOff_){
								tailShort = true;
								//i'll leave the headless and tailless removal to it's own function rather than doing it here
							}
							if (tailShort) {
								n->on_ = false;
							}
						}
					} else if (1 == n->headCount() && 1 != n->tailCount()) {
						//one head

						//next should only have a head count of 1 for it to be here
						std::set<uint32_t> readsEntering(
								keptHead.front()->inReadNamesIdx_.begin(),
								keptHead.front()->inReadNamesIdx_.end());
						//and to reach here tailEdges_ should only be one as well
						std::set<uint32_t> readsLeaving(
								n->inReadNamesIdx_.begin(),
								n->inReadNamesIdx_.end());
						std::vector<uint32_t> readsThatEnteredLeaving;
						std::set_intersection(
								readsEntering.begin(), readsEntering.end(),
								readsLeaving.begin(),  readsLeaving.end(),
								std::back_inserter(readsThatEnteredLeaving));
#if defined(PATHWEAVERSUPERDEBUG)
						{ //if (debug_) {// || true ) {
							std::cout << "Checking For head count 1 uid: " << n->nodeUid_ << " " << n->nameUid_ << std::endl;;
							std::cout << "\tReadsEnteringNumber    : "
									<< readsEntering.size() << std::endl;
							std::cout << "\tReadsLeavingNumber     : "
									<< readsLeaving.size() << std::endl;
							std::cout << "\tReadsLeavingThatEntered: "
									<< readsThatEnteredLeaving.size() << std::endl;
							std::cout << njh::bashCT::red << "Breaking, only 1 headed and no tailed" << njh::bashCT::reset << std::endl;
							std::cout << std::endl;
						}
#endif
						//breaking for one headed, not sufficient evidence that these two belong together, if we had bridging to a next node then we would have evidence but
						//without that we don't know if this is just connected because of low coverage in a multicopy sample
						auto head = n->getFirstOnHeadEdge();
						head->on_ = false;
						if(n->tailless() && pars.throwAwayConservedAddNodesDuringDisentaglement_){
							//it's only here if we added this node already to something else,
							//if we are turning off it's only remaining heads and tails then we might as well turn of the node as well
							n->on_ = false;
						}else{
							//if the head and/or tail nodes are shorter than the short tip cut offs, then throw this away anyways
							auto headNode = head->head_.lock();
							bool headShort = false;
//							if(pars.trimShortTips_ &&  headNode->headless() && headNode->k_.size() - klen_ < pars.shortTipNumber_ && headNode->inReadNamesIdx_.size() < pars.shortTipCutOff_){
//								headShort= true;
//								headNode->on_ = false;
//							}
//							if(pars.removeHeadlessTaillessAfterDisentaglement_ && headNode->headless() && headNode->k_.length() <=pars.headlessTaillessCutOff_){
//								headShort = true;
//								//i'll leave the headless and tailless removal to it's own function rather than doing it here
//							}
							//the above doesn't take into account the fact that some of the nodes that are going to be become headless or tailless are already disconected from other nodes
							//should be more precise about this but can't take the time right now njh 2020-03-26
							if(pars.trimShortTips_ && headNode->k_.size() - klen_ < pars.shortTipNumber_ && headNode->inReadNamesIdx_.size() < pars.shortTipCutOff_){
								headShort= true;
								headNode->on_ = false;
							}
							if(pars.removeHeadlessTaillessAfterDisentaglement_ && headNode->k_.length() <=pars.headlessTaillessCutOff_){
								headShort = true;
								//i'll leave the headless and tailless removal to it's own function rather than doing it here
							}
							if (headShort) {
								n->on_ = false;
							}
						}
					}
				} else{
					n->on_ = false;
				}
			}
		}
		//remove off nodes;
		removeOffNodes();
		//remove any connectors below cut off
		turnOffEdgesBelowCutOff(connectorCutOff);
		removeOffEdges();

	} //for (const auto & n : nodesToProcess)
	if(nodesSplit){
		resetNodePositions();
	}
//	if(!conservative){
//		exit(1);
//	}
	return nodesSplit;
}




} //namespace njhseq
