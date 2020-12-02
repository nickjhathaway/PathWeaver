/*
 * KmerPathwayGraphDev_copyGraph.cpp
 *
 *  Created on: Jan 8, 2019
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




KmerPathwayGraphDev KmerPathwayGraphDev::copyGraph() const{
	KmerPathwayGraphDev outGraph(klen_, occurenceCutOff_);
//	outGraph.debug_ = debug_;
//	outGraph.verbose_ = verbose_;

	outGraph.allowableErrorForHPIndexCollapse_ = allowableErrorForHPIndexCollapse_;
	outGraph.homopolymerIndelCollapseFreqMultiplier_ = homopolymerIndelCollapseFreqMultiplier_;
	outGraph.bridgingCutOff_ = bridgingCutOff_;

	outGraph.numThreads_ = numThreads_;

	outGraph.kCounts_ = kCounts_;
	outGraph.nodePositions_ = nodePositions_;
	outGraph.readNames_ = readNames_;
	//outGraph.throwAwayConservedAddNodesDuringDisentaglement_ = throwAwayConservedAddNodesDuringDisentaglement_;

	outGraph.nodes_ = std::vector<std::shared_ptr<node>>(nodes_.size());
	std::vector<uint64_t> nodePositions(nodes_.size());
	njh::iota<uint64_t>(nodePositions, 0);
	njh::concurrent::LockableQueue<uint64_t> positionsQueue(nodePositions);
	std::function<void()> addNodes = [this,&outGraph,&positionsQueue](){
		uint64_t pos = std::numeric_limits<uint64_t>::max();
		while(positionsQueue.getVal(pos)){
			const auto & n = nodes_[pos];
			auto cpNode = std::make_shared<node>(n->k_, n->cnt_, n->kLen_);
			cpNode->group_ = n->group_;
			cpNode->inReadNamesIdx_ = n->inReadNamesIdx_;
			cpNode->on_ = n->on_;
			cpNode->nodeUid_ = n->nodeUid_;
			cpNode->nameUid_ = n->nameUid_;
			cpNode->visitCount_ = n->visitCount_;
			outGraph.nodes_[pos] = cpNode;
		}
	};
	njh::concurrent::runVoidFunctionThreaded(addNodes, numThreads_);

//	//first copy nodes without the edges
//	for(const auto & n : nodes_){
//		auto cpNode = std::make_shared<node>(n->k_, n->cnt_, n->kLen_);
//		cpNode->group_ = n->group_;
//		cpNode->inReadNamesIdx_ = n->inReadNamesIdx_;
//		cpNode->on_ = n->on_;
//		cpNode->uid_ = n->uid_;
//		cpNode->visitCount_ = n->visitCount_;
//		outGraph.nodes_.emplace_back(cpNode);
//	}

	//now copy in edges

	for(const auto nodePos : iter::range(nodes_.size())){
		//copy in heads

		for(const auto & head : nodes_[nodePos]->headEdges_){
			bool found = false;
			for(const auto & otherHead : outGraph.nodes_[nodePos]->headEdges_){
				if(otherHead->createUidWithNameUid() == head->createUidWithNameUid()){
					found = true;
					break;
				}
			}
			if(!found){
				auto headNodePos = head->head_.lock()->nodeUid_;
				auto addingHead = std::make_shared<KmerPathwayGraphDev::edge>(
						outGraph.nodes_[headNodePos],
						outGraph.nodes_[head->tail_.lock()->nodeUid_],
						head->cnt_,
						head->inReadNamesIdx_);
				addingHead->on_ = head->on_;
				//add
				outGraph.nodes_[nodePos]->addHead(addingHead);
				//head node means this node is the tail, so add to the head node to it's tail nodes
				outGraph.nodes_[headNodePos]->addTail(addingHead);
			}
		}
		//copy in tails
		for(const auto & tail : nodes_[nodePos]->tailEdges_){
			bool found = false;
			for(const auto & otherTail : outGraph.nodes_[nodePos]->tailEdges_){
				if(otherTail->createUidWithNameUid() == tail->createUidWithNameUid()){
					found = true;
					break;
				}
			}
			if(!found){
				auto tailNodePos = tail->tail_.lock()->nodeUid_;
				auto addingTail = std::make_shared<KmerPathwayGraphDev::edge>(
						outGraph.nodes_[tail->head_.lock()->nodeUid_],
						outGraph.nodes_[tailNodePos],
						tail->cnt_,
						tail->inReadNamesIdx_);
				addingTail->on_ = tail->on_;
				//add
				outGraph.nodes_[nodePos]->addTail(addingTail);
				//tail node means this node is the head, so add to the tail node to it's head nodes
				outGraph.nodes_[tailNodePos]->addHead(addingTail);
			}
		}
	}

	return outGraph;
}


KmerPathwayGraphDev::PossibleHap::PossibleHap():seqBase_("", ""){

}

KmerPathwayGraphDev::PossibleHap::PossibleHap(const node & firstNode):seqBase_("", ""){
	initWithFirstContig(firstNode);
}



KmerPathwayGraphDev::PossibleHap::PossibleHap(const seqInfo & seqBase,
		const Bed6RecordCore & firstContigLoc):
		seqBase_(seqBase), contigLocations_({firstContigLoc}){

}

void KmerPathwayGraphDev::PossibleHap::initWithFirstContig(const node & nextNode){
	std::string addingSeq = nextNode.k_;
	contigLocations_.emplace_back(seqBase_.name_,
			seqBase_.seq_.size(),
			seqBase_.seq_.size() + addingSeq.size(),
			nextNode.nameUid_,
			addingSeq.size(),
			'+');
	seqBase_.append(addingSeq);
}

void KmerPathwayGraphDev::PossibleHap::addContig(const node & nextNode, const std::string & edgeUID){
	std::string addingSeq = nextNode.k_;
	if("" != seqBase_.seq_){
		//already added before, this means we need to take only suffix of the next node
		addingSeq = nextNode.k_.substr(nextNode.kLen_ - 1);
	}
	contigLocations_.emplace_back(seqBase_.name_,
			seqBase_.seq_.size(),
			seqBase_.seq_.size() + addingSeq.size(),
			nextNode.nameUid_,
			addingSeq.size(),
			'+');
	seqBase_.append(addingSeq);

	edgesUids_.emplace(edgeUID);
}

void growPossibleHap(const KmerPathwayGraphDev::node & node,
		std::vector<KmerPathwayGraphDev::PossibleHap> currentHaps,
		std::vector<KmerPathwayGraphDev::PossibleHap> & finalHaps){
	if(node.tailless()){
		addOtherVec(finalHaps, currentHaps);
	}else{
		for(auto & currentHap : currentHaps){
			bool didAdd = false;
			for(const auto & tail : node.tailEdges_){
				if(tail->on_ && !njh::in(tail->createUidWithNameUid(), currentHap.edgesUids_)){
					didAdd = true;
					auto nextNode = tail->tail_.lock();
					currentHap.addContig(*nextNode, tail->createUidWithNameUid());
					growPossibleHap(*nextNode, currentHaps, finalHaps);
				}
			}
			if(!didAdd){
				//if ends in a loop, didn't grow anymore, need to add
				addOtherVec(finalHaps, currentHaps);
			}
		}
	}
}

std::vector<KmerPathwayGraphDev::PossibleHap> KmerPathwayGraphDev::generateAllPossibleHaps() const{
	std::vector<KmerPathwayGraphDev::PossibleHap> ret;
	std::vector<std::shared_ptr<node>> headlessNodes;
	for(const auto & n : nodes_){
		if(n->on_ && n->headless()){
			headlessNodes.push_back(n);
		}
	}

	for(const auto & headless : headlessNodes){
		std::vector<KmerPathwayGraphDev::PossibleHap> currentHaps{KmerPathwayGraphDev::PossibleHap(*headless)};
		growPossibleHap(*headless, currentHaps, ret);
	}

	return ret;
}


void KmerPathwayGraphDev::turnOffLowestEdges(){
	for(const auto & n : nodes_){
		if(n->tailCount() > 1){
			uint32_t maxTailReadCount = 0;
			for(const auto & t : n->tailEdges_){
				if(t->on_){
					if(t->inReadNamesIdx_.size() > maxTailReadCount){
						maxTailReadCount = t->inReadNamesIdx_.size();
					}
				}
			}
			for(const auto & t : n->tailEdges_){
				if(t->on_ && t->inReadNamesIdx_.size() < maxTailReadCount){
					t->on_ = false;
				}
			}
		}
		if(n->headCount() > 1){
			uint32_t maxHeadReadCount = 0;
			for(const auto & h : n->headEdges_){
				if(h->on_){
					if(h->inReadNamesIdx_.size() > maxHeadReadCount){
						maxHeadReadCount = h->inReadNamesIdx_.size();
					}
				}
			}
			for(const auto & h : n->headEdges_){
				if(h->on_ && h->inReadNamesIdx_.size() < maxHeadReadCount){
					h->on_ = false;
				}
			}
		}
	}
}



}  // namespace njhseq
