/*
 * ContigsCompareGraph.cpp
 *
 *  Created on: Nov 5, 2019
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


#include "ContigsCompareGraph.hpp"


namespace njhseq {

ContigsCompareGraph::node::node(const std::string & k, uint32_t cnt, uint32_t kLen) :
		k_(k),uid_(k), cnt_(cnt), kLen_(kLen) {
}
void ContigsCompareGraph::node::resetVisitCount(){
	visitCount_ = 0;
}
void ContigsCompareGraph::node::addHead(const std::shared_ptr<edge> & e){
	headEdges_.push_back(e);
}
void ContigsCompareGraph::node::addTail(const std::shared_ptr<edge> & e){
	tailEdges_.push_back(e);
}
bool ContigsCompareGraph::node::headless() const {
	return 0 == headCount();
}
uint32_t ContigsCompareGraph::node::headCount() const {
	uint32_t headCount = 0;
	for (const auto & head : headEdges_) {
		if (head->on_) {
			++headCount;
		}
	}
	return headCount;
}
std::shared_ptr<ContigsCompareGraph::edge> ContigsCompareGraph::node::getFirstOnHeadEdge() const {
	for (const auto & head : headEdges_) {
		if (head->on_) {
			return head;
		}
	}
	return nullptr;
}
std::shared_ptr<ContigsCompareGraph::edge> ContigsCompareGraph::node::getLastOnHeadEdge() const {
	for (const auto & head : iter::reversed(headEdges_)) {
		if (head->on_) {
			return head;
		}
	}
	return nullptr;
}
bool ContigsCompareGraph::node::tailless() const{
	return 0 == tailCount();
}
uint32_t ContigsCompareGraph::node::tailCount() const{
	uint32_t tailCount = 0;
	for(const auto & tail : tailEdges_){
		if(tail->on_){
			++tailCount;
		}
	}
	return tailCount;
}
std::shared_ptr<ContigsCompareGraph::edge> ContigsCompareGraph::node::getFirstOnTailEdge() const {
	for (const auto & tail : tailEdges_) {
		if (tail->on_) {
			return tail;
		}
	}
	return nullptr;
}
std::shared_ptr<ContigsCompareGraph::edge> ContigsCompareGraph::node::getLastOnTailEdge() const {
	for (const auto & tail : iter::reversed(tailEdges_)) {
		if (tail->on_) {
			return tail;
		}
	}
	return nullptr;
}
std::shared_ptr<ContigsCompareGraph::node> ContigsCompareGraph::node::getHeadNode(const std::string & uid)const{
	std::shared_ptr<node> ret = nullptr;
	for(const auto & head : headEdges_){
		if(uid == head->head_.lock()->uid_){
			ret = head->head_.lock();
			break;
		}
	}
	return ret;
}
std::shared_ptr<ContigsCompareGraph::node> ContigsCompareGraph::node::getTailNode(const std::string & uid)const{
	std::shared_ptr<node> ret = nullptr;
	for(const auto & tail : tailEdges_){
		if(uid == tail->tail_.lock()->uid_){
			ret = tail->tail_.lock();
			break;
		}
	}
	return ret;
}


} // namespace njhseq
