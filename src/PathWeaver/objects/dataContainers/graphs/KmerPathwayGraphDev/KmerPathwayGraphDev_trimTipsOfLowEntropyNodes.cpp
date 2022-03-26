/*
 * KmerPathwayGraphDev_trimTipsOfLowEntropyNodes.cpp
 *
 *  Created on: Mar 2, 2020
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


bool KmerPathwayGraphDev::trimEdgesNodeTipsWithLowEntropy(const readVecTrimmer::TrimEdgesByLowEntropyPars & pars){
	std::vector<std::shared_ptr<node>> nodesToProcess;
	for(const auto & n : nodes_){
			if(n->on_ && (n->tailless() || n->headless())){
			nodesToProcess.emplace_back(n);
		}
	}
	if(!nodesToProcess.empty()){
		bool modifiedNodes = false;
		bool trimmedNodes = false;
		for(const auto & n : nodesToProcess){

			auto trimPos = readVecTrimmer::determineTrimPostionsByLowEntropy(n->k_, pars);
			if(0 == trimPos.start_ && n->k_.size() == trimPos.end_){
				//nothing to trim
				continue;
			}
			if (n->tailless() && n->headless()) {
				if (trimPos.start_ < trimPos.end_) {
					if (trimPos.end_ - trimPos.start_ < klen_) {
						//it all failed, turn off the node
						n->on_ = false;
						n->turnOffAllEdges(); //probably don't need
						modifiedNodes = true;
					} else {
						trimmedNodes = true;
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << "n->k_.size()  : "<< n->k_.size() << std::endl;
							std::cout << "trimPos.start_: "<< trimPos.start_ << std::endl;
							std::cout << "trimPos.end_  : "<< trimPos.end_ << std::endl;
						}
#endif
						n->k_ = n->k_.substr(trimPos.start_, trimPos.end_); //trim the edges
					}
				} else {
					//it all failed, turn off the node
					n->on_ = false;
					n->turnOffAllEdges(); //probably don't need
					modifiedNodes = true;
				}
			} else if (n->tailless()) {
				if(trimPos.end_ < n->k_.size()){
					if(trimPos.end_ < klen_){
						//too much failed, turn off the node
						n->on_ = false;
						n->turnOffAllEdges();
						modifiedNodes = true;
					}else{
						trimmedNodes = true;
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << "n->k_.size()  : "<< n->k_.size() << std::endl;
							std::cout << "trimPos.start_: "<< trimPos.start_ << std::endl;
							std::cout << "trimPos.end_  : "<< trimPos.end_ << std::endl;
						}
#endif
						n->k_ = n->k_.substr(0, trimPos.end_); //trim off the back
					}
				} //else it's tail is fine
			} else if (n->headless()) {
				if(trimPos.start_ > 0){
					if(n->k_.size() - trimPos.start_ < klen_){
						//too much failed, turn off the node
						n->on_ = false;
						n->turnOffAllEdges();
						modifiedNodes = true;
					}else{
						trimmedNodes = true;
#if defined(PATHWEAVERSUPERDEBUG)
						{
							std::cout << __FILE__ << " " << __LINE__ << std::endl;
							std::cout << "n->k_.size()  : "<< n->k_.size() << std::endl;
							std::cout << "trimPos.start_: "<< trimPos.start_ << std::endl;
							std::cout << "trimPos.end_  : "<< trimPos.end_ << std::endl;
						}
#endif
						n->k_ = n->k_.substr(trimPos.start_);//trim off the front
					}
				}//else it's front is fine
			}
#if defined(PATHWEAVERSUPERDEBUG)
			{
				std::cout << std::endl;
			}
#endif
		}

		if(modifiedNodes){
			removeOffNodes();
			//removeOffEdges(); //remove off edges gets run in removeOffNodes
			return true;
		}
		if(trimmedNodes){
			return true;
		}
	}
	return false;
}





bool KmerPathwayGraphDev::trimTipsOfLowEntropyNodes(double entropyCutOff){
	std::vector<std::shared_ptr<node>> nodesToProcess;
	for(const auto & n : nodes_){
		if(n->on_ && (n->tailless() || n->headless())){
			nodesToProcess.emplace_back(n);
		}
	}
	if(!nodesToProcess.empty()){
		bool modifiedNodes = false;
		for(const auto & n : nodesToProcess){
			kmerInfo kInfo(n->k_, 2, false);
			if(kInfo.computeKmerEntropy() < entropyCutOff){
				n->on_ = false;
				n->turnOffAllEdges();
				modifiedNodes = true;
			}
		}
		if(modifiedNodes){
			removeOffNodes();
			//removeOffEdges(); //remove off edges gets run in removeOffNodes
			return true;
		}
	}

	return false;
}





}  // namespace njhseq
