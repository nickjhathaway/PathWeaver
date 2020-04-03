

/*
 * TandemRepeatUtils.cpp
 *
 *  Created on: Oct 17, 2019
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


#include "TandemRepeatUtils.hpp"


namespace njhseq {


StartStopTanPos::StartStopTanPos(uint32_t start, uint32_t stop) :
		start_(start), stop_(stop){
}

uint32_t StartStopTanPos::length() const {
	return stop_ - start_;
}



TandemRepeatPlusAdjustedSize::TandemRepeatPlusAdjustedSize(const TandemRepeat & repeat,
		uint32_t adjustedPartialSize) :
		repeat_(repeat), adjustedPartialSize_(adjustedPartialSize) {
}


std::vector<TandemRepeatPlusAdjustedSize> readInTandemsWithPartInfo(const bfs::path & tandemInfoFnp){
	std::vector<TandemRepeatPlusAdjustedSize> allTandems;
	if(bfs::exists(tandemInfoFnp)){
		TableReader tandemInfoReader(TableIOOpts::genTabFileIn(tandemInfoFnp, true));
		VecStr row;
		while(tandemInfoReader.getNextRow(row)){
			//repeat	size	adjustedPartialSize
			auto repeat = row[tandemInfoReader.header_.getColPos("repeat")];
			uint32_t size = njh::StrToNumConverter::stoToNum<uint32_t>(row[tandemInfoReader.header_.getColPos("size")]);
			uint32_t adjustedPartialSize = njh::StrToNumConverter::stoToNum<uint32_t>(row[tandemInfoReader.header_.getColPos("adjustedPartialSize")]);
			allTandems.emplace_back(TandemRepeat(repeat, size/repeat.size(), size * 2, 0, size), adjustedPartialSize);
		}
	}
	return allTandems;
}


void ProcessedDeterminedTandemRepeatsResults::writeOutInfos(const bfs::path & currentKmerDirectory, bool overWrite){
	OutOptions finalTandemsOutOpts(njh::files::make_path(currentKmerDirectory, "finalTandemSeqs.txt"));
	finalTandemsOutOpts.overWriteFile_ = overWrite;
	OutOptions finalTandemsCountOutOpts(njh::files::make_path(currentKmerDirectory, "finalTandemSeqsCounts.txt"));
	finalTandemsCountOutOpts.overWriteFile_ = overWrite;
	OutOptions finalTandemsAltOutOpts(njh::files::make_path(currentKmerDirectory, "finalTandemSeqsAlts.txt"));
	finalTandemsAltOutOpts.overWriteFile_ = overWrite;

	{
		OutputStream finalTandemsOut(finalTandemsOutOpts);
		OutputStream finalTandemsCountsOut(finalTandemsCountOutOpts);
		finalTandemsCountsOut << "tandem\tcount\n" << std::endl;
		njh::sort(finalTandemSeqs);
		finalTandemsOut << njh::conToStr(finalTandemSeqs, "\n") << std::endl;
		OutputStream finalTandemsAltsOut(finalTandemsAltOutOpts);
		for(const auto & final : finalTandemSeqs){
			finalTandemsCountsOut << final << '\t' << tandemsCounts[final] << std::endl;
			finalTandemsAltsOut << final << "\t" << njh::conToStr(tandemsAlts[final], ",") << std::endl;
		}
	}
}



ProcessedDeterminedTandemRepeatsResults processGlobalTandems(
		const std::vector<TandemRepeatPlusAdjustedSize> & allTandems,
		uint32_t currentKLen,
		const HaploPathFinder::PathFinderCorePars & extractionPars){


	auto tandemOccurenceCutOff = extractionPars.finalTandemOccurenceCutOff_;
	if(extractionPars.finalTandemOccurenceCutOff_ >=6){
		tandemOccurenceCutOff = extractionPars.finalTandemOccurenceCutOff_/2;
	}

	ProcessedDeterminedTandemRepeatsResults ret;
	//		uint32_t maxTandemRepeatFullSize = 0;
	//		TandemRepeatPlusAdjustedSize maxTandemRepeat{TandemRepeat("", 0, 0, 0, 0	), 0};
	for(const auto & t : allTandems){
//			if(t.adjustedPartialSize_ < currentKLen){
//				continue;
//			}
//			std::cout << "t.repeat_.getSize(): " << t.repeat_.getSize() << std::endl;
//			std::cout << "currentKLen: " << currentKLen << std::endl;
		if(t.repeat_.getSize() < currentKLen){
			continue;
		}
		++ret.tandemsCounts[t.repeat_.repeat_];
//			if(t.adjustedPartialSize_ > maxTandemRepeatFullSize){
//				maxTandemRepeatFullSize = t.adjustedPartialSize_;
//				maxTandemRepeat = t;
//			}
	}
//		exit(1);
	VecStr tandemSeqs = njh::getVecOfMapKeys(ret.tandemsCounts);
	njh::sort(tandemSeqs, [&ret](const std::string & str1, const std::string & str2){
		if(ret.tandemsCounts[str1] == ret.tandemsCounts[str2]){
			return str1 > str2;
		}else{
			return ret.tandemsCounts[str1] > ret.tandemsCounts[str2];
		}
	});
	if (extractionPars.debug) {
		std::cout << "Tandem Counts" << std::endl;
		for (const auto & t : tandemSeqs) {
			std::cout << t << ": " << ret.tandemsCounts[t] << std::endl;
		}
	}
//	for (const auto & t : tandemSeqs) {
//		std::cout << t << ": " << tandemsCounts[t] << std::endl;
//	}
	//uint32_t finalTandemOccurenceCutOff = 2;
	//uint32_t finalTandemOccurenceCutOff = extractionPars.kOccurenceCutOff;
	for(const auto & t : tandemSeqs){
		bool alreadyHave = false;
		for(const auto & otherT : ret.finalTandemSeqs){
			if(checkTwoRotatingStrings(t, otherT, 0).size() > 0){
				alreadyHave = true;
				ret.tandemsCountsAdjusted[otherT] += ret.tandemsCounts[t];
				break;
			}
		}
//		if(!alreadyHave && ret.tandemsCounts[t] > tandemOccurenceCutOff){
//			ret.finalTandemSeqs.emplace_back(t);
//			ret.tandemsCountsAdjusted[t] += ret.tandemsCounts[t];
//		}
		if(!alreadyHave){
			ret.finalTandemSeqs.emplace_back(t);
			ret.tandemsCountsAdjusted[t] += ret.tandemsCounts[t];
		}
	}
	//filter final tandems based on final adjusted tandems
	if (extractionPars.debug) {
		std::cout << std::endl;
		std::cout << "final " << std::endl;
		for(const auto & t : ret.finalTandemSeqs){
			std::cout << t << ": " << ret.tandemsCounts[t] << std::endl;
		}
		std::cout << "final adjusted counts " << std::endl;
		for(const auto & t : ret.finalTandemSeqs){
			std::cout << t << ": " << ret.tandemsCountsAdjusted[t] << std::endl;
		}
	}
	if(!extractionPars.collapseFinalDeterminedTandems_){
		ret.finalTandemSeqs.clear();
		for(const auto & tan : ret.tandemsCountsAdjusted){
			if(tan.second > tandemOccurenceCutOff){
				ret.finalTandemSeqs.emplace_back(tan.first);
			}
		}
	}
	if (extractionPars.debug) {
		std::cout << "After occurence filter"<< std::endl;
		std::cout << "final " << std::endl;
		for(const auto & t : ret.finalTandemSeqs){
			std::cout << t << ": " << ret.tandemsCounts[t] << std::endl;
		}
		std::cout << "final adjusted counts " << std::endl;
		for(const auto & t : ret.finalTandemSeqs){
			std::cout << t << ": " << ret.tandemsCountsAdjusted[t] << std::endl;
		}
	}

	if(extractionPars.collapseFinalDeterminedTandems_){
		VecStr toBeFinalTandems;
		//sort so the smaller tandmes are first, then we determine if the smaller tandems are simply multiples of the same tandem
		//system is not perfect right now, will mostly work best for di repeats more work should be done NJH 2018/12/17
		std::unordered_map<std::string, std::set<std::string>> tandemsAlts;
		for(const auto & t : allTandems){
			if (t.repeat_.getSize() < currentKLen) {
				continue;
			}
			for(const auto & finalT : ret.finalTandemSeqs){
				if(checkTwoRotatingStrings(t.repeat_.repeat_, finalT, 0).size() > 0){
					if(t.repeat_.repeat_ != finalT){
						tandemsAlts[finalT].emplace(t.repeat_.repeat_);
					}
					break;
				}
			}
		}
		njh::sort(ret.finalTandemSeqs, [&ret](const std::string & str1, const std::string & str2){
			if(str1.size() == str2.size()){
				if(ret.tandemsCounts[str1] == ret.tandemsCounts[str2]){
					return str1 < str2;
				}else{
					return ret.tandemsCounts[str1] > ret.tandemsCounts[str2];
				}
			}else{
				return str1.size() < str2.size();
			}
		});
		for(const auto & tand : ret.finalTandemSeqs){
			bool add = true;
			for(const auto & otherTands : toBeFinalTandems){
				auto tr = aligner::findTandemRepeatOfStrInSequence(tand, otherTands, 2, -2, -7, 2 * (tand.size()/otherTands.size() ) );
				//tr.outPutInfoFormated(std::cout, tand);
				if(tr.getSize() == tand.size()){
					add = false;
					ret.tandemsCountsAdjustedAfterCollapse[otherTands] += ret.tandemsCountsAdjusted[tand];
					break;
				}
				//also check alternative tandems;
				if(add){
					for(const auto & altOtherTands : tandemsAlts[otherTands]){
						auto tr = aligner::findTandemRepeatOfStrInSequence(tand, altOtherTands, 2, -2, -7, 2 * (tand.size()/altOtherTands.size() ) );
						//tr.outPutInfoFormated(std::cout, tand);
						if(tr.getSize() == tand.size()){
							add = false;
							ret.tandemsCountsAdjustedAfterCollapse[otherTands] += ret.tandemsCountsAdjusted[tand];
							break;
						}
					}
				}
				if(!add){
					break;
				}
			}
			if(add){
				toBeFinalTandems.emplace_back(tand);
				ret.tandemsCountsAdjustedAfterCollapse[tand] += ret.tandemsCountsAdjusted[tand];
			}
		}
		//ret.finalTandemSeqs = toBeFinalTandems;
		ret.finalTandemSeqs.clear();
		for(const auto & tan : 	ret.tandemsCountsAdjustedAfterCollapse ){
			if(tan.second > tandemOccurenceCutOff){
				ret.finalTandemSeqs.emplace_back(tan.first);
			}
		}
		njh::sort(ret.finalTandemSeqs, [&ret](const std::string & str1, const std::string & str2){
			if(str1.size() == str2.size()){
				if(ret.tandemsCountsAdjustedAfterCollapse[str1] == ret.tandemsCountsAdjustedAfterCollapse[str2]){
					return str1 < str2;
				}else{
					return ret.tandemsCountsAdjustedAfterCollapse[str1] > ret.tandemsCountsAdjustedAfterCollapse[str2];
				}
			}else{
				return str1.size() < str2.size();
			}
		});
		if (extractionPars.debug) {
			std::cout << std::endl;
			std::cout << "final after collapse " << std::endl;
			for(const auto & t : ret.finalTandemSeqs){
				std::cout << t << ": " << ret.tandemsCounts[t] << std::endl;
			}

			std::cout << "final after collapse adjusted " << std::endl;
			for(const auto & t : ret.finalTandemSeqs){
				std::cout << t << ": " << ret.tandemsCountsAdjustedAfterCollapse[t] << std::endl;
			}
		}
		//exit(1);
	}

	for(const auto & t : allTandems){
		//			if(t.adjustedPartialSize_ < currentKLen){
		//				continue;
		//			}
		if (t.repeat_.getSize() < currentKLen) {
			continue;
		}
		for(const auto & finalT : ret.finalTandemSeqs){
			if(checkTwoRotatingStrings(t.repeat_.repeat_, finalT, 0).size() > 0){
				if(t.repeat_.repeat_ != finalT){
					ret.tandemsAlts[finalT].emplace(t.repeat_.repeat_);
				}
				if(t.adjustedPartialSize_ > ret.maxTandemSizesForEachTandem[finalT] ){
//						std::cout << finalT << " " << t.repeat_.getSize() << " " << t.adjustedPartialSize_ << std::endl;
//						std::cout << std::endl;
					ret.maxTandemSizesForEachTandem[finalT] = t.adjustedPartialSize_;
					ret.maxTandemEachTandem[finalT] = std::make_shared<TandemRepeatPlusAdjustedSize>(t);
				}
				break;
			}
		}
	}


	for(const auto & t : ret.finalTandemSeqs){
		ret.tandemMots.emplace(t, motif(t));
		for(const auto & altT : ret.tandemsAlts[t]){
			ret.tandemsAltMots[t].emplace(altT, motif(altT));
		}
	}


//	std::cout << "alts: " << std::endl;
//	for(const auto & tands : tandemsAlts){
//		std::cout << tands.first << std::endl;
//		std::cout << "\t" << njh::conToStr(tands.second, ",") << std::endl;
//	}
//	for(const auto & max : maxTandemEachTandem){
//		std::cout << max.first << std::endl;
//		std::cout << "\t" << "max.second->startPos_           : " << max.second->repeat_.startPos_ << std::endl;
//		std::cout << "\t" << "max.second->stopPos_            : " << max.second->repeat_.stopPos_ << std::endl;
//		std::cout << "\t" << "max.second->alignScore_         : " << max.second->repeat_.alignScore_ << std::endl;
//		std::cout << "\t" << "max.second->numberOfRepeats_    : " << max.second->repeat_.numberOfRepeats_ << std::endl;
//		std::cout << "\t" << "max.second->repeat_             : " << max.second->repeat_.repeat_ << std::endl;
//		std::cout << "\t" << "max.second->getSize()           : " << max.second->repeat_.getSize() << std::endl;
//		std::cout << "\t" << "max.second->adjustedPartialSize_: " << max.second->adjustedPartialSize_ << std::endl;
//	}
	return ret;
}





}  // namespace njhseq
