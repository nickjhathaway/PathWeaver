
// Created on 2015/01/13
// main.cpp

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


#include "PathWeaverPrograms.h"
#include <njhcpp/progutils/oneRing.hpp>


namespace njhseq {

class PathWeaverRunner: public njh::progutils::OneRing {
public:
	PathWeaverRunner();
};
PathWeaverRunner::PathWeaverRunner() :
		njh::progutils::OneRing(
				{
					addRing<WeaverRunner>(),
					addRing<CarmenServerRunner>(),
				},//
				{ }, "PathWeaver", "1", "0", "0-dev") {
}
} //namespace njhseq

int main(int argc, char* argv[]) {
	try {
		njhseq::PathWeaverRunner runner;
		return runner.run(argc, argv);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
