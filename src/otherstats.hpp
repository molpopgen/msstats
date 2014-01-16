#ifndef __OTHERSTATS_HPP__
#define __OTHERSTATS_HPP__
#include <Sequence/SimData.hpp>
#include <utility>

std::pair<double,double> RozasR(const Sequence::SimData & matrix, 
				const double & thetapi,
				const unsigned & segsites );
unsigned Rm_MG( const Sequence::SimData & matrix, unsigned segsites, unsigned nhaps );
#endif
