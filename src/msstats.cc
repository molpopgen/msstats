/* 
   msstats - read data from ms via stdin, calculate common summary statistics

   Copyright (C) 2002-2007 Kevin Thornton

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  

*/

#include <iostream>
#include <vector>
#include <cstdio>
#include <numeric>
#include <utility>
#include <Sequence/SimParams.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/Recombination.hpp>
#include <otherstats.hpp>

using namespace std;
using namespace Sequence;

void calcstats(const SimData & d, const unsigned & mincount);

int main(int argc, char *argv[]) 
{
  std::vector<int> config;
  bool multipop = false;
  int mincount = 1;
  for(int arg = 1 ; arg < argc ; ++arg)
    {
      if( string(argv[arg]) == "-I" )
	{
	  multipop = true;
	  int npop = atoi(argv[++arg]);
	  for( int i=0;i<npop;++i )
	    {
	      config.push_back(atoi(argv[++arg]));
	    }
	}
      else if (string(argv[arg]) == "-m")
	{
	  mincount = atoi(argv[++arg]);
	}
    }
  int total = std::accumulate(config.begin(),config.end(),0,plus<int>());
  SimParams p;
  p.fromfile(stdin);
  SimData d;

#if __GNUG__ && __GNUC__ >= 3
  std::ios_base::sync_with_stdio(true);
#endif
  if(multipop)
    {
      std::cout << "rep\tpop\t";
    }
  std::cout << "S\t"
	    << "n1\t"
	    << "next\t"
	    << "theta\t"
	    << "pi\t"
	    << "thetaH\t"
	    << "Hprime\t"
	    << "tajd\t"
	    << "fulif\t"
	    << "fulid\t"
	    << "fulifs\t"
	    << "fulids\t"
	    << "rm\t"
	    << "rmmg\t"
	    << "nhaps\t"
	    << "hdiv\t"
	    << "wallb\t"
	    << "wallq\t"
	    << "rosasrf\t"
	    << "rosasru\t"
	    << "zns"
	    << endl;
  int rv;
  int rep=0;

  while( (rv=d.fromfile(stdin)) != EOF )
    {
      if(!multipop)
	{
	  calcstats(d,mincount);
	}
      else
	{
	  /*
	  if(d.size() != total)
	    {
	      std::cerr << "oh crap\n";
	      exit(10);
	    }
	  */
	  int sum = 0;
	  for(int i = 0 ; i < config.size() ; ++i)
	    {
	      SimData d2;
	      if(!d.empty())
		d2.assign(&*d.pbegin(),d.numsites(),
			  &d[sum],config[i]);
	      sum += config[i];
	      RemoveInvariantColumns(&d2);
	      cout << rep << '\t' << i << '\t';
	      calcstats(d2,mincount);
	    }
	}
      ++rep;
    } 
}

void calcstats(const SimData & d, const unsigned & mincount)
{
  PolySIM P(&d);

  cout << P.NumPoly()   << '\t' 
       << P.NumSingletons() << '\t'
       << P.NumExternalMutations() << '\t'
       << P.ThetaW()    << '\t' 
       << P.ThetaPi()   << '\t'
       << P.ThetaH()    << '\t' 
       << P.Hprime()    << '\t'
       << P.TajimasD()  << '\t'
       << P.FuLiF()     << '\t'
       << P.FuLiD()     << '\t'
       << P.FuLiFStar() << '\t'
       << P.FuLiDStar() << '\t';
  unsigned rm = P.Minrec(),nhaps=P.DandVK();
  cout << ((rm != SEQMAXUNSIGNED) ? rm : strtod("NAN",NULL)) << '\t'
       << Rm_MG(d,P.NumPoly(),nhaps) << '\t'
       << nhaps    << '\t'
       << P.DandVH()    << '\t'
       << P.WallsB()    << '\t'
       << P.WallsQ()    << '\t';
  pair<double,double> r2 = RozasR(d,P.ThetaPi(),P.NumPoly());
  cout << r2.first << '\t'
       << r2.second << '\t';
  if(d.numsites() > 1)
    {
      unsigned site1=0,site2=1;
      vector<double> LDSTATS(6);
      bool allskipped=true;
      double zns=0.;
      unsigned npairsLD=0;
      while( Recombination::Disequilibrium(&d,LDSTATS,&site1,&site2,false,0,mincount))
	{
	  if( !LDSTATS[LDSTATS.size()-1] )  //if site pair NOT skipped!
	    {
	      allskipped=false;
	      zns += LDSTATS[2];
	      ++npairsLD;
	    }
	}
      zns /= double(npairsLD);
      if( allskipped )
	{
	  cout << strtod("NAN",NULL) ;//<< "\n";
	}
      else
	{
	  cout << zns ;//<< '\n';
	}
	/*
	  vector< vector<double> > ld = P.Disequilibrium(mincount);
	  if(ld.empty()) //no sites made the mincount cutoff
	  {
	  cout << "nan\n";
	  }
	  else
	  {
	  double zns = 0.;
	  for(unsigned i=0;i<ld.size();++i)
	  {
	  zns += ld[i][2];
	  }
	  zns /= double(ld.size());
	  cout << zns ;
	  }
	*/
    }
  else 
    {
      cout << strtod("NAN",NULL);
    }
  cout << endl;
}
