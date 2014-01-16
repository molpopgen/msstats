#include <otherstats.hpp>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace Sequence;

std::pair<double,double> RozasR(const SimData & matrix, const double & thetapi, const unsigned & segsites ) 
/*
  Ramon-Osnin & Rozas' (2005, MBE 19: 2092) R2 and R2E statistics, returned in that order in a pair.
*/
{
  if(segsites == 0) return std::make_pair(strtod("NAN",NULL),strtod("NAN",NULL));
  unsigned int nsam = matrix.size();
  unsigned int sites = matrix.numsites();
  vector<short> derived(nsam,0),folded(nsam,0);
  unsigned i=0;
  for(SimData::const_site_iterator itr = matrix.sbegin() ; itr != matrix.send() ; ++itr,++i)
    {
      size_t c = count(itr->second.begin(),itr->second.end(),'1');
      if ( c == 1 ) //goes in to both the derived and folded vector
	{
	  string::size_type sst = itr->second.find('1');
	  derived[sst]++;
	  folded[sst]++;
	}
      else if( c == nsam-1 )
	{
	  string::size_type sst = itr->second.find('0');
	  folded[sst]++;
	}
    }
  double sum_f=0.,sum_u=0.; 
  for(unsigned i=0;i<nsam;++i)
    {
      sum_f += pow(( folded[i] - ( thetapi / 2. )), 2);   
      sum_u += pow(( derived[i] - ( thetapi / 2. )), 2);
    }
  return make_pair( sqrt(sum_f/double(nsam))/double(segsites), sqrt(sum_u/double(nsam))/double(segsites) );
}

unsigned Rm_MG( const Sequence::SimData & matrix, unsigned segsites, unsigned nhaps )
/*
  Returns a simple bound on the minumum number of recombination events, calculated according
  to equation 4 of Myers & Griffiths paper
*/
{
  bool anc = ( std::find(matrix.begin(),matrix.end(),std::string(segsites,'0')) != matrix.end() );
  return (nhaps > segsites) ? ( (anc) ? nhaps-segsites-1 : nhaps-segsites ) : 0;
}
