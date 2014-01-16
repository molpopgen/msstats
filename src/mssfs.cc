#include <Sequence/SimData.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace Sequence;

int main(int argc, char ** argv)
{
  int argn = 1;

  vector<unsigned> sfs;

  SimData d;
  int rv;
  unsigned nreps=0;
  while( (rv = d.fromfile(stdin)) != EOF )
    {
      ++nreps;
      if( sfs.empty() && ! d.empty() )
	{
	  sfs = vector<unsigned>(d.size()-1,0u);
	}
      for(SimData::const_site_iterator i = d.sbegin() ; 
	  i != d.send() ; ++i)
	{
	  unsigned c = count(i->second.begin(),
			     i->second.end(),'1');
	  sfs[c-1]++;
	}
    }
  for(unsigned i=0;i<sfs.size();++i)
    {
      cout << double(sfs[i])/double(nreps) << endl;
    }
}
