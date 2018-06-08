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
#include <sstream>
#include <vector>
#include <cstdio>
#include <numeric>
#include <utility>
#include <limits>
#include <Sequence/SimParams.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/Recombination.hpp>
#include <Sequence/FST.hpp>
#include <otherstats.hpp>
#include <tbb/parallel_invoke.h>
#include <tbb/task_scheduler_init.h>
using namespace std;
using namespace Sequence;

std::string calcstats(const SimData &d, const unsigned &mincount);
std::string process_input(const SimData &d, const std::vector<int> &config,
                          const int mincount, const bool multipop, int &rep);

int
main(int argc, char *argv[])
{
    std::vector<int> config;
    bool multipop = false;
    int mincount = 1;
    int nthreads = -1;
    for (int arg = 1; arg < argc; ++arg)
        {
            if (string(argv[arg]) == "-I")
                {
                    multipop = true;
                    int npop = atoi(argv[++arg]);
                    for (int i = 0; i < npop; ++i)
                        {
                            config.push_back(atoi(argv[++arg]));
                        }
                }
            else if (string(argv[arg]) == "-m")
                {
                    mincount = atoi(argv[++arg]);
                }
            else if (string(argv[arg]) == "-t")
                {
                    nthreads = atoi(argv[++arg]);
                }
        }
    tbb::task_scheduler_init init(nthreads);
    unsigned configSUM = accumulate(config.begin(), config.end(), 0u);
    SimParams p;
    p.fromfile(stdin);
    SimData d;

#if __GNUG__ && __GNUC__ >= 3
    std::ios_base::sync_with_stdio(true);
#endif
    if (multipop)
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
              << "zns";
    // Add header info for FST
    if (!config.empty())
        {
            cout << '\t';
            for (unsigned i = 0; i < config.size() - 1; ++i)
                {
                    for (unsigned j = i + 1; j < config.size(); ++j)
                        {
                            std::cout << "hsm" << i << j << '\t';
                        }
                }
        }
    std::cout << endl;

    int rv;
    int rep = 0;

    while ((rv = d.fromfile(stdin)) != EOF)
        {
            if (!config.empty() && d.size() != configSUM)
                {
                    cerr << "error: input sample size does not equal sum of "
                            "deme sizes\n";
                    exit(10);
                }
            std::cout << process_input(d, config, mincount, multipop, rep);
        }
}

std::string
calcstats(const SimData &d, const unsigned &mincount)
{
    PolySIM P(&d);
    std::ostringstream buffer;
    unsigned npoly, nsingle, nexternal, rm, nhaps;
    double thetaw, thetapi, thetah, hprime, tajd, fulif, fulid, fulifstar,
        fulidstar, dvh, wallb, wallq;
    tbb::parallel_invoke(
        [&npoly, &P] { npoly = P.NumPoly(); },
        [&nsingle, &P] { nsingle = P.NumSingletons(); },
        [&nexternal, &P] { nexternal = P.NumExternalMutations(); },
        [&thetapi, &P] { thetapi = P.ThetaPi(); },
        [&thetaw, &P] { thetaw = P.ThetaW(); },
        [&thetah, &P] { thetah = P.ThetaH(); },
        [&hprime, &P] { hprime = P.Hprime(); },
        [&tajd, &P] { tajd = P.TajimasD(); },
        [&fulid, &P] { fulid = P.FuLiD(); });
    tbb::parallel_invoke(
        [&fulif, &P] { fulif = P.FuLiF(); },
        [&fulidstar, &P] { fulidstar = P.FuLiDStar(); },
        [&fulifstar, &P] { fulifstar = P.FuLiFStar(); },
        [&rm, &P] { rm = P.Minrec(); }, [&nhaps, &P] { nhaps = P.DandVK(); },
        [&dvh, &P] { dvh = P.DandVH(); }, [&wallb, &P] { wallb = P.WallsB(); },
        [&wallq, &P] { wallq = P.WallsQ(); });
    buffer << npoly << '\t' << nsingle << '\t' << nexternal << '\t' << thetaw
           << '\t' << thetapi << '\t' << thetah << '\t' << hprime << '\t'
           << tajd << '\t' << fulif << '\t' << fulid << '\t' << fulifstar
           << '\t' << fulidstar << '\t';
    buffer << ((rm != SEQMAXUNSIGNED)
                   ? rm
                   : std::numeric_limits<double>::quiet_NaN())
           << '\t' << Rm_MG(d, P.NumPoly(), nhaps) << '\t' << nhaps << '\t'
           << dvh << '\t' << wallb << '\t' << wallq << '\t';
    pair<double, double> r2 = RozasR(d, P.ThetaPi(), P.NumPoly());
    buffer << r2.first << '\t' << r2.second << '\t';
    if (d.numsites() > 1)
        {
            /*
                        unsigned site1 = 0, site2 = 1;
            vector<double> LDSTATS(6);
            bool allskipped = true;
            double zns = 0.;
            unsigned npairsLD = 0;
            while (Recombination::Disequilibrium(&d, LDSTATS, &site1, &site2,
                                                 false, 0, mincount))
                {
                    if (!LDSTATS[LDSTATS.size()
                                 - 1]) // if site pair NOT skipped!
                        {
                            allskipped = false;
                            zns += LDSTATS[2];
                            ++npairsLD;
                        }
                }
                        */
            double zns = 0.;
            unsigned npairsLD = 0;
            bool allskipped = true;
            auto LD = P.Disequilibrium(mincount);
            for (auto &&LDi : LD)
                {
                    if (!LDi.skipped) // not skipped
                        {
                            allskipped = false;
                            zns += LDi.rsq;
                            ++npairsLD;
                        }
                }
            zns /= static_cast<double>(npairsLD);
            if (allskipped)
                {
                    buffer << std::numeric_limits<double>::quiet_NaN(); //<< "\n";
                }
            else
                {
                    buffer << zns; //<< '\n';
                }
        }
    else
        {
            buffer << strtod("NAN", NULL);
        }
    return buffer.str();
}

std::string
process_input(const SimData &d, const std::vector<int> &config,
              const int mincount, const bool multipop, int &rep)
{
    std::ostringstream buffer;
    if (!multipop)
        {
            buffer << calcstats(d, mincount);
            buffer << endl;
        }
    else
        {
            int sum = 0;

            // Calculate FST
            ostringstream fstout;
            if (!config.empty())
                {
                    for (unsigned i = 0; i < config.size() - 1; ++i)
                        {
                            for (unsigned j = i + 1; j < config.size(); ++j)
                                {
                                    std::vector<polymorphicSite> demeData;
                                    unsigned offset1 = std::accumulate(
                                        config.begin(), config.begin() + i,
                                        0u);
                                    unsigned offset2 = std::accumulate(
                                        config.begin(), config.begin() + j,
                                        0u);

                                    for (SimData::const_site_iterator sitr
                                         = d.sbegin();
                                         sitr != d.send(); ++sitr)
                                        {
                                            string data(
                                                string(sitr->second.begin()
                                                           + offset1,
                                                       sitr->second.begin()
                                                           + offset1
                                                           + config[i])
                                                + string(sitr->second.begin()
                                                             + offset2,
                                                         sitr->second.begin()
                                                             + offset2
                                                             + config[j]));
                                            demeData.push_back(
                                                make_pair(sitr->first, data));
                                        }
                                    SimData demeij(demeData.begin(),
                                                   demeData.end());
                                    unsigned configij[2];
                                    configij[0] = config[i];
                                    configij[1] = config[j];
                                    FST fst(&demeij, 2, configij);
                                    fstout << fst.HSM() << '\t';
                                }
                        }
                }
            for (std::size_t i = 0; i < config.size(); ++i)
                {
                    SimData d2;
                    if (!d.empty())
                        {
                            d2.assign(d.GetPositions(),
                                      std::vector<string>(d.begin() + sum,
                                                          d.begin() + sum
                                                              + config[i]));
                        }
                    removeInvariantPos(d2);
                    buffer << rep << '\t' << i << '\t';
                    buffer << calcstats(d2, mincount);
                    if (!fstout.str().empty())
                        {
                            buffer << '\t' << fstout.str();
                        }
                    buffer << endl;
                }
        }
    ++rep;
    return buffer.str();
}
