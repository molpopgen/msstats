	msstats - read data from ms via stdin, calculate common summary statistics



  Copyright (C) 2002 Kevin Thornton

  msstats is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Comments are welcome.

	- Kevin Thornton <krthornt@uci.edu>

#Usage

[coalesent simulation] | msstats | gzip > stats.gz

"coalescent simulation is assumed to be any program that prints biallelic marker output in the format of Dick Hudson's ms (http://home.uchicago.edu/~rhudson1)

This format looks like this:

//<br>
segsites: 2<br>
positions: 0.1 0.9<br>
01<br>
10<br>

The above is for 2 "SNPs" in a sample of size n = 2

#Installation

./configure

make

sudo make install

##If dependent libraries are in non-standard locations.  For example, "/opt":

CXXFLAGS=-I/opt/include LDFLAGS="$LDFLAGS -L/opt/lib" ./configure

make 

sudo make install

##Installing somewhere other than /usr/local/bin

./configure --prefix=/path/to/where/you/want/it

For example,

./configure --prefix=$HOME

make 

make install

will result in msstats being in ~/bin.

##Notes on calculations for metapopulations.

If the input data contain ordered samples from multiple demes, you may pass that info to _msstats_ as follows:

-I D n0 n1 ... nD,

where D is the number of demes and n0 is the sample size in the first deme, etc.

When -I is used, summary statistics will be calculated for each deme within each replicate.  The columns "rep" and "pop" will tell you which deme in which replicate each statistic corresponds to.

The program will exit with an error if the sum of deme sizes does not equal the input sample size read from STDIN.

When the -I option is used, _msstats_ will report Hudson, Slatkin, and Maddison's Fst statistic.  The reference for this statistic is:

Hudson, Slatkin and Maddison (1992) Estimation of levels of gene flow from population data. Genetics 132:583-589

For all pairwise comparisions amongst the D demes, you will see output columns with headers hsmij, which is the Fst statistic calculated between demes _i_ and _j_.  These numbers will be repeated for each deme within each replicate.  Thus, to get the actual distribution of Fst for a simulation, you should condition on the line of results for just one population.  For example:

ms 30 100 -t 10 -I 3 10 10 10 1 | ./src/msstats -I 3 10 10 10 > output<br>

Then, using __R__,

> x=read.table("output",header=T)<br>
> mean(x$hsm01[x$pop==0])<br>
[1] 0.3511204<br>

#The output

##For a single population

There are 21 columns in the output:

1. S = the number of "segregating sites", aka mutations. ([Watterson, 1975](http://www.ncbi.nlm.nih.gov/pubmed/1145509))
2. n1 = the number of singletons in the data.  This is the number of mutations at frequency 1 and n-1 in the sample.
3. next = the number of "external mutations" (sensu [Fu](http://www.ncbi.nlm.nih.gov/pubmed/7482370)) = the number of derived singletons.
4. theta = [Watterson's](http://www.ncbi.nlm.nih.gov/pubmed/1145509) estimate of theta
5. pi = "sum of site heterogzygosity" = sum of 2pq over the S sites. (Nei, [Tajima](http://www.genetics.org/content/105/2/437.abstract), others).
6. thetaH = [Fay and Wu's](http://www.genetics.org/content/155/3/1405.abstract) estimator of theta.  Their H statistic is pi - thetaH, hence no column for it.
7. Hprime = [Zeng et al.'s](http://www.genetics.org/content/174/3/1431.abstract) normalized Fay and Wu's H.  
8. tajd = [Tajima's D](http://www.genetics.org/content/123/3/585.abstract)
9. fulif = [Fu and Li's](http://www.genetics.org/content/133/3/693.abstract) F
10. fulid = [Fu and Li's](http://www.genetics.org/content/133/3/693.abstract) D
11. fulifs = [Fu and Li's](http://www.genetics.org/content/133/3/693.abstract) F-star
12. fulifds = [Fu and Li's](http://www.genetics.org/content/133/3/693.abstract) D-star
13. rm = [Hudson and Kaplan's](http://www.genetics.org/content/111/1/147.abstract) Rmin 
14. rmmg = [Myers and Griffiths](http://www.ncbi.nlm.nih.gov/pubmed/12586723) simplest lowest bound on Rmin
15. nhap = Number of distinct haplotypes. [Depaulis and Veuille](http://mbe.oxfordjournals.org/content/15/12/1788)
16. hdiv = haplotype diversity. [Depaulis and Veuille](http://mbe.oxfordjournals.org/content/15/12/1788)
17. wallB = [Jeff Wall's](https://www.cambridge.org/core/journals/genetics-research/article/div-classtitlerecombination-and-the-power-of-statistical-tests-of-neutralitydiv/F66762884D13A7215EF59D45189D9E4B) B
18. wallQ = [Jeff Wall's](https://www.cambridge.org/core/journals/genetics-research/article/div-classtitlerecombination-and-the-power-of-statistical-tests-of-neutralitydiv/F66762884D13A7215EF59D45189D9E4B) Q
19. rosarf = [Ramos-Onsins and Rozas](http://mbe.oxfordjournals.org/content/19/12/2092) Rf statistic
20. rosarf = [Ramos-Onsins and Rozas](http://mbe.oxfordjournals.org/content/19/12/2092) Ru statistic
21. zns = [Kelly's](http://www.genetics.org/content/146/3/1197) Zns = average pairwise r-squared.
