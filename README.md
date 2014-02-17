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