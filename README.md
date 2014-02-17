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
