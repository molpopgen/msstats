bin_PROGRAMS=msstats mssfs
msstats_SOURCES=msstats.cc otherstats.cc otherstats.hpp 
mssfs_SOURCES=mssfs.cc
AM_CXXFLAGS=-DNDEBUG -O3 -Wall -W
AM_CPPFLAGS=-I$(prefix)/include
#  uncomment the following if msstats requires the math library
#msstats_LDADD=-lm

#EXTRA_DIST=msstats.lsm.in #NAME.texinfo msstats.spec.in

#  if you write a self-test script named `chk', uncomment the
#  following and add `chk' to the EXTRA_DIST list
#TESTS=chk

#  build and install the .info pages
#info_TEXINFOS = msstats.texinfo
#msstats_TEXINFOS = gpl.texinfo

#  install the man pages
man_MANS=msstats.1
