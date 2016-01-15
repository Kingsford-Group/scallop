bin_PROGRAMS=scallop

# add htslib library
samtools = $(SAMTOOLS)
hts = $(HTSLIB)
hts_headers = $(hts)/htslib

scallop_CPPFLAGS = -I$(hts) -I$(hts_headers)
scallop_LDFLAGS = -L$(samtools) -L$(hts) -pthread
scallop_LDADD = -lbam -lhts -lz

scallop_SOURCES = main.cc \
				  scallop.h scallop.cc \
				  bundle.h bundle.cc \
				  config.h config.cc \
				  hit.h hit.cc \
				  boundary.h boundary.cc \
				  region.h region.cc \
				  bridge.h bridge.cc \
				  sgraph.h sgraph.cc \
				  bbase.h bbase.cc \
				  dgraph.h dgraph.cc \
				  imap.h imap.cc \
				  binomial.h binomial.cc \
				  draw.h draw.cc
