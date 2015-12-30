bin_PROGRAMS=scallop

# add htslib library
samtools = $(SAMTOOLS)
hts = $(HTSLIB)
hts_headers = $(hts)/htslib

scallop_CPPFLAGS = -I$(hts) -I$(hts_headers)
scallop_LDFLAGS = -L$(samtools) -L$(hts)
scallop_LDADD = -lbam -lhts -lz -lpthread

scallop_SOURCES = main.cc \
				  scallop.h scallop.cc \
				  bundle.h bundle.cc \
				  common.h common.cc \
				  config.h config.cc \
				  hit.h hit.cc \
				  boundary.h boundary.cc
