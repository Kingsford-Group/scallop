bin_PROGRAMS=scallop

# add htslib library
samtools = /home/mingfus/Data/Tools/samtools/samtools-1.2
hts = $(samtools)/htslib-1.2.1
hts_headers = $(hts)/htslib

scallop_CPPFLAGS = -I$(hts) -I$(hts_headers) 
scallop_LDFLAGS = -L$(samtools) -L$(hts)
scallop_LDADD = -lbam -lhts -lz

scallop_SOURCES = main.cc
