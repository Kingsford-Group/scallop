bin_PROGRAMS=hayer2gtf

GTFDIR = ../../../lib/gtf/src
UTILDIR = ../../../lib/util/src

hayer2gtf_CPPFLAGS = -I$(GTFDIR) -I$(UTILDIR)
hayer2gtf_LDFLAGS = -L$(GTFDIR) -L$(UTILDIR)
hayer2gtf_LDADD = -lgtf -lutil

hayer2gtf_SOURCES = hayer.h hayer.cc \
					main.cc
