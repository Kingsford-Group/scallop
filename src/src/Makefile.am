bin_PROGRAMS=scallop

# add htslib library
htslib = $(HTSLIB)

scallop_CPPFLAGS = -I$(htslib) -I$(GUROBI_HOME)/include
scallop_LDFLAGS = -L$(htslib) -L$(GUROBI_HOME)/lib -pthread
scallop_LDADD = -lhts -lz -lgurobi_c++

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
				  draw.h draw.cc \
				  lpsolver.h lpsolver.cc
