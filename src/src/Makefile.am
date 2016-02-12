bin_PROGRAMS=scallop

scallop_CPPFLAGS = -I$(HTSLIB) -I$(GUROBI_HOME)/include -I$(BOOST_HOME)
scallop_LDFLAGS = -L$(HTSLIB) -L$(GUROBI_HOME)/lib -pthread
scallop_LDADD = -lhts -lz -lgurobi_c++ -lgurobi

scallop_SOURCES = main.cc \
				  dgraph.h dgraph.cc \
				  scallop.h scallop.cc \
				  bundle.h bundle.cc \
				  config.h config.cc \
				  hit.h hit.cc \
				  boundary.h boundary.cc \
				  region.h region.cc \
				  bridge.h bridge.cc \
				  sgraph.h sgraph.cc \
				  bbase.h bbase.cc \
				  imap.h imap.cc \
				  binomial.h binomial.cc \
				  draw.h draw.cc \
				  lpsolver.h lpsolver.cc \
				  path.h path.cc \
				  fheap.h fheap.cc
