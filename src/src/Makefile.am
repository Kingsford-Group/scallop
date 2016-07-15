bin_PROGRAMS=scallop

GTFDIR = ../../lib/gtf/src
UTILDIR = ../../lib/util/src
GRAPHDIR = ../../lib/graph/src

scallop_CPPFLAGS = -pg -std=c++11 -I$(HTSLIB) -I$(GUROBI_HOME)/include -I$(BOOST_HOME) -I$(GTFDIR) -I$(GRAPHDIR) -I$(UTILDIR)
scallop_LDFLAGS = -pg -L$(HTSLIB) -L$(GUROBI_HOME)/lib -L$(GTFDIR) -L$(GRAPHDIR) -L$(UTILDIR) -pthread
scallop_LDADD = -lhts -lz -lgurobi_c++ -lgurobi -lgtf -lgraph -lutil

scallop_SOURCES = splice_graph.h splice_graph.cc \
				  nested_graph.h nested_graph.cc \
				  interval_map.h interval_map.cc \
				  binomial.h binomial.cc \
				  fibonacci_heap.h fibonacci_heap.cc \
				  subsetsum.h subsetsum.cc \
				  subsetsum1.h subsetsum1.cc \
				  disjoint_sets.h disjoint_sets.cc \
				  config.h config.cc \
				  hit.h hit.cc \
				  boundary.h boundary.cc \
				  partial_exon.h partial_exon.cc \
				  region.h region.cc \
				  slope.h slope.cc \
				  junction.h junction.cc \
				  bundle_base.h bundle_base.cc \
				  bundle.h bundle.cc \
				  bundle_extend.h bundle_extend.cc \
				  path.h path.cc \
				  smoother.h smoother.cc \
				  equation.h equation.cc \
				  scallop1.h scallop1.cc \
				  scallop2.h scallop2.cc \
				  gtf.h gtf.cc \
				  assembler.h assembler.cc \
				  main.cc
