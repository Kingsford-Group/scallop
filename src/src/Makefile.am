bin_PROGRAMS=scallop

GTFDIR = ../../lib/gtf/src
UTILDIR = ../../lib/util/src
GRAPHDIR = ../../lib/graph/src

scallop_CPPFLAGS = -std=c++11 -I$(HTSLIB) -I$(GUROBI_HOME)/include -I$(BOOST_HOME) -I$(GTFDIR) -I$(GRAPHDIR) -I$(UTILDIR)
scallop_LDFLAGS = -L$(HTSLIB) -L$(GUROBI_HOME)/lib -L$(GTFDIR) -L$(GRAPHDIR) -L$(UTILDIR) -pthread -Wl,-rpath,$(GUROBI_HOME)/lib
scallop_LDADD = -lhts -lz -lgurobi_c++ -lgurobi65 -lgtf -lgraph -lutil

scallop_SOURCES = splice_graph.h splice_graph.cc \
				  super_graph.h super_graph.cc \
				  sgraph_compare.h sgraph_compare.cc \
				  vertex_info.h vertex_info.cc \
				  edge_info.h edge_info.cc \
				  interval_map.h interval_map.cc \
				  subsetsum.h subsetsum.cc \
				  config.h config.cc \
				  hit.h hit.cc \
				  partial_exon.h partial_exon.cc \
				  hyper_set.h hyper_set.cc \
				  router.h router.cc \
				  region.h region.cc \
				  junction.h junction.cc \
				  bundle_base.h bundle_base.cc \
				  bundle.h bundle.cc \
				  path.h path.cc \
				  smoother.h smoother.cc \
				  equation.h equation.cc \
				  gtf.h gtf.cc \
				  scallop.h scallop.cc \
				  assembler.h assembler.cc \
				  estimator.h estimator.cc \
				  filter.h filter.cc \
				  main.cc
