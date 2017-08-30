bin_PROGRAMS = scallop

GTFDIR = ../lib/gtf
UTILDIR = ../lib/util
GRAPHDIR = ../lib/graph

scallop_CPPFLAGS = -std=c++11 -I$(GTFDIR) -I$(GRAPHDIR) -I$(UTILDIR)
scallop_LDFLAGS = -pthread -L$(GTFDIR) -L$(GRAPHDIR) -L$(UTILDIR)
scallop_LDADD = -lgtf -lgraph -lutil

scallop_SOURCES = splice_graph.h splice_graph.cc \
				  super_graph.h super_graph.cc \
				  sgraph_compare.h sgraph_compare.cc \
				  vertex_info.h vertex_info.cc \
				  edge_info.h edge_info.cc \
				  interval_map.h interval_map.cc \
				  config.h config.cc \
				  hit.h hit.cc \
				  partial_exon.h partial_exon.cc \
				  hyper_set.h hyper_set.cc \
				  subsetsum.h subsetsum.cc \
				  router.h router.cc \
				  region.h region.cc \
				  junction.h junction.cc \
				  bundle_base.h bundle_base.cc \
				  bundle.h bundle.cc \
				  path.h path.cc \
				  equation.h equation.cc \
				  gtf.h gtf.cc \
				  scallop.h scallop.cc \
				  previewer.h previewer.cc \
				  assembler.h assembler.cc \
				  filter.h filter.cc \
				  main.cc
