bin_PROGRAMS=scallop

scallop_CPPFLAGS = -std=c++11 -I$(HTSLIB) -I$(GUROBI_HOME)/include -I$(BOOST_HOME)
scallop_LDFLAGS = -L$(HTSLIB) -L$(GUROBI_HOME)/lib -pthread
scallop_LDADD = -lhts -lz -lgurobi_c++ -lgurobi

scallop_SOURCES = edge_base.h edge_base.cc \
				  vertex_base.h vertex_base.cc \
				  graph_base.h graph_base.cc \
				  directed_graph.h directed_graph.cc \
				  undirected_graph.h undirected_graph.cc \
				  splice_graph.h splice_graph.cc \
				  nested_graph.h nested_graph.cc \
				  interval_map.h interval_map.cc \
				  binomial.h binomial.cc \
				  fibonacci_heap.h fibonacci_heap.cc \
				  subsetsum.h subsetsum.cc \
				  algebra.h algebra.cc \
				  disjoint_sets.h disjoint_sets.cc \
				  config.h config.cc \
				  hit.h hit.cc \
				  boundary.h boundary.cc \
				  region.h region.cc \
				  junction.h junction.cc \
				  bundle_base.h bundle_base.cc \
				  bundle.h bundle.cc \
				  draw.h draw.cc \
				  path.h path.cc \
				  smoother.h smoother.cc \
				  assembler.h assembler.cc \
				  stringtie.h stringtie.cc \
				  scallop.h scallop.cc \
				  gtf_exon.h gtf_exon.cc \
				  gtf_gene.h gtf_gene.cc \
				  path_space.h path_space.cc \
				  manager.h manager.cc \
				  main.cc
