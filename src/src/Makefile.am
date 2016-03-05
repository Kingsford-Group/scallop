bin_PROGRAMS=scallop

scallop_CPPFLAGS = -I$(HTSLIB) -I$(GUROBI_HOME)/include -I$(BOOST_HOME)
scallop_LDFLAGS = -L$(HTSLIB) -L$(GUROBI_HOME)/lib -pthread
scallop_LDADD = -lhts -lz -lgurobi_c++ -lgurobi

scallop_SOURCES = main.cc \
				  manager.h manager.cc \
				  bundle.h bundle.cc \
				  config.h config.cc \
				  hit.h hit.cc \
				  boundary.h boundary.cc \
				  region.h region.cc \
				  junction.h junction.cc \
				  splice_graph.h splice_graph.cc \
				  bundle_base.h bundle_base.cc \
				  interval_map.h interval_map.cc \
				  binomial.h binomial.cc \
				  draw.h draw.cc \
				  path.h path.cc \
				  fibonacci_heap.h fibonacci_heap.cc \
				  smoother.h smoother.cc \
				  assembler.h assembler.cc \
				  stringtie.h stringtie.cc \
				  scallop1.h scallop1.cc \
				  scallop2.h scallop2.cc \
				  scallop3.h scallop3.cc \
				  gtf_exon.h gtf_exon.cc \
				  gtf_gene.h gtf_gene.cc \
				  subsetsum.h subsetsum.cc \
				  decomposer.h decomposer.cc \
				  algebra.h algebra.cc 
