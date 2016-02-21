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
				  lpsolver.h lpsolver.cc \
				  path.h path.cc \
				  fibonacci_heap.h fibonacci_heap.cc \
				  assembler.h assembler.cc \
				  stringtie.h stringtie.cc \
				  scallop.h scallop.cc \
				  gtf_exon.h gtf_exon.cc \
				  gtf_gene.h gtf_gene.cc \
				  shannon.h shannon.cc
