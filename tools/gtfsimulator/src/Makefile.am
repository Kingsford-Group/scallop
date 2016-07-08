bin_PROGRAMS=gtfsimulator

gtfsimulator_CPPFLAGS = -I../../../lib/gtf/src -I../../../lib/util/src
gtfsimulator_LDFLAGS = -L../../../lib/gtf/src -L../../../lib/util/src
gtfsimulator_LDADD = -lgtf -lutil

gtfsimulator_SOURCES = main.cc \
					   simulator.h simulator.cc 
