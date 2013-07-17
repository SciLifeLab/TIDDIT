ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS=FindTranslocations FindTranslocations_debug


FindTranslocations_SOURCES=FindTranslocations.cpp \
				dataStructures/Translocations.cpp
	

	

FindTranslocations_debug_SOURCES=$(FindTranslocations_SOURCES)

SAMTOOLS=samtools/
INCLUDES = -I$(SAMTOOLS) 

FindTranslocations_CXXFLAGS= -O3 $(BOOST_CPPFLAGS)
FindTranslocations_debug_CXXFLAGS= -O0 -g3 $(BOOST_CPPFLAGS)
FindTranslocations_CFLAGS= -O3 
FindTranslocations_debug_CFLAGS= -O0 -g3 


FindTranslocations_LDFLAGS = -L$(SAMTOOLS) -lbam $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_THREAD_LIBS)  $(BOOST_IOSTREAMS_LIBS)  -lz -lbz2
FindTranslocations_LIBS = -L$(SAMTOOLS) -lbam $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_THREAD_LIBS) $(BOOST_IOSTREAMS_LIBS)  -lz -lbz2
FindTranslocations_debug_LDFLAGS = $(FindTranslocations_LDFLAGS)
FindTranslocations_debug_LIBS = $(FindTranslocations_LIBS)

