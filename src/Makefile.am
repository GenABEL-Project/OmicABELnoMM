# NOTE: This is not a recursive makefile. This file is to be included
# in the top Makefile.am

bin_PROGRAMS = omicabelnomm omicout

oanomm_cpp = src/AIOwrapper.cpp src/Algorithm.cpp src/Utility.cpp
oanomm_headers = $(oanomm_cpp:.cpp=.h) src/Definitions.h

omicabelnomm_SOURCES = $(oanomm_headers) $(oanomm_cpp) src/main.cpp
omicabelnomm_CPPFLAGS = $(AM_CPPFLAGS)
if USE_MPI
omicabelnomm_CPPFLAGS += -DUSE_MPI
endif

omicout_SOURCES = OmicOUT/main.cpp
omicout_CPPFLAGS = $(AM_CPPFLAGS)