bin_PROGRAMS = omicabelnomm

oanomm_cpp = src/AIOwrapper.cpp src/Algorithm.cpp src/Utility.cpp
oanomm_headers = $(oanomm_cpp:.cpp=.h) src/Definitions.h

omicabelnomm_SOURCES = $(oanomm_headers) $(oanomm_cpp) src/main.cpp


## Stuff needed for the checks in the tests/ directory
check_PROGRAMS = tests/test_omicabelnomm
#bin_PROGRAMS += tests/test_omicabelnomm
tests_test_omicabelnomm_SOURCES =		\
 $(oanomm_headers)				\
 $(oanomm_cpp)					\
 tests/test.cpp

TESTS = $(check_PROGRAMS)
CLEANFILES = tempAL.bin tempB.bin tempY.bin tempAR.bin


## Stuff needed for documentation in the doc/ directory
dist_doc_DATA = doc/howtocompile.txt doc/UserGuide.tex \
COPYING INSTALL ChangeLog AUTHORS

## Specify where to store the examples data
exampledir = $(pkgdatadir)/examples
dist_example_DATA = examples/CreateData.R examples/XL.fvd	\
	examples/XL.fvi  examples/XR.fvd  examples/XR.fvi	\
	examples/Y.fvd  examples/Y.fvi
