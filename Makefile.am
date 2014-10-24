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


## Specify where to store the examples data
exampledir = $(pkgdatadir)/examples
dist_example_DATA = examples/CreateData.R examples/XL.fvd	\
	examples/XL.fvi  examples/XR.fvd  examples/XR.fvi	\
	examples/Y.fvd  examples/Y.fvi


## Stuff needed for documentation in the doc/ directory
dist_doc_DATA = doc/howtocompile.txt doc/UserGuide.tex	\
 COPYING INSTALL ChangeLog AUTHORS

## Build the PDF documentation if building of the LaTeX docs is
## enabled via ./configure.
if BUILD_latexdoc
if HAVE_PDFLATEX
DOCDIR = doc/
MANNAME = $(DOCDIR)/UserGuide
MANPDF = $(MANNAME).pdf
MANTEXSRC = $(MANNAME).tex
MANAUX = $(MANNAME).aux
doc_DATA = $(MANPDF)
PDFLATEXOPTS = --output-directory=$(DOCDIR)

CLEANFILES += $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
 $(MANNAME).toc $(MANNAME).idx $(MANNAME).ilg $(MANNAME).ind $(MANAUX) .btmp
endif
endif


# Several make rules to generate the PDF from the LaTeX source
.aux.pdf:
	@echo === Making PDF: $@ from $^ ===
	$(PDFLATEX) $(PDFLATEXOPTS) $(srcdir)/$(MANTEXSRC)
	@while ( grep "Rerun to "			\
		$(MANNAME).log  ); do			\
		echo '** Re-running LaTeX **';		\
		$(PDFLATEX) $(PDFLATEXOPTS) $(srcdir)/$(MANTEXSRC);	\
	done
	if [ -f $(MANNAME).idx ]; then \
		echo === Making index ===; \
		makeindex $(MANNAME); \
	fi
	@echo === Making final PDF ===
	$(PDFLATEX) $(PDFLATEXOPTS) $(srcdir)/$(MANTEXSRC)

.tex.aux:
	@echo === Making $@ file from $^, we are in $(shell pwd) ===
	$(PDFLATEX) $(PDFLATEXOPTS) $(srcdir)/$(MANTEXSRC)
## Look for citations. Make sure grep never returns an error code.
	@grep "^\\\\citation" $(MANAUX) > .btmp.new || true
## If the citations are not changed, don't do anything. Otherwise replace
## the .btmp file to make sure BibTeX will be run.
	@if ( diff .btmp.new .btmp  > /dev/null ); then \
		rm .btmp.new; \
	else \
		mv .btmp.new .btmp; \
	fi

# A target needed to keep track of the nr. of LaTeX runs
.btmp:
