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
cleanfiles_tests = examples/results/dosages__dis.txt	\
 examples/results/dosages_add__dis.txt			\
 examples/results/dosages_dom__dis.txt			\
 examples/results/dosages_excl__dis.txt			\
 examples/results/dosages_excl_add__dis.txt		\
 examples/results/dosages_excl_dom__dis.txt		\
 examples/results/dosages_excl_myadd__dis.txt		\
 examples/results/dosages_excl_mylin__dis.txt		\
 examples/results/dosages_excl_res__dis.txt		\
 examples/results/dosages_myadd__dis.txt		\
 examples/results/dosages_mylin__dis.txt		\
 examples/results/dosages_res__dis.txt			\
 examples/results/multi_inter_011_dis.txt		\
 examples/results/multi_inter_0_dis.txt			\
 examples/results/multi_inter_10_dis.txt		\
 examples/results/multi_inter_11_dis.txt		\
 examples/results/multi_inter_1RR_dis.txt		\
 examples/results/multi_inter_1_dis.txt			\
 examples/results/multi_inter_R0_dis.txt		\
 examples/results/multi_inter_R1_dis.txt		\
 examples/results/multi_inter_R_dis.txt			\
 examples/results/multi_inter_add_011_dis.txt		\
 examples/results/multi_inter_add_0_dis.txt		\
 examples/results/multi_inter_add_10_dis.txt		\
 examples/results/multi_inter_add_11_dis.txt		\
 examples/results/multi_inter_add_1RR_dis.txt		\
 examples/results/multi_inter_add_1_dis.txt		\
 examples/results/multi_inter_add_R0_dis.txt		\
 examples/results/multi_inter_add_R1_dis.txt		\
 examples/results/multi_inter_add_R_dis.txt		\
 examples/results/multi_inter_dom_011_dis.txt		\
 examples/results/multi_inter_dom_0_dis.txt		\
 examples/results/multi_inter_dom_10_dis.txt		\
 examples/results/multi_inter_dom_11_dis.txt		\
 examples/results/multi_inter_dom_1RR_dis.txt		\
 examples/results/multi_inter_dom_1_dis.txt		\
 examples/results/multi_inter_dom_R0_dis.txt		\
 examples/results/multi_inter_dom_R1_dis.txt		\
 examples/results/multi_inter_dom_R_dis.txt		\
 examples/results/multi_inter_keep_R0_dis.txt		\
 examples/results/multi_inter_keep_RR_dis.txt		\
 examples/results/multi_inter_keep_R_dis.txt		\
 examples/results/multi_inter_keep_add_R0_dis.txt	\
 examples/results/multi_inter_keep_add_RR_dis.txt	\
 examples/results/multi_inter_keep_add_R_dis.txt	\
 examples/results/multi_inter_keep_dom_R0_dis.txt	\
 examples/results/multi_inter_keep_dom_RR_dis.txt	\
 examples/results/multi_inter_keep_dom_R_dis.txt	\
 examples/results/multi_inter_keep_myadd_R0_dis.txt	\
 examples/results/multi_inter_keep_myadd_RR_dis.txt	\
 examples/results/multi_inter_keep_myadd_R_dis.txt	\
 examples/results/multi_inter_keep_mylin_R0_dis.txt	\
 examples/results/multi_inter_keep_mylin_RR_dis.txt	\
 examples/results/multi_inter_keep_mylin_R_dis.txt	\
 examples/results/multi_inter_keep_res_R0_dis.txt	\
 examples/results/multi_inter_keep_res_RR_dis.txt	\
 examples/results/multi_inter_keep_res_R_dis.txt	\
 examples/results/multi_inter_myadd_011_dis.txt		\
 examples/results/multi_inter_myadd_0_dis.txt		\
 examples/results/multi_inter_myadd_10_dis.txt		\
 examples/results/multi_inter_myadd_11_dis.txt		\
 examples/results/multi_inter_myadd_1RR_dis.txt		\
 examples/results/multi_inter_myadd_1_dis.txt		\
 examples/results/multi_inter_myadd_R0_dis.txt		\
 examples/results/multi_inter_myadd_R1_dis.txt		\
 examples/results/multi_inter_myadd_R_dis.txt		\
 examples/results/multi_inter_mylin_011_dis.txt		\
 examples/results/multi_inter_mylin_0_dis.txt		\
 examples/results/multi_inter_mylin_10_dis.txt		\
 examples/results/multi_inter_mylin_11_dis.txt		\
 examples/results/multi_inter_mylin_1RR_dis.txt		\
 examples/results/multi_inter_mylin_1_dis.txt		\
 examples/results/multi_inter_mylin_R0_dis.txt		\
 examples/results/multi_inter_mylin_R1_dis.txt		\
 examples/results/multi_inter_mylin_R_dis.txt		\
 examples/results/multi_inter_res_011_dis.txt		\
 examples/results/multi_inter_res_0_dis.txt		\
 examples/results/multi_inter_res_10_dis.txt		\
 examples/results/multi_inter_res_11_dis.txt		\
 examples/results/multi_inter_res_1RR_dis.txt		\
 examples/results/multi_inter_res_1_dis.txt		\
 examples/results/multi_inter_res_R0_dis.txt		\
 examples/results/multi_inter_res_R1_dis.txt		\
 examples/results/multi_inter_res_R_dis.txt		\
 examples/results/normal_dis.txt			\
 examples/results/single_inter_011_dis.txt		\
 examples/results/single_inter_0_dis.txt		\
 examples/results/single_inter_10_dis.txt		\
 examples/results/single_inter_11_dis.txt		\
 examples/results/single_inter_1RR_dis.txt		\
 examples/results/single_inter_1_dis.txt		\
 examples/results/single_inter_R0_dis.txt		\
 examples/results/single_inter_R1_dis.txt		\
 examples/results/single_inter_R_dis.txt		\
 examples/results/single_inter_add_011_dis.txt		\
 examples/results/single_inter_add_0_dis.txt		\
 examples/results/single_inter_add_10_dis.txt		\
 examples/results/single_inter_add_11_dis.txt		\
 examples/results/single_inter_add_1RR_dis.txt		\
 examples/results/single_inter_add_1_dis.txt		\
 examples/results/single_inter_add_R0_dis.txt		\
 examples/results/single_inter_add_R1_dis.txt		\
 examples/results/single_inter_add_R_dis.txt		\
 examples/results/single_inter_dom_011_dis.txt		\
 examples/results/single_inter_dom_0_dis.txt		\
 examples/results/single_inter_dom_10_dis.txt		\
 examples/results/single_inter_dom_11_dis.txt		\
 examples/results/single_inter_dom_1RR_dis.txt		\
 examples/results/single_inter_dom_1_dis.txt		\
 examples/results/single_inter_dom_R0_dis.txt		\
 examples/results/single_inter_dom_R1_dis.txt		\
 examples/results/single_inter_dom_R_dis.txt		\
 examples/results/single_inter_keep_R0_dis.txt		\
 examples/results/single_inter_keep_RR_dis.txt		\
 examples/results/single_inter_keep_R_dis.txt		\
 examples/results/single_inter_keep_add_R0_dis.txt	\
 examples/results/single_inter_keep_add_RR_dis.txt	\
 examples/results/single_inter_keep_add_R_dis.txt	\
 examples/results/single_inter_keep_dom_R0_dis.txt	\
 examples/results/single_inter_keep_dom_RR_dis.txt	\
 examples/results/single_inter_keep_dom_R_dis.txt	\
 examples/results/single_inter_keep_myadd_R0_dis.txt	\
 examples/results/single_inter_keep_myadd_RR_dis.txt	\
 examples/results/single_inter_keep_myadd_R_dis.txt	\
 examples/results/single_inter_keep_mylin_R0_dis.txt	\
 examples/results/single_inter_keep_mylin_RR_dis.txt	\
 examples/results/single_inter_keep_mylin_R_dis.txt	\
 examples/results/single_inter_keep_res_R0_dis.txt	\
 examples/results/single_inter_keep_res_RR_dis.txt	\
 examples/results/single_inter_keep_res_R_dis.txt	\
 examples/results/single_inter_myadd_011_dis.txt	\
 examples/results/single_inter_myadd_0_dis.txt		\
 examples/results/single_inter_myadd_10_dis.txt		\
 examples/results/single_inter_myadd_11_dis.txt		\
 examples/results/single_inter_myadd_1RR_dis.txt	\
 examples/results/single_inter_myadd_1_dis.txt		\
 examples/results/single_inter_myadd_R0_dis.txt		\
 examples/results/single_inter_myadd_R1_dis.txt		\
 examples/results/single_inter_myadd_R_dis.txt		\
 examples/results/single_inter_mylin_011_dis.txt	\
 examples/results/single_inter_mylin_0_dis.txt		\
 examples/results/single_inter_mylin_10_dis.txt		\
 examples/results/single_inter_mylin_11_dis.txt		\
 examples/results/single_inter_mylin_1RR_dis.txt	\
 examples/results/single_inter_mylin_1_dis.txt		\
 examples/results/single_inter_mylin_R0_dis.txt		\
 examples/results/single_inter_mylin_R1_dis.txt		\
 examples/results/single_inter_mylin_R_dis.txt		\
 examples/results/single_inter_res_011_dis.txt		\
 examples/results/single_inter_res_0_dis.txt		\
 examples/results/single_inter_res_10_dis.txt		\
 examples/results/single_inter_res_11_dis.txt		\
 examples/results/single_inter_res_1RR_dis.txt		\
 examples/results/single_inter_res_1_dis.txt		\
 examples/results/single_inter_res_R0_dis.txt		\
 examples/results/single_inter_res_R1_dis.txt		\
 examples/results/single_inter_res_R_dis.txt


CLEANFILES = tempAL.bin tempB.bin tempY.bin tempAR.bin	\
 $(cleanfiles_tests)


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
	@echo === Making $@ file from $^ ===
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
