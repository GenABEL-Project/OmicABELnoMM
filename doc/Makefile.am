# NOTE: This is not a recursive makefile. This file is to be included
# in the top Makefile.am

## Stuff needed for documentation in the doc/ directory
dist_doc_DATA = doc/howtocompile.txt doc/UserGuide.tex	\
 COPYING INSTALL ChangeLog AUTHORS

## Build the PDF documentation if building of the LaTeX docs is
## enabled via ./configure.
if BUILD_latexdoc
if HAVE_PDFLATEX
DOCDIR = doc
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

SUFFIXES = .aux .tex .pdf
# Several make rules to generate the PDF from the LaTeX source
.aux.pdf:
	@echo === Making PDF: $@ from $^ ===
	$(MKDIR_P) $(DOCDIR)
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
	$(MKDIR_P) $(DOCDIR)
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
