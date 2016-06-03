# NOTE: This is not a recursive makefile. This file is to be included
# in the top Makefile.am

exmpldir = examples

examplescripts = $(exmpldir)/CreateData.R	\
 $(exmpldir)/create_Y.R				\
 $(exmpldir)/create_XL.R				\
 $(exmpldir)/create_XR.R				\
 $(exmpldir)/settings.R


examplefiles = $(exmpldir)/exclude_individuals.txt	\
 $(exmpldir)/dosages_1.txt				\
 $(exmpldir)/dosages_2.txt

## Specify where to store the example data
exampledir = $(pkgdatadir)/examples
dist_example_DATA = $(examplefiles) $(examplescripts)

cleanfiles_examples=

## These files are large (well, XR.fvd is), so we don't distribute
## them separately. These are used by the tests (see Makefile.am in
## the tests directory)
largeexamplefiles =	\
 $(exmpldir)/XL.fvd	\
 $(exmpldir)/XL.fvi	\
 $(exmpldir)/XR.fvd	\
 $(exmpldir)/XR.fvi	\
 $(exmpldir)/Y.fvd	\
 $(exmpldir)/Y.fvi


## Define how to create the example filevector files
## Note that $(<F) is a GNU make feature. It contains the directory
## part and the file-within-directory part of the first prerequisite.
$(exmpldir)/%.fvi $(exmpldir)/%.fvd: $(exmpldir)/create_%.R $(exmpldir)/settings.R
	cd $(exmpldir) && R --vanilla --no-save < $(<F)

## An explicit rule for XR because it depends on XL and Y
$(exmpldir)/XR.fvi $(exmpldir)/XR.fvi: $(exmpldir)/create_XR.R $(exmpldir)/settings.R $(exmpldir)/XL.fvd $(exmpldir)/XL.fvi $(exmpldir)/Y.fvd $(exmpldir)/Y.fvi


## These are small data files for interaction examples. These are
## distributed.
intexmpldir = $(exmpldir)/interactions
intexamplefiles =		\
 $(intexmpldir)/INT011.fvd	\
 $(intexmpldir)/INT011.fvi	\
 $(intexmpldir)/INT01R.fvd	\
 $(intexmpldir)/INT01R.fvi	\
 $(intexmpldir)/INT11.fvd	\
 $(intexmpldir)/INT11.fvi	\
 $(intexmpldir)/INTR.fvd	\
 $(intexmpldir)/INTR.fvi	\
 $(intexmpldir)/INT0.fvd	\
 $(intexmpldir)/INT0.fvi	\
 $(intexmpldir)/INTR0.fvd	\
 $(intexmpldir)/INTR0.fvi	\
 $(intexmpldir)/INT01.fvd	\
 $(intexmpldir)/INT01.fvi	\
 $(intexmpldir)/INT1.fvd	\
 $(intexmpldir)/INT1.fvi	\
 $(intexmpldir)/INTRR.fvd	\
 $(intexmpldir)/INTRR.fvi	\
 $(intexmpldir)/INT1RR.fvd	\
 $(intexmpldir)/INT1RR.fvi	\
 $(intexmpldir)/INT10.fvd	\
 $(intexmpldir)/INT10.fvi	\
 $(intexmpldir)/INTR1.fvd	\
 $(intexmpldir)/INTR1.fvi	\
 $(intexmpldir)/XL.fvd		\
 $(intexmpldir)/XL.fvi		\
 $(intexmpldir)/XR.fvd		\
 $(intexmpldir)/XR.fvi		\
 $(intexmpldir)/Y.fvd		\
 $(intexmpldir)/Y.fvi

intexampledir = $(pkgdatadir)/examples/interactions
dist_intexample_DATA = $(intexamplefiles)

CLEANFILES += $(cleanfiles_examples)
