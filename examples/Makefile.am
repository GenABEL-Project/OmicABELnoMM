# NOTE: This is not a recursive makefile. This file is to be included
# in the top Makefile.am

exmpldir = examples

examplescripts = CreateData.R

examplefiles = $(exmpldir)/exclude_individuals.txt	\
 $(exmpldir)/dosages_1.txt				\
 $(exmpldir)/dosages_2.txt

largeexamplefiles = $(exmpldir)/XL.fvd		\
 $(exmpldir)/XL.fvi				\
 $(exmpldir)/XR.fvd				\
 $(exmpldir)/XR.fvi				\
 $(exmpldir)/Y.fvd				\
 $(exmpldir)/Y.fvi

intexmpldir = examples/interactions
intexamplefiles = $(intexmpldir)/INT011.fvd	\
 $(intexmpldir)/INT01R.fvi			\
 $(intexmpldir)/INT11.fvd			\
 $(intexmpldir)/INT1RR.fvi			\
 $(intexmpldir)/INTR.fvd			\
 $(intexmpldir)/XL.fvi				\
 $(intexmpldir)/INT011.fvi			\
 $(intexmpldir)/INT0.fvd			\
 $(intexmpldir)/INT11.fvi			\
 $(intexmpldir)/INTR0.fvd			\
 $(intexmpldir)/INTR.fvi			\
 $(intexmpldir)/XR.fvd				\
 $(intexmpldir)/INT01.fvd			\
 $(intexmpldir)/INT0.fvi			\
 $(intexmpldir)/INT1.fvd			\
 $(intexmpldir)/INTR0.fvi			\
 $(intexmpldir)/INTRR.fvd			\
 $(intexmpldir)/XR.fvi				\
 $(intexmpldir)/INT01.fvi			\
 $(intexmpldir)/INT10.fvd			\
 $(intexmpldir)/INT1.fvi			\
 $(intexmpldir)/INTR1.fvd			\
 $(intexmpldir)/INTRR.fvi			\
 $(intexmpldir)/Y.fvd				\
 $(intexmpldir)/INT01R.fvd			\
 $(intexmpldir)/INT10.fvi			\
 $(intexmpldir)/INT1RR.fvd			\
 $(intexmpldir)/INTR1.fvi			\
 $(intexmpldir)/XL.fvd				\
 $(intexmpldir)/Y.fvi


## Specify where to store the example data
exampledir = $(pkgdatadir)/examples
dist_example_DATA = $(examplefiles)

if BUILD_examples
largeexamplefiles: $(exmpldir)/$(examplescripts) $(examplefiles)
	cd $(exmpldir) && R --vanilla --no-save < $(examplescripts)

dist_example_DATA += $(largeexamplefiles)
endif

intexampledir = $(pkgdatadir)/examples/interactions
dist_intexample_DATA = $(intexamplefiles)

cleanfiles_examples = $(largeexamplefiles)

CLEANFILES += $(cleanfiles_examples)
