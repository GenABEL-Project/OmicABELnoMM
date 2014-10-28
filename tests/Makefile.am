# NOTE: This is not a recursive makefile. This file is to be included
# in the top Makefile.am



## Stuff needed for the checks in the tests/ directory
check_PROGRAMS = tests/test_omicabelnomm
if USE_MPI
tests_test_omicabelnomm_CXXFLAGS = -DUSE_MPI $(AM_CXXFLAGS)
endif

# List the sources for the test binary. Also include the empty
# directory for the results, so that it gets included in the tar-ball
# when running 'make dist'. Otherwise running 'make check' fails if
# run by the user (since the test program doesn't check for the
# directory's existence).
tests_test_omicabelnomm_SOURCES =		\
 $(oanomm_headers)				\
 $(oanomm_cpp)					\
 examples/results				\
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


CLEANFILES += $(cleanfiles_tests)
