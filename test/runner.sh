#! /bin/sh
#
# Run tests directly from the test directory without perl-tap/prove

# add meson builddir
PATH=../builddir:$PATH

bash ../test/t/00_region_and_target_handling.t
echo $?
exit
bash ../test/t/01_call_variants.t
bash ../test/t/02_multi_bam.t
bash ../test/t/03_reference_bases.t
