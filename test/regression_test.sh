#!/usr/bin/env bash
#
# Simple regression test for freebayes
# Runs freebayes on test data and compares output to baseline
#
# Usage:
#   ./regression_test.sh               # Run regression test
#   ./regression_test.sh --update      # Update baseline (after verifying changes are correct)
#

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DATA_DIR="${SCRIPT_DIR}/data"
BASELINE_DIR="${BASELINE_DIR:-${SCRIPT_DIR}/regression_baseline}"
OUTPUT_DIR="${SCRIPT_DIR}/regression_output"

# Create directories
mkdir -p "${BASELINE_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Find freebayes executable (FREEBAYES env var overrides auto-detection)
if [ -z "${FREEBAYES}" ]; then
    if [ -f "${SCRIPT_DIR}/../build/freebayes" ]; then
        FREEBAYES="${SCRIPT_DIR}/../build/freebayes"
    elif command -v freebayes &> /dev/null; then
        FREEBAYES="freebayes"
    else
        echo -e "${RED}ERROR: freebayes not found${NC}"
        echo "Build it first: cd build && ninja"
        exit 1
    fi
fi

echo "Using freebayes: ${FREEBAYES}"

# Test cases
run_test() {
    local test_name="$1"
    local ref="$2"
    local bam="$3"
    local extra_args="$4"

    echo ""
    echo "Running test: ${test_name}"
    echo "  binary:    ${FREEBAYES}"
    echo "  reference: ${ref}"
    echo "  bam:       ${bam}"
    echo "  args:      ${extra_args:-<none>}"
    echo "  command:   ${FREEBAYES} -f ${ref} ${extra_args} ${bam}"

    local baseline="${BASELINE_DIR}/${test_name}.vcf"
    local output="${OUTPUT_DIR}/${test_name}.vcf"

    # Run freebayes
    ${FREEBAYES} -f "${ref}" ${extra_args} "${bam}" > "${output}" 2>/dev/null || {
        echo -e "${RED}FAILED: freebayes crashed${NC}"
        return 1
    }

    # Update mode: save as new baseline
    if [ "$UPDATE_BASELINE" = "true" ]; then
        cp "${output}" "${baseline}"
        echo -e "${YELLOW}Updated baseline${NC}"
        return 0
    fi

    # Check if baseline exists
    if [ ! -f "${baseline}" ]; then
        echo -e "${YELLOW}No baseline found, creating one${NC}"
        cp "${output}" "${baseline}"
        return 0
    fi

    # Print baseline file path for debugging
    echo "Baseline file: ${baseline}"

    # Compare output to baseline, ignoring volatile headers (date, path, commandline)
    filter_vcf() { grep -v '^##fileDate\|^##commandline\|^##reference'; }
    if diff -u <(filter_vcf < "${baseline}") <(filter_vcf < "${output}") > "${OUTPUT_DIR}/${test_name}.diff"; then
        echo -e "${GREEN}PASSED${NC}"
        rm "${OUTPUT_DIR}/${test_name}.diff"
        return 0
    else
        echo -e "${RED}FAILED: Output differs from baseline${NC}"
        echo "Diff saved to: ${OUTPUT_DIR}/${test_name}.diff"
        echo ""
        echo "First 20 lines of diff:"
        head -20 "${OUTPUT_DIR}/${test_name}.diff"
        return 1
    fi
}

# Parse command line
UPDATE_BASELINE=false
if [ "$1" = "--update" ]; then
    UPDATE_BASELINE=true
    echo -e "${YELLOW}UPDATE MODE: Will update baselines${NC}"
fi

# Run test suite
echo "========================================"
echo "FreeBayes Regression Test Suite"
echo "========================================"

FAILED=0
PASSED=0

# Test 1: Basic variant calling on tiny test
if run_test "basic" \
    "${TEST_DATA_DIR}/test.ref" \
    "${TEST_DATA_DIR}/test.bam" \
    ""; then
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi

# Test 2: With region specified
if run_test "region" \
    "${TEST_DATA_DIR}/test.ref" \
    "${TEST_DATA_DIR}/test.bam" \
    "-r ref:1-11"; then
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi

# Test 3: Different parameters
if run_test "min_alt_frac" \
    "${TEST_DATA_DIR}/test.ref" \
    "${TEST_DATA_DIR}/test.bam" \
    "-F 0.1 -C 1"; then
    PASSED=$((PASSED + 1))
else
    FAILED=$((FAILED + 1))
fi

# NA12878 chr22 tests (require tiny/ data)
TINY_REF="${SCRIPT_DIR}/tiny/q.fa"
TINY_BAM="${SCRIPT_DIR}/tiny/NA12878.chr22.tiny.bam"

run_tiny_test() {
    local name="$1"; shift
    if run_test "$name" "${TINY_REF}" "${TINY_BAM}" "$*"; then
        PASSED=$((PASSED + 1))
    else
        FAILED=$((FAILED + 1))
    fi
}

if [ -f "${TINY_REF}" ] && [ -f "${TINY_BAM}" ]; then

    # Regions
    run_tiny_test "tiny_chr22"              "-r q:1-1000"
    run_tiny_test "tiny_chr22_region2"      "-r q:1000-5000"
    run_tiny_test "tiny_chr22_region3"      "-r q:5000-10000"
    run_tiny_test "tiny_chr22_full"         ""

    # GVCF output (exercises non-variant block logic including cross-chromosome boundary fix)
    run_tiny_test "tiny_chr22_gvcf"         "-r q:1-1000 --gvcf"
    run_tiny_test "tiny_chr22_gvcf_nochunk" "-r q:1-1000 --gvcf --gvcf-dont-use-chunk true"

    # Quality filtering
    run_tiny_test "tiny_chr22_min_mapq"     "-r q:1-1000 --min-mapping-quality 20"
    run_tiny_test "tiny_chr22_min_baseq"    "-r q:1-1000 --min-base-quality 20"

    # Coverage limit
    run_tiny_test "tiny_chr22_limit_cov"    "-r q:1-1000 --limit-coverage 20"

    # Ploidy
    run_tiny_test "tiny_chr22_ploidy1"      "-r q:1-1000 --ploidy 1"

    # Haplotype length
    run_tiny_test "tiny_chr22_hap_len_0"    "-r q:1-1000 --haplotype-length 0"

fi

# Summary
echo ""
echo "========================================"
echo "Summary:"
echo "  Passed: ${PASSED}"
echo "  Failed: ${FAILED}"

if [ ${FAILED} -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed${NC}"
    exit 1
fi
