#!/bin/bash -eu
# OSS-Fuzz build script for pfloat.
#
# Compiles each `fuzz/fuzz_targets/*.rs` binary with the
# libFuzzer sanitizer and copies the resulting executables to
# $OUT/. ADR-0013 records the layout and the
# no-checked-in-corpus policy.

cd $SRC/pfloat/fuzz

cargo +nightly fuzz build -O --debug-assertions

FUZZ_TARGETS=(parse arith exp_log_family trig hyperbolic specials fmt)

for target in "${FUZZ_TARGETS[@]}"; do
  cp ../target/x86_64-unknown-linux-gnu/release/"$target" "$OUT"/
done
