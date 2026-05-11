# OSS-Fuzz Dockerfile for pfloat.
#
# Copied to `projects/pfloat/Dockerfile` in the upstream
# google/oss-fuzz PR. ADR-0013 records the layout decision.

FROM gcr.io/oss-fuzz-base/base-builder-rust

RUN git clone --depth 1 https://github.com/ixmatus/pfloat.git pfloat

WORKDIR pfloat

COPY build.sh $SRC/
