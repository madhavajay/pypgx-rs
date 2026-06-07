# syntax=docker/dockerfile:1
#
# pypgx-rs — self-contained pharmacogenomics image: the Rust CLI + beagle-rs
# binary + the full pypgx-bundle reference data (1KGP phasing panels + CNV models
# converted to native Rust weights) all baked in. No compiling or data fetching
# needed at run time:
#
#   docker run --rm -v "$PWD":/data ghcr.io/madhavajay/pypgx-rs \
#       run-ngs-pipeline --vcf /data/sample.vcf.gz --assembly GRCh38 --output /data/out
#
# Multi-arch (linux/amd64 + linux/arm64): the Rust binaries build per target
# (under QEMU for the non-native arch in CI); the reference bundle is
# arch-independent and built once on the native build platform.

# ---- Stage 1: build the Rust binaries (pypgx CLI + beagle-rs) per target arch ----
FROM rust:1-bookworm AS rust-builder
WORKDIR /src
# Submodules (repos/beagle-rs, repos/samtools-rs) must be checked out in the
# build context before `docker build` (the CI does this).
COPY . .
RUN cargo build --release --features beagle --bin pypgx \
 && cargo build --release -p beagle-rs-cli --manifest-path repos/beagle-rs/Cargo.toml \
 && mkdir -p /out \
 && cp target/release/pypgx /out/pypgx \
 && cp repos/beagle-rs/target/release/beagle-rs /out/beagle-rs

# ---- Stage 2: assemble the reference bundle (arch-independent → build once) ----
# Pinned to $BUILDPLATFORM so the heavy fetch + sklearn-unpickle conversion runs
# natively once, not once per target arch under emulation.
FROM --platform=$BUILDPLATFORM python:3.10-slim-bookworm AS bundle-builder
ARG PYPGX_BUNDLE_VERSION=0.26.0
WORKDIR /build
RUN apt-get update \
 && apt-get install -y --no-install-recommends git bash ca-certificates \
 && rm -rf /var/lib/apt/lists/*
# Only scikit-learn + numpy are needed to unpickle the bundle's Model[CNV] (a
# pure sklearn object) and dump its weights — no pypgx/fuc. manylinux wheels, so
# no build toolchain required. Pinned to match the env that wrote the pickles.
RUN pip install --no-cache-dir numpy==2.2.6 scikit-learn==1.7.2
COPY tools/convert_cnv_models_all.py tools/fetch-bundle.sh /build/tools/
RUN bash /build/tools/fetch-bundle.sh /opt/pypgx-bundle

# ---- Stage 3: runtime ----
FROM debian:bookworm-slim AS runtime
LABEL org.opencontainers.image.source="https://github.com/madhavajay/pypgx-rs" \
      org.opencontainers.image.description="pypgx-rs CLI + beagle-rs + baked pypgx-bundle reference data"
# tabix: region-slices the input VCF in `run-ngs-pipeline`.
RUN apt-get update \
 && apt-get install -y --no-install-recommends tabix \
 && rm -rf /var/lib/apt/lists/*
COPY --from=rust-builder /out/pypgx     /usr/local/bin/pypgx
COPY --from=rust-builder /out/beagle-rs /usr/local/bin/beagle-rs
COPY --from=bundle-builder /opt/pypgx-bundle /opt/pypgx-bundle
ENV PYPGX_BUNDLE=/opt/pypgx-bundle \
    BEAGLE_RS_BIN=/usr/local/bin/beagle-rs
ENTRYPOINT ["pypgx"]
CMD ["--help"]
