# Runtime image for rust-varlap, built by .github/workflows/release.yml.
#
# The binary is compiled on native GitHub runners and copied in, so the image
# needs only the runtime libraries, not a Rust toolchain. The build context must
# contain <arch>/rust-varlap for each target platform (amd64, arm64).
FROM ubuntu:24.04

ARG TARGETARCH

# Runtime libraries that rust-htslib links against
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        libcurl4t64 \
        zlib1g \
        libbz2-1.0 \
        liblzma5 \
    && rm -rf /var/lib/apt/lists/*

COPY ${TARGETARCH}/rust-varlap /usr/local/bin/rust-varlap
RUN chmod +x /usr/local/bin/rust-varlap

ENV LC_ALL=C
WORKDIR /data
CMD ["rust-varlap", "--help"]
