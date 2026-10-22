# Installation

rust-varlap can be installed in several ways. Pick whichever suits your environment:

| Method | Best for |
|---|---|
| [Pre-built binary](#pre-built-binaries) | Quickest start on Linux or macOS |
| [Docker](#docker) | Containerised workflows, Nextflow / nf-core |
| [Apptainer / Singularity](#apptainer--singularity) | HPC clusters without Docker |
| [Bioconda](#bioconda) | Conda-based environments |
| [crates.io](#cratesio) | Rust users with a toolchain already installed |
| [Build from source](#build-from-source) | Development, or unsupported platforms |

The examples below use a `VERSION` variable. Set it once and reuse it:

```bash
VERSION="v0.1.0-alpha.5"
```

---

## Pre-built binaries

Binaries for Linux and macOS are available from the
[Releases page](https://github.com/jchchiu/rust-varlap/releases).

{{#tabs global="os" }}

{{#tab name="Linux x86_64" }}

```bash
wget https://github.com/jchchiu/rust-varlap/releases/download/${VERSION}/rust-varlap_${VERSION}_linux-x86_64
chmod +x rust-varlap_${VERSION}_linux-x86_64

# Optional: move to somewhere on your PATH
mv rust-varlap_${VERSION}_linux-x86_64 ~/.local/bin/rust-varlap

rust-varlap --version
```

{{#endtab }}

{{#tab name="Linux aarch64" }}

```bash
wget https://github.com/jchchiu/rust-varlap/releases/download/${VERSION}/rust-varlap_${VERSION}_linux-aarch64
chmod +x rust-varlap_${VERSION}_linux-aarch64

# Optional: move to somewhere on your PATH
mv rust-varlap_${VERSION}_linux-aarch64 ~/.local/bin/rust-varlap

rust-varlap --version
```

{{#endtab }}

{{#tab name="macOS Apple Silicon" }}

```bash
curl -fsSLO https://github.com/jchchiu/rust-varlap/releases/download/${VERSION}/rust-varlap_${VERSION}_macos-aarch64
chmod +x rust-varlap_${VERSION}_macos-aarch64

# Optional: move to somewhere on your PATH
mkdir -p ~/.local/bin
mv rust-varlap_${VERSION}_macos-aarch64 ~/.local/bin/rust-varlap

rust-varlap --version
```

{{#endtab }}

{{#tab name="macOS Intel" }}

```bash
curl -fsSLO https://github.com/jchchiu/rust-varlap/releases/download/${VERSION}/rust-varlap_${VERSION}_macos-x86_64
chmod +x rust-varlap_${VERSION}_macos-x86_64

# Optional: move to somewhere on your PATH
mkdir -p ~/.local/bin
mv rust-varlap_${VERSION}_macos-x86_64 ~/.local/bin/rust-varlap

rust-varlap --version
```

{{#endtab }}

{{#endtabs }}

> **Note:** While rust-varlap is in alpha, GitHub's `releases/latest` URL does not
> resolve to pre-releases, so always pin a specific `VERSION`.

---

## Docker

Run rust-varlap without any local installation:

```bash
docker run --rm -v $(pwd):/data ghcr.io/jchchiu/rust-varlap:${VERSION} \
  rust-varlap --version
```

Mount the directory containing your data at `/data` and refer to your input
files with `/data/...` paths.

---

## Apptainer / Singularity

Two options are available: pull the image from the GitHub Container Registry,
or download a pre-built `.sif` from the Releases page.

{{#tabs }}

{{#tab name="Pull from registry" }}

```bash
apptainer pull rust-varlap_${VERSION}.sif docker://ghcr.io/jchchiu/rust-varlap:${VERSION}
```

{{#endtab }}

{{#tab name="Download .sif" }}

Download the `.sif` image from the [Releases page](https://github.com/jchchiu/rust-varlap/releases).

```bash
wget https://github.com/jchchiu/rust-varlap/releases/download/${VERSION}/rust-varlap_${VERSION}.sif
```

{{#endtab }}

{{#endtabs }}

**Run:**

```bash
# Directly via the runscript
./rust-varlap_${VERSION}.sif --version

# Or explicitly with apptainer
apptainer run rust-varlap_${VERSION}.sif --version

# Bind directories when data lives outside your home folder
apptainer exec \
  -B /scratch/$USER:/scratch/$USER \
  rust-varlap_${VERSION}.sif \
  rust-varlap --version
```

> **Note:** An image pulled from the registry may not include a runscript. If
> `./rust-varlap_${VERSION}.sif` does not work, use `apptainer exec` as shown above.

---

## Bioconda

> **Coming soon.** A Bioconda recipe is planned. Once available, it will also
> provide a BioContainers image for use in Nextflow / nf-core pipelines.

```bash
conda install -c bioconda rust-varlap
```

---

## crates.io

If you already have a Rust toolchain, you can install from [crates.io](https://crates.io/crates/rust-varlap).
Pre-release versions must be requested explicitly:

```bash
cargo install rust-varlap --version 0.1.0-alpha.5 --locked
```

This places the binary in `~/.cargo/bin/rust-varlap`, which is on your `PATH` after a
standard rustup install. You will need the [system dependencies](#build-from-source)
listed below, as rust-varlap builds against [rust-htslib](https://github.com/rust-bio/rust-htslib).

---

## Build from source

Building from source requires Rust ≥ 1.70, a C toolchain, and a few system
libraries needed by rust-htslib for SAM/BAM/CRAM I/O.

### 1. Install system dependencies

{{#tabs global="os-build" }}

{{#tab name="Ubuntu / Debian" }}

```bash
sudo apt install cmake zlib1g-dev libbz2-dev liblzma-dev \
  libcurl4-openssl-dev libssl-dev clang
```

{{#endtab }}

{{#tab name="Fedora / RHEL" }}

```bash
sudo dnf install cmake zlib-devel bzip2-devel xz-devel \
  libcurl-devel openssl-devel clang
```

{{#endtab }}

{{#tab name="macOS" }}

```bash
brew install cmake zlib bzip2 xz curl openssl
```

{{#endtab }}

{{#endtabs }}

If Rust is not installed yet:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### 2. Build

{{#tabs }}

{{#tab name="cargo install (git tag)" }}

```bash
cargo install \
  --git https://github.com/jchchiu/rust-varlap \
  --tag ${VERSION} \
  --locked
```

The binary is placed in `~/.cargo/bin/rust-varlap`.

{{#endtab }}

{{#tab name="Git clone" }}

```bash
git clone https://github.com/jchchiu/rust-varlap.git
cd rust-varlap
git checkout ${VERSION}   # optional: build a specific release
cargo build --release

# Binary is at target/release/rust-varlap
./target/release/rust-varlap --version
```

Copy the binary to a directory on your `PATH` for convenient access.

{{#endtab }}

{{#endtabs }}

---

## Verify installation

```bash
rust-varlap --version
rust-varlap --help
```

If the binary is not on your `PATH`, use the full path instead, for example
`./target/release/rust-varlap --version`.