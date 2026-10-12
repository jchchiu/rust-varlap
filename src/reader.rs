use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context, Result};

use flate2::read::MultiGzDecoder;

pub enum MaybeGzipReader {
    Plain(File),
    Gzip(MultiGzDecoder<File>),
}

impl Read for MaybeGzipReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            MaybeGzipReader::Plain(file) => file.read(buf),
            MaybeGzipReader::Gzip(decoder) => decoder.read(buf),
        }
    }
}

pub fn open_variant_input(path: &Path) -> Result<MaybeGzipReader> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open variants file '{}'", path.display()))?;

    if is_gzip(&mut file)? {
        Ok(MaybeGzipReader::Gzip(MultiGzDecoder::new(file)))
    } else {
        Ok(MaybeGzipReader::Plain(file))
    }
}

// USING GZIP HANDLER AS PER RUSTQC
// https://github.com/seqeralabs/RustQC/blob/main/src/io.rs
// REWRITE MYSELF LATER

// https://cseweb.ucsd.edu/classes/sp22/cse223B-a/tribbler/flate2/read/struct.MultiGzDecoder.html
// Need to use multigzdecoder

/// Gzip magic bytes: the first two bytes of any gzip-compressed file.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

fn is_gzip(file: &mut File) -> Result<bool> {
    let mut magic = [0; 2];
    let bytes = file
        .read(&mut magic)
        .with_context(|| format!("Failed to read first two bytes from '{:?}'", file))?;

    // Go back to the start of the file
    file.seek(SeekFrom::Start(0))
        .with_context(|| format!("Failed to reseek start of file '{:?}'", file))?;

    Ok(bytes >= 2 && magic == GZIP_MAGIC)
}