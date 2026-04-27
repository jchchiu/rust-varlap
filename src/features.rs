use std::rc::Rc;
use std::cmp;

use serde::Serialize;
use rust_htslib::bam::Record;
use rust_htslib::bam::record::{Aux, Cigar};

use crate::variant::VarType;

#[derive(Debug, Clone)]
pub enum LocusFeatures {
    Snv(LocusFeaturesSnv),
    Indel(LocusFeaturesIndel),
}

impl LocusFeatures {
    pub fn normalized_row(&self) -> NormalizedLocusFeaturesRow {
        match self {
            LocusFeatures::Snv(f) => f.common.normalized_row(),
            LocusFeatures::Indel(f) => f.common.normalized_row(),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct CommonLocusFeatures {
    pub ref_read_features: ReadFeatures,
    pub alt_read_features: ReadFeatures,
    pub all_read_features: ReadFeatures,
}

impl CommonLocusFeatures {
    pub fn normalized_row(&self) -> NormalizedLocusFeaturesRow {
        let r = self.ref_read_features.normalized();
        let a = self.alt_read_features.normalized();
        let all = self.all_read_features.normalized();

        NormalizedLocusFeaturesRow {
            ref_nm: r.nm,
            ref_base_qual: r.base_qual,
            ref_map_qual: r.map_qual,
            ref_align_len: r.align_len,
            ref_clipping: r.clipping,
            ref_indel: r.indel,
            ref_forward_strand: r.forward_strand,
            ref_reverse_strand: r.reverse_strand,
            ref_supplementary: r.supplementary,
            ref_normalised_read_position: r.normalised_read_position,

            alt_nm: a.nm,
            alt_base_qual: a.base_qual,
            alt_map_qual: a.map_qual,
            alt_align_len: a.align_len,
            alt_clipping: a.clipping,
            alt_indel: a.indel,
            alt_forward_strand: a.forward_strand,
            alt_reverse_strand: a.reverse_strand,
            alt_supplementary: a.supplementary,
            alt_normalised_read_position: a.normalised_read_position,

            all_nm: all.nm,
            all_base_qual: all.base_qual,
            all_map_qual: all.map_qual,
            all_align_len: all.align_len,
            all_clipping: all.clipping,
            all_indel: all.indel,
            all_forward_strand: all.forward_strand,
            all_reverse_strand: all.reverse_strand,
            all_supplementary: all.supplementary,
            all_normalised_read_position: all.normalised_read_position,
        }
    }    
}

#[derive(Debug, Clone, Default)]
pub struct LocusFeaturesSnv {
    pub base_counts: AlleleCountsSnv,
    pub common: CommonLocusFeatures,
}

impl LocusFeaturesSnv {
    pub fn count(
        &mut self,
        read: &Rc<Record>,
        base: Option<u8>,
        refr: char,
        alt: char,
        query_pos: Option<u32>,
    ) {
        if let Some(base_u8) = base {
            let base_char = base_u8 as char;

            self.base_counts.count(base_char);

            if base_char == refr {
                self.common.ref_read_features.count(read, query_pos);
            } else if base_char == alt {
                self.common.alt_read_features.count(read, query_pos);
            }
        }

        self.common.all_read_features.count(&read, query_pos);
    }
}

#[derive(Debug, Clone)]
pub struct IndelEvent {
    pub indel_type: VarType,
    pub start: u64,
    pub end: u64,
    pub bases: String,
}

#[derive(Debug, Clone, Default)]
pub struct LocusFeaturesIndel {
    pub overlapping_indels_count: u64,
    pub common: CommonLocusFeatures,
}

impl LocusFeaturesIndel {
    pub fn count(
        &mut self,
        read: &Rc<Record>,
        refr: &str,
        alt: &str,
        ref_pos: u64,
        query_pos: Option<u32>,
        indel_type: &VarType,
    ) {
        let start = self.get_indel_start_coord(ref_pos, refr, alt);
        let size = refr.len().abs_diff(alt.len()) as u64;
        let end = start + size - 1;

        let overlapping_indels = self.indels_overlapping_variant(&read, start, end);
        self.overlapping_indels_count += overlapping_indels.len() as u64;

        let mut read_supports_alt = false;

        let bases = match indel_type {
            VarType::Del => Some(String::new()),
            VarType::Ins => Some(alt[1..].to_string()),
            _ => None,
        };

        for event in overlapping_indels.iter() {
            if event.indel_type == *indel_type {
                if (event.start == start) && (event.end == end) {
                    if matches!(indel_type, VarType::Del) ||           // FIX UNWRAP HERE
                    (matches!(indel_type, VarType::Ins) && event.bases == bases.clone().unwrap()) {
                        read_supports_alt = true;
                        break;
                    }
                }
            }
        }

        let mut read_supports_ref = false;
        if overlapping_indels.is_empty() && let Some(qpos) = query_pos {
            let seq_bytes = read.seq().as_bytes();
            let read_bases = match indel_type {
                VarType::Ins => {
                    Some((seq_bytes[qpos as usize] as char)
                    .to_string()
                    .to_ascii_uppercase())
                },
                VarType::Del => {
                    Some(String::from_utf8(
                        seq_bytes[qpos as usize .. (qpos + (size as u32) + 1) as usize].to_vec()
                    )
                    .unwrap()
                    .to_ascii_uppercase())
                },
                _ => None,
            };
            
            // FIX UNWRAP HERE
            if refr == read_bases.unwrap() {
                read_supports_ref = true;
            }
        }

        if read_supports_ref {
            self.common.ref_read_features.count(read, query_pos);
        } else if read_supports_alt {
            self.common.alt_read_features.count(read, query_pos);
        }

        self.common.all_read_features.count(&read, query_pos);
    }

    // get the genome coordinates of where an INDEL variant will actually
    // start, as opposed to the location of where the variant is reported
    // the starting position must take into account the context bases that
    // are given when the variant is reported in the VCF file.
    // Note: if we assume normalized biallelic variants, could just do (pos + 1)
    pub fn get_indel_start_coord(&self, pos: u64, refr: &str, alt: &str) -> u64 {
        let shortest_len = cmp::min(refr.len(), alt.len());
        pos + shortest_len as u64
    }

    // True if the intervals of two indels overlap
    pub fn interval_overlaps(&self, start1: u64, end1: u64, start2: u64, end2: u64) -> bool {
        !((end1 < start2) || (start1 > end2))
    }

    // Determine the allele in the read at the locus of an INDEL variant
    pub fn indels_overlapping_variant(
        &self,
        read: &Rc<Record>,
        var_start: u64,
        var_end: u64,
    ) -> Vec<IndelEvent> {
        let mut read_pos: u32 = 0;
        let mut ref_pos = read.pos() as u32;
        let mut result = Vec::new();

        // See https://samtools.github.io/hts-specs/SAMv1.pdf page 8 for how CIGAR consumes
        for c in read.cigar().iter() {
            match *c {
                // Consumes both reference and query
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    ref_pos += len;
                    read_pos += len;
                }
                // Only consumes query
                Cigar::Ins(len) => {
                    let this_start = ref_pos as u64;
                    let this_end = this_start + len as u64 - 1;

                    if self.interval_overlaps(var_start, var_end, this_start, this_end) {
                        let seq_bytes = read.seq().as_bytes();
                        let inserted_bases = String::from_utf8(
                            seq_bytes[read_pos as usize .. (read_pos + len) as usize].to_vec()
                        )
                        .unwrap()
                        .to_ascii_uppercase();

                        result.push(IndelEvent {
                            indel_type: VarType::Ins,
                            start: this_start,
                            end: this_end,
                            bases: inserted_bases,
                        });
                    }

                    read_pos += len;
                }
                // Only consumes reference
                Cigar::Del(len) => {
                    let this_start = ref_pos as u64;
                    let this_end = this_start + len as u64 - 1;

                    if self.interval_overlaps(var_start, var_end, this_start, this_end) {
                        result.push(IndelEvent {
                            indel_type: VarType::Del,
                            start: this_start,
                            end: this_end,
                            bases: String::new(),
                        });
                    }

                    ref_pos += len;
                }
                // Only consumes reference
                Cigar::RefSkip(len) => {
                    ref_pos += len;
                }
                // Only consumes query
                Cigar::SoftClip(len) => {
                    read_pos += len;
                }
                // Consumes neither reference nor query
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }        

        result
    }

}

impl LocusFeaturesIndel {
    pub fn stats(&self) -> AlleleCountsIndelStats {
        let depth = self.common.all_read_features.num_reads;
        let ref_count = self.common.ref_read_features.num_reads;
        let alt_count = self.common.alt_read_features.num_reads;
        let other_count = depth - ref_count - alt_count;
        let alt_vaf = if depth > 0 {
            alt_count as f64 / depth as f64
        } else {
            0.0
        };

        AlleleCountsIndelStats {
            depth,
            ref_count,
            alt_count,
            other_count,
            alt_vaf,
            overlapping_indels_count: self.overlapping_indels_count,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct AlleleCountsIndelStats {
    pub depth: u32,
    pub ref_count: u32,
    pub alt_count: u32,
    pub other_count: u32,
    pub alt_vaf: f64,
    pub overlapping_indels_count: u64,
}

#[derive(Debug, Clone, Default)]
pub struct AlleleCountsSnv {
    pub a: u32,
    pub c: u32,
    pub g: u32,
    pub t: u32,
    pub n: u32,
}

impl AlleleCountsSnv {
    pub fn count(&mut self, base: char) {
        match base {
            'A' => self.a += 1,
            'C' => self.c += 1,
            'G' => self.g += 1,
            'T' => self.t += 1,
            'N' => self.n += 1,
            _ => eprintln!("Warning: Base does not match: {}", base),
        }
    }

    pub fn count_for_base(&self, base: char) -> u32 {
        match base {
            'A' => self.a,
            'C' => self.c,
            'G' => self.g,
            'T' => self.t,
            'N' => self.n,
            _ => 0,
        }
    }

    // fn depth(&self) -> u32 {
    //     self.a + self.c + self.g + self.t + self.n
    // }

    pub fn stats(&self, refr: char, alt: char) -> AlleleCountsSnvStats {
        let depth = self.a + self.c + self.g + self.t + self.n;
        let ref_count = self.count_for_base(refr);
        let alt_count = self.count_for_base(alt);
        let alt_vaf = if depth > 0 {
            alt_count as f64 / depth as f64
        } else {
            0.0
        };

        AlleleCountsSnvStats {
            depth,
            ref_count,
            alt_count,
            alt_vaf,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct AlleleCountsSnvStats {
    pub depth: u32,
    pub ref_count: u32,
    pub alt_count: u32,
    pub alt_vaf: f64,
}

#[derive(Debug, Clone, Default)]
pub struct ReadFeatures {
    pub nm: u32,
    pub base_qual: u32,
    pub map_qual: u32,
    pub align_len: u32,
    pub clipping: u32,
    pub indel: u32,
    pub forward_strand: u32,
    pub reverse_strand: u32,
    pub supplementary: u32,
    pub normalised_read_position: f64, 
    pub num_reads: u32,
}

impl ReadFeatures {
    pub fn count(
        &mut self,
        read: &Rc<Record>,
        query_pos: Option<u32>,
    ) {
        // Instead of counting num of reads, can create a function that sums foward and reverse strand?
        self.num_reads += 1;
        let query_len = read.seq_len() as usize;

        if let Some(qpos) = query_pos {
            if query_len > 0 {
                self.normalised_read_position += qpos as f64 / query_len as f64;
            }

            let qpos_usize = qpos as usize;
            if qpos_usize < read.qual().len() {
                let pos_qual = read.qual()[qpos_usize];
                self.base_qual += pos_qual as u32;
            }
        }

        self.align_len += self.query_alignment_length(&read) as u32;
        self.map_qual += read.mapq() as u32;

        for c in read.cigar().iter() {
            match *c {
                Cigar::Ins(len) | Cigar::Del(len) => self.indel += len as u32,
                Cigar::SoftClip(len) | Cigar::HardClip(len) => self.clipping += len as u32,
                _ => {}
            }
        }

        if let Ok(aux) = read.aux(b"NM") {
            match aux {
                Aux::I8(v) => self.nm += v as u32,
                Aux::U8(v) => self.nm += v as u32,
                Aux::I16(v) => self.nm += v as u32,
                Aux::U16(v) => self.nm += v as u32,
                Aux::I32(v) => self.nm += v as u32,
                Aux::U32(v) => self.nm += v as u32,
                _ => {}
            }
        }

        if read.is_reverse() {
            self.reverse_strand += 1;
        } else {
            self.forward_strand += 1;
        }
        if read.is_supplementary() {
            self.supplementary += 1;
        }

    }

    pub fn query_alignment_length(&self, record: &Rc<Record>) -> u32 {
        let mut len = 0;
        for c in record.cigar().iter() {
            match *c {
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | Cigar::Ins(l) => len += l,
                _ => {}
            }
        }
        len
    }

    pub fn normalized(&self) -> NormalizedReadFeatures {
        if self.num_reads > 0 {
            let n = self.num_reads as f64;
            NormalizedReadFeatures {
                nm: Some(self.nm as f64 / n),
                base_qual: Some(self.base_qual as f64 / n),
                map_qual: Some(self.map_qual as f64 / n),
                align_len: Some(self.align_len as f64 / n),
                clipping: Some(self.clipping as f64 / n),
                indel: Some(self.indel as f64 / n),
                forward_strand: Some(self.forward_strand as f64 / n),
                reverse_strand: Some(self.reverse_strand as f64 / n),
                supplementary: Some(self.supplementary as f64 / n),
                normalised_read_position: Some(self.normalised_read_position / n),
            }
        } else {
            NormalizedReadFeatures::default()
        }
    }
}

#[derive(Debug, Clone, Default, serde::Serialize)]
pub struct NormalizedLocusFeaturesRow {
    pub ref_nm: Option<f64>,
    pub ref_base_qual: Option<f64>,
    pub ref_map_qual: Option<f64>,
    pub ref_align_len: Option<f64>,
    pub ref_clipping: Option<f64>,
    pub ref_indel: Option<f64>,
    pub ref_forward_strand: Option<f64>,
    pub ref_reverse_strand: Option<f64>,
    pub ref_supplementary: Option<f64>,
    pub ref_normalised_read_position: Option<f64>,

    pub alt_nm: Option<f64>,
    pub alt_base_qual: Option<f64>,
    pub alt_map_qual: Option<f64>,
    pub alt_align_len: Option<f64>,
    pub alt_clipping: Option<f64>,
    pub alt_indel: Option<f64>,
    pub alt_forward_strand: Option<f64>,
    pub alt_reverse_strand: Option<f64>,
    pub alt_supplementary: Option<f64>,
    pub alt_normalised_read_position: Option<f64>,

    pub all_nm: Option<f64>,
    pub all_base_qual: Option<f64>,
    pub all_map_qual: Option<f64>,
    pub all_align_len: Option<f64>,
    pub all_clipping: Option<f64>,
    pub all_indel: Option<f64>,
    pub all_forward_strand: Option<f64>,
    pub all_reverse_strand: Option<f64>,
    pub all_supplementary: Option<f64>,
    pub all_normalised_read_position: Option<f64>,
}

// Use option or not? Because in python script if there are no values it defaults to ''
// Compared to default which would initialize it as 0.0
#[derive(Debug, Clone, Copy, Default, serde::Serialize)]
pub struct NormalizedReadFeatures {
    pub nm: Option<f64>,
    pub base_qual: Option<f64>,
    pub map_qual: Option<f64>,
    pub align_len: Option<f64>,
    pub clipping: Option<f64>,
    pub indel: Option<f64>,
    pub forward_strand: Option<f64>,
    pub reverse_strand: Option<f64>,
    pub supplementary: Option<f64>,
    pub normalised_read_position: Option<f64>,
}