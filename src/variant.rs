// use std::fs::File;
// use std::io::{BufRead, BufReader};
// use std::collections::VecDeque;
use clap::ValueEnum;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum VarClass {
    Snv,
    Indel,
}
