use std::collections::HashMap;
use std::io::{BufRead, Write};

use timsrust::converters::ConvertableDomain;
use timsrust::readers::{FrameReader, MetadataReader};

const AA: [(u8, f64); 20] = [
    (b'G', 57.02146),
    (b'A', 71.03711),
    (b'S', 87.03203),
    (b'P', 97.05276),
    (b'V', 99.06841),
    (b'T', 101.04768),
    (b'C', 103.00919),
    (b'L', 113.08406),
    (b'I', 113.08406),
    (b'N', 114.04293),
    (b'D', 115.02694),
    (b'Q', 128.05858),
    (b'K', 128.09496),
    (b'E', 129.04259),
    (b'M', 131.04049),
    (b'H', 137.05891),
    (b'F', 147.06841),
    (b'R', 156.10111),
    (b'Y', 163.06333),
    (b'W', 186.07931),
];
const P: f64 = 1.007_276_466_88;
const H2O: f64 = 18.010_564_686_3;

fn residues(pep: &str) -> Vec<f64> {
    let m: HashMap<u8, f64> = AA.into_iter().collect();
    let core = pep.split('.').nth(1).unwrap_or(pep);
    let b = core.as_bytes();
    let mut r: Vec<f64> = Vec::new();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'+' || b[i] == b'-' {
            let j = (i + 1..b.len())
                .find(|&k| !(b[k].is_ascii_digit() || b[k] == b'.'))
                .unwrap_or(b.len());
            let d: f64 = core[i..j].parse().unwrap_or(0.0);
            if let Some(l) = r.last_mut() {
                *l += d;
            }
            i = j;
            continue;
        }
        r.push(m.get(&b[i]).copied().unwrap_or(0.0));
        i += 1;
    }
    r
}

fn ions(r: &[f64], z: u8) -> Vec<f64> {
    let tot: f64 = r.iter().sum();
    let mut out = Vec::new();
    let mut pre = 0.0;
    for k in 0..r.len().saturating_sub(1) {
        pre += r[k];
        for c in 1..=(if z >= 3 { 2 } else { 1 }) {
            let cf = c as f64;
            if k + 1 >= 3 {
                out.push((pre + cf * P) / cf);
            }
            if r.len() - k - 1 >= 3 {
                out.push((tot - pre + H2O + cf * P) / cf);
            }
        }
    }
    out
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (d, rows, out) = (&a[1], &a[2], &a[3]);
    let frames = FrameReader::new(d).expect("frames");
    let meta = MetadataReader::new(d).expect("metadata");
    let mut pos: HashMap<usize, usize> = HashMap::new();
    for k in 0..frames.len() {
        pos.insert(frames.get_frame_without_coordinates(k).unwrap().index, k);
    }
    let mzc = meta.mz_converter;
    let mut req: HashMap<usize, Vec<(String, usize, usize, f64, Vec<f64>)>> = HashMap::new();
    for (n, line) in std::io::BufReader::new(std::fs::File::open(rows).unwrap())
        .lines()
        .enumerate()
    {
        let line = line.unwrap();
        if n == 0 {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let (id, frame, lo, hi, pmz, z, pep) = (
            f[0].to_string(),
            f[1].parse::<usize>().unwrap(),
            f[2].parse::<usize>().unwrap(),
            f[3].parse::<usize>().unwrap(),
            f[4].parse::<f64>().unwrap(),
            f[5].parse::<u8>().unwrap(),
            f[6],
        );
        let ion = ions(&residues(pep), z);
        req.entry(frame).or_default().push((id, lo, hi, pmz, ion));
    }
    let mut o = std::io::BufWriter::new(std::fs::File::create(out).unwrap());
    writeln!(o, "SpecId\tn_ions_seen\tim_spread\tim_top_frac").unwrap();
    let ppm = 15e-6;
    for (frame, list) in req {
        let Some(&p) = pos.get(&frame) else { continue };
        let fr = frames.get(p).unwrap();
        let nscan = fr.scan_offsets.len().saturating_sub(1);
        for (id, lo, hi, _pmz, ion) in list {
            let mut cents: Vec<(f64, f64)> = Vec::new();
            for &mz in &ion {
                let (mut w, mut s) = (0.0f64, 0.0f64);
                for scan in lo..hi.min(nscan) {
                    for k in fr.scan_offsets[scan]..fr.scan_offsets[scan + 1] {
                        let m = mzc.convert(fr.tof_indices[k] as f64);
                        if (m - mz).abs() <= mz * ppm {
                            let it = fr.intensities[k] as f64;
                            w += it;
                            s += it * scan as f64;
                        }
                    }
                }
                if w > 0.0 {
                    cents.push((s / w, w));
                }
            }
            let (spread, top) = if cents.len() >= 3 {
                let tw: f64 = cents.iter().map(|c| c.1).sum();
                let mean = cents.iter().map(|c| c.0 * c.1).sum::<f64>() / tw;
                let sd = (cents.iter().map(|c| c.1 * (c.0 - mean).powi(2)).sum::<f64>() / tw).sqrt();
                let near = cents
                    .iter()
                    .filter(|c| (c.0 - mean).abs() <= 10.0)
                    .map(|c| c.1)
                    .sum::<f64>()
                    / tw;
                (sd, near)
            } else {
                (f64::NAN, f64::NAN)
            };
            writeln!(o, "{id}\t{}\t{spread:.2}\t{top:.3}", cents.len()).unwrap();
        }
    }
}
