use std::collections::HashMap;
use std::io::{BufRead, Write};

const AA: [(u8, f64); 20] = [
    (b'G', 57.02146),
    (b'A', 71.03711),
    (b'S', 87.03203),
    (b'P', 97.05276),
    (b'V', 99.06841),
    (b'T', 101.04768),
    (b'C', 160.03065),
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

fn ions(seq: &str) -> Vec<f64> {
    let m: HashMap<u8, f64> = AA.into_iter().collect();
    let mut r: Vec<f64> = Vec::new();
    let s = seq.as_bytes();
    let mut i = 0;
    while i < s.len() {
        if s[i] == b'(' {
            let j = seq[i..].find(')').map_or(s.len(), |k| i + k + 1);
            if seq[i..j].contains(":35") {
                if let Some(l) = r.last_mut() {
                    *l += 15.994915;
                }
            }
            i = j;
            continue;
        }
        r.push(m.get(&s[i]).copied().unwrap_or(0.0));
        i += 1;
    }
    let tot: f64 = r.iter().sum();
    let mut out = Vec::new();
    let mut pre = 0.0;
    for k in 0..r.len().saturating_sub(1) {
        pre += r[k];
        if k + 1 >= 3 {
            out.push(pre + P);
        }
        if r.len() - k - 1 >= 3 {
            out.push(tot - pre + H2O + P);
        }
    }
    out
}

fn hits(peaks: &[(f64, f32)], mz: &[f64], ppm: f64) -> usize {
    mz.iter()
        .filter(|&&x| {
            let t = x * ppm * 1e-6;
            let i = peaks.partition_point(|p| p.0 < x - t);
            peaks.get(i).is_some_and(|p| p.0 <= x + t)
        })
        .count()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (d, targets, out) = (&a[1], &a[2], &a[3]);
    let reader = input::TimsTofReader::open(d).expect("open .d");
    let spectra: Vec<input::Spectrum> = reader.filter_map(Result::ok).collect();
    eprintln!("{} spectra", spectra.len());
    let mut o = std::io::BufWriter::new(std::fs::File::create(out).unwrap());
    writeln!(o, "key\tn_ions\tbest_win_hits7\tbest_win_hits15\tbest_win_hits30\tn_spectra\tmean_peaks").unwrap();
    for (n, line) in std::io::BufReader::new(std::fs::File::open(targets).unwrap())
        .lines()
        .enumerate()
    {
        let line = line.unwrap();
        if n == 0 {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let (key, seq, mz, rt) = (f[0], f[1], f[2].parse::<f64>().unwrap(), f[3].parse::<f64>().unwrap() * 60.0);
        let ion = ions(seq);
        let near: Vec<&input::Spectrum> = spectra
            .iter()
            .filter(|s| {
                let c = s.isolation_target_mz.unwrap_or(0.0);
                let h = s.isolation_lower_offset.unwrap_or(0.0);
                (mz - c).abs() <= h && s.rt_seconds.is_some_and(|t| (t - rt).abs() <= 3.0)
            })
            .collect();
        let best = |ppm| near.iter().map(|s| hits(&s.peaks, &ion, ppm)).max().unwrap_or(0);
        let mp = near.iter().map(|s| s.peaks.len()).sum::<usize>() as f64 / near.len().max(1) as f64;
        writeln!(o, "{key}\t{}\t{}\t{}\t{}\t{}\t{mp:.0}", ion.len(), best(7.0), best(15.0), best(30.0), near.len()).unwrap();
    }
}
