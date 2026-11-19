use std::collections::HashMap;
use std::io::{BufRead, Write};

use rayon::prelude::*;
use timsrust::converters::ConvertableDomain;
use timsrust::readers::{FrameReader, MetadataReader};
use timsrust::MSLevel;

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
const C13: f64 = 1.003_354_835;
const K: usize = 5;
const BIN: usize = 8;
const TOP: usize = 8;
const PPM: f64 = 15e-6;
const COLS: [&str; 13] = [
    "pg_corr",
    "pg_clean",
    "pg_coapex",
    "pg_peak_frac",
    "pg_n_apex",
    "pg_log_int",
    "pg_rt_off",
    "pg_im_corr",
    "pg_rt_corr",
    "pg_ms1_corr",
    "pg_iso_corr",
    "pg_ms1_log",
    "pg_apex_scan",
];

struct Slab {
    tof: Vec<u32>,
    scan: Vec<u16>,
    int: Vec<u32>,
}

fn slab(scan_offsets: &[usize], tof: &[u32], int: &[u32]) -> Slab {
    let mut v: Vec<(u32, u16, u32)> = Vec::with_capacity(tof.len());
    for s in 0..scan_offsets.len().saturating_sub(1) {
        for k in scan_offsets[s]..scan_offsets[s + 1] {
            v.push((tof[k], s as u16, int[k]));
        }
    }
    v.sort_unstable_by_key(|p| p.0);
    Slab {
        tof: v.iter().map(|p| p.0).collect(),
        scan: v.iter().map(|p| p.1).collect(),
        int: v.iter().map(|p| p.2).collect(),
    }
}

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

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if a.len() < 3 {
        return 0.0;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma) * (x - ma);
        sbb += (y - mb) * (y - mb);
    }
    if saa <= 0.0 || sbb <= 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}

fn smooth(m: &[f64], t: usize, b: usize) -> Vec<f64> {
    let mut o = vec![0.0; m.len()];
    for i in 0..t {
        for j in 0..b {
            let mut s = 0.0;
            for di in i.saturating_sub(1)..(i + 2).min(t) {
                for dj in j.saturating_sub(1)..(j + 2).min(b) {
                    s += m[di * b + dj];
                }
            }
            o[i * b + j] = s;
        }
    }
    o
}

fn argmax(m: &[f64]) -> usize {
    (0..m.len())
        .max_by(|&x, &y| m[x].total_cmp(&m[y]))
        .unwrap_or(0)
}

struct Req {
    id: String,
    frame: usize,
    lo: usize,
    hi: usize,
    mz: f64,
    z: u8,
    ions: Vec<f64>,
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (d, rows, out) = (&a[1], &a[2], &a[3]);
    let frames = FrameReader::new(d).expect("frames");
    let meta = MetadataReader::new(d).expect("metadata");
    let mzc = meta.mz_converter;
    let n = frames.len();
    let t0 = std::time::Instant::now();
    let loaded: Vec<(usize, f64, bool, u8, Slab)> = (0..n)
        .into_par_iter()
        .map(|k| {
            let f = frames.get(k).expect("frame");
            let s = slab(&f.scan_offsets, &f.tof_indices, &f.intensities);
            (f.index, f.rt_in_seconds, f.ms_level == MSLevel::MS1, f.window_group, s)
        })
        .collect();
    eprintln!("loaded {n} frames in {:.1}s", t0.elapsed().as_secs_f64());
    let mut pos: HashMap<usize, usize> = HashMap::new();
    let mut groups: HashMap<u8, Vec<usize>> = HashMap::new();
    let mut parent: Vec<Option<usize>> = vec![None; n];
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| loaded[x].1.total_cmp(&loaded[y].1));
    let mut last = None;
    for &k in &order {
        pos.insert(loaded[k].0, k);
        if loaded[k].2 {
            last = Some(k);
        } else {
            parent[k] = last;
            groups.entry(loaded[k].3).or_default().push(k);
        }
    }
    let mut slot: HashMap<usize, (u8, usize)> = HashMap::new();
    for (g, v) in &groups {
        for (p, &k) in v.iter().enumerate() {
            slot.insert(k, (*g, p));
        }
    }
    let mut reqs: Vec<Req> = Vec::new();
    for (i, line) in std::io::BufReader::new(std::fs::File::open(rows).unwrap())
        .lines()
        .enumerate()
    {
        let line = line.unwrap();
        if i == 0 {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let z = f[5].parse::<u8>().unwrap();
        reqs.push(Req {
            id: f[0].to_string(),
            frame: f[1].parse().unwrap(),
            lo: f[2].parse().unwrap(),
            hi: f[3].parse().unwrap(),
            mz: f[4].parse().unwrap(),
            z,
            ions: ions(&residues(f[6]), z),
        });
    }
    let tofs = |mz: f64| -> (u32, u32) {
        (
            mzc.invert(mz * (1.0 - PPM)).floor().max(0.0) as u32,
            mzc.invert(mz * (1.0 + PPM)).ceil() as u32,
        )
    };
    let xim = |s: &Slab, mz: f64, lo: usize, hi: usize, nb: usize, row: &mut [f64]| {
        let (a, b) = tofs(mz);
        let i0 = s.tof.partition_point(|&t| t < a);
        for k in i0..s.tof.len() {
            if s.tof[k] > b {
                break;
            }
            let sc = s.scan[k] as usize;
            if sc >= lo && sc < hi {
                row[((sc - lo) / BIN).min(nb - 1)] += s.int[k] as f64;
            }
        }
    };
    let t1 = std::time::Instant::now();
    let feats: Vec<[f64; COLS.len()]> = reqs
        .par_iter()
        .map(|r| {
            let mut f = [0.0; COLS.len()];
            let Some(&k) = pos.get(&r.frame) else {
                return f;
            };
            let Some(&(g, p)) = slot.get(&k) else {
                return f;
            };
            let win = &groups[&g];
            let (c0, c1) = (p.saturating_sub(K), (p + K + 1).min(win.len()));
            let t = c1 - c0;
            let nb = (r.hi - r.lo).div_ceil(BIN).max(1);
            let mats: Vec<Vec<f64>> = r
                .ions
                .iter()
                .map(|&mz| {
                    let mut m = vec![0.0; t * nb];
                    for c in 0..t {
                        let s = &loaded[win[c0 + c]].4;
                        xim(s, mz, r.lo, r.hi, nb, &mut m[c * nb..(c + 1) * nb]);
                    }
                    smooth(&m, t, nb)
                })
                .collect();
            let tot: Vec<f64> = mats.iter().map(|m| m.iter().sum()).collect();
            let mut top: Vec<usize> = (0..mats.len()).filter(|&i| tot[i] > 0.0).collect();
            top.sort_by(|&x, &y| tot[y].total_cmp(&tot[x]));
            top.truncate(TOP);
            if top.len() < 3 {
                f[6] = K as f64;
                return f;
            }
            let norm: Vec<Vec<f64>> = top
                .iter()
                .map(|&i| {
                    let mx = mats[i].iter().cloned().fold(0.0, f64::max);
                    mats[i].iter().map(|v| v / mx).collect()
                })
                .collect();
            let cons: Vec<f64> = (0..t * nb)
                .map(|x| {
                    let mut v: Vec<f64> = norm.iter().map(|m| m[x]).collect();
                    v.sort_by(f64::total_cmp);
                    let n = v.len();
                    if n % 2 == 1 {
                        v[n / 2]
                    } else {
                        (v[n / 2 - 1] + v[n / 2]) / 2.0
                    }
                })
                .collect();
            let ap = argmax(&cons);
            let (at, ab) = (ap / nb, ap % nb);
            let near = |x: usize, dt: usize, db: usize| {
                (x / nb).abs_diff(at) <= dt && (x % nb).abs_diff(ab) <= db
            };
            let corrs: Vec<f64> = top.iter().map(|&i| pearson(&mats[i], &cons)).collect();
            f[0] = corrs.iter().sum::<f64>() / corrs.len() as f64;
            f[1] = corrs.iter().filter(|&&c| c >= 0.6).count() as f64 / corrs.len() as f64;
            f[2] = top
                .iter()
                .filter(|&&i| near(argmax(&mats[i]), 1, 2))
                .count() as f64
                / top.len() as f64;
            let (mut inside, mut all) = (0.0, 0.0);
            for &i in &top {
                for (x, v) in mats[i].iter().enumerate() {
                    all += v;
                    if near(x, 1, 3) {
                        inside += v;
                    }
                }
            }
            f[3] = inside / all.max(1.0);
            f[4] = mats
                .iter()
                .filter(|m| m.iter().enumerate().any(|(x, v)| *v > 0.0 && near(x, 1, 2)))
                .count() as f64
                / mats.len().max(1) as f64;
            f[5] = (1.0 + inside).ln();
            f[6] = (c0 + at).abs_diff(p) as f64;
            let im_prof = |m: &[f64]| -> Vec<f64> {
                (0..nb)
                    .map(|j| {
                        (at.saturating_sub(1)..(at + 2).min(t))
                            .map(|i| m[i * nb + j])
                            .sum()
                    })
                    .collect()
            };
            let rt_prof = |m: &[f64]| -> Vec<f64> {
                (0..t)
                    .map(|i| {
                        (ab.saturating_sub(2)..(ab + 3).min(nb))
                            .map(|j| m[i * nb + j])
                            .sum()
                    })
                    .collect()
            };
            let (ci, cr) = (im_prof(&cons), rt_prof(&cons));
            f[7] = top
                .iter()
                .map(|&i| pearson(&im_prof(&mats[i]), &ci))
                .sum::<f64>()
                / top.len() as f64;
            f[8] = top
                .iter()
                .map(|&i| pearson(&rt_prof(&mats[i]), &cr))
                .sum::<f64>()
                / top.len() as f64;
            let ms1 = |mz: f64| -> Vec<f64> {
                let mut m = vec![0.0; t * nb];
                for c in 0..t {
                    if let Some(q) = parent[win[c0 + c]] {
                        xim(&loaded[q].4, mz, r.lo, r.hi, nb, &mut m[c * nb..(c + 1) * nb]);
                    }
                }
                smooth(&m, t, nb)
            };
            let m0 = ms1(r.mz);
            let m1 = ms1(r.mz + C13 / r.z as f64);
            f[9] = pearson(&m0, &cons);
            f[10] = pearson(&m0, &m1);
            f[11] = (1.0
                + m0.iter()
                    .enumerate()
                    .filter(|(x, _)| near(*x, 1, 3))
                    .map(|(_, v)| v)
                    .sum::<f64>())
            .ln();
            f[12] = (r.lo + ab * BIN + BIN / 2) as f64;
            f
        })
        .collect();
    eprintln!(
        "scored {} rows in {:.1}s",
        reqs.len(),
        t1.elapsed().as_secs_f64()
    );
    let mut o = std::io::BufWriter::new(std::fs::File::create(out).unwrap());
    writeln!(o, "SpecId\t{}", COLS.join("\t")).unwrap();
    for (r, f) in reqs.iter().zip(&feats) {
        let v: Vec<String> = f.iter().map(|x| format!("{x:.4}")).collect();
        writeln!(o, "{}\t{}", r.id, v.join("\t")).unwrap();
    }
}
