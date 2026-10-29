# Two-sided rigorous optimality bounds for the PGLib-OPF v23.07 typical cases

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23752138.svg)](https://doi.org/10.5281/zenodo.23752138)

**Status: computer-assisted proof candidate, not peer reviewed.**

Author: Puyi Sun (ORCID [0009-0004-1644-5300](https://orcid.org/0009-0004-1644-5300)).
AI systems were used extensively, for writing code, running computations and reviewing; see [AI_DISCLOSURE.md](AI_DISCLOSURE.md).

## The claim

For each of the 26 typical-operating-condition (TYP) cases of PGLib-OPF v23.07, with the AC optimal power flow model of
that release and the exact data of its case files, this repository contains two certificates:

- an **upper certificate**: a small box proven (Krawczyk interval-Newton test, interval arithmetic with outward rounding)
  to contain an exact solution of the AC power-flow equations at which every operational limit holds, so **OPT ≤ U**;
- a **lower certificate**: a sum-of-squares / Lagrangian dual certificate, checked in interval arithmetic with every PSD
  correction and every residual priced, so **OPT ≥ L**.

OPT is the global optimum of the nonconvex AC-OPF problem. 13 of the 26 cases have a relative gap of at most 1e-6, i.e. the
global optimum is pinned to about six significant digits.

| case | buses | L (lower bound, $/h) | U (upper bound, $/h) | relative gap (U−L)/U | lower-bound certificate |
|---|---:|---:|---:|---:|---|
| case3_lmbd | 3 | 5812.642962 | 5812.643005 | 7.33e-09 | sparse SOS, mixed order 1/2 |
| case5_pjm | 5 | 17551.890711 | 17551.891094 | 2.17e-08 | sparse SOS, mixed order 1/2 |
| case14_ieee | 14 | 2178.080425 | 2178.080449 | 1.02e-08 | chordal Shor SDP (Lagrangian dual) |
| case24_ieee_rts | 24 | 63352.202444 | 63352.203674 | 1.94e-08 | chordal Shor SDP (Lagrangian dual) |
| case30_as | 30 | 803.127218 | 803.127313 | 1.17e-07 | chordal Shor SDP (Lagrangian dual) |
| case30_ieee | 30 | 8208.515455 | 8208.515549 | 1.13e-08 | chordal Shor SDP (Lagrangian dual) |
| case39_epri | 39 | 138415.554411 | 138415.563671 | 6.69e-08 | sparse SOS, mixed order 1/2 |
| case57_ieee | 57 | 37589.336873 | 37589.338384 | 4.02e-08 | sparse SOS, mixed order 1/2 |
| case60_c | 60 | 92693.623915 | 92693.670117 | 4.98e-07 | sparse SOS, mixed order 1/2 |
| case73_ieee_rts | 73 | 189764.081348 | 189764.084939 | 1.89e-08 | chordal Shor SDP (Lagrangian dual) |
| case89_pegase | 89 | 106968.658190 | 107285.675532 | 2.95e-03 | chordal Shor SDP (Lagrangian dual) |
| case118_ieee | 118 | 97213.585027 | 97213.607809 | 2.34e-07 | sparse SOS, mixed order 1/2 |
| case162_ieee_dtc | 162 | 106512.858738 | 108075.648889 | 1.45e-02 | sparse SOS, mixed order 1/2 |
| case179_goc | 179 | 754128.038771 | 754266.421037 | 1.83e-04 | sparse SOS, mixed order 1/2 |
| case197_snem | 197 | 1.501257 | 1.501747 | 3.25e-04 | chordal Shor SDP (Lagrangian dual) |
| case200_activ | 200 | 27557.570830 | 27557.571313 | 1.75e-08 | chordal Shor SDP (Lagrangian dual) |
| case240_pserc | 240 | 3283231.188888 | 3329670.118665 | 1.39e-02 | sparse SOS, mixed order 1/2 |
| case300_ieee | 300 | 565219.433978 | 565220.008801 | 1.02e-06 | sparse SOS, mixed order 1/2 |
| case500_goc | 500 | 454945.910696 | 454945.987308 | 1.68e-07 | sparse SOS, mixed order 1/2 |
| case588_sdet | 588 | 311343.651681 | 313139.784598 | 5.74e-03 | sparse SOS, mixed order 1/2 |
| case793_goc | 793 | 258581.712068 | 260197.852143 | 6.21e-03 | sparse SOS, mixed order 1/2 |
| case1354_pegase | 1354 | 1257170.185668 | 1258844.004024 | 1.33e-03 | sparse SOS, mixed order 1/2 |
| case1803_snem | 1803 | 94825.222708 | 98334.547531 | 3.57e-02 | chordal Shor SDP (Lagrangian dual) |
| case1888_rte | 1888 | 1375653.865114 | 1402531.010635 | 1.92e-02 | sparse SOS, mixed order 1/2 |
| case1951_rte | 1951 | 2085565.057043 | 2085581.868859 | 8.06e-06 | sparse SOS, mixed order 1/2 |
| case2000_goc | 2000 | 973063.826000 | 973432.478973 | 3.79e-04 | sparse SOS, mixed order 1/2 |

L is rounded down and U up to six decimals for display; the exact binary64 endpoints are the `L` / `U` fields of the
certificates (also in `cases.json`), and a checker accepts a certificate only if it proves the declared value. The gap
column is computed from the full-precision declared endpoints, so for the tightest cases it can differ from a gap
recomputed from the six-decimal display. Model, certificate formats and checkers are
described in [METHOD.md](METHOD.md).

## What is not claimed

- No external peer review has taken place.
- No priority claim. Certified (rigorous) SDP **lower** bounds for AC-OPF were published by Oustry, D'Ambrosio, Liberti and
  Ruiz (EPSR 212, 2022, with numerical upper bounds); global-optimality gaps from convex relaxations and bound tightening by
  Gopinath et al. (EPSR 189, 2020) and by sparse moment-SOS hierarchies (Wang, Magron, Lasserre: CS-TSSOS) were reported with
  floating-point solvers; interval verification of power-flow solutions is classical. What this repository adds is a fixed,
  machine-replayable artifact in which **both** ends of the interval are verified with directed rounding, on all 26 TYP
  cases up to 2000 buses, re-checked by checkers written separately from the certificate generators.
- The statements are about the PGLib v23.07 TYP model exactly as written (apparent-power thermal limits, angle-difference
  limits). They do not cover other variants (e.g. current-magnitude limits), other PGLib versions or real-world operation.
- Several gaps are not small: case1803_snem 3.6 %, case1888_rte 1.9 %, case162_ieee_dtc 1.4 %, case240_pserc 1.4 %.

## How to verify

```
python3 -m pip install numpy scipy mpmath
python3 replay.py case14_ieee            # one case, a few seconds
python3 replay.py --jobs 8               # all 52 endpoints; the largest upper-bound checks take 20-40 min each
```

`replay.py` first checks every file against `MANIFEST.sha256`, then runs one checker per endpoint and writes
`replay_results.csv`. A single endpoint by hand (exit code 0 = accepted):

```
python3 checkers/independent_ub_krawczyk_r12.py cases/pglib_opf_case14_ieee.m certificates/upper/case14_ieee_upper.json
python3 checkers/independent_chordal_r13.py cases/pglib_opf_case14_ieee.m certificates/lower/case14_ieee_lower_chordal.json
```

Each checker reads the raw `.m` file, rebuilds the model itself, compares the SHA-256 of the case file with the one stored
in the certificate, and imports nothing from the programs that produced the certificates.

Provenance note: the three files in `checkers/` are byte-identical copies of the checker versions that produced the
recorded acceptances (their SHA-256 are in `MANIFEST.sha256`). Their docstrings are kept unedited and mention the private
research repository in which they were developed and, in one case, an earlier state of testing; those remarks are historical.

## Files

| path | content |
|---|---|
| `cases/` | the 26 PGLib-OPF v23.07 TYP case files, unmodified |
| `certificates/upper/`, `certificates/lower/` | one upper and one lower certificate per case (JSON, binary64 values as hex floats) |
| `checkers/` | independent checkers (Python; numpy, scipy, mpmath) |
| `cases.json` | which certificate and checker belong to which case, declared L and U |
| `replay.py` | runs all checks |
| `replay_records/` | a complete replay of all 52 endpoints on a second machine (the author's computer) |
| `paper/` | short paper (PDF and LaTeX source) |
| `MANIFEST.sha256` | SHA-256 of every file |

## Licence

Code: MIT ([LICENSE](LICENSE)). Certificates, tables and text: CC BY 4.0. PGLib-OPF case files: CC BY 4.0 by the PGLib-OPF
authors, unmodified copies of tag v23.07. Details in [LICENSE-DATA.md](LICENSE-DATA.md).

## Citation

Archived on Zenodo: https://doi.org/10.5281/zenodo.23752138 (v1.0.1; mathematical content identical to v1.0). See [CITATION.cff](CITATION.cff). Data source: S. Babaeinejadsarookolaee et al., *The Power Grid Library for Benchmarking AC
Optimal Power Flow Algorithms*, arXiv:1908.02788.
