#!/usr/bin/env python3
"""Re-check the certificates of this repository with the independent checkers in checkers/.

usage:  python3 replay.py [CASE ...] [--jobs N] [--only {U,L}]
Needs Python >= 3.9 with numpy, scipy and mpmath.  Steps:
  1. every file listed in MANIFEST.sha256 must match its SHA-256;
  2. cases.json must list exactly the 26 cases, and for every selected case the declared L / U and the case-file SHA-256
     recorded in cases.json must equal the values stored in the certificates themselves;
  3. one checker process per endpoint (U: Krawczyk upper-bound checker; L: the lower-bound checker named in cases.json).
An endpoint is ACCEPTED iff the checker exits with code 0 (lower-bound checkers must also print a line starting with ACCEPT).
Writes replay_results.csv and replay_logs/.  Exit status: 0 if every requested endpoint is accepted, 1 if any is not,
2 for a usage or consistency error (unknown case, bad option, manifest or metadata mismatch)."""
import argparse, csv, concurrent.futures as cf, hashlib, json, os, re, subprocess, sys, time

here = os.path.dirname(os.path.abspath(__file__)); os.chdir(here)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
ap = argparse.ArgumentParser(description='Replay the TOP-035 certificates with the independent checkers.')
ap.add_argument('cases', nargs='*', help='case names (default: all)')
ap.add_argument('--jobs', type=int, default=1, help='parallel checker processes (default 1)')
ap.add_argument('--only', choices=['U', 'L'], help='check only upper (U) or only lower (L) certificates')
a = ap.parse_args()
if a.jobs < 1: ap.error('--jobs must be >= 1')


def fail(msg): print('ERROR:', msg); sys.exit(2)


bad = []
for line in open('MANIFEST.sha256'):
    h, p = line.rstrip('\n').split('  ', 1)
    if not os.path.exists(p) or sha(p) != h: bad.append(p)
if bad: fail(f'MANIFEST mismatch: {bad}')
print('MANIFEST.sha256: all files match')

allc = json.load(open('cases.json'))
if len(allc) != 26 or len({c['case'] for c in allc}) != 26: fail('cases.json must list exactly 26 distinct cases')
known = {c['case']: c for c in allc}
unknown = [n for n in a.cases if n not in known]
if unknown: fail(f'unknown case(s): {unknown}; known: {sorted(known)}')
cases = [known[n] for n in a.cases] if a.cases else allc
for c in cases:                                   # metadata must agree with the certificates
    U = json.load(open(c['upper_certificate'])); L = json.load(open(c['lower_certificate'])); h = sha(c['case_file'])
    if U['U'] != c['U_declared'] or L['L'] != c['L_declared']: fail(f"{c['case']}: declared values in cases.json differ from the certificates")
    if U['sha256'] != h or L['sha256'] != h: fail(f"{c['case']}: case-file SHA-256 differs from the one recorded in a certificate")
print(f'cases.json: {len(cases)} case(s) consistent with their certificates')

PAT = {'independent_ub_krawczyk_r12.py': r'^U (\S+) declared', 'independent_sparse.py': r'L_outward (\S+)',
       'independent_chordal_r13.py': r'^RIGOROUS_L (\S+)'}
work = []
for c in cases:
    if a.only in (None, 'U'): work.append((c['case'], 'U', c['upper_checker'], c['upper_certificate'], c['U_declared'], c['case_file']))
    if a.only in (None, 'L'): work.append((c['case'], 'L', c['lower_checker'], c['lower_certificate'], c['L_declared'], c['case_file']))
os.makedirs('replay_logs', exist_ok=True)


def run(w):
    case, end, ck, cert, declared, mfile = w; t = time.time()
    p = subprocess.run([sys.executable, '-I', '-u', ck, mfile, cert], capture_output=True, text=True)
    out = p.stdout + p.stderr; open(f'replay_logs/{case}_{end}.log', 'w').write(out)
    m = re.search(PAT[os.path.basename(ck)], out, re.M)
    ok = p.returncode == 0 and (end == 'U' or re.search(r'^ACCEPT', out, re.M) is not None)
    return dict(case=case, endpoint=end, declared=repr(declared), checker_value=m.group(1) if m else '',
                verdict='ACCEPT' if ok else 'NOT ACCEPTED', exit_code=p.returncode, seconds=round(time.time() - t, 1))


res = []
with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
    for r in ex.map(run, sorted(work, key=lambda w: -os.path.getsize(w[5]))):
        res.append(r)
        print(f"{r['case']:18s} {r['endpoint']}  {r['verdict']:12s} declared {r['declared']:>24s}  checker {r['checker_value']:>24s}  [{r['seconds']}s]", flush=True)
with open('replay_results.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
n_ok = sum(r['verdict'] == 'ACCEPT' for r in res)
print(f'{n_ok}/{len(res)} endpoints accepted; details in replay_results.csv and replay_logs/')
sys.exit(0 if n_ok == len(res) else 1)
