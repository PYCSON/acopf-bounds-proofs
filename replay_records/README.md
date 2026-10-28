# Replay records

`author_computer_2026-10-28.csv`: all 52 endpoints re-checked on the author's own computer (WSL2 Ubuntu,
Intel i9-14900HX, Python 3.14, numpy 2.5.4, scipy 1.18.1, mpmath 1.3.0), 2026-10-28, with checker, certificate and case
files byte-identical to those in this repository (SHA-256 compared). `checker_value` is the bound the checker proved
(upper: ≤ declared U; lower: ≥ declared L); `log_sha256` identifies the full checker log, kept in the private research record.
The last digits of `checker_value` can differ between machines because floating point is used to *propose* PSD shifts and
preconditioners; every accepted value is a valid bound and each one is on the safe side of the declared value.
