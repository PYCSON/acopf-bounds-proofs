# AI assistance — who did what

This work was carried out with AI systems under the direction of the human author, Puyi Sun, who chose the problem, set
the acceptance rules, decided what to keep or discard, operated the machines and is responsible for this release.
**The human author did not personally verify the mathematics line by line.** Correctness rests on the certificates and on
the checkers in `checkers/`, which anyone can re-run (`replay.py`).

| part | written by | run by |
|---|---|---|
| Certificate-generating programs (AC-OPF model, Krawczyk upper certificates, SDP / sparse SOS lower certificates, outward rounding) | Claude (Anthropic); one native SDP-assembly module by Codex (OpenAI) | Claude (Linux sandbox) for most certificates; Codex (Windows machine) for some lower certificates |
| Independent checkers in `checkers/` | GPT (OpenAI ChatGPT), without access to or imports from the generator code | Claude on Linux for all 52 endpoints; GPT in its sandbox for the small cases; Codex on Windows for many endpoints; the author's own computer (WSL2 Ubuntu) for all 52 endpoints, see replay_records/ |
| Review of the checkers' source | Claude (each checker version reviewed before its results were used) | — |
| Mathematical review of certificate formats | GPT and Claude, cross-reviewing each other | — |
| Tamper tests (deliberately corrupted certificates that must be rejected) | Claude and GPT | Claude; GPT |
| README, METHOD, replay script | Claude | — |
| Paper (`paper/`) | GPT; reviewed by Claude before release | — |

"Second-party verified" in this repository means: both endpoints of a case were re-proved by the GPT-authored checkers,
which were written independently of the generators. It does **not** mean external human peer review.
